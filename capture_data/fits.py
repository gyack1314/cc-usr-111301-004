"""FITS二进制表采集：按2880字节块遍历HDU，逐字节留痕，不猜测跨度。

只依赖标准库。主HDU（NAXIS=0）与BINTABLE扩展按FITS标准计算跨度；
不支持图像像素解码，不支持的扩展保留原始字节区间但不继续猜测。
所有读取经范围检查的读取器，记录每个字节区间的来源。
"""
import math
import re
import struct
from fractions import Fraction

CARD_SIZE = 80
BLOCK_SIZE = 2880

# BINTABLE原子类型：格式字母 -> (字节大小, struct字符, 数值类别)
# 数值类别：'int'整数原位型，'float'浮点原位型，'logical'逻辑，'char'字符
_ATOMS = {
    'L': (1, 'c', 'logical'),
    'B': (1, 'B', 'int'),
    'I': (2, 'h', 'int'),
    'J': (4, 'i', 'int'),
    'K': (8, 'q', 'int'),
    'E': (4, 'f', 'float'),
    'D': (8, 'd', 'float'),
    'A': (1, 's', 'char'),
}
_DESC_ATOMS = {'P': (4, 'i'), 'Q': (8, 'q')}
_TDIM_RE = re.compile(r'^\s*\((.*)\)\s*$')


class FITSError(ValueError):
    """结构性违约：头卡非法、跨度无法确定或必需值不符。"""


def _pad_to_block(size):
    return (size + BLOCK_SIZE - 1) // BLOCK_SIZE * BLOCK_SIZE


class _Reader:
    """范围检查读取器：任何越界立即失败，并记录实际字节来源。"""

    def __init__(self, data):
        self.data = data
        self.calls = []

    def read(self, start, length):
        if type(start) is not int or type(length) is not int:
            raise FITSError('偏移与长度必须为整数')
        if start < 0 or length < 0 or start + length > len(self.data):
            raise FITSError(f'字节范围越界: start={start} length={length} file={len(self.data)}')
        self.calls.append({'start': start, 'length': length})
        return self.data[start:start + length]

    def __len__(self):
        return len(self.data)


# ---------------------------------------------------------------- 头卡解析

class Card:
    """一张80字节卡片的解析结果与原始来源。"""

    __slots__ = ('index', 'keyword', 'value', 'raw', 'source', 'is_commentary')

    def __init__(self, index, keyword, value, raw, source, is_commentary):
        self.index = index
        self.keyword = keyword
        self.value = value
        self.raw = raw
        self.source = source
        self.is_commentary = is_commentary


def _is_printable_ascii(b):
    return all(0x20 <= c <= 0x7E for c in b)


def _parse_value_field(body):
    """解析卡片第10..79列。body为69字节bytes（已保证可打印ASCII）。

    字符串单引号内两个单引号转义一个字面引号；只有引号外的'/'是注释起点。
    返回 (value, kind)，kind ∈ {'string','logical','number','bare'}。
    """
    text = body.decode('ascii')
    stripped = text.lstrip(' ')
    if stripped.startswith("'"):
        # 字符串从首个引号开始，前导空白不属于字符串内容
        lead = len(text) - len(stripped)
        chars = []
        i = lead + 1
        while i < len(text):
            ch = text[i]
            if ch == "'":
                # 连续两个单引号 -> 一个字面单引号；否则字符串闭合
                if i + 1 < len(text) and text[i + 1] == "'":
                    chars.append("'")
                    i += 2
                    continue
                tail = text[i + 1:].strip()
                if tail and not tail.startswith('/'):
                    raise FITSError('闭合引号后只允许空白或注释')
                return ''.join(chars), 'string'
            chars.append(ch)
            i += 1
        raise FITSError("字符串卡片缺少闭合单引号")
    # 非字符串：引号外的第一个'/'开始注释
    slash = text.find('/')
    token = text if slash < 0 else text[:slash]
    token = token.strip()
    if token == '':
        return None, 'bare'
    if token == 'T':
        return True, 'logical'
    if token == 'F':
        return False, 'logical'
    # 数字：FITS允许D/d作指数；保留精确有理值
    candidate = token.replace('D', 'E').replace('d', 'e')
    try:
        return Fraction(candidate), 'number'
    except ValueError:
        raise FITSError(f'无法解析的卡片值: {token!r}')


def parse_header(reader, header_start):
    """从 header_start 解析一个头块，返回 (cards, header_length, padding)。

    头块由若干2880字节块组成，END卡片后到块尾为空格padding。
    """
    if header_start % BLOCK_SIZE != 0:
        raise FITSError(f'头起点不在2880字节边界: {header_start}')
    cards = []
    end_seen = False
    block_count = 0
    end_card_index = None
    while True:
        block = reader.read(header_start + block_count * BLOCK_SIZE, BLOCK_SIZE)
        if not _is_printable_ascii(block):
            raise FITSError(f'头块含非ASCII字节 @ {header_start + block_count * BLOCK_SIZE}')
        for k in range(BLOCK_SIZE // CARD_SIZE):
            card_bytes = block[k * CARD_SIZE:(k + 1) * CARD_SIZE]
            card_index = block_count * (BLOCK_SIZE // CARD_SIZE) + k
            source = {'start': header_start + card_index * CARD_SIZE,
                      'length': CARD_SIZE, 'hdu_relative': card_index * CARD_SIZE}
            if end_seen:
                # END之后只允许全空格卡片（头padding）
                if card_bytes != b' ' * CARD_SIZE:
                    raise FITSError('END卡片后存在非空格头padding')
                continue
            if not _is_printable_ascii(card_bytes):
                raise FITSError('卡片含非ASCII字节')
            keyword = card_bytes[:8].decode('ascii').strip()
            if keyword == '':
                # 空白卡：保留，不取值
                cards.append(Card(card_index, None, None, card_bytes, source, True))
                continue
            if keyword in ('COMMENT', 'HISTORY'):
                # 重复的COMMENT/HISTORY逐条保留，不合并、不丢弃
                cards.append(Card(card_index, keyword,
                                  card_bytes[8:].decode('ascii').rstrip(),
                                  card_bytes, source, True))
                continue
            if keyword == 'END':
                if card_bytes[8:].strip():
                    raise FITSError('END卡片不允许带值')
                end_seen = True
                end_card_index = card_index
                cards.append(Card(card_index, 'END', None, card_bytes, source, False))
                continue
            # 普通键值卡：第9、10列须为 '= '
            if card_bytes[8:10] != b'= ':
                raise FITSError(f'卡片 {keyword!r} 缺少值指示符 "= "')
            value, kind = _parse_value_field(card_bytes[10:])
            cards.append(Card(card_index, keyword, (value, kind), card_bytes,
                              source, False))
        block_count += 1
        if end_seen:
            break
        if header_start + block_count * BLOCK_SIZE > len(reader):
            raise FITSError('头块未遇到END即到达文件尾')
    header_length = block_count * BLOCK_SIZE
    padding = header_length - (end_card_index + 1) * CARD_SIZE
    return cards, header_length, padding


def _keyed_cards(cards):
    """普通键值卡 -> {keyword: Card}；重复键（非注释类）即违约。"""
    keyed = {}
    for c in cards:
        if c.is_commentary or c.keyword == 'END':
            continue
        if c.keyword in keyed:
            raise FITSError(f'关键字 {c.keyword!r} 在同一头中重复出现')
        keyed[c.keyword] = c
    return keyed


def _require(keyed, keyword):
    c = keyed.get(keyword)
    if c is None:
        raise FITSError(f'缺少必需卡片 {keyword}')
    return c


def _int_value(card):
    value, kind = card.value
    if kind != 'number' or value.denominator != 1:
        raise FITSError(f'卡片 {card.keyword} 必须为整数')
    return int(value)


def _bool_value(card):
    value, kind = card.value
    if kind != 'logical':
        raise FITSError(f'卡片 {card.keyword} 必须为逻辑值 T/F')
    return value


def _string_value(card):
    value, kind = card.value
    if kind != 'string':
        raise FITSError(f'卡片 {card.keyword} 必须为字符串')
    return value


# ---------------------------------------------------------------- HDU遍历

class HDU:
    """一个HDU的结构骨架与原始字节区间（原始HDU始终保留）。"""

    __slots__ = ('index', 'kind', 'header_start', 'header_length',
                 'header_padding', 'data_start', 'data_size', 'data_padding',
                 'total_span', 'cards', 'keyed', 'raw_source')

    def __init__(self, **kw):
        for name in self.__slots__:
            setattr(self, name, kw.get(name))


def walk_hdus(data):
    """严格按2880字节块顺序遍历全部HDU。

    绝不扫描字节找下一个头：下一个HDU起点完全由本HDU跨度算术决定，
    因此END后的头padding与数据区尾部padding都不会被误当成新头。
    """
    reader = _Reader(data)
    hdus = []
    offset = 0
    index = 0
    while offset < len(reader):
        if len(reader) - offset < BLOCK_SIZE:
            raise FITSError(f'文件尾部剩余 {len(reader) - offset} 字节，不足一个头块')
        cards, header_length, header_padding = parse_header(reader, offset)
        keyed = _keyed_cards(cards)
        data_start = offset + header_length

        if index == 0:
            simple = _require(keyed, 'SIMPLE')
            if not _bool_value(simple):
                raise FITSError('主HDU的SIMPLE必须为T')
        else:
            xt = _require(keyed, 'XTENSION')
            if _string_value(xt) != 'BINTABLE':
                # 不支持的扩展：必需卡仍须合法才能算出跨度，保留原区间，
                # 但绝不靠猜长度继续。
                hdu = _build_unsupported(index, offset, header_length,
                                         header_padding, data_start, cards, keyed,
                                         reader, xt.value[0])
                hdus.append(hdu)
                offset += hdu.total_span
                index += 1
                continue

        bitpix_c = _require(keyed, 'BITPIX')
        bitpix = _int_value(bitpix_c)
        if bitpix not in (8, 16, 32, 64, -32, -64):
            raise FITSError(f'BITPIX非法: {bitpix}')
        naxis = _int_value(_require(keyed, 'NAXIS'))
        if naxis < 0:
            raise FITSError('NAXIS不能为负')
        axis_sizes = []
        for n in range(1, naxis + 1):
            axis_sizes.append(_int_value(_require(keyed, f'NAXIS{n}')))
            if axis_sizes[-1] < 0:
                raise FITSError(f'NAXIS{n}不能为负')
        pcount = _int_value(keyed['PCOUNT']) if 'PCOUNT' in keyed else 0
        gcount = _int_value(keyed['GCOUNT']) if 'GCOUNT' in keyed else 1
        if pcount < 0 or gcount < 0:
            raise FITSError('PCOUNT/GCOUNT不能为负')

        is_bintable = index > 0 and _string_value(keyed['XTENSION']) == 'BINTABLE'
        if is_bintable:
            if naxis != 2:
                raise FITSError('BINTABLE的NAXIS必须为2')
            if bitpix != 8:
                raise FITSError('BINTABLE的BITPIX必须为8')
            if gcount != 1:
                raise FITSError('BINTABLE的GCOUNT必须为1')

        # FITS数据跨度：|BITPIX|/8 * 积(NAXISn) * GCOUNT + PCOUNT
        if naxis == 0:
            data_size = 0
            if pcount != 0 or gcount != 1:
                raise FITSError('NAXIS=0的HDU要求PCOUNT=0、GCOUNT=1')
        else:
            element_count = gcount
            for a in axis_sizes:
                element_count = _checked_mul(element_count, a, 'NAXIS乘积')
            data_size = _checked_mul(abs(bitpix) // 8, element_count, '数据区大小')
            data_size = _checked_add(data_size, pcount, 'PCOUNT')

        # 数据块向2880取整；零长数据区不占数据块，但整个HDU至少有一个头块。
        # 下一个HDU起点由该跨度算术决定，绝不扫描字节找头。
        data_padding = (_pad_to_block(data_size) - data_size
                        if data_size > 0 else 0)
        total_span = header_length + data_size + data_padding
        # 先读数据区（含padding）确认整体不越界，原始HDU区间随之保留
        reader.read(offset, total_span)
        hdu = HDU(
            index=index,
            kind='bintable' if is_bintable else ('primary' if index == 0 else 'image'),
            header_start=offset, header_length=header_length,
            header_padding=header_padding, data_start=data_start,
            data_size=data_size, data_padding=data_padding,
            total_span=total_span, cards=cards, keyed=keyed,
            raw_source={'start': offset, 'length': total_span})
        hdus.append(hdu)
        offset += total_span
        index += 1
    if offset != len(reader):
        raise FITSError('HDU跨度合计与文件长度不一致')
    return hdus, reader


def _build_unsupported(index, offset, header_length, header_padding,
                       data_start, cards, keyed, reader, ext_name):
    """非BINTABLE扩展：仅当跨度可由标准必需卡确定时保留原区间。"""
    bitpix = _int_value(_require(keyed, 'BITPIX'))
    if bitpix not in (8, 16, 32, 64, -32, -64):
        raise FITSError(f'BITPIX非法: {bitpix}')
    naxis = _int_value(_require(keyed, 'NAXIS'))
    axis_sizes = [_int_value(_require(keyed, f'NAXIS{n}')) for n in range(1, naxis + 1)]
    pcount = _int_value(keyed['PCOUNT']) if 'PCOUNT' in keyed else 0
    gcount = _int_value(keyed['GCOUNT']) if 'GCOUNT' in keyed else 1
    if naxis == 0:
        data_size = 0
    else:
        element_count = gcount
        for a in axis_sizes:
            element_count = _checked_mul(element_count, a, 'NAXIS乘积')
        data_size = _checked_mul(abs(bitpix) // 8, element_count, '数据区大小')
        data_size = _checked_add(data_size, pcount, 'PCOUNT')
    data_padding = (_pad_to_block(data_size) - data_size
                    if data_size > 0 else 0)
    total_span = header_length + data_size + data_padding
    reader.read(offset, total_span)
    return HDU(index=index, kind=f'unsupported:{ext_name}',
               header_start=offset, header_length=header_length,
               header_padding=header_padding, data_start=data_start,
               data_size=data_size, data_padding=data_padding,
               total_span=total_span, cards=cards, keyed=keyed,
               raw_source={'start': offset, 'length': total_span})


def _checked_mul(a, b, what):
    if a < 0 or b < 0:
        raise FITSError(f'{what}含负值')
    result = a * b
    if b and result // b != a:
        raise FITSError(f'{what}乘法溢出')
    return result


def _checked_add(a, b, what):
    if a < 0 or b < 0:
        raise FITSError(f'{what}含负值')
    return a + b


# ---------------------------------------------------------------- 列布局

class Column:
    """一列的TFORM布局、换算规则与字节来源。"""

    __slots__ = ('number', 'name', 'format_letter', 'repeat', 'atom_size',
                 'atom_code', 'value_kind', 'is_descriptor', 'descriptor_size',
                 'offset', 'width', 'tnull', 'tscal', 'tzero', 'tdim',
                 'form_card_source', 'selected')

    def __init__(self, **kw):
        for name in self.__slots__:
            setattr(self, name, kw.get(name))


def _parse_tform(text):
    """解析TFORM值，返回 (repeat, letter, heap_atom)。

    普通列形如 '10J'、'20A'、'1D'；变长列形如 '1PD'/'1QD'，
    P/Q后必须跟一个堆元素原子字母。不支持变长字符串（PA/QA）。
    """
    s = text.strip()
    m = re.match(r'^(\d+)([A-Za-z])([A-Za-z]?)$', s)
    if not m:
        raise FITSError(f'TFORM无法解析: {text!r}')
    repeat_s, letter, tail = m.groups()
    repeat = int(repeat_s)
    if repeat < 1:
        raise FITSError('TFORM重复数必须为正')
    if letter in ('X', 'C', 'M'):
        raise FITSError(f'本项目不支持的列类型 {letter}（bit/复数列）')
    if letter in _DESC_ATOMS:
        if repeat != 1:
            raise FITSError('变长描述符列的TFORM重复数必须为1')
        if tail == '':
            raise FITSError(f'变长列 {text!r} 缺少堆元素原子类型')
        if tail not in _ATOMS or tail == 'A':
            raise FITSError(f'变长列 {text!r} 的堆原子类型不支持: {tail}（不支持变长字符串）')
        return repeat, letter, tail
    if tail:
        raise FITSError(f'TFORM含无法识别的尾部: {text!r}')
    if letter not in _ATOMS:
        raise FITSError(f'不支持的列类型: {letter}')
    return repeat, letter, None


def build_layout(hdu):
    """从BINTABLE头卡构建列布局，并核对NAXIS1行宽与列布局一致。"""
    keyed = hdu.keyed
    tfields = _int_value(_require(keyed, 'TFIELDS'))
    if tfields < 0:
        raise FITSError('TFIELDS不能为负')
    naxis1 = _int_value(_require(keyed, 'NAXIS1'))
    naxis2 = _int_value(_require(keyed, 'NAXIS2'))
    if naxis1 < 0 or naxis2 < 0:
        raise FITSError('NAXIS1/NAXIS2不能为负')

    columns = []
    cursor = 0
    names = set()
    for n in range(1, tfields + 1):
        # TFIELDS要求 1..TFIELDS 的TTYPE/TFORM连续存在，缺卡即结构违约
        if f'TTYPE{n}' not in keyed:
            raise FITSError(f'缺少列 {n} 的TTYPE{n}卡片')
        if f'TFORM{n}' not in keyed:
            raise FITSError(f'缺少列 {n} 的TFORM{n}卡片')
        ttype_card = _require(keyed, f'TTYPE{n}')
        name = _string_value(ttype_card)
        if name in names:
            raise FITSError(f'列名重复: {name!r}')
        names.add(name)
        form_card = _require(keyed, f'TFORM{n}')
        form_text = _string_value(form_card)
        repeat, letter, heap_atom = _parse_tform(form_text)

        is_desc = letter in _DESC_ATOMS
        if is_desc:
            desc_width = _DESC_ATOMS[letter][0]
            atom_size, atom_code, value_kind = _ATOMS[heap_atom]
            width = desc_width * 2
            col = Column(number=n, name=name, format_letter=letter,
                         repeat=repeat, atom_size=atom_size, atom_code=atom_code,
                         value_kind=value_kind,
                         is_descriptor=True, descriptor_size=desc_width,
                         offset=cursor, width=width,
                         form_card_source=dict(form_card.source))
        else:
            atom_size, atom_code, value_kind = _ATOMS[letter]
            width = _checked_mul(atom_size, repeat, f'列 {name} 宽度')
            col = Column(number=n, name=name, format_letter=letter,
                         repeat=repeat, atom_size=atom_size, atom_code=atom_code,
                         value_kind=value_kind,
                         is_descriptor=False, descriptor_size=None,
                         offset=cursor, width=width,
                         form_card_source=dict(form_card.source))

        # TNULL（只允许整数列）
        if f'TNULL{n}' in keyed:
            tnull_card = keyed[f'TNULL{n}']
            if col.value_kind != 'int':
                raise FITSError(f'TNULL{n} 只允许用于整数列，列 {name} 不是')
            col.tnull = _int_value(tnull_card)
        else:
            col.tnull = None
        # TSCAL/TZERO：整数列要求精确有理；缺省 1 / 0
        tscal = Fraction(1)
        tzero = Fraction(0)
        if f'TSCAL{n}' in keyed:
            c = keyed[f'TSCAL{n}']
            if c.value[1] != 'number':
                raise FITSError(f'TSCAL{n} 必须为数值')
            tscal = c.value[0]
        if f'TZERO{n}' in keyed:
            c = keyed[f'TZERO{n}']
            if c.value[1] != 'number':
                raise FITSError(f'TZERO{n} 必须为数值')
            tzero = c.value[0]
        col.tscal = tscal
        col.tzero = tzero
        # TDIM：只对固定长度数组生效；维度积必须等于重复数
        col.tdim = None
        if f'TDIM{n}' in keyed:
            dim_text = _string_value(keyed[f'TDIM{n}'])
            m2 = _TDIM_RE.match(dim_text)
            if not m2:
                raise FITSError(f'TDIM{n} 格式非法: {dim_text!r}')
            inner = m2.group(1).strip()
            dims = tuple(int(x.strip()) for x in inner.split(',')) if inner else ()
            if any(d < 1 for d in dims):
                raise FITSError(f'TDIM{n} 维度必须为正')
            if not is_desc:
                product = 1
                for d in dims:
                    product = _checked_mul(product, d, 'TDIM维度积')
                if product != repeat:
                    raise FITSError(
                        f'TDIM{n} 维度积 {product} 与TFORM重复数 {repeat} 不一致')
            col.tdim = dims  # FITS次序：第一维变化最快
        columns.append(col)
        cursor = _checked_add(cursor, width, f'列 {name} 偏移')

    # FITS允许NAXIS1 >= 列宽合计：差值是每行的未用填充字节，
    # 属于行的一部分，绝不能被当成下一列或下一行的数据。
    if cursor > naxis1:
        raise FITSError(
            f'列布局宽度合计 {cursor} 超过NAXIS1行宽 {naxis1}')
    row_padding = naxis1 - cursor
    row_area = _checked_mul(naxis1, naxis2, '固定行区')
    return columns, naxis1, naxis2, row_area, row_padding


# ---------------------------------------------------------------- 值解码

def _fraction_to_number(frac):
    """整数换算的精确有理值：整数给int，否则给分子/分母证据。"""
    if frac.denominator == 1:
        return int(frac), 'integer'
    return {'numerator': frac.numerator, 'denominator': frac.denominator,
            'decimal': f'{frac.numerator / frac.denominator:.17g}'}, 'rational'


def _decode_int(buf, code):
    (raw,) = struct.unpack('>' + code, buf)
    return raw


def _decode_float(buf, code):
    (raw,) = struct.unpack('>' + code, buf)
    if math.isnan(raw):
        return None, True  # 浮点NaN作为缺测，绝不转成字符串"nan"
    return raw, False


def _decode_logical(buf):
    b = buf
    if b == b'T':
        return True, False
    if b == b'F':
        return False, False
    if b == b' ' or b == b'\x00':
        return None, True  # 空缺
    raise FITSError(f'L列非法逻辑字节: {b!r}')


def _scaled_int(raw, col):
    """整数：先在原值域判TNULL空，再做TSCAL*x+TZERO精确换算。"""
    missing = col.tnull is not None and raw == col.tnull
    result = {'raw_integer': raw, 'missing': missing,
              'scaled': None, 'scaled_kind': None}
    if missing:
        return result
    physical = col.tscal * raw + col.tzero
    value, kind = _fraction_to_number(physical)
    result['scaled'] = value
    result['scaled_kind'] = kind
    return result


def _scaled_float(raw, missing, col):
    """浮点：保留原位型；有穷数再给换算值，缺测不换算。"""
    result = {'raw_float': raw, 'missing': missing,
              'scaled': None, 'scaled_kind': None}
    if missing:
        return result
    physical = float(col.tscal) * raw + float(col.tzero)
    if math.isnan(physical):
        result['missing'] = True
        return result
    result['scaled'] = physical
    result['scaled_kind'] = 'float'
    return result


def _source(reader_relative_start, length, hdu_data_start):
    """给出文件绝对位置与HDU数据区相对位置两套坐标。"""
    return {'start': hdu_data_start + reader_relative_start,
            'length': length,
            'hdu_data_relative': reader_relative_start}


# ---------------------------------------------------------------- 采集

def collect(data, hdu_index, columns, max_array_elements=None):
    """采集指定BINTABLE HDU的请求列，冻结选择并逐字段留字节来源。

    参数
    ----
    data: 完整FITS字节
    hdu_index: 目标HDU序号（主HDU为0）
    columns: 请求列名序列，进入时立即冻结为元组
    max_array_elements: 变长数组元素数上界（结构检查用）
    """
    requested = tuple(columns)
    if len(set(requested)) != len(requested):
        raise FITSError('请求列存在重复')
    hdus, reader = walk_hdus(data)
    if not 0 <= hdu_index < len(hdus):
        raise FITSError(f'HDU序号 {hdu_index} 越界，共 {len(hdus)} 个HDU')
    target = hdus[hdu_index]
    if target.kind != 'bintable':
        raise FITSError(f'HDU {hdu_index} 不是BINTABLE（类型 {target.kind}），不采集')

    layout_cols, naxis1, naxis2, row_area, row_padding = build_layout(target)
    by_name = {c.name: c for c in layout_cols}
    for name in requested:
        if name not in by_name:
            raise FITSError(f'请求列 {name!r} 不存在于HDU {hdu_index}')
    for c in layout_cols:
        c.selected = c.name in requested

    keyed = target.keyed
    pcount = _int_value(keyed['PCOUNT'])
    # THEAP缺省紧接固定行区；显式THEAP与行区之间的间隙属于PCOUNT
    if 'THEAP' in keyed:
        theap = _int_value(keyed['THEAP'])
        theap_explicit = True
    else:
        theap = row_area
        theap_explicit = False
    if theap < row_area:
        raise FITSError(f'THEAP={theap} 早于固定行区结尾 {row_area}')
    gap_size = theap - row_area
    heap_size = pcount - gap_size
    if heap_size < 0:
        raise FITSError('PCOUNT不足以容纳THEAP间隙与堆')
    if theap > row_area + pcount:
        raise FITSError(
            f'THEAP={theap} 超出数据区末尾 {row_area + pcount}')
    data_end = row_area + pcount
    if data_end != target.data_size:
        raise FITSError('行区+PCOUNT与HDU数据跨度不一致')

    issues = []
    layout_evidence = _layout_evidence(target, layout_cols, naxis1, naxis2,
                                       row_area, row_padding, pcount, theap,
                                       theap_explicit, gap_size, heap_size)

    # 描述符先全部结构性核验（无论该列是否被选），但不读未选堆列的负载
    descriptor_map = {}  # (col_number, row) -> 描述符证据
    for col in layout_cols:
        if not col.is_descriptor:
            continue
        for row in range(naxis2):
            desc = _read_descriptor(reader, target, col, row, naxis1)
            count, offset = desc['count'], desc['offset']
            negative = count < 0 or offset < 0
            # 负值不进入算术：先标记，窗口长度按0处理避免乘法/加法抛错
            if negative:
                elem_bytes = 0
                window_end = offset  # 必然 < 0 或配合负count判负
                within = False
            else:
                elem_bytes = _checked_mul(count, col.atom_size, '堆元素字节数')
                window_end = _checked_add(offset, elem_bytes, '堆窗口')
                within = window_end <= heap_size
            desc['heap_relative_window'] = {'start': offset, 'length': elem_bytes}
            desc['data_relative_window'] = {
                'start': theap + offset if offset >= 0 else offset,
                'length': elem_bytes}
            desc['within_heap'] = within
            if (not negative and max_array_elements is not None
                    and count > max_array_elements):
                desc['over_max_elements'] = True
                issues.append({'kind': 'array_too_large', 'column': col.name,
                               'row': row, 'count': count,
                               'limit': max_array_elements})
            if not within:
                issues.append({'kind': 'heap_out_of_range', 'column': col.name,
                               'row': row, 'offset': offset, 'count': count,
                               'heap_size': heap_size,
                               'negative': negative})
            descriptor_map[(col.number, row)] = desc

    # 行采集：只对请求列解码负载
    records = []
    for row in range(naxis2):
        record = {'row_index': row,
                  'row_source': _source(row * naxis1, naxis1, target.data_start),
                  'fields': {}}
        for col in layout_cols:
            if not col.selected:
                continue
            if col.is_descriptor:
                record['fields'][col.name] = _collect_descriptor(
                    reader, target, col, row, naxis1, theap, heap_size,
                    descriptor_map[(col.number, row)])
            else:
                record['fields'][col.name] = _collect_fixed(
                    reader, target, col, row, naxis1)
        records.append(record)

    return {
        'hdu_index': hdu_index,
        'selected_columns': list(requested),
        'max_array_elements': max_array_elements,
        'hdu': _hdu_evidence(target),
        'layout': layout_evidence,
        'records': records,
        'issues': issues,
        'structure_valid': not issues,
        'read_calls': list(reader.calls),
        'raw_hdus': [{'hdu_index': h.index, 'kind': h.kind,
                      'source': dict(h.raw_source)} for h in hdus],
    }


def _read_descriptor(reader, hdu, col, row, naxis1):
    """读P/Q描述符（两个有符号整数），校验负值；不触碰堆负载。"""
    cell = col.offset
    desc_start = row * naxis1 + cell
    code = 'i' if col.descriptor_size == 4 else 'q'
    buf = reader.read(hdu.data_start + desc_start, col.descriptor_size * 2)
    count, offset = struct.unpack('>' + code * 2, buf)
    evidence = {
        'column': col.name, 'row': row,
        'descriptor_kind': 'P' if col.descriptor_size == 4 else 'Q',
        'count': count, 'offset': offset,
        'descriptor_source': _source(desc_start, col.descriptor_size * 2,
                                     hdu.data_start),
        'count_source': _source(desc_start, col.descriptor_size, hdu.data_start),
        'offset_source': _source(desc_start + col.descriptor_size,
                                 col.descriptor_size, hdu.data_start),
        'negative_count': count < 0,
        'negative_offset': offset < 0,
    }
    if count < 0 or offset < 0:
        evidence['within_heap'] = False
    return evidence


def _collect_descriptor(reader, hdu, col, row, naxis1, theap, heap_size, desc):
    """组装一个变长字段：描述符证据 + 堆切片证据（共享不报错）。"""
    out = {
        'column': col.name, 'row': row,
        'value_kind': col.value_kind,
        'descriptor': {k: desc[k] for k in (
            'descriptor_kind', 'count', 'offset', 'descriptor_source',
            'count_source', 'offset_source', 'negative_count',
            'negative_offset', 'heap_relative_window',
            'data_relative_window', 'within_heap')},
        'elements': None,
        'heap_slice': None,
        'readable': False,
    }
    if desc['negative_count'] or desc['negative_offset'] or not desc['within_heap']:
        # 越界/负值：不给负载，绝不跨HDU读取
        out['error'] = 'heap_out_of_range'
        return out
    count, offset = desc['count'], desc['offset']
    elem_bytes = col.atom_size * count
    heap_abs = hdu.data_start + theap + offset
    payload = reader.read(heap_abs, elem_bytes)
    out['readable'] = True
    out['heap_slice'] = {
        'start': heap_abs, 'length': elem_bytes,
        'hdu_heap_relative': offset,
        'hdu_data_relative': theap + offset,
        'raw_hex': payload.hex(),
    }
    elements = []
    for k in range(count):
        buf = payload[k * col.atom_size:(k + 1) * col.atom_size]
        elem_src = {'start': heap_abs + k * col.atom_size,
                    'length': col.atom_size,
                    'hdu_heap_relative': offset + k * col.atom_size}
        elements.append(_decode_atom(col, buf, elem_src))
    out['elements'] = elements
    return out


def _collect_fixed(reader, hdu, col, row, naxis1):
    """固定长度列：按重复数逐元素解码，每个元素带独立字节来源。

    A字符列是特例：rA是一个长度r的单字符串单元（repeat是字符数，
    不是元素数），整体解码为一个元素并保留原字节。
    """
    cell_start = row * naxis1 + col.offset
    buf = reader.read(hdu.data_start + cell_start, col.width)
    out = {
        'column': col.name, 'row': row,
        'value_kind': col.value_kind,
        'format_letter': col.format_letter,
        'repeat': col.repeat,
        'tdim': col.tdim,
        'cell_source': _source(cell_start, col.width, hdu.data_start),
        'raw_hex': buf.hex(),
        'elements': [],
        'value': None,
    }
    if col.value_kind == 'char':
        elem_src = _source(cell_start, col.width, hdu.data_start)
        element = _decode_atom(col, buf, elem_src)
        out['elements'] = [element]
        out['value'] = element
        return out
    for k in range(col.repeat):
        piece = buf[k * col.atom_size:(k + 1) * col.atom_size]
        elem_src = _source(cell_start + k * col.atom_size, col.atom_size,
                           hdu.data_start)
        out['elements'].append(_decode_atom(col, piece, elem_src))
    if col.repeat == 1:
        out['value'] = out['elements'][0]
    else:
        out['value'] = out['elements']
    return out


def _decode_atom(col, buf, source):
    """按列原子类型解码一个元素，同时保留原位型与换算值。"""
    base = {'source': source, 'raw_hex': buf.hex(), 'missing': False}
    if col.value_kind == 'int':
        raw = _decode_int(buf, col.atom_code)
        base.update(_scaled_int(raw, col))
    elif col.value_kind == 'float':
        raw, missing = _decode_float(buf, col.atom_code)
        base['raw_hex'] = buf.hex()
        base.update(_scaled_float(raw, missing, col))
    elif col.value_kind == 'logical':
        value, missing = _decode_logical(buf)
        base['raw_logical'] = None if missing else value
        base['missing'] = missing
    elif col.value_kind == 'char':
        # A列：保留原始字节，按约定去右空格给出字符串；不删内部/前导空格
        text = buf.decode('ascii')
        base['raw_text'] = text
        base['string'] = text.rstrip(' ')
        base['trailing_spaces'] = len(text) - len(base['string'])
    else:  # pragma: no cover - 布局阶段已拦截
        raise FITSError(f'未知原子类型: {col.value_kind}')
    return base


def _layout_evidence(hdu, columns, naxis1, naxis2, row_area, row_padding,
                     pcount, theap, theap_explicit, gap_size, heap_size):
    return {
        'naxis1': naxis1, 'naxis2': naxis2,
        'row_area_size': row_area,
        'row_end_padding_per_row': row_padding,
        'columns': [{
            'number': c.number, 'name': c.name,
            'tform_repeat': c.repeat, 'format_letter': c.format_letter,
            'is_variable': c.is_descriptor,
            'row_offset': c.offset, 'width': c.width,
            'atom_size': c.atom_size, 'value_kind': c.value_kind,
            'tdim': c.tdim, 'selected': c.selected,
            'has_tnull': c.tnull is not None,
            'tscal': str(c.tscal), 'tzero': str(c.tzero),
            'form_card_source': c.form_card_source,
        } for c in columns],
        'pcount': pcount,
        'theap': theap, 'theap_explicit': theap_explicit,
        'theap_gap_size': gap_size, 'heap_size': heap_size,
        'layout_width_matches_naxis1':
            sum(c.width for c in columns) + row_padding == naxis1,
    }


def _hdu_evidence(hdu):
    return {
        'index': hdu.index, 'kind': hdu.kind,
        'header_start': hdu.header_start,
        'header_length': hdu.header_length,
        'header_padding_after_end': hdu.header_padding,
        'data_start': hdu.data_start,
        'data_size': hdu.data_size,
        'data_trailing_padding': hdu.data_padding,
        'total_span': hdu.total_span,
        'raw_source': dict(hdu.raw_source),
        'commentary': [
            {'card_index': c.index, 'keyword': c.keyword, 'text': c.value,
             'source': c.source}
            for c in hdu.cards if c.is_commentary and c.keyword is not None
        ],
    }
