"""独立FITS夹具构造器：只用原始bytes/struct手工拼卡与数据。

与 capture_data.fits 完全不共享解析代码，作为测试中的独立原值表来源。
所有数值大端；头80字节卡片，头块与数据块各自向2880取整。
"""
import struct

BLOCK = 2880
CARD = 80


def card(keyword, value=None, comment=None):
    """构造一张80字节卡，字符串/整数/逻辑/注释卡各自按FITS排版。"""
    if keyword in ('COMMENT', 'HISTORY'):
        body = (keyword.ljust(8) + str(value))[:CARD]
        return body.ljust(CARD).encode('ascii')
    if keyword == 'END':
        return b'END'.ljust(CARD)
    if keyword == '':
        return b' ' * CARD
    head = keyword.ljust(8) + '= '
    if isinstance(value, bool):
        v = 'T' if value else 'F'
        field = v.rjust(20)
    elif isinstance(value, int):
        field = str(value).rjust(20)
    elif isinstance(value, float):
        # 数值卡：不带引号，保留足够精度
        field = repr(value).rjust(20)
    elif isinstance(value, str):
        field = ("'" + value + "'").ljust(20)
    else:
        raise TypeError(f'不支持的卡片值类型: {type(value)}')
    text = head + field
    if comment:
        text = (' ' + text + ' /' + comment)[:CARD]
    return text.ljust(CARD).encode('ascii')


def header_block(cards):
    raw = b''.join(cards) + card('END')
    if len(raw) % BLOCK:
        raw += b' ' * (BLOCK - len(raw) % BLOCK)
    return raw


def pad_data(data):
    if len(data) % BLOCK:
        data += b'\x00' * (BLOCK - len(data) % BLOCK)
    return data


def primary_header():
    return header_block([
        card('SIMPLE', True),
        card('BITPIX', 8),
        card('NAXIS', 0),
        card('EXTEND', True),
    ])


def bintable(header_cards, data_bytes):
    """主HDU + BINTABLE扩展拼成完整文件。"""
    return primary_header() + header_block(header_cards) + pad_data(data_bytes)


def std_table_cards(naxis1, naxis2, pcount, theap=None, extra=(),
                    spectrum_tform='1PD'):
    """标准三列表 ID(J) SPECTRUM(PD/QD) QUALITY(I) 的头卡。"""
    cards = [
        card('XTENSION', 'BINTABLE'),
        card('BITPIX', 8),
        card('NAXIS', 2),
        card('NAXIS1', naxis1),
        card('NAXIS2', naxis2),
        card('PCOUNT', pcount),
        card('GCOUNT', 1),
        card('TFIELDS', 3),
        card('TTYPE1', 'ID'),
        card('TFORM1', '1J'),
        card('TTYPE2', 'SPECTRUM'),
        card('TFORM2', spectrum_tform),
        card('TTYPE3', 'QUALITY'),
        card('TFORM3', '1I'),
        card('TNULL3', -999),
        card('TSCAL3', 2),
        card('TZERO3', 1),
    ]
    if theap is not None:
        cards.append(card('THEAP', theap))
    cards.extend(extra)
    return cards


# ------------------------------------------------------------- 各场景夹具

def build_shared_heap():
    """两行描述符共享同一堆片段（P描述符，非默认THEAP，带间隙）。"""
    naxis1, naxis2, pcount, theap = 14, 2, 24, 36
    cards = std_table_cards(naxis1, naxis2, pcount, theap,
                            extra=[card('OBJECT', 'orbit/a')])
    rows = b''
    rows += struct.pack('>i', 1) + struct.pack('>ii', 2, 0) + struct.pack('>h', -999)
    rows += struct.pack('>i', 2) + struct.pack('>ii', 2, 0) + struct.pack('>h', 3)
    gap = b'\xAA' * (theap - naxis1 * naxis2)
    heap = struct.pack('>dd', 1.5, 2.5)
    heap += b'\x00' * (pcount - len(gap) - len(heap))
    return bintable(cards, rows + gap + heap), {
        'naxis1': naxis1, 'naxis2': naxis2, 'pcount': pcount, 'theap': theap,
        'rows': [
            {'id': 1, 'count': 2, 'offset': 0, 'quality_raw': -999,
             'quality_missing': True, 'quality_scaled': None},
            {'id': 2, 'count': 2, 'offset': 0, 'quality_raw': 3,
             'quality_missing': False, 'quality_scaled': 7},
        ],
        'heap_values': [1.5, 2.5],
    }


def build_q_out_of_range():
    """Q描述符（i64对），第二行offset越界。"""
    naxis1, naxis2, pcount, theap = 22, 2, 24, 52
    cards = std_table_cards(naxis1, naxis2, pcount, theap,
                            spectrum_tform='1QD')
    rows = b''
    rows += struct.pack('>i', 1) + struct.pack('>qq', 2, 0) + struct.pack('>h', -999)
    rows += struct.pack('>i', 2) + struct.pack('>qq', 2, 10000) + struct.pack('>h', 3)
    gap = b'\xBB' * (theap - naxis1 * naxis2)
    heap = struct.pack('>dd', 1.5, 2.5)
    heap += b'\x00' * (pcount - len(gap) - len(heap))
    return bintable(cards, rows + gap + heap), {
        'naxis1': naxis1, 'theap': theap,
        'bad_row': 1, 'bad_offset': 10000, 'bad_count': 2,
        'heap_size': pcount - (theap - naxis1 * naxis2),
    }


def build_slash_string():
    """A列数据内含斜杠与字面单引号（数据区字节，按约定去右空格）。"""
    label = b"vega/lyra's tar"   # 15字节，内含 / 与字面单引号
    assert len(label) == 15
    cell0 = label + b' ' * 2     # 17字节单元，2个右空格
    cell1 = b'spaced label' + b' ' * 5  # 17字节单元，5个右空格
    assert len(cell0) == len(cell1) == 17
    naxis1, naxis2 = 4 + 17, 2
    cards = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', naxis1), card('NAXIS2', naxis2),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 2),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
        # 列名本身含斜杠：引号内斜杠不得被当注释起点
        ("TTYPE2  = 'orbit/label'").ljust(80).encode('ascii'),
        card('TFORM2', '17A'),
    ]
    rows = b''
    rows += struct.pack('>i', 7) + cell0
    rows += struct.pack('>i', 8) + cell1
    return bintable(cards, rows), {
        'naxis1': naxis1,
        'column_name': 'orbit/label',
        'expected': [
            (7, "vega/lyra's tar", 2),
            (8, 'spaced label', 5),
        ],
    }


def build_row_padding():
    """NAXIS1大于列宽合计：行尾padding字节非零，且数据区尾部padding非零。

    列宽合计5，NAXIS1=7，每行2字节行尾填充；验证遍历不把padding当数据。
    """
    naxis1, naxis2, pcount = 7, 2, 6
    cards = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', naxis1), card('NAXIS2', naxis2),
        card('PCOUNT', pcount), card('GCOUNT', 1), card('TFIELDS', 2),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
        card('TTYPE2', 'FLAG'), card('TFORM2', '1L'),
    ]
    rows = b''
    rows += struct.pack('>i', 11) + b'T' + b'\xDE\xDE'
    rows += struct.pack('>i', 12) + b'F' + b'\xDE\xDE'
    heap = b'\xCC' * pcount
    return bintable(cards, rows + heap), {
        'naxis1': naxis1, 'naxis2': naxis2, 'pcount': pcount,
        'row_padding_per_row': 2,
        'expected': [(11, True), (12, False)],
    }


def build_type_garden():
    """覆盖 L B I J K E D A、TDIM固定数组、浮点NaN、整数TNULL、有理缩放。"""
    # 列：FLAG 1L(1) B 1B(1) I 1I(2) J 2J(8, TDIM=(2,1)) K 1K(8)
    #     E 1E(4, NaN) D 1D(8) A 4A(4) => 宽 1+1+2+8+8+4+8+4 = 36
    naxis1 = 36
    cards = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', naxis1), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 8),
        card('TTYPE1', 'FLAG'), card('TFORM1', '1L'),
        card('TTYPE2', 'BYTE'), card('TFORM2', '1B'),
        card('TTYPE3', 'SHORT'), card('TFORM3', '1I'),
        card('TNULL3', -1), card('TSCAL3', 1), card('TZERO3', 0),
        card('TTYPE4', 'VEC'), card('TFORM4', '2J'), card('TDIM4', '(2,1)'),
        card('TTYPE5', 'LONG'), card('TFORM5', '1K'),
        card('TSCAL5', 2), card('TZERO5', 1),
        card('TTYPE6', 'REAL'), card('TFORM6', '1E'),
        card('TTYPE7', 'DOUB'), card('TFORM7', '1D'),
        card('TSCAL7', 2),
        card('TTYPE8', 'TEXT'), card('TFORM8', '4A'),
    ]
    row = b''
    row += b' '                                   # FLAG 空缺
    row += struct.pack('>B', 200)                 # BYTE
    row += struct.pack('>h', -1)                  # SHORT = TNULL
    row += struct.pack('>ii', 5, -7)              # VEC 两元素
    row += struct.pack('>q', 3)                   # LONG -> 2*3+1=7
    row += struct.pack('>f', float('nan'))        # REAL NaN 缺测
    row += struct.pack('>d', 2.5)                 # DOUB -> 5.0
    row += b'abc '                                # TEXT -> 'abc'
    return bintable(cards, row), {
        'naxis1': naxis1,
        'expected': {
            'flag_missing': True,
            'byte_raw': 200,
            'short_missing': True,
            'vec_raw': [5, -7], 'tdim': (2, 1),
            'long_raw': 3, 'long_scaled': 7,
            'real_missing': True,
            'doub_raw': 2.5, 'doub_scaled': 5.0,
            'text': 'abc', 'text_trailing': 1,
        },
    }


def build_rational_scale():
    """整数TSCAL=1/2 TZERO=1/4：换算为非整数有理值，须保留分子分母。"""
    cards = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'R'), card('TFORM1', '1J'),
        card('TSCAL1', 0.5), card('TZERO1', 0.25),
    ]
    row = struct.pack('>i', 3)  # 3/2 + 1/4 = 7/4
    return bintable(cards, row), {'num': 7, 'den': 4}


def build_commentary():
    """重复COMMENT/HISTORY与空白卡不丢失。"""
    cards = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
        card('COMMENT', 'first note'),
        card(''),
        card('COMMENT', 'second note'),
        card('HISTORY', 'version one'),
        card('HISTORY', 'version two'),
    ]
    return bintable(cards, struct.pack('>i', 1)), {
        'comments': ['first note', 'second note'],
        'history': ['version one', 'version two'],
    }


def build_multihdu():
    """主HDU + 不支持的IMAGE扩展 + BINTABLE：图像像素不解，区间保留。"""
    image_header = header_block([
        card('XTENSION', 'IMAGE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 3), card('PCOUNT', 0),
        card('GCOUNT', 1),
    ])
    image_data = pad_data(bytes(range(12)))
    table_header = header_block([
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
    ])
    table_data = pad_data(struct.pack('>i', 42))
    return (primary_header() + image_header + image_data + table_header
            + table_data), {'table_hdu': 2, 'image_hdu': 1, 'table_id': 42}


def build_bad_identity_cards():
    """若干结构性坏头：重复关键卡、BITPIX非法、NAXIS1与布局不符。"""
    dup = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
        card('NAXIS1', 5),
    ]
    bad_bitpix = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 12), card('NAXIS', 2),
        card('NAXIS1', 4), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
    ]
    # NAXIS1小于列宽合计(4)：真正的结构矛盾（行内只允许多余填充，不允许不足）
    bad_width = [
        card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2),
        card('NAXIS1', 3), card('NAXIS2', 1),
        card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 1),
        card('TTYPE1', 'ID'), card('TFORM1', '1J'),
    ]
    return {
        'duplicate_key': bintable(dup, struct.pack('>i', 1)),
        'bad_bitpix': bintable(bad_bitpix, struct.pack('>i', 1)),
        'bad_width': bintable(bad_width, struct.pack('>i', 1)),
    }
