"""FITS表采集测试：期望值取自独立构造的原值表（tests/fits_fixtures.py）。

原值表由手工bytes/struct拼出，与 capture_data.fits 不共享任何解析代码；
采集器输出只与该独立期望逐字段对照，并核对每个字节来源。
"""
import struct
import unittest

from capture_data.fits import collect, walk_hdus, FITSError
from tests import fits_fixtures as fx


def _read_desc_raw(data, start, length, kind):
    """测试侧独立解码描述符，不经过采集器。"""
    code = 'i' if kind == 'P' else 'q'
    return struct.unpack('>' + code * 2, data[start:start + length])


class HDUWalkTests(unittest.TestCase):
    def test_primary_and_bintable_spans(self):
        data, expect = fx.build_shared_heap()
        hdus, reader = walk_hdus(data)
        self.assertEqual([h.index for h in hdus], [0, 1])
        primary, table = hdus
        # 主HDU NAXIS=0：仅一个2880头块，无数据块
        self.assertEqual(primary.total_span, 2880)
        self.assertEqual(primary.data_size, 0)
        self.assertEqual(primary.data_padding, 0)
        # BINTABLE：头2880 + 数据52 + 尾部padding到2880
        self.assertEqual(table.header_start, 2880)
        self.assertEqual(table.data_start, 5760)
        self.assertEqual(table.data_size, 52)
        self.assertEqual(table.data_padding, 2880 - 52)
        self.assertEqual(table.total_span, 5760)
        # 跨度合计恰好覆盖文件，padding从未被当成下一个头
        self.assertEqual(primary.total_span + table.total_span, len(data))

    def test_unsupported_extension_kept_but_pixels_not_decoded(self):
        data, expect = fx.build_multihdu()
        hdus, _ = walk_hdus(data)
        self.assertEqual([h.kind for h in hdus],
                         ['primary', 'unsupported:IMAGE', 'bintable'])
        image = hdus[1]
        # 12字节像素 -> 数据块2880；原始区间保留
        self.assertEqual(image.data_size, 12)
        self.assertEqual(image.total_span, 2880 + 2880)
        # 采集仍能到达后面的BINTABLE（不靠猜长度）
        result = collect(data, expect['table_hdu'], ['ID'])
        self.assertEqual(
            result['records'][0]['fields']['ID']['elements'][0]['raw_integer'],
            expect['table_id'])

    def test_header_padding_after_end_is_measured(self):
        data, _ = fx.build_shared_heap()
        hdus, _ = walk_hdus(data)
        # 20张卡(含END) -> END后还有16张空格卡 = 1280字节头padding
        self.assertEqual(hdus[1].header_padding, 1280)

    def test_structural_bad_headers_rejected(self):
        for name, blob in fx.build_bad_identity_cards().items():
            with self.subTest(name=name):
                with self.assertRaises(FITSError):
                    # bad_width只有进入列布局才能识别，因此走完整采集
                    collect(blob, 1, ['ID'])

    def test_truncated_tail_rejected(self):
        data, _ = fx.build_shared_heap()
        with self.assertRaises(FITSError):
            walk_hdus(data[:-100])


class HeaderSemanticsTests(unittest.TestCase):
    def test_duplicate_comment_history_preserved(self):
        data, expect = fx.build_commentary()
        result = collect(data, 1, ['ID'])
        commentary = result['hdu']['commentary']
        comments = [c['text'] for c in commentary if c['keyword'] == 'COMMENT']
        history = [c['text'] for c in commentary if c['keyword'] == 'HISTORY']
        self.assertEqual(comments, expect['comments'])
        self.assertEqual(history, expect['history'])
        # 每条注释带来源且卡片序号不同
        sources = {c['source']['start'] for c in commentary}
        self.assertEqual(len(sources), len(commentary))

    def test_slash_inside_string_is_not_comment(self):
        data, expect = fx.build_slash_string()
        col_name = expect['column_name']  # 'orbit/label'
        result = collect(data, 1, ['ID', col_name])
        # 列名中的斜杠被完整保留，证明引号内斜杠不是注释
        label_col = next(c for c in result['layout']['columns']
                         if c['name'] == col_name)
        self.assertEqual(label_col['name'], col_name)
        for row, (exp_id, exp_label, exp_trail) in zip(
                result['records'], expect['expected']):
            label = row['fields'][col_name]['elements'][0]
            idv = row['fields']['ID']['elements'][0]['raw_integer']
            self.assertEqual(idv, exp_id)
            self.assertEqual(label['string'], exp_label)
            self.assertEqual(label['trailing_spaces'], exp_trail)
            # 原字节始终保留
            self.assertEqual(bytes.fromhex(label['raw_hex']),
                             label['raw_text'].encode('ascii'))

    def test_escaped_quote_in_string(self):
        # 独立小夹具：头卡里字面单引号按FITS写成两个单引号
        cards = [
            fx.card('XTENSION', 'BINTABLE'), fx.card('BITPIX', 8),
            fx.card('NAXIS', 2), fx.card('NAXIS1', 5), fx.card('NAXIS2', 1),
            fx.card('PCOUNT', 0), fx.card('GCOUNT', 1), fx.card('TFIELDS', 1),
            b"TTYPE1  = 'a''b'".ljust(80), fx.card('TFORM1', '5A'),
        ]
        blob = fx.bintable(cards, b"a'b  ")
        # 列名中的两个单引号应解析为一个字面引号
        result = collect(blob, 1, ["a'b"])
        self.assertEqual(
            result['records'][0]['fields']["a'b"]['elements'][0]['string'],
            "a'b")
        self.assertEqual(
            result['records'][0]['fields']["a'b"]['elements'][0][
                'trailing_spaces'], 2)


class ColumnTypeTests(unittest.TestCase):
    def test_all_fixed_types_tdim_nan_tnull_rational(self):
        data, expect = fx.build_type_garden()
        result = collect(data, 1,
                         ['FLAG', 'BYTE', 'SHORT', 'VEC', 'LONG',
                          'REAL', 'DOUB', 'TEXT'])
        fields = result['records'][0]['fields']
        e = expect['expected']

        flag = fields['FLAG']['elements'][0]
        self.assertTrue(flag['missing'])
        self.assertIsNone(flag['raw_logical'])

        self.assertEqual(fields['BYTE']['elements'][0]['raw_integer'],
                         e['byte_raw'])

        short = fields['SHORT']['elements'][0]
        # TNULL在原值域先判空，空值不缩放
        self.assertTrue(short['missing'])
        self.assertEqual(short['raw_integer'], -1)
        self.assertIsNone(short['scaled'])

        vec = fields['VEC']
        self.assertEqual([x['raw_integer'] for x in vec['elements']],
                         e['vec_raw'])
        self.assertEqual(vec['tdim'], e['tdim'])  # 第一维最快
        self.assertEqual(vec['cell_source']['length'], 8)

        longv = fields['LONG']['elements'][0]
        self.assertEqual(longv['raw_integer'], e['long_raw'])
        self.assertEqual(longv['scaled'], e['long_scaled'])

        real = fields['REAL']['elements'][0]
        # NaN是缺测标记，绝不输出"nan"字符串
        self.assertTrue(real['missing'])
        self.assertIsNone(real['raw_float'])
        self.assertIsNone(real['scaled'])

        doub = fields['DOUB']['elements'][0]
        self.assertEqual(doub['raw_float'], e['doub_raw'])
        self.assertEqual(doub['scaled'], e['doub_scaled'])

        text = fields['TEXT']['elements'][0]
        self.assertEqual(text['string'], e['text'])
        self.assertEqual(text['trailing_spaces'], e['text_trailing'])
        self.assertEqual(text['raw_text'], 'abc ')

    def test_rational_integer_scaling_keeps_numerator_denominator(self):
        data, expect = fx.build_rational_scale()
        result = collect(data, 1, ['R'])
        elem = result['records'][0]['fields']['R']['elements'][0]
        self.assertEqual(elem['raw_integer'], 3)
        self.assertEqual(elem['scaled_kind'], 'rational')
        self.assertEqual(elem['scaled']['numerator'], expect['num'])
        self.assertEqual(elem['scaled']['denominator'], expect['den'])

    def test_logical_illegal_byte_rejected(self):
        cards = [
            fx.card('XTENSION', 'BINTABLE'), fx.card('BITPIX', 8),
            fx.card('NAXIS', 2), fx.card('NAXIS1', 1), fx.card('NAXIS2', 1),
            fx.card('PCOUNT', 0), fx.card('GCOUNT', 1), fx.card('TFIELDS', 1),
            fx.card('TTYPE1', 'L'), fx.card('TFORM1', '1L'),
        ]
        blob = fx.bintable(cards, b'X')
        with self.assertRaises(FITSError):
            collect(blob, 1, ['L'])


class DescriptorHeapTests(unittest.TestCase):
    def test_two_rows_share_heap_slice_not_flagged_duplicate(self):
        data, expect = fx.build_shared_heap()
        result = collect(data, 1, ['ID', 'SPECTRUM', 'QUALITY'],
                         max_array_elements=1000)
        self.assertTrue(result['structure_valid'])
        self.assertEqual(result['issues'], [])
        s0 = result['records'][0]['fields']['SPECTRUM']
        s1 = result['records'][1]['fields']['SPECTRUM']
        # 两行的堆切片证据指向完全相同的文件区间
        self.assertEqual(s0['heap_slice']['start'], s1['heap_slice']['start'])
        self.assertEqual(s0['heap_slice']['length'],
                         s1['heap_slice']['length'])
        # 元素值来自独立原值表
        for field in (s0, s1):
            vals = [e['raw_float'] for e in field['elements']]
            self.assertEqual(vals, expect['heap_values'])
        # 描述符自身的字节来源不同（在各自行内）
        self.assertNotEqual(
            s0['descriptor']['descriptor_source']['start'],
            s1['descriptor']['descriptor_source']['start'])

    def test_q_descriptor_offset_out_of_range_no_cross_hdu_read(self):
        data, expect = fx.build_q_out_of_range()
        result = collect(data, 1, ['ID', 'SPECTRUM', 'QUALITY'],
                         max_array_elements=1000)
        self.assertFalse(result['structure_valid'])
        issue = result['issues'][0]
        self.assertEqual(issue['kind'], 'heap_out_of_range')
        self.assertEqual(issue['row'], expect['bad_row'])
        self.assertEqual(issue['offset'], expect['bad_offset'])
        bad = result['records'][1]['fields']['SPECTRUM']
        self.assertFalse(bad['readable'])
        self.assertIsNone(bad['heap_slice'])
        self.assertEqual(bad['error'], 'heap_out_of_range')
        # 越界窗口的绝对位置：数据区起点 + THEAP + 10000
        layout = result['layout']
        bad_abs_start = (result['hdu']['data_start'] + layout['theap']
                         + expect['bad_offset'])
        # 任何读取调用都不得与越界窗口发生区间重叠
        bad_window_end = bad_abs_start + expect['bad_count'] * 8
        for call in result['read_calls']:
            call_end = call['start'] + call['length']
            self.assertFalse(
                call['start'] < bad_window_end and call_end > bad_abs_start,
                f'读取进入越界堆窗口: {call}')
        # 所有读取都不越过本HDU原始区间（绝不跨HDU）
        hdu_end = (result['hdu']['raw_source']['start']
                   + result['hdu']['raw_source']['length'])
        for call in result['read_calls']:
            self.assertLessEqual(call['start'] + call['length'], hdu_end)
        # 第一行合法描述符仍可读
        good = result['records'][0]['fields']['SPECTRUM']
        self.assertTrue(good['readable'])
        self.assertEqual(good['descriptor']['descriptor_kind'], 'Q')

    def test_non_default_theap_gap_belongs_to_pcount(self):
        data, expect = fx.build_shared_heap()
        result = collect(data, 1, ['SPECTRUM'])
        layout = result['layout']
        self.assertTrue(layout['theap_explicit'])
        self.assertEqual(layout['theap'], expect['theap'])
        self.assertEqual(layout['theap_gap_size'], 8)
        self.assertEqual(layout['heap_size'], 16)
        # 间隙+堆 == PCOUNT
        self.assertEqual(layout['theap_gap_size'] + layout['heap_size'],
                         layout['pcount'])
        # 堆切片绝对位置 = 数据区起点 + THEAP + offset
        sl = result['records'][0]['fields']['SPECTRUM']['heap_slice']
        self.assertEqual(sl['start'],
                         result['hdu']['data_start'] + expect['theap'])

    def test_default_theap_follows_row_area(self):
        # 无THEAP卡：默认THEAP=NAXIS1*NAXIS2，间隙为0
        cards = [
            fx.card('XTENSION', 'BINTABLE'), fx.card('BITPIX', 8),
            fx.card('NAXIS', 2), fx.card('NAXIS1', 8), fx.card('NAXIS2', 1),
            fx.card('PCOUNT', 16), fx.card('GCOUNT', 1), fx.card('TFIELDS', 1),
            fx.card('TTYPE1', 'S'), fx.card('TFORM1', '1PD'),
        ]
        rows = struct.pack('>ii', 2, 0)
        heap = struct.pack('>dd', 9.0, 10.0)
        blob = fx.bintable(cards, rows + heap)
        result = collect(blob, 1, ['S'])
        self.assertFalse(result['layout']['theap_explicit'])
        self.assertEqual(result['layout']['theap'], 8)
        self.assertEqual(result['layout']['theap_gap_size'], 0)
        self.assertTrue(result['structure_valid'])

    def test_max_array_elements_limit_flagged(self):
        data, _ = fx.build_shared_heap()
        result = collect(data, 1, ['SPECTRUM'], max_array_elements=1)
        kinds = {i['kind'] for i in result['issues']}
        self.assertIn('array_too_large', kinds)


class RowPaddingTests(unittest.TestCase):
    def test_row_end_padding_and_data_trailing_padding(self):
        data, expect = fx.build_row_padding()
        result = collect(data, 1, ['ID', 'FLAG'])
        self.assertEqual(result['layout']['naxis1'], expect['naxis1'])
        for row, (exp_id, exp_flag) in zip(result['records'],
                                           expect['expected']):
            self.assertEqual(
                row['fields']['ID']['elements'][0]['raw_integer'], exp_id)
            self.assertEqual(
                row['fields']['FLAG']['elements'][0]['raw_logical'],
                exp_flag)
            # 行来源宽度是含padding的完整NAXIS1
            self.assertEqual(row['row_source']['length'], expect['naxis1'])
        # FLAG单元来源在padding之前，绝不把0xDE当逻辑字节
        flag_src = result['records'][0]['fields']['FLAG']['elements'][0]['source']
        self.assertEqual(
            flag_src['start'],
            result['hdu']['data_start'] + 4)
        # 数据区尾部padding独立记录
        self.assertEqual(result['hdu']['data_trailing_padding'],
                         2880 - (expect['naxis1'] * expect['naxis2']
                                 + expect['pcount']))


class SelectionAndProvenanceTests(unittest.TestCase):
    def test_column_selection_frozen_unselected_heap_not_read(self):
        data, _ = fx.build_shared_heap()
        result = collect(data, 1, ['ID'])
        # 未选堆列：描述符结构仍核验，但堆负载读取不存在
        heap_payload_calls = [
            c for c in result['read_calls']
            if c['start'] >= result['hdu']['data_start'] + result['layout']['theap']
        ]
        self.assertEqual(heap_payload_calls, [])
        # 记录里只有选中的列
        self.assertEqual(
            set(result['records'][0]['fields']), {'ID'})
        # 布局里SPECTRUM标记为未选
        spec = next(c for c in result['layout']['columns']
                    if c['name'] == 'SPECTRUM')
        self.assertFalse(spec['selected'])

    def test_unknown_column_and_bad_hdu_rejected(self):
        data, _ = fx.build_shared_heap()
        with self.assertRaises(FITSError):
            collect(data, 1, ['NOPE'])
        with self.assertRaises(FITSError):
            collect(data, 9, ['ID'])
        with self.assertRaises(FITSError):
            collect(data, 0, ['ID'])  # 主HDU不是BINTABLE

    def test_every_field_source_matches_raw_bytes(self):
        data, expect = fx.build_shared_heap()
        result = collect(data, 1, ['ID', 'SPECTRUM', 'QUALITY'])
        data_start = result['hdu']['data_start']
        for row in result['records']:
            for name, field in row['fields'].items():
                if 'cell_source' in field:
                    src = field['cell_source']
                    self.assertEqual(
                        data[src['start']:src['start'] + src['length']].hex(),
                        field['raw_hex'])
                    self.assertEqual(src['start'],
                                     data_start + src['hdu_data_relative'])
                else:
                    d = field['descriptor']
                    ds = d['descriptor_source']
                    # 描述符字节必须落在该字段的行内单元区间
                    self.assertEqual(
                        ds['hdu_data_relative'],
                        row['row_source']['hdu_data_relative']
                        + result['layout']['columns'][1]['row_offset'])
                    # 来源切片与采集器给出的描述符原值一致（独立重读源文件）
                    count, offset = _read_desc_raw(
                        data, ds['start'], ds['length'],
                        d['descriptor_kind'])
                    self.assertEqual(count, d['count'])
                    self.assertEqual(offset, d['offset'])
                    if field['readable']:
                        hs = field['heap_slice']
                        self.assertEqual(
                            data[hs['start']:hs['start'] + hs['length']].hex(),
                            hs['raw_hex'])

    def test_raw_hdu_intervals_cover_whole_file(self):
        data, _ = fx.build_shared_heap()
        result = collect(data, 1, ['ID'])
        intervals = [(h['source']['start'], h['source']['length'])
                     for h in result['raw_hdus']]
        self.assertEqual(intervals[0][0], 0)
        total = sum(length for _, length in intervals)
        self.assertEqual(total, len(data))


class RealFileOracleTests(unittest.TestCase):
    """用独立硬编码期望核对三个真实资产（期望不来自采集器）。"""

    def test_catalog_shared_p_descriptor(self):
        with open('data/catalog.fits', 'rb') as f:
            data = f.read()
        result = collect(data, 1, ['ID', 'SPECTRUM', 'QUALITY'],
                         max_array_elements=1000)
        self.assertTrue(result['structure_valid'])
        # 独立期望（手工struct解码得到，见测试说明）
        self.assertEqual(
            [r['fields']['ID']['elements'][0]['raw_integer']
             for r in result['records']], [1, 2])
        qualities = [r['fields']['QUALITY']['elements'][0]
                     for r in result['records']]
        self.assertTrue(qualities[0]['missing'])      # -999 判空
        self.assertEqual(qualities[1]['scaled'], 7)   # 2*3+1
        # 两行共享同一堆切片
        s0 = result['records'][0]['fields']['SPECTRUM']['heap_slice']
        s1 = result['records'][1]['fields']['SPECTRUM']['heap_slice']
        self.assertEqual(s0, s1)
        vals = [e['raw_float']
                for e in result['records'][0]['fields']['SPECTRUM']['elements']]
        self.assertEqual(vals, [1.5, 2.5])

    def test_wide_q_descriptor_and_explicit_theap(self):
        with open('data/wide-descriptor.fits', 'rb') as f:
            data = f.read()
        result = collect(data, 1, ['ID', 'SPECTRUM'])
        self.assertEqual(result['layout']['theap'], 52)
        self.assertEqual(
            result['records'][0]['fields']['SPECTRUM']
            ['descriptor']['descriptor_kind'], 'Q')

    def test_outside_heap_flagged(self):
        with open('data/outside-heap.fits', 'rb') as f:
            data = f.read()
        result = collect(data, 1, ['ID', 'SPECTRUM', 'QUALITY'])
        self.assertFalse(result['structure_valid'])
        self.assertEqual(result['issues'][0]['kind'], 'heap_out_of_range')


if __name__ == '__main__':
    unittest.main()