"""额外对抗边界：负值描述符、部分共享、非零padding、选择冻结一致性。

期望仍由独立夹具与手工struct给出，不引用采集器输出作为期望。
"""
import struct
import unittest

from capture_data.fits import collect, walk_hdus, FITSError
from tests import fits_fixtures as fx


class NegativeDescriptorTests(unittest.TestCase):
    def _build(self, count, offset):
        naxis1, naxis2, pcount, theap = 14, 1, 16, 14
        cards = fx.std_table_cards(naxis1, naxis2, pcount, theap)
        # 只留1行：把QUALITY的TNULL等卡保留无妨
        row = (struct.pack('>i', 1) + struct.pack('>ii', count, offset)
               + struct.pack('>h', 0))
        heap = struct.pack('>dd', 1.0, 2.0)
        return fx.bintable(cards, row + heap)

    def test_negative_count_flagged_not_read(self):
        data = self._build(-1, 0)
        result = collect(data, 1, ['SPECTRUM'])
        self.assertFalse(result['structure_valid'])
        kind = result['issues'][0]['kind']
        self.assertEqual(kind, 'heap_out_of_range')
        field = result['records'][0]['fields']['SPECTRUM']
        self.assertTrue(field['descriptor']['negative_count'])
        self.assertFalse(field['readable'])
        self.assertIsNone(field['heap_slice'])

    def test_negative_offset_flagged_not_read(self):
        data = self._build(2, -8)
        result = collect(data, 1, ['SPECTRUM'])
        self.assertFalse(result['structure_valid'])
        field = result['records'][0]['fields']['SPECTRUM']
        self.assertTrue(field['descriptor']['negative_offset'])
        self.assertFalse(field['readable'])

    def test_zero_count_is_valid_empty_array(self):
        data = self._build(0, 0)
        result = collect(data, 1, ['SPECTRUM'])
        self.assertTrue(result['structure_valid'], result['issues'])
        field = result['records'][0]['fields']['SPECTRUM']
        self.assertTrue(field['readable'])
        self.assertEqual(field['elements'], [])
        self.assertEqual(field['heap_slice']['length'], 0)


class PartialShareTests(unittest.TestCase):
    def test_partial_overlap_is_legal(self):
        # 三个描述符：[0,16) [8,16) [0,8)，两两部分重叠，全部合法
        naxis1, naxis2, pcount, theap = 14, 3, 24, 42
        cards = fx.std_table_cards(naxis1, naxis2, pcount, theap)
        rows = b''
        rows += struct.pack('>i', 1) + struct.pack('>ii', 2, 0) + struct.pack('>h', 0)
        rows += struct.pack('>i', 2) + struct.pack('>ii', 1, 8) + struct.pack('>h', 0)
        rows += struct.pack('>i', 3) + struct.pack('>ii', 1, 0) + struct.pack('>h', 0)
        gap = b'\x11' * (theap - naxis1 * naxis2)
        heap = struct.pack('>dd', 4.0, 5.0)
        heap += b'\x00' * (pcount - len(gap) - len(heap))
        data = fx.bintable(cards, rows + gap + heap)
        result = collect(data, 1, ['SPECTRUM'])
        self.assertTrue(result['structure_valid'], result['issues'])
        vals = [
            [e['raw_float'] for e in r['fields']['SPECTRUM']['elements']]
            for r in result['records']
        ]
        # 独立期望：直接按堆字节struct解释
        self.assertEqual(vals, [[4.0, 5.0], [5.0], [4.0]])


class PaddingRobustnessTests(unittest.TestCase):
    def test_nonzero_data_trailing_padding_never_parsed_as_header(self):
        # 构造尾部padding非空格/非零的文件，遍历仍精确按跨度结束
        data, _ = fx.build_shared_heap()
        # 数据区尾部padding从 data_start+52 到下一块边界
        # 写入伪XTENSION字样，验证不会被当成新HDU头
        poisoned = bytearray(data)
        pad_pos = 5760 + 52
        poisoned[pad_pos:pad_pos + 8] = b'XTENSION'
        hdus, _ = walk_hdus(bytes(poisoned))
        self.assertEqual(len(hdus), 2)
        self.assertEqual(hdus[1].data_padding, 2880 - 52)

    def test_nonzero_header_padding_rejected_if_not_spaces(self):
        # END之后的头padding若不是全空格，必须报错而不是容忍
        data, _ = fx.build_shared_heap()
        poisoned = bytearray(data)
        # END在卡片19（索引），其后第一张padding卡的开头
        end_card_end = 2880 + 20 * 80
        poisoned[end_card_end] = ord('X')
        with self.assertRaises(FITSError):
            walk_hdus(bytes(poisoned))


class SelectionFreezeTests(unittest.TestCase):
    def test_selection_tuple_independent_of_caller_mutation(self):
        data, _ = fx.build_shared_heap()
        chosen = ['ID', 'SPECTRUM']
        result = collect(data, 1, chosen)
        self.assertEqual(result['selected_columns'], ['ID', 'SPECTRUM'])
        # 调用后修改原列表不影响已冻结选择
        chosen.append('QUALITY')
        self.assertEqual(result['selected_columns'], ['ID', 'SPECTRUM'])

    def test_duplicate_requested_columns_rejected(self):
        data, _ = fx.build_shared_heap()
        with self.assertRaises(FITSError):
            collect(data, 1, ['ID', 'ID'])

    def test_unselected_descriptor_still_structurally_checked(self):
        # 越界文件中不选SPECTRUM：结构问题仍须上报，不能宣称整文件有效
        with open('data/outside-heap.fits', 'rb') as f:
            data = f.read()
        result = collect(data, 1, ['ID', 'QUALITY'])
        self.assertFalse(result['structure_valid'])
        self.assertEqual(result['issues'][0]['kind'], 'heap_out_of_range')
        # 但堆负载确实没有被读取
        heap_calls = [
            c for c in result['read_calls']
            if c['start'] >= result['hdu']['data_start']
            + result['layout']['theap']
        ]
        self.assertEqual(heap_calls, [])

    def test_theap_beyond_data_area_rejected(self):
        # THEAP指向行区+PCOUNT之后：结构性硬错误
        naxis1, naxis2, pcount = 8, 1, 8
        cards = [
            fx.card('XTENSION', 'BINTABLE'), fx.card('BITPIX', 8),
            fx.card('NAXIS', 2), fx.card('NAXIS1', naxis1),
            fx.card('NAXIS2', naxis2), fx.card('PCOUNT', pcount),
            fx.card('GCOUNT', 1), fx.card('TFIELDS', 1),
            fx.card('TTYPE1', 'S'), fx.card('TFORM1', '1PD'),
            fx.card('THEAP', naxis1 * naxis2 + pcount + 4),
        ]
        blob = fx.bintable(cards, struct.pack('>ii', 1, 0) + b'\x00' * 8)
        with self.assertRaises(FITSError):
            collect(blob, 1, ['S'])


if __name__ == '__main__':
    unittest.main()
