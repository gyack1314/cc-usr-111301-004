# FITS观测目录

文件首HDU是NAXIS=0主头，下一HDU为BINTABLE。卡片80字节，头和数据区各自向2880取整。BITPIX/NAXIS/PCOUNT/GCOUNT的跨度遵循FITS，BINTABLE数据跨度=NAXIS1*NAXIS2+PCOUNT。只接收ASCII头；单引号字符串内两个单引号表示字面引号，只有引号外的斜杠才是注释起点。

原表含ID J、SPECTRUM PD/QD、QUALITY I。变长描述符的两数是元素个数和相对THEAP的字节offset，P为两个i32，Q为两个i64，皆非负。PCOUNT含固定行区到THEAP间隙和堆；堆可重叠共享。TFORM的P/Q限定元素原子类型，本项目不支持变长字符串、bit列X或复杂数C/M；数组最大长度括号只作上界不能替代实际count。TDIM如存在，其维度积等固定重复数，第一维变化最快。

TNULL在整数原域判断，之后TSCAL*x+TZERO；整数换算要求精确有理值，浮点另保留原位型并给有穷数或缺测标记。L只接受T/F/空缺，A保留原切片，表行没有跨HDU读取权限。材料包括合法共享堆、Q描述符和结构上合法但实际越界的描述符。初始工具不解HDU。

## 采集器（capture_data.fits）

`walk_hdus`按80字节卡片与2880字节块顺序遍历，下一个HDU起点完全由本HDU跨度（头块 + 向2880取整的数据块）算术决定，绝不扫描字节找头；END后的头padding必须全空格，数据区尾部padding只记录不解析。主HDU（NAXIS=0）占且仅占一个头块；支持BINTABLE，图像/其他扩展保留原始字节区间但不解像素，跨度无法确定即报错，不靠猜长度继续。关键卡（SIMPLE/XTENSION/BITPIX/NAXIS/NAXISn/PCOUNT/GCOUNT/TFIELDS/TTYPEn/TFORMn等）缺失、重复或必需值不符均抛`FITSError`。字符串卡内两个单引号转义一个字面引号，只有引号外斜杠是注释；重复COMMENT/HISTORY逐条保留。

`collect(data, hdu_index, columns, max_array_elements=None)`在入口冻结列选择（元组、查重），返回逐字段记录、描述符与堆切片证据、布局证据、问题清单和全部读取调用。列类型覆盖L/B/I/J/K/E/D/A及P/Q变长描述符，数值一律大端。TFORM重复数与列偏移累加结果必须不超过NAXIS1，差值记为每行尾padding；TDIM仅对定长数组生效，维度积等于重复数，按FITS第一维最快次序保留。THEAP缺省等于行区大小，显式间隙计入PCOUNT；描述符count/offset为负或窗口越界记`heap_out_of_range`问题且不读负载，多个描述符共享同一堆片段合法且不重复报错。整数TNULL先在原值域判空再做有理换算（非整数给分子/分母），浮点NaN作为缺测（JSON输出allow_nan=False），A列去右空格但保留原字节。未选堆列不读取其负载，但其描述符仍参与结构检查，因此存在结构问题时`structure_valid`为false，不能宣称整文件数据有效。
