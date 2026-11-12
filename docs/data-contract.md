# FITS观测目录

文件首HDU是NAXIS=0主头，下一HDU为BINTABLE。卡片80字节，头和数据区各自向2880取整。BITPIX/NAXIS/PCOUNT/GCOUNT的跨度遵循FITS，BINTABLE数据跨度=NAXIS1*NAXIS2+PCOUNT。只接收ASCII头；单引号字符串内两个单引号表示字面引号，只有引号外的斜杠才是注释起点。

原表含ID J、SPECTRUM PD/QD、QUALITY I。变长描述符的两数是元素个数和相对THEAP的字节offset，P为两个i32，Q为两个i64，皆非负。PCOUNT含固定行区到THEAP间隙和堆；堆可重叠共享。TFORM的P/Q限定元素原子类型，本项目不支持变长字符串、bit列X或复杂数C/M；数组最大长度括号只作上界不能替代实际count。TDIM如存在，其维度积等固定重复数，第一维变化最快。

TNULL在整数原域判断，之后TSCAL*x+TZERO；整数换算要求精确有理值，浮点另保留原位型并给有穷数或缺测标记。L只接受T/F/空缺，A保留原切片，表行没有跨HDU读取权限。材料包括合法共享堆、Q描述符和结构上合法但实际越界的描述符。初始工具不解HDU。
