# 空间观测FITS表的堆描述符与列值采集

此项目保存科学资料交换的原始文件、清单及冻结HTTP记录。`capture_data.io`读取JSON、核对资产身份，提供固定字节范围读取和顺序响应重放；`capture_data.domain`只核对资料类别、路径和数量，`capture_data.cli`输出概况。[资料约定](docs/data-contract.md)说明时间、字段与格式范围。所有站点、观测和机构资料均为合成，不连接航天器、接收机或真实研究系统。

Python3.11以上，无第三方Python依赖。项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q capture_data tests
python -m capture_data.cli inspect data/normal.json
```

`inspect`不解文件格式或选择HTTP版本，只显示资产身份及数量；`invalid-identity.json`有预期资产摘要错误，其余清单有合法来源身份但内容可含需领域操作识别的异常。`Source.read_at`记录请求范围，`Replay.next`记录实际method、URL与headers；它不替调用者完成协议或缓存策略。

Dockerfile固定Linux基础镜像，`docker build -t science-material .`构建；`docker run --rm science-material`读取概况，`docker run --rm science-material python -m unittest discover -s tests -v`运行已有检查。资产均可离线读取，不需要下载数据或真实账号。
