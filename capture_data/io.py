"""读取冻结清单及验证资产字节身份。"""
import hashlib,json
from pathlib import Path

def load(path):
    raw=Path(path).read_bytes();v=json.loads(raw.decode('utf8'),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('非有限JSON常量')))
    if not isinstance(v,dict):raise ValueError('清单必须为对象')
    return v,hashlib.sha256(raw).hexdigest()

def asset(root,ref):
    base=Path(root).resolve();p=(base/ref['path']).resolve()
    if not p.is_relative_to(base):raise ValueError('资产路径越界')
    raw=p.read_bytes()
    if type(ref['bytes']) is not int or ref['bytes']<0:raise ValueError('资产长度非法')
    if len(raw)!=ref['bytes'] or hashlib.sha256(raw).hexdigest()!=ref['sha256']:raise ValueError('资产身份不符')
    return raw

class Source:
    """固定版本字节范围读取器，记录实际请求。"""
    def __init__(self,root,ref):self.raw=asset(root,ref);self.sha256=ref['sha256'];self.calls=[]
    def read_at(self,start,length):
        if type(start) is not int or type(length) is not int or start<0 or length<0 or start+length>len(self.raw):raise ValueError('范围越界')
        self.calls.append({'start':start,'length':length});return self.raw[start:start+length]

class Replay:
    """按序返回冻结HTTP记录，不完成采集策略。"""
    def __init__(self,records):self.records=records;self.position=0;self.requests=[]
    def next(self,method,url,headers=None):
        if self.position>=len(self.records):raise EOFError('没有后续响应')
        row=self.records[self.position];expected=row['request']
        if (method,url)!=(expected['method'],expected['url']):raise ValueError('请求不匹配')
        self.requests.append({'method':method,'url':url,'headers':dict(headers or {})});self.position+=1
        return row['response']
