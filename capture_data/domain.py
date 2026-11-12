"""清点空间观测FITS表的堆描述符与列值采集原始资产，不解释内容格式。"""
from .io import asset

def summarize(v,root):
    if v['kind']!='fits-table':raise ValueError('资料类别不匹配')
    rows=v['assets']
    if not isinstance(rows,list) or not rows:raise ValueError('资产清单为空')
    paths=[x['path'] for x in rows]
    if len(paths)!=len(set(paths)):raise ValueError('重复资产路径')
    total=sum(len(asset(root,r)) for r in rows)
    return {'kind':v['kind'],'asset_count':len(rows),'raw_bytes':total,'http_records':len(v.get('calls',[]))}
