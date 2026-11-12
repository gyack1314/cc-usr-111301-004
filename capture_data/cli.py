"""输出原始清单、资产数量与FITS表采集结果。"""
import argparse
import json
from pathlib import Path

from .io import load
from .domain import summarize
from .fits import collect


def main():
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['inspect', 'collect'])
    p.add_argument('file')
    p.add_argument('--hdu', type=int, default=None,
                   help='目标HDU序号；缺省取清单hdu_index')
    p.add_argument('--columns', default=None,
                   help='逗号分隔的请求列；缺省取清单columns')
    p.add_argument('--max-array-elements', type=int, default=None,
                   dest='max_array_elements')
    a = p.parse_args()

    v, digest = load(a.file)
    root = Path(a.file).parent
    if a.command == 'inspect':
        print(json.dumps(
            {'manifest_sha256': digest,
             'summary': summarize(v, root)},
            ensure_ascii=False, indent=2))
        return

    # 采集入口：选择在开始时冻结（collect内部再冻结一次并查重）
    ref = v.get('source') or v['assets'][0]
    from .io import asset
    raw = asset(root, ref)
    hdu_index = a.hdu if a.hdu is not None else v['hdu_index']
    if a.columns is not None:
        columns = tuple(c.strip() for c in a.columns.split(',') if c.strip())
    else:
        columns = tuple(v['columns'])
    limit = a.max_array_elements
    if limit is None:
        limit = v.get('max_array_elements')
    result = collect(bytes(raw), hdu_index, columns,
                     max_array_elements=limit)
    print(json.dumps(
        {'manifest_sha256': digest, 'collection': result},
        ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
