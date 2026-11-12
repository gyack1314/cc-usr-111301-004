"""已有资料读取与版本字节核验。"""
import hashlib,json,tempfile,unittest
from pathlib import Path
from capture_data.io import load,asset,Source,Replay
from capture_data.domain import summarize

class MaterialTests(unittest.TestCase):
    def test_normal(self):
        v,d=load('data/normal.json');r=summarize(v,Path('data'))
        self.assertEqual(len(d),64);self.assertEqual(r['asset_count'],len(v['assets']));self.assertGreater(r['raw_bytes'],0)
    def test_identity_rejected(self):
        with self.assertRaises(ValueError):summarize(load('data/invalid-identity.json')[0],Path('data'))
    def test_boundary_sources(self):
        for p in Path('data').glob('*.json'):
            if p.stem=='invalid-identity':continue
            with self.subTest(name=p.name):self.assertIsInstance(summarize(load(p)[0],p.parent),dict)
    def test_read_at(self):
        v,_=load('data/normal.json');s=Source('data',v['assets'][0]);self.assertEqual(s.read_at(0,1),s.raw[:1]);self.assertEqual(s.calls,[{'start':0,'length':1}])
        for pair in ((-1,1),(0,len(s.raw)+1),(True,1)):
            with self.assertRaises(ValueError):s.read_at(*pair)
    def test_asset_path_and_json(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_bytes(b'abc');r={'path':'x','bytes':3,'sha256':hashlib.sha256(b'abc').hexdigest()};self.assertEqual(asset(d,r),b'abc')
            with self.assertRaises(ValueError):asset(d,{**r,'path':'../x'})
            for raw in ('[]','null','{"a":NaN}'):
                p.write_text(raw,encoding='utf8')
                with self.assertRaises(ValueError):load(p)
    def test_replay_position(self):
        r=Replay([{'request':{'method':'GET','url':'https://local.invalid/a'},'response':{'status':200}}])
        with self.assertRaises(ValueError):r.next('POST','https://local.invalid/a')
        self.assertEqual(r.position,0);self.assertEqual(r.next('GET','https://local.invalid/a')['status'],200)
        with self.assertRaises(EOFError):r.next('GET','https://local.invalid/a')
