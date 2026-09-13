import csv, hashlib, json, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class NormalizeRankingTests(unittest.TestCase):
 def invoke(self, text):
  root=Path(tempfile.mkdtemp(prefix='ranking-gate-'))
  manifest=root/'manifest.csv';manifest.write_text('compound_id,canonical_smiles\nB,CC\nA,CCC\n',encoding='utf-8')
  meta=dict(pdb_id='TEST',consensus_id='C001',tier='T2',outputs=dict(molecule_manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest()))
  (root/'input.json').write_text(json.dumps(meta));(root/'ranked.txt').write_text(text)
  proc=subprocess.run([sys.executable,str(ROOT/'scripts/12_normalize_drugclip_results.py'),'--ranked',str(root/'ranked.txt'),'--manifest',str(manifest),'--output',str(root/'out.csv'),'--pdb-id','TEST','--consensus-id','C001','--tier','T2','--input-run',str(root/'input.json')],capture_output=True,text=True)
  return proc,root
 def test_ties_use_stable_id(self):
  proc,root=self.invoke('CC\t1\nCCC\t1\n');self.assertEqual(proc.returncode,0,proc.stderr)
  with (root/'out.csv').open(encoding='utf-8-sig') as f:self.assertEqual([r['compound_id'] for r in csv.DictReader(f)],['A','B'])
 def test_missing_duplicate_nan_rejected(self):
  for text in ('CC\t1\n','CC\t1\nCC\t1\n','CC\tnan\nCCC\t1\n'):
   proc,root=self.invoke(text);self.assertNotEqual(proc.returncode,0);self.assertFalse((root/'out.csv').exists())
if __name__=='__main__':unittest.main()
