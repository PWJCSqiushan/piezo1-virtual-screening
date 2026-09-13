import csv, io, json, sys, tempfile, unittest, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
try:
 from piezo_vs.candidate_audit import audit_zip, digest, validate_ranking
 HAS_RDKIT=True
except ImportError:
 HAS_RDKIT=False


def csv_bytes(rows):
 f=io.StringIO(newline='');w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows);return f.getvalue().encode()

@unittest.skipUnless(HAS_RDKIT,'RDKit required; use the DrugCLIP environment')
class CandidateAuditTests(unittest.TestCase):
 def fixture(self, kind='valid'):
  root=Path(tempfile.mkdtemp(prefix='candidate-audit-'))
  rows=[]; files={}
  for task,smiles in [('TEST_C001','CCO'),('TEST_C002','CCO' if kind!='conflict' else 'CCN')]:
   row=dict(rank='1',compound_id='A',entity_id='E1',smiles=smiles,drugclip_score='1',vina_docking_score='',pdb_id='TEST',consensus_id=task.split('_')[1],job_id=task,library_name='finite')
   if kind=='nan':row['drugclip_score']='nan'
   if kind=='rank':row['rank']='2'
   raw=dict(entity_id=row['entity_id'],mol_id=row['compound_id'],library_name='finite',smiles=smiles,drugclip_score=row['drugclip_score'],vina_docking_score='')
   meta=dict(pdb_id='TEST',consensus_id=row['consensus_id'],server=dict(job_id=task,original_export_file='original.csv'))
   prefix='tasks/'+task+'/'
   files.update({prefix+'ranked_compounds.csv':csv_bytes([row]),prefix+'molecule_manifest.csv':csv_bytes([row]),prefix+'input_run.json':json.dumps(meta).encode(),prefix+'original.csv':csv_bytes([raw])})
  manifest='\n'.join(digest(b)+'  '+n for n,b in files.items())
  if kind=='hash': files[next(iter(files))]+=b'\n'
  files['SHA256SUMS.txt']=manifest.encode()
  archive=root/'handoff.zip'
  with zipfile.ZipFile(archive,'w') as z:
   for n,b in files.items(): z.writestr(n,b)
  return root,archive
 def test_traceable_dedup_and_eligibility(self):
  root,archive=self.fixture();s=audit_zip(archive,root/'out')
  self.assertEqual((s['record_count'],s['source_id_count'],s['canonical_structure_count']),(2,1,1))
  self.assertEqual(s['exploratory_pool_status'],'eligible')
  self.assertEqual(s['platform_pipeline_reproduction'],'not_reproduced')
  with (root/'out/source_mapping.csv').open(encoding='utf-8-sig') as f:
   self.assertEqual(len(list(csv.DictReader(f))),2)
 def test_hash_conflict_nan_rank_block_candidate_export(self):
  for kind in ('hash','conflict','nan','rank'):
   with self.subTest(kind=kind):
    root,archive=self.fixture(kind);s=audit_zip(archive,root/'out')
    self.assertEqual(s['exploratory_pool_status'],'blocked')
    self.assertFalse((root/'out/candidate_master.csv').exists())
    self.assertTrue((root/'out/source_records.csv').exists())
 def test_no_overwrite(self):
  root,archive=self.fixture();audit_zip(archive,root/'out')
  with self.assertRaises(FileExistsError):audit_zip(archive,root/'out')
 def test_duplicate_and_reverse_score_rejected(self):
  for rows in ([dict(rank=1,compound_id='A',drugclip_score=1),dict(rank=2,compound_id='A',drugclip_score=0)], [dict(rank=1,compound_id='A',drugclip_score=0),dict(rank=2,compound_id='B',drugclip_score=1)]):
   with self.assertRaises(ValueError):validate_ranking(rows)
if __name__=='__main__':unittest.main()
