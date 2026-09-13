import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from piezo_vs.docking import read_candidate_list
from piezo_vs.candidate_selection import write_csv_rows
from piezo_vs.docking_validation import inspect_poses
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT=True
except ImportError:
    HAS_RDKIT=False


class FormalGateTests(unittest.TestCase):
    def test_unchecked_batch_cannot_start(self):
        from piezo_vs.docking_validation import check_technical_evidence
        with self.assertRaisesRegex(ValueError, '3 successful'):
            check_technical_evidence([],dict(purpose='formal_screening'),'receptor','binary',{'A'})

    def test_pending_review_or_assessment_blocks_input(self):
        directory=Path(tempfile.mkdtemp(prefix='formal-gate-test-'))
        row=dict(compound_id='A',canonical_smiles='CC',site_direction='C006_C007',pdb_id='8YEZ',consensus_id='C006',
                 approved='true',evaluation_complete='true',review_status='approved')
        path=directory/'candidates.csv'
        for changed in (dict(approved='false'),dict(evaluation_complete='false'),dict(review_status='pending')):
            write_csv_rows(path,[dict(row,**changed)])
            with self.assertRaises(ValueError):read_candidate_list(path,'8YEZ','C006',purpose='formal_screening')


@unittest.skipUnless(HAS_RDKIT,'RDKit required for real SDF validation')
class PoseTests(unittest.TestCase):
    def fixture(self,mode='valid'):
        directory=Path(tempfile.mkdtemp(prefix='gnina-pose-test-'))
        candidates=[dict(compound_id=str(i),canonical_smiles=s,site_direction='C006_C007',pdb_id='8YEZ',consensus_id='C006',preparation_status='prepared')
                    for i,s in enumerate(('CCO','CCN','CCC'))]
        path=directory/'docked.sdf'; writer=Chem.SDWriter(str(path))
        for row in (candidates[:2] if mode=='missing' else candidates):
            mol=Chem.AddHs(Chem.MolFromSmiles(row['canonical_smiles']));AllChem.EmbedMolecule(mol,randomSeed=1)
            mol.SetProp('_Name',row['compound_id'])
            for k,v in row.items():mol.SetProp(k,v)
            for k,v in [('minimizedAffinity','-5'),('CNNscore','0.5'),('CNNaffinity','4')]:mol.SetProp(k,v)
            if mode=='nan':mol.SetProp('CNNscore','nan')
            if mode=='wrong_site':mol.SetProp('consensus_id','C016')
            writer.write(mol)
        writer.close()
        if mode=='corrupt':
            with path.open('a') as handle:handle.write('bad record\n$$$$\n')
        return path,candidates

    def test_valid_output_reconciles_all_three(self):
        path,candidates=self.fixture();poses,best,ledger,issues=inspect_poses(path,candidates,'fixture')
        self.assertEqual((len(poses),len(best),len(ledger),issues),(3,3,3,[]))
        self.assertTrue(all(r['run_id']=='fixture' for r in poses))

    def test_missing_corrupt_nonfinite_wrong_site_are_not_success(self):
        for mode in ('missing','corrupt','nan','wrong_site'):
            with self.subTest(mode=mode):
                path,candidates=self.fixture(mode);_,_,ledger,issues=inspect_poses(path,candidates,'fixture')
                self.assertTrue(issues)
                self.assertEqual(len(ledger),3)


if __name__=='__main__':unittest.main()
