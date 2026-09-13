from __future__ import annotations
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from piezo_vs.candidate_selection import (
    assess_candidate_rows, select_candidates, validate_ranking_rows, validate_master_rows,
    maxmin_select, _external_result_state, write_external_input_bundle, read_csv_rows,
    write_csv_rows, sha256_file,
)
from piezo_vs.candidate_workflow import import_predictions


def bit(index):
    return ('0' * index + '1').ljust(2048, '0')


class SelectionTests(unittest.TestCase):
    def test_nan_wrong_rank_ties_and_missing_universe_block(self):
        base = [dict(compound_id='A', rank=1, drugclip_score=2), dict(compound_id='B', rank=2, drugclip_score=1)]
        self.assertEqual(len(validate_ranking_rows(base, 'C016', master_ids={'A','B'})), 2)
        for rows in ([dict(base[0],drugclip_score='nan'),base[1]], [base[0],dict(base[1],rank=3)],
                     [base[0],dict(base[1],drugclip_score=3)], [base[0]],
                     [dict(base[0],compound_id='B'),dict(base[1],compound_id='A',drugclip_score=2)]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                validate_ranking_rows(rows,'C016',master_ids={'A','B'})

    def test_review_alias_conflict_is_rejected(self):
        from piezo_vs.candidate_selection import _review_decision
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            _review_decision(dict(approval_status='pending',approved='true'))

    def test_copied_ranking_cannot_change_site(self):
        from piezo_vs.candidate_workflow import verified_ranking
        directory=Path(tempfile.mkdtemp(prefix='ranking-site-test-'))
        path=directory/'ranked.csv'
        write_csv_rows(path,[dict(compound_id='A',rank=1,drugclip_score=1)])
        path.with_suffix('.csv.run.json').write_text(json.dumps(dict(pdb_id='8ZU3',consensus_id='C016',tier='T2',output_sha256=sha256_file(path))))
        with self.assertRaisesRegex(ValueError,'site/tier mismatch'):
            verified_ranking(path,directory/'master.csv','8YEZ','C006')

    def test_duplicate_structure_rejected(self):
        with self.assertRaises(ValueError):
            validate_master_rows([dict(compound_id='A',canonical_smiles='CC'),dict(compound_id='B',canonical_smiles='CC')])

    def test_status_only_and_pending_predictions_never_present(self):
        for row in ({'status':'pending'}, {'status':'completed'}, {'status':'success','comment':'done'},
                    {'status':'success','prediction_payload_json':'{}'}, {'status':'failed','value':'1'}):
            self.assertNotEqual(_external_result_state(row),'present')

    def test_maxmin_distance_and_determinism(self):
        rows = [dict(compound_id=cid,rank=n,drugclip_score=4-n,morgan_fingerprint=fp)
                for cid,n,fp in [('A',1,bit(0)),('B',2,bit(0)),('C',3,bit(1))]]
        self.assertEqual([r['compound_id'] for r in maxmin_select(rows,2)],['A','C'])
        self.assertEqual(maxmin_select(rows,2),maxmin_select(list(reversed(rows)),2))

    def test_nonempty_twenty_unique_and_partial(self):
        rows=[]
        for i in range(30):
            rows.append(dict(compound_id=f'M{i:02}',canonical_smiles='C'*(i+1),eligible_for_selection=True,
                assessment_state='eligible',human_approval_status='approved',swissadme_result_state='present',
                admetlab3_result_state='present',protox3_result_state='present',rank_C006_C007=i+1,
                drugclip_score_C006_C007=30-i,rank_C016=i+1,drugclip_score_C016=30-i,morgan_fingerprint=bit(i)))
        chosen, summary=select_candidates(rows)
        again, _=select_candidates(list(reversed(rows)))
        self.assertEqual(len({r['compound_id'] for r in chosen}),20)
        self.assertEqual([r['compound_id'] for r in chosen],[r['compound_id'] for r in again])
        self.assertEqual(summary['status'],'complete')
        _, partial=select_candidates(rows[:12])
        self.assertEqual((partial['status'],partial['selected_count']),('partial',12))
        for row in rows: row.update(eligible_for_selection='false',assessment_state='pending')
        self.assertEqual(select_candidates(rows)[0],[])

    def test_swissadme_batches_preserve_order_and_size(self):
        rows=[dict(compound_id=f'M{i}',canonical_smiles='C'*(i+1)) for i in range(342)]
        directory=Path(tempfile.mkdtemp(prefix='swiss-batch-test-'))
        result=write_external_input_bundle(rows,directory)
        self.assertEqual([len(Path(p).read_text().splitlines()) for p in result['swissadme_batches']],[200,142])
        self.assertEqual(Path(result['swissadme_batches'][1]).read_text().splitlines()[0].split()[-1],'M200')

    def test_approval_requires_exact_evidence(self):
        master=[dict(compound_id='A',canonical_smiles='CC')]
        descriptor=dict(rdkit_status='parsed',morgan_fingerprint=bit(0),structure_flags='')
        external={t:[dict(compound_id='A',canonical_smiles='CC',result_status='completed',
                    prediction_payload_json='{"value": "0"}',raw_file_sha256='a'*64,raw_row_number=1)]
                  for t in ('swissadme','admetlab3','protox3')}
        with patch('piezo_vs.candidate_selection.assess_smiles',return_value=descriptor):
            pending=assess_candidate_rows(master,external_rows_by_tool=external)[0]
            review=dict(compound_id='A',canonical_smiles='CC',approval_status='approved',reviewer='fixture',
                        reviewed_at_utc='2026-01-01T00:00:00Z',assessment_comment='synthetic fixture',
                        assessment_evidence_sha256=pending['assessment_evidence_sha256'])
            accepted=assess_candidate_rows(master,external_rows_by_tool=external,review_rows=[review])[0]
            self.assertTrue(accepted['eligible_for_selection'])
            external['protox3'][0]['prediction_payload_json']='{"value":"1"}'
            with self.assertRaisesRegex(ValueError,'changed assessment'):
                assess_candidate_rows(master,external_rows_by_tool=external,review_rows=[review])

    def test_raw_import_hash_nan_unknown_ids(self):
        directory=Path(tempfile.mkdtemp(prefix='prediction-import-test-'))
        master=directory/'master.csv'
        write_csv_rows(master,[dict(compound_id='A',canonical_smiles='CC')])
        raw=directory/'raw.csv'; write_csv_rows(raw,[dict(Name='A',prediction='0')])
        spec=dict(path='raw.csv',sha256=sha256_file(raw),id_column='Name',prediction_columns=['prediction'],prediction_types={'prediction':'number'})
        manifest=dict(tool='protox3',version='synthetic-test',retrieved_at_utc='2026-01-01T00:00:00Z',
                      input_master_sha256=sha256_file(master),files=[spec])
        path=directory/'manifest.json';path.write_text(json.dumps(manifest))
        self.assertEqual(import_predictions(path,master,directory/'archive')[1][0]['result_status'],'completed')
        write_csv_rows(raw,[dict(Name='A',prediction='nan')])
        with self.assertRaisesRegex(ValueError,'hash mismatch'):import_predictions(path,master,directory/'archive2')
        spec['sha256']=sha256_file(raw);path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'Nonfinite'):import_predictions(path,master,directory/'archive3')
        for value in ('N/A', 'error', '-'):
            write_csv_rows(raw,[dict(Name='A',prediction=value)]);spec['sha256']=sha256_file(raw);path.write_text(json.dumps(manifest))
            self.assertEqual(import_predictions(path,master,directory/'placeholder')[1][0]['result_status'],'failed')
        write_csv_rows(raw,[dict(Name='A',prediction='1',status='failed')]);spec['sha256']=sha256_file(raw);path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,'status columns'):import_predictions(path,master,directory/'ignored_status')


if __name__=='__main__': unittest.main()
