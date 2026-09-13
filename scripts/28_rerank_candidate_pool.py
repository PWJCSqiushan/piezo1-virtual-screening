from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from piezo_vs.io_utils import sha256_file
from piezo_vs.candidate_audit import write_csv, csv_rows, validate_ranking


def main():
    p=argparse.ArgumentParser(description='Rerank an audited finite pool with shared conformers and deduplicated identical pockets.')
    p.add_argument('--batch-dir',type=Path,required=True); a=p.parse_args()
    batch=a.batch_dir.resolve(); audit=json.loads((batch/'audit/audit_summary.json').read_text())
    if audit['exploratory_pool_status']!='eligible': raise ValueError('Audit does not allow this pool')
    compounds=batch/'audit/candidate_master.csv'
    # Bind the actual exported master to the audit's structure inventory.
    checked=batch/'audit/candidate_master_audited.csv'
    if sha256_file(compounds)!=sha256_file(checked): raise ValueError('Candidate master changed after audit')
    shared=batch/'inputs_v2/8YEZ_C006_two_plus'
    def run(script,*args):
        subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,args)],cwd=ROOT,check=True)
    aliases=[]; computed={}; statuses=[]
    for policy in ('two_plus','all_supporting_tools'):
      for pdb,site in [('8YEZ','C006'),('8YEZ','C007'),('8ZU3','C016'),('8YFC','C016'),('9VMX','C016')]:
        key=f'{pdb}_{site}_{policy}'; inp=batch/'inputs_v2'/key
        if not (inp/'input_run.json').exists():
            args=['--pdb-id',pdb,'--consensus-id',site,'--compounds',compounds,'--output-dir',inp,'--residue-policy',policy]
            if (shared/'input_run.json').exists(): args+=['--shared-molecule-inputs',shared]
            run('10_prepare_drugclip_inputs.py',*args)
        meta=json.loads((inp/'input_run.json').read_text())
        if meta['source_files']['compounds_sha256'] != sha256_file(compounds):
            raise ValueError(f'{key}: stale candidate inputs')
        for filename, hashkey in [('mols.lmdb','mols_lmdb_sha256'),('pocket.lmdb','pocket_lmdb_sha256'),('molecule_manifest.csv','molecule_manifest_sha256'),('standardized_pocket.pdb','standardized_pocket_pdb_sha256')]:
            if sha256_file(inp/filename) != meta['outputs'][hashkey]:
                raise ValueError(f'{key}: input hash mismatch: {filename}')
        for name in ('structure', 'consensus'):
            if sha256_file(ROOT/meta['source_files'][name]) != meta['source_files'][name+'_sha256']:
                raise ValueError(f'{key}: source {name} changed')
        if meta['rejected_compound_count'] or meta['accepted_compound_count']!=audit['canonical_structure_count']:
            raise ValueError(f'{key}: incomplete conformer preparation; inspect rejected_compounds.csv')
        if sha256_file(inp/'mols.lmdb')!=json.loads((shared/'input_run.json').read_text())['outputs']['mols_lmdb_sha256']:
            raise ValueError('Shared conformer mismatch')
        atom_lines=[l[12:16]+l[17:20]+l[30:54]+l[76:78] for l in (inp/'standardized_pocket.pdb').read_text().splitlines() if l.startswith('ATOM')]
        import hashlib
        signature=hashlib.sha256('\n'.join(atom_lines).encode()).hexdigest()
        representative=computed.get(signature,key)
        aliases.append(dict(pdb_id=pdb,consensus_id=site,residue_policy=policy,pocket_geometry_sha256=signature,calculation_key=representative,computed_separately=(representative==key),independent_biological_evidence=False))
        write_csv(batch/'pocket_calculation_mapping.csv',aliases)
        if signature in computed: continue
        computed[signature]=key
        out=batch/'rerank'/key
        status=dict(key=key,status='running');statuses.append(status)
        (batch/'rerank_status.json').write_text(json.dumps(statuses,indent=2))
        try:
            if not (out/'ranked_compounds.csv').exists():
                run('11_run_drugclip.py','--pdb-id',pdb,'--consensus-id',site,'--tier',meta['tier'],'--mol-lmdb',inp/'mols.lmdb','--pocket-lmdb',inp/'pocket.lmdb','--checkpoint',ROOT/'tools/drugclip/official/checkpoint_best.pt','--output-dir',out,'--no-fp16')
                run('12_normalize_drugclip_results.py','--ranked',out/'embeddings/ranked_compounds.txt','--manifest',inp/'molecule_manifest.csv','--output',out/'ranked_compounds.csv','--pdb-id',pdb,'--consensus-id',site,'--tier',meta['tier'],'--input-run',inp/'input_run.json')
            rows=csv_rows((out/'ranked_compounds.csv').read_bytes())
            validate_ranking(rows)
            if {r['compound_id'] for r in rows} != {r['compound_id'] for r in csv_rows(compounds.read_bytes())}:
                raise ValueError(f'{key}: incomplete ranked universe')
            normalized=json.loads((out/'ranked_compounds.csv.run.json').read_text())
            if normalized['output_sha256'] != sha256_file(out/'ranked_compounds.csv'):
                raise ValueError(f'{key}: modified normalized ranking')
            runmeta=json.loads((out/'run.json').read_text())
            if runmeta['status']!='success' or any(runmeta['input_sha256'][k]!=sha256_file(p) for k,p in [('mols_lmdb',inp/'mols.lmdb'),('pocket_lmdb',inp/'pocket.lmdb'),('checkpoint',ROOT/'tools/drugclip/official/checkpoint_best.pt')]):
                raise ValueError(f'{key}: run provenance mismatch')
            status['status']='computed';status['ranked_sha256']=sha256_file(out/'ranked_compounds.csv')
        except Exception as e:
            status.update(status='failed',error=str(e));raise
        finally: (batch/'rerank_status.json').write_text(json.dumps(statuses,indent=2))
    for policy in ('two_plus','all_supporting_tools'):
        run('19_compare_drugclip_runs.py','--first',batch/f'rerank/8YEZ_C006_{policy}/ranked_compounds.csv','--second',batch/f'rerank/8YEZ_C007_{policy}/ranked_compounds.csv','--top-k',20,'--output',batch/f'comparison_{policy}.json')
    return 0
if __name__=='__main__': sys.exit(main())
