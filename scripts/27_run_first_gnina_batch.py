from __future__ import annotations
import argparse
import json
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from piezo_vs.candidate_selection import read_csv_rows, sha256_file


def main():
    p=argparse.ArgumentParser(description="Run reviewed technical checks, then the explicit first GNINA batch.")
    p.add_argument('--assessment-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--receptor-review',type=Path,required=True)
    p.add_argument('--gnina',type=Path,required=True)
    a=p.parse_args()
    out=a.output_dir.resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError('Refusing to overwrite a batch')
    out.mkdir(parents=True,exist_ok=True)
    summary=dict(status='not_started',requested_count=20,computed_count=0,runs=[])
    try:
        assessment=a.assessment_dir.resolve()
        selection=json.loads((assessment/'selection_manifest.json').read_text())
        if sha256_file(assessment/'selected_candidates.csv')!=selection['selection_sha256']:
            raise ValueError('Selected candidate file changed')
        candidates=read_csv_rows(assessment/'selected_candidates.csv')
        summary['selected_count']=len(candidates)
        if len(candidates)<3:
            raise ValueError('Fewer than 3 evaluated, approved candidates; Step 8 and review must finish first')
        review=json.loads(a.receptor_review.read_text(encoding='utf-8-sig'))
        evidence=[]
        def run(script,*args):
            command=[sys.executable,str(ROOT/'scripts'/script),*map(str,args)]
            result=subprocess.run(command,capture_output=True,text=True)
            tag=str(len(summary['runs'])+1).zfill(2)
            (out/(tag+'.stdout.log')).write_text(result.stdout,encoding='utf-8')
            (out/(tag+'.stderr.log')).write_text(result.stderr,encoding='utf-8')
            summary['runs'].append(dict(command=command,returncode=result.returncode))
            if result.returncode: raise RuntimeError('Stage failed: '+script+'; see '+tag+' logs')
        for stage,prefix in [('technical_check','technical'),('batch','selected')]:
            for pdb,site in [('8YEZ','C006'),('8ZU3','C016')]:
                key=pdb+'_'+site
                candidate_file=assessment/(prefix+'_'+key+'.csv')
                if not read_csv_rows(candidate_file): continue
                r=review.get(key,{})
                if r.get('status')!='approved' or any(not r.get(k) for k in ('reviewer','reviewed_at_utc','note','heteroatom_policy')):
                    raise ValueError('Receptor environment review missing for '+key)
                inp=out/(prefix+'_'+key+'_inputs');result_dir=out/(prefix+'_'+key+'_results')
                run('15_prepare_gnina_inputs.py','--pdb-id',pdb,'--consensus-id',site,'--candidate-list',candidate_file,
                    '--purpose','formal_screening','--output-dir',inp,'--heteroatom-policy',r['heteroatom_policy'],
                    '--receptor-review-approved','--receptor-review-note',r['reviewer']+' '+r['reviewed_at_utc']+' '+r['note'])
                flags=[]
                if stage=='batch':
                    for path in evidence: flags+=['--technical-check-run',path]
                run('16_run_gnina.py','--gnina',a.gnina.resolve(),'--input-dir',inp,'--output-dir',result_dir,'--stage',stage,*flags)
                if stage=='technical_check': evidence.append(result_dir/'run.json')
                else:
                    result=json.loads((result_dir/'run.json').read_text())
                    summary['computed_count']+=result['computed_count']
        summary['status']='complete' if summary['computed_count']==20 else 'partial'
    except Exception as exc:
        summary['status']='partial' if summary['computed_count'] else 'not_started' if not summary['runs'] else 'failed'
        summary['reason']=str(exc)
    (out/'batch_status.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    return 0 if summary['status']=='complete' else 2


if __name__=='__main__':raise SystemExit(main())
