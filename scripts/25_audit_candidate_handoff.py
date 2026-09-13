from pathlib import Path
import argparse, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from piezo_vs.candidate_audit import audit_zip
if __name__=='__main__':
    p=argparse.ArgumentParser(description='Audit a finite web-returned candidate pool without executing or extracting handoff code.')
    p.add_argument('--handoff-zip',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--structure-dir',type=Path,default=ROOT/'data/raw/structures')
    a=p.parse_args(); result=audit_zip(a.handoff_zip,a.output_dir,a.structure_dir)
    print(json.dumps(result,ensure_ascii=False,indent=2)); sys.exit(0 if result['exploratory_pool_status']=='eligible' else 2)
