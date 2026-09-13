from __future__ import annotations
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from piezo_vs.candidate_selection import sha256_file, write_csv_rows
from piezo_vs.docking_validation import validate_prepared_inputs, inspect_poses, check_technical_evidence


def main():
    p = argparse.ArgumentParser(description='Run pinned GNINA with explicit candidates and complete input/output reconciliation.')
    p.add_argument('--gnina', type=Path, required=True)
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--exhaustiveness', type=int, default=8)
    p.add_argument('--num-modes', type=int, default=9)
    p.add_argument('--cpu', type=int, default=8)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--gpu-id', type=int, default=0)
    p.add_argument('--no-gpu', action='store_true')
    p.add_argument('--stage', choices=('technical_check', 'batch'), default='technical_check')
    p.add_argument('--technical-check-run', type=Path, action='append', default=[])
    a = p.parse_args()
    if min(a.exhaustiveness, a.num_modes, a.cpu) < 1:
        raise ValueError('Positive exhaustiveness, num-modes, cpu required')
    output = a.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Refusing to overwrite a GNINA run')
    output.mkdir(parents=True, exist_ok=True)
    meta = {'schema_version': 3, 'run_id': output.name, 'status': 'failed', 'stage': a.stage,
            'started_at_utc': datetime.now(timezone.utc).isoformat()}
    try:
        config, candidates = validate_prepared_inputs(a.input_dir)
        binary = a.gnina.resolve()
        pinned = json.loads((ROOT/'config/tool_sources.json').read_text())['tools']['gnina']
        binary_hash = sha256_file(binary)
        if binary_hash != pinned['binary_sha256'] or binary.stat().st_size != int(pinned['binary_size_bytes']):
            raise ValueError('GNINA binary size/hash mismatch')
        receptor, ligands = a.input_dir/'receptor.pdb', a.input_dir/'ligands.sdf'
        receptor_hash = sha256_file(receptor)
        if a.stage == 'technical_check' and len(candidates) > 3:
            raise ValueError('Technical checks contain at most 3 molecules; use reviewed batches after checks pass')
        if a.stage == 'batch':
            check_technical_evidence(a.technical_check_run, config, receptor_hash, binary_hash, {r['compound_id'] for r in candidates})
        version = subprocess.run([str(binary), '--version'], capture_output=True, text=True, check=True)
        command = [str(binary), '--receptor', str(receptor.resolve()), '--ligand', str(ligands.resolve())]
        for key in ('center_x','center_y','center_z','size_x','size_y','size_z'):
            command += ['--'+key, str(config['box'][key])]
        command += ['--out', str(output/'docked.sdf'), '--cnn_scoring', 'rescore', '--pose_sort_order', 'CNNscore',
                    '--exhaustiveness', str(a.exhaustiveness), '--num_modes', str(a.num_modes), '--cpu', str(a.cpu), '--seed', str(a.seed)]
        command += ['--no_gpu'] if a.no_gpu else ['--device', str(a.gpu_id)]
        meta.update(status='running', purpose=config['purpose'], pdb_id=config['pdb_id'], consensus_id=config['consensus_id'],
                    box=config['box'], candidate_ids=[r['compound_id'] for r in candidates], requested_count=len(candidates),
                    gnina_sha256=binary_hash, gnina_version=(version.stdout or version.stderr).strip(), command=command,
                    input_sha256={'receptor_pdb':receptor_hash,'ligands_sdf':sha256_file(ligands),
                                  'docking_config':sha256_file(a.input_dir/'docking_config.json'),
                                  'input_run':sha256_file(a.input_dir/'input_run.json')},
                    technical_evidence=[{'path':str(p.resolve()),'sha256':sha256_file(p)} for p in a.technical_check_run],
                    scientific_status='computational_docking_prediction_only')
        (output/'run.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        (output/'stdout.log').write_text(process.stdout,encoding='utf-8')
        (output/'stderr.log').write_text(process.stderr,encoding='utf-8')
        poses, best, ledger, issues = inspect_poses(output/'docked.sdf', candidates, output.name)
        if process.returncode:
            issues.append('gnina_nonzero_exit:'+str(process.returncode))
        base_fields = list(candidates[0]) + ['run_id','pose_index','minimizedAffinity','CNNscore','CNNaffinity','interpretation']
        write_csv_rows(output/'all_poses.csv',poses,base_fields)
        write_csv_rows(output/'best_poses.csv',best,base_fields)
        write_csv_rows(output/'candidate_results.csv',ledger)
        meta.update(status='success' if not issues else 'partial' if best else 'failed', returncode=process.returncode,
                    issues=issues, computed_count=len(best), outputs={
                        'docked_sdf_sha256':sha256_file(output/'docked.sdf') if (output/'docked.sdf').exists() else None,
                        'candidate_results_sha256':sha256_file(output/'candidate_results.csv'),
                        'all_poses_csv_sha256':sha256_file(output/'all_poses.csv'),
                        'best_poses_csv_sha256':sha256_file(output/'best_poses.csv'),
                        'pose_count':len(poses),'compound_count':len(best)})
    except Exception as exc:
        meta.update(status='failed', error=str(exc))
    finally:
        meta['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        (output/'run.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta,indent=2))
    return 0 if meta['status']=='success' else 2


if __name__ == '__main__':
    raise SystemExit(main())
