from __future__ import annotations
import csv, hashlib, io, json, math, zipfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from rdkit import Chem, rdBase


def digest(data):
    return hashlib.sha256(data).hexdigest()


def csv_rows(data):
    return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'), newline='')))


def write_csv(path, rows, fields=None):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def validate_ranking(rows):
    ids, scores = set(), []
    for expected, row in enumerate(rows, 1):
        cid = row.get('compound_id', '')
        if not cid or cid in ids:
            raise ValueError('missing_or_duplicate_id')
        ids.add(cid)
        if int(row['rank']) != expected:
            raise ValueError('non_contiguous_rank')
        score = float(row['drugclip_score'])
        if not math.isfinite(score):
            raise ValueError('nonfinite_score')
        scores.append(score)
    if any(a < b for a,b in zip(scores, scores[1:])):
        raise ValueError('non_descending_scores')


def audit_zip(archive, output, structure_dir=None):
    archive, output = Path(archive), Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Refusing to overwrite nonempty audit output')
    output.mkdir(parents=True, exist_ok=True)
    errors, warnings, hashes, tasks, records = [], [], [], [], []
    identities, entities = defaultdict(set), defaultdict(set)
    pocket_mapping = []
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        if len(names) != len(set(names)):
            raise ValueError('duplicate_zip_members')
        for name in names:
            p = PurePosixPath(name.replace('\\', '/'))
            if p.is_absolute() or '..' in p.parts or ':' in name:
                raise ValueError('unsafe_zip_path')
        covered = set()
        for line in z.read('SHA256SUMS.txt').decode('utf-8-sig').splitlines():
            if not line.strip(): continue
            expected, name = line.split(maxsplit=1); name = name.lstrip('*')
            if name in covered: errors.append('duplicate_hash_entry:'+name)
            covered.add(name)
            actual = digest(z.read(name)) if name in names else ''
            ok = actual == expected.lower()
            hashes.append(dict(path=name,expected_sha256=expected,actual_sha256=actual,passed=ok))
            if not ok: errors.append('hash_mismatch:'+name)
        task_names = sorted({n.split('/')[1] for n in names if n.startswith('tasks/') and n.endswith('/ranked_compounds.csv')})
        if not task_names: errors.append('no_tasks')
        for task in task_names:
            prefix = 'tasks/'+task+'/'
            before = len(errors)
            required = [prefix+s for s in ('ranked_compounds.csv','molecule_manifest.csv','input_run.json')]
            for name in required:
                if name not in covered: errors.append('unhashed_required_file:'+name)
            rows = csv_rows(z.read(required[0])); manifest = csv_rows(z.read(required[1])); meta = json.loads(z.read(required[2]))
            try: validate_ranking(rows)
            except (KeyError,TypeError,ValueError) as e: errors.append(task+':'+str(e))
            original_name = prefix+meta['server']['original_export_file']
            if original_name not in covered: errors.append('unhashed_original:'+original_name)
            original = csv_rows(z.read(original_name))
            if len(rows)!=len(original) or len(rows)!=len(manifest): errors.append(task+':row_count_mismatch')
            missing = Counter()
            for i,row in enumerate(rows):
                for field in ('compound_id','entity_id','smiles','drugclip_score','pdb_id','consensus_id','job_id','library_name'):
                    if not row.get(field): errors.append(f'{task}:{i+1}:missing_{field}')
                for field,value in row.items():
                    if value in ('',None): missing[field]+=1
                if row.get('pdb_id') != meta['pdb_id'] or row.get('consensus_id') != meta['consensus_id'] or row.get('job_id') != meta['server']['job_id']:
                    errors.append(f'{task}:{i+1}:task_identity_mismatch')
                if i<len(original):
                    for key,rawkey in [('compound_id','mol_id'),('entity_id','entity_id'),('smiles','smiles'),('drugclip_score','drugclip_score'),('vina_docking_score','vina_docking_score'),('library_name','library_name')]:
                        if row.get(key,'') != original[i].get(rawkey,''): errors.append(f'{task}:{i+1}:original_mismatch_{key}')
                if i<len(manifest):
                    for key in ('compound_id','smiles','entity_id','job_id','pdb_id','consensus_id'):
                        if row.get(key,'') != manifest[i].get(key,''): errors.append(f'{task}:{i+1}:manifest_mismatch_{key}')
                mol = Chem.MolFromSmiles(row.get('smiles',''))
                canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) if mol else ''
                if not canonical: errors.append(f'{task}:{i+1}:invalid_smiles')
                identities[row.get('compound_id','')].add(canonical)
                entities[row.get('entity_id','')].add(canonical)
                records.append(dict(row, raw_original_row_json=json.dumps(original[i] if i<len(original) else {},ensure_ascii=False,sort_keys=True), task_key=task, source_row=i+1, source_file=original_name, canonical_smiles=canonical, structure_sha256=digest(canonical.encode()), normalization_changed=(canonical!=row.get('smiles')), score_kind='unknown'))
            receptor_status = 'not_checked'
            if structure_dir:
                receptor = Path(structure_dir)/(meta['pdb_id']+'.pdb')
                uploaded = prefix+'uploaded_original/'+meta['upload']['filename']
                if uploaded not in names: uploaded = prefix+meta['upload']['filename']
                if uploaded not in covered: errors.append(task+':unhashed_upload')
                if uploaded not in names or not receptor.is_file() or digest(z.read(uploaded))!=digest(receptor.read_bytes()):
                    errors.append(task+':receptor_mismatch'); receptor_status='failed'
                else:
                    receptor_status='passed'
                    def atom_key(line):
                        return (line[30:54], line[76:78].strip())
                    receptor_atoms = defaultdict(list)
                    for line in receptor.read_text(encoding='ascii',errors='replace').splitlines():
                        if line.startswith(('ATOM  ', 'HETATM')):
                            receptor_atoms[atom_key(line)].append(line)
                    pocket_names = [n for n in names if n.startswith(prefix+'processed_pocket/') and n.endswith('.pdb')]
                    if len(pocket_names)!=1:
                        errors.append(task+':missing_or_ambiguous_downloaded_pocket')
                    for pocket_name in pocket_names:
                        if pocket_name not in covered: errors.append(task+':unhashed_downloaded_pocket')
                        for line in z.read(pocket_name).decode('ascii',errors='replace').splitlines():
                            if not line.startswith(('ATOM  ','HETATM')): continue
                            matches = receptor_atoms[atom_key(line)]
                            if len(matches)!=1: errors.append(task+':unmapped_or_ambiguous_pocket_atom:'+line[6:11].strip())
                            match = matches[0] if len(matches)==1 else ''
                            pocket_mapping.append(dict(task_key=task, downloaded_pocket=pocket_name,
                                pocket_serial=line[6:11].strip(), mapped_count=len(matches), receptor_chain=match[21:22],
                                receptor_residue=match[22:27].strip(), receptor_atom_name=match[12:16].strip()))
            tasks.append(dict(task_key=task,pdb_id=meta['pdb_id'],consensus_id=meta['consensus_id'],record_count=len(rows),receptor_status=receptor_status,record_checks='passed' if before==len(errors) else 'failed',missing_fields=json.dumps(missing),source_rank_sha256=digest(z.read(required[0]))))
        for cid,structures in identities.items():
            if len(structures)>1: errors.append('compound_id_structure_conflict:'+cid)
        for eid,structures in entities.items():
            if len(structures)>1: errors.append('entity_id_structure_conflict:'+eid)
        uncovered = sorted(set(n for n in names if not n.endswith('/') and n!='SHA256SUMS.txt') - covered)
        if uncovered: warnings.append('uncovered_nonrequired_files:'+str(len(uncovered)))
        structures = defaultdict(list)
        for r in records: structures[r['canonical_smiles']].append(r)
        master = []
        for smiles, sources in structures.items():
            cid = 'MOL_'+digest(smiles.encode())[:16]
            for r in sources: r['pool_compound_id']=cid
            master.append(dict(compound_id=cid,canonical_smiles=smiles,structure_sha256=digest(smiles.encode()),source='finite_web_returned_pool',source_ids=';'.join(sorted({r['compound_id'] for r in sources})),source_tasks=';'.join(sorted({r['task_key'] for r in sources})),source_record_count=len(sources)))
        master.sort(key=lambda r:r['compound_id'])
        if pocket_mapping: write_csv(output/'downloaded_pocket_atom_mapping.csv', pocket_mapping)
        write_csv(output/'hash_checks.csv', hashes)
        write_csv(output/'task_table.csv', tasks)
        write_csv(output/'source_records.csv', records)
        write_csv(output/'candidate_master_audited.csv',master)
        write_csv(output/'source_mapping.csv',[{k:r[k] for k in ('task_key','source_row','compound_id','entity_id','pool_compound_id','canonical_smiles','normalization_changed','structure_sha256','source_file')} for r in records])
        integrity = all(h['passed'] for h in hashes) and bool(hashes)
        summary = dict(schema_version=1,archive_sha256=digest(archive.read_bytes()),archive_path=str(archive.resolve()),hash_count=len(hashes),uncovered_file_count=len(uncovered),downloaded_pocket_atom_count=len(pocket_mapping),downloaded_pocket_unmapped_count=sum(r['mapped_count']!=1 for r in pocket_mapping),record_count=len(records),source_id_count=len(identities),canonical_structure_count=len(master),file_integrity='passed' if integrity else 'failed',exploratory_pool_status='eligible' if not errors else 'blocked',platform_pipeline_reproduction='not_reproduced',formal_gnina_status='pending_step8_and_human_review',normalization='RDKit canonical isomeric SMILES; no desalting, neutralization or tautomer merging',rdkit_version=rdBase.rdkitVersion,errors=errors,warnings=warnings,web_score_kind='unknown',score_threshold_applied=False)
        (output/'audit_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
        if not errors: write_csv(output/'candidate_master.csv',master)
        return summary
