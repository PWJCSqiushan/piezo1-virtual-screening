from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from piezo_vs.io_utils import sha256_file, write_json
from piezo_vs.pubchem import fetch_property_batch


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch a deterministic PubChem CID slice for engineering-scale DrugCLIP validation. "
            "This is not a disease-specific or formal screening library."
        )
    )
    parser.add_argument("--start-cid", type=int, default=1)
    parser.add_argument("--count", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "runs" / "validation" / "pubchem_cid_1_128.csv",
    )
    args = parser.parse_args()
    if args.start_cid < 1 or args.count < 1:
        raise SystemExit("start-cid and count must be positive")
    if not 1 <= args.batch_size <= 100:
        raise SystemExit("batch-size must be between 1 and 100")

    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    end_cid = args.start_cid + args.count
    records: list[dict[str, str]] = []
    request_urls: list[str] = []
    for first in range(args.start_cid, end_cid, args.batch_size):
        cids = list(range(first, min(first + args.batch_size, end_cid)))
        rows, url = fetch_property_batch(cids)
        request_urls.append(url)
        for row in rows:
            cid = row["CID"]
            records.append(
                {
                    "compound_id": f"pubchem_cid_{cid}",
                    "smiles": row["SMILES"],
                    "source": "PubChem PUG REST",
                    "source_id": f"CID:{cid}",
                    "source_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                    "title": row.get("Title", ""),
                    "molecular_formula": row.get("MolecularFormula", ""),
                    "molecular_weight": row.get("MolecularWeight", ""),
                }
            )

    if not records:
        raise RuntimeError("PubChem returned no usable CID/SMILES rows")
    fields = list(records[0])
    with output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    metadata = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "engineering_scale_smoke_test_not_formal_candidate_library",
        "selection_rule": f"consecutive PubChem CIDs {args.start_cid} through {end_cid - 1}",
        "requested_cid_count": args.count,
        "returned_row_count": len(records),
        "request_urls": request_urls,
        "output_sha256": sha256_file(output),
        "usage_note": (
            "This deterministic slice validates data preparation and model throughput only. "
            "Its ranking must not be reported as a PIEZO1 or COPD candidate conclusion."
        ),
    }
    write_json(output.with_suffix(output.suffix + ".run.json"), metadata)
    print(
        f"PUBCHEM_SMOKE_LIBRARY_OK rows={len(records)} "
        f"cid_range={args.start_cid}-{end_cid - 1} output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
