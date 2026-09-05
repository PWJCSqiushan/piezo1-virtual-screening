from __future__ import annotations

import csv
import io
import time
import urllib.error
import urllib.request

PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    "{cids}/property/SMILES,Title,MolecularFormula,MolecularWeight/CSV"
)


def parse_property_csv(text: str) -> list[dict[str, str]]:
    rows = []
    for row in csv.DictReader(io.StringIO(text)):
        cid = (row.get("CID") or "").strip()
        smiles = (row.get("SMILES") or "").strip()
        if cid and smiles:
            rows.append({(key or "").strip(): (value or "").strip() for key, value in row.items()})
    return rows


def fetch_property_batch(
    cids: list[int], *, timeout: float = 60.0, retries: int = 3
) -> tuple[list[dict[str, str]], str]:
    if not cids:
        return [], ""
    url = PUBCHEM_PROPERTY_URL.format(cids=",".join(str(cid) for cid in cids))
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "piezo1-vs/0.3 academic-engineering-validation"},
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                charset = response.headers.get_content_charset() or "utf-8"
                return parse_property_csv(response.read().decode(charset)), url
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(1.0 * (attempt + 1))
    assert last_error is not None
    raise last_error
