from __future__ import annotations

import argparse
import hashlib
import os
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    resolved_destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as handle:
        for member in handle.getmembers():
            member_path = (destination / member.name).resolve()
            if resolved_destination not in member_path.parents and member_path != resolved_destination:
                raise RuntimeError(f"Unsafe path in P2Rank archive: {member.name}")
        handle.extractall(destination)


def main() -> int:
    parser = argparse.ArgumentParser(description="Download and verify a local P2Rank release.")
    parser.add_argument("--version", default="2.5.1")
    args = parser.parse_args()

    version = args.version
    cache_dir = ROOT / "tools" / "cache"
    tools_dir = ROOT / "tools" / "p2rank"
    archive = cache_dir / f"p2rank_{version}.tar.gz"
    destination = tools_dir / f"p2rank_{version}"
    url = f"https://github.com/rdk/p2rank/releases/download/{version}/p2rank_{version}.tar.gz"
    expected = "d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274" if version == "2.5.1" else ""

    cache_dir.mkdir(parents=True, exist_ok=True)
    tools_dir.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "piezo1-vs/0.2"})
        partial = archive.with_suffix(archive.suffix + ".part")
        with urllib.request.urlopen(request, timeout=300) as response, partial.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        partial.replace(archive)

    actual = sha256(archive)
    if expected and actual != expected:
        raise RuntimeError(f"P2Rank archive SHA256 mismatch: expected {expected}, got {actual}")
    if not destination.exists():
        safe_extract(archive, tools_dir)
    launcher = destination / ("prank.bat" if os.name == "nt" else "prank")
    if not launcher.exists():
        raise FileNotFoundError(f"P2Rank launcher was not found after extraction: {launcher}")
    if os.name != "nt":
        launcher.chmod(launcher.stat().st_mode | 0o111)
    print(f"P2Rank {version}: {destination}")
    print(f"Archive SHA256: {actual}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
