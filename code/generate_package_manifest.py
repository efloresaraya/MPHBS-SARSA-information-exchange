from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "FILE_MANIFEST_SHA256.txt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


files = sorted(
    path for path in ROOT.rglob("*")
    if path.is_file() and path != OUTPUT and path.name != ".DS_Store"
)

with OUTPUT.open("w", encoding="utf-8", newline="\n") as handle:
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        handle.write(f"{sha256(path)}  {relative}\n")

print(f"Wrote {OUTPUT} with {len(files)} entries.")
