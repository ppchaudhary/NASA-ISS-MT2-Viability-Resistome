#!/usr/bin/env python3
"""Stage processed input files for reproducible downstream analyses."""

from pathlib import Path
import json
import shutil
import hashlib

ROOT = Path(__file__).resolve().parents[1]

mapping = json.loads(
    (ROOT / "input_staging_manifest.json").read_text()
)

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

for target, source in mapping.items():
    src = ROOT / source
    dst = ROOT / target

    if not src.is_file():
        raise FileNotFoundError(src)

    if dst.is_symlink():
        print("EXISTING LINK:", target)
        continue

    if dst.exists():
        if sha256(src) != sha256(dst):
            raise RuntimeError(
                f"Existing file differs from reference: {target}"
            )
        print("VERIFIED:", target)
        continue

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    assert sha256(src) == sha256(dst)

    print("STAGED:", target)

print("Input staging complete.")
