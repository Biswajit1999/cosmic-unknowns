"""Verify saved research products and local document links. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest_path = ROOT/"results/study-001/sha256.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for relative, expected in manifest.items():
        path = ROOT/relative
        assert path.is_file(), f"Missing artifact: {relative}"
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, f"Hash mismatch: {relative}"
    summary = json.loads((ROOT/"results/study-001/summary.json").read_text(encoding="utf-8"))
    assert summary["author"] == "Biswajit Jana"
    with (ROOT/"results/study-001/spectra.csv").open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 525
    assert max(abs(float(r["T"])+float(r["R"])-1) for r in rows) < 1e-8
    documents = list(ROOT.glob("*.md"))
    for folder in ("reports", "notes", "references", "results"):
        documents.extend((ROOT/folder).rglob("*.md"))
    count = 0
    for path in documents:
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if target.startswith(("https:", "http:", "#", "mailto:")):
                continue
            target = target.split("#", 1)[0]
            assert (path.parent/target).exists(), f"Broken local link: {path.name} -> {target}"
            count += 1
    print(f"Verified {len(manifest)} artifact hashes, 525 spectrum rows, and {count} local links.")


if __name__ == "__main__":
    main()
