"""Portable line endings and artifact hashes. Author: Biswajit Jana."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def write_manifest():
    folders = (ROOT/"results/study-001", ROOT/"figures/study-001")
    paths = [p for folder in folders for p in sorted(folder.glob("*"))
             if p.name != "sha256.json" and p.suffix in (".csv", ".json", ".svg", ".png")]
    for path in paths:
        if path.suffix != ".png":
            path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n"))
    hashes = {str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in paths}
    (folders[0]/"sha256.json").write_text(json.dumps(hashes, indent=2)+"\n", encoding="utf-8", newline="\n")
    return hashes


if __name__ == "__main__":
    print(f"Recorded {len(write_manifest())} portable artifact hashes.")
