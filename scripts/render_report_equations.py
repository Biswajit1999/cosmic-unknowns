"""Render display mathematics for the PDF. Author: Biswajit Jana."""
import hashlib
import json
from pathlib import Path
import re
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["mathtext.fontset"] = "stix"
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"reports/007-calibration-limited-identification.md"
OUT = ROOT/"tmp/pdfs/equations"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    equations = re.findall(r"\$\$\s*(.*?)\s*\$\$", SOURCE.read_text(encoding="utf-8"), flags=re.S)
    manifest = []
    for i, equation in enumerate(equations, 1):
        clean = " ".join(equation.split())
        path = OUT/f"equation-{i:02d}.png"
        mathtext.math_to_image("$"+clean+"$", path, prop=FontProperties(family="STIXGeneral", size=13), dpi=300, format="png", color="black")
        manifest.append({"number": i, "latex": clean, "path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    (OUT/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(f"Rendered {len(manifest)} equations")


if __name__ == "__main__":
    main()
