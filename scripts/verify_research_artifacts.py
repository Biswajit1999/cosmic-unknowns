"""Verify saved research products and local document links. Author: Biswajit Jana."""
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = {}
    for manifest_path in sorted((ROOT/"results").glob("study-*/sha256.json")):
        manifest.update(json.loads(manifest_path.read_text(encoding="utf-8")))
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
    verify_inference()
    camera = json.loads((ROOT/"results/study-003/summary.json").read_text())
    assert camera["author"] == "Biswajit Jana"
    assert camera["max_angle_quadrature_error"] < 1e-8
    assert camera["uniform_radiance_max_error"] < 1e-12
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


def verify_inference():
    import numpy as np
    from transmission_inference import pt_transmission, deviance, wilson
    folder = ROOT/"results/study-002"
    def read_csv(name):
        with (folder/name).open() as f:
            return list(csv.DictReader(f))
    spectra = read_csv("spectra.csv")
    trials = read_csv("trials.csv")
    examples = read_csv("example-counts.csv")
    config = json.loads((folder/"design.json").read_text())
    assert len(trials) == 12000 and len(spectra) == len(examples) == 48
    rng = np.random.default_rng(config["seed"])
    grouped = {}
    for r in trials:
        s = [v for v in spectra if v["design"] == r["design"]]
        k = np.array([float(v["K"]) for v in s])
        p = np.array([float(v["T_ellis"] if r["group"] == "ellis_power" else v["T_pt_projection"]) for v in s])
        n = int(r["N"])
        y = np.array([int(r[f"count_{i:02d}"]) for i in range(24)])
        assert np.array_equal(y, rng.binomial(n, p)), "Saved seed/count mismatch"
        q = pt_transmission(k, [float(r["logC"]), float(r["logb"])])
        assert abs(deviance(y, n, q)-float(r["deviance"])) < 1e-7
        grouped.setdefault((r["design"], n, r["group"]), []).append(float(r["deviance"]))
    for r in read_csv("rates.csv"):
        d, n, g = r["design"], int(r["N"]), r["group"]
        threshold = float(np.quantile(grouped[(d, n, "calibration")], .95, method="higher"))
        rejected = int(sum(np.array(grouped[(d, n, g)]) > threshold))
        assert rejected == int(r["rejected"])
        np.testing.assert_allclose(wilson(rejected, 500), [float(r["wilson_low"]), float(r["wilson_high"])])
    for r in examples:
        n, y = int(r["N"]), int(r["transmitted"])
        assert y+int(r["reflected"]) == n
        np.testing.assert_allclose(wilson(y, n), [float(r["wilson_low"]), float(r["wilson_high"])])
    audit = json.loads((folder/"independent-audit.json").read_text())
    assert audit["original_non_success_audited"] == sum(r["optimizer_success"] == "0" for r in trials)
    assert audit["alternative_non_success"] == 0
    print("Replayed 288,000 synthetic counts; verified all deviances, rates, intervals and fit-warning audit.")


if __name__ == "__main__":
    main()
