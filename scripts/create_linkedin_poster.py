"""Render a shareable scientific poster from the archived studies.

Author: Biswajit Jana. No new simulations or illustrative image generation.
Run from the repository root: .venv/Scripts/python scripts/create_linkedin_poster.py
"""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outreach" / "linkedin"
OUT.mkdir(parents=True, exist_ok=True)
RATE_PATH = ROOT / "results/study-005/rates.csv"
SUMMARY_PATH = ROOT / "results/study-004/summary.json"
with RATE_PATH.open(newline="", encoding="utf-8") as stream:
    rows = list(csv.DictReader(stream))
summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
delta = summary["designs"]["broad_band"]["achievable_error_budget"]

BG, INK, MUTED = "#F7F5EF", "#142B3A", "#4B5F69"
TEAL, ORANGE, RULE = "#007E87", "#B65421", "#CFD5D1"
plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                     "svg.hashsalt": "biswajit-jana-linkedin-poster-v1"})
fig = plt.figure(figsize=(10, 12.5), dpi=180, facecolor=BG)

def text(x, y, value, size=12, color=INK, weight="normal", **kw):
    return fig.text(x, y, value, fontsize=size, color=color, weight=weight,
                    va="top", **kw)

def rule(y):
    fig.add_artist(Line2D([.07, .93], [y, y], transform=fig.transFigure,
                         color=RULE, linewidth=.8))

text(.07, .957, "COSMIC UNKNOWNS  /  RESEARCH NOTE 01", 10, TEAL, "bold")
text(.93, .957, "NOV 2026", 10, MUTED, ha="right")
text(.065, .912, "Would a wormhole\nsignal be unique?", 40, weight="bold",
     linespacing=1.08)
text(.07, .796, "A scalar-wave benchmark for future measurements", 14, MUTED)
text(.07, .763, "Biswajit Jana", 14, weight="bold")
text(.93, .759, "PUBLIC WORKING RESEARCH", 9, MUTED, ha="right")
rule(.728)

for x, number, label in [(.07, "150,000", "synthetic experiments"),
                         (.40, "25", "count budgets per band"),
                         (.73, "24", "frequency channels")]:
    text(x, .709, number, 27, weight="bold")
    text(x, .673, label, 10, MUTED)
rule(.647)
text(.07, .628, "01  /  IDEAL DETECTOR", 10, TEAL, "bold")
text(.07, .604, "More counts can separate the models.", 19, weight="bold")

ax = fig.add_axes([.135, .342, .79, .231], facecolor=BG)
for design, label, color, marker in [
    ("broad_band", "Broad band", TEAL, "o"),
    ("low_band", "Low band", ORANGE, "s"),
]:
    data = sorted([r for r in rows if r["design"] == design and
                   r["group"] == "ellis_power"], key=lambda r: int(r["N"]))
    assert len(data) == 25
    n = np.array([int(r["N"]) for r in data])
    rate = np.array([float(r["rate"]) for r in data])
    lo = np.array([float(r["wilson_low"]) for r in data])
    hi = np.array([float(r["wilson_high"]) for r in data])
    ax.errorbar(n, rate, yerr=[rate-lo, hi-rate], color=color, marker=marker,
                markersize=3.8, linewidth=1.55, elinewidth=1, capsize=2,
                label=label)
ax.set_xscale("log")
ax.set_ylim(-.025, 1.045)
ax.set_xlim(85, 120000)
ax.set_yticks([0, .25, .5, .75, 1])
ax.yaxis.set_major_formatter(PercentFormatter(1))
ax.set_xlabel("Incident quanta per channel", fontsize=11, color=INK, labelpad=7)
ax.set_ylabel("Comparator rejection rate", fontsize=11, color=INK, labelpad=9)
ax.tick_params(colors=MUTED, labelsize=10, length=3)
ax.spines[["top", "right"]].set_visible(False)
for side in ("bottom", "left"):
    ax.spines[side].set_color(RULE)
ax.grid(axis="y", color=RULE, linewidth=.6, alpha=.8)
ax.set_axisbelow(True)
ax.legend(loc="center right", frameon=False, fontsize=11, labelcolor=INK)
ax.annotate("92.1% at 3,162", xy=(3162, .921), xytext=(500, .74),
            fontsize=11, weight="bold", color=TEAL,
            arrowprops={"arrowstyle": "-", "color": TEAL, "lw": .8})
text(.07, .294, "Ellis-model data vs fitted sech² barrier · 1,000 evaluation trials per point", 9.4, MUTED)
text(.07, .277, "Bars: pointwise 95% Wilson intervals, conditional on a calibrated threshold.", 9.4, MUTED)
rule(.252)
text(.07, .235, "02  /  DETECTOR AMBIGUITY", 10, TEAL, "bold")
text(.07, .212, f"≈{100 * delta:.4f}%", 30, weight="bold")
text(.405, .210, "An allowed error budget can make a\ncompeting model match all 24 broad-band\nrecorded probabilities.", 12.5, linespacing=1.4)
text(.07, .154, "Independent false-positive and false-negative errors per channel.\nConstructive numerical witness; the smallest possible budget is not certified.",
     10, MUTED, linespacing=1.5)
rule(.111)
text(.07, .094, "Future instruments need calibration as well as sensitivity.", 13, weight="bold")
text(.07, .066, "READ THE REPORT + EXPLORE THE DATA", 9, TEAL, "bold")
text(.07, .047, "github.com/Biswajit1999/cosmic-unknowns", 11)
text(.07, .022, "MODEL STUDY · NO OBSERVATIONAL DETECTION · SPECIALIST REVIEW PENDING", 8, MUTED)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for artist in fig.texts:
    box = artist.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= 1800 and box.y1 <= 2250, artist.get_text()
png = OUT / "wormhole-signal-poster-1800x2250.png"
svg = OUT / "wormhole-signal-poster.svg"
fig.savefig(png, dpi=180, facecolor=BG, metadata={"Author": "Biswajit Jana"})
fig.savefig(svg, facecolor=BG, metadata={"Creator": "Biswajit Jana", "Date": "2026-11-10"})
svg.write_bytes(("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n").encode("utf-8"))
plt.close(fig)
manifest = {
    "author": "Biswajit Jana", "dimensions_px": [1800, 2250],
    "data_type": "synthetic-model analysis; no observational data",
    "files": {str(p.relative_to(ROOT)).replace("\\", "/"):
              hashlib.sha256(p.read_bytes()).hexdigest()
              for p in [RATE_PATH, SUMMARY_PATH, png, svg]},
}
(OUT / "manifest.json").write_bytes((json.dumps(manifest, indent=2) + "\n").encode("utf-8"))
print(f"Created {png}\nCreated {svg}")
