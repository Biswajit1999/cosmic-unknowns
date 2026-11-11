# Research poster

**Author: Biswajit Jana**

- [Portrait PNG, 1800 × 2250](wormhole-signal-poster-1800x2250.png): ready to attach to a LinkedIn post.
- [Editable vector SVG](wormhole-signal-poster.svg): the same composition with text preserved.
- [Suggested first-person caption and alt text](post-caption.md).

The poster is drawn directly from archived numerical results, using Matplotlib. It contains no generated astronomical imagery or observational claims. Its two results assume different detector models: ideal counting in the chart, and flexible independent classification errors in the matching construction. The error allowance is achievable for the reported witness, not a certified global minimum. The poster is a summary, not a substitute for the full uncertainty analysis and prior-work discussion in the [report](../../reports/007-calibration-limited-identification.md).

Sources: [study 005 rates](../../results/study-005/rates.csv) and [study 004 summary](../../results/study-004/summary.json). The figure omits the separate false-positive and threshold-sensitivity panels to keep the social poster readable; both remain available in the full report and data. The quoted 150,000 experiments include calibration, null-check and power trials across both bands.

Rebuild from the repository root:

```powershell
.\.venv\Scripts\python.exe scripts/create_linkedin_poster.py
```

The [manifest](manifest.json) records source and output hashes. No LinkedIn post has been submitted automatically.
