# Study 008: independent validation data

**Author: Biswajit Jana. All count records are synthetic.**

[Scientific report](../../reports/010-validation-restores-identification.md) · [Current literature context](../../references/validation-and-observability-2026-11-22.md)

- `design.json`: all 29 cases, spectra, direction, parameter domain, confidence budget, environment versions and source hashes. Python source hashes normalize CRLF to LF for portability; data hashes use raw bytes.
- `case-00.csv.gz` through `case-28.csv.gz`: 1,000 experiments each. Columns include group, replicate, drift interval, 24 science counts `y0..23`, reference false-negative/positive counts `rf0..23` and `rg0..23`, validation counts `vf0..23` and `vg0..23`, then four procedure statuses and node counts.
- `rates.csv`: 232 rows with rejection counts, pointwise 95% Wilson intervals, transfer-validity flags and incident budgets. Procedures share samples; rows are not independent experiments.
- `coverage.csv`: 58 group/case rows for reference-to-validation drift-interval coverage and descriptive median bounds. This does not measure validation-to-science transfer.
- `bound-histogram.csv`: 81 fixed bins, including overflow, for the specified representative upper-bound distribution. Fractions and pointwise Wilson endpoints are saved even for bins outside the displayed range.
- `summary.json`, `sha256.json`, `verification.json`: totals, artifact hashes and independent replay checks. Verification output is not self-hashed.

M and V are per known class per channel; N is per science channel. B=N+2M+2V is the per-channel incident cost. Multiply by 24 for total cost. The plotted validation fraction is 2V/(B-2M), excluding fixed reference cost. There are 3,480,000 count slots, including 96,000 structural zeros from V=0 cases.

Reference-to-validation drift is 0.006722145420158068 in absolute probability. Main cases have zero validation-to-science residual; stress cases have 0.001, 0.003 or 0.005. Probability 0.003 means 0.3 percentage points. The guarded procedure assumes a fixed residual bound of 0.003; it does not learn that residual from the reference/validation counts.

Only `excluded` counts as rejection. Inconclusive computation is non-rejection. The domain remains the bounded generic barrier family from earlier studies. No astrophysical detection or stability claim follows from these counts.

```powershell
.\.venv\Scripts\python.exe scripts/run_validation_study.py --repetitions 500 --workers 4
.\.venv\Scripts\python.exe scripts/verify_validation_study.py
.\.venv\Scripts\python.exe -m unittest discover -s tests
```

Install dependencies from `requirements-lock.txt` for the recorded package versions. The design records Python and numerical-library versions. PNG/SVG figures are in `figures/study-008`; they can be recreated from the saved CSV/count records through the runner's `plots` function without repeating simulations.
