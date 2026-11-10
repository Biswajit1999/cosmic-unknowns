# Formal report

**Author: Biswajit Jana**

`calibration-limited-identification.pdf` is generated from [the canonical Markdown manuscript](../../reports/007-calibration-limited-identification.md). `report-manifest.json` records the source/PDF hashes, page count and figure/equation counts. The PDF is a public working report, not a peer-reviewed publication or an external submission.

After generating studies 004 and 005, run:

```text
python scripts/render_report_equations.py
python scripts/build_scientific_report.py
```

The first command requires the numerical dependencies. The second additionally requires `reportlab`, `pypdf` and Pillow. The tested PDF runtime used reportlab 4.4.9 and pypdf 6.21.4; these came from the bundled document runtime. To build in a single local environment, install the numerical requirements plus `requirements-report.txt`. The builder discovers the DejaVu fonts supplied with matplotlib.

Equation PNGs and PDF review renders are intermediate files under ignored `tmp/pdfs/`. Scientific plot files and source data remain versioned. The final PDF is rendered with Poppler for visual review; metadata and extracted text are checked programmatically. Font and pagination details can vary across library versions; the manifest identifies the delivered build.
