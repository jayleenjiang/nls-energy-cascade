# Autocorrelation sector follow-up

Analysis-only re-evaluation of the saved three-mode autocorrelations. No
simulator is called by these scripts.

- Frozen specification and input hashes: `PROTOCOL.md`
- Analysis: `analyze_sectors.py`
- Machine-readable results: `results/`
- Figure source: `make_sector_figure.py`
- LaTeX/PDF report: `sector_followup_report.tex` and
  `sector_followup_report.pdf`
- Claim summary: `FINAL_VERDICT.md`

Run from the parent experiment directory:

```bash
python3 sector_followup/analyze_sectors.py
python3 sector_followup/make_sector_figure.py
```
