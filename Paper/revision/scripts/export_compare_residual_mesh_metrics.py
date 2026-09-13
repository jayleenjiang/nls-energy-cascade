#!/usr/bin/env python3
"""Export quantitative metrics for the MATLAB-style LTE residual mesh figure.

The manuscript uses ``report_assets/compare_residual_mesh.pdf`` as a
structural diagnostic generated from the convention in
``/Users/jayleenjiang/Documents/MATLAB/lte/compare_residual.m``.  This helper
recomputes the numerical quantities behind that exact plotting convention so
the figure is accompanied by source-traced residual magnitudes rather than only
visual evidence.

The convention intentionally differs from the fixed weighted-core estimator in
``source_trace_metrics.json``:

* the fit is the MATLAB ``polyfit(log(pn(mask)), log(qn(mask)), 1)`` with
  ``mask=(Q>50)&(P>50)`` over the full histogram domain;
* the displayed surface masks ``Q<50`` or ``P==0``;
* the slice index is ``theta=20`` in MATLAB's one-based indexing.

The resulting metrics are descriptive diagnostics for the mesh figure, not a
replacement for the main LTE slope and even/odd residual tables.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
OUT_JSON = REVISION / "report_assets/compare_residual_mesh_metrics.json"
OUT_MD = REVISION / "report_assets/compare_residual_mesh_metrics.md"
OUT_TXT = REVISION / "report_assets/compare_residual_mesh_metrics.txt"

NB_EXPECTED = 80
THETA_MATLAB_INDEX = 20


@dataclass(frozen=True)
class MeshCase:
    label: str
    n: int
    j: int
    q_hist: str
    p_hist: str


CASES = [
    MeshCase(
        label="n15",
        n=15,
        j=4,
        q_hist="lte/n15 data/n15_j4.hist",
        p_hist="lte/n15 data/n15_eq_j4.hist",
    ),
    MeshCase(
        label="n25",
        n=25,
        j=6,
        q_hist="lte/n25 data/n25_j6.hist",
        p_hist="lte/n25 data/n25_eq_j6.hist",
    ),
    MeshCase(
        label="n50",
        n=50,
        j=12,
        q_hist="lte/n50 data/simd_n50_j12.hist",
        p_hist="lte/n50 data/simd_n50_eq_j12.hist",
    ),
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    return {
        "path": rel(path),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else None,
        "sha256": sha256_file(path) if path.exists() else None,
    }


def load_dense_counts(path: Path) -> tuple[dict[str, Any], np.ndarray]:
    """Load a dense ``ia ib it count`` histogram as a MATLAB-compatible vector."""

    meta: dict[str, Any] = {}
    data_start = 0
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if parts[0] == "NB":
            meta["NB"] = int(parts[1])
        elif parts[0] in {"I_LO", "I_HI", "TH_LO", "TH_HI"}:
            meta[parts[0]] = float(parts[1])
        elif parts[0] in {"TOTAL", "OVERFLOW"}:
            meta[parts[0]] = int(parts[1])
        elif len(parts) == 4:
            data_start = i
            break
    nb = int(meta["NB"])
    counts = np.zeros(nb * nb * nb, dtype=np.float64)
    cursor = 0
    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) != 4:
            continue
        counts[cursor] = float(parts[3])
        cursor += 1
    if cursor != nb * nb * nb:
        raise ValueError(f"{path}: expected {nb**3} dense rows, found {cursor}")
    return meta, counts


def slice_stats(values: np.ndarray) -> dict[str, Any]:
    if values.size == 0:
        raise ValueError("Cannot summarize an empty residual slice")
    abs_values = np.abs(values)
    return {
        "bins": int(values.size),
        "rms": float(np.sqrt(np.mean(values * values))),
        "mean_abs": float(np.mean(abs_values)),
        "median": float(np.median(values)),
        "q95_abs": float(np.quantile(abs_values, 0.95)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
    }


def compute_case(case: MeshCase) -> dict[str, Any]:
    q_path = ROOT / case.q_hist
    p_path = ROOT / case.p_hist
    q_meta, q_counts = load_dense_counts(q_path)
    p_meta, p_counts = load_dense_counts(p_path)
    if q_meta["NB"] != p_meta["NB"]:
        raise ValueError(f"{case.label}: mismatched NB values")
    nb = int(q_meta["NB"])
    if nb != NB_EXPECTED:
        raise ValueError(f"{case.label}: expected NB={NB_EXPECTED}, found {nb}")
    ntot = nb * nb * nb

    qn = q_counts / float(np.sum(q_counts))
    pn = p_counts / float(np.sum(p_counts))
    fit_mask = (q_counts > 50) & (p_counts > 50)
    slope, intercept = np.polyfit(np.log(pn[fit_mask]), np.log(qn[fit_mask]), 1)

    residual = np.full(ntot, np.nan, dtype=float)
    valid = (q_counts > 0) & (p_counts > 0)
    residual[valid] = np.log(qn[valid]) - (slope * np.log(pn[valid]) + intercept)
    residual[(q_counts < 50) | (p_counts == 0)] = np.nan

    theta_zero_based = THETA_MATLAB_INDEX - 1
    slice_indices = np.arange(theta_zero_based, ntot, nb, dtype=np.int64)
    slice_values = residual[slice_indices]
    finite = np.isfinite(slice_values)

    ia_centers = (np.floor(slice_indices / (nb * nb)) + 0.5) * 4.0 / nb
    ib_centers = (np.mod(np.floor(slice_indices / nb), nb) + 0.5) * 4.0 / nb
    core = finite & (ia_centers < 2.5) & (ib_centers < 2.5)

    theta_center = float(
        float(q_meta["TH_LO"])
        + (theta_zero_based + 0.5) * (float(q_meta["TH_HI"]) - float(q_meta["TH_LO"])) / nb
    )
    finite_values = slice_values[finite]
    core_values = slice_values[core]
    return {
        "label": case.label,
        "n": case.n,
        "pair_index_j": case.j,
        "a_over_n": case.j / case.n,
        "q_hist": file_record(q_path),
        "p_hist": file_record(p_path),
        "fit": {
            "slope_x": float(slope),
            "intercept": float(intercept),
            "fit_mask_bins": int(np.sum(fit_mask)),
            "fit_mask": "Q>50 and P>50 over the full 80^3 histogram",
            "fit_method": "unweighted polyfit(log(pn), log(qn), 1), matching compare_residual.m",
        },
        "slice": {
            "theta_matlab_index": THETA_MATLAB_INDEX,
            "theta_center": theta_center,
            "display_mask": "Q>=50 and P>0 on the theta slice, matching residual(Q<50 | P==0)=NaN",
            "display": slice_stats(finite_values),
            "core_IaIb_lt_2p5": slice_stats(core_values),
        },
        "manuscript_rounding": {
            "slope_3dp": f"{slope:.3f}",
            "display_rms_3dp": f"{slice_stats(finite_values)['rms']:.3f}",
            "core_rms_3dp": f"{slice_stats(core_values)['rms']:.3f}",
            "display_mean_abs_3dp": f"{slice_stats(finite_values)['mean_abs']:.3f}",
        },
    }


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# LTE residual mesh metrics",
        "",
        f"Generated: `{payload['generated_at_utc']}`",
        "",
        "These metrics reproduce the plotting convention of",
        "`/Users/jayleenjiang/Documents/MATLAB/lte/compare_residual.m` for",
        "`report_assets/compare_residual_mesh.pdf`.  They are descriptive",
        "mesh-slice diagnostics and are not replacements for the fixed",
        "weighted-core LTE estimators in `source_trace_metrics.json`.",
        "",
        "## Method",
        "",
        f"- Histogram grid: `NB={NB_EXPECTED}`.",
        f"- MATLAB theta index: `{THETA_MATLAB_INDEX}`.",
        "- Fit mask: `Q>50` and `P>50` over the full histogram.",
        "- Surface display mask: `Q>=50` and `P>0` on the displayed slice.",
        "- Reported norms are unweighted over displayed mesh bins.",
        "",
        "## Summary",
        "",
        "| case | pair | slope x | displayed bins | slice RMS | slice mean abs | core RMS | residual range |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in payload["rows"]:
        display = row["slice"]["display"]
        core = row["slice"]["core_IaIb_lt_2p5"]
        lines.append(
            f"| `{row['label']}` | {row['pair_index_j']} | "
            f"{row['fit']['slope_x']:.6f} | {display['bins']} | "
            f"{display['rms']:.6f} | {display['mean_abs']:.6f} | "
            f"{core['rms']:.6f} | [{display['min']:.6f}, {display['max']:.6f}] |"
        )
    path.write_text("\n".join(lines) + "\n")


def write_txt(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "Generated by exact Python reproduction of /Users/jayleenjiang/Documents/MATLAB/lte/compare_residual.m",
        "NB=80, theta=20, mask=(Q>50)&(P>50), display mask=(Q<50)|(P==0), fit=polyfit(log(pn),log(qn),1)",
    ]
    for row in payload["rows"]:
        display = row["slice"]["display"]
        core = row["slice"]["core_IaIb_lt_2p5"]
        lines.append(
            f"{row['label']} j={row['pair_index_j']} "
            f"slope={row['fit']['slope_x']:.8f} "
            f"mask_bins={row['fit']['fit_mask_bins']} "
            f"display_bins={display['bins']} "
            f"slice_rms={display['rms']:.8f} "
            f"slice_mean_abs={display['mean_abs']:.8f} "
            f"core_rms={core['rms']:.8f} "
            f"resid_range=[{display['min']:.8f},{display['max']:.8f}]"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    rows = [compute_case(case) for case in CASES]
    payload = {
        "generated_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "figure": file_record(REVISION / "report_assets/compare_residual_mesh.pdf"),
        "matlab_convention": {
            "source": "/Users/jayleenjiang/Documents/MATLAB/lte/compare_residual.m",
            "NB": NB_EXPECTED,
            "theta_matlab_index": THETA_MATLAB_INDEX,
            "fit": "polyfit(log(pn(mask)), log(qn(mask)), 1)",
            "fit_mask": "(Q > 50) & (P > 50)",
            "display_mask": "resid(Q < 50 | P == 0) = NaN",
        },
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_markdown(OUT_MD, payload)
    write_txt(OUT_TXT, payload)
    print(json.dumps({"wrote": [rel(OUT_JSON), rel(OUT_MD), rel(OUT_TXT)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
