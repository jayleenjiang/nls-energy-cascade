#!/usr/bin/env python3
"""Export source-trace metrics used by the manuscript integrity audit.

This script is intentionally conservative:

* LTE table entries are recomputed from archived histogram/profile files.
* Short-chain neural-network numbers are extracted from archived notebook
  outputs and paired with source/image hashes.  The script does not attempt to
  rerun TensorFlow training or inference.
* Notebook cells that contain diagnostic code but no archived output are
  reported explicitly as unverified, so the manuscript can avoid relying on
  them as quantitative evidence.
"""

from __future__ import annotations

import hashlib
import json
import re
import struct
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import generate_manuscript_figures as manuscript_figures


ROOT = Path(__file__).resolve().parents[3]
REVISION = Path(__file__).resolve().parents[1]
OUT_PATH = REVISION / "source_trace_metrics.json"


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": rel(path),
        "exists": path.exists(),
    }
    if path.exists():
        record["bytes"] = path.stat().st_size
        record["sha256"] = sha256_file(path)
    return record


def png_dimensions(path: Path) -> tuple[int, int] | None:
    if not path.exists() or path.suffix.lower() != ".png":
        return None
    with path.open("rb") as f:
        header = f.read(24)
    if len(header) >= 24 and header[:8] == b"\x89PNG\r\n\x1a\n":
        return struct.unpack(">II", header[16:24])
    return None


def image_record(path: Path) -> dict[str, Any]:
    record = file_record(path)
    dims = png_dimensions(path)
    if dims:
        record["width_px"], record["height_px"] = dims
    return record


def read_profile(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            values[int(parts[0])] = float(parts[1])
    if not values:
        raise ValueError(f"No profile values in {path}")
    return values


def pair_kinetic_temperature(ness_profile: Path, eq_profile: Path, j: int) -> dict[str, float]:
    ness = read_profile(ness_profile)
    eq = read_profile(eq_profile)
    numerator = 0.5 * (ness[j] + ness[j + 1])
    denominator = 0.5 * (eq[j] + eq[j + 1])
    return {
        "pair_mean_ness_action": numerator,
        "pair_mean_eq_action": denominator,
        "T_kin": 6.0 * numerator / denominator,
    }


def rounded_row(fit: dict[str, Any], tkin: float) -> dict[str, Any]:
    x = float(fit["slope"])
    r2 = float(fit["weighted_r2"])
    return {
        "x_4dp": f"{x:.4f}",
        "weighted_r2_4dp": f"{r2:.4f}",
        "six_over_x_2dp": f"{6.0 / x:.2f}",
        "T_kin_2dp": f"{tkin:.2f}",
    }


def export_lte_table() -> dict[str, Any]:
    cases = [
        {
            "n": 25,
            "j": 6,
            "a_over_n": 0.24,
            "q": ROOT / "lte/n25 data/n25_j6.hist",
            "p": ROOT / "lte/n25 data/n25_eq_j6.hist",
            "profile": ROOT / "experiments/lte/test_profile.txt",
            "eq_profile": ROOT / "experiments/lte/eigen/eigen_eq/eq_T6_profile.txt",
        },
        {
            "n": 25,
            "j": 12,
            "a_over_n": 0.48,
            "q": ROOT / "lte/n25 data/n25_j12.hist",
            "p": ROOT / "lte/n25 data/n25_eq_j12.hist",
            "profile": ROOT / "experiments/lte/test_profile.txt",
            "eq_profile": ROOT / "experiments/lte/eigen/eigen_eq/eq_T6_profile.txt",
        },
        {
            "n": 25,
            "j": 18,
            "a_over_n": 0.72,
            "q": ROOT / "lte/n25 data/n25_j18.hist",
            "p": ROOT / "lte/n25 data/n25_eq_j18.hist",
            "profile": ROOT / "experiments/lte/test_profile.txt",
            "eq_profile": ROOT / "experiments/lte/eigen/eigen_eq/eq_T6_profile.txt",
        },
        {
            "n": 50,
            "j": 12,
            "a_over_n": 0.24,
            "q": ROOT / "lte/n50 data/simd_n50_j12.hist",
            "p": ROOT / "lte/n50 data/simd_n50_eq_j12.hist",
            "profile": ROOT / "lte/n50 data/simd_n50_profile.txt",
            "eq_profile": ROOT / "lte/n50 data/simd_n50_eq_profile.txt",
        },
        {
            "n": 50,
            "j": 24,
            "a_over_n": 0.48,
            "q": ROOT / "lte/n50 data/simd_n50_j24.hist",
            "p": ROOT / "lte/n50 data/simd_n50_eq_j24.hist",
            "profile": ROOT / "lte/n50 data/simd_n50_profile.txt",
            "eq_profile": ROOT / "lte/n50 data/simd_n50_eq_profile.txt",
        },
        {
            "n": 50,
            "j": 36,
            "a_over_n": 0.72,
            "q": ROOT / "lte/n50 data/simd_n50_j36.hist",
            "p": ROOT / "lte/n50 data/simd_n50_eq_j36.hist",
            "profile": ROOT / "lte/n50 data/simd_n50_profile.txt",
            "eq_profile": ROOT / "lte/n50 data/simd_n50_eq_profile.txt",
        },
        {
            "n": 100,
            "j": 24,
            "a_over_n": 0.24,
            "q": ROOT / "lte/n100 data/n100_j24.hist",
            "p": ROOT / "lte/n100 data/n100_eq_j24.hist",
            "profile": ROOT / "lte/n100 data/n100_profile.txt",
            "eq_profile": ROOT / "lte/n100 data/n100_eq_profile.txt",
        },
        {
            "n": 100,
            "j": 48,
            "a_over_n": 0.48,
            "q": ROOT / "lte/n100 data/n100_j48.hist",
            "p": ROOT / "lte/n100 data/n100_eq_j48.hist",
            "profile": ROOT / "lte/n100 data/n100_profile.txt",
            "eq_profile": ROOT / "lte/n100 data/n100_eq_profile.txt",
        },
        {
            "n": 100,
            "j": 72,
            "a_over_n": 0.72,
            "q": ROOT / "lte/n100 data/n100_j72.hist",
            "p": ROOT / "lte/n100 data/n100_eq_j72.hist",
            "profile": ROOT / "lte/n100 data/n100_profile.txt",
            "eq_profile": ROOT / "lte/n100 data/n100_eq_profile.txt",
        },
    ]

    rows: list[dict[str, Any]] = []
    for case in cases:
        q_hist = manuscript_figures.load_hist(case["q"])
        p_hist = manuscript_figures.load_hist(case["p"])
        fit = manuscript_figures.fit_rescaled_equilibrium(q_hist, p_hist, min_count=1, i_max=2.5)
        tkin = pair_kinetic_temperature(case["profile"], case["eq_profile"], int(case["j"]))
        row = {
            "n": case["n"],
            "pair_index_j": case["j"],
            "a_over_n": case["a_over_n"],
            "q_hist": file_record(case["q"]),
            "p_hist": file_record(case["p"]),
            "profile": file_record(case["profile"]),
            "eq_profile": file_record(case["eq_profile"]),
            "fit": {
                "slope_x": float(fit["slope"]),
                "intercept": float(fit["intercept"]),
                "weighted_r2": float(fit["weighted_r2"]),
                "six_over_x": 6.0 / float(fit["slope"]),
                "bins_used": int(fit["bins_used"]),
                "mask": "q_count>=1 and p_count>=1 and I_a,I_b<2.5",
                "weights": "q_count",
            },
            "kinetic_temperature": tkin,
        }
        row["manuscript_rounding"] = rounded_row(fit, float(tkin["T_kin"]))
        rows.append(row)

    control_cases = []
    for j in [6, 12, 18]:
        q = ROOT / f"lte/T712:8:9 data/simd_eq_T712_j{j}.hist"
        p = ROOT / f"lte/n25 data/n25_eq_j{j}.hist"
        fit = manuscript_figures.fit_rescaled_equilibrium(
            manuscript_figures.load_hist(q),
            manuscript_figures.load_hist(p),
            min_count=1,
            i_max=2.5,
        )
        control_cases.append(
            {
                "label": "p_7.12_vs_p_6",
                "j": j,
                "q_hist": file_record(q),
                "p_hist": file_record(p),
                "slope_x": float(fit["slope"]),
                "weighted_r2": float(fit["weighted_r2"]),
                "six_over_x": 6.0 / float(fit["slope"]),
                "bins_used": int(fit["bins_used"]),
            }
        )

    q = ROOT / "lte/n25 data/n25_j12.hist"
    p = ROOT / "lte/T712:8:9 data/simd_eq_T712_j12.hist"
    fit = manuscript_figures.fit_rescaled_equilibrium(
        manuscript_figures.load_hist(q),
        manuscript_figures.load_hist(p),
        min_count=1,
        i_max=2.5,
    )
    control_cases.append(
        {
            "label": "NESS_vs_p_7.12_midchain",
            "j": 12,
            "q_hist": file_record(q),
            "p_hist": file_record(p),
            "slope_x": float(fit["slope"]),
            "weighted_r2": float(fit["weighted_r2"]),
            "six_over_x": 6.0 / float(fit["slope"]),
            "bins_used": int(fit["bins_used"]),
        }
    )

    p712 = [c for c in control_cases if c["label"] == "p_7.12_vs_p_6"]
    return {
        "description": "Source-traced recomputation of Table tab:lte from archived histograms and profiles.",
        "algorithm": {
            "fit": "weighted least squares of log q = c + x log p over core bins",
            "mask": "q_count>=1 and p_count>=1 and I_a,I_b<2.5",
            "weights": "q_count",
            "T_kin": "6 times pair-averaged NESS action divided by pair-averaged equilibrium action",
        },
        "rows": rows,
        "controls": {
            "sitewise": control_cases,
            "p_7.12_vs_p_6_mean_slope": sum(c["slope_x"] for c in p712) / len(p712),
            "p_7.12_vs_p_6_mean_weighted_r2": sum(c["weighted_r2"] for c in p712) / len(p712),
        },
    }


def load_notebook(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def output_text(cell: dict[str, Any]) -> str:
    parts: list[str] = []
    for out in cell.get("outputs", []):
        if "text" in out:
            text = out["text"]
            parts.append("".join(text) if isinstance(text, list) else str(text))
        data = out.get("data", {})
        if "text/plain" in data:
            text = data["text/plain"]
            parts.append("".join(text) if isinstance(text, list) else str(text))
    return "\n".join(parts)


def source_text(cell: dict[str, Any]) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else str(source)


def regex_float(pattern: str, text: str) -> float | None:
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def regex_int(pattern: str, text: str) -> int | None:
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


def export_short_chain() -> dict[str, Any]:
    fp_nb_path = ROOT / "KDE/4:15_NN/FKE_5d_NLS.ipynb"
    eigen_nb_path = ROOT / "KDE/FKE_eigen.ipynb"
    fp_nb = load_notebook(fp_nb_path)
    eigen_nb = load_notebook(eigen_nb_path)

    fp_cells = fp_nb["cells"]
    eigen_cells = eigen_nb["cells"]

    load_text = output_text(fp_cells[2])
    eval_text = output_text(fp_cells[8])
    residual_text = output_text(fp_cells[9])
    symmetry_text = output_text(fp_cells[11])
    eq_model_text = output_text(fp_cells[19])
    eq_validation_text = output_text(fp_cells[20])

    eq_validation = []
    for match in re.finditer(
        r"I=([0-9.]+): mean err=([0-9.]+)%, log ratio=\[([-0-9.]+),([-0-9.]+)\]",
        eq_validation_text,
    ):
        eq_validation.append(
            {
                "I": float(match.group(1)),
                "mean_relative_error_percent": float(match.group(2)),
                "log_ratio_min": float(match.group(3)),
                "log_ratio_max": float(match.group(4)),
            }
        )

    eigen_load_text = output_text(eigen_cells[1])
    eigen_eval_text = output_text(eigen_cells[7])
    eigen_residual_text = output_text(eigen_cells[8])

    unexecuted = []
    for name, idx in [
        ("phase_locking_peak_table", 14),
        ("angular_width_ratio", 15),
        ("middle_mode_current_balance", 16),
    ]:
        cell = fp_cells[idx]
        unexecuted.append(
            {
                "diagnostic": name,
                "notebook": rel(fp_nb_path),
                "cell_index": idx,
                "archived_output_present": bool(cell.get("outputs")),
                "source_sha256": hashlib.sha256(source_text(cell).encode("utf-8")).hexdigest(),
                "status": "unverified_no_archived_output" if not cell.get("outputs") else "archived_output_present",
            }
        )

    images = {
        "eq_validation": image_record(REVISION / "eq_validation.png"),
        "neq_density": image_record(REVISION / "neq_density.png"),
        "symmetry_breaking": image_record(REVISION / "symmetry_breaking.png"),
        "Q1_slices": image_record(REVISION / "Q1_slices.png"),
    }

    return {
        "description": "Archived notebook-output metrics and figure hashes for the short-chain Fokker--Planck and eigenfunction sections.",
        "tensorflow_rerun": {
            "attempted": False,
            "reason": "Current local Python environments do not provide TensorFlow; this export records archived notebook outputs instead of rerunning models.",
        },
        "source_files": {
            "fokker_planck_notebook": file_record(fp_nb_path),
            "fokker_planck_boxes": file_record(ROOT / "KDE/4:15_NN/NLS_FP_boxes.txt"),
            "fokker_planck_density": file_record(ROOT / "KDE/4:15_NN/NLS_FP_density.txt"),
            "equilibrium_boxes": file_record(ROOT / "KDE/4:15_NN/eq/NLS_FP_boxes.txt"),
            "equilibrium_density": file_record(ROOT / "KDE/4:15_NN/eq/NLS_FP_density.txt"),
            "nonequilibrium_model": file_record(ROOT / "KDE/4:15_NN/h5_files/final.keras"),
            "equilibrium_model": file_record(ROOT / "KDE/4:15_NN/h5_files_eq/final.keras"),
            "eigen_notebook": file_record(eigen_nb_path),
            "eigen_X": file_record(ROOT / "KDE/backward_NLS_X.txt"),
            "eigen_Q1": file_record(ROOT / "KDE/backward_NLS_Q1.txt"),
            "eigen_model": file_record(ROOT / "KDE/h5_files_eigen/final.keras"),
        },
        "figures": images,
        "fokker_planck_training_archive": {
            "data_loss_points": regex_int(r"Data loss points: ([0-9]+)", load_text),
            "zero_density_boxes": regex_int(r"Zero-density boxes: ([0-9]+)", load_text),
            "log_density_mean": regex_float(r"Log density: mean=([-0-9.]+)", load_text),
            "log_density_std": regex_float(r"std=([-0-9.]+)", load_text),
            "pde_collocation_points": regex_int(r"PDE collocation: ([0-9]+)", load_text),
            "rmse_log_scaled": regex_float(r"RMSE \(log-scaled\): ([0-9.eE+-]+)", eval_text),
            "rmse_density": regex_float(r"RMSE \(density\):\s+([0-9.eE+-]+)", eval_text),
            "fp_residual_mean_abs_over_p": regex_float(r"mean \|residual/p\|: ([0-9.eE+-]+)", residual_text),
            "fp_residual_median": regex_float(r"median:\s+([0-9.eE+-]+)", residual_text),
            "fp_residual_90th_percentile": regex_float(r"90th percentile:\s+([0-9.eE+-]+)", residual_text),
        },
        "equilibrium_validation": {
            "notebook": rel(fp_nb_path),
            "cell_index": 20,
            "eq_model_loaded": {
                "mu": regex_float(r"mu=([-0-9.]+)", eq_model_text),
                "sigma": regex_float(r"sig=([-0-9.]+)", eq_model_text),
            },
            "rows": eq_validation,
        },
        "symmetry_breaking_archive": {
            "notebook": rel(fp_nb_path),
            "cell_index": 11,
            "max_asymmetry_percent_unmasked": regex_float(r"Max asymmetry: ([0-9.]+)%", symmetry_text),
            "mean_asymmetry_percent_unmasked": regex_float(r"Mean asymmetry: ([0-9.]+)%", symmetry_text),
            "manuscript_use": "qualitative_only; unmasked extrema are sensitive to low-density regions",
        },
        "eigen_surrogate_archive": {
            "notebook": rel(eigen_nb_path),
            "data_points": regex_int(r"Data points: ([0-9]+)", eigen_load_text),
            "Q1_min": regex_float(r"Q1 range: \[([-0-9.]+),", eigen_load_text),
            "Q1_max": regex_float(r"Q1 range: \[-?[0-9.]+, ([0-9.]+)\]", eigen_load_text),
            "Q1_mean": regex_float(r"Q1 mean=([-0-9.]+)", eigen_load_text),
            "Q1_std": regex_float(r"std=([-0-9.]+)", eigen_load_text),
            "collocation_points": regex_int(r"Collocation points: ([0-9]+)", eigen_load_text),
            "rmse": regex_float(r"RMSE = ([0-9.]+)", eigen_eval_text),
            "relative_rmse": regex_float(r"Relative RMSE = ([0-9.]+)", eigen_eval_text),
            "normalized_pde_residual_median": regex_float(r"median = ([0-9.]+)", eigen_residual_text),
            "normalized_pde_residual_mean": regex_float(r"mean\s+= ([0-9.]+)", eigen_residual_text),
            "normalized_pde_residual_90th": regex_float(r"90th\s+= ([0-9.]+)", eigen_residual_text),
        },
        "unverified_notebook_diagnostics": unexecuted,
    }


def main() -> None:
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "lte_table": export_lte_table(),
        "short_chain": export_short_chain(),
    }
    OUT_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"wrote": rel(OUT_PATH), "bytes": OUT_PATH.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
