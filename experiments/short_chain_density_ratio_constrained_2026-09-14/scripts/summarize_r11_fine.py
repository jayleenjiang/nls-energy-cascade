#!/usr/bin/env python3
"""Create auditable R11 diagnostic, coefficient, and density-shift tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import tensorflow as tf

from double_calibrated_transport_ratio_model import load_bundle as load_double
from triple_calibrated_transport_ratio_model import load_bundle as load_triple
from v4_common import evaluation_roles as v4_roles, load_rows as v4_load
from v6_common import evaluation_roles as v6_roles, load_rows as v6_load

ROOT = Path(__file__).resolve().parent.parent
SEEDS = (8201, 8202, 8203)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def deterministic_indices(length: int, count: int, salt: str) -> np.ndarray:
    seed = int(hashlib.sha256(salt.encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed).choice(length, min(length, count), replace=False)


def values(model, state: np.ndarray, batch: int = 8192) -> np.ndarray:
    return np.concatenate([
        model.log_ratio(tf.constant(state[start:start + batch], tf.float64)).numpy()
        for start in range(0, len(state), batch)
    ])


def main() -> None:
    output = ROOT / "recovery_r11_fine" / "analysis"
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    tf.keras.backend.set_floatx("float64")

    sets = (
        ("fine_two_moment", ROOT / "formal_v7_fine/driven_test"),
        ("fine_three_moment", ROOT / "recovery_r11_fine/driven_test"),
        ("coarse_three_moment", ROOT / "recovery_r10_coarse/driven_test"),
    )
    diagnostic_rows = []
    for protocol, directory in sets:
        agreement = read_csv(directory / "seed_agreement.csv")
        maximum_rms = max(float(row["centered_log_ratio_rms"]) for row in agreement)
        for row in read_csv(directory / "model_summary.csv"):
            diagnostic_rows.append({
                "protocol": protocol,
                "training_base_seed": row["training_base_seed"],
                "reported_model_seed": row["model_seed"],
                "ess_fraction": row["ess_fraction"],
                "fp_abs_median": row["fp_abs_median"],
                "fp_abs_p90": row["fp_abs_p90"],
                "max_abs_moment_z": row["max_abs_moment_z"],
                "max_marginal_tv": row["max_marginal_tv"],
                "maximum_pairwise_centered_log_ratio_rms": maximum_rms,
                "per_model_pass": row["per_model_pass"],
            })
    write_csv(output / "diagnostic_comparison.csv", diagnostic_rows)

    coefficient_rows = []
    for calibration, template in (
        ("fine_two_moment", ROOT / "recovery_r9/factor_0.625/models/seed_{seed}"),
        ("fine_three_moment", ROOT / "recovery_r11_fine/models/seed_{seed}"),
        ("coarse_three_moment", ROOT / "recovery_r10_coarse/models/seed_{seed}"),
    ):
        for seed in SEEDS:
            directory = Path(str(template).format(seed=seed))
            config = json.loads((directory / "model_config.json").read_text())
            alpha = np.asarray(config["calibration_alpha"], np.float64)
            scale = np.asarray(config["calibration_scale"], np.float64)
            unshrunk = np.asarray(config["unshrunk_calibration_alpha"], np.float64)
            for index, observable in enumerate(config["calibration_observables"]):
                coefficient_rows.append({
                    "calibration": calibration,
                    "training_base_seed": seed,
                    "observable": observable,
                    "unshrunk_standardized_alpha": unshrunk[index],
                    "shrink_factor": config["calibration_shrink_factor"],
                    "shrunk_standardized_alpha": alpha[index],
                    "observable_center": config["calibration_center"][index],
                    "observable_scale": scale[index],
                    "shrunk_physical_coefficient": alpha[index] / scale[index],
                })
    write_csv(output / "calibration_coefficients.csv", coefficient_rows)

    old_models = [load_double(ROOT / f"recovery_r9/factor_0.625/models/seed_{s}") for s in SEEDS]
    new_models = [load_triple(ROOT / f"recovery_r11_fine/models/seed_{s}") for s in SEEDS]
    old_logz = {
        int(row["training_base_seed"]): float(row["log_normalizer"])
        for row in read_csv(ROOT / "formal_v7_fine/driven_test/model_summary.csv")
    }
    new_logz = {
        int(row["training_base_seed"]): float(row["log_normalizer"])
        for row in read_csv(ROOT / "recovery_r11_fine/driven_test/model_summary.csv")
    }

    shift_rows = []
    for support in ("fine", "coarse"):
        if support == "fine":
            target_rows, _ = v6_roles("driven", "fine", "test")
            state, _ = v6_load(target_rows)
        else:
            target_rows, _ = v4_roles("driven", "coarse", "test")
            state, _ = v4_load(target_rows)
        state = state[deterministic_indices(len(state), 50_000, f"r11-density-shift|{support}")]
        old_values = np.asarray([values(model, state) for model in old_models])
        new_values = np.asarray([values(model, state) for model in new_models])
        raw_shift = new_values - old_values
        centered_shift = raw_shift - raw_shift.mean(axis=1, keepdims=True)
        for index, seed in enumerate(SEEDS):
            if support == "fine":
                normalized_shift = raw_shift[index] - (new_logz[seed] - old_logz[seed])
                normalized_mean = float(normalized_shift.mean())
                normalized_rms = float(np.sqrt(np.mean(normalized_shift ** 2)))
            else:
                normalized_mean = float("nan")
                normalized_rms = float("nan")
            absolute = np.abs(centered_shift[index])
            shift_rows.append({
                "support": support,
                "comparison": f"seed_{seed}",
                "sample_count": len(state),
                "centered_shift_mean": centered_shift[index].mean(),
                "centered_shift_rms": np.sqrt(np.mean(centered_shift[index] ** 2)),
                "centered_abs_shift_p90": np.quantile(absolute, 0.90),
                "centered_abs_shift_p99": np.quantile(absolute, 0.99),
                "centered_abs_shift_max": absolute.max(),
                "normalized_log_density_shift_mean": normalized_mean,
                "normalized_log_density_shift_rms": normalized_rms,
            })
        old_ensemble = old_values.mean(axis=0)
        new_ensemble = new_values.mean(axis=0)
        ensemble_shift = new_ensemble - old_ensemble
        ensemble_centered = ensemble_shift - ensemble_shift.mean()
        absolute = np.abs(ensemble_centered)
        if support == "fine":
            delta_logz = np.mean([new_logz[s] - old_logz[s] for s in SEEDS])
            normalized_shift = ensemble_shift - delta_logz
            normalized_mean = float(normalized_shift.mean())
            normalized_rms = float(np.sqrt(np.mean(normalized_shift ** 2)))
        else:
            normalized_mean = float("nan")
            normalized_rms = float("nan")
        shift_rows.append({
            "support": support,
            "comparison": "ensemble_mean",
            "sample_count": len(state),
            "centered_shift_mean": ensemble_centered.mean(),
            "centered_shift_rms": np.sqrt(np.mean(ensemble_centered ** 2)),
            "centered_abs_shift_p90": np.quantile(absolute, 0.90),
            "centered_abs_shift_p99": np.quantile(absolute, 0.99),
            "centered_abs_shift_max": absolute.max(),
            "normalized_log_density_shift_mean": normalized_mean,
            "normalized_log_density_shift_rms": normalized_rms,
        })
    write_csv(output / "density_shift.csv", shift_rows)

    verdict = {
        "fine_validation_pass": json.loads((ROOT / "recovery_r11_fine/driven_validation/verdict.json").read_text())["complete_pass"],
        "fine_test_pass": json.loads((ROOT / "recovery_r11_fine/driven_test/verdict.json").read_text())["complete_pass"],
        "equilibrium_validation_pass": json.loads((ROOT / "recovery_r11_fine/equilibrium_validation/verdict.json").read_text())["complete_pass"],
        "equilibrium_test_pass": json.loads((ROOT / "recovery_r11_fine/equilibrium_test/verdict.json").read_text())["complete_pass"],
        "three_moment_timestep_pass": json.loads((ROOT / "recovery_r11_fine/timestep_comparison/verdict.json").read_text())["complete_pass"],
    }
    verdict["complete_pass"] = bool(all(verdict.values()))
    (output / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")


if __name__ == "__main__":
    main()

