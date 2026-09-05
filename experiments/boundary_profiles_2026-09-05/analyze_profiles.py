#!/usr/bin/env python3
"""Merge fixed-seed profile runs and report frozen diagnostics."""

import csv
import hashlib
import json
import math
from statistics import NormalDist
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
OUT = ROOT / "analysis"
CURATED = ROOT / "curated_results" / "profiles"
BC_EQUATIONS = {
    "BC1": ("b_I=2*gamma*(2*T-F)", "b_phi=gamma*P"),
    "BC2": ("b_I=2*gamma*(2*T-F)", "b_phi=0"),
    "BC3": ("b_I=2*gamma*(2*T-I^2)", "b_phi=0"),
    "BC3b": ("b_I=2*gamma*(2*T-I^2)", "b_phi=gamma*P"),
}


def fnum(value):
    if value is None or value == "":
        return math.nan
    return float(value)


def read_commented_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def read_summary(path):
    with path.open() as handle:
        return next(csv.DictReader(handle))


def pooled(parts):
    """Pool (n, mean, se) sufficient statistics across independent groups."""
    parts = [(n, m, se) for n, m, se in parts if n > 1 and math.isfinite(m) and math.isfinite(se)]
    if not parts:
        return math.nan, math.nan, 0
    total = sum(n for n, _, _ in parts)
    mean = sum(n * m for n, m, _ in parts) / total
    ss = 0.0
    for n, m, se in parts:
        variance = se * se * n
        ss += (n - 1) * variance + n * (m - mean) ** 2
    variance = ss / (total - 1) if total > 1 else math.nan
    return mean, math.sqrt(max(variance, 0.0) / total), total


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_matrix():
    with (ROOT / "RUN_MATRIX.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    groups = defaultdict(list)
    for row in rows:
        key = (row["bc_label"], row["temp_label"], int(row["n"]))
        groups[key].append(row)
    return groups


def merge_table(paths, count_by_rep, value_pairs, key_fields):
    tables = [read_commented_csv(path) for path in paths]
    by_key = []
    for table in tables:
        by_key.append({tuple(row[k] for k in key_fields): row for row in table})
    keys = list(by_key[0])
    merged = []
    for key in keys:
        row = dict(zip(key_fields, key))
        for mean_field, se_field in value_pairs:
            parts = []
            for rep, table in enumerate(by_key):
                source = table[key]
                parts.append((count_by_rep[rep], fnum(source[mean_field]), fnum(source[se_field])))
            mean, se, _ = pooled(parts)
            row[mean_field] = mean
            row[se_field] = se
        merged.append(row)
    return merged


def source_prefix(row):
    return RAW / f'{row["bc_label"]}_{row["temp_label"]}_n{row["n"]}_rep{row["replicate"]}'


def midpoint(profile, n, mean_field, se_field, is_bond=False):
    # Protocol averages the two central objects for an even-sized object list.
    # With no cross-site covariance saved, use the sharp conservative upper
    # bound (se1+se2)/2 for the SE of a two-object average.
    size = n - 1 if is_bond else n
    indices = [size // 2] if size % 2 else [size // 2 - 1, size // 2]
    vals = [fnum(profile[i][mean_field]) for i in indices]
    ses = [fnum(profile[i][se_field]) for i in indices]
    return float(np.mean(vals)), float(sum(ses) / len(ses)), ";".join(str(i + 1) for i in indices)


def relative_checkpoint_change(checkpoints, n, field, site_index):
    keyed = defaultdict(dict)
    for row in checkpoints:
        keyed[fnum(row["checkpoint_time"])][int(row["j"]) - 1] = row
    times = sorted(keyed)
    if len(times) < 2:
        return math.nan, math.nan
    previous = fnum(keyed[times[-2]][site_index][field])
    final = fnum(keyed[times[-1]][site_index][field])
    denom = max(abs(final), 1e-15)
    return (final - previous) / denom, final - previous


def linear_fit(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or not np.all(np.isfinite(y)):
        return math.nan, math.nan, math.nan
    design = np.column_stack([np.ones(len(x)), x])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan
    slope_se = math.nan
    if len(x) > 2:
        cov = np.linalg.inv(design.T @ design) * ss_res / (len(x) - 2)
        slope_se = math.sqrt(max(cov[1, 1], 0.0))
    return float(coef[1]), slope_se, r2


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CURATED.mkdir(parents=True, exist_ok=True)
    groups = load_matrix()
    run_rows = []
    equilibrium_rows = []
    all_profiles = {}
    missing = []

    for key, reps in sorted(groups.items()):
        bc, temp, n = key
        reps = sorted(reps, key=lambda r: int(r["replicate"]))
        prefixes = [source_prefix(row) for row in reps]
        required = [Path(str(prefix) + suffix) for prefix in prefixes for suffix in
                    ("_profile.csv", "_burnin_checkpoints.csv", "_checkpoints.csv",
                     "_trajectory_diagnostics.csv", "_summary.csv")]
        absent = [str(path) for path in required if not path.exists()]
        if absent:
            missing.extend(absent)
            continue
        summaries = [read_summary(Path(str(prefix) + "_summary.csv")) for prefix in prefixes]
        counts = [int(s["valid_trajectories"]) for s in summaries]
        profiles = merge_table(
            [Path(str(p) + "_profile.csv") for p in prefixes], counts,
            [("mean_I", "se_mean_I"), ("mean_sin_theta", "se_mean_sin_theta")], ["j"])
        burns = merge_table(
            [Path(str(p) + "_burnin_checkpoints.csv") for p in prefixes], counts,
            [("mean_I", "se_mean_I"), ("mean_sin_theta", "se_mean_sin_theta")],
            ["checkpoint_time", "j"])
        checks = merge_table(
            [Path(str(p) + "_checkpoints.csv") for p in prefixes], counts,
            [("mean_I", "se_mean_I"), ("mean_sin_theta", "se_mean_sin_theta")],
            ["checkpoint_time", "j"])
        all_profiles[key] = profiles

        seeds = [int(s["seed"]) for s in summaries]
        valid = sum(counts)
        discarded = sum(int(s["discarded_trajectories"]) for s in summaries)
        nonfinite = sum(int(s["nonfinite_trajectories"]) for s in summaries)
        projections = sum(int(s["projection_count"]) for s in summaries)
        prefix_out = CURATED / f"{bc}_{temp}_n{n}"
        action_eq, phase_eq = BC_EQUATIONS[bc]
        with Path(str(prefix_out) + "_profile.csv").open("w", newline="") as handle:
            handle.write(f"# bc_label={bc}\n# {action_eq}\n# {phase_eq}\n")
            handle.write("# sigma_I=2*sqrt(2*gamma*T*I)\n# sigma_phi=sqrt(2*gamma*T/I)\n")
            handle.write(f"# n={n} T1={reps[0]['T1']} Tn={reps[0]['Tn']} gamma=0.1 dt_max=0.0005 ")
            handle.write(f"burnin={reps[0]['burnin']} measure={reps[0]['measure']} ")
            handle.write(f"seeds={';'.join(map(str,seeds))} samples={valid} ")
            handle.write(f"nonfinite_trajectories={nonfinite} discarded_trajectories={discarded} ")
            handle.write(f"projection_count={projections}\n")
            writer = csv.DictWriter(handle, fieldnames=["j", "mean_I", "se_mean_I",
                                                        "mean_sin_theta", "se_mean_sin_theta"])
            writer.writeheader()
            writer.writerows(profiles)
        write_csv(Path(str(prefix_out) + "_burnin_checkpoints.csv"), burns,
                  ["checkpoint_time", "j", "mean_I", "se_mean_I",
                   "mean_sin_theta", "se_mean_sin_theta"])
        write_csv(Path(str(prefix_out) + "_measurement_checkpoints.csv"), checks,
                  ["checkpoint_time", "j", "mean_I", "se_mean_I",
                   "mean_sin_theta", "se_mean_sin_theta"])

        # Replicate discrepancy in units of independent-replicate SE.
        rep_profiles = [read_commented_csv(Path(str(p) + "_profile.csv")) for p in prefixes]
        z_values = []
        for j in range(n):
            denom = math.hypot(fnum(rep_profiles[0][j]["se_mean_I"]),
                               fnum(rep_profiles[1][j]["se_mean_I"]))
            if denom > 0:
                z_values.append(abs(fnum(rep_profiles[0][j]["mean_I"]) -
                                    fnum(rep_profiles[1][j]["mean_I"])) / denom)
        for j in range(n - 1):
            denom = math.hypot(fnum(rep_profiles[0][j]["se_mean_sin_theta"]),
                               fnum(rep_profiles[1][j]["se_mean_sin_theta"]))
            if denom > 0:
                z_values.append(abs(fnum(rep_profiles[0][j]["mean_sin_theta"]) -
                                    fnum(rep_profiles[1][j]["mean_sin_theta"])) / denom)

        mid_I, mid_I_se, mid_sites = midpoint(profiles, n, "mean_I", "se_mean_I")
        mid_sin, mid_sin_se, mid_bonds = midpoint(
            profiles, n, "mean_sin_theta", "se_mean_sin_theta", True)
        left_change, _ = relative_checkpoint_change(checks, n, "mean_I", 0)
        right_change, _ = relative_checkpoint_change(checks, n, "mean_I", n - 1)
        max_left = max(fnum(r["max_I_left"]) for p in prefixes
                       for r in read_commented_csv(Path(str(p) + "_trajectory_diagnostics.csv"))
                       if r["valid"] == "1")
        max_right = max(fnum(r["max_I_right"]) for p in prefixes
                        for r in read_commented_csv(Path(str(p) + "_trajectory_diagnostics.csv"))
                        if r["valid"] == "1")
        run_rows.append({
            "bc_label": bc, "temp_label": temp, "T1": reps[0]["T1"], "Tn": reps[0]["Tn"],
            "n": n, "burnin": reps[0]["burnin"], "measure": reps[0]["measure"],
            "seeds": ";".join(map(str, seeds)), "valid_trajectories": valid,
            "nonfinite_trajectories": nonfinite, "discarded_trajectories": discarded,
            "projection_count": projections, "mean_I_left": profiles[0]["mean_I"],
            "se_I_left": profiles[0]["se_mean_I"], "mean_I_right": profiles[-1]["mean_I"],
            "se_I_right": profiles[-1]["se_mean_I"], "mean_I_mid": mid_I,
            "se_I_mid_conservative": mid_I_se, "mid_sites": mid_sites,
            "mean_sin_mid": mid_sin, "se_sin_mid_conservative": mid_sin_se,
            "mid_bonds": mid_bonds, "max_seed_z": max(z_values) if z_values else math.nan,
            "median_seed_z": float(np.median(z_values)) if z_values else math.nan,
            "left_last_quarter_relative_change": left_change,
            "right_last_quarter_relative_change": right_change,
            "max_I_left_over_trajectories": max_left,
            "max_I_right_over_trajectories": max_right,
            "unstable_nonfinite": nonfinite > 0,
        })

        if temp == "T6_T6":
            means = np.array([fnum(r["mean_I"]) for r in profiles])
            ses = np.array([fnum(r["se_mean_I"]) for r in profiles])
            weighted_mean = float(np.sum(means / ses**2) / np.sum(1.0 / ses**2))
            action_z = np.abs((means - weighted_mean) / ses)
            sine_z = np.array([abs(fnum(r["mean_sin_theta"])) /
                               fnum(r["se_mean_sin_theta"]) for r in profiles[:-1]])
            bonf = NormalDist().inv_cdf(1.0 - 0.05 / (2.0 * len(sine_z)))
            equilibrium_rows.append({
                "bc_label": bc, "n": n, "valid_trajectories": valid,
                "mean_I_spatial": float(np.mean(means)), "I_range": float(np.ptp(means)),
                "relative_I_range": float(np.ptp(means) / np.mean(means)),
                "max_action_z_from_weighted_constant": float(np.max(action_z)),
                "sine_pointwise_95_failures": int(np.sum(sine_z > 1.96)),
                "sine_max_abs_z": float(np.max(sine_z)),
                "sine_bonferroni_threshold_approx": bonf,
                "sine_simultaneous_pass": bool(np.max(sine_z) <= bonf),
            })

    if missing:
        (OUT / "missing_outputs.txt").write_text("\n".join(missing) + "\n")
        print(f"incomplete: {len(missing)} required files missing")
        raise SystemExit(3)

    run_fields = list(run_rows[0])
    write_csv(OUT / "logical_run_summary.csv", run_rows, run_fields)
    write_csv(OUT / "equilibrium_checks.csv", equilibrium_rows, list(equilibrium_rows[0]))

    scaling_rows = []
    grouped = defaultdict(list)
    for row in run_rows:
        grouped[(row["bc_label"], row["temp_label"])].append(row)
    for (bc, temp), rows in sorted(grouped.items()):
        rows.sort(key=lambda r: int(r["n"]))
        ns = [int(r["n"]) for r in rows]
        mids = [fnum(r["mean_I_mid"]) for r in rows]
        sins = [fnum(r["mean_sin_mid"]) for r in rows]
        slope_I, slope_I_se, r2_I = linear_fit(np.log(ns), np.log(mids))
        # Sine can change sign; report n-dependence only when one nonzero sign is shared.
        if all(s != 0 for s in sins) and len({math.copysign(1, s) for s in sins}) == 1:
            slope_s, slope_s_se, r2_s = linear_fit(np.log(ns), np.log(np.abs(sins)))
        else:
            slope_s = slope_s_se = r2_s = math.nan
        scaling_rows.append({
            "bc_label": bc, "temp_label": temp, "n_values": ";".join(map(str, ns)),
            "mid_I_values": ";".join(f"{x:.17g}" for x in mids),
            "loglog_slope_mid_I": slope_I, "slope_se_three_point": slope_I_se,
            "R2_mid_I": r2_I, "reference_slope": -0.5 if bc in ("BC1", "BC2") else 0.0,
            "mid_sin_values": ";".join(f"{x:.17g}" for x in sins),
            "loglog_slope_abs_mid_sin": slope_s, "sine_slope_se": slope_s_se,
            "R2_abs_mid_sin": r2_s,
        })
    write_csv(OUT / "midpoint_scaling.csv", scaling_rows, list(scaling_rows[0]))

    manifest = []
    for path in sorted([p for p in ROOT.rglob("*") if p.is_file() and
                        not any(part == ".git" for part in p.parts)]):
        if path.name == "FILE_HASHES.csv":
            continue
        manifest.append({"sha256": sha256(path), "bytes": path.stat().st_size,
                         "path": str(path.relative_to(ROOT))})
    write_csv(OUT / "FILE_HASHES.csv", manifest, ["sha256", "bytes", "path"])
    print(f"complete: {len(groups)} logical runs, {sum(int(r['valid_trajectories']) for r in run_rows)} valid trajectories")


if __name__ == "__main__":
    main()
