#!/usr/bin/env python3
"""Apply the frozen Phase-III multi-k gates and build final products."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
PRIMARY = ROOT / "analysis" / "primary"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=1.0e-12)


def one(rows: list[dict[str, str]], description: str, **criteria: object) -> dict[str, str]:
    matches = []
    for row in rows:
        passed = True
        for key, expected in criteria.items():
            actual = row[key]
            if isinstance(expected, float):
                passed = passed and close(float(actual), expected)
            elif isinstance(expected, int):
                passed = passed and int(actual) == expected
            else:
                passed = passed and actual == expected
        if passed:
            matches.append(row)
    if len(matches) != 1:
        raise ValueError(f"expected one {description}, found {len(matches)}")
    return matches[0]


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    frozen = read_csv(ROOT / "FROZEN_GRID.csv")
    pairs = read_csv(PRIMARY / "paired_gc_residuals.csv")
    aggregates = read_csv(PRIMARY / "controlled_aggregate.csv")
    time_rows = read_csv(PRIMARY / "time_convergence_pairs.csv")
    population_members = read_csv(PRIMARY / "population_convergence_members.csv")
    population_pairs = read_csv(PRIMARY / "population_convergence_pairs.csv")
    phase2 = read_csv(
        ROOT.parent / "long_chain_ft_phase2_2026-08-29" / "FINAL_SUMMARY.csv"
    )
    control_members = read_csv(
        ROOT / "analysis" / "n40_outer_comparisons" / "control_members.csv"
    )
    control_pairs = read_csv(
        ROOT / "analysis" / "n40_outer_comparisons" / "control_pairs.csv"
    )

    result_rows: list[dict[str, object]] = []
    chain_status: dict[int, bool] = {}
    for n in (10, 20, 30, 40):
        baseline = one(phase2, "Phase-II baseline", n=n)
        baseline_pass = baseline["final_status"] == "controlled_consistency"
        chain_pass = baseline_pass
        result_rows.append(
            {
                "n": n,
                "role": "phase2_baseline",
                "k": baseline["k"],
                "one_minus_k": baseline["one_minus_k"],
                "clone_count": baseline["clone_count"],
                "horizon": baseline["horizon"],
                "residual": baseline["gc_residual"],
                "ci_low": baseline["gc_ci_low"],
                "ci_high": baseline["gc_ci_high"],
                "support_gate": baseline["support_gate"],
                "symmetry_gate": "pass" if baseline_pass else "fail",
                "time_gate": baseline["time_gate"],
                "population_gate": baseline["population_gate"],
                "numerical_control_gate": (
                    "pass" if baseline["timestep_gate"] == "pass"
                    and baseline["selection_gate"] == "pass" else "fail"
                ),
                "final_status": "pass" if baseline_pass else "fail",
            }
        )
        selected = [row for row in frozen if int(row["n"]) == n]
        if len(selected) != 2:
            raise ValueError(f"expected two selected pairs at n={n}")
        for grid in selected:
            k = float(grid["k_low"])
            clones = int(grid["clone_count"])
            horizon = float(grid["horizon"])
            pair = one(
                pairs, "primary pair", n=n, clone_count=clones,
                horizon=horizon, selection_time=2.0, dt=0.0005, k=k,
            )
            time = one(
                time_rows, "time comparison", n=n, clone_count=clones,
                selection_time=2.0, dt=0.0005, k=k,
            )
            support_pass = int(pair["paired_support_gate"]) == 1
            symmetry_pass = int(pair["paired_gc_gate"]) == 1
            time_pass = int(time["time_convergence_gate"]) == 1
            population_pass = True
            numerical_pass = True
            if grid["role"] == "outer":
                lower = clones // 2
                pair_pop = one(
                    population_pairs, "pair population comparison", n=n,
                    horizon=horizon, selection_time=2.0, dt=0.0005, k=k,
                    lower_clone_count=lower, upper_clone_count=clones,
                )
                member_low = one(
                    population_members, "low-k population comparison", n=n,
                    horizon=horizon, selection_time=2.0, dt=0.0005, k=k,
                    lower_clone_count=lower, upper_clone_count=clones,
                )
                member_high = one(
                    population_members, "high-k population comparison", n=n,
                    horizon=horizon, selection_time=2.0, dt=0.0005,
                    k=1.0-k, lower_clone_count=lower, upper_clone_count=clones,
                )
                population_pass = (
                    int(pair_pop["pair_population_gate"]) == 1
                    and int(member_low["population_gate"]) == 1
                    and int(member_high["population_gate"]) == 1
                )
                if n == 40:
                    numerical_pass = (
                        len(control_members) == 4
                        and all(int(row["member_convergence_gate"]) == 1
                                for row in control_members)
                        and len(control_pairs) == 2
                        and all(int(row["paired_control_gate"]) == 1
                                for row in control_pairs)
                    )
            pair_pass = (
                support_pass and symmetry_pass and time_pass
                and population_pass and numerical_pass
            )
            chain_pass = chain_pass and pair_pass
            result_rows.append(
                {
                    "n": n,
                    "role": grid["role"],
                    "k": k,
                    "one_minus_k": 1.0-k,
                    "clone_count": clones,
                    "horizon": horizon,
                    "residual": pair["residual"],
                    "ci_low": pair["residual_ci_low"],
                    "ci_high": pair["residual_ci_high"],
                    "support_gate": "pass" if support_pass else "fail",
                    "symmetry_gate": "pass" if symmetry_pass else "fail",
                    "time_gate": "pass" if time_pass else "fail",
                    "population_gate": (
                        "pass" if population_pass else "fail"
                    ) if grid["role"] == "outer" else "not_required",
                    "numerical_control_gate": (
                        "pass" if numerical_pass else "fail"
                    ) if n == 40 and grid["role"] == "outer" else "reused_phase2",
                    "final_status": "pass" if pair_pass else "fail",
                }
            )
        chain_status[n] = chain_pass

    write_csv(ROOT / "FINAL_SUMMARY.csv", result_rows)

    colors = {10: "tab:red", 20: "tab:orange", 30: "tab:green", 40: "tab:blue"}
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for n in (10, 20, 30, 40):
        rows = [row for row in result_rows if int(row["n"]) == n]
        xs = [float(row["k"]) for row in rows]
        ys = [float(row["residual"]) for row in rows]
        low = [float(row["ci_low"]) for row in rows]
        high = [float(row["ci_high"]) for row in rows]
        order = sorted(range(len(xs)), key=xs.__getitem__)
        ax.errorbar(
            [xs[i] for i in order], [ys[i] for i in order],
            yerr=[
                [ys[i]-low[i] for i in order],
                [high[i]-ys[i] for i in order],
            ], marker="o", capsize=3, color=colors[n], label=f"n={n}",
        )
    ax.axhline(0.0, color="black", ls="--", lw=1)
    ax.set_xlabel(r"lower complementary tilt $k$")
    ax.set_ylabel(r"$\psi_n(k)-\psi_n(1-k)$")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(ROOT / "multik_gc_residuals.png", dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.8), sharex=True)
    for n, ax in zip((10, 20, 30, 40), axes.flat):
        grids = [row for row in frozen if int(row["n"]) == n]
        ks = []
        for grid in grids:
            ks.extend([float(grid["k_low"]), float(grid["k_high"])])
        baseline = one(phase2, "Phase-II baseline plot row", n=n)
        ks.extend([float(baseline["k"]), float(baseline["one_minus_k"])])
        points=[]
        for k in sorted(set(ks)):
            grid_match = [g for g in grids if close(float(g["k_low"]),k)
                          or close(float(g["k_high"]),k)]
            if grid_match:
                g=grid_match[0]
                clones = int(g["clone_count"])
                horizon = float(g["horizon"])
            else:
                clones = int(baseline["clone_count"])
                horizon = float(baseline["horizon"])
            agg=one(
                aggregates,"SCGF aggregate",n=n,
                clone_count=clones, horizon=horizon, selection_time=2.0,
                dt=0.0005,k=k,
            )
            points.append((k,float(agg["mean_scgf"]),float(agg["run_scgf_se"])))
        ax.errorbar([p[0] for p in points],[p[1] for p in points],
                    yerr=[p[2] for p in points],marker="o",capsize=3,
                    color=colors[n])
        ax.set_title(f"n={n}")
        ax.grid(alpha=0.25)
    fig.supxlabel(r"tilt $k$")
    fig.supylabel(r"SCGF $\psi_n(k)$")
    fig.tight_layout()
    fig.savefig(ROOT / "multik_scgf.png", dpi=240)
    plt.close(fig)

    passed = all(chain_status.values())
    audit = {
        "chain_status": {str(n): "pass" if value else "fail"
                         for n, value in chain_status.items()},
        "all_chains_pass": passed,
        "pairs_reported": len(result_rows),
        "new_pairs": sum(row["role"] != "phase2_baseline" for row in result_rows),
    }
    (ROOT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")

    lines = [
        "# Phase-III multi-k validation report", "",
        f"Overall frozen-gate status: **{'PASS' if passed else 'FAIL'}**", "",
        "| n | pair role | pair | residual | 95% CI | support | time | population | numerical | final |",
        "|---:|:---|:---:|---:|:---:|:---:|:---:|:---:|:---:|:---:|",
    ]
    for row in result_rows:
        lines.append(
            f"| {row['n']} | {row['role']} | ({float(row['k']):g},"
            f"{float(row['one_minus_k']):g}) | {float(row['residual']):+.6f} "
            f"| [{float(row['ci_low']):+.6f},{float(row['ci_high']):+.6f}] "
            f"| {row['support_gate']} | {row['time_gate']} | "
            f"{row['population_gate']} | {row['numerical_control_gate']} "
            f"| {row['final_status']} |"
        )
    lines.extend([
        "", "## Claim boundary", "",
        "A pass supports numerical consistency with the GC SCGF symmetry at",
        "three resolved complementary tilt pairs for each tested finite chain.",
        "It does not prove exact symmetry for every tilt, uniformity in chain",
        "length, or the infinite-chain limit.", "",
    ])
    (ROOT / "VALIDATION_REPORT.md").write_text("\n".join(lines))
    verdict = [
        "# Final Phase-III verdict", "",
        f"Status: **{'PASS' if passed else 'FAIL'}**", "",
        "See `VALIDATION_REPORT.md`, `FINAL_SUMMARY.csv`, and the two multi-k",
        "figures for the frozen-gate result and its scope.", "",
    ]
    (ROOT / "FINAL_VERDICT.md").write_text("\n".join(verdict))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
