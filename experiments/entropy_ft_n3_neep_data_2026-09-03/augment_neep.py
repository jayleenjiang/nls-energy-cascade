#!/usr/bin/env python3
"""Append periodic encodings to the unchanged n=3 transition CSV."""

from __future__ import annotations

import csv
import math
import sys


RAW_HEADER = [
    "stream_id", "block_id", "q_left", "q_right", "delta_energy",
    "entropy_medium", "entropy_rate", "action_current",
    "energy_balance_error", "start_I1", "start_I2", "start_I3",
    "start_theta1", "start_theta3", "end_I1", "end_I2", "end_I3",
    "end_theta1", "end_theta3",
]
TRIG_HEADER = [
    "start_theta1_cos", "start_theta1_sin",
    "start_theta3_cos", "start_theta3_sin",
    "end_theta1_cos", "end_theta1_sin",
    "end_theta3_cos", "end_theta3_sin",
]


def main() -> int:
    reader = csv.reader(sys.stdin)
    writer = csv.writer(sys.stdout, lineterminator="\n")
    try:
        header = next(reader)
    except StopIteration:
        raise RuntimeError("empty input")
    if header != RAW_HEADER:
        raise RuntimeError(f"unexpected raw header: {header}")
    writer.writerow(RAW_HEADER + TRIG_HEADER)
    count = 0
    for row in reader:
        if len(row) != len(RAW_HEADER):
            raise RuntimeError(f"row {count + 2} has {len(row)} fields")
        angles = [float(row[index]) for index in (12, 13, 17, 18)]
        trig: list[str] = []
        for angle in angles:
            trig.extend((format(math.cos(angle), ".17g"),
                         format(math.sin(angle), ".17g")))
        writer.writerow(row + trig)
        count += 1
    print(f"augmented_rows={count}", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
