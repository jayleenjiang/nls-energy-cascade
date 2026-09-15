#!/usr/bin/env python3
"""Formal-V3 regular trainer using combined train+development streams."""

import train_development
from common_v2 import rows


def combined_roles(phase, timestep):
    suffix = "dt2p5e-4" if timestep == "fine" else "dt1e-3"
    equilibrium = f"equilibrium_{suffix}"
    if phase == "driven":
        return (
            rows(f"driven_{suffix}", "train") + rows(f"driven_{suffix}", "validation"),
            rows(equilibrium, "train") + rows(equilibrium, "validation"),
        )
    raise ValueError("formal V3 trains only the driven ratio")


train_development.formal_training_roles = combined_roles
train_development.SCHEDULES.update({
    "pde100_m20_huber": {"pde": 1.0, "moment": 20.0, "loss": "huber"},
    "pde300_m40_huber": {"pde": 3.0, "moment": 40.0, "loss": "huber"},
    "pde2000_m80_huber": {"pde": 20.0, "moment": 80.0, "loss": "huber"},
    "pde5000_m160_huber": {"pde": 50.0, "moment": 160.0, "loss": "huber"},
})

if __name__ == "__main__":
    train_development.main()
