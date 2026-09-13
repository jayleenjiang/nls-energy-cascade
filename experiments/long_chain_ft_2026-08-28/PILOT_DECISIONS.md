# Pilot decisions for long-chain GC sampling

This file records proposal choices separately from the fluctuation-symmetry
test, so the sampler is not tuned to make the GC residual look small.

## Initial frozen settings

- observable: `Sigma_R = -(1/T_R-1/T_L) Q_R`;
- temperatures: `T_L=10`, `T_R=2`;
- timestep: `dt=0.0005` for the primary series;
- initial selection interval: `0.5`;
- fixed resampling at every selection event;
- primary nontrivial tilt pair: `(k,1-k)=(0.3,0.7)`;
- proposal control scale for replacement production: `0.50`.

## Why scale 0.25 was rejected

The first four-seed `n=10`, `N_c=1024`, `t=40` diagnostic had a minimum
selection-weight ESS of `97.84` at `k=0.7`, below the predeclared floor
`0.1 N_c = 102.4`.  Its mean cumulative paired residual was `-0.01333` at
`t=40` (95% Welch interval `[-0.02218,-0.00447]`), so it was not accepted.
The late-half residual was smaller (`-0.00666`) but this did not override the
failed cumulative and support gates.

## Support-only control refinement

At `n=10`, `N_c=1024`, `k=0.7`, `t=40`, identical pilot seed `87101` gave:

| control scale | minimum weight ESS | final unique roots | final root-weight ESS | diagnosis |
|---:|---:|---:|---:|:---|
| 0.50 | 615.03 | 29 | 19.19 | retained |
| 0.75 | 535.06 | 13 | 1.51 | rejected: severe state/genealogy collapse |

The earlier scale `0.25` four-seed diagnostic had substantially lower and less
stable high-tilt ESS, including the failure above.  Scale `0.50` was selected
from these support diagnostics only.  The corresponding SCGF values were not
used as a tuning objective.

All rows in the old `production_n10/N1024` directory are therefore treated as
diagnostic pilot output.  Accepted results must come from a replacement series
using the frozen scale `0.50` and must independently pass every gate in
`PREDECLARED_GATES.md`.

## Selection-interval amendment

Longer `n=10` pilots showed that resampling every `0.5` time unit caused rapid
genealogical collapse even when the per-selection weight ESS was healthy.
No-resampling and adaptive-resampling pilots retained weights for too long and
instead developed very low cumulative ESS at the high tilt.  A support-only
comparison therefore selected `selection_time=2.0` for replacement production.
At this interval the controlled estimator also reproduced the independent
million-block direct estimate at `t=20` within uncertainty.  The paired GC
residual was not used as the selection criterion.

The production settings after this amendment are:

- `dt=0.0005`;
- `selection_time=2.0`;
- fixed resampling at each selection event;
- exact finite-step Gaussian likelihood correction;
- `gauge_shift=0.1` (right-bath entropy-current gauge);
- `control_scale=0.50`;
- four independent seeds at each `(n,k,N_c)`;
- population doubling from `N_c=512` to `N_c=1024` (and the existing
  `N_c=1024` to `2048` comparison for `n=10`).

## Resolved tilt interval

The pair `(0.3,0.7)` remains the strongest target.  It passes the production
audit for `n=10`, but the `k=0.7` member is support/finite-time limited for
longer chains.  This is recorded as a failed strong-pair diagnostic, not
discarded.  The predeclared secondary pair `(0.4,0.6)` is therefore used for
`n=20,30,40`, conditional on the same support, independent-seed, time, and
population gates.  This narrows the numerical claim; it does not retroactively
change the original strongest target.

## n=40 late-time failure and remedial extension

The first completed `n=40`, `N_c=1024`, `t=80` production series is retained
as a failed diagnostic.  Its full-window intervals include zero at `t=60` and
`t=80`, and its `t=60` to `t=80` time-comparison gate passes.  However, the
late-half residual is `0.01527` at `t=60` and `0.01013` at `t=80`, above the
frozen absolute tolerance `0.01`.  No threshold or rounding convention is
changed.

The remedial plan was fixed before further seeds were inspected: rerun four
new independent seed pairs at `t=120` for both `N_c=512` and `N_c=1024`, using
seed bases `89600` and `89700`.  All other settings remain frozen.  The prior
failed series is not pooled with this extension.

Two old `N_c=512` run-3 files were automatically evicted by macOS/iCloud and
became local `dataless` placeholders.  Their already generated aggregate audit
is retained, but they are not used as a substitute for raw records.  The same
deterministic seeds (`88605`, `88606`) are rerun into a separate recovery
directory and compared against the prior aggregate values; no original file is
deleted or overwritten.
