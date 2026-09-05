# n=3 direct finite-time entropy experiment

This prospective experiment extends the direct bath-heat sampler to the
five-dimensional `n=3` reduced state.  Read `PROTOCOL.md` before interpreting
any output.  `MODEL_COMPATIBILITY_AUDIT.md` explains why the saved Section-4
`(2,8)` NN density cannot be used for the requested `(10,2)` total-entropy
calculation.

Frozen production command:

```bash
/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_n3_direct_2026-08-31/bin/entropy_ft_n3 \
  sample_n3 10 2 3 8 500 20 7813 0.0005 2026083133 8 \
  /Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_n3_direct_2026-08-31/raw/n3 1
```

The AC-aware launcher is `run_when_ac.sh`.  It waits without altering the
frozen parameters and executes `run_production.sh` once external power is
available.
