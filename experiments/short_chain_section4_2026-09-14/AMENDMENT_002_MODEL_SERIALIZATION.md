# Amendment 002: Keras serialization-only repair

Date: 2026-09-14

This amendment was written after the first two frozen equilibrium fits had
finished their numerical optimization and validation calculations but failed
while serializing the trained model.  Both failed directories and their logs
are preserved with a `_failed_serialization_20260914` suffix.

## Failure

`NormalizedNLS3Density` is a subclassed Keras model.  Training calls its
component methods directly, so Keras did not mark the outer model as built.
The final `save_weights` call therefore raised:

```
ValueError: You are saving a model that has not yet been built.
```

The pre-repair SHA-256 of `stationary_density/normalized_density.py` is
`aa2b0560493b3d5fed1d2c8ddec304ba21fbad7c1a7ac4bc514cf4afded8c179`.

## Frozen repair

Immediately before `save_weights`, call the existing model `call` path once
on the fixed dummy state `(1,1,1,0,0)`.  This creates no data-dependent value,
does not alter a trainable parameter, and changes no optimizer, random seed,
architecture, loss, epoch, early-stopping, split, validation metric, or gate.
Then save weights exactly as originally specified.

Post-repair source SHA-256:
`735ffcb2d9a2af5898d8aad8d390080c50e70d73cd223d4b84f49bf26023ba94`.
A serialization smoke test first instantiated the conditional-layer weights
without calling the outer model (reproducing `built == false`), then called the
repaired `save_bundle`; the HDF5 weights were written, `built == true`, and
every pre-existing trainable array remained bitwise unchanged.

Because the in-memory fitted weights were lost when each failed process
exited, the two failed seeds must be rerun from the beginning with their exact
original commands and seeds.  They are not replaced by their previously
written metrics, and the complete 18-run matrix remains mandatory.  No test
stream has been read by the density branch.
