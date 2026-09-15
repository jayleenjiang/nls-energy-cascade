# Fresh fine-timestep held-out trajectory protocol for formal V7

Frozen on 2026-09-14 after R9 selected calibration factor 0.625 using only
opened development data and before any V6 state exists.

Use the unchanged validated binary `NLS_stationary_state5_f64`, SHA-256
`58c22764eb864f156bc1df4a0190a0ac1851ff2505fda5fe07b9fde2e4469c82`,
and source SHA-256
`0d7000b628257311515e8edd43492a348c9bdbc3903c95c06bff5c718bf1f8a2`.

Generate two independent 64-stream cases with gamma 0.1, dt 2.5e-4, burn-in
500, stationary duration 1000, snapshot spacing 0.01, four SIMD batches and
four OpenMP threads:

- driven fine `(T1,T3,seed)=(2,8,2026091701)`;
- equilibrium fine `(T1,T3,seed)=(5,5,2026091702)`.

Rank streams by SHA-256 of `formal-v6|case|stream_id` before reading values;
the first 32 are validation and the remaining 32 are test.  Test remains blind
until the exact selected R9 bundles pass all unchanged validation gates.

Every stream file must have shape 100001-by-5, finite values and strictly
positive actions.  All projection, floor, zero-radius, midpoint-failure and
non-finite counters must be zero.  The integrity audit may establish only these
conditions and may not compute physical statistics.
