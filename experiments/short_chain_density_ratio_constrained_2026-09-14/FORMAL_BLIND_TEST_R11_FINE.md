# R11 fine three-moment final test authorization

Frozen after both R11 validation verdicts passed and before any R11 test
statistic was read.

- driven validation verdict SHA-256:
  `09392334b1a7cecba9098ac62e7aa23b7aebb883384c7f17aa5b8c6cf8c9f09c`;
- driven validation edges SHA-256:
  `b0b65d6e32126bfa2f8f691282900143e152c3fab9101183d8a2436e0c150830`;
- equilibrium validation verdict SHA-256:
  `f1b7cac596d919cdcd08dbe50aad625674d8951dbfdcea81379e6c56d59e1dbf`;
- equilibrium validation edges SHA-256:
  `1c13e837e1e9ad208770705437c2373f0e80be7b68ba416bff8b422ddfcc9720`.

The models, three calibration observables, 0.625 shrink factor, trajectory
assignments, edges, normalizer treatment, bootstrap, and gates are exactly
those frozen in `RECOVERY_PROTOCOL_R11_FINE_THREE_MOMENT.md`.

`run_recovery_r11_fine_test.sh` is the one additional and final authorized
R11 reading of the V6 test physical statistics.  The test result cannot be
used to change or retry the model.  Both PASS and FAIL outcomes are final for
this protocol and must be retained.

