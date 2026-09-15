# Formal V7 fine-timestep blind-test authorization

Frozen on 2026-09-14 after both V6 validation verdicts passed completely.

- driven validation verdict SHA-256:
  `15f252647c1ca927c0d26b96dfa1ea109c62db7400a58b6c6bf21ab9b475b0cf`;
- driven edges SHA-256:
  `b0b65d6e32126bfa2f8f691282900143e152c3fab9101183d8a2436e0c150830`;
- equilibrium validation verdict SHA-256:
  `b20baac548dbaf826086e34cbe78eba6b7b43a3877595a68b55530948fafe418`;
- equilibrium edges SHA-256:
  `1c13e837e1e9ad208770705437c2373f0e80be7b68ba416bff8b422ddfcc9720`.

The exact R9 factor-0.625 bundles and all original gates are unchanged.
`run_formal_v7_test.sh` is the only authorized reading of V6 test physical
statistics.  The resulting verdict is final for the fine-timestep model and
cannot be used for another adjustment.
