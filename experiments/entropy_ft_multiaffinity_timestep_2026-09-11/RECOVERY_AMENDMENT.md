# Recovery amendment: missing fine-timestep runs

The original four-process production pipeline exited with status 1 at
2026-09-12 14:38:32 EDT.  The two `dt=2.5e-4` cases had already completed and
their 4,000,256-row block files were preserved.  The two `dt=1.25e-4`
processes exited simultaneously before writing any CSV output; their original
logs were empty.  No kernel memory-pressure or kill message was found in the
available system log around the exit, so the external interruption cause is
unresolved.

The recovery does not change any scientific or statistical parameter.  It
reruns only the two missing cases with the exact frozen commands, seeds,
source, binary, and two-thread allocation.  Existing completed block files are
read-only inputs to the final integrity audit and analysis.  The initial
`pipeline.exitcode=1`, empty logs (renamed `*.initial_failed.log`), original
PID table, and original runtime table are retained as provenance.

Recovery cases:

- `weak_dt1p25e4`: `(T_L,T_R)=(6.5,5.5)`, `dt=0.000125`, seed `2026091102`.
- `moderate_dt1p25e4`: `(T_L,T_R)=(8,4)`, `dt=0.000125`, seed `2026091104`.

The recovery launcher retains the original AC-only stop/resume policy.  After
both missing runs finish, it audits all four frozen cases and runs the original
analysis/report pipeline without changing bins, support gates, bootstrap, or
fit definitions.
