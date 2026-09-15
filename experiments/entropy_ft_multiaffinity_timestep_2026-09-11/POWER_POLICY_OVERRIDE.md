# Power-policy override

At 2026-09-13 00:30 UTC the user explicitly authorized the two active
fine-timestep recovery processes to continue regardless of AC/battery state.
This changes only the operational pause policy.  It does not change the model,
integrator, timestep, temperatures, burn-in, sample count, seeds, estimator,
or analysis gates.

The already-running recovery simulator PIDs are resumed in place.  The helper
`keep_recovery_running.sh` checks those recorded PIDs and sends `SIGCONT` only
when either is found in a stopped state.  It cannot launch a simulator and
exits when both recorded processes are gone.  Events are appended to
`unconditional_continuation.log`.
