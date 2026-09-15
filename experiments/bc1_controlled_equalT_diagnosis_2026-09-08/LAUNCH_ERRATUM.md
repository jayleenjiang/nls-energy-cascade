# Launch erratum

The first detached launch on 2026-09-08 started the first eight `n=25` process
commands and then the detached process tree disappeared before producing any
simulation artifact. The eight per-run logs were all zero bytes; no profile,
metadata, final-state, or summary file existed. Those empty logs are preserved
under `provenance/failed_detached_launch_1/`.

This was classified as a launcher/session-lifetime failure, not a numerical
result. The identical frozen pipeline was restarted in a persistent PTY with
the same source, binary, parameters, seeds, and gates. No scientific parameter
was changed and no completed batch was overwritten.
