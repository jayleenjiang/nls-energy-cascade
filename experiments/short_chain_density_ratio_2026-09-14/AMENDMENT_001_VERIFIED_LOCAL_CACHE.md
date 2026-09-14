# Amendment 001: verified local read-only cache

The first equilibrium blind evaluation produced no output because macOS file
coordination blocked indefinitely while reading
`equilibrium_dt2p5e-4/stream_005.f64` from the repository's Documents tree.
An independent SHA-256 read of the same file blocked as well.  This occurred
after equilibrium fitting and before any blind metric was calculated.

The already delivered, losslessly compressed Section-4 data package was
integrity-tested and its `raw/direct_state_f64` subtree was extracted to:

`/Users/jayleenjiang/NLS_section4_density_cache/short_chain_section4_2026-09-14`

All 256 extracted files were checked against the frozen
`provenance/DIRECT_STATE_SHA256.tsv` manifest.  Every byte count and SHA-256
matches; see `provenance/CACHE_AUDIT.json`.  The loader now accepts the
environment variable `NLS_SECTION4_RAW_CACHE` to resolve those identical bytes
outside the synchronized Documents tree.  No state, stream, split, seed,
model, optimizer, gate, or statistic changes.  All subsequent commands set
this variable explicitly.  The failed no-output blind attempt remains recorded
in `logs/equilibrium_blind_first_io_block.log`.

