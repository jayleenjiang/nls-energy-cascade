# Revision roadmap

## Gate A: canonical model and reproducible simulation

- Write one authoritative SDE specification in both full and reduced
  coordinates.
- Add a unit/sanity test showing equilibrium agreement with
  `exp(-H/(2T))`.
- Repair the flux simulator and record deterministic seeds and run metadata.
- Run short equilibrium and NESS smoke tests before expensive production runs.

## Gate B: long-chain evidence

- Produce converged action profiles for multiple chain lengths.
- Recompute current scaling with valid uncertainty and timestep sensitivity.
- Analyze finite-time averaged-current distributions at several window lengths.
- Recompute LTE distances with block/bootstrap uncertainty and fixed masks.
- Decide whether entropy production is completed or removed from the paper.

## Gate C: short-chain evidence

- Audit the 5D neural Fokker--Planck loss, normalization, validation slices,
  and held-out error metrics.
- Recompute symmetry-breaking and phase-locking statistics.
- Resolve or substantially narrow the generator/eigenfunction claims.

## Gate D: manuscript

- Add abstract, conclusion, limitations, methods table, reproducibility
  statement, and data/code availability.
- Replace every placeholder with a source-traceable figure.
- Move exploratory or weakly supported claims to a clearly labeled discussion.
- Standardize terminology: action/mass current versus Hamiltonian energy.
- Convert references to a maintainable BibTeX database with verified metadata.

## Gate E: final quality

- Compile cleanly with no unresolved references or overfull boxes.
- Visually inspect every PDF page.
- Run a fresh citation/data/claim integrity audit.
- Perform an independent reviewer-style methodology and presentation review.
- Commit and push only coherent, reproducible milestones.

