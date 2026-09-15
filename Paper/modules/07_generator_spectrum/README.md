# Module 07: backward-generator relaxation and eigenfunction surrogate

## Manuscript coverage

Eigenfunction of the generator and slow-relaxation diagnostics.

## Experiment

- Backward Monte Carlo evolves ensembles from 4096 initial points and records
  the observable E_x[cos(theta1(X_t))].
- Exponential fits produce an observable-dependent decay rate and amplitudes.
- A neural surrogate is trained for the corresponding amplitude field Q.

## Key limitation

The median retained rate depends strongly on the fitting window: it is about
-1.60 on [0,2], -0.934 on [0,5], and -0.676 on [0.5,5].  Therefore the current
paper can call -0.934 a diagnostic relaxation rate for cos(theta1), but not a
resolved spectral gap or long-chain eigenvalue.

The Q surrogate is also only moderately accurate in the archived rerun.  Its
figures show qualitative structure rather than a high-accuracy eigenfunction.

## Code/data/output chain

- Simulation: cpp/backward/NLS_backward.cpp.
- Initial-point/data helpers: python/data_gen.
- Neural solver: NN notebooks/FKE_eigen.ipynb.
- Raw curves and amplitudes: KDE/backward_NLS_X.txt,
  KDE/4:23_eigen/NLS_backward_Y_train.txt, KDE/backward_NLS_Q1.txt.
- Fit-window audit: Paper/revision/eigen_fit_sensitivity.json.
- Figures: eigenvalue_scatter.png and Q1_slices.png.

## Status

This section is LOCAL_ENHANCED relative to the advisor-facing draft.  Present
the window-sensitivity limitation before migrating it.
