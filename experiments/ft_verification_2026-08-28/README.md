# Fluctuation-theorem verification synthesis

This directory synthesizes two independent positive controls:

1. an exact finite-step forward/reverse path-probability experiment for the
   split Cartesian integrator;
2. a short-time `n=2` driven-NESS total-entropy experiment with an independently
   validated stationary endpoint density.

The report source is `report/ft_verification_report.tex`.  Large raw trajectory
files remain in the source experiment directories and are ignored by Git; all
audits, summaries, figures, source code, and exact commands are versioned.

Accepted scope: numerical verification of finite-time total-entropy
fluctuation relations for `n=2`.  The long-chain asymptotic medium-entropy
Gallavotti--Cohen symmetry remains unresolved.
