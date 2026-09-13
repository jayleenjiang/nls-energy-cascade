# Final verdict: n=3 direct finite-time entropy sampling

## Integrity

PASS after the documented analysis-only roundoff-gate repair.  The simulation
completed 8/8 batches and wrote 1,000,064 finite, ordered rows.  There were
zero midpoint failures.  Entropy and first-law identities recompute within a
scale-aware binary64/CSV roundoff bound, and consecutive reduced endpoints
match exactly after parsing.  The raw block SHA-256 is
`4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f`.

## Direct negative-tail result

| t | evaluated blocks | negative medium-entropy blocks | probability | resolved by frozen gate |
|---:|---:|---:|---:|:---:|
| 20 | 1,000,064 | 41 | 4.0997376e-5 | no |
| 40 | 499,968 | 0 | 0 | no |
| 60 | 333,312 | 0 | 0 | no |
| 80 | 249,984 | 0 | 0 | no |
| 100 | 199,936 | 0 | 0 | no |
| 120 | 166,656 | 0 | 0 | no |
| 140 | 142,848 | 0 | 0 | no |
| 160 | 124,928 | 0 | 0 | no |
| 180 | 111,104 | 0 | 0 | no |
| 200 | 99,968 | 0 | 0 | no |

At `t=20`, the 95% independent-stream bootstrap interval for the probability
is `[2.8998144e-5,5.3996544e-5]`.  There are zero symmetric FD bin pairs with
at least 20 raw counts on each side.  Therefore a detailed-symmetry slope,
confidence interval, and R-squared are **not estimable by direct sampling**.
No fit window was selected.

## First law

At `t=20`, the residual `Q_left+Q_right-Delta E` has mean `-1.4168990e-5`,
standard deviation `2.9704513e-4`, and RMS `2.9738271e-4`.  The residual RMS is
`2.2442619e-6` of the RMS left-bath heat.  Complete horizon-by-horizon raw
distribution summaries are in `analysis/first_law_residuals.csv`.

## Medium-only exponential average

At `t=20`, `log mean exp(-Sigma_m)=-5.7735605` with stream-bootstrap interval
`[-8.0681734,-4.8380775]`.  Its exponential-weight ESS is only `2.2880` out of
1,000,064 and the largest sample contributes `61.55%` of the weight.  At longer
horizons the ESS remains approximately one to five.  These are unresolved
rare-event estimates and medium entropy is not total entropy.

## Total entropy and cloning

The requested total-entropy test is not evaluated.  The saved Section-4 NN is
for `(T1,T3)=(2,8)`, not `(10,2)`, and its documented log-density/FP errors are
not adequate for endpoint entropy at the requested precision.  It was not
applied to these data.

Direct sampling therefore fails the support-feasibility test at `n=3`.
Cloning or another rare-event method may now be considered explicitly, but no
cloning run was started or substituted in this experiment.  The data neither
verify nor refute the fluctuation theorem.
