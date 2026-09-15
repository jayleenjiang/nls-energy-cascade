# Current-scaling fit-window sensitivity

## Material Passport

- Artifact type: finite-size robustness analysis
- Model version: `gibbs-canonical-v1`
- Verification status: VERIFIED against canonical summary/sample files
- Scope: fit-window sensitivity only; not an asymptotic theorem

## Fit-window results

| label | chain lengths | exponent | 95% bootstrap CI | R^2 |
|---|---:|---:|---:|---:|
| primary n=10--40 | 10,20,30,40 | -1.85008 | [-1.87032, -1.83081] | 0.99801 |
| with n=50 | 10,20,30,40,50 | -1.89449 | [-1.91717, -1.87295] | 0.99761 |
| with n=50,60 | 10,20,30,40,50,60 | -1.92956 | [-1.95424, -1.90603] | 0.99739 |
| tail n=20--60 | 20,30,40,50,60 | -2.05926 | [-2.10697, -2.01242] | 0.99933 |
| leave out n=10 | 20,30,40,50,60 | -2.05926 | [-2.10798, -2.01274] | 0.99933 |
| leave out n=20 | 10,30,40,50,60 | -1.91232 | [-1.93624, -1.88943] | 0.99799 |
| leave out n=30 | 10,20,40,50,60 | -1.92964 | [-1.95427, -1.90553] | 0.99816 |
| leave out n=40 | 10,20,30,50,60 | -1.93303 | [-1.96062, -1.90657] | 0.99734 |
| leave out n=50 | 10,20,30,40,60 | -1.92214 | [-1.95249, -1.89317] | 0.99703 |
| leave out n=60 | 10,20,30,40,50 | -1.89449 | [-1.91685, -1.87324] | 0.99761 |

## Adjacent local slopes

| interval | local exponent |
|---:|---:|
| 10--20 | -1.71976 |
| 20--30 | -1.92619 |
| 30--40 | -2.10925 |
| 40--50 | -2.12473 |
| 50--60 | -2.17771 |

## Interpretation

All tested three-or-more-point fit windows remain faster than Fourier scaling, i.e. their fitted current exponent is more negative than `-1`.  The spread across windows is treated as finite-size sensitivity, separate from the bootstrap Monte Carlo intervals.
