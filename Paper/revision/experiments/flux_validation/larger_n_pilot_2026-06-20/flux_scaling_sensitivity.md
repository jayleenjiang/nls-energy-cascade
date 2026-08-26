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
| tail n=20--50 | 20,30,40,50 | -2.03265 | [-2.07868, -1.98781] | 0.99935 |
| leave out n=10 | 20,30,40,50 | -2.03265 | [-2.07989, -1.98740] | 0.99935 |
| leave out n=20 | 10,30,40,50 | -1.88002 | [-1.90094, -1.86019] | 0.99853 |
| leave out n=30 | 10,20,40,50 | -1.90049 | [-1.92438, -1.87797] | 0.99823 |
| leave out n=40 | 10,20,30,50 | -1.89151 | [-1.91870, -1.86571] | 0.99722 |
| leave out n=50 | 10,20,30,40 | -1.85008 | [-1.86995, -1.83060] | 0.99801 |

## Adjacent local slopes

| interval | local exponent |
|---:|---:|
| 10--20 | -1.71976 |
| 20--30 | -1.92619 |
| 30--40 | -2.10925 |
| 40--50 | -2.12473 |

## Interpretation

All tested three-or-more-point fit windows remain faster than Fourier scaling, i.e. their fitted current exponent is more negative than `-1`.  The spread across windows is treated as finite-size sensitivity, separate from the bootstrap Monte Carlo intervals.
