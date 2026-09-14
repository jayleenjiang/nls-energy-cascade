# Historical input-integrity failure

The frozen input audit completed before any density fit and returned **FAIL**
for one of the four historical cases.

## Raw findings

| Case | Rows | Minimum reconstructed action | Gate |
|---|---:|---:|---|
| driven, `dt=1e-3` | 6,400,064 | `1.7881393432617188e-07` | PASS |
| driven, `dt=2.5e-4` | 6,400,064 | `1.1548399925231934e-07` | PASS |
| equilibrium, `dt=1e-3` | 6,400,064 | `5.960464477539063e-08` | PASS |
| equilibrium, `dt=2.5e-4` | 6,400,064 | `0` | FAIL |

The two exact cancellations in the failed case are:

| Stream | Snapshot row | Stored sum | Stored difference | Failed reconstruction |
|---:|---:|---:|---:|---|
| 2 | 40629 | `0.8611177206039429` | `0.8611177206039429` | `I3=0` |
| 34 | 46751 | `11.036346435546875` | `-11.036346435546875` | `I1=0` |

Both streams belong to the density-training split. In the independent mode
split, stream 2 is test and stream 34 is train.

## Source-level cause

The controlled simulation evolves Cartesian state in float64 without an
action projection or floor. Its historical writer casts `I1+I3` and `I1-I3`
separately to float32. Near a vanishing endpoint action those two rounded
numbers can become equal in magnitude, so their later sum/difference loses the
small action exactly. The run metadata gives a true minimum sampled action of
`1.1709150955926485e-08`, zero radius events, zero projection/floor events and
zero midpoint/non-finite failures.

## Consequence

No density or mode fit is permitted from the historical representation under
the frozen strict-positivity gate. No row is dropped, clipped or imputed. The
predeclared repair in `AMENDMENT_001_DIRECT_STATE_OUTPUT.md` replays the same
trajectories and writes the five reduced coordinates directly in float64.

