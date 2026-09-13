# Analysis implementation history

The scientific data and frozen parameters were unchanged throughout.

1. Analysis v1 stopped after the T=6 KDE returned zero density for 52 endpoint
   pairs.  A diagnostic confirmed the points were inside the coordinate grid;
   the finite four-bandwidth convolution support was the cause.  The failed
   output and log are preserved locally as `analysis_failed_v1*`.
2. The analysis was corrected to report unsupported points as a failed KDE
   gate and to forbid total-entropy statistics rather than extrapolating or
   adding a density floor.  Analysis v2 then stopped on a heterogeneous CSV
   field-name serialization error after computing the driven stability audit.
   Its output and log are preserved locally as `analysis_failed_v2*`.
3. The CSV writer was corrected to use the stable union of row fields.  The
   final analysis completed without changing any KDE, bootstrap, histogram,
support, or FT parameter.
