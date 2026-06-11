#!/usr/bin/env python3
"""Convert existing sparse .hist files to dense format (all NB^3 bins, zeros included).
Usage: python3 sparse_to_dense.py file1.hist [file2.hist ...]
Writes <name>_dense.hist alongside each input. Header lines are preserved."""
import sys
import numpy as np

NB = 80

for fn in sys.argv[1:]:
    header, H = [], np.zeros((NB, NB, NB), dtype=np.uint64)
    with open(fn) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            tok = s.split()
            if s.startswith('#') or len(tok) != 4 or not tok[0].isdigit():
                header.append(line.rstrip('\n'))
                continue
            i, j, k, c = tok
            H[int(i), int(j), int(k)] = int(c)
    out = fn.rsplit('.hist', 1)[0] + '_dense.hist'
    with open(out, 'w') as f:
        for h in header:
            if 'nonzero only' in h:
                h = '# format: ia ib it count (dense: all NB^3 bins, zeros included)'
            f.write(h + '\n')
        for ia in range(NB):
            for ib in range(NB):
                for it in range(NB):
                    f.write(f'{ia} {ib} {it} {H[ia, ib, it]}\n')
    print(f'{fn} -> {out}  (total counts {H.sum()})')