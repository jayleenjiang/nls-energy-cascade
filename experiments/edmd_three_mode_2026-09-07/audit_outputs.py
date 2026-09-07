#!/usr/bin/env python3
"""Integrity and provenance audit for the EDMD analysis artifacts."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


CASES = ('driven_dt1e-3','driven_dt2p5e-4','equilibrium_dt1e-3','equilibrium_dt2p5e-4')


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()


def nrows(path):
    with path.open() as f:
        return sum(1 for _ in f)-1


def main():
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,required=True); a=p.parse_args()
    root=a.root.resolve(); analysis=root/'analysis'; checks=[]
    for case in CASES:
        spec=analysis/f'{case}_spectra.csv'
        boot=analysis/f'{case}_bootstrap_raw.csv'
        gram=analysis/f'{case}_gram_diagnostics.csv'
        checks.append({'file':str(spec),'kind':'spectra','rows':nrows(spec),'pass':nrows(spec)>0})
        with boot.open() as f: br=list(csv.DictReader(f))
        modes=len(set(int(r['reference_mode']) for r in br))
        expected=500*modes
        checks.append({'file':str(boot),'kind':'bootstrap','rows':len(br),'expected':expected,
                       'matched':sum(int(r['matched']) for r in br),'pass':len(br)==expected})
        with gram.open() as f: gr=list(csv.DictReader(f))
        checks.append({'file':str(gram),'kind':'trivial','max_error':max(float(r['trivial_error']) for r in gr),
                       'pass':all(float(r['trivial_error'])<=1e-8 for r in gr)})
        for tag in ('A','B_tau0.10','B_tau0.25','B_tau0.50','B_tau1.00'):
            path=analysis/'matrices'/f'{case}_D3_{tag}.csv.gz'
            with gzip.open(path,'rt') as f:
                first=f.readline().rstrip().split(',')
                rows=1+sum(1 for _ in f)
            checks.append({'file':str(path),'kind':'matrix','rows':rows,'columns':len(first),
                           'pass':rows==500 and len(first)==500})
    summary=analysis/'leading_modes_summary.csv'
    checks.append({'file':str(summary),'kind':'summary','rows':nrows(summary),'pass':nrows(summary)>0})
    required=[root/'PROTOCOL.md',root/'run_edmd.py',root/'summarize_edmd.py',root/'make_figures.py',
              analysis/'gate_verdict.json',analysis/'reportable_eigenfunction_slices.csv',
              root/'figures'/'edmd_lag_convergence.pdf',root/'figures'/'edmd_observable_weights.pdf',
              root/'figures'/'edmd_reportable_eigenfunctions.pdf',
              root/'report'/'edmd_report.tex',root/'report'/'edmd_report.pdf',
              root/'FINAL_VERDICT.md',root/'VALIDATION_REPORT.md',root/'DELIVERY_CONTENTS.md']
    checks.extend({'file':str(x),'kind':'required','pass':x.exists()} for x in required)
    files=[]
    for path in sorted(set(Path(c['file']) for c in checks if Path(c['file']).exists())):
        files.append({'path':str(path.relative_to(root)),'sha256':digest(path),'bytes':path.stat().st_size})
    payload={'all_pass':all(c['pass'] for c in checks),'checks':checks,'artifacts':files,
             'input_source_sha256':'58c871882f1180701a7aae637a8f4be113315a72fb9d7815562f5fc90faf7976',
             'input_binary_sha256':'f7cd820bb8f716a8790f688db060259da5657e4f20fa09100546e66efb696b3f',
             'input_launch_commit':'899a5d345a98e4c8d5605a6a5a166a08bb86adcd'}
    prov=root/'provenance'; prov.mkdir(exist_ok=True)
    with (prov/'integrity_audit.json').open('w') as f: json.dump(payload,f,indent=2)
    with (prov/'analysis_outputs.sha256').open('w') as f:
        for x in files: f.write(f"{x['sha256']}  {x['path']}\n")
    print(json.dumps({'all_pass':payload['all_pass'],'checks':len(checks),'artifacts':len(files)},indent=2))
    raise SystemExit(0 if payload['all_pass'] else 1)


if __name__=='__main__': main()
