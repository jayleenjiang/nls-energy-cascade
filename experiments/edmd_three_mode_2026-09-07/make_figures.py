#!/usr/bin/env python3
"""Publication figures and validated eigenfunction slices for the EDMD audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CASES = ('driven_dt1e-3','driven_dt2p5e-4','equilibrium_dt1e-3','equilibrium_dt2p5e-4')
TITLES = {
    'driven_dt1e-3': r'driven, $\Delta t=10^{-3}$',
    'driven_dt2p5e-4': r'driven, $\Delta t=2.5\times10^{-4}$',
    'equilibrium_dt1e-3': r'equilibrium, $\Delta t=10^{-3}$',
    'equilibrium_dt2p5e-4': r'equilibrium, $\Delta t=2.5\times10^{-4}$'}


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def load_edmd_module(root):
    path = root / 'run_edmd.py'
    spec = importlib.util.spec_from_file_location('edmd_run', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules['edmd_run'] = module
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    analysis = root / 'analysis'
    figures = root / 'figures'
    figures.mkdir(exist_ok=True)

    plt.rcParams.update({
        'font.size': 8.5, 'axes.labelsize': 9, 'axes.titlesize': 9,
        'legend.fontsize': 7, 'figure.dpi': 160, 'savefig.dpi': 300,
        'axes.spines.top': False, 'axes.spines.right': False})

    # Lag convergence for all frozen reference candidates.
    lag = rows(analysis / 'lag_convergence.csv')
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True)
    for ax, case in zip(axes.flat, CASES):
        case_rows = [r for r in lag if r['case'] == case]
        modes = sorted(set(int(r['mode']) for r in case_rows))
        for mode in modes[:6]:
            rr = sorted([r for r in case_rows if int(r['mode']) == mode], key=lambda r: float(r['tau']))
            x = np.array([float(r['tau']) for r in rr])
            y = np.array([float(r['lambda_real']) for r in rr])
            ok = rr[0]['pass'] == 'True'
            ax.plot(x, y, marker='o', ms=3, lw=1, alpha=0.95 if ok else 0.45,
                    label=f'mode {mode}' + (' (lag pass)' if ok else ''))
        ax.axhline(0, color='0.65', lw=.7)
        ax.set_xscale('log')
        ax.set_title(TITLES[case])
        ax.set_ylabel(r'$\mathrm{Re}\,\lambda$')
        ax.legend(ncol=2, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel(r'EDMD lag $\tau$')
    fig.tight_layout()
    fig.savefig(figures / 'edmd_lag_convergence.png')
    fig.savefig(figures / 'edmd_lag_convergence.pdf')
    plt.close(fig)

    # Observable modal weights for the fine-timestep driven case.
    wr = rows(analysis / 'driven_dt2p5e-4_observable_modal_weights.csv')
    observables = ('I2','I1_plus_I3','cos_theta3','cos_even_sum')
    modes = sorted(set(int(r['reference_mode']) for r in wr))[:4]
    W = np.full((len(observables), len(modes)), np.nan)
    for i, obs in enumerate(observables):
        for j, mode in enumerate(modes):
            hit = [r for r in wr if r['observable'] == obs and int(r['reference_mode']) == mode]
            if hit:
                W[i, j] = float(hit[0]['weight_real'])
    fig, ax = plt.subplots(figsize=(5.0, 2.6))
    lim = max(.3, float(np.nanmax(np.abs(W))))
    im = ax.imshow(W, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
    ax.set_xticks(range(len(modes)), [f'mode {m}' for m in modes])
    ax.set_yticks(range(len(observables)),
                  [r'$I_2$',r'$I_1+I_3$',r'$\cos\theta_3$',r'$\cos\theta_1+\cos\theta_3$'])
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            ax.text(j, i, f'{W[i,j]:.2f}', ha='center', va='center', fontsize=7)
    ax.set_title('Driven fine-step EDMD modal covariance weights')
    fig.colorbar(im, ax=ax, label='signed fraction of reconstructed variance')
    fig.tight_layout()
    fig.savefig(figures / 'edmd_observable_weights.png')
    fig.savefig(figures / 'edmd_observable_weights.pdf')
    plt.close(fig)

    # Re-solve the primary matrices and save slices for every final reportable mode.
    edmd = load_edmd_module(root)
    summary = rows(analysis / 'leading_modes_summary.csv')
    selected = [r for r in summary if r['reportable'] == 'True'
                and r['case'] in ('driven_dt2p5e-4','equilibrium_dt2p5e-4')]
    slice_rows = []
    plot_payload = []
    angle = np.linspace(-math.pi, math.pi, 81)
    for case in ('driven_dt2p5e-4','equilibrium_dt2p5e-4'):
        with gzip.open(analysis / 'matrices' / f'{case}_D3_A.csv.gz', 'rt') as f:
            A = np.loadtxt(f, delimiter=',')
        with gzip.open(analysis / 'matrices' / f'{case}_D3_B_tau0.50.csv.gz', 'rt') as f:
            B = np.loadtxt(f, delimiter=',')
        fit = edmd.solve_edmd(A, B, .5, 1e-10, vectors=True)
        spec_rows = rows(analysis / f'{case}_spectra.csv')
        refs = []
        for r in spec_rows:
            if r['dictionary']=='D3' and float(r['tau'])==.5 and float(r['cutoff'])==1e-10:
                lam=complex(float(r['lambda_real']),float(r['lambda_imag']))
                mu=complex(float(r['mu_real']),float(r['mu_imag']))
                if -3<=lam.real<=-.05 and -1e-8<=lam.imag<=12 and 0<abs(mu)<=1.05:
                    refs.append(r)
                if len(refs)==12: break
        standard = rows(analysis / f'{case}_standardization.csv')
        mean = np.array([float(r['log_mean']) for r in standard])
        sd = np.array([float(r['log_sd']) for r in standard])
        for srow in [r for r in selected if r['case']==case]:
            mode = int(srow['mode'])
            erow = refs[mode]
            eig = int(erow['eigen_index'])
            lam = complex(float(erow['lambda_real']),float(erow['lambda_imag']))
            coeff = fit.coeff[:,eig].copy()
            pivot=int(np.argmax(np.abs(coeff)))
            coeff *= np.exp(-1j*np.angle(coeff[pivot]))
            maps=[]
            for action in (1.,2.,4.):
                vals=np.empty((len(angle),len(angle)))
                for a,th1 in enumerate(angle):
                    raw=np.zeros((len(angle),edmd.NCOLS))
                    raw[:,0]=np.cos(th1); raw[:,1]=np.sin(th1)
                    raw[:,2]=np.cos(angle); raw[:,3]=np.sin(angle)
                    raw[:,5]=action; raw[:,6]=2*action; raw[:,7]=0
                    raw[:,8]=raw[:,0]+raw[:,2]
                    value=edmd.features(raw,mean,sd)@coeff
                    vals[a]=value.real
                    for th3,val in zip(angle,value):
                        slice_rows.append((case,mode,lam.real,lam.imag,action,th1,th3,val.real,val.imag))
                maps.append(vals)
            plot_payload.append((case,mode,lam,maps))
    edmd.write_csv(analysis / 'reportable_eigenfunction_slices.csv',
                   ('case','mode','lambda_real','lambda_imag','I_equal','theta1','theta3','phi_real','phi_imag'),
                   slice_rows)

    if plot_payload:
        fig, axes = plt.subplots(len(plot_payload), 3, figsize=(7.2, 2.15*len(plot_payload)),
                                 squeeze=False, layout='constrained')
        for row,(case,mode,lam,maps) in enumerate(plot_payload):
            vmax=max(np.max(np.abs(m)) for m in maps)
            for col,(action,m) in enumerate(zip((1,2,4),maps)):
                ax=axes[row,col]
                im=ax.imshow(m,origin='lower',extent=(-math.pi,math.pi,-math.pi,math.pi),
                             cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto')
                if row == 0:
                    ax.set_title(rf'$I_1=I_2=I_3={action:g}$')
                ax.set_xlabel(r'$\theta_3$')
                if col==0:
                    bath = 'driven' if case.startswith('driven') else 'equilibrium'
                    ax.set_ylabel(bath + f' mode {mode}\n' +
                                  rf'$\lambda={lam.real:.3f}$' + '\n' + r'$\theta_1$')
            fig.colorbar(im,ax=axes[row,:].tolist(),shrink=.70,pad=.025,label=r'$\mathrm{Re}\,\phi$')
        fig.savefig(figures / 'edmd_reportable_eigenfunctions.png')
        fig.savefig(figures / 'edmd_reportable_eigenfunctions.pdf')
        plt.close(fig)


if __name__ == '__main__':
    main()
