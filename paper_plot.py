"""Reproduce the purple/green HMF figure from released catalog histograms."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from hmf import MassFunction


def plot_release(samples_file, output):
    data = np.load(samples_file)
    edges, z_left = data['log_mass_edges'], data['z_left']
    mass = 10**(edges[:-1] + 0.05)
    lcdm = FlatLambdaCDM(H0=67.74, Om0=0.3089, Ob0=0.0486, Tcmb0=2.725)
    volumes = (lcdm.comoving_volume(z_left + 0.05) - lcdm.comoving_volume(z_left)).value * 0.6774**3
    conversion = volumes[None, :, None] * 0.1 * np.log(10) * mass[None, None, :]
    raw = data['raw_counts'] / conversion
    calibrated = data['calibrated_counts'] / conversion
    means = [raw.mean(axis=0), calibrated.mean(axis=0)]
    errors = [raw.std(axis=0), calibrated.std(axis=0)]
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['cmr10'],
                         'mathtext.fontset': 'cm', 'axes.formatter.use_mathtext': True,
                         'font.size': 14, 'font.weight': 'bold', 'axes.labelweight': 'bold'})
    fig, axes = plt.subplots(2, 4, figsize=(11, 6), sharex='col', sharey='row')
    for i, left in enumerate(z_left):
        common = dict(Mmin=edges[0]+0.05, Mmax=16.0+1e-6, dlog10m=0.1,
                      z=float(left+0.025), cosmo_model=lcdm, sigma_8=0.8159, n=0.9667,
                      transfer_model='CAMB', transfer_params={'extrapolate_with_eh': True})
        watson = MassFunction(hmf_model='Watson_FoF', **common)
        tinker = MassFunction(hmf_model='Tinker08', **common)
        if not np.allclose(watson.m, mass, rtol=1e-10):
            raise ValueError('HMF and histogram bin centers differ')
        for row in range(2):
            ax = axes[row, i]
            denominator = watson.dndm if row else np.ones_like(mass)
            ax.plot(mass, tinker.dndm/denominator, '--', color='#7F7F7F', label='Tinker08')
            for mean, error, color, label in zip(means, errors, ['#B48EAD', '#A3BE8C'], ['RFOF', 'Calibrated']):
                ax.errorbar(mass, mean[i]/denominator, yerr=error[i]/denominator, color=color, label=label, lw=1.5)
            ax.set_xscale('log')
            ax.set_xlim(10**edges[0], 10**15.5)
            ax.tick_params(which='both', direction='in')
            if row:
                ax.set_ylim(0, 2.75)
                ax.set_yticks(np.arange(0, 3, 0.5))
                ax.axhline(1, color='#C7C7C7', alpha=0.5, zorder=0)
                ax.set_xlabel(r'Mass [$M_\odot/h$]')
            else:
                ax.set_yscale('log')
                ax.set_ylim(1e-30, 1e-6)
                ax.set_yticks(10.**np.arange(-30, -9, 5))
    axes[0, 0].legend(frameon=False, fontsize=12)
    axes[0, 0].set_ylabel(r'$dn/dm\,[h^4{\rm Mpc}^{-3}M_\odot^{-1}]$')
    axes[1, 0].set_ylabel(r'$dn/dn_{\rm Watson\,FoF}$')
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['.png', '.pdf']:
        fig.savefig(output.with_suffix(suffix), dpi=180)
    plt.close(fig)
    return output.with_suffix('.png')
