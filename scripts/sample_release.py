"""Read released catalogs without modifying them; export small validation samples.

Run on Perlmutter with the NERSC Python module. Only numpy and h5py are needed.
"""
import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np

M0 = 1.9487254007799495
SIM_MIN = np.log10(320 * 1e10 * M0)
EDGES = np.arange(SIM_MIN, 16.0 + 0.05, 0.1)
Z_LEFT = np.arange(4, dtype=float)


def descending_bound(dataset, value, strict=False):
    """First row with z <= value (or z < value when strict)."""
    lo, hi = 0, len(dataset)
    while lo < hi:
        mid = (lo + hi) // 2
        if dataset[mid] < value if strict else dataset[mid] <= value:
            hi = mid
        else:
            lo = mid + 1
    return lo


def stratified_sample(f, per_bin=32, chunk_size=1_000_000, seed=0, min_log_mass=None):
    rng = np.random.default_rng(seed)
    kept = np.empty((0, 5), dtype=np.float64)
    for start in range(0, len(f['redshift']), chunk_size):
        end = start + chunk_size
        z = f['redshift'][start:end].astype(float)
        lm = f['log10_m_rfof'][start:end].astype(float)
        target = np.log10(f['halo_mass_m200m'][start:end])
        if min_log_mass is not None:
            select = lm >= min_log_mass
            z, lm, target = z[select], lm[select], target[select]
            if not len(z):
                continue
        keys = np.floor(z / 0.05).astype(int) * 100 + np.clip(
            np.floor((lm - SIM_MIN) / 0.05).astype(int), 0, 99)
        priority = rng.random(len(z))
        order = np.lexsort((priority, keys))
        ordered_keys = keys[order]
        first = np.r_[0, np.flatnonzero(np.diff(ordered_keys)) + 1]
        rank = np.arange(len(order)) - np.repeat(first, np.diff(np.r_[first, len(order)]))
        take = order[rank < per_bin]
        new = np.column_stack((keys[take], priority[take], z[take], lm[take], target[take]))
        merged = np.concatenate((kept, new))
        order = np.lexsort((merged[:, 1], merged[:, 0]))
        merged = merged[order]
        first = np.r_[0, np.flatnonzero(np.diff(merged[:, 0])) + 1]
        rank = np.arange(len(merged)) - np.repeat(first, np.diff(np.r_[first, len(merged)]))
        kept = merged[rank < per_bin]
    return kept[:, 2:]


def histograms(f):
    raw = np.zeros((4, len(EDGES) - 1), dtype=np.int64)
    calibrated = raw.copy()
    slices = []
    z = f['redshift']
    probes = z[np.linspace(0, len(z) - 1, 1000, dtype=int)]
    if not np.all(np.diff(probes) <= 0):
        raise ValueError('Expected catalog sorted by decreasing redshift')
    for i, left in enumerate(Z_LEFT):
        a = descending_bound(z, left + 0.05)
        b = descending_bound(z, left, strict=True)
        slices.append([a, b])
        for start in range(a, b, 1_000_000):
            end = min(start + 1_000_000, b)
            zi = z[start:end]
            lm = f['log10_m_rfof'][start:end]
            valid = (zi >= left) & (zi <= left + 0.05) & (lm > SIM_MIN)
            raw[i] += np.histogram(lm[valid], bins=EDGES)[0]
            calibrated[i] += np.histogram(
                np.log10(f['halo_mass_m200m'][start:end][valid]), bins=EDGES)[0]
    return raw, calibrated, slices


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog-dir', type=Path, default=Path('/global/cfs/cdirs/cmb/halfdome/full_res/halos'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sample-seeds', type=int, nargs='+', default=[100, 102])
    p.add_argument('--per-bin', type=int, default=32)
    p.add_argument('--min-sample-log-mass', type=float)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    records, arrays, raw, calibrated = [], {}, [], []
    seeds = list(range(100, 122, 2))
    for seed in seeds:
        path = args.catalog_dir / f'lightcone_{seed}.hdf5'
        stat = path.stat()
        with h5py.File(path, 'r') as f:
            r, c, slices = histograms(f)
            raw.append(r)
            calibrated.append(c)
            if seed in args.sample_seeds:
                sample = stratified_sample(f, args.per_bin, seed=seed, min_log_mass=args.min_sample_log_mass)
                arrays[f'sample_{seed}'] = sample
                print(f'seed {seed}: {len(sample)} stratified pairs', flush=True)
            records.append({'seed': seed, 'path': str(path), 'size': stat.st_size,
                            'mtime_ns': stat.st_mtime_ns, 'rows': len(f['redshift']),
                            'plot_slices': slices})
        print(f'seed {seed}: histograms complete, {time.monotonic()-start:.1f}s', flush=True)
    np.savez_compressed(args.output / 'release_samples.npz', **arrays,
                        seeds=seeds, log_mass_edges=EDGES, z_left=Z_LEFT,
                        raw_counts=raw, calibrated_counts=calibrated)
    metadata = {'catalogs': records, 'sample_columns': ['redshift', 'log10_m_rfof', 'log10_m200m'],
                'sample_per_bin': args.per_bin, 'particle_mass': M0 * 1e10,
                'elapsed_seconds': time.monotonic()-start}
    (args.output / 'release_samples.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'elapsed': metadata['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
