#!/usr/bin/env python
"""Technical comparisons of the smoke-test attention outputs (no truth set is read, nothing is scored).

  python -I analysis/scfm_smoke_compare.py ~/scfm_models/smoke_out
Determinism (run1 vs run2, same inputs), observed vs expression-permuted-across-cells, and two different cell
subsets of the same size (seed 1 vs 2) as the yardstick for how much a change of input moves the mean attention.
"""
import sys, os, json
import numpy as np
from scipy.stats import spearmanr

d = os.path.expanduser(sys.argv[1])
rng = np.random.default_rng(0)


def top(M, k):
    M = M.copy(); np.fill_diagonal(M, -1); return set(np.argsort(M.ravel())[-k:])


for model, pre in (("scPRINT medium-v1.5", "attn_"), ("scGPT whole-human", "scgpt_attn_")):
    L = lambda t: np.load(f"{d}/{pre}{t}.npz", allow_pickle=True)["A"]
    a, b, p, s2 = L("run1"), L("run2"), L("perm"), L("seed2")
    G, nL = a.shape[-1], a.shape[0]
    idx = rng.choice(G * G, 200000, replace=False)
    m = lambda A: A.mean((0, 1))
    print(f"== {model}: A shape {a.shape}, NaN {int(np.isnan(a).sum())}, row sums [{a.sum(-1).min():.6f}, {a.sum(-1).max():.6f}]")
    print(f"determinism run1 vs run2: max abs diff {np.abs(a - b).max():.3g}, bitwise identical {np.array_equal(a, b)}")

    def cmp(x, y, name):
        print(f"{name}: relL1(all-layer-head mean) {np.abs(m(x) - m(y)).sum() / m(x).sum():.4f}, "
              f"Spearman {spearmanr(m(x).ravel()[idx], m(y).ravel()[idx])[0]:.4f}, top-1000 overlap {len(top(m(x), 1000) & top(m(y), 1000)) / 1000:.3f}")
        print("   per-layer relL1:", [round(float(np.abs(x[l] - y[l]).sum() / x[l].sum()), 3) for l in range(nL)])
    cmp(a, s2, "different 64-cell subset (seed1 vs seed2)")
    cmp(a, p, "observed vs permuted across cells         ")
    print(f"observed vs uniform 1/G: relL1 {np.abs(a - 1 / G).sum() / a.sum():.3f}")
