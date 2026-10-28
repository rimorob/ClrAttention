#!/usr/bin/env python3
"""GRNBoost2 baseline (Moerman et al. 2019). Runs arboreto's own per-target
routine pieces, with the GRNBoost2 profile (regressor_type="GBM", SGBM_KWARGS,
early-stopping window 25, importances x n_trees, target removed from its
regressors), over a multiprocessing pool instead of dask (arboreto's dask layer
fails with current dask releases). Every gene is a candidate regulator.
Input: genes x samples TSV (first column = gene names).

--early_stop (decision of 27 Oct 2026, docs/plans/grnboost2_diagnostic.md):
  legacy   (default) arboreto's stopping rule with the out-of-bag improvement it was
           written for: the OOB loss of the SAME out-of-bag samples before and
           after each stage, as scikit-learn < 1.3 computed it (tools/grnboost2_variants.py,
           LegacyOOBMonitor). This is GRNBoost2 as published.
  sklearn  arboreto's stock EarlyStopMonitor on the installed scikit-learn. Under
           scikit-learn >= 1.3 this is the faulty behaviour (oob_improvement_ compares
           two different OOB subsets; most targets stop after 25-36 trees). Under
           scikit-learn < 1.3 it equals legacy. Used for the superseded runs.
  none     no early stopping (5,000 trees per target); sensitivity only.
    ~/venvs/clr-baselines/bin/python tools/grnboost2.py in.tsv out.tsv --workers 24 --seed 1 [--early_stop legacy]
"""
import argparse
import multiprocessing as mp
import os
import sys
import time
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from arboreto.core import (infer_partial_network, clean, to_links_df, SGBM_KWARGS,
                           EARLY_STOP_WINDOW_LENGTH)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grnboost2_variants import LegacyOOBMonitor

_M = None; _NAMES = None; _SEED = None; _MODE = None

def _init(M, names, seed, mode):
    global _M, _NAMES, _SEED, _MODE
    _M, _NAMES, _SEED, _MODE = M, names, seed, mode

def _one(j):
    name = _NAMES[j]; y = _M[:, j]
    if _MODE == "legacy":    # same steps as infer_partial_network: clean -> fit -> to_links_df
        Xc, nc = clean(_M, _NAMES, name)
        reg = GradientBoostingRegressor(random_state=_SEED, **SGBM_KWARGS)
        reg.fit(Xc, y, monitor=LegacyOOBMonitor(EARLY_STOP_WINDOW_LENGTH))
        return to_links_df("GBM", SGBM_KWARGS, reg, nc, name), float(len(reg.estimators_))
    win = EARLY_STOP_WINDOW_LENGTH if _MODE == "sklearn" else 10 ** 9
    df, meta = infer_partial_network("GBM", SGBM_KWARGS, _M, _NAMES, name, y, include_meta=True,
                                     early_stop_window_length=win, seed=_SEED)
    return df, float(meta["n_estimators"].iloc[0]) if len(meta) else np.nan

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--early_stop", choices=["legacy", "sklearn", "none"], default="legacy")
    a = ap.parse_args()
    t0 = time.time()
    X = pd.read_csv(a.inp, sep="\t", index_col=0)
    M = np.ascontiguousarray(X.T.values, dtype=np.float32)     # samples x genes
    names = list(X.index)
    with mp.get_context("fork").Pool(a.workers, initializer=_init, initargs=(M, names, a.seed, a.early_stop)) as pool:
        parts = pool.map(_one, range(len(names)), chunksize=4)
    net = pd.concat([p[0] for p in parts], ignore_index=True)
    net.to_csv(a.out, sep="\t", index=False)
    nt = np.array([p[1] for p in parts]); nt = nt[np.isfinite(nt)]
    q = (f"; trees per target median {np.median(nt):.0f} (IQR {np.percentile(nt, 25):.0f}-{np.percentile(nt, 75):.0f})"
         if len(nt) else "")
    import sklearn
    print(f"grnboost2[early_stop={a.early_stop}, scikit-learn {sklearn.__version__}] seed {a.seed}: "
          f"{M.shape[1]} genes x {M.shape[0]} samples -> {len(net)} edges, {time.time() - t0:.1f} s{q}")
