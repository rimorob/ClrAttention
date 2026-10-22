#!/usr/bin/env python3
"""GRNBoost2 diagnostic variants (analysis/grnboost2_check.R). Input/output as
tools/grnboost2.py: genes x samples TSV in, TSV [TF, target, importance] out.

  --mode legacy_es  arboreto's GRNBoost2 profile (GBM, SGBM_KWARGS, window 25,
                    importances x n_trees, target removed from its regressors)
                    but with the out-of-bag early-stopping signal computed as
                    scikit-learn < 1.3 did and arboreto assumes: the OOB loss of
                    the SAME out-of-bag samples before and after each stage.
                    scikit-learn >= 1.3 defines oob_improvement_[i] as
                    oob_scores_[i-1] - oob_scores_[i], two different random 10%
                    OOB subsets, so arboreto's 25-round window mean telescopes to
                    loss(subset lo-1) - loss(subset hi): sampling noise, and the
                    monitor stops at a near-random early round.
  --mode no_es      arboreto profile with early stopping disabled (window longer
                    than n_estimators = 5000; importances x 5000, a constant).
  --mode dask       arboreto.algo.grnboost2 itself on a local dask
                    Client(processes=False); one patch: arboreto calls
                    from_delayed([]) for the unused meta frame, which current
                    dask rejects, so that call returns None when empty.
  --mode wrapper    in-process replica of tools/grnboost2.py (stock monitor).
Every gene is a candidate regulator in all modes (tf_names = all genes).
    ~/venvs/clr-baselines/bin/python tools/grnboost2_variants.py in.tsv out.tsv --mode legacy_es --workers 8 --seed 1
"""
import argparse
import multiprocessing as mp
import time
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from arboreto.core import (infer_partial_network, clean, to_links_df, SGBM_KWARGS,
                           EARLY_STOP_WINDOW_LENGTH)

class LegacyOOBMonitor:
    """arboreto's EarlyStopMonitor rule (stop when the mean OOB improvement over
    the last `window` stages is < 0), on same-sample OOB improvements."""
    def __init__(self, window=EARLY_STOP_WINDOW_LENGTH):
        self.window = window; self.imp = []; self.prev = None

    def __call__(self, i, est, env):
        oob = ~env["sample_mask"]; y = env["y"]; w = env["sample_weight"]
        f_new = np.asarray(env["raw_predictions"]).ravel()
        if self.prev is None:     # stage 0: undo this stage's tree
            f_old = f_new - est.learning_rate * est.estimators_[i, 0].predict(env["X"])
        else:
            f_old = self.prev
        ww = w[oob]
        lo = np.sum(ww * (y[oob] - f_old[oob]) ** 2) / np.sum(ww)
        ln = np.sum(ww * (y[oob] - f_new[oob]) ** 2) / np.sum(ww)
        self.imp.append(lo - ln); self.prev = f_new.copy()
        if i >= self.window - 1:
            return np.mean(self.imp[-self.window:]) < 0
        return False

_M = None; _NAMES = None; _SEED = None; _MODE = None

def _init(M, names, seed, mode):
    global _M, _NAMES, _SEED, _MODE
    _M, _NAMES, _SEED, _MODE = M, names, seed, mode

def _one(j):
    name = _NAMES[j]; y = _M[:, j]
    if _MODE == "wrapper":
        df, meta = infer_partial_network("GBM", SGBM_KWARGS, _M, _NAMES, name, y, include_meta=True,
                                         early_stop_window_length=EARLY_STOP_WINDOW_LENGTH, seed=_SEED)
        return df, float(meta["n_estimators"].iloc[0]) if len(meta) else np.nan
    if _MODE == "no_es":
        df, meta = infer_partial_network("GBM", SGBM_KWARGS, _M, _NAMES, name, y, include_meta=True,
                                         early_stop_window_length=10 ** 9, seed=_SEED)
        return df, float(meta["n_estimators"].iloc[0]) if len(meta) else np.nan
    # legacy_es: same steps as infer_partial_network (clean -> fit -> to_links_df)
    Xc, nc = clean(_M, _NAMES, name)
    reg = GradientBoostingRegressor(random_state=_SEED, **SGBM_KWARGS)
    reg.fit(Xc, y, monitor=LegacyOOBMonitor(EARLY_STOP_WINDOW_LENGTH))
    return to_links_df("GBM", SGBM_KWARGS, reg, nc, name), float(len(reg.estimators_))

def run_pool(M, names, seed, mode, workers):
    with mp.get_context("fork").Pool(workers, initializer=_init, initargs=(M, names, seed, mode)) as pool:
        parts = pool.map(_one, range(len(names)), chunksize=4)
    net = pd.concat([p[0] for p in parts], ignore_index=True)
    return net, np.array([p[1] for p in parts])

def run_dask(X, seed, workers):
    from distributed import Client
    import arboreto.core as ac
    from arboreto.algo import grnboost2
    fd = ac.from_delayed
    ac.from_delayed = lambda dfs, meta: fd(dfs, meta=meta) if len(dfs) else None
    df = pd.DataFrame(np.ascontiguousarray(X.T.values, dtype=np.float64), columns=list(X.index))
    c = Client(processes=False, n_workers=1, threads_per_worker=workers, dashboard_address=None)
    try:
        net = grnboost2(df, client_or_address=c, seed=seed)
    finally:
        c.close()
    return net.reset_index(drop=True), np.array([np.nan])

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--mode", required=True, choices=["wrapper", "legacy_es", "no_es", "dask"])
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    X = pd.read_csv(a.inp, sep="\t", index_col=0)
    t0 = time.time()
    if a.mode == "dask":
        net, ntree = run_dask(X, a.seed, a.workers)
    else:
        M = np.ascontiguousarray(X.T.values, dtype=np.float32)     # samples x genes, as tools/grnboost2.py
        net, ntree = run_pool(M, list(X.index), a.seed, a.mode, a.workers)
    net[["TF", "target", "importance"]].to_csv(a.out, sep="\t", index=False)
    nt = ntree[np.isfinite(ntree)]
    q = (f"; trees per target median {np.median(nt):.0f} (IQR {np.percentile(nt, 25):.0f}-{np.percentile(nt, 75):.0f})"
         if len(nt) else "")
    print(f"grnboost2[{a.mode}] seed {a.seed}: {X.shape[0]} genes x {X.shape[1]} samples -> {len(net)} edges, "
          f"{net['target'].nunique()} targets, {time.time() - t0:.1f} s{q}")
