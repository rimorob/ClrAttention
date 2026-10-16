#!/usr/bin/env python3
"""GRNBoost2 baseline (Moerman et al. 2019). Runs arboreto's own per-target
routine, arboreto.core.infer_partial_network, with the GRNBoost2 profile
(regressor_type="GBM", SGBM_KWARGS, early stopping) -- exactly what
arboreto.algo.grnboost2 schedules -- over a multiprocessing pool instead of
dask (arboreto's dask layer fails with current dask releases). Every gene is a
candidate regulator. Input: genes x samples TSV (first column = gene names).
    ~/venvs/clr-baselines/bin/python tools/grnboost2.py in.tsv out.tsv --workers 24 --seed 1
"""
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
from arboreto.core import infer_partial_network, SGBM_KWARGS, EARLY_STOP_WINDOW_LENGTH

_M = None; _NAMES = None; _SEED = None

def _init(M, names, seed):
    global _M, _NAMES, _SEED
    _M, _NAMES, _SEED = M, names, seed

def _one(j):
    return infer_partial_network(regressor_type="GBM", regressor_kwargs=SGBM_KWARGS,
                                 tf_matrix=_M, tf_matrix_gene_names=_NAMES,
                                 target_gene_name=_NAMES[j], target_gene_expression=_M[:, j],
                                 include_meta=False, early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
                                 seed=_SEED)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    X = pd.read_csv(a.inp, sep="\t", index_col=0)
    M = np.ascontiguousarray(X.T.values, dtype=np.float32)     # samples x genes
    names = list(X.index)
    with mp.get_context("fork").Pool(a.workers, initializer=_init, initargs=(M, names, a.seed)) as pool:
        parts = pool.map(_one, range(len(names)), chunksize=4)
    net = pd.concat(parts, ignore_index=True)
    net.to_csv(a.out, sep="\t", index=False)
    print(f"grnboost2: {M.shape[1]} genes x {M.shape[0]} samples -> {len(net)} edges")
