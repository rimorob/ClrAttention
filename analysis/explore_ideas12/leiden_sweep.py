#!/usr/bin/env python3
"""EXPLORATORY (docs/plans/explore_ideas12_plan.md, idea 1). Leiden sweep over a log-resolution grid on half-data attention graphs,
the rewired-operator null, nVI reproducibility weights and consensus co-membership counts.
Usage: ~/venvs/scfm/bin/python -I leiden_sweep.py --graphs DIR --out DIR --G N [--splits 6] [--seeds 5] [--workers 8]
Inputs  DIR/{full,s<k>h<1|2>}.tsv : i<TAB>j<TAB>w, 0-based, i<j, undirected weighted graph W=(A+A')/2.
Outputs OUT/counts_{obs,null,full}_g<idx>.u8 (GxG uint8, number of runs co-clustered), OUT/weights.csv, OUT/memberships.npz,
        OUT/summary.csv"""
import argparse, os, sys, time, random, itertools
import numpy as np, igraph as ig, leidenalg as la
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument("--graphs"); ap.add_argument("--out"); ap.add_argument("--G", type=int)
ap.add_argument("--splits", type=int, default=6); ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--workers", type=int, default=8); ap.add_argument("--kmin", type=int, default=-6); ap.add_argument("--kmax", type=int, default=6)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
GAMMAS = [2.0 ** k for k in range(a.kmin, a.kmax + 1)]; NG = len(GAMMAS); G = a.G
SEED0 = 20261200; REWIRE0 = 20261300
t00 = time.time()
def log(*x): print("[%6.0fs]" % (time.time() - t00), *x, flush=True)

def load(name):
    d = np.loadtxt(os.path.join(a.graphs, name + ".tsv"), dtype=np.float64, ndmin=2)
    g = ig.Graph(n=G, edges=list(zip(d[:, 0].astype(int).tolist(), d[:, 1].astype(int).tolist()))); g.es["weight"] = d[:, 2].tolist(); return g

GR = {}                                            # (kind, name) -> graph, shared by fork
GR[("full", "full")] = load("full")
for k in range(1, a.splits + 1):
    for h in (1, 2):
        nm = "s%dh%d" % (k, h); g = load(nm); GR[("obs", nm)] = g
        gn = g.copy(); random.seed(REWIRE0 + 100 * k + 10 * h)                    # degree-preserving edge swaps, weights permuted across edges
        gn.rewire(n=10 * gn.ecount(), allowed_edge_types="simple")
        w = np.array(g.es["weight"]); np.random.default_rng(REWIRE0 + 100 * k + 10 * h).shuffle(w); gn.es["weight"] = w.tolist(); GR[("null", nm)] = gn
log("graphs loaded/rewired:", len(GR), "edges per half graph ~", GR[("obs", "s1h1")].ecount())

def task(t):
    kind, nm, gi, s = t
    p = la.find_partition(GR[(kind, nm)], la.RBConfigurationVertexPartition, weights="weight", resolution_parameter=GAMMAS[gi], n_iterations=-1, seed=SEED0 + s)
    return t, np.asarray(p.membership, dtype=np.int32)

tasks = [("full", "full", gi, s) for gi in range(NG) for s in range(1, a.seeds + 1)]
for kind in ("obs", "null"):
    for k in range(1, a.splits + 1):
        for h in (1, 2):
            tasks += [(kind, "s%dh%d" % (k, h), gi, s) for gi in range(NG) for s in range(1, a.seeds + 1)]
log("tasks:", len(tasks))
res = {}
with Pool(a.workers) as pool:
    for n, (t, m) in enumerate(pool.imap_unordered(task, tasks, chunksize=2)):
        res[t] = m
        if (n + 1) % 100 == 0: log("done", n + 1, "of", len(tasks))
np.savez_compressed(os.path.join(a.out, "memberships.npz"), **{"%s|%s|%d|%d" % t: m for t, m in res.items()})

def nvi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True)
    kx, ky = xi.max() + 1, yi.max() + 1
    c = np.bincount(xi.astype(np.int64) * ky + yi, minlength=kx * ky).reshape(kx, ky).astype(np.float64) / len(x)
    px, py = c.sum(1), c.sum(0); nz = c > 0
    vi = -(c[nz] * (np.log(c[nz] / np.repeat(px, ky).reshape(kx, ky)[nz]) + np.log(c[nz] / np.tile(py, kx).reshape(kx, ky)[nz]))).sum()
    return vi / np.log(len(x))

rows = []
for gi in range(NG):
    r = {"gamma_index": gi, "gamma": GAMMAS[gi]}
    for kind in ("obs", "null"):
        v = [nvi(res[(kind, "s%dh1" % k, gi, s)], res[(kind, "s%dh2" % k, gi, s)]) for k in range(1, a.splits + 1) for s in range(1, a.seeds + 1)]
        r["nvi_" + kind] = float(np.mean(v)); r["nvi_%s_sd" % kind] = float(np.std(v, ddof=1))
        ncl = [res[(kind, "s%dh%d" % (k, h), gi, s)].max() + 1 for k in range(1, a.splits + 1) for h in (1, 2) for s in range(1, a.seeds + 1)]
        big = [np.bincount(res[(kind, "s%dh%d" % (k, h), gi, s)]).max() for k in range(1, a.splits + 1) for h in (1, 2) for s in range(1, a.seeds + 1)]
        r["n_comm_" + kind] = float(np.mean(ncl)); r["largest_" + kind] = float(np.mean(big))
    ncl = [res[("full", "full", gi, s)].max() + 1 for s in range(1, a.seeds + 1)]; big = [np.bincount(res[("full", "full", gi, s)]).max() for s in range(1, a.seeds + 1)]
    r["n_comm_full"] = float(np.mean(ncl)); r["largest_full"] = float(np.mean(big))
    # singleton share, mean over obs runs
    r["singleton_frac_obs"] = float(np.mean([(np.bincount(res[("obs", "s%dh%d" % (k, h), gi, s)]) == 1).sum() / G for k in range(1, a.splits + 1) for h in (1, 2) for s in range(1, a.seeds + 1)]))
    rows.append(r)
import csv
exc = np.maximum(0.0, np.array([r["nvi_null"] - r["nvi_obs"] for r in rows]))
wR = exc / exc.sum() if exc.sum() > 0 else np.full(NG, 1.0 / NG)
gm = np.array(GAMMAS)
for i, r in enumerate(rows): r["excess"] = float(exc[i]); r["w_R"] = float(wR[i]); r["w_U"] = 1.0 / NG; r["w_G"] = float(gm[i] / gm.sum())
with open(os.path.join(a.out, "weights.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
log("weights written; sum excess = %.4f%s" % (exc.sum(), "" if exc.sum() > 0 else " (all zero: uniform fallback)"))

for kind, nm_list, nruns in (("obs", ["s%dh%d" % (k, h) for k in range(1, a.splits + 1) for h in (1, 2)], a.splits * 2 * a.seeds),
                             ("null", ["s%dh%d" % (k, h) for k in range(1, a.splits + 1) for h in (1, 2)], a.splits * 2 * a.seeds),
                             ("full", ["full"], a.seeds)):
    for gi in range(NG):
        cnt = np.zeros((G, G), dtype=np.uint8)
        for nm in nm_list:
            for s in range(1, a.seeds + 1):
                m = res[(kind, nm, gi, s)]; cnt += (m[:, None] == m[None, :])
        assert cnt.max() <= nruns
        cnt.tofile(os.path.join(a.out, "counts_%s_g%02d.u8" % (kind, gi)))
    log("counts written:", kind, "(runs per gamma %d)" % nruns)
log("done")
