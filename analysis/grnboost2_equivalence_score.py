#!/usr/bin/env python3
"""Score the equivalence gate (tools/grnboost2_equivalence.sh): per arm and seed, directed AUPRC ratio on the
BEELINE case's three primary truths (candidates = truth regulators x all genes; same candidate set as
beeline_baselines.R's score_truth; average precision from sklearn, ties handled by sklearn's step rule), plus
agreement of the edge weights between arms (Spearman) and the trees-per-target lines from the log.
    ~/venvs/clr-baselines/bin/python analysis/grnboost2_equivalence_score.py [case_dir_name=mDC_500 species=mouse ds=mDC]"""
import sys, glob, os
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score
from scipy.stats import spearmanr
case = sys.argv[1] if len(sys.argv) > 1 else "mDC_500"; sp = sys.argv[2] if len(sys.argv) > 2 else "mouse"
out = "results/grnboost2_equivalence"
genes = list(pd.read_csv(f"results/beeline_baselines/{case}/expression.tsv", sep="\t", usecols=[0]).iloc[:, 0])
up = {g.upper(): i for i, g in enumerate(genes)}; G = len(genes)
d = f"data/beeline/Networks/{sp}"
tr = {"celltype_ChIP": f"{d}/mDC-ChIP-seq-network.csv" if case.startswith("mDC") else None,
      "nonspecific_ChIP": f"{d}/Non-Specific-ChIP-seq-network.csv", "STRING": f"{d}/STRING-network.csv"}
def mat(f):
    e = pd.read_csv(f, sep="\t"); M = np.zeros((G, G))
    ia = e["TF"].map(lambda x: up.get(x.upper())); ib = e["target"].map(lambda x: up.get(x.upper()))
    ok = ia.notna() & ib.notna(); M[ia[ok].astype(int), ib[ok].astype(int)] = e["importance"][ok]; return M
res = []; arms = {}
for f in sorted(glob.glob(f"{out}/*_s[0-9].tsv")):
    arm, s = os.path.basename(f)[:-4].rsplit("_s", 1); arms[(arm, int(s))] = mat(f)
for tn, tf in tr.items():
    gt = pd.read_csv(tf); gi = gt["Gene1"].str.upper().map(up); gj = gt["Gene2"].str.upper().map(up)
    ok = gi.notna() & gj.notna() & (gi != gj); gi = gi[ok].astype(int).values; gj = gj[ok].astype(int).values
    regs = np.unique(gi); P = np.zeros((G, G), bool); P[gi, gj] = True
    C = np.zeros((G, G), bool); C[regs, :] = True; np.fill_diagonal(C, False)
    y = P[C].astype(int); dens = y.mean()
    for (arm, s), M in arms.items():
        res.append(dict(truth=tn, arm=arm, seed=s, ratio=average_precision_score(y, M[C]) / dens, n_edges=nz if (nz := int((M[C] > 0).sum())) else 0))
r = pd.DataFrame(res); pd.set_option("display.width", 200)
print(r.groupby(["truth", "arm"]).ratio.agg(["mean", "std", "min", "max"]).round(4))
print("\nnon-zero candidate edges (mean over truths/seeds):"); print(r.groupby("arm").n_edges.mean().round(0))
print("\nSpearman of edge weights between arms, same seed (all gene pairs):")
for s in (1, 2, 3):
    for a, b in (("legacy", "stock_pool"), ("stock_dask", "stock_pool"), ("legacy", "stock_dask")):
        if (a, s) in arms and (b, s) in arms:
            print(f"  seed {s} {a} vs {b}: rho = {spearmanr(arms[(a, s)].ravel(), arms[(b, s)].ravel())[0]:.4f}")
print("\nlog:"); print(open(f"{out}/log.txt").read())
r.to_csv(f"{out}/scores.csv", index=False)
