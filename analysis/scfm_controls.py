#!/usr/bin/env python
"""Step 3 controls that need no GPU: (1) the data-free gene-embedding cosine prior of scPRINT medium-v1.5 and (2) the truth-free
control-distance QC table. NO TRUTH FILE IS READ.

  python -I analysis/scfm_controls.py esm2cos --out ~/scfm_models/step3/esm2_cos
  python -I analysis/scfm_controls.py qc --step3 ~/scfm_models/step3 --out results/scfm_step3/control_distances.csv

esm2cos: for each of the 14 cases, the cosine matrix (G x G, float32, diagonal included) between the rows of the checkpoint's
`gene_encoder.embedding.weight` (the ESM2-derived 256-d gene embeddings stored in the checkpoint) for the arm-A genes (gene-map status ok).
qc: relative L1 / Spearman (off-diagonal) of the all-layer all-head mean R0 between the zero-shot run and each control
(permutation seeds 1001-1003, random weights seeds 1-3 in two versions), and between the controls, per case, plus the distance of
R0 to the uniform matrix 1/G. These are descriptive QC numbers; nothing is scored.
"""
import argparse, glob, hashlib, json, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = [f"{d}_{n}" for d in ("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L") for n in (500, 1000)]


def offdiag(a, b):
    from scipy.stats import rankdata
    G = a.shape[0]
    m = ~np.eye(G, dtype=bool)
    x, y = a[m].astype(np.float64), b[m].astype(np.float64)
    return float(np.abs(x - y).sum() / (0.5 * (x.sum() + y.sum()))), float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def esm2cos(out):
    import torch
    ck = torch.load(os.path.expanduser("~/scfm_models/scprint/medium-v1.5.ckpt"), map_location="cpu", weights_only=False)
    W = ck["state_dict"]["gene_encoder.embedding.weight"].float().numpy()
    vocab = list(ck["hyper_parameters"]["genes"])
    vidx = {g: i for i, g in enumerate(vocab)}
    os.makedirs(out, exist_ok=True)
    rows = []
    for c in CASES:
        gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{c}.tsv", sep="\t")
        u = gm[gm.status == "ok"].reset_index(drop=True)
        E = W[[vidx[e] for e in u.ensembl]]
        E = E / np.linalg.norm(E, axis=1, keepdims=True)
        C = (E @ E.T).astype(np.float32)
        p = f"{out}/{c}.esm2cos.npz"
        np.savez(p, cos=C, genes=np.array(u.symbol), ensembl=np.array(u.ensembl))
        rows.append(dict(case=c, G=len(u), file=os.path.basename(p), sha256=hashlib.sha256(open(p, "rb").read()).hexdigest()))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{out}/esm2cos_manifest.csv", index=False)


def qc(step3, out):
    rows = []
    for c in CASES:
        base = np.load(f"{step3}/scprint_medium_A/{c}.nets.npz")["R0"]
        G = base.shape[0]
        ctrl = {}
        for s in (1001, 1002, 1003):
            f = f"{step3}/ctrl_perm_s{s}/{c}.nets.npz"
            if os.path.exists(f):
                ctrl[f"perm{s}"] = np.load(f)["R0"]
        for v in ("rand_esm", "rand_emb"):
            for s in (1, 2, 3):
                f = f"{step3}/ctrl_{v}_s{s}/{c}.nets.npz"
                if os.path.exists(f):
                    ctrl[f"{v}{s}"] = np.load(f)["R0"]
        uni = np.full_like(base, 1.0 / G)
        r = dict(case=c, G=G)
        r["relL1_vs_uniform"], r["spearman_vs_uniform"] = float(np.abs(base - uni).sum() / (0.5 * (base.sum() + uni.sum()))), float("nan")
        for k, v in ctrl.items():
            r[f"{k}_relL1"], r[f"{k}_spearman"] = offdiag(base, v)
        rows.append(r)
        print(c, {k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()}, flush=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["esm2cos", "qc"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--step3", default=os.path.expanduser("~/scfm_models/step3"))
    a = ap.parse_args()
    esm2cos(os.path.expanduser(a.out)) if a.mode == "esm2cos" else qc(os.path.expanduser(a.step3), os.path.expanduser(a.out))
