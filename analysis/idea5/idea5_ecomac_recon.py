#!/usr/bin/env python
"""POST HOC, label-free (plan section 10): held-out masked-gene reconstruction of the M3D-trained refit models on the
EcoMAC arrays without an M3D counterpart (subset B, ComBat by platform, results/ecomac_nonoverlap/B/expression.tsv),
shared genes only (attention softmax renormalised over the shared genes; signs from the M3D training data).
No truth is read. Usage: python -I analysis/idea5/idea5_ecomac_recon.py --prep results/idea5/prep --models results/idea5/models --out results/idea5/ecomac_recon.csv"""
import argparse, json, csv, os, re
import numpy as np, torch
p = argparse.ArgumentParser()
p.add_argument("--prep", required=True); p.add_argument("--models", required=True); p.add_argument("--out", required=True)
p.add_argument("--expr", default="results/ecomac_nonoverlap/B/expression.tsv"); p.add_argument("--device", type=int, default=1)
a = p.parse_args()
dev = torch.device(f"cuda:{a.device}" if torch.cuda.is_available() else "cpu")
meta = json.load(open(f"{a.prep}/meta.json")); G, N = meta["G"], meta["N"]
genes = [l.strip() for l in open(f"{a.prep}/genes.txt")]
X = np.fromfile(f"{a.prep}/X.f64", dtype=np.float64).reshape(G, N, order="F")
Zf = np.fromfile(f"{a.prep}/Z_full.f64", dtype=np.float64).reshape(G, G, order="F").astype(np.float32)
with open(a.expr) as f:
    hdr = f.readline(); names = []; rows = []
    for line in f:
        s = line.rstrip("\n").split("\t"); names.append(s[0]); rows.append(np.array(s[1:], dtype=np.float64))
E = np.vstack(rows); pos = {g: i for i, g in enumerate(names)}
sh = [i for i, g in enumerate(genes) if g in pos]; ei = [pos[genes[i]] for i in sh]
Es = E[ei]; Es = (Es - Es.mean(1, keepdims=True)) / Es.std(1, ddof=1, keepdims=True)
Xm = X[sh]; Xm = (Xm - Xm.mean(1, keepdims=True)) / Xm.std(1, ddof=1, keepdims=True)
S = len(sh); print(f"shared genes {S} of {G}; EcoMAC arrays {E.shape[1]}", flush=True)
shi = torch.tensor(sh, device=dev)
Xe = torch.from_numpy(Es.astype(np.float32)).to(dev)
Sg = torch.where(torch.from_numpy((Xm @ Xm.T / (N - 1)).astype(np.float32)).to(dev) >= 0, 1.0, -1.0)
eye = torch.eye(S, dtype=torch.bool, device=dev)
Zs = (torch.from_numpy(Zf).to(dev)[shi][:, shi]) / meta["tau0"]["full"]
def read_D(tag):
    D = np.fromfile(f"{a.models}/{tag}.D.f32", dtype=np.float32).reshape(G, G)   # python wrote D.T: flat index j*G+i -> [j, i]
    return torch.from_numpy(D.T.copy()).to(dev)[shi][:, shi]
def recon(L):
    A = torch.softmax(L.masked_fill(eye, float("-inf")), 1)
    xh = (A * Sg) @ Xe
    mse = float(((xh - Xe) ** 2).mean()); return mse, 1 - mse / float((Xe ** 2).mean())
rows = []
mse, r2 = recon(Zs); rows.append(("CLR_step0_full", "clr", "", "", mse, r2))
for f in sorted(os.listdir(a.models)):
    m = re.match(r"^(clr|scratch)_d(\d+)_s(\d+)_full\.json$", f) or re.match(r"^ctrl_(scratch)_d(\d+)_s(\d+)\.json$", f)
    if not m: continue
    tag = f[:-5]; init = m.group(1); D = read_D(tag)
    L = D + (torch.from_numpy(np.fromfile(f"{a.models}/{tag}.gamma.f64", dtype=np.float64).astype(np.float32)).to(dev)[shi].exp()[:, None] * Zs if init == "clr" else 0)
    mse, r2 = recon(L); rows.append((tag, init, m.group(2), m.group(3), mse, r2))
with open(a.out, "w") as f:
    f.write("model,init,d,seed,mse,r2\n")
    for r in rows: f.write(",".join(str(x) for x in r) + "\n")
print(f"wrote {a.out}")
