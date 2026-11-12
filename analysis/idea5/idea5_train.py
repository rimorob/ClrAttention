#!/usr/bin/env python
"""EXPLORATORY (docs/plans/explore_idea5_plan.md). One attention layer over genes, trained self-supervised
(leave-one-gene-out masked reconstruction) on the M3D compendium. No truth is read here.

Model (plan section 2):  L_ij = a_i Z_ij / tau0 + D_ij  (j != i),  D = Q K' / sqrt(d) row-centred over j != i,
A = softmax_j L_ij (diagonal masked), xhat_i(s) = sum_j A_ij s_ij x_j(s), s_ij = sign(Pearson r_ij) on training samples,
loss = mean over genes and samples of (x - xhat)^2 on standardised expression.

Stages:  gate  (export step-0 CLR-initialised models for the equivalence gate)
         lr    (tuning seed 0, lr grid, both directions, per config; writes lr_selection.csv)
         final (seeds 1-5: both directions + refit on all samples + untrained scratch control)
Usage: python -I analysis/idea5/idea5_train.py --prep results/idea5/prep --out results/idea5/models --stage final [--device 1]
"""
import argparse, json, math, os, time, csv
import numpy as np
import torch

p = argparse.ArgumentParser()
p.add_argument("--prep", required=True); p.add_argument("--out", required=True)
p.add_argument("--stage", required=True, choices=["gate", "lr", "final"])
p.add_argument("--device", type=int, default=1)
p.add_argument("--dims", default="32,128"); p.add_argument("--inits", default="clr,scratch")
p.add_argument("--seeds", default="1,2,3,4,5"); p.add_argument("--lrs", default="0.001,0.003,0.01")
p.add_argument("--max_epochs", type=int, default=600); p.add_argument("--patience", type=int, default=30)
p.add_argument("--batch", type=int, default=64); p.add_argument("--warmup", type=int, default=10)
p.add_argument("--min_rel", type=float, default=1e-4)
a = p.parse_args()

dev = torch.device(f"cuda:{a.device}" if torch.cuda.is_available() and a.device >= 0 else "cpu")
meta = json.load(open(f"{a.prep}/meta.json")); G, N = meta["G"], meta["N"]
X_all = np.fromfile(f"{a.prep}/X.f64", dtype=np.float64).reshape(G, N, order="F")
split = np.array([int(r["half"]) for r in csv.DictReader(open(f"{a.prep}/split.csv"))])
Zs = {}
def get_Z(name):
    if name not in Zs:
        Zs[name] = torch.from_numpy(np.fromfile(f"{a.prep}/Z_{name}.f64", dtype=np.float64).reshape(G, G, order="F").astype(np.float32)).to(dev)
    return Zs[name]
eye = torch.eye(G, dtype=torch.bool, device=dev)
os.makedirs(a.out, exist_ok=True)


def standardise(tr_idx, other_idx):
    Xt = X_all[:, tr_idx]
    mu = Xt.mean(1, keepdims=True); sd = Xt.std(1, ddof=1, keepdims=True); sd[sd <= 0] = 1
    f = lambda M: torch.from_numpy(((M - mu) / sd).astype(np.float32)).to(dev)
    return f(Xt), (f(X_all[:, other_idx]) if len(other_idx) else None)


class Layer(torch.nn.Module):
    def __init__(self, init, d, seed, Zb):
        super().__init__()
        g = torch.Generator().manual_seed(seed)
        q = torch.randn(G, d, generator=g)
        k = torch.zeros(G, d) if init == "clr" else torch.randn(G, d, generator=g)
        self.q = torch.nn.Parameter(q.to(dev)); self.k = torch.nn.Parameter(k.to(dev))
        self.init, self.d, self.Zb = init, d, Zb
        self.gamma = torch.nn.Parameter(torch.zeros(G, device=dev)) if init == "clr" else None

    def D(self):
        D = (self.q @ self.k.T) / math.sqrt(self.d)
        return D - ((D.sum(1) - D.diagonal()) / (G - 1))[:, None]

    def forward(self, Xb, Sg):
        L = self.D()
        if self.init == "clr": L = L + torch.exp(self.gamma)[:, None] * self.Zb
        A = torch.softmax(L.masked_fill(eye, float("-inf")), 1)
        return (A * Sg) @ Xb


def evaluate(model, Xv, Sg):
    with torch.no_grad():
        mse = float(((model(Xv, Sg) - Xv) ** 2).mean()); base = float((Xv ** 2).mean())
    return mse, 1 - mse / base


def export(model, tag, tau0):
    with torch.no_grad():
        D = model.D(); D.fill_diagonal_(0)
        D.T.contiguous().cpu().numpy().astype(np.float32).tofile(f"{a.out}/{tag}.D.f32")
        if model.init == "clr": model.gamma.detach().cpu().numpy().astype(np.float64).tofile(f"{a.out}/{tag}.gamma.f64")


def run(init, d, seed, mode, lr, tag, epochs_fixed=None, export_step0_only=False):
    """mode: 'h1'/'h2' (train on that half, validate on the other), 'full' (all samples, fixed epochs)."""
    t0 = time.time(); torch.cuda.reset_peak_memory_stats() if dev.type == "cuda" else None
    if mode == "full": tr = np.arange(N); va = np.array([], dtype=int)
    else: h = int(mode[1]); tr = np.where(split == h)[0]; va = np.where(split != h)[0]
    zname = "full" if mode == "full" else mode
    Xt, Xv = standardise(tr, va)
    corr = (Xt @ Xt.T) / (Xt.shape[1] - 1)
    Sg = torch.where(corr >= 0, 1.0, -1.0).to(torch.float32); del corr
    Zb = (get_Z(zname) / meta["tau0"][zname]) if init == "clr" else None
    torch.manual_seed(seed)
    model = Layer(init, d, seed, Zb).to(dev)
    gen = torch.Generator().manual_seed(seed * 1000 + {"h1": 1, "h2": 2, "full": 3}[mode])
    params = [model.q, model.k] + ([model.gamma] if model.gamma is not None else [])
    opt = torch.optim.Adam(params, lr=lr, betas=(0.9, 0.999), eps=1e-8, weight_decay=0)
    curve = []; step = 0
    if Xv is not None:
        v0, r0 = evaluate(model, Xv, Sg); best, best_ep = v0, 0
        best_state = [x.detach().clone().cpu() for x in params]
    else: v0 = r0 = float("nan")
    ntr = Xt.shape[1]; n_ep = 0; stale = 0
    max_ep = epochs_fixed if epochs_fixed is not None else a.max_epochs
    if export_step0_only: max_ep = 0
    for ep in range(1, max_ep + 1):
        perm = torch.randperm(ntr, generator=gen).to(dev); tot = 0.0
        for b in range(0, ntr, a.batch):
            idx = perm[b:b + a.batch]; Xb = Xt[:, idx]
            for gp in opt.param_groups: gp["lr"] = lr * min(1.0, (step + 1) / a.warmup)
            loss = ((model(Xb, Sg) - Xb) ** 2).mean()
            opt.zero_grad(set_to_none=True); loss.backward(); opt.step(); step += 1
            tot += float(loss) * len(idx)
        n_ep = ep
        if Xv is not None:
            v, r2 = evaluate(model, Xv, Sg); curve.append((ep, tot / ntr, v, r2))
            if v < best * (1 - a.min_rel):
                best, best_ep, stale = v, ep, 0
                best_state = [x.detach().clone().cpu() for x in params]
            else: stale += 1
            if stale >= a.patience: break
        else: curve.append((ep, tot / ntr, float("nan"), float("nan")))
    if Xv is not None:
        with torch.no_grad():
            for x, s in zip(params, best_state): x.copy_(s.to(dev))
        vb, rb = evaluate(model, Xv, Sg)
    else: best_ep, vb, rb = n_ep, float("nan"), float("nan")
    export(model, tag, meta["tau0"][zname])
    with open(f"{a.out}/{tag}_curve.csv", "w") as f:
        f.write("epoch,train_mse,val_mse,val_r2\n")
        if Xv is not None: f.write(f"0,,{v0:.8g},{r0:.8g}\n")
        for r in curve: f.write(",".join(f"{x:.8g}" if i else str(int(x)) for i, x in enumerate(r)) + "\n")
    info = dict(tag=tag, init=init, d=d, seed=seed, mode=mode, lr=lr, best_epoch=int(best_ep), epochs_run=int(n_ep), steps=int(step),
                val_mse_step0=v0, val_r2_step0=r0, val_mse_best=vb, val_r2_best=rb, seconds=round(time.time() - t0, 1),
                peak_gpu_mb=round(torch.cuda.max_memory_allocated() / 2**20) if dev.type == "cuda" else None)
    json.dump(info, open(f"{a.out}/{tag}.json", "w"))
    del model, Xt, Xv, Sg, opt
    if dev.type == "cuda": torch.cuda.empty_cache()
    return info


def done(tag): return os.path.exists(f"{a.out}/{tag}.json")
def log(*x): print(time.strftime("[%H:%M:%S]"), *x, flush=True)
dims = [int(x) for x in a.dims.split(",")]; inits = a.inits.split(","); seeds = [int(x) for x in a.seeds.split(",")]
log(f"device {dev}, G={G}, N={N}, stage {a.stage}, torch {torch.__version__}")

if a.stage == "gate":
    for m in ["h1", "h2", "full"]:
        run("clr", 32, 1, m, 0.001, f"gate_clr_d32_{m}", export_step0_only=True)
    log("gate exports written")
elif a.stage == "lr":
    lrs = [float(x) for x in a.lrs.split(",")]; rows = []
    for init in inits:
        for d in dims:
            for lr in lrs:
                res = []
                for m in ["h1", "h2"]:
                    tag = f"lr_{init}_d{d}_lr{lr:g}_{m}"
                    if not done(tag): info = run(init, d, 0, m, lr, tag); log(tag, f"best_ep {info['best_epoch']} val {info['val_mse_best']:.5f} step0 {info['val_mse_step0']:.5f} {info['seconds']}s")
                    res.append(json.load(open(f"{a.out}/{tag}.json")))
                rows.append((init, d, lr, np.mean([r["val_mse_best"] for r in res]), res[0]["best_epoch"], res[1]["best_epoch"]))
    with open(f"{a.out}/lr_selection.csv", "w") as f:
        f.write("init,d,lr,mean_best_val_mse,best_epoch_h1,best_epoch_h2,selected\n")
        for init in inits:
            for d in dims:
                cand = [r for r in rows if r[0] == init and r[1] == d]
                sel = min(cand, key=lambda r: (round(r[3], 12), r[2]))[2]
                for r in cand: f.write(f"{r[0]},{r[1]},{r[2]:g},{r[3]:.8g},{r[4]},{r[5]},{int(r[2] == sel)}\n")
    log("lr selection written")
else:
    sel = {}
    for r in csv.DictReader(open(f"{a.out}/lr_selection.csv")):
        if r["selected"] == "1": sel[(r["init"], int(r["d"]))] = float(r["lr"])
    for seed in seeds:
        for init in inits:
            for d in dims:
                lr = sel[(init, d)]; ep = []
                for m in ["h1", "h2"]:
                    tag = f"{init}_d{d}_s{seed}_{m}"
                    if not done(tag): info = run(init, d, seed, m, lr, tag); log(tag, f"lr {lr:g} best_ep {info['best_epoch']} run {info['epochs_run']} R2 {info['val_r2_step0']:.4f}->{info['val_r2_best']:.4f} {info['seconds']}s {info['peak_gpu_mb']}MB")
                    ep.append(json.load(open(f"{a.out}/{tag}.json"))["best_epoch"])
                E = int(math.floor(np.mean(ep) + 0.5)); tag = f"{init}_d{d}_s{seed}_full"
                if not done(tag): info = run(init, d, seed, "full", lr, tag, epochs_fixed=E); log(tag, f"refit E={E} {info['seconds']}s")
                if init == "scratch":    # untrained control: the step-0 layer of this seed (data-independent)
                    tag = f"ctrl_scratch_d{d}_s{seed}"
                    if not done(tag): run(init, d, seed, "full", lr, tag, export_step0_only=True)
    log("final stage done")
