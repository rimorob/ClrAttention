#!/usr/bin/env python
"""Step 3.10 (descriptive arm, E10/E17): continued self-supervised training of scGPT whole-human on a case's own cells.
NO TRUTH FILE IS READ. No labels of any kind are used.

Run in ~/venvs/scfm (scgpt 0.2.4), python -I, one CUDA device, never device 0.
  python -I analysis/scfm_finetune_scgpt.py --cases hESC_500,hESC_1000,hHep_500,hHep_1000 --epochs 10 --out ~/scfm_models/step5_ft/scgpt_ckpt --device 1

scGPT 0.2.4 has no training entry point that runs without data-bank / wandb machinery; the package's own self-supervised fine-tuning
step is `scgpt.trainer.train` (gene-expression prediction GEP + masked-value cell-embedding prediction GEPC, both masked MSE of
`scgpt.loss.masked_mse_loss`, as in the package tutorials). This script performs exactly that step (same `model(...)` call, same loss
composition, same `scgpt.tokenizer.random_mask_value` masking, same gradient clipping at 1.0), with micro-batch accumulation so that
a 6 GB GPU can hold the 12-layer model at 500-1,400 tokens; see docs/preregistration_scfm.md E17 for the settings and deviations.
Inputs are those of the zero-shot extraction (analysis/scfm_extract_scgpt.py): arm A, all cells, the case genes covered by the scGPT
vocabulary, normalize_total 1e4 + log1p + 51-bin per-cell quantile binning (np.random.seed(20261021) before the binning), `<cls>` first
with value 0. Writes <out>/<case>.pt (fp32 state dict, deleted after the extraction) and <out>/<case>.json.
"""
import argparse, hashlib, json, os, sys, time
import numpy as np, pandas as pd, torch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
STEP2 = os.path.expanduser("~/scfm_models/step2/inputs")
SEED = 20261021
MASK_RATIO, MASK_VALUE, PAD_VALUE = 0.4, -1, -2   # scGPT tutorial defaults (mask_ratio 0.4, mask_value -1, pad_value -2)
LR, SCHED_GAMMA, CLIP = 1e-4, 0.9, 1.0            # checkpoint args.json lr 1e-4; tutorial StepLR(1, 0.9) per epoch; trainer.train clip 1.0


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def flat(sd):
    return {k: v.detach().float().cpu() for k, v in sd.items() if v.dtype.is_floating_point}


def rel_change(a, b):
    num = sum(float(((b[k] - a[k]) ** 2).sum()) for k in a) ** 0.5
    den = sum(float((a[k] ** 2).sum()) for k in a) ** 0.5
    return num / den


def step_loss(model, ids, vals, target, masked, crit):
    """The loss composition of scgpt.trainer.train with GEP=True, GEPC=True, CLS=False, ESC=False, DAR=False, explicit_zero_prob=False."""
    out = model(ids, vals, src_key_padding_mask=torch.zeros_like(ids, dtype=torch.bool), batch_labels=None, CLS=False, MVC=True, ECS=False)
    l_gep = crit(out["mlm_output"], target, masked)
    l_gepc = crit(out["mvc_output"], target, masked)
    return l_gep + l_gepc, l_gep, l_gepc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True)
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch", type=int, default=32, help="effective batch size (checkpoint args.json batch_size 32)")
    ap.add_argument("--micro-bs", type=int, default=0, help="cells per forward pass; 0 = auto by token count")
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--ncells-limit", type=int, default=0, help="TEST ONLY")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert a.device != 0
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
    dev = torch.device(f"cuda:{a.device}")
    os.makedirs(a.out, exist_ok=True)
    from scfm_smoke_scgpt import load_scgpt, SHA, D
    from scgpt.preprocess import binning
    from scgpt.tokenizer import random_mask_value
    from scgpt.loss import masked_mse_loss
    import scgpt
    assert all(sha256(f"{D}/{f}") == h for f, h in SHA.items())
    nb = 51

    for case in a.cases.split(","):
        base = f"{a.out}/{case}"
        if os.path.exists(base + ".json") and os.path.exists(base + ".pt"):
            print("skip (done)", case, flush=True)
            continue
        t0 = time.time()
        torch.manual_seed(SEED)
        np.random.seed(SEED)
        torch.cuda.manual_seed_all(SEED)
        model, vocab, args, info = load_scgpt(torch.device("cpu"))
        z = np.load(f"{STEP2}/{case}.npz", allow_pickle=True)
        gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t")
        assert list(gm.symbol) == list(z["genes"]) and case.startswith("h")
        sel = gm.in_scgpt.values & gm.scgpt_token.notna().values
        use = gm[sel].reset_index(drop=True)
        x = np.asarray(z["tpm_beeline"][sel], dtype=np.float64).T
        if a.ncells_limit:
            x = x[:a.ncells_limit]
        n, G = x.shape
        norm = np.log1p(x / x.sum(1, keepdims=True) * 1e4)
        np.random.seed(SEED)                                   # immediately before the binning (E4)
        binned = np.stack([binning(r, nb) for r in norm]).astype(np.float32)
        gids = torch.tensor([vocab["<cls>"]] + [vocab[t] for t in use.scgpt_token])
        vals_full = torch.cat([torch.zeros(n, 1), torch.tensor(binned)], 1)          # <cls> value 0
        ids_full = gids.unsqueeze(0).expand(n, -1)
        micro = a.micro_bs or (4 if G <= 700 else 2 if G <= 1000 else 1)
        micro = min(micro, a.batch)
        accum = a.batch // micro
        assert micro * accum == a.batch
        w0 = flat(model.state_dict())
        qk0 = {k: v for k, v in w0.items() if "in_proj_weight" in k}
        model.to(dev)
        opt = torch.optim.Adam(model.parameters(), lr=LR, eps=1e-8)        # tutorials: Adam, lr 1e-4 (eps 1e-8 without AMP)
        sched = torch.optim.lr_scheduler.StepLR(opt, 1, gamma=SCHED_GAMMA)
        crit = masked_mse_loss
        rng = np.random.default_rng(SEED)

        @torch.no_grad()
        def probe(nc=min(n, 64)):
            model.eval()
            np.random.seed(1)
            m = random_mask_value(vals_full[:nc], mask_ratio=MASK_RATIO, mask_value=MASK_VALUE, pad_value=PAD_VALUE)
            tot = 0.0
            for s in range(0, nc, micro):
                ids = ids_full[s:s + micro].to(dev)
                v, t = m[s:s + micro].to(dev), vals_full[s:s + micro].to(dev)
                l, _, _ = step_loss(model, ids, v, t, v.eq(MASK_VALUE), crit)
                tot += float(l) * ids.shape[0]
            return tot / nc

        loss_before = probe()
        torch.cuda.reset_peak_memory_stats(dev)
        t1 = time.time()
        steps, epoch_loss, epoch_gep, epoch_gepc = [], [], [], []
        np.random.seed(SEED)
        for ep in range(a.epochs):
            model.train()
            # the package tutorials re-draw the random masking of the whole training set at the start of every epoch
            masked_all = random_mask_value(vals_full, mask_ratio=MASK_RATIO, mask_value=MASK_VALUE, pad_value=PAD_VALUE)
            perm = rng.permutation(n)
            tl, tg, tc, nb_ = 0.0, 0.0, 0.0, 0
            for s in range(0, n, a.batch):
                idx = perm[s:s + a.batch]
                opt.zero_grad(set_to_none=True)
                bl = 0.0
                for m0 in range(0, len(idx), micro):
                    mi = idx[m0:m0 + micro]
                    ids = ids_full[mi].to(dev)
                    v, t = masked_all[mi].to(dev), vals_full[mi].to(dev)
                    l, lg, lc = step_loss(model, ids, v, t, v.eq(MASK_VALUE), crit)
                    (l * len(mi) / len(idx)).backward()
                    bl += float(l) * len(mi) / len(idx)
                    tg += float(lg) * len(mi)
                    tc += float(lc) * len(mi)
                torch.nn.utils.clip_grad_norm_(model.parameters(), CLIP, error_if_nonfinite=True)
                opt.step()
                steps.append(dict(epoch=ep, step=len(steps), loss=bl, lr=opt.param_groups[0]["lr"]))
                tl += bl * len(idx)
                nb_ += len(idx)
            sched.step()
            epoch_loss.append(tl / nb_)
            epoch_gep.append(tg / nb_)
            epoch_gepc.append(tc / nb_)
            print(f"{case} epoch {ep + 1}/{a.epochs} loss {epoch_loss[-1]:.4f} gep {epoch_gep[-1]:.4f} gepc {epoch_gepc[-1]:.4f} ({time.time() - t1:.0f}s)", flush=True)
        torch.cuda.synchronize(dev)
        train_s = time.time() - t1
        loss_after = probe()
        sd = {k: v.detach().cpu() for k, v in model.state_dict().items()}
        w1 = flat(sd)
        qk1 = {k: v for k, v in w1.items() if "in_proj_weight" in k}
        torch.save(sd, base + ".pt")
        R = dict(case=case, n_cells=int(n), n_genes=int(G), scgpt=scgpt.__version__, torch=torch.__version__, device=torch.cuda.get_device_name(dev),
                 seed=SEED, epochs=a.epochs, effective_batch=a.batch, micro_bs=micro, accum=accum, lr=LR, scheduler=f"StepLR(1, {SCHED_GAMMA}) per epoch",
                 optimizer="Adam eps 1e-8", grad_clip=CLIP, precision="fp32", mask_ratio=MASK_RATIO, mask_value=MASK_VALUE, pad_value=PAD_VALUE,
                 losses="GEP + GEPC (masked MSE); CLS, ECS, DAR, CCE off; explicit_zero_prob False", model_dropout=args["dropout"],
                 optimizer_steps=len(steps), epoch_mean_loss=epoch_loss, epoch_mean_gep=epoch_gep, epoch_mean_gepc=epoch_gepc,
                 probe_loss_before=loss_before, probe_loss_after=loss_after, rel_weight_change_all=rel_change(w0, w1),
                 rel_weight_change_inproj=rel_change(qk0, qk1), gpu_peak_alloc_MiB=torch.cuda.max_memory_allocated(dev) / 2**20,
                 train_s=train_s, total_s=time.time() - t0, state_sha256=sha256(base + ".pt"), state_bytes=os.path.getsize(base + ".pt"),
                 step_log=steps[::max(1, len(steps) // 200)], any_nonfinite_loss=bool(not np.isfinite([s["loss"] for s in steps]).all()))
        json.dump(R, open(base + ".json", "w"), indent=1, default=str)
        print(json.dumps({k: R[k] for k in ("case", "n_cells", "n_genes", "micro_bs", "optimizer_steps", "epoch_mean_loss", "probe_loss_before",
                                              "probe_loss_after", "rel_weight_change_all", "gpu_peak_alloc_MiB", "train_s", "any_nonfinite_loss")}), flush=True)
        del model, opt, w0, w1, sd
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
