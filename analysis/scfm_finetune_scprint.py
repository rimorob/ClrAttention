#!/usr/bin/env python
"""Step 3.10 (descriptive arm, E10/E17): continued self-supervised training of scPRINT medium-v1.5 on a case's own cells.
NO TRUTH FILE IS READ. No labels of any kind are used.

Route: scprint 2.2.0 in ~/venvs/scfm_alt (python -I), cwd = ~/scfm_models/lamin_instance, one CUDA device (never device 0).
  python -I analysis/scfm_finetune_scprint.py --cases hESC_500 --epochs 10 --out ~/scfm_models/step5_ft/scprint_ckpt --device 1

The package's own training code is used unchanged: `scPrint._full_training` / `_compute_loss` (losses), `configure_optimizers`
(AdamW, lr 1e-4 = hparams.lr, weight decay 0.01), `optimizer_step` (500-step linear warm-up), driven by a lightning Trainer.
The loss flags are set by the package's own `scprint.trainer.TrainingMode` callback with its DEFAULTS (denoising forward with
noise 0.6, bottleneck "generate" task, no masking, no full forward, class_embd_diss_scale 0.1, ZINB expression loss), except
`do_cls=False` (the classification loss needs cell-type / tissue / disease labels, which these data do not have and which would
be truths) and `lr_reduce_monitor=None` (the ReduceLROnPlateau schedule and the EarlyStopping callback of the package CLI watch a
validation loss; there is no validation split here, see docs/preregistration_scfm.md E17).
Inputs are exactly those of the zero-shot extraction (`case_inputs(case, "A")`): arm A, all cells, all case genes covered by the
checkpoint vocabulary, case-gene order. Writes <out>/<case>.pt (fp32 state dict, deleted after the extraction) and <out>/<case>.json.
"""
import argparse, copy, hashlib, json, os, sys, time
import numpy as np, torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scfm_extract_scprint as X

SEED = 20261021   # E10


class CellData(torch.utils.data.Dataset):
    def __init__(self, x, gpos):
        self.x = torch.tensor(x, dtype=torch.float32)
        self.g = gpos.reshape(-1)
        self.depth = self.x.sum(1)           # per-cell total over the case genes (the Collator rule used for the zero-shot depth tokens)

    def __len__(self):
        return self.x.shape[0]

    def __getitem__(self, i):
        return {"x": self.x[i], "genes": self.g, "depth": self.depth[i]}


def param_vec(model, keys=None):
    return {k: v.detach().float().cpu().clone() for k, v in model.state_dict().items() if (keys is None or keys(k)) and v.dtype.is_floating_point}


def rel_change(a, b):
    num = sum(float(((b[k] - a[k]) ** 2).sum()) for k in a) ** 0.5
    den = sum(float((a[k] ** 2).sum()) for k in a) ** 0.5
    return num / den


@torch.no_grad()
def probe_loss(model, ds, dev, n=64, bs=8):
    """Package training loss (same flags) on the first n cells, fixed torch seed (denoising noise is random)."""
    model.eval()
    torch.manual_seed(1)
    tot, cnt = 0.0, 0
    for s in range(0, n, bs):
        b = {k: v.to(dev) for k, v in torch.utils.data.default_collate([ds[i] for i in range(s, min(s + bs, n))]).items()}
        loss, _ = model._full_training(batch=b, do_denoise=model.do_denoise, noise=model.noise, do_next_tp=False, do_cce=False, cce_temp=model.cce_temp,
                                       do_ecs=False, do_mvc=False, do_adv_cls=False, do_adv_batch=False, do_cls=False,
                                       do_generate=model.do_generate, run_full_forward=model.run_full_forward, mask_ratio=model.mask_ratio)
        tot += float(loss) * b["x"].shape[0]
        cnt += b["x"].shape[0]
    return tot / cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True)
    ap.add_argument("--ckpt-name", default="medium-v1.5", choices=list(X.CKPTS))
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch", type=int, default=8, help="effective batch size (E17: 8, as the batch size of E10)")
    ap.add_argument("--micro-bs", default="auto", help="cells per forward pass; 'auto' = 8 if G <= 700, 4 if G <= 900, else 2 (6 GB GPU); "
                    "accumulation = batch / micro-bs, so the optimizer sees batch cells per step whatever the micro-batch")
    ap.add_argument("--clip", type=float, default=100.0, help="gradient_clip_val of the package config base.yml")
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--ncells-limit", type=int, default=0, help="TEST ONLY: first N cells")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert a.device != 0, "device 0 is the desktop GPU: never used"
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
    dev = torch.device(f"cuda:{a.device}")
    os.makedirs(a.out, exist_ok=True)
    import lightning as L, scprint
    from lightning.pytorch.callbacks import Callback
    from scprint.trainer import TrainingMode
    ckpt = os.path.expanduser(X.CKPTS[a.ckpt_name][0])
    assert X.sha256(ckpt) == X.CKPTS[a.ckpt_name][1], "checkpoint hash mismatch"

    class Rec(Callback):
        def __init__(self, accum):
            self.accum = accum
            self.steps, self.epoch_loss = [], []
            self._acc = []

        def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
            v = float(outputs["loss"] if isinstance(outputs, dict) else outputs) * self.accum   # lightning returns loss / accumulate_grad_batches
            self._acc.append(v)
            self.steps.append(dict(global_step=trainer.global_step, epoch=trainer.current_epoch, loss=v, lr=trainer.optimizers[0].param_groups[0]["lr"]))

        def on_train_epoch_end(self, trainer, pl_module):
            self.epoch_loss.append(float(np.mean(self._acc)))
            self._acc = []

    for case in a.cases.split(","):
        base = f"{a.out}/{case}"
        if os.path.exists(base + ".json") and os.path.exists(base + ".pt"):
            print("skip (done)", case, flush=True)
            continue
        t0 = time.time()
        L.seed_everything(SEED, workers=True)
        torch.manual_seed(SEED)
        np.random.seed(SEED)
        model, info = X.load_model(ckpt, "none", 0, torch.device("cpu"))
        vidx = {g: i for i, g in enumerate(model.genes)}
        use, x, cells, note = X.case_inputs(case, "A")
        if a.ncells_limit:
            x, cells = x[:a.ncells_limit], cells[:a.ncells_limit]
        gp = torch.tensor([[vidx[e] for e in use.ensembl]], dtype=torch.long)
        ds = CellData(x, gp)
        n, G = x.shape
        micro = (8 if G <= 700 else 4 if G <= 900 else 2) if a.micro_bs == "auto" else int(a.micro_bs)
        micro = min(micro, a.batch)
        accum = a.batch // micro
        assert micro * accum == a.batch
        gen = torch.Generator().manual_seed(SEED)
        dl = torch.utils.data.DataLoader(ds, batch_size=micro, shuffle=True, generator=gen, num_workers=0, drop_last=False)
        # the package's own TrainingMode defaults, then the two documented deviations
        tm = TrainingMode()
        tm.do_cls = False
        tm.lr_reduce_monitor = None
        tm.setup(None, model)
        flags = {k: getattr(model, k) for k in ("do_denoise", "noise", "do_generate", "do_cls", "do_cce", "do_ecs", "do_mvc", "do_adv_cls", "do_adv_batch",
                                                 "run_full_forward", "mask_ratio", "class_scale", "class_embd_diss_scale", "zinb_and_mse", "var_context_length",
                                                 "warmup_duration", "lr", "optim", "weight_decay", "lr_reduce_monitor", "fused_adam")}
        w0 = param_vec(model)
        qk0 = {k: v for k, v in w0.items() if "Wqkv" in k}
        model.to(dev)
        loss_before = probe_loss(model, ds, dev)
        rec = Rec(accum)
        tr = L.Trainer(max_epochs=a.epochs, accelerator="gpu", devices=[a.device], precision="32-true", gradient_clip_val=a.clip,
                       accumulate_grad_batches=accum, logger=False, enable_checkpointing=False, enable_progress_bar=False,
                       enable_model_summary=False, callbacks=[rec], num_sanity_val_steps=0, limit_val_batches=0, deterministic=False,
                       log_every_n_steps=10 ** 9)
        torch.cuda.reset_peak_memory_stats(dev)
        t1 = time.time()
        tr.fit(model, train_dataloaders=dl)
        torch.cuda.synchronize(dev)
        train_s = time.time() - t1
        model.to(dev).eval()
        loss_after = probe_loss(model, ds, dev)
        w1 = param_vec(model)
        qk1 = {k: v for k, v in w1.items() if "Wqkv" in k}
        sd = {k: v.detach().cpu() for k, v in model.state_dict().items()}
        torch.save(sd, base + ".pt")
        nsteps = int(tr.global_step)
        R = dict(case=case, n_cells=int(n), n_genes=int(G), scprint=scprint.__version__, torch=torch.__version__, device=torch.cuda.get_device_name(dev),
                 ckpt=ckpt, ckpt_sha256=X.CKPTS[a.ckpt_name][1], seed=SEED, epochs=a.epochs, micro_bs=micro, accum=accum, effective_batch=a.batch,
                 grad_clip=a.clip, precision="32-true", optimizer_steps=nsteps, train_flags=flags, final_lr=rec.steps[-1]["lr"],
                 epoch_mean_loss=rec.epoch_loss, probe_loss_before=loss_before, probe_loss_after=loss_after,
                 rel_weight_change_all=rel_change(w0, w1), rel_weight_change_Wqkv=rel_change(qk0, qk1),
                 gpu_peak_alloc_MiB=torch.cuda.max_memory_allocated(dev) / 2**20, train_s=train_s, total_s=time.time() - t0,
                 state_sha256=X.sha256(base + ".pt"), state_bytes=os.path.getsize(base + ".pt"), step_log=rec.steps[::max(1, len(rec.steps) // 200)],
                 any_nonfinite_loss=bool(not np.isfinite([s["loss"] for s in rec.steps]).all()))
        json.dump(R, open(base + ".json", "w"), indent=1, default=str)
        print(json.dumps({k: R[k] for k in ("case", "n_cells", "n_genes", "optimizer_steps", "epoch_mean_loss", "probe_loss_before", "probe_loss_after",
                                              "rel_weight_change_all", "rel_weight_change_Wqkv", "gpu_peak_alloc_MiB", "train_s", "any_nonfinite_loss")}), flush=True)
        del model, tr, ds, dl, w0, w1, sd
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
