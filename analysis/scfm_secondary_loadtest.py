#!/usr/bin/env python
"""Step 3.2(c) gate: load test of scPRINT large-v1 and small-v1 under scprint 2.2.0 (strict load_from_checkpoint), architecture
summary, vocabulary comparison with medium-v1.5, and a 16-cell forward on hESC_500 (arm A). No truth. Run in ~/venvs/scfm_alt.
  python -I analysis/scfm_secondary_loadtest.py --out ~/scfm_models/step3/_gates/secondary_loadtest.json --device 1
"""
import argparse, json, os, sys, time, traceback
import numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scfm_extract_scprint as ex

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--device", type=int, default=1)
a = ap.parse_args()
assert a.device != 0
dev = torch.device(f"cuda:{a.device}")
import scprint
R = {"scprint": scprint.__version__, "note": "depth token fed (row sum over the case genes) for depth_atinput models"}
med = torch.load(os.path.expanduser(ex.CKPTS["medium-v1.5"][0]), map_location="cpu", weights_only=False)["hyper_parameters"]["genes"]
use, x, cells, _ = ex.case_inputs("hESC_500", "A")
for name in ("large-v1", "small-v1"):
    r = {}
    ck = os.path.expanduser(ex.CKPTS[name][0])
    try:
        r["sha256_ok"] = ex.sha256(ck) == ex.CKPTS[name][1]
        t = time.time()
        model, info = ex.load_model(ck, "none", 0, dev)
        r["load_s"] = time.time() - t
        hp = torch.load(ck, map_location="cpu", weights_only=False)["hyper_parameters"]
        r.update(strict_load="ok", nlayers=model.nlayers, nheads=model.transformer.blocks[0].mixer.num_heads, d_model=model.d_model,
                 depth_atinput=bool(model.depth_atinput), n_prefix_tokens=int(model.cell_embs_count),
                 vocab_size=len(model.genes), vocab_identical_to_medium=list(model.genes) == list(med),
                 n_params=int(sum(p.numel() for p in model.parameters())), classes=list(model.classes),
                 hparams_other={k: str(hp.get(k)) for k in ("expr_emb_style", "normalization", "gene_pos_enc", "transformer", "cell_specific_blocks", "compress_class_dim", "d_hid", "mvc_decoder")})
        vidx = {g: i for i, g in enumerate(model.genes)}
        missing = [e for e in use.ensembl if e not in vidx]
        r["hESC_500_genes_missing_from_vocab"] = len(missing)
        keep = [i for i, e in enumerate(use.ensembl) if e in vidx]
        gp = torch.tensor([[vidx[use.ensembl[i]] for i in keep]], dtype=torch.long)
        t = time.time()
        As, _, mass, hid = ex.run(model, gp, torch.tensor(x[:16][:, keep], dtype=torch.float32), dev, 8)
        A = (As / 16).float().cpu().numpy()
        r.update(forward_16_cells_s=time.time() - t, forward_shape=list(A.shape), forward_nan=int(np.isnan(A).sum()),
                 forward_row_sum=[float(A.sum(-1).min()), float(A.sum(-1).max())],
                 gpu_peak_alloc_MiB=torch.cuda.max_memory_allocated(dev) / 2**20)
        del model, As
    except Exception as e:
        r["error"] = traceback.format_exc()[-2500:]
    R[name] = r
    torch.cuda.empty_cache()
json.dump(R, open(a.out, "w"), indent=1, default=str)
print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "hparams_other"} if isinstance(v, dict) else v) for k, v in R.items()}, indent=1, default=str))
