#!/usr/bin/env python
"""Step 3 production extraction of scGPT whole-human attention and cosine readouts. NO TRUTH FILE IS READ.

Run in ~/venvs/scfm (scgpt 0.2.4), python -I, one CUDA device, never device 0.
  python -I analysis/scfm_extract_scgpt.py --cases hESC_500,hESC_1000,hHep_500,hHep_1000 --out ~/scfm_models/step3/scgpt_A --device 1
  python -I analysis/scfm_extract_scgpt.py --reproduce-step1 ~/scfm_models/smoke_out --out ... (unit check, hESC_1000, 64 cells)
  python -I analysis/scfm_extract_scgpt.py --cases hESC_1000 --binning-seed 1 --tag seed1 --out ... (stability rerun)

Per case (arm A only: scGPT's preprocessing is invariant to a per-cell factor, E6): all cells; genes = the case genes covered by
the scGPT vocabulary (gene-map column in_scgpt, token scgpt_token), zeros included; input 2^y - 1 of the BEELINE matrix;
normalize_total 1e4 per cell, log1p, per-cell quantile binning into 51 bins over the case genes (zeros stay bin 0); `<cls>` first
with value 0. np.random.seed(binning-seed) is set immediately before the binning of each case (E4; default 20261021).
Writes <case>.json, <case>.nets.npz with
  attn_last  (heads, G, G)  last-layer attention per head, cell-mean, <cls> removed and rows renormalised (softmax over gene keys)
  attn_last_mean  (G, G)    head mean of attn_last: the preregistered scGPT attention readout
  tok_cos / tok_cos_ln      cosine of the gene-token embeddings (raw embedding weights / after the encoder LayerNorm); data-free
  hid_cos                   final-layer hidden-state cosine between gene tokens, averaged over cells
and in the json the <cls> mass and relative row entropy for every layer x head.
"""
import argparse, hashlib, json, os, sys, time
import numpy as np, pandas as pd, torch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
STEP2 = os.path.expanduser("~/scfm_models/step2/inputs")
PROD_SEED = 20261021


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


@torch.no_grad()
def run(model, gene_ids, cls_id, vals, device, bs):
    layers = model.transformer_encoder.layers
    nL, nH = len(layers), layers[0].self_attn.num_heads
    G = gene_ids.shape[0]
    A = torch.zeros(nL, nH, G, G, dtype=torch.float64, device=device)
    cmass = torch.zeros(nL, nH, dtype=torch.float64, device=device)
    H = torch.zeros(G, G, dtype=torch.float64, device=device)
    n = vals.shape[0]
    chk = {}
    for s in range(0, n, bs):
        vb = vals[s:s + bs].to(device)
        B = vb.shape[0]
        ids = torch.cat([torch.full((B, 1), cls_id), gene_ids.unsqueeze(0).expand(B, -1)], 1).to(device)
        v = torch.cat([torch.zeros(B, 1, device=device), vb], 1)
        x = model.encoder(ids) + model.value_encoder(v)
        for l, lay in enumerate(layers):
            sa = lay.self_attn
            W, b = sa.in_proj_weight, sa.in_proj_bias
            d = x.shape[-1]
            q = torch.nn.functional.linear(x, W[:d], b[:d]).view(B, -1, nH, d // nH)
            k = torch.nn.functional.linear(x, W[d:2 * d], b[d:2 * d]).view(B, -1, nH, d // nH)
            sc = torch.einsum("bthd,bshd->bhts", q, k * (d // nH) ** -0.5)[:, :, 1:, :]
            lse_all = torch.logsumexp(sc, -1)
            lse_gene = torch.logsumexp(sc[..., 1:], -1)
            cmass[l] += (1 - torch.exp(lse_gene - lse_all)).double().mean(2).sum(0)
            att = torch.softmax(sc[..., 1:], -1)
            A[l] += att.double().sum(0)
            if s == 0 and l == 0:
                _, w = sa(x, x, x, need_weights=True, average_attn_weights=False)
                full = torch.softmax(torch.einsum("bthd,bshd->bhts", q, k * (d // nH) ** -0.5), -1)
                chk["max_abs_diff_vs_torch_MHA_layer0"] = float((w - full).abs().max())
            x = lay(x)
            del sc, att
        hn = torch.nn.functional.normalize(x[:, 1:].float(), dim=-1)           # final-layer hidden states of the gene tokens
        H += torch.einsum("bgd,bhd->gh", hn.double(), hn.double())
    return A / n, cmass / n, H / n, chk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="hESC_500,hESC_1000,hHep_500,hHep_1000")
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--bs", type=int, default=4)
    ap.add_argument("--binning-seed", type=int, default=PROD_SEED)
    ap.add_argument("--tag", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ncells-limit", type=int, default=0, help="TEST ONLY")
    ap.add_argument("--load-state-dir", default="", help="E17 fine-tuning arm: per case, replace the model weights by <dir>/<case>.pt (fp32 state dict "
                    "written by analysis/scfm_finetune_scgpt.py) before the unchanged readout; empty = zero-shot (default, no effect)")
    ap.add_argument("--reproduce-step1", default="")
    a = ap.parse_args()
    assert a.device != 0
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
    dev = torch.device(f"cuda:{a.device}")
    os.makedirs(a.out, exist_ok=True)
    from scfm_smoke_scgpt import load_scgpt, SHA, D
    from scgpt.preprocess import binning
    R = {"torch": torch.__version__, "device": torch.cuda.get_device_name(dev), "args": vars(a),
         "sha256_ok": {f: sha256(f"{D}/{f}") == h for f, h in SHA.items()}}
    assert all(R["sha256_ok"].values())
    import scgpt
    R["scgpt"] = scgpt.__version__
    model, vocab, args, info = load_scgpt(dev)
    R["load_info"] = info
    nb = args["n_bins"]

    def bin_cells(tpm, seed):
        norm = np.log1p(tpm / tpm.sum(1, keepdims=True) * 1e4)
        np.random.seed(seed)                                # immediately before binning (E4)
        return np.stack([binning(r, nb) for r in norm]).astype(np.float32)

    if a.reproduce_step1:
        sm = os.path.expanduser(a.reproduce_step1)
        z = np.load(f"{sm}/scgpt_attn_run1.npz", allow_pickle=True)
        Aref, genes = z["A"], list(z["genes"])
        case = pd.read_csv(ROOT + "/results/beeline_baselines/hESC_1000/expression.tsv", sep="\t", index_col=0)
        cells = np.sort(np.random.default_rng(1).choice(case.shape[1], 64, replace=False))
        tpm = 2.0 ** case.loc[genes].iloc[:, cells].values.T.astype(np.float64) - 1.0
        binned = bin_cells(tpm, 1)
        A, cm, H, chk = run(model, torch.tensor([vocab[g] for g in genes]), vocab["<cls>"], torch.tensor(binned), dev, 4)
        A = A.float().cpu().numpy()
        d = np.abs(A - Aref)
        R["reproduce_step1"] = dict(shape=list(A.shape), max_abs_diff=float(d.max()), bitwise_identical=bool(np.array_equal(A, Aref)),
                                    rel_L1=float(d.sum() / np.abs(Aref).sum()), check=chk)
        json.dump(R, open(f"{a.out}/reproduce_step1_scgpt.json", "w"), indent=1, default=str)
        print(json.dumps(R["reproduce_step1"], indent=1))
        return

    tok = model.encoder.embedding.weight.detach()
    manifest = []
    for case in a.cases.split(","):
        tc = time.time()
        ft = {}
        if a.load_state_dir:
            sp = os.path.join(os.path.expanduser(a.load_state_dir), case + ".pt")
            ft = dict(ft_state_file=sp, ft_state_sha256=sha256(sp))
            model.load_state_dict(torch.load(sp, map_location="cpu"), strict=True)
            model.eval().to(dev)
            tok = model.encoder.embedding.weight.detach()
        z = np.load(f"{STEP2}/{case}.npz", allow_pickle=True)
        gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t")
        assert list(gm.symbol) == list(z["genes"]) and case.startswith("h")
        sel = gm.in_scgpt.values & gm.scgpt_token.notna().values
        use = gm[sel].reset_index(drop=True)
        assert all(t in vocab for t in use.scgpt_token)
        x = np.asarray(z["tpm_beeline"][sel], dtype=np.float64).T
        if a.ncells_limit:
            x = x[:a.ncells_limit]
        n, G = x.shape
        binned = bin_cells(x, a.binning_seed)
        ids = torch.tensor([vocab[t] for t in use.scgpt_token])
        torch.cuda.reset_peak_memory_stats(dev)
        t1 = time.time()
        A, cm, H, chk = run(model, ids, vocab["<cls>"], torch.tensor(binned), dev, a.bs)
        torch.cuda.synchronize(dev)
        infer_s = time.time() - t1
        A, cm = A.float().cpu().numpy(), cm.float().cpu().numpy()
        ent = -(A * np.log(np.clip(A, 1e-30, None))).sum(-1).mean(-1) / np.log(G)
        al = A[-1]
        with torch.no_grad():
            e_raw = torch.nn.functional.normalize(tok[ids.to(dev)].float(), dim=-1)
            e_ln = torch.nn.functional.normalize(model.encoder(ids.unsqueeze(0).to(dev))[0].float(), dim=-1)
            tok_cos, tok_cos_ln = (e_raw @ e_raw.T).cpu().numpy(), (e_ln @ e_ln.T).cpu().numpy()
        M = dict(case=case, n_cells=int(n), n_genes=int(G), binning_seed=a.binning_seed, tag=a.tag, infer_s=infer_s, **ft,
                 gpu_peak_alloc_MiB=torch.cuda.max_memory_allocated(dev) / 2**20, check=chk,
                 attn_nan=int(np.isnan(A).sum()), attn_inf=int(np.isinf(A).sum()), row_sum_min=float(A.sum(-1).min()), row_sum_max=float(A.sum(-1).max()),
                 cls_token_mass=np.round(cm, 5).tolist(), rel_entropy=np.round(ent, 4).tolist(),
                 hid_cos_nan=int(np.isnan(H.cpu().numpy()).sum()), binned_max=float(binned.max()), binned_frac_zero=float((binned == 0).mean()),
                 attn_all_layers_head_mean_sha256=hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest())
        base = f"{a.out}/{case}{('_' + a.tag) if a.tag else ''}"
        np.savez(base + ".nets.npz", attn_last=al.astype(np.float32), attn_last_mean=al.mean(0).astype(np.float32),
                 attn_all_mean=A.mean((0, 1)).astype(np.float32), tok_cos=tok_cos.astype(np.float32), tok_cos_ln=tok_cos_ln.astype(np.float32),
                 hid_cos=H.float().cpu().numpy(), genes=np.array(use.symbol), tokens=np.array(use.scgpt_token))
        M["nets_file_sha256"] = sha256(base + ".nets.npz")
        M["total_s"] = time.time() - tc
        json.dump(M, open(base + ".json", "w"), indent=1, default=str)
        manifest.append(dict(case=case, n=n, G=G, infer_s=round(infer_s, 1), nan=M["attn_nan"], cls_mass_last_layer_max=float(cm[-1].max())))
        print(json.dumps(manifest[-1]), flush=True)
        del A
        torch.cuda.empty_cache()
    R["manifest"], R["total_s"] = manifest, None
    json.dump(R, open(f"{a.out}/run{('_' + a.tag) if a.tag else ''}.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
