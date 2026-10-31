#!/usr/bin/env python
"""Step 1 smoke test of the secondary comparator: scGPT whole-human attention on hESC.

TECHNICAL SANITY ONLY. Never reads a truth network, scores nothing (docs/preregistration_scfm.md).

  python -I analysis/scfm_smoke_scgpt.py --device 1 --ncells 64 --out OUTDIR [--permute] [--tag NAME]

Attention readout, as in the preregistration: every layer and head, per cell, from Q and K of the
layer input (post-norm encoder, no flash-attn), softmax over gene keys only (= strip the <cls> token and
renormalise rows, exact and free of 0/0), averaged over cells. Input: the case genes (HGNC symbols that are in
the scGPT vocabulary), zeros included, TPM back-transformed from BEELINE's log2(TPM+1), per-cell
quantile binning into 51 bins (as scgpt.preprocess.Preprocessor(binning=51); zeros stay bin 0), <cls> first.
"""
import argparse, json, os, time, hashlib
import numpy as np, pandas as pd, torch

D = os.path.expanduser("~/scfm_models/scgpt_human")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPR = ROOT + "/results/beeline_baselines/hESC_1000/expression.tsv"
SHA = {"best_model.pt": "6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922",
       "vocab.json": "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f",
       "args.json": "c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def load_scgpt(device):
    import warnings
    warnings.filterwarnings("ignore")
    import torchtext
    torchtext.disable_torchtext_deprecation_warning()
    from scgpt.model import TransformerModel
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    from scgpt.utils import load_pretrained
    a = json.load(open(f"{D}/args.json"))
    vocab = GeneVocab.from_file(f"{D}/vocab.json")
    for s in (a["pad_token"], "<cls>", "<eoc>"):
        assert s in vocab
    vocab.set_default_index(vocab["<pad>"])
    m = TransformerModel(len(vocab), a["embsize"], a["nheads"], a["d_hid"], a["nlayers"], nlayers_cls=a["n_layers_cls"], n_cls=1,
                         vocab=vocab, dropout=a["dropout"], pad_token=a["pad_token"], pad_value=a["pad_value"], do_mvc=a["MVC"],
                         do_dab=False, use_batch_labels=False, domain_spec_batchnorm=False, n_input_bins=a["n_bins"],
                         input_emb_style=a["input_emb_style"], ecs_threshold=0.8, explicit_zero_prob=False,
                         use_fast_transformer=False, pre_norm=False)
    sd = torch.load(f"{D}/best_model.pt", map_location="cpu")
    before = {k: v.clone() for k, v in m.state_dict().items()}
    m = load_pretrained(m, sd, strict=False, verbose=False)
    after = m.state_dict()
    sd2 = {k.replace("Wqkv.", "in_proj_"): v for k, v in sd.items()}
    info = {"ckpt_keys_not_in_model": sorted(k for k in sd2 if k not in after),
            "model_keys_not_in_ckpt": sorted(k for k in after if k not in sd2),
            "n_params_unchanged_by_load": int(sum(torch.equal(before[k], after[k]) for k in after if k in sd2))}
    return m.eval().to(device), vocab, a, info


@torch.no_grad()
def run(model, gene_ids, cls_id, vals, device, bs):
    """vals: cells x G binned values. Returns mean (L, H, G, G) attention (<cls> stripped, rows renormalised), cls mass."""
    layers = model.transformer_encoder.layers
    nL, nH = len(layers), layers[0].self_attn.num_heads
    G = gene_ids.shape[0]
    A = torch.zeros(nL, nH, G, G, dtype=torch.float64, device=device)
    cmass = torch.zeros(nL, nH, dtype=torch.float64, device=device)
    n = vals.shape[0]
    chk = {}
    for s in range(0, n, bs):
        vb = vals[s:s + bs].to(device)
        B = vb.shape[0]
        ids = torch.cat([torch.full((B, 1), cls_id), gene_ids.unsqueeze(0).expand(B, -1)], 1).to(device)
        v = torch.cat([torch.zeros(B, 1, device=device), vb], 1)                      # <cls> value 0 as in scgpt.tokenizer
        x = model.encoder(ids) + model.value_encoder(v)                               # (B, 1+G, 512); no DSBN / batch labels
        for l, lay in enumerate(layers):
            sa = lay.self_attn
            W, b = sa.in_proj_weight, sa.in_proj_bias
            d = x.shape[-1]
            q = torch.nn.functional.linear(x, W[:d], b[:d]).view(B, -1, nH, d // nH)
            k = torch.nn.functional.linear(x, W[d:2 * d], b[d:2 * d]).view(B, -1, nH, d // nH)
            sc = torch.einsum("bthd,bshd->bhts", q, k * (d // nH) ** -0.5)[:, :, 1:, :]   # gene queries
            lse_all = torch.logsumexp(sc, -1)
            lse_gene = torch.logsumexp(sc[..., 1:], -1)
            cmass[l] += (1 - torch.exp(lse_gene - lse_all)).double().mean(2).sum(0)
            att = torch.softmax(sc[..., 1:], -1)
            A[l] += att.double().sum(0)
            if s == 0 and l == 0:                                                        # check against torch's own MHA weights
                _, w = sa(x, x, x, need_weights=True, average_attn_weights=False)
                full = torch.softmax(torch.einsum("bthd,bshd->bhts", q, k * (d // nH) ** -0.5), -1)
                chk["max_abs_diff_vs_torch_MHA_layer0"] = float((w - full).abs().max())
            x = lay(x)                                                                    # the layer's own forward (no mask: all tokens real)
            del sc, att
    A /= n
    cmass /= n
    return A.float().cpu().numpy(), cmass.float().cpu().numpy(), chk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--ncells", type=int, default=64)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bs", type=int, default=4)
    ap.add_argument("--permute", action="store_true")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
    dev = torch.device(f"cuda:{a.device}")
    t0 = time.time()
    M = {"tag": a.tag, "device": torch.cuda.get_device_name(dev), "torch": torch.__version__, "ncells": a.ncells, "seed": a.seed,
         "permute": a.permute, "sha256_ok": {f: sha256(f"{D}/{f}") == h for f, h in SHA.items()}}
    assert all(M["sha256_ok"].values())
    model, vocab, args, info = load_scgpt(dev)
    M["load_s"], M["load_info"] = time.time() - t0, info
    M["nlayers"], M["nheads"], M["vocab_size"] = len(model.transformer_encoder.layers), model.transformer_encoder.layers[0].self_attn.num_heads, len(vocab)
    case = pd.read_csv(EXPR, sep="\t", index_col=0)
    syms = list(case.index)
    inv = [s for s in syms if s in vocab]
    M["case_genes"], M["genes_in_vocab"], M["genes_dropped"] = len(syms), len(inv), len(syms) - len(inv)
    M["dropped_examples"] = [s for s in syms if s not in vocab][:15]
    gene_ids = torch.tensor([vocab[s] for s in inv])
    rng = np.random.default_rng(a.seed)
    cells = np.sort(rng.choice(case.shape[1], a.ncells, replace=False))
    lg = case.loc[inv].iloc[:, cells].values.T.astype(np.float64)
    tpm = 2.0 ** lg - 1.0
    if a.permute:
        prng = np.random.default_rng(1000 + a.seed)
        tpm = np.stack([prng.permutation(tpm[:, j]) for j in range(tpm.shape[1])], axis=1)
    from scgpt.preprocess import binning
    norm = np.log1p(tpm / tpm.sum(1, keepdims=True) * 1e4)                  # normalize_total 1e4 then log1p (binning is rank-based, so monotone)
    np.random.seed(a.seed)                                                  # scgpt's _digitize breaks bin-edge ties at random (side='both'): seed it
    binned = np.stack([binning(r, args["n_bins"]) for r in norm]).astype(np.float32)
    M["binned_max"], M["binned_frac_zero"] = float(binned.max()), float((binned == 0).mean())
    torch.cuda.reset_peak_memory_stats(dev)
    t1 = time.time()
    A, cmass, chk = run(model, gene_ids, vocab["<cls>"], torch.tensor(binned), dev, a.bs)
    torch.cuda.synchronize(dev)
    M["infer_s"], M["check"] = time.time() - t1, chk
    M["gpu_peak_alloc_MiB"], M["gpu_peak_reserved_MiB"] = torch.cuda.max_memory_allocated(dev) / 2**20, torch.cuda.max_memory_reserved(dev) / 2**20
    G = A.shape[-1]
    M["A_shape"], M["A_nan"], M["A_inf"] = list(A.shape), int(np.isnan(A).sum()), int(np.isinf(A).sum())
    M["A_row_sum_minmax"] = [float(A.sum(-1).min()), float(A.sum(-1).max())]
    M["cls_token_mass_layer_x_head"] = np.round(cmass, 4).tolist()
    M["A_rel_entropy_layer_x_head"] = np.round(-(A * np.log(np.clip(A, 1e-30, None))).sum(-1).mean(-1) / np.log(G), 3).tolist()
    np.savez_compressed(os.path.join(a.out, f"scgpt_attn_{a.tag}.npz"), A=A, genes=np.array(inv))
    M["total_s"] = time.time() - t0
    json.dump(M, open(os.path.join(a.out, f"scgpt_metrics_{a.tag}.json"), "w"), indent=1, default=str)
    print(json.dumps({k: M[k] for k in ("tag", "device", "genes_in_vocab", "genes_dropped", "A_shape", "A_nan", "infer_s", "gpu_peak_alloc_MiB", "total_s", "check")}, indent=1))


if __name__ == "__main__":
    main()
