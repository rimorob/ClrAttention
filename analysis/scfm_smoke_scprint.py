#!/usr/bin/env python
"""Step 1 smoke test of the learned-attention comparator: scPRINT medium-v1.5 attention on hESC.

TECHNICAL SANITY ONLY. This script never reads a truth network and scores nothing
(docs/preregistration_scfm.md forbids scoring before the preregistration is committed).

Run (see docs/scfm_smoke_test.md for the exact environment and the systemd-run wrapper):
  python -I analysis/scfm_smoke_scprint.py --device 1 --ncells 64 --out OUTDIR [--permute] [--tag NAME]

What it does
  1. Loads the scPRINT medium-v1.5 checkpoint with a shim (see load_scprint): the installed scprint 1.0.13
     does not match the checkpoint's architecture (decoder takes depth as a concatenated feature; no
     depth token at the input), so only the encoder side is used and the depth token is not fed.
  2. Maps the case gene symbols to Ensembl IDs (Ensembl release 110 table via bionty) and to the
     checkpoint vocabulary; reports dropped genes.
  3. Runs the encoder in fp32 on N cells and computes the attention of every layer and head per cell
     from the captured Q and K (softmax over all tokens), removes the 9 special (cell/class) tokens
     (rows and columns), renormalises each row to sum 1, averages over cells.
  4. Writes the (layers x heads x G x G) mean attention, metrics JSON and the gene mapping.
"""
import argparse, json, os, time, inspect, hashlib
import numpy as np, pandas as pd, torch

CKPT = os.path.expanduser("~/scfm_models/scprint/medium-v1.5.ckpt")
CKPT_SHA256 = "96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPR = ROOT + "/results/beeline_baselines/hESC_1000/expression.tsv"        # case matrix (genes x cells), log2(TPM+1); gene list only is used from the case
RAW = ROOT + "/data/beeline/BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv"
ENS_TSV = os.path.expanduser("~/scfm_models/lamin_instance/ensembl110_human_genes.tsv")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def load_scprint(device):
    """Build scPrint from the checkpoint without lightning's CLI instantiator and without the depth parts.

    Shims (all logged in docs/scfm_smoke_test.md):
      - drop checkpoint hparams the installed class does not accept (compress_class_dim, cell_specific_blocks,
        depth_atinput, label_counts); the first three are all falsy/None, label_counts duplicates classes;
      - transformer='normal' (no flash-attn on Pascal), precpt_gene_emb=None (the weights are in the state dict);
      - drop expr_decoder.fc.0.weight (ckpt: 256x257 = depth concatenated at the decoder; installed: 256x256);
        the decoder is not used here;
      - the checkpoint has no depth_encoder weights: depth is never fed (full_depth=None), so the input prefix is
        the 9 class tokens (base + 8 classes), exactly as in the checkpoint's training.
    """
    from scprint import scPrint
    from scprint.model.flash_attn.flashformer import FlashTransformerEncoder
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    hp = {k: v for k, v in ck["hyper_parameters"].items() if k != "_instantiator"}
    ok = set(inspect.signature(scPrint.__init__).parameters) | set(inspect.signature(FlashTransformerEncoder.__init__).parameters)
    dropped = sorted(k for k in hp if k not in ok)
    hp = {k: v for k, v in hp.items() if k in ok}
    hp.update(precpt_gene_emb=None, transformer="normal")
    m = scPrint(**hp)
    sd = {k: v for k, v in ck["state_dict"].items() if k != "expr_decoder.fc.0.weight"}
    r = m.load_state_dict(sd, strict=False)
    info = dict(dropped_hparams=dropped, missing_keys=list(r.missing_keys), unexpected_keys=list(r.unexpected_keys),
                label_decoders=ck["hyper_parameters"]["label_decoders"], classes=list(m.classes))
    return m.eval().to(device), info


def ensembl_table():
    if os.path.exists(ENS_TSV):
        return pd.read_csv(ENS_TSV, sep="\t")
    import bionty as bt
    df = bt.Gene.public(organism="human").df()
    df.to_csv(ENS_TSV, sep="\t", index=False)
    return df


def map_genes(symbols, vocab):
    """Exact HGNC-symbol match to Ensembl release 110, then to the checkpoint vocabulary."""
    ens = ensembl_table()
    ens = ens[ens.ensembl_gene_id.str.startswith("ENSG")]
    vset = set(vocab)
    by_sym = ens.groupby("symbol")["ensembl_gene_id"].apply(list).to_dict()
    rows = []
    for s in symbols:
        ids = by_sym.get(s, [])
        inv = [i for i in ids if i in vset]
        if not ids:
            st, e = "no_ensembl_symbol", None
        elif not inv:
            st, e = "not_in_scprint_vocab", ids[0]
        elif len(inv) > 1:
            st, e = "ambiguous", None
        else:
            st, e = "ok", inv[0]
        rows.append((s, e, st))
    df = pd.DataFrame(rows, columns=["symbol", "ensembl", "status"])
    dup = df[df.status == "ok"].duplicated("ensembl", keep=False)       # two symbols to one Ensembl ID
    df.loc[dup[dup].index, "status"] = "duplicate_ensembl"
    return df


@torch.no_grad()
def run(model, gene_pos, x, device, bs, return_hidden=False):
    """x: cells x G (log2(1+x) is applied inside the model: normalization='log', so pass back-transformed TPM).

    Returns mean (layers, heads, G, G) of row-renormalised, special-token-stripped attention, plus the
    hidden states of the prefix tokens (for the organism-head check)."""
    nL, nH = model.nlayers, model.transformer.blocks[0].mixer.num_heads
    npre = model.cell_embs_count - 1                          # 9 prefix tokens, no depth token
    G = gene_pos.shape[1]
    A = torch.zeros(nL, nH, G, G, dtype=torch.float64, device=device)
    hid_pre = []
    special_mass = torch.zeros(nL, nH, dtype=torch.float64, device=device)
    n = x.shape[0]
    for s in range(0, n, bs):
        xb = x[s:s + bs].to(device)
        gb = gene_pos.expand(xb.shape[0], -1).to(device)
        enc = model._encoder(gb, xb, None, None, None)         # full_depth=None: no depth token
        assert enc.shape[1] == npre + G, enc.shape
        out, qkvs = model.transformer(enc, return_qkv=list(range(nL)))
        hid_pre.append(out[:, :npre].float().cpu())
        for l, qkv in enumerate(qkvs):
            q, k = qkv[:, :, 0].float(), qkv[:, :, 1].float()    # (B, S, H, D)
            sc = torch.einsum("bthd,bshd->bhts", q, k * q.shape[-1] ** -0.5)
            sc = sc[:, :, npre:, :]                                # queries: gene tokens only
            lse_all = torch.logsumexp(sc, dim=-1)                  # log normaliser over all tokens (specials + genes)
            lse_gene = torch.logsumexp(sc[..., npre:], dim=-1)     # log normaliser over gene keys only
            special_mass[l] += (1 - torch.exp(lse_gene - lse_all)).double().mean(2).sum(0)   # (H,) mean mass on special tokens
            # softmax over gene keys only == strip the special tokens and renormalise each row (exact; no 0/0 when a
            # row puts all its mass on the special tokens, which happens in the last layer)
            att = torch.softmax(sc[..., npre:], dim=-1)           # (B, H, G, G)
            A[l] += att.double().sum(0)
            del sc, att
        del qkvs, out, enc
    A /= n
    special_mass /= n
    return A.float().cpu().numpy(), torch.cat(hid_pre), special_mass.float().cpu().numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--ncells", type=int, default=64)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--bs", type=int, default=8)
    ap.add_argument("--permute", action="store_true", help="permute each gene independently across the selected cells")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
    dev = torch.device(f"cuda:{a.device}")
    torch.use_deterministic_algorithms(False)
    t0 = time.time()
    M = {"tag": a.tag, "device": torch.cuda.get_device_name(dev), "torch": torch.__version__, "ncells": a.ncells,
         "seed": a.seed, "permute": a.permute, "ckpt_sha256": sha256(CKPT)}
    assert M["ckpt_sha256"] == CKPT_SHA256, "checkpoint hash mismatch"
    model, info = load_scprint(dev)
    M["load_s"] = time.time() - t0
    M["shim"] = {k: info[k] for k in ("dropped_hparams", "missing_keys", "unexpected_keys")}
    M["nlayers"], M["nheads"], M["d_model"] = model.nlayers, model.transformer.blocks[0].mixer.num_heads, 256
    M["n_prefix_tokens"] = model.cell_embs_count - 1
    vocab = list(model.genes)
    M["vocab_size"] = len(vocab)

    # ---- data: case gene set (symbols), cells ----
    case = pd.read_csv(EXPR, sep="\t", index_col=0)
    raw = pd.read_csv(RAW, index_col=0)
    M["raw_shape_genes_cells"] = list(raw.shape)
    v = raw.values
    M["raw_max"], M["raw_frac_zero"] = float(v.max()), float((v == 0).mean())
    tpm_sum = (2.0 ** v - 1).sum(0)
    M["tpm_colsum_median_if_log2p1"] = float(np.median(tpm_sum))
    M["tpm_colsum_range_if_log2p1"] = [float(tpm_sum.min()), float(tpm_sum.max())]
    M["case_shape_genes_cells"] = list(case.shape)
    mp = map_genes(list(case.index), vocab)
    M["gene_status"] = mp.status.value_counts().to_dict()
    mp.to_csv(os.path.join(a.out, "gene_mapping.tsv"), sep="\t", index=False)
    use = mp[mp.status == "ok"].reset_index(drop=True)
    G = len(use)
    M["n_genes_used"] = G
    gene_pos = torch.tensor([[vocab.index(e) for e in use.ensembl]], dtype=torch.long)
    rng = np.random.default_rng(a.seed)
    cells = np.sort(rng.choice(case.shape[1], a.ncells, replace=False))
    M["cells"] = [str(c) for c in case.columns[cells]]
    lg = case.loc[use.symbol].iloc[:, cells].values.T.astype(np.float64)       # cells x G, log2(TPM+1)
    xs = 2.0 ** lg - 1.0                                                         # back to TPM; the model re-applies log2(1+x)
    if a.permute:
        prng = np.random.default_rng(1000 + a.seed)
        xs = np.stack([prng.permutation(xs[:, j]) for j in range(G)], axis=1)
    x = torch.tensor(xs, dtype=torch.float32)
    M["input_tpm_colsum_median"] = float(np.median(xs.sum(1)))

    # ---- inference ----
    torch.cuda.reset_peak_memory_stats(dev)
    t1 = time.time()
    A, hid, smass = run(model, gene_pos, x, dev, a.bs)
    torch.cuda.synchronize(dev)
    M["infer_s"] = time.time() - t1
    M["gpu_peak_alloc_MiB"] = torch.cuda.max_memory_allocated(dev) / 2**20
    M["gpu_peak_reserved_MiB"] = torch.cuda.max_memory_reserved(dev) / 2**20
    M["special_token_mass_layer_x_head"] = np.round(smass, 4).tolist()
    M["A_shape"] = list(A.shape)
    M["A_nan"], M["A_inf"] = int(np.isnan(A).sum()), int(np.isinf(A).sum())
    M["A_row_sum_minmax"] = [float(A.sum(-1).min()), float(A.sum(-1).max())]
    M["A_mean_selfattn_times_G"] = float(np.mean([np.diagonal(A[l, h]).mean() * G for l in range(A.shape[0]) for h in range(A.shape[1])]))  # 1 = uniform
    M["A_layer_head_mean_max_entry"] = float(A.mean((0, 1)).max())
    # per-layer-head mean row entropy relative to uniform (sink / hub check)
    ent = -(A * np.log(np.clip(A, 1e-30, None))).sum(-1).mean(-1) / np.log(G)
    M["A_rel_entropy_layer_x_head"] = np.round(ent, 3).tolist()
    # key (column) concentration of the all-layer all-head mean: top key's share of a row
    Am = A.mean((0, 1))
    M["Amean_top_key_col_mean"] = float(Am.mean(0).max() * G)
    # ---- organism-head check (functional sanity of the loaded encoder; not a truth) ----
    from scprint.model.model import scPrint  # noqa
    ld = info["label_decoders"]
    cls = info["classes"]
    heads = {}
    for i, c in enumerate(cls):
        dec = model.cls_decoders[c].to("cpu").eval()
        for off, name in ((1, "idx_1+i"), (2, "idx_2+i")):
            if off + i >= hid.shape[1]:
                continue
            with torch.no_grad():
                p = dec(hid[:, off + i]).argmax(-1).numpy()
            lab = pd.Series([ld[c].get(int(k), str(int(k))) for k in p]).value_counts().head(3).to_dict()
            heads.setdefault(c, {})[name] = lab
    M["head_predictions_top3_by_token_offset"] = heads
    np.savez_compressed(os.path.join(a.out, f"attn_{a.tag}.npz"), A=A.astype(np.float32), genes=np.array(use.symbol))
    M["total_s"] = time.time() - t0
    with open(os.path.join(a.out, f"metrics_{a.tag}.json"), "w") as f:
        json.dump(M, f, indent=1, default=str)
    print(json.dumps({k: M[k] for k in ("tag", "device", "n_genes_used", "A_shape", "A_nan", "infer_s", "gpu_peak_alloc_MiB", "total_s")}, indent=1))


if __name__ == "__main__":
    main()
