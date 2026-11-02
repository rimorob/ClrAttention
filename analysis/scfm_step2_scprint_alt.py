#!/usr/bin/env python
"""scFM comparator step 2.4: does a released scprint version match the medium-v1.5 checkpoint?

Runs in the SEPARATE venv ~/venvs/scfm_alt (scprint 2.2.0 --no-deps + simpler-flash 1.0.7; torch/lightning/etc. are borrowed
read-only from ~/venvs/scfm through a .pth file; ~/venvs/scfm is not modified).
  python -I analysis/scfm_step2_scprint_alt.py --device 1 --smoke ~/scfm_models/smoke_out --out results/scfm_step2

1. Build scprint.scPrint 2.2.0 from the checkpoint's own hyper_parameters (only precpt_gene_emb -> None, the path is on the authors'
   cluster, and transformer -> "normal" = no flash-attn) and load the state dict with strict=True.
2. Re-run the 64 hESC cells of the step-1 smoke test (same cells, same genes, same input scale) and compare the layer x head
   attention (special tokens stripped, rows renormalised) with the step-1 shim result attn_run1.npz.
No truth file is read.
"""
import argparse, inspect, json, os
import numpy as np, pandas as pd, torch

CKPT = os.path.expanduser("~/scfm_models/scprint/medium-v1.5.ckpt")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--smoke", default=os.path.expanduser("~/scfm_models/smoke_out"))
    ap.add_argument("--out", default=ROOT + "/results/scfm_step2")
    a = ap.parse_args()
    import scprint
    from scprint import scPrint
    R = {"scprint_version": scprint.__version__, "torch": torch.__version__}
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    hp = {k: v for k, v in ck["hyper_parameters"].items() if k != "_instantiator"}
    sig = set(inspect.signature(scPrint.__init__).parameters)
    R["hparams_not_in_init_signature (go to **flash_attention_kwargs)"] = sorted(k for k in hp if k not in sig)
    R["depth_atinput_in_ckpt"] = hp.get("depth_atinput")
    hp.update(precpt_gene_emb=None, transformer="normal")
    # the package's own route: Lightning load_from_checkpoint (runs scPrint.on_load_checkpoint, which renames
    # gene_encoder.embedding.weight -> embeddings.weight and drops the training-only grad_reverse_discriminator keys), strict
    try:
        m = scPrint.load_from_checkpoint(CKPT, map_location="cpu", precpt_gene_emb=None, transformer="normal", strict=True)
        R["strict_load"] = "ok via scPrint.load_from_checkpoint(strict=True) with only precpt_gene_emb=None, transformer='normal' overridden"
    except Exception as e:
        R["load_from_checkpoint_error"] = (type(e).__name__ + ": " + str(e))[:1500]
        m = scPrint(**hp)
        try:
            m.on_load_checkpoint(ck)
        except Exception as e2:
            R["on_load_checkpoint_error"] = str(e2)[:300]
        r = m.load_state_dict(ck["state_dict"], strict=False)
        R["strict_load"] = "load_from_checkpoint failed; manual hook + strict=False"
        R["missing_keys"], R["unexpected_keys"] = list(r.missing_keys), list(r.unexpected_keys)
    sd = ck["state_dict"]
    R["cell_embs_count"] = int(m.cell_embs_count)
    R["has_depth_encoder"] = hasattr(m, "depth_encoder")
    R["expr_decoder_fc0_weight_shape_ckpt"] = list(sd["expr_decoder.fc.0.weight"].shape)
    R["expr_decoder_fc0_weight_shape_model"] = list(m.expr_decoder.fc[0].weight.shape) if hasattr(m.expr_decoder, "fc") else None
    dev = torch.device(f"cuda:{a.device}")
    m = m.eval().to(dev)

    # ---- parity with the step-1 shim on the same 64 cells ----
    z = np.load(os.path.join(a.smoke, "attn_run1.npz"), allow_pickle=True)
    A1, genes = z["A"], list(z["genes"])
    M1 = json.load(open(os.path.join(a.smoke, "metrics_run1.json")))
    gm = pd.read_csv(os.path.join(a.smoke, "gene_mapping.tsv"), sep="\t").set_index("symbol")
    vocab = list(m.genes)
    gene_pos = torch.tensor([[vocab.index(gm.loc[g, "ensembl"]) for g in genes]], dtype=torch.long)
    case = pd.read_csv(ROOT + "/results/beeline_baselines/hESC_1000/expression.tsv", sep="\t", index_col=0)
    cells = M1["cells"]
    x = torch.tensor((2.0 ** case.loc[genes, cells].values.T.astype(np.float64) - 1.0), dtype=torch.float32)
    nL, nH = m.nlayers, m.transformer.blocks[0].mixer.num_heads
    npre = m.cell_embs_count - (1 if m.depth_atinput else 0)           # tokens before the genes
    G = len(genes)
    A = torch.zeros(nL, nH, G, G, dtype=torch.float64, device=dev)
    with torch.no_grad():
        for s in range(0, x.shape[0], 8):
            xb = x[s:s + 8].to(dev)
            gb = gene_pos.expand(xb.shape[0], -1).to(dev)
            enc = m._encoder(gb, xb, None, None, None)                  # no depth token fed: req_depth=None
            R["encoder_input_tokens"] = int(enc.shape[1]); R["expected_prefix_tokens"] = int(npre)
            assert enc.shape[1] == npre + G, (enc.shape, npre, G)
            out, qkvs = m.transformer(enc, return_qkv=list(range(nL)))
            for l, qkv in enumerate(qkvs):
                q, k = qkv[:, :, 0].float(), qkv[:, :, 1].float()
                sc = torch.einsum("bthd,bshd->bhts", q, k * q.shape[-1] ** -0.5)[:, :, npre:, npre:]
                A[l] += torch.softmax(sc, -1).double().sum(0)
    A = (A / x.shape[0]).float().cpu().numpy()
    R["A_shape"] = list(A.shape)
    d = np.abs(A - A1)
    R["parity_vs_step1_shim"] = dict(max_abs_diff=float(d.max()), rel_L1=float(d.sum() / np.abs(A1).sum()),
                                     mean_abs_diff_over_mean_entry=float(d.mean() / A1.mean()))
    R["parity_per_layer_relL1"] = [float(d[l].sum() / np.abs(A1[l]).sum()) for l in range(nL)]
    # sensitivity check: how large is the difference if a depth token were (wrongly) prepended? not run; the checkpoint has no depth_encoder weight.
    os.makedirs(a.out, exist_ok=True)
    json.dump(R, open(os.path.join(a.out, "scprint_alt_check.json"), "w"), indent=1)
    print(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
