#!/usr/bin/env python
"""Step 3 production extraction of scPRINT attention (R0, R0', R2 and diagnostics). NO TRUTH FILE IS READ.

Preregistration: docs/preregistration_scfm.md (change log 2026-11-01, entries E1-E12).
Route (E3): scprint 2.2.0 `scPrint.load_from_checkpoint(ckpt, precpt_gene_emb=None, transformer="normal", strict=True)`,
run in ~/venvs/scfm_alt (python -I). Device: one CUDA device, never device 0 (enforced).

  python -I analysis/scfm_extract_scprint.py --cases all --arm A --out ~/scfm_models/step3/scprint_medium_A --device 1
  python -I analysis/scfm_extract_scprint.py --cases hESC_1000 --reproduce-step1 --ncells-check 64 ...   (unit check)

Inputs: ~/scfm_models/step2/inputs/<case>.npz (arm A: `tpm_beeline` = 2^y - 1 of the BEELINE matrix; arm B: GEO `counts`,
used only on genes with `counts_ok`) and results/scfm_step2/gene_map/<case>.tsv (rows with status == ok). All cells, case-gene
order = the gene-map order. The encoder applies log2(1 + x) itself (normalization="log"); no depth token (depth_atinput False).

Per case, writes (OUT = --out):
  <case>.attn.npy     (layers, heads, G, G) float32 cell-mean of the row-renormalised, special-token-stripped attention
                      (query x key; softmax over gene keys only = strip the special tokens and renormalise, exact)
                      [only with --store-attn, default for zero-shot runs]
  <case>.nets.npz     R0 (mean over all layers and heads), R0p (R0': heads with cell-mean special-token mass < 0.5, or the
                      4 lowest-mass heads if fewer than 4 pass), R2 (last layer, head mean), genes, ensembl, head set
  <case>.json         diagnostics: special-token mass and relative row entropy per layer x head, R0' head set, organism-head
                      accuracy and class-head predictions, NaN / row-sum checks, split-half stability, timings, SHA-256 of files
Controls (--control): perm = expression permuted independently per gene across cells (seed --seed);
rand_esm = random-weight model of the same architecture (seed --seed), ESM2 gene embeddings kept;
rand_emb = random weights including the gene embeddings. Controls store R0 only (--store-attn off).
"""
import argparse, hashlib, json, os, sys, time
import numpy as np, pandas as pd, torch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEP2 = os.path.expanduser("~/scfm_models/step2/inputs")
CKPTS = {"medium-v1.5": ("~/scfm_models/scprint/medium-v1.5.ckpt", "96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8"),
         "large-v1": ("~/scfm_models/scprint/large-v1.ckpt", "c0c53989dfc0121cb7502652f643b491e43538e5bd616b9334c504ac6be03ef9"),
         "small-v1": ("~/scfm_models/scprint/small-v1.ckpt", "269bc14be52fd0264631facecef400cd95bb1e03f36c23ff400ba02d7c9cd090")}
CASES = [f"{d}_{n}" for d in ("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L") for n in (500, 1000)]
COUNT_CASES = [c for c in CASES if c.split("_")[0] in ("hESC", "mHSC-E", "mHSC-GM", "mHSC-L")]
SPLIT_SEED = 20261021          # E10: split-half seed
MASS_THR, MIN_HEADS = 0.5, 4   # E7: R0' rule


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def array_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def load_model(ckpt, control, seed, device):
    """Returns (model, info). Zero-shot: the package's own load_from_checkpoint, strict (E3)."""
    from scprint import scPrint
    m = scPrint.load_from_checkpoint(ckpt, map_location="cpu", precpt_gene_emb=None, transformer="normal", strict=True)
    info = {"route": "scprint.scPrint.load_from_checkpoint(strict=True)"}
    if control in ("rand_esm", "rand_emb"):
        import inspect
        ck = torch.load(ckpt, map_location="cpu", weights_only=False)
        hp = {k: v for k, v in ck["hyper_parameters"].items() if k != "_instantiator"}
        hp.update(precpt_gene_emb=None, transformer="normal")
        torch.manual_seed(seed)
        r = scPrint(**hp)                       # package-native random initialisation (utils._init_weights)
        if control == "rand_esm":
            r.gene_encoder.load_state_dict(m.gene_encoder.state_dict())   # keep the pretrained (ESM2-derived) gene embeddings
        info.update(route=f"random weights ({control}), seed {seed}, scPrint(**ckpt hparams) + torch.manual_seed",
                    gene_emb_std_pretrained=float(m.gene_encoder.embeddings.weight.std()),
                    gene_emb_std_model=float(r.gene_encoder.embeddings.weight.std()))
        del m
        m = r
        del ck
    return m.eval().to(device), info


@torch.no_grad()
def run(model, gene_pos, x, device, bs, half_mask=None, want_hidden=True):
    """x: cells x G (log2(1+x) applied inside the model). Same arithmetic as analysis/scfm_smoke_scprint.py::run, plus the
    half-sample sums. Returns Asum (float64, L,H,G,G on device, SUM over cells), A1sum (cells with half_mask True), mass sums,
    hidden prefix states."""
    nL, nH = model.nlayers, model.transformer.blocks[0].mixer.num_heads
    npre = model.cell_embs_count - (1 if model.depth_atinput else 0)   # tokens before the genes (9: base + 8 class tokens)
    G = gene_pos.shape[1]
    A = torch.zeros(nL, nH, G, G, dtype=torch.float64, device=device)
    A1 = torch.zeros_like(A) if half_mask is not None else None
    mass = torch.zeros(nL, nH, dtype=torch.float64, device=device)
    hid = []
    n = x.shape[0]
    hm = torch.tensor(half_mask) if half_mask is not None else None
    for s in range(0, n, bs):
        xb = x[s:s + bs].to(device)
        gb = gene_pos.expand(xb.shape[0], -1).to(device)
        enc = model._encoder(gb, xb, None, None, None)
        assert enc.shape[1] == npre + G, enc.shape
        out, qkvs = model.transformer(enc, return_qkv=list(range(nL)))
        if want_hidden:
            hid.append(out[:, :npre].float().cpu())
        sel = hm[s:s + bs].to(device) if hm is not None else None
        for l, qkv in enumerate(qkvs):
            q, k = qkv[:, :, 0].float(), qkv[:, :, 1].float()
            sc = torch.einsum("bthd,bshd->bhts", q, k * q.shape[-1] ** -0.5)[:, :, npre:, :]
            lse_all = torch.logsumexp(sc, dim=-1)
            lse_gene = torch.logsumexp(sc[..., npre:], dim=-1)
            mass[l] += (1 - torch.exp(lse_gene - lse_all)).double().mean(2).sum(0)
            att = torch.softmax(sc[..., npre:], dim=-1).double()
            A[l] += att.sum(0)
            if A1 is not None and bool(sel.any()):
                A1[l] += att[sel].sum(0)
            del sc, att
        del qkvs, out, enc
    return A, A1, mass, (torch.cat(hid) if want_hidden else None)


def head_set(mass):
    """E7: heads with cell-mean special-token mass < 0.5; if fewer than 4, the 4 lowest-mass heads. Returns sorted (l, h) list."""
    flat = mass.reshape(-1)
    idx = np.where(flat < MASS_THR)[0]
    rule = "mass<0.5"
    if len(idx) < MIN_HEADS:
        idx = np.argsort(flat, kind="stable")[:MIN_HEADS]
        rule = "4 lowest mass (fewer than 4 below 0.5)"
    idx = np.sort(idx)
    H = mass.shape[1]
    return [(int(i // H), int(i % H)) for i in idx], rule


def offdiag_stats(a, b):
    """relL1 and Spearman between two G x G networks on the off-diagonal (self pairs excluded, as they are scored 0)."""
    from scipy.stats import rankdata
    G = a.shape[0]
    m = ~np.eye(G, dtype=bool)
    x, y = a[m].astype(np.float64), b[m].astype(np.float64)
    rel = float(np.abs(x - y).sum() / (0.5 * (x.sum() + y.sum())))
    rx, ry = rankdata(x), rankdata(y)
    return rel, float(np.corrcoef(rx, ry)[0, 1])


def case_inputs(case, arm):
    z = np.load(f"{STEP2}/{case}.npz", allow_pickle=True)
    gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t")
    assert list(gm.symbol) == list(z["genes"])
    ok = (gm.status == "ok").values
    note = {}
    if arm == "A":
        X = z["tpm_beeline"]
    else:
        assert z["counts_ok"].any(), f"{case} has no counts (arm B)"
        okc = z["counts_ok"].astype(bool)
        note["genes_ok_without_counts"] = int((ok & ~okc).sum())
        ok = ok & okc
        X = z["counts"]
    assert bool(z["cells_ok"].all())
    use = gm[ok].reset_index(drop=True)
    x = np.asarray(X[ok], dtype=np.float64).T          # cells x G
    assert np.isfinite(x).all() and (x >= 0).all()
    return use, x, list(z["cells"]), note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="all", help="'all', 'counts' (the 8 arm-B cases), or comma list")
    ap.add_argument("--arm", default="A", choices=["A", "B"])
    ap.add_argument("--ckpt-name", default="medium-v1.5", choices=list(CKPTS))
    ap.add_argument("--control", default="none", choices=["none", "perm", "rand_esm", "rand_emb"])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", type=int, default=1)
    ap.add_argument("--bs", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--store-attn", choices=["auto", "yes", "no"], default="auto")
    ap.add_argument("--ncells-limit", type=int, default=0, help="TEST ONLY: first N cells")
    ap.add_argument("--reproduce-step1", default="", help="step-1 smoke dir: reproduce attn_run1.npz exactly (hESC_1000, 64 cells) and exit")
    a = ap.parse_args()
    assert a.device != 0, "device 0 is the desktop GPU: never used"
    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "2")))
    dev = torch.device(f"cuda:{a.device}")
    torch.use_deterministic_algorithms(False)
    os.makedirs(a.out, exist_ok=True)
    ckpt = os.path.expanduser(CKPTS[a.ckpt_name][0])
    t0 = time.time()
    import scprint
    R = {"scprint": scprint.__version__, "torch": torch.__version__, "device": torch.cuda.get_device_name(dev), "args": vars(a),
         "ckpt": ckpt, "ckpt_sha256": sha256(ckpt)}
    assert R["ckpt_sha256"] == CKPTS[a.ckpt_name][1], "checkpoint hash mismatch"
    model, info = load_model(ckpt, a.control, a.seed, dev)
    R.update(info, load_s=time.time() - t0)
    vocab = list(model.genes)
    vidx = {g: i for i, g in enumerate(vocab)}
    ck_hp = torch.load(ckpt, map_location="cpu", weights_only=False)["hyper_parameters"]
    ldec, classes = ck_hp["label_decoders"], list(model.classes)
    R["nlayers"], R["nheads"] = model.nlayers, model.transformer.blocks[0].mixer.num_heads
    R["n_prefix_tokens"] = model.cell_embs_count - (1 if model.depth_atinput else 0)
    R["depth_atinput"] = bool(model.depth_atinput)

    if a.reproduce_step1:
        sm = os.path.expanduser(a.reproduce_step1)
        z = np.load(f"{sm}/attn_run1.npz", allow_pickle=True)
        A1ref, genes = z["A"], list(z["genes"])
        M1 = json.load(open(f"{sm}/metrics_run1.json"))
        gmap = pd.read_csv(f"{sm}/gene_mapping.tsv", sep="\t").set_index("symbol")
        case = pd.read_csv(ROOT + "/results/beeline_baselines/hESC_1000/expression.tsv", sep="\t", index_col=0)
        x = torch.tensor(2.0 ** case.loc[genes, M1["cells"]].values.T.astype(np.float64) - 1.0, dtype=torch.float32)
        gp = torch.tensor([[vidx[gmap.loc[g, "ensembl"]] for g in genes]], dtype=torch.long)
        As, _, mass, _ = run(model, gp, x, dev, 8)
        A = (As / x.shape[0]).float().cpu().numpy()
        d = np.abs(A - A1ref)
        # same input from the step-2 npz (tpm_beeline), same cells and genes
        z2 = np.load(f"{STEP2}/hESC_1000.npz", allow_pickle=True)
        gi = {g: i for i, g in enumerate(z2["genes"])}
        ci = {c: i for i, c in enumerate(z2["cells"])}
        x2 = z2["tpm_beeline"][[gi[g] for g in genes]][:, [ci[c] for c in M1["cells"]]].T
        R["reproduce_step1"] = dict(shape=list(A.shape), max_abs_diff=float(d.max()), bitwise_identical=bool(np.array_equal(A, A1ref)),
                                    rel_L1=float(d.sum() / np.abs(A1ref).sum()),
                                    mass_diff_vs_step1=float(np.abs(np.round(mass.cpu().numpy() / x.shape[0], 4) - np.array(M1["special_token_mass_layer_x_head"])).max()),
                                    step2_npz_input_max_abs_diff_vs_step1_input=float(np.abs(x2.astype(np.float64) - x.numpy().astype(np.float64)).max()),
                                    step2_npz_input_identical=bool(np.array_equal(x2.astype(np.float32), x.numpy())))
        json.dump(R, open(f"{a.out}/reproduce_step1.json", "w"), indent=1, default=str)
        print(json.dumps(R["reproduce_step1"], indent=1))
        return

    cases = CASES if a.cases == "all" else COUNT_CASES if a.cases == "counts" else a.cases.split(",")
    if a.arm == "B":
        assert all(c in COUNT_CASES for c in cases), "arm B exists only for hESC and mHSC-E/GM/L"
    store = {"auto": a.control == "none", "yes": True, "no": False}[a.store_attn]
    exp_org = lambda c: "NCBITaxon:9606" if c.startswith("h") else "NCBITaxon:10090"
    manifest = []
    for case in cases:
        tc = time.time()
        use, x, cells, note = case_inputs(case, a.arm)
        if a.ncells_limit:
            x, cells = x[:a.ncells_limit], cells[:a.ncells_limit]
        n, G = x.shape
        if a.control == "perm":
            prng = np.random.default_rng(a.seed)
            x = np.stack([prng.permutation(x[:, j]) for j in range(G)], axis=1)
        gp = torch.tensor([[vidx[e] for e in use.ensembl]], dtype=torch.long)
        sprng = np.random.default_rng(SPLIT_SEED)
        half = np.zeros(n, dtype=bool)
        half[sprng.permutation(n)[: n // 2]] = True
        torch.cuda.reset_peak_memory_stats(dev)
        t1 = time.time()
        As, A1s, mass, hid = run(model, gp, torch.tensor(x, dtype=torch.float32), dev, a.bs, half_mask=half)
        torch.cuda.synchronize(dev)
        infer_s = time.time() - t1
        A = (As / n).float().cpu().numpy()                                  # (L,H,G,G)
        mass = (mass / n).float().cpu().numpy()
        A1 = (A1s / half.sum()).float().cpu().numpy()
        A2 = ((As - A1s) / (n - half.sum())).float().cpu().numpy()
        del As, A1s
        nL, nH = A.shape[:2]
        hs, rule = head_set(mass)
        R0 = A.mean((0, 1))
        R0p = np.mean([A[l, h] for l, h in hs], axis=0)
        R2 = A[-1].mean(0)
        ent = -(A * np.log(np.clip(A, 1e-30, None))).sum(-1).mean(-1) / np.log(G)
        M = dict(case=case, arm=a.arm, control=a.control, seed=a.seed, ckpt=a.ckpt_name, n_cells=int(n), n_genes=int(G), **note,
                 infer_s=infer_s, gpu_peak_alloc_MiB=torch.cuda.max_memory_allocated(dev) / 2**20,
                 attn_shape=list(A.shape), attn_nan=int(np.isnan(A).sum()), attn_inf=int(np.isinf(A).sum()),
                 row_sum_min=float(A.sum(-1).min()), row_sum_max=float(A.sum(-1).max()),
                 special_token_mass=np.round(mass, 5).tolist(), rel_entropy=np.round(ent, 4).tolist(),
                 r0p_heads=hs, r0p_rule=rule, r0p_n_heads=len(hs), r0p_mass_threshold=MASS_THR,
                 attn_array_sha256=array_sha(A))
        # split-half stability (cells split once, seed 20261021)
        sh = {}
        for name, f in (("R0", lambda T: T.mean((0, 1))), ("R0p", lambda T: np.mean([T[l, h] for l, h in hs], axis=0)), ("R2", lambda T: T[-1].mean(0))):
            r, s = offdiag_stats(f(A1), f(A2))
            sh[name] = dict(relL1=r, spearman=s)
        M["split_half"] = sh
        # class-head predictions (functional check of the input scale; token index 1+i, as in the step-1 smoke test)
        org = exp_org(case)
        preds = {}
        with torch.no_grad():
            for i, c in enumerate(classes):
                dec = model.cls_decoders[c].to("cpu").eval()
                p = dec(hid[:, 1 + i]).argmax(-1).numpy()
                lab = pd.Series([ldec[c].get(int(k), str(int(k))) for k in p]).value_counts()
                preds[c] = lab.head(3).to_dict()
                if c == "organism_ontology_term_id":
                    M["organism_expected"], M["organism_accuracy"] = org, float(lab.get(org, 0) / n)
        M["class_head_top3"] = preds
        # files
        base = f"{a.out}/{case}"
        if store:
            np.save(base + ".attn.npy", A)
            M["attn_file_sha256"] = sha256(base + ".attn.npy")
        np.savez(base + ".nets.npz", R0=R0.astype(np.float32), R0p=R0p.astype(np.float32), R2=R2.astype(np.float32),
                 genes=np.array(use.symbol), ensembl=np.array(use.ensembl), r0p_heads=np.array(hs))
        M["nets_file_sha256"] = sha256(base + ".nets.npz")
        M["total_s"] = time.time() - tc
        json.dump(M, open(base + ".json", "w"), indent=1, default=str)
        manifest.append(dict(case=case, n=n, G=G, infer_s=round(infer_s, 1), nan=M["attn_nan"], org_acc=M.get("organism_accuracy"),
                             r0p_n=len(hs), sh_R0=round(sh["R0"]["relL1"], 4)))
        print(json.dumps(manifest[-1]), flush=True)
        del A, A1, A2, hid
        torch.cuda.empty_cache()
    R["total_s"] = time.time() - t0
    R["manifest"] = manifest
    json.dump(R, open(f"{a.out}/run_{a.arm}_{a.control}{a.seed}.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
