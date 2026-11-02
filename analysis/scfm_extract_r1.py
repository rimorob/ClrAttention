#!/usr/bin/env python
"""Step 3.5 / 3.2(d): R1 = the scPRINT package's own gene-network readout (scprint.tasks.GNInfer, scprint 2.2.0). NO TRUTH FILE IS READ.

Run in ~/venvs/scfm_alt inside the separate lamin instance (the job sets LAMIN_SETTINGS_DIR=~/scfm_models/lamin_r1/home and cds to
~/scfm_models/lamin_r1/inst; see analysis/scfm_r1_setup.py):
  CUDA_VISIBLE_DEVICES=1 python -I analysis/scfm_extract_r1.py --cases hESC_500 --arm A --out ~/scfm_models/step3/r1_A --device 1

Settings (E3, logged in <case>.r1.json as `gninfer_kwargs`): how="given", genes = covered case genes (Ensembl ids),
num_genes = max_len = len(genes), max_cells = 0 (all cells), drop_unexpressed = False, layer = None (all layers), head_agg = "mean",
preprocess = "softmax", filtration = "none", apc = False, symmetrize = False, known_grn = None, precision = "32", dtype = float32,
forward_mode = "none", doplot = False, comp_attn default; one constant cell-type label so one network is produced per case.
Unlisted GNInfer arguments: batch_size (8, as the R0 runs), num_workers 0. The adata holds every vocabulary gene of the organism
(zeros except the case genes), as scdataloader's Collator requires. The network is the package's output: softmax(Q-bar K-bar^T / sqrt(d))
per layer and head from the cell-mean Q and K, special tokens' rows and columns dropped WITHOUT renormalising (the package's rule), mean over
heads and layers; rows are queries (same convention as A). Stored: R1 (G x G float32) in the case-gene order of the gene map.
"""
import argparse, hashlib, json, os, sys, time, traceback
import numpy as np, pandas as pd, torch, scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scfm_extract_scprint as ex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="all", help="'all', 'counts', or comma list")
    ap.add_argument("--arm", default="A", choices=["A", "B"])
    ap.add_argument("--device", type=int, default=1, help="ignored: see CUDA_VISIBLE_DEVICES")
    ap.add_argument("--bs", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ncells-limit", type=int, default=0, help="TEST ONLY")
    a = ap.parse_args()
    # GNInfer moves its batches to the bare device string "cuda" (= the current device), so this script must run with
    # CUDA_VISIBLE_DEVICES=1: physical device 1 (GTX 1060) is then the only visible card, logical cuda:0. Physical device 0 is never visible.
    assert os.environ.get("CUDA_VISIBLE_DEVICES") == "1", "set CUDA_VISIBLE_DEVICES=1"
    assert torch.cuda.device_count() == 1 and "1060" in torch.cuda.get_device_name(0), torch.cuda.get_device_name(0)
    a.device = 0
    assert "lamin_r1" in os.environ.get("LAMIN_SETTINGS_DIR", ""), "run in the separate lamin instance"
    os.makedirs(a.out, exist_ok=True)
    dev = torch.device(f"cuda:{a.device}")
    import anndata, scprint
    from scprint.tasks import GNInfer
    ckpt = os.path.expanduser(ex.CKPTS["medium-v1.5"][0])
    assert ex.sha256(ckpt) == ex.CKPTS["medium-v1.5"][1]
    model, info = ex.load_model(ckpt, "none", 0, dev)
    genes_all = list(model.genes)
    gpos = {g: i for i, g in enumerate(genes_all)}
    cases = ex.CASES if a.cases == "all" else ex.COUNT_CASES if a.cases == "counts" else a.cases.split(",")
    manifest = []
    for case in cases:
        t0 = time.time()
        use, x, cells, note = ex.case_inputs(case, a.arm)
        if a.ncells_limit:
            x, cells = x[:a.ncells_limit], cells[:a.ncells_limit]
        n, G = x.shape
        org = "NCBITaxon:9606" if case.startswith("h") else "NCBITaxon:10090"
        pref = "ENSG" if org == "NCBITaxon:9606" else "ENSMUSG"
        orggenes = [g for g in genes_all if g.startswith(pref)]
        col = {g: i for i, g in enumerate(orggenes)}
        idx = np.array([col[e] for e in use.ensembl])
        Xd = sp.csr_matrix((x.astype(np.float32)))
        X = sp.csr_matrix((Xd.data, idx[Xd.indices], Xd.indptr), shape=(n, len(orggenes)))
        X.sort_indices()
        sym = pd.Series(orggenes, index=orggenes)
        sym.loc[list(use.ensembl)] = list(use.symbol)
        ad = anndata.AnnData(X=X, obs=pd.DataFrame({"organism_ontology_term_id": org, "cell_type": "all"}, index=[str(c) for c in cells]),
                             var=pd.DataFrame({"symbol": sym.values}, index=orggenes))
        kw = dict(how="given", genes=list(use.ensembl), num_genes=G, max_len=G, max_cells=0, drop_unexpressed=False, layer=None,
                  head_agg="mean", preprocess="softmax", filtration="none", apc=False, symmetrize=False, known_grn=None,
                  precision="32", dtype=torch.float32, forward_mode="none", doplot=False, batch_size=a.bs, num_workers=0,
                  cell_type_col="cell_type")
        gn = GNInfer(**kw)
        torch.cuda.reset_peak_memory_stats(dev)
        M = dict(case=case, arm=a.arm, n_cells=int(n), n_genes=int(G), scprint=scprint.__version__, gninfer_kwargs={k: str(v) for k, v in kw.items()}, **note)
        try:
            grn = gn(model, ad, cell_type=None)
            torch.cuda.synchronize(dev)
            M["gninfer_s"] = time.time() - t0
            M["gpu_peak_alloc_MiB"] = torch.cuda.max_memory_allocated(dev) / 2**20
            W = grn.varp["GRN"] if "GRN" in grn.varp else grn.grn
            W = W.toarray() if sp.issparse(W) else np.asarray(W)
            out_ids = list(grn.var.index)
            M["output_shape"] = list(W.shape)
            M["n_out_genes"] = len(out_ids)
            # reorder to the case-gene (gene-map) order; R1 must cover exactly the case genes
            pos = {g: i for i, g in enumerate(out_ids)}
            miss = [e for e in use.ensembl if e not in pos]
            M["case_genes_missing_from_output"] = len(miss)
            assert not miss and len(out_ids) == G, (len(miss), len(out_ids), G)
            o = np.array([pos[e] for e in use.ensembl])
            R1 = W[np.ix_(o, o)].astype(np.float32)
            M.update(nan=int(np.isnan(R1).sum()), inf=int(np.isinf(R1).sum()), row_sum_min=float(R1.sum(1).min()), row_sum_max=float(R1.sum(1).max()),
                     min=float(R1.min()), max=float(R1.max()), frac_zero=float((R1 == 0).mean()),
                     output_gene_order_equals_vocab_order=bool([e for e in genes_all if e in set(use.ensembl)] == out_ids))
            np.savez(f"{a.out}/{case}.r1.npz", R1=R1, genes=np.array(use.symbol), ensembl=np.array(use.ensembl))
            M["file_sha256"] = ex.sha256(f"{a.out}/{case}.r1.npz")
            # cross-check with R0: rank correlation between R1 and the all-head mean of the same case, when available (QC only)
        except Exception:
            M["error"] = traceback.format_exc()[-3000:]
        M["total_s"] = time.time() - t0
        json.dump(M, open(f"{a.out}/{case}.r1.json", "w"), indent=1, default=str)
        manifest.append({k: M.get(k) for k in ("case", "n_cells", "n_genes", "gninfer_s", "nan", "row_sum_min", "row_sum_max", "error")})
        print(json.dumps(manifest[-1], default=str)[:600], flush=True)
        torch.cuda.empty_cache()
    json.dump(dict(manifest=manifest), open(f"{a.out}/run_r1_{a.arm}.json", "w"), indent=1, default=str)
    if any(m.get("error") for m in manifest):
        sys.exit(1)


if __name__ == "__main__":
    main()
