#!/usr/bin/env python
"""K562 replication stratum, step 3.8 (2026-11-02): the preregistered gene set and the five cell subsamples, and the extractor inputs.

Preregistration (Data): "K562 replication stratum (Replogle et al. 2022, data/perturbseq/K562_essential_raw_singlecell_01.h5ad). Cells: non-targeting
control cells only, 5 seeded subsamples of 1,024 cells (seeds 1-5). Genes: CollecTRI TFs with >= 10 targets that are detected in >= 5% of control
cells, plus the 1,000 most variable detected genes."  PI decision 2026-11-02: CollecTRI may be read ONLY to count targets per TF.

What this script reads (and nothing else):
  * data/perturbseq/K562_essential_raw_singlecell_01.h5ad: obs `gene` (to find the non-targeting cells) and `cell_barcode`, var `gene_name`,
    `gene_id`, and the X rows of the NON-TARGETING CONTROL CELLS only (no perturbed cell's expression is read);
  * data/human_truth/collectri.tsv: ONLY the two columns source_genesymbol and target_genesymbol, aggregated to the number of distinct targets
    per source symbol (a regulon size). No target gene name is inspected, printed or stored. No scoring.
  * no DoRothEA, no knockdown-response (bulk) file, no other truth.

Rules fixed here (the preregistration leaves them open; they use no truth other than the CollecTRI regulon sizes):
  detected   : nonzero UMI count in >= 5% of ALL non-targeting control cells (not per subsample)
  TF part    : source_genesymbol with >= 10 distinct CollecTRI targets (counted over the whole CollecTRI table, not restricted to the gene set),
               whose symbol equals a deposit gene_name exactly (case-sensitive) and that is detected
  variable   : the 1,000 detected genes with the largest variance of log1p(count / cell total * median cell total) over all control cells
               (cell total = sum over the 8,563 matrix genes, as analysis/perturbseq.R); ties broken by gene order
  gene set   : union of the TF part and the variable part (as BEELINE's "TFs + N"), deposit order
  subsamples : k = 1..5: numpy.random.default_rng(k).choice(control cell indices in file order, 1024, replace=False), sorted
  gene ids   : mapped to scPRINT/scGPT with the E5 rule (analysis/scfm_step2_genemap.py, human)
  model input: raw UMI counts of the subsample (each model applies its own normalization); extractor case names hK562_s1..5 (the "h" prefix
               selects the human organism in the extractors); `tpm_beeline` = counts, `counts_ok` all True.
Outputs: results/scfm_k562/{geneset.tsv,subsamples.tsv,geneset_meta.json}, ~/scfm_models/step2/inputs/hK562_s{k}.npz,
results/scfm_step2/gene_map/hK562_s{k}.tsv (same columns as the BEELINE gene maps).

  ~/venvs/scfm/bin/python -I analysis/scfm_k562_prep.py
"""
import hashlib, json, os, sys
import numpy as np, pandas as pd, h5py

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H5 = f"{ROOT}/data/perturbseq/K562_essential_raw_singlecell_01.h5ad"
COLLECTRI = f"{ROOT}/data/human_truth/collectri.tsv"
OUT = f"{ROOT}/results/scfm_k562"
RAW = os.path.expanduser("~/scfm_models/raw")
INP = os.path.expanduser("~/scfm_models/step2/inputs")
MIN_TARGETS, MIN_DETECT, N_VAR, N_CELLS = 10, 0.05, 1000, 1024


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def cat(f, grp, field):
    codes = f[f"{grp}/{field}"][:]
    cats = f[f"{grp}/__categories/{field}"][:]
    return np.array([c.decode() if isinstance(c, bytes) else c for c in cats])[codes]


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(f"{ROOT}/results/scfm_step2/gene_map", exist_ok=True)
    f = h5py.File(H5, "r")
    pert = cat(f, "obs", "gene")
    barcode = np.array([b.decode() if isinstance(b, bytes) else b for b in f["obs/cell_barcode"][:]])
    sym = cat(f, "var", "gene_name")
    gid = np.array([b.decode() if isinstance(b, bytes) else b for b in f["var/gene_id"][:]])
    ctrl = np.where(pert == "non-targeting")[0]
    print(f"{len(pert)} cells, {len(sym)} genes, {len(ctrl)} non-targeting control cells", flush=True)
    C = f["X"][ctrl.tolist(), :]                 # control cells x genes, raw UMI (file order); perturbed cells are not read
    f.close()
    assert C.shape == (len(ctrl), len(sym)) and np.all(C >= 0) and np.all(C == np.round(C))
    det = (C > 0).mean(0)
    tot = C.sum(1)
    L = np.log1p(C / tot[:, None] * np.median(tot))
    var = L.var(0, ddof=1)
    detected = det >= MIN_DETECT
    # CollecTRI: regulon sizes only (distinct targets per source symbol); nothing else is kept
    col = pd.read_csv(COLLECTRI, sep="\t", usecols=["source_genesymbol", "target_genesymbol"], dtype=str)
    n_rows = len(col)
    size = col.dropna().drop_duplicates().groupby("source_genesymbol").size()
    meta = dict(collectri_sha256=sha256(COLLECTRI), collectri_rows=int(n_rows), collectri_sources=int(len(size)),
                collectri_sources_ge10=int((size >= MIN_TARGETS).sum()))
    del col
    dup = pd.Series(sym).duplicated(keep=False).values
    tf_ok = np.array([(s in size.index) and size[s] >= MIN_TARGETS for s in sym]) & detected & ~dup
    order = np.argsort(-np.where(detected & ~dup, var, -np.inf), kind="stable")
    hv = np.zeros(len(sym), bool)
    hv[order[:N_VAR]] = True
    assert hv.sum() == N_VAR and detected[hv].all()
    use = tf_ok | hv
    gs = pd.DataFrame(dict(symbol=sym, deposit_gene_id=gid, detect_frac=det, lognorm_var=var,
                           n_collectri_targets=[int(size.get(s, 0)) for s in sym], is_tf_part=tf_ok, in_variable_part=hv))[use].reset_index(drop=True)
    gs.to_csv(f"{OUT}/geneset.tsv", sep="\t", index=False)
    meta.update(n_cells_total=int(len(pert)), n_genes_matrix=int(len(sym)), n_control_cells=int(len(ctrl)), n_detected_ge5pct=int(detected.sum()),
                n_duplicate_symbols_in_matrix=int(dup.sum()), tf_part=int(tf_ok.sum()), variable_part=int(hv.sum()), overlap=int((tf_ok & hv).sum()),
                gene_set_size=int(use.sum()), median_control_umi=float(np.median(tot)), h5ad_bytes=os.path.getsize(H5))
    # subsamples
    sub = {}
    rows = []
    for k in range(1, 6):
        idx = np.sort(np.random.default_rng(k).choice(len(ctrl), N_CELLS, replace=False))
        sub[k] = idx
        for j in idx:
            rows.append((k, barcode[ctrl[j]], int(ctrl[j])))
    pd.DataFrame(rows, columns=["subsample", "cell_barcode", "file_row"]).to_csv(f"{OUT}/subsamples.tsv", sep="\t", index=False)
    # gene mapping (E5 rule, human)
    import scfm_step2_genemap as gm
    work = os.path.expanduser("~/scfm_models/step2")
    vocab = gm.scprint_vocab(work)
    vh = {x for x in vocab if x.startswith("ENSG")}
    sgv = set(json.load(open(gm.SCGPT_VOCAB)))
    HM, hg, hgnc = gm.human_mapper(RAW, vh)
    tfs = set(pd.read_csv(f"{ROOT}/data/beeline/human-tfs.csv").TF.str.upper())
    rows = []
    for s in gs.symbol:
        i, route, cur = HM.resolve(s)
        st = route if i is None else ("ok" if i in vh else "not_in_scprint_vocab")
        tok, tr = "", ""
        if s in sgv:
            tok, tr = s, "case_symbol"
        elif cur and cur in sgv:
            tok, tr = cur, "current_symbol"
        rows.append((s, i or "", route if i else "", cur or "", st, tok, tr))
    df = pd.DataFrame(rows, columns=["symbol", "ensembl", "route", "current_symbol", "status", "scgpt_token", "scgpt_route"])
    ok = df.status == "ok"
    exact = df.route.str.contains("exact_primary|tiebreak_|curated_id")
    dupid = df[ok].duplicated("ensembl", keep=False)
    for idx in df[ok][dupid].index:
        grp = df[ok & (df.ensembl == df.ensembl[idx])]
        ex = grp[exact[grp.index]]
        if len(ex) == 1 and idx != ex.index[0]:
            df.loc[idx, "status"] = "duplicate_id"
        elif len(ex) != 1:
            df.loc[idx, "status"] = "duplicate_id"
    df["in_scprint"] = df.status == "ok"
    dd = df[df.scgpt_token != ""].duplicated("scgpt_token", keep="first")
    df.loc[dd[dd].index, "scgpt_token"] = ""
    df["in_scgpt"] = df.scgpt_token != ""
    df["is_tf"] = df.symbol.str.upper().isin(tfs)
    # agreement with the deposit's own Ensembl ids (information only)
    agree = (df.ensembl.values == gs.deposit_gene_id.values) & df.in_scprint.values
    meta.update(scprint_covered=int(df.in_scprint.sum()), scgpt_covered=int(df.in_scgpt.sum()),
                intersection=int((df.in_scprint & df.in_scgpt).sum()), tf_part_covered=int((df.in_scprint & gs.is_tf_part).sum()),
                mapped_equal_to_deposit_gene_id=int(agree.sum()),
                drop_reasons=df[~df.in_scprint].status.value_counts().to_dict())
    os.makedirs(INP, exist_ok=True)
    cols = np.where(use)[0]
    for k in range(1, 6):
        case = f"hK562_s{k}"
        df.to_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t", index=False)
        X = C[sub[k]][:, cols].T.astype(np.float64)      # genes x cells, raw UMI counts
        np.savez(f"{INP}/{case}.npz", genes=np.array(gs.symbol, dtype=object), cells=np.array(barcode[ctrl[sub[k]]], dtype=object),
                 tpm_beeline=X, counts=X, counts_ok=np.ones(len(gs), bool), cells_ok=np.ones(N_CELLS, bool))
    meta["geneset_sha256"] = sha256(f"{OUT}/geneset.tsv")
    meta["subsamples_sha256"] = sha256(f"{OUT}/subsamples.tsv")
    meta["gene_map_sha256"] = sha256(f"{ROOT}/results/scfm_step2/gene_map/hK562_s1.tsv")
    json.dump(meta, open(f"{OUT}/geneset_meta.json", "w"), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
