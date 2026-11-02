#!/usr/bin/env python
"""scFM comparator step 2.2: raw-count / linear-scale inputs for the 7 BEELINE datasets (14 cases).

No truth file is read. For each dataset it
  1. checks that results/beeline_baselines/<ds>_{500,1000}/expression.tsv equals the corresponding block of the BEELINE
     ExpressionData.csv (the cell set is all cells; the gene set is the case set),
  2. aligns the GEO deposit (by cell name, by gene symbol or Ensembl ID) to the BEELINE cells and genes,
  3. runs the transform check: BEELINE value y vs the GEO linear scale x. The BEELINE value is log2(1 + x * s_c) when the
     ratio (2^y - 1)/x is constant within every cell (s_c = cell scale); the script reports the fraction of nonzero
     entries within 1e-3 of the cell's median ratio, plus the rank-1 (gene + cell) residual of log ratios,
  4. writes the model inputs for each case to <work>/inputs/<case>.npz:
       genes (case symbols), cells, tpm_beeline (genes x cells, 2^y - 1 = BEELINE's back-transform, the preregistered
       primary input), geo_linear (GEO linear-scale values on the case genes where the deposit has any, else NaN rows),
       counts (GEO raw counts on the case genes where the deposit has counts, else NaN), counts_ok (row has counts),
       cells_ok (cell reproduces BEELINE's values from the GEO values up to its scale; geo_linear/counts are NaN in other columns).
     Rows (genes) and columns (cells) are used only if verified: a cell is verified when >= 99% of its nonzero entries have
     (2^y - 1)/x within 1e-3 of the cell's median ratio and the zero pattern agrees; a gene row likewise on verified cells.
Run:  python -I analysis/scfm_step2_counts.py --raw ~/scfm_models/raw --work ~/scfm_models/step2 --out results/scfm_step2
"""
import argparse, gzip, json, os, re, time
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BEE = ROOT + "/results/beeline_baselines"
EXPR = ROOT + "/data/beeline/BEELINE-data/inputs/scRNA-Seq"
DATASETS = ["hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L"]


def gtf_tables(path):
    """gene_id -> upper-case name, and transcript_id (no version) -> gene_id."""
    nm, tx = {}, {}
    rg, rn, rt = (re.compile(k + r' "([^"]+)"') for k in ("gene_id", "gene_name", "transcript_id"))
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln[0] == "#":
                continue
            p = ln.split("\t", 8)
            if p[2] == "gene":
                n = rn.search(p[8])
                nm[rg.search(p[8]).group(1)] = n.group(1) if n else ""
            elif p[2] == "transcript":
                tx[rt.search(p[8]).group(1)] = rg.search(p[8]).group(1)
    return nm, tx


def transform_check(lin, ref):
    """lin = 2^y - 1 (BEELINE), ref = GEO values (same shape, NaN allowed). Returns a dict of diagnostics."""
    ok = ~np.isnan(ref).any(1)
    lin, ref = lin[ok], ref[ok]
    both = (lin > 0) & (ref > 0)
    out = dict(n_genes=int(ok.sum()), zero_pattern_agree=float(((lin > 0) == (ref > 0)).mean()))
    with np.errstate(all="ignore"):
        r = np.where(both, lin / np.where(both, ref, 1), np.nan)
    cm = np.nanmedian(r, 0)
    close = np.abs(r / cm - 1) < 1e-3
    out["frac_entries_const_ratio_within_cell(1e-3)"] = float(np.nansum(close & both) / both.sum())
    out["cell_scale_s_c_quantiles(min,25,50,75,max)"] = [float(v) for v in np.nanpercentile(cm, [0, 25, 50, 75, 100])]
    m = both & (ref >= 1)
    L = np.where(m, np.log(np.where(m, lin, 1) / np.where(m, ref, 1)), np.nan)
    gm = np.nanmean(L, 1, keepdims=True)
    cmn = np.nanmean(L - gm, 0, keepdims=True)
    out["rank1_log_ratio_resid_sd"] = float(np.nanstd(L - gm - cmn))
    out["gene_effect_sd"] = float(np.nanstd(gm))
    out["cell_effect_sd"] = float(np.nanstd(cmn))
    out["max_abs_log2_diff_if_scaled_by_cell_median"] = float(np.nanmax(np.abs(np.log2(1 + lin) - np.log2(1 + ref * cm))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.expanduser("~/scfm_models/raw"))
    ap.add_argument("--work", default=os.path.expanduser("~/scfm_models/step2"))
    ap.add_argument("--out", default=ROOT + "/results/scfm_step2")
    a = ap.parse_args()
    os.makedirs(a.work + "/inputs", exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    nm, tx = gtf_tables(a.raw + "/ensembl110_mouse_gtf/Mus_musculus.GRCm39.110.gtf.gz")
    report, rows = {}, []
    for ds in DATASETS:
        t0 = time.time()
        b = pd.read_csv(f"{EXPR}/{ds}/ExpressionData.csv", index_col=0)
        info = dict(beeline_genes=int(b.shape[0]), beeline_cells=int(b.shape[1]))
        cases = {}
        for n in (500, 1000):
            c = pd.read_csv(f"{BEE}/{ds}_{n}/expression.tsv", sep="\t", index_col=0)
            same_cells = list(c.columns) == list(b.columns)
            bb = b.loc[c.index, c.columns]
            info[f"case{n}_equals_beeline_block"] = bool(same_cells and np.allclose(bb.values, c.values, atol=1e-9, rtol=0))
            info[f"case{n}_genes"] = int(c.shape[0])
            info[f"case{n}_zero_range_genes_dropped_vs_beeline_gene_ids"] = None
            cases[n] = c
        lin_all = 2.0 ** b.values.astype(np.float64) - 1.0
        info["beeline_max_y"] = float(b.values.max())
        info["beeline_frac_zero"] = float((b.values == 0).mean())
        info["beeline_linear_colsum_median"] = float(np.median(lin_all.sum(0)))
        geo_lin = None      # DataFrame genes(BEELINE symbol) x cells, GEO linear scale
        counts = None       # DataFrame, raw counts
        src = {}
        if ds == "hESC":
            g = pd.read_csv(f"{a.raw}/geo_GSE75748/GSE75748_sc_time_course_ec.csv.gz", index_col=0)
            g.columns = [x.replace(".", "_") for x in g.columns]
            counts = g
            geo_lin = g
            src = dict(accession="GSE75748", file="GSE75748_sc_time_course_ec.csv.gz", kind="RSEM expected counts (non-integer, 2 decimals); not UMI")
        elif ds == "hHep":
            g = pd.read_csv(f"{a.raw}/geo_GSE81252/GSE81252_data.cast.log2.lineage.csv.gz", index_col=0).drop(columns="experiment").T
            geo_lin = 2.0 ** g.astype(np.float64) - 1.0
            src = dict(accession="GSE81252", file="GSE81252_data.cast.log2.lineage.csv.gz", kind="log2 values only (no counts in the deposit); BEELINE values are these values")
            info["beeline_equals_geo_log2_exactly"] = bool(np.array_equal(g.loc[b.index.intersection(g.index), b.columns].values, b.loc[b.index.intersection(g.index)].values))
        elif ds == "mDC":
            g = pd.read_csv(f"{a.raw}/geo_GSE48968/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz", sep="\t", index_col=0)
            g.index = [i.upper() for i in g.index]
            g = g[~g.index.duplicated(keep=False)]
            geo_lin = g
            src = dict(accession="GSE48968", file="GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz", kind="RSEM TPM only (no counts in the deposit)")
        elif ds == "mESC":
            g = pd.read_csv(f"{a.raw}/geo_GSE98664_mergedGTF/GSE98664_tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course.txt.gz", sep="\t", index_col=0)
            g = g[b.columns.intersection(g.columns)]
            tg = {}
            rn, rt = re.compile(r'gene_name "([^"]+)"'), re.compile(r'transcript_id "([^"]+)"')
            with gzip.open(f"{a.raw}/geo_GSE98664_mergedGTF/GSE98664_mergedGTF.gtf.gz", "rt") as f:
                for ln in f:
                    if ln[0] != "#":
                        p = ln.split("\t", 8)
                        if len(p) > 8 and p[2] == "transcript":
                            tg[rt.search(p[8]).group(1)] = rn.search(p[8]).group(1)
            gid = pd.Series([tg.get(i) for i in g.index], index=g.index)
            info["mESC_transcripts"] = int(len(g))
            info["mESC_transcripts_mapped_to_gene_name"] = int(gid.notna().sum())
            gg = g[gid.notna()].groupby(gid[gid.notna()].values).sum()
            geo_lin = gg
            src = dict(accession="GSE98664", file="GSE98664_tpm_sailfish_mergedGTF_RamDA_mESC_differentiation_time_course.txt.gz (+ GSE98664_mergedGTF.gtf.gz for transcript->gene_name)", kind="sailfish transcript TPM (mergedGTF), summed to gene_name; no counts for the RamDA time course in the deposit (its count files cover 10pg and C1 cells only)")
        else:
            g = pd.read_csv(f"{a.raw}/geo_GSE81682/GSE81682_HTSeq_counts.txt.gz", sep="\t", index_col=0)
            g.columns = [x.replace("LT-HSC", "LT_HSC") for x in g.columns]
            g = g[[i.startswith("ENSMUSG") for i in g.index]]
            # symbol -> ID by exact (upper-case) Ensembl 110 name; unique only; ID must be in the deposit
            byname = {}
            for i in g.index:
                n = nm.get(i, "")
                if n:
                    byname.setdefault(n.upper(), []).append(i)
            gm = pd.read_csv(f"{a.out}/gene_map/{ds}_1000.tsv", sep="\t")
            step1 = dict(zip(gm.symbol, gm.ensembl.fillna("")))
            gm5 = pd.read_csv(f"{a.out}/gene_map/{ds}_500.tsv", sep="\t")
            step1.update({k: v for k, v in zip(gm5.symbol, gm5.ensembl.fillna("")) if v and not step1.get(k)})
            ids, route = {}, {}
            for s in b.index:
                e = step1.get(s, "")
                if e and e in g.index:
                    ids[s], route[s] = e, "step1_id"
                elif len(byname.get(s.upper(), [])) == 1:
                    ids[s], route[s] = byname[s.upper()][0], "gtf_name"
            gsel = g.loc[[ids[s] for s in ids]]
            gsel.index = list(ids)
            counts = gsel
            geo_lin = gsel
            info["mHSC_gene_route"] = pd.Series(route).value_counts().to_dict()
            src = dict(accession="GSE81682", file="GSE81682_HTSeq_counts.txt.gz", kind="HTSeq integer counts, Ensembl IDs")
        info["source"] = src
        cells_ok = [c for c in b.columns if c in geo_lin.columns]
        info["cells_matched"] = len(cells_ok)
        genes_ok = [s for s in b.index if s in geo_lin.index]
        info["genes_matched_all_beeline_genes"] = len(genes_ok)
        for attempt in (0, 1):
            ref = geo_lin.loc[genes_ok, cells_ok].values.astype(np.float64)
            lin = (2.0 ** b.loc[genes_ok, cells_ok].values.astype(np.float64)) - 1.0
            info["transform_check"] = transform_check(lin, ref) if len(genes_ok) and len(cells_ok) else None
            # per-gene verification: the gene row of the deposit must reproduce the BEELINE row up to the cell scale
            with np.errstate(all="ignore"):
                rr = np.where((lin > 0) & (ref > 0), lin / np.where(ref > 0, ref, 1), np.nan)
            cmed = np.nanmedian(rr, 0)
            cons = np.abs(rr / cmed - 1) < 1e-3
            # pass 1: gene rows that mostly agree (>= 50% of nonzero entries consistent, >= 80% zero-pattern agreement); cells are judged on these
            n1 = ((lin > 0) & (ref > 0)).sum(1)
            g1 = (np.where(n1 > 0, cons.sum(1) / np.maximum(n1, 1), 0) >= 0.5) & (((lin > 0) == (ref > 0)).mean(1) >= 0.8)
            ncell = (~np.isnan(rr[g1])).sum(0)
            cfr = np.where(ncell > 0, cons[g1].sum(0) / np.maximum(ncell, 1), 1.0)
            cell_ok = (cfr >= 0.99) & (((lin[g1] > 0) == (ref[g1] > 0)).mean(0) >= 0.99)
            info["cells_verified_reproduce_beeline"] = int(cell_ok.sum())
            info["cells_failing_verification"] = [str(c) for c in np.array(cells_ok)[~cell_ok]]
            # gene rows are verified on the verified cells only
            nz = ((lin > 0) & (ref > 0))[:, cell_ok].sum(1)
            fr = np.where(nz > 0, (cons & (lin > 0) & (ref > 0))[:, cell_ok].sum(1) / np.maximum(nz, 1), 1.0)
            za = ((lin > 0) == (ref > 0))[:, cell_ok].mean(1)
            gok = (fr >= 0.99) & (za >= 0.99)
            verified = set(np.array(genes_ok)[gok])
            bad_cells = set(np.array(cells_ok)[~cell_ok])
            info["genes_verified_row_matches_beeline"] = len(verified)
            info["genes_matched_but_failing_verification"] = len(genes_ok) - len(verified)
            info["failing_examples"] = [str(x) for x in np.array(genes_ok)[~gok][:10]]
            if attempt == 0 and (~cell_ok).any() and len(set(geo_lin.columns) - set(cells_ok)) > 0:
                # BEELINE cell labels that do not reproduce from the same-named GEO column: look for the GEO column whose
                # zero pattern equals the BEELINE cell's (>= 99.9% agreement, unique, not itself a verified BEELINE cell)
                info["cells_failing_name_match"] = [str(c) for c in np.array(cells_ok)[~cell_ok]]
                gv = geo_lin.loc[genes_ok]
                used = set(np.array(cells_ok)[cell_ok])
                cmap = {}
                for j in np.where(~cell_ok)[0]:
                    zp = ((gv.values > 0) == (lin[:, j] > 0)[:, None]).mean(0)
                    top = np.argsort(-zp)[:2]
                    if zp[top[0]] >= 0.999 and (len(top) < 2 or zp[top[1]] < 0.999) and gv.columns[top[0]] not in used:
                        cmap[cells_ok[j]] = gv.columns[top[0]]
                info["cells_rematched_by_content"] = {k: v for k, v in cmap.items()}
                if not cmap:
                    break
                newcols = {c: cmap.get(c, c) for c in cells_ok}
                geo_lin = geo_lin[[newcols[c] for c in cells_ok]].copy()
                geo_lin.columns = cells_ok
                if counts is not None:
                    counts = counts[[newcols[c] for c in cells_ok]].copy()
                    counts.columns = cells_ok
                continue
            break
        if counts is not None:
            info["counts_integer_valued"] = bool(np.all(counts.loc[genes_ok, cells_ok].values == np.round(counts.loc[genes_ok, cells_ok].values)))
            info["counts_colsum_median_matched_genes"] = float(np.median(counts.loc[genes_ok, cells_ok].values.sum(0)))
        # write case inputs
        for n, c in cases.items():
            gs = list(c.index)
            tpm_b = (2.0 ** c.values.astype(np.float64) - 1.0).astype(np.float32)
            have = [s in verified for s in gs]
            gl = np.full(c.shape, np.nan, np.float32)
            if any(have) and len(cells_ok) == c.shape[1]:
                gl[np.array(have)] = geo_lin.loc[[s for s, h in zip(gs, have) if h], list(c.columns)].values.astype(np.float32)
            cn = np.full(c.shape, np.nan, np.float32)
            hc = np.zeros(len(gs), bool)
            if counts is not None and len(cells_ok) == c.shape[1]:
                hc = np.array([s in counts.index and s in verified for s in gs])
                if hc.any():
                    cn[hc] = counts.loc[[s for s, h in zip(gs, hc) if h], list(c.columns)].values.astype(np.float32)
            cok = np.array([x not in bad_cells for x in c.columns])
            gl[:, ~cok] = np.nan
            cn[:, ~cok] = np.nan
            np.savez_compressed(f"{a.work}/inputs/{ds}_{n}.npz", genes=np.array(gs), cells=np.array(c.columns), tpm_beeline=tpm_b,
                                geo_linear=gl, counts=cn, counts_ok=hc, cells_ok=cok)
            rows.append(dict(case=f"{ds}_{n}", genes=len(gs), cells=c.shape[1], cells_matched_to_geo=len(cells_ok), cells_verified=int(cok.sum()),
                             genes_in_geo_linear_verified=int(np.sum(have)), genes_with_counts=int(hc.sum()), counts_available=bool(counts is not None),
                             geo_accession=src["accession"], geo_kind=src["kind"]))
        report[ds] = info
        print(ds, "done", round(time.time() - t0, 1), "s", flush=True)
    json.dump(report, open(a.out + "/counts_alignment.json", "w"), indent=1, default=str)
    pd.DataFrame(rows).to_csv(a.out + "/counts_alignment.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 60):
        print(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
