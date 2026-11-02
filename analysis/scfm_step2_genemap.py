#!/usr/bin/env python
"""scFM comparator step 2.1: gene-symbol -> primary-assembly Ensembl ID mapping for the 14 BEELINE cases.

No truth file is read (BEELINE gene lists come from results/beeline_baselines/<case>/expression.tsv; the TF lists
data/beeline/{human,mouse}-tfs.csv define the gene sets, they are not truths).

Run (python -I; reads downloaded reference files, see docs/scfm_step2.md):
  python -I analysis/scfm_step2_genemap.py --raw ~/scfm_models/raw --out results/scfm_step2 --work ~/scfm_models/step2

MAPPING RULE (documented here and in docs/scfm_step2.md)
Human (Ensembl release 110 GRCh38 GTF, HGNC complete set downloaded 2026-11-01):
  primary = gene on chromosomes 1-22, X, Y, MT; everything else (patch, haplotype, scaffold) is "non-primary".
  1. exact_primary   : exactly one primary-chromosome gene has gene_name == symbol -> that ID.
  2. tiebreak        : several primary genes share the name -> narrow by (a) membership in the scPRINT checkpoint vocabulary,
                       (b) the HGNC/MGI-curated Ensembl ID, (c) biotype protein_coding, (d) the longest genomic span (the full locus, not a partial duplicate model); applied in this order, each only if non-empty;
                       unique -> that ID, else ambiguous.
  3. hgnc_id         : no primary gene_name match, but the symbol is an HGNC-approved symbol whose HGNC Ensembl ID is a primary gene.
  4. prev_symbol     : symbol is a unique HGNC previous symbol -> resolve the current approved symbol by 1-3.
  5. alias_symbol    : symbol is a unique HGNC alias -> resolve the current approved symbol by 1-3.
  6. excel_date      : "1-Mar"-style Excel mangling -> MARCH/SEPT/DEC symbols -> 1-5. (reported separately)
  Not mapped: spike_in (ERCC-), patch_only (a gene_name exists only on non-primary contigs and HGNC has no primary ID),
  ambiguous (several candidate IDs / several current symbols), no_match.
Mouse (Ensembl release 110 GRCm39 GTF, MGI MRK_List2 + MRK_ENSEMBL, downloaded 2026-11-01): the same ladder with
  primary = 1-19, X, Y, MT, matching case-insensitively (the mDC and mHSC BEELINE files use upper-case mouse symbols
  such as QK, PARK2, CLEC1B), MGI symbol->Ensembl (MRK_ENSEMBL) in place of HGNC IDs, MGI synonyms in place of
  prev/alias symbols.
Then: ID must be in the scPRINT checkpoint vocabulary of that organism (ENSG*/ENSMUSG*), else "not_in_scprint_vocab".
Two symbols resolving to the same ID: exact-route symbol is kept; if neither is exact both are dropped ("duplicate_id").
scGPT (human cases only): token = the case symbol if in the 60,697-symbol vocabulary, else the resolved current symbol.
"""
import argparse, gzip, os, re, json, hashlib, collections
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BEE = ROOT + "/results/beeline_baselines"
DATASETS = {"hESC": "human", "hHep": "human", "mDC": "mouse", "mESC": "mouse", "mHSC-E": "mouse", "mHSC-GM": "mouse", "mHSC-L": "mouse"}
HUMAN_PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
MOUSE_PRIMARY = {str(i) for i in range(1, 20)} | {"X", "Y", "MT"}
CKPT = os.path.expanduser("~/scfm_models/scprint/medium-v1.5.ckpt")
SCGPT_VOCAB = os.path.expanduser("~/scfm_models/scgpt_human/vocab.json")
MONTH = {"mar": "MARCH", "sep": "SEPT", "dec": "DEC"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_gtf_genes(path):
    rows = []
    rx = {k: re.compile(k + r' "([^"]*)"') for k in ("gene_id", "gene_name", "gene_biotype")}
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln[0] == "#":
                continue
            p = ln.split("\t", 8)
            if p[2] != "gene":
                continue
            a = p[8]
            g = {k: (r.search(a).group(1) if r.search(a) else "") for k, r in rx.items()}
            rows.append((g["gene_id"], g["gene_name"], g["gene_biotype"], p[0], int(p[4]) - int(p[3]) + 1))
    return pd.DataFrame(rows, columns=["gene_id", "gene_name", "biotype", "chrom", "span"])


def scprint_vocab(work):
    p = os.path.join(work, "scprint_vocab_ids.txt")
    if not os.path.exists(p):
        import torch
        ck = torch.load(CKPT, map_location="cpu", weights_only=False)
        os.makedirs(work, exist_ok=True)
        with open(p, "w") as f:
            f.write("\n".join(ck["hyper_parameters"]["genes"]) + "\n")
    ids = [x.strip() for x in open(p) if x.strip()]
    return ids


class Mapper:
    """Resolve symbols to a primary-assembly Ensembl ID. `upper` = case-insensitive matching (mouse)."""

    def __init__(self, gtf, primary, current, ens_of_current, synonyms, upper, vocab=()):
        # gtf: gene table; current: set of current official symbols; ens_of_current: official symbol -> list of IDs (HGNC/MGI)
        self.upper = upper
        self.vocab = set(vocab)
        k = (lambda s: s.upper()) if upper else (lambda s: s)
        self.k = k
        gtf = gtf.copy()
        gtf["primary"] = gtf.chrom.isin(primary)
        self.chrom = dict(zip(gtf.gene_id, gtf.chrom))
        self.biotype = dict(zip(gtf.gene_id, gtf.biotype))
        self.span = dict(zip(gtf.gene_id, gtf.span))
        self.name_of = dict(zip(gtf.gene_id, gtf.gene_name))
        self.prim_by_name = collections.defaultdict(set)
        self.any_by_name = collections.defaultdict(set)
        for gid, nm, pr in zip(gtf.gene_id, gtf.gene_name, gtf.primary):
            if nm:
                self.any_by_name[k(nm)].add(gid)
                if pr:
                    self.prim_by_name[k(nm)].add(gid)
        self.prim_ids = set(gtf.gene_id[gtf.primary])
        self.current = {k(s): s for s in current}
        self.ens_cur = {k(s): [i for i in v if i] for s, v in ens_of_current.items()}
        self.syn = collections.defaultdict(set)         # kind -> key -> set(current symbols)
        for kind, d in synonyms.items():
            m = collections.defaultdict(set)
            for cur, lst in d.items():
                for s in lst:
                    if k(s) != k(cur):
                        m[k(s)].add(cur)
            self.syn[kind] = m

    def direct(self, s):
        """Steps 1-3 for a (possibly already current) symbol. Returns (id, route, current_symbol) or (None, reason, None)."""
        ks = self.k(s)
        c = self.prim_by_name.get(ks, set())
        cur_ens = [i for i in self.ens_cur.get(ks, []) if i in self.prim_ids]
        if len(c) == 1:
            return next(iter(c)), "exact_primary", self.current.get(ks, s)
        if len(c) > 1 or len(cur_ens) > 1:
            return self.tiebreak(set(c) | set(cur_ens) if not c else set(c), set(cur_ens), ks, s)
        if len(cur_ens) == 1:
            return cur_ens[0], "curated_id", self.current.get(ks, s)
        if ks in self.any_by_name:
            return None, "patch_only", None
        return None, "no_match", None

    def tiebreak(self, cand, cur_ens, ks, s):
        """Several primary-chromosome genes: narrow by (1) the checkpoint vocabulary, (2) the curated HGNC/MGI ID,
        (3) biotype protein_coding; each narrowing is applied only if it leaves a non-empty set; unique -> that ID."""
        steps = (("vocab", lambda i: i in self.vocab), ("curated", lambda i: i in cur_ens),
                 ("protein_coding", lambda i: self.biotype.get(i) == "protein_coding"),
                 ("longest_span", lambda i: self.span.get(i) == max(self.span.get(j, 0) for j in cand)))
        used = []
        for nm, f in steps:
            if len(cand) == 1:
                break
            sub = {i for i in cand if f(i)}
            if sub:
                if sub != cand:
                    used.append(nm)
                cand = sub
        if len(cand) == 1:
            return next(iter(cand)), "tiebreak_" + "+".join(used or ["unique"]), self.current.get(ks, s)
        return None, "ambiguous_primary_ids", None

    def resolve(self, s):
        if s.startswith("ERCC-"):
            return None, "spike_in", None
        m = re.fullmatch(r"(\d+)-(Mar|Sep|Dec)", s)
        s2, pre = s, ""
        if m:
            s2, pre = MONTH[m.group(2).lower()] + m.group(1), "excel_date+"
        i, route, cur = self.direct(s2)
        if i is not None:
            return i, pre + route, cur
        if route in ("ambiguous_primary_ids",):
            return None, route, None
        for kind in ("prev", "alias"):
            cands = self.syn[kind].get(self.k(s2), set())
            if not cands:
                continue
            ok = {}
            for cs in cands:
                j, r2, c2 = self.direct(cs)
                if j is not None:
                    ok[j] = c2
            if len(ok) == 1:
                j = next(iter(ok))
                return j, pre + kind + "_symbol", ok[j]
            if len(ok) > 1 or len(cands) > 1:
                return None, "ambiguous_" + kind + "_symbol", None
        return None, route, None


def human_mapper(raw, vocab):
    gtf = read_gtf_genes(raw + "/ensembl110_human_gtf/Homo_sapiens.GRCh38.110.gtf.gz")
    h = pd.read_csv(raw + "/hgnc/hgnc_complete_set.txt", sep="\t", dtype=str).fillna("")
    app = h[h.status == "Approved"]
    cur = list(app.symbol)
    ens = {s: [e] for s, e in zip(app.symbol, app.ensembl_gene_id)}
    spl = lambda x: [t for t in x.split("|") if t]
    syn = {"prev": {s: spl(p) for s, p in zip(app.symbol, app.prev_symbol)},
           "alias": {s: spl(p) for s, p in zip(app.symbol, app.alias_symbol)}}
    return Mapper(gtf, HUMAN_PRIMARY, cur, ens, syn, upper=False, vocab=vocab), gtf, h


def mouse_mapper(raw, vocab):
    gtf = read_gtf_genes(raw + "/ensembl110_mouse_gtf/Mus_musculus.GRCm39.110.gtf.gz")
    lst = pd.read_csv(raw + "/mgi/MRK_List2.rpt", sep="\t", dtype=str).fillna("")
    lst = lst[(lst["Status"] == "O") & (lst["Marker Type"] == "Gene")]
    en = pd.read_csv(raw + "/mgi/MRK_ENSEMBL.rpt", sep="\t", header=None, dtype=str).fillna("")
    en = en.rename(columns={1: "symbol", 5: "ensembl"})
    ens_by = en.groupby("symbol")["ensembl"].apply(lambda x: sorted(set(i for i in x if i.startswith("ENSMUSG"))))
    cur = list(lst["Marker Symbol"])
    ens = {s: ens_by.get(s, []) for s in cur}
    syn = {"prev": {}, "alias": {s: [t for t in x.split("|") if t] for s, x in zip(lst["Marker Symbol"], lst["Marker Synonyms (pipe-separated)"])}}
    return Mapper(gtf, MOUSE_PRIMARY, cur, ens, syn, upper=True, vocab=vocab), gtf, lst


def case_genes(case):
    with open(f"{BEE}/{case}/expression.tsv") as f:
        f.readline()
        return [ln.split("\t", 1)[0] for ln in f]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.expanduser("~/scfm_models/raw"))
    ap.add_argument("--out", default=ROOT + "/results/scfm_step2")
    ap.add_argument("--work", default=os.path.expanduser("~/scfm_models/step2"))
    a = ap.parse_args()
    os.makedirs(a.out + "/gene_map", exist_ok=True)
    vocab = scprint_vocab(a.work)
    vh = {x for x in vocab if x.startswith("ENSG")}
    vm = {x for x in vocab if x.startswith("ENSMUSG")}
    sg = json.load(open(SCGPT_VOCAB))
    sgv = set(sg) if isinstance(sg, dict) else set(sg)
    HM, hg, hgnc = human_mapper(a.raw, vh)
    MM, mg, mgi = mouse_mapper(a.raw, vm)
    tfs = {"human": set(pd.read_csv(ROOT + "/data/beeline/human-tfs.csv").TF.str.upper()), "mouse": set(pd.read_csv(ROOT + "/data/beeline/mouse-tfs.csv").TF.str.upper())}
    summ, dropped = [], []
    for ds, org in DATASETS.items():
        M, vv = (HM, vh) if org == "human" else (MM, vm)
        for n in (500, 1000):
            case = f"{ds}_{n}"
            syms = case_genes(case)
            rows = []
            for s in syms:
                i, route, cur = M.resolve(s)
                st = route if i is None else ("ok" if i in vv else "not_in_scprint_vocab")
                tok, tr = "", ""
                if org == "human":
                    if s in sgv:
                        tok, tr = s, "case_symbol"
                    elif cur and cur in sgv:
                        tok, tr = cur, "current_symbol"
                rows.append((s, i or "", route if i else "", cur or "", st, tok, tr))
            df = pd.DataFrame(rows, columns=["symbol", "ensembl", "route", "current_symbol", "status", "scgpt_token", "scgpt_route"])
            # duplicates: two symbols -> one ID
            ok = df.status == "ok"
            exact = df.route.str.contains("exact_primary|tiebreak_|curated_id")
            dup = df[ok].duplicated("ensembl", keep=False)
            for idx in df[ok][dup].index:
                grp = df[ok & (df.ensembl == df.ensembl[idx])]
                ex = grp[exact[grp.index]]
                if len(ex) == 1 and idx != ex.index[0]:
                    df.loc[idx, "status"] = "duplicate_id"
                elif len(ex) != 1:
                    df.loc[idx, "status"] = "duplicate_id"
            df["in_scprint"] = df.status == "ok"
            if org == "human":
                dd = df[df.scgpt_token != ""].duplicated("scgpt_token", keep="first")
                df.loc[dd[dd].index, "scgpt_token"] = ""
                df["in_scgpt"] = df.scgpt_token != ""
            else:
                df["in_scgpt"] = False
            df["is_tf"] = df.symbol.str.upper().isin(tfs[org])
            df.to_csv(f"{a.out}/gene_map/{case}.tsv", sep="\t", index=False)
            n_g = len(df)
            r = dict(case=case, organism=org, genes=n_g, scprint_covered=int(df.in_scprint.sum()), tfs=int(df.is_tf.sum()),
                     tfs_scprint=int((df.is_tf & df.in_scprint).sum()))
            for stn, c in df.status.value_counts().items():
                if stn != "ok":
                    r["drop_" + stn] = int(c)
            for rt, c in df[df.in_scprint].route.str.replace(r"^excel_date\+", "", regex=True).value_counts().items():
                r["route_" + rt] = int(c)
            if org == "human":
                r.update(scgpt_covered=int(df.in_scgpt.sum()), intersection=int((df.in_scprint & df.in_scgpt).sum()),
                         tfs_scgpt=int((df.is_tf & df.in_scgpt).sum()), tfs_intersection=int((df.is_tf & df.in_scprint & df.in_scgpt).sum()))
            r["frac_scprint"] = round(r["scprint_covered"] / n_g, 4)
            summ.append(r)
            for _, x in df[~df.in_scprint].iterrows():
                dropped.append((case, x.symbol, x.status, x.is_tf, x.ensembl, x.in_scgpt))
    S = pd.DataFrame(summ).fillna(0)
    for c in S.columns:
        if S[c].dtype == float and c != "frac_scprint":
            S[c] = S[c].astype(int)
    S.to_csv(a.out + "/coverage_by_case.csv", index=False)
    pd.DataFrame(dropped, columns=["case", "symbol", "reason", "is_tf", "ensembl_candidate", "in_scgpt"]).to_csv(a.out + "/dropped_genes.csv", index=False)
    meta = dict(ensembl_human_gtf_sha256=sha256(a.raw + "/ensembl110_human_gtf/Homo_sapiens.GRCh38.110.gtf.gz"),
                ensembl_mouse_gtf_sha256=sha256(a.raw + "/ensembl110_mouse_gtf/Mus_musculus.GRCm39.110.gtf.gz"),
                hgnc_sha256=sha256(a.raw + "/hgnc/hgnc_complete_set.txt"), mgi_list2_sha256=sha256(a.raw + "/mgi/MRK_List2.rpt"),
                mgi_ensembl_sha256=sha256(a.raw + "/mgi/MRK_ENSEMBL.rpt"), scprint_vocab_human=len(vh), scprint_vocab_mouse=len(vm),
                scprint_vocab_other=len(vocab) - len(vh) - len(vm), scprint_ckpt_sha256=sha256(CKPT), scgpt_vocab_size=len(sgv),
                vocab_nonprimary_human=int(len(vh - set(hg.gene_id[hg.chrom.isin(HUMAN_PRIMARY)]))),
                vocab_nonprimary_mouse=int(len(vm - set(mg.gene_id[mg.chrom.isin(MOUSE_PRIMARY)]))))
    json.dump(meta, open(a.out + "/genemap_meta.json", "w"), indent=1)
    with pd.option_context("display.width", 250, "display.max_columns", 50):
        print(S)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
