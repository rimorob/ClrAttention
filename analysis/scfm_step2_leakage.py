#!/usr/bin/env python
"""scFM comparator step 2.3: pretraining-corpus leakage audit for the 5 BEELINE accessions + the Replogle K562 deposit.

Metadata only; no expression data and no truth file is read.
Run: python -I analysis/scfm_step2_leakage.py --raw ~/scfm_models/raw --out results/scfm_step2

Evidence sources (all public, fetched 2026-11-01 into ~/scfm_models/raw/):
  cxg_api/collections.json        CELLxGENE Discover curation API, all PUBLIC collections (current state)
  cxg_api/h5ads_<release>.xml     S3 listing of the census h5ad files of the LTS releases 2023-05-15, 2023-07-25, 2023-12-15
                                  (one <dataset_id>.h5ad per dataset in that release)
  scbasecount_meta/*.parquet      scBaseCount sample_metadata (SRX accessions), releases 2025-02-25 and 2026-01-12
  sra_runinfo/<PRJ>.csv           SRA run tables of the BioProjects behind each GEO series (SRX <-> accession)
Model corpora:
  scPRINT v1 (medium-v1.5): CELLxGENE census LTS 2023-12-15 (Kalfon et al. 2025, Data availability; the scPRINT README adds that
      "some cells, datasets were dropped due to low quality data and some were randomly removed" for validation/test).
  scGPT whole-human: "33 million normal human cells" from CELLxGENE (README); the census snapshot is NOT stated in any
      source reachable here; all three LTS releases are checked, 2023-05-15 (earliest) as the likely upper bound on age.
      Only normal (non-diseased) human cells, so mouse and K562 (cancer line) cells are out of its corpus by construction.
  scPRINT-2 (optional arm): scBaseCount (10x droplet data keyed by SRX), checked for completeness.
"""
import argparse, json, os, re
import pandas as pd

STUDIES = [
    dict(ds="hESC", gse="GSE75748", prj="PRJNA305280", doi="10.1186/s13059-016-1033-x", author="Chu", year=2016, organism="human",
         tech="SMARTer/C1 full-length plate", kw=["embryonic stem", "hESC"]),
    dict(ds="hHep", gse="GSE81252", prj="PRJNA322308", doi="10.1038/nature22796", author="Camp", year=2017, organism="human",
         tech="SMART-seq2 / C1 plate", kw=["liver bud"]),
    dict(ds="mDC", gse="GSE48968", prj="PRJNA212439", doi="10.1038/nature13437", author="Shalek", year=2014, organism="mouse",
         tech="SMART-seq plate", kw=["dendritic", "paracrine"]),
    dict(ds="mESC", gse="GSE98664", prj="PRJNA385832", doi="10.1038/s41467-018-02866-0", author="Hayashi", year=2018, organism="mouse",
         tech="RamDA-seq plate", kw=["RamDA", "recursive splicing"]),
    dict(ds="mHSC", gse="GSE81682", prj="PRJNA322317", doi="10.1182/blood-2016-05-716480", author="Nestorowa", year=2016, organism="mouse",
         tech="SMART-seq2 plate", kw=["blood stem", "hematopoietic stem and progenitor"]),
    dict(ds="K562", gse="(SRA only)", prj="PRJNA831566", doi="10.1016/j.cell.2022.05.013", author="Replogle", year=2022, organism="human (K562 cancer line)",
         tech="10x 3' Perturb-seq", kw=["Perturb-seq"]),
]
RELEASES = ["2023-05-15", "2023-07-25", "2023-12-15"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.expanduser("~/scfm_models/raw"))
    ap.add_argument("--out", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/results/scfm_step2")
    a = ap.parse_args()
    C = json.load(open(a.raw + "/cxg_api/collections.json"))
    rel = {r: set(re.findall(r"cell-census/%s/h5ads/([0-9a-f-]{36})\.h5ad" % r, open(a.raw + f"/cxg_api/h5ads_{r}.xml").read())) for r in RELEASES}
    cur_ids = {d["dataset_id"] for c in C for d in c["datasets"]}
    unresolved = {r: len(rel[r] - cur_ids) for r in RELEASES}
    print("census datasets per release:", {r: len(v) for r, v in rel.items()}, "| not in any current public collection:", unresolved)

    def text(c):
        pm = c.get("publisher_metadata") or {}
        au = " ".join(f"{x.get('family','')}" for x in pm.get("authors", []))
        return " ".join([c.get("name") or "", c.get("description") or "", c.get("doi") or "", au,
                         " ".join(l.get("link_url") or "" for l in c.get("links", [])), pm.get("journal") or ""])

    hits, rows = [], []
    for s in STUDIES:
        found = []
        for c in C:
            t = text(c)
            ev = []
            for key in (s["gse"], s["prj"]):
                if key.startswith(("GSE", "PRJ")) and key in t:
                    ev.append("accession " + key)
            if s["doi"].lower() in (c.get("doi") or "").lower() or s["doi"].lower() in t.lower():
                ev.append("doi")
            pm = c.get("publisher_metadata") or {}
            fam = [x.get("family", "") for x in pm.get("authors", [])[:1]]
            if fam and fam[0] == s["author"] and pm.get("published_year") == s["year"]:
                ev.append("first author+year")
            if ev:
                dsets = [d["dataset_id"] for d in c["datasets"]]
                found.append(dict(study=s["ds"], collection_id=c["collection_id"], name=c["name"], doi=c.get("doi"), evidence="; ".join(ev),
                                  n_datasets=len(dsets), organisms=",".join(sorted({o["label"] for d in c["datasets"] for o in d["organism"]})),
                                  **{f"datasets_in_census_{r}": sum(d in rel[r] for d in dsets) for r in RELEASES}))
        hits += found
        # exact-accession/doi evidence is strong; author+year alone is weak
        strong = [f for f in found if "accession" in f["evidence"] or "doi" in f["evidence"]]
        rows.append(dict(study=s["ds"], accession=s["gse"], bioproject=s["prj"], doi=s["doi"], organism=s["organism"], technology=s["tech"],
                         n_cxg_collections_strong=len(strong), n_cxg_collections_weak=len(found) - len(strong),
                         in_census_2023_05_15=sum(f[f"datasets_in_census_2023-05-15"] for f in strong),
                         in_census_2023_07_25=sum(f[f"datasets_in_census_2023-07-25"] for f in strong),
                         in_census_2023_12_15=sum(f[f"datasets_in_census_2023-12-15"] for f in strong)))
    H = pd.DataFrame(hits)
    S = pd.DataFrame(rows)

    # scBaseCount (SRX keyed)
    sb, sbrows = [], []
    for r in ("2025-02-25", "2026-01-12"):
        for sp in ("Homo_sapiens", "Mus_musculus"):
            m = pd.read_parquet(a.raw + f"/scbasecount_meta/{r}_{sp}_sample_metadata.parquet")
            m["release"], m["species"] = r, sp
            sb.append(m)
    SB = pd.concat(sb)
    lp = SB.groupby(["release", "species", "lib_prep"]).size().rename("n_samples").reset_index()
    for s in STUDIES:
        ri = pd.read_csv(a.raw + f"/sra_runinfo/{s['prj']}.csv")
        ri = ri[ri.Run != "Run"]
        srx = set(ri["Experiment"].dropna())
        for r in ("2025-02-25", "2026-01-12"):
            sub = SB[SB.release == r]
            hit = sub[sub.srx_accession.isin(srx)]
            sbrows.append(dict(study=s["ds"], bioproject=s["prj"], release=r, n_srx_in_bioproject=len(srx), n_srx_in_scbasecount=hit.srx_accession.nunique(),
                               lib_prep_of_hits=";".join(f"{k}:{v}" for k, v in hit.lib_prep.value_counts().items())))
        S.loc[S.study == s["ds"], "n_srx_bioproject"] = len(srx)
    SBR = pd.DataFrame(sbrows)
    for r in ("2025-02-25", "2026-01-12"):
        S[f"scbasecount_{r}_srx_hits"] = [int(SBR[(SBR.study == x) & (SBR.release == r)].n_srx_in_scbasecount.iloc[0]) for x in S.study]

    # verdict table (model x study)
    uns = json.load(open(a.raw + "/cxg_api/unresolved_census_datasets_meta.json"))
    ver = []
    for _, r in S.iterrows():
        human = r.organism.startswith("human")
        cxg_hits = int(r.n_cxg_collections_strong)
        base = ("no CELLxGENE collection cites the accession, BioProject or DOI, and none has this first author and year; "
                f"{sum(unresolved.values())} census dataset-release entries ({len(uns)} datasets) have no current public collection and were identified by title in the h5ad (none is related). "
                "Re-analysis inside an atlas that does not cite the study cannot be excluded from metadata.")
        # scGPT
        if not human:
            v = ("no", "scGPT corpus is human-only (mouse study not in it by construction); not used for mouse cases")
        elif r.study == "K562":
            v = ("no", "scGPT whole-human is normal-cell only and Replogle K562 (cancer line, CRISPRi) is not in any census release: " + base)
        else:
            v = ("no (moderate confidence)", base + " Census snapshot used by scGPT is not stated in reachable sources; all three LTS releases checked.")
        ver.append(dict(model="scGPT whole-human", study=r.study, accession=r.accession, status=v[0], evidence=v[1]))
        ver.append(dict(model="scPRINT v1 medium-v1.5", study=r.study, accession=r.accession, status="no (moderate confidence)" if cxg_hits == 0 else "YES",
                        evidence="corpus = CELLxGENE census LTS 2023-12-15 (Kalfon 2025); " + base))
        sb = SBR[SBR.study == r.study]
        n25 = int(sb[sb.release == "2025-02-25"].n_srx_in_scbasecount.iloc[0]); n26 = int(sb[sb.release == "2026-01-12"].n_srx_in_scbasecount.iloc[0])
        tot = int(sb.n_srx_in_bioproject.iloc[0])
        st = "no" if n25 == 0 else "YES"
        ev = f"scBaseCount SRX overlap with BioProject {r.bioproject}: {n25}/{tot} (2025-02-25), {n26}/{tot} (2026-01-12)"
        if n25 == 0 and n26 > 0:
            st = "no for the 2025-02-25 snapshot; YES in 2026-01-12 (unknown which snapshot scPRINT-2 used)"
        ver.append(dict(model="scPRINT-2 (optional)", study=r.study, accession=r.accession, status=st, evidence=ev))
    pd.DataFrame(ver).to_csv(a.out + "/leakage_audit.csv", index=False)
    os.makedirs(a.out, exist_ok=True)
    H.to_csv(a.out + "/leakage_cxg_collection_hits.csv", index=False)
    S.to_csv(a.out + "/leakage_audit_raw.csv", index=False)
    SBR.to_csv(a.out + "/leakage_scbasecount.csv", index=False)
    lp.to_csv(a.out + "/scbasecount_lib_prep_counts.csv", index=False)
    # positive control of the membership logic: datasets of collections first published before 2023-01-01 should mostly be in the 2023-12-15 census
    old = [d["dataset_id"] for c in C if (c.get("published_at") or "9")[:4] < "2023" for d in c["datasets"]]
    ctrl = {r: round(sum(d in rel[r] for d in old) / max(len(old), 1), 3) for r in RELEASES}
    print("control: fraction of datasets from collections published before 2023 present in each release:", ctrl, "n =", len(old))
    json.dump(dict(control_fraction_pre2023_datasets_in_release=ctrl, control_n=len(old), census_datasets={r: len(v) for r, v in rel.items()}, census_datasets_not_in_current_public_collections=unresolved,
                   n_public_collections=len(C), n_public_datasets=len(cur_ids)), open(a.out + "/leakage_meta.json", "w"), indent=1)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 70):
        print(S.T); print(len(H), "weak/strong CELLxGENE collection hits"); print(SBR); print(lp); print(pd.DataFrame(ver)[["model", "study", "status"]])


if __name__ == "__main__":
    main()
