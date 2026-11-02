#!/usr/bin/env python
"""Step 3 truth-free QC summary and SHA-256 manifest of every network / attention file. NO TRUTH FILE IS READ.
  python -I analysis/scfm_step3_summary.py --step3 ~/scfm_models/step3 --out results/scfm_step3
Writes (in --out): sha256_manifest.csv, qc_scprint_A.csv, r0p_head_sets.csv, a_vs_b.csv, r1_vs_r0.csv, scgpt_qc.csv, sink_mass_A.csv, summary.json
"""
import argparse, glob, hashlib, json, os
import numpy as np, pandas as pd
from scipy.stats import rankdata

CASES = [f"{d}_{n}" for d in ("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L") for n in (500, 1000)]
BC = [c for c in CASES if c.split("_")[0] in ("hESC", "mHSC-E", "mHSC-GM", "mHSC-L")]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def offd(a, b):
    G = a.shape[0]
    m = ~np.eye(G, dtype=bool)
    x, y = a[m].astype(np.float64), b[m].astype(np.float64)
    return float(np.abs(x - y).sum() / (0.5 * (x.sum() + y.sum()))), float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step3", default=os.path.expanduser("~/scfm_models/step3"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S, O = os.path.expanduser(a.step3), a.out
    os.makedirs(O, exist_ok=True)
    # 1. manifest
    rows = []
    for d in sorted(glob.glob(f"{S}/*")):
        if os.path.basename(d).startswith("_"):
            continue
        for f in sorted(glob.glob(f"{d}/*.np[yz]")):
            rows.append(dict(path=os.path.relpath(f, S), bytes=os.path.getsize(f), sha256=sha(f)))
    pd.DataFrame(rows).to_csv(f"{O}/sha256_manifest.csv", index=False)
    # 2. scPRINT A / B QC
    qa, hs_rows, mass_rows = [], [], []
    J = lambda arm, c: json.load(open(f"{S}/scprint_medium_{arm}/{c}.json"))
    for c in CASES:
        m = J("A", c)
        r = dict(case=c, n_cells=m["n_cells"], G=m["n_genes"], org_acc_A=m["organism_accuracy"], infer_s=round(m["infer_s"], 1),
                 r0p_n_heads=m["r0p_n_heads"], r0p_rule=m["r0p_rule"])
        for k in ("R0", "R0p", "R2"):
            r[f"split_{k}_relL1"], r[f"split_{k}_spearman"] = m["split_half"][k]["relL1"], m["split_half"][k]["spearman"]
        mass = np.array(m["special_token_mass"])
        r.update(mass_min=mass.min(), mass_median=float(np.median(mass)), mass_max=mass.max(), heads_mass_gt_0p9=int((mass > 0.9).sum()),
                 heads_mass_lt_0p5=int((mass < 0.5).sum()))
        if c in BC:
            mb = J("B", c)
            r["org_acc_B"] = mb["organism_accuracy"]
        qa.append(r)
        hs_rows.append(dict(case=c, n_heads=m["r0p_n_heads"], rule=m["r0p_rule"], heads_layer_head_1based=";".join(f"L{l + 1}H{h + 1}" for l, h in m["r0p_heads"])))
        for l in range(mass.shape[0]):
            for h in range(mass.shape[1]):
                mass_rows.append(dict(case=c, layer=l + 1, head=h + 1, special_mass=mass[l, h], rel_entropy=m["rel_entropy"][l][h]))
    pd.DataFrame(qa).to_csv(f"{O}/qc_scprint_A.csv", index=False)
    pd.DataFrame(hs_rows).to_csv(f"{O}/r0p_head_sets.csv", index=False)
    pd.DataFrame(mass_rows).to_csv(f"{O}/sink_mass_A.csv", index=False)
    # 3. A vs B (common genes, R0)
    ab = []
    for c in BC:
        na, nb = np.load(f"{S}/scprint_medium_A/{c}.nets.npz", allow_pickle=True), np.load(f"{S}/scprint_medium_B/{c}.nets.npz", allow_pickle=True)
        ea, eb = list(na["ensembl"]), list(nb["ensembl"])
        com = [e for e in ea if e in set(eb)]
        ia, ib = [ea.index(e) for e in com], [eb.index(e) for e in com]
        r = dict(case=c, G_A=len(ea), G_B=len(eb), G_common=len(com))
        for k in ("R0", "R0p", "R2"):
            A_, B_ = na[k][np.ix_(ia, ia)], nb[k][np.ix_(ib, ib)]
            # renormalise rows over the common genes only if G differs (it differs by at most one gene)
            r[f"{k}_relL1"], r[f"{k}_spearman"] = offd(A_, B_)
        ab.append(r)
    pd.DataFrame(ab).to_csv(f"{O}/a_vs_b.csv", index=False)
    # 4. R1 vs R0
    rr = []
    for c in CASES:
        r1 = np.load(f"{S}/r1_A/{c}.r1.npz", allow_pickle=True)["R1"]
        n0 = np.load(f"{S}/scprint_medium_A/{c}.nets.npz", allow_pickle=True)
        assert list(np.load(f"{S}/r1_A/{c}.r1.npz", allow_pickle=True)["ensembl"]) == list(n0["ensembl"])
        r = dict(case=c, G=r1.shape[0])
        for k in ("R0", "R0p", "R2"):
            r[f"R1_vs_{k}_spearman"] = offd(r1, n0[k])[1]
        rr.append(r)
    pd.DataFrame(rr).to_csv(f"{O}/r1_vs_r0.csv", index=False)
    # 5. scGPT
    sg = []
    for c in ("hESC_500", "hESC_1000", "hHep_500", "hHep_1000"):
        m = json.load(open(f"{S}/scgpt_A/{c}.json"))
        r = dict(case=c, n_cells=m["n_cells"], G=m["n_genes"], seed=m["binning_seed"], infer_s=round(m["infer_s"], 1),
                 cls_mass_last_layer_max=max(m["cls_token_mass"][-1]), cls_mass_max_any_layer=float(np.max(m["cls_token_mass"])))
        z = np.load(f"{S}/scgpt_A/{c}.nets.npz", allow_pickle=True)
        r["tok_cos_vs_tok_cos_ln_spearman"] = offd(z["tok_cos"], z["tok_cos_ln"])[1]
        r["attn_last_vs_hid_cos_spearman"] = offd(z["attn_last_mean"], z["hid_cos"])[1]
        if c == "hESC_1000" and os.path.exists(f"{S}/scgpt_A/hESC_1000_seed1.nets.npz"):
            z1 = np.load(f"{S}/scgpt_A/hESC_1000_seed1.nets.npz", allow_pickle=True)
            for k in ("attn_last_mean", "attn_all_mean", "hid_cos"):
                r[f"seed1_vs_prod_{k}_relL1"], r[f"seed1_vs_prod_{k}_spearman"] = offd(z[k], z1[k])
        sg.append(r)
    pd.DataFrame(sg).to_csv(f"{O}/scgpt_qc.csv", index=False)
    print(pd.DataFrame(qa).round(3).to_string())
    print(pd.DataFrame(ab).round(3).to_string())
    print(pd.DataFrame(rr).round(3).to_string())
    print(pd.DataFrame(sg).round(4).to_string())
    print(len(rows), "files hashed")


if __name__ == "__main__":
    main()
