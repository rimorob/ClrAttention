#!/usr/bin/env python3
"""Label-free input facts for the DREAM5 design (docs/plans/dream5_design.md).

Reads only the training data (data/dream5/synapse/net{1,3,4}, data/dream5/synapse/original,
and the organisers' login-free archive). Never touches data/dream5/synapse/scoring/ or any
gold standard. Output: stdout (saved to docs/dream5_selections/input_facts.txt).
"""
import collections
import numpy as np
import pandas as pd

R = "data/dream5/"
S = R + "synapse/"
P = R + "public_gnw/DREAM5_NetworkInferenceChallenge_AlternativeDataFormats/"
ORIG = {3: ("Ecoli/ecoli_gene_names.tsv", "Ecoli/ecoli_tf_names.tsv"),
        4: ("Scerevisiae/yeast_gene_names.tsv", "Scerevisiae/yeast_tf_names.tsv")}


def lines(p):
    return [l.strip() for l in open(p) if l.strip()]


for n in (1, 3, 4):
    d = f"{S}net{n}/net{n}_"
    X = pd.read_csv(d + "expression_data.tsv", sep="\t")
    A = X.values.astype(float)
    cf = pd.read_csv(d + "chip_features.tsv", sep="\t", keep_default_na=False, dtype=str)
    gid = pd.read_csv(d + "gene_ids.tsv", sep="\t", keep_default_na=False)
    gid.columns = ["id", "name"]
    tf = lines(d + "transcription_factors.tsv")
    nm = dict(zip(gid.id, gid.name))
    nt, ng = len(tf), X.shape[1]
    print(f"=== net{n}")
    print(f"matrix chips x genes {A.shape}; NaN {int(np.isnan(A).sum())} (fraction {np.isnan(A).mean():.3g}); "
          f"inf {int(np.isinf(A).sum())}; range {A.min():.4g}..{A.max():.4g}; mean {A.mean():.4g}; "
          f"zero-variance genes {int((A.std(0) == 0).sum())}; duplicate gene columns {int(X.columns.duplicated().sum())}")
    print(f"gene_ids rows {len(gid)}; ids equal expression header: {list(gid.id) == list(X.columns)}")
    print(f"TFs {nt}; all in gene list: {all(t in nm for t in tf)}; TF ids are G1..G{nt}: "
          f"{tf == [f'G{i}' for i in range(1, nt + 1)]}; TF fraction of genes {nt / ng:.3f}")
    dec = gid.name.str.startswith("decoy")
    dtf = sum(nm[t].startswith("decoy") for t in tf)
    print(f"decoy genes (random expression, name 'decoy*') {int(dec.sum())}; decoy TFs {dtf}; "
          f"real genes {ng - int(dec.sum())}; real TFs {nt - dtf}")
    print(f"candidate TF x gene pairs {nt * ng}; without self-pairs {nt * ng - nt}")
    if n in ORIG:
        og, ot = lines(S + "original/" + ORIG[n][0]), lines(S + "original/" + ORIG[n][1])
        real = set(gid.name[~dec]); tfn = {nm[t] for t in tf}
        print(f"original (pre-anonymisation) genes {len(og)}, all in challenge: {set(og) <= real}, "
              f"challenge real genes not in original {len(real - set(og))}; original TFs {len(ot)}, "
              f"in challenge TF list {len(set(ot) & tfn)}; TF-list entries not in original TF list "
              f"(non-decoy) {len({x for x in tfn if not x.startswith('decoy')} - set(ot))}")
    r = cf["Repeat"].astype(int).values
    grp = np.cumsum(r == 1)
    sizes = pd.Series(grp).value_counts()
    ser = cf["#Experiment"].astype(int).values
    per = pd.Series(grp).groupby(ser).nunique()
    print(f"chip features columns {list(cf.columns)}; rows {len(cf)}")
    print(f"'#Experiment' = coarse series index: {len(set(ser))} values (non-decreasing {bool((np.diff(ser) >= 0).all())}); "
          f"replicate group = run starting at Repeat==1: {grp.max()} groups; Repeat is 1..k within every group: "
          f"{all(list(r[grp == k]) == list(range(1, int((grp == k).sum()) + 1)) for k in range(1, grp.max() + 1))}")
    print(f"chips per replicate group {dict(sorted(sizes.value_counts().to_dict().items()))}; mean {sizes.mean():.3f}; "
          f"chips in multi-chip groups {int((sizes[sizes > 1]).sum())} ({(sizes[sizes > 1]).sum() / len(cf):.3f})")
    print(f"replicate groups per series: median {per.median():.0f}, max {per.max()}; series with >1 group {int((per > 1).sum())}")
    for c in cf.columns[1:7]:
        v = cf[c]; k = v[v != "NA"]
        print(f"  {c}: non-NA chips {len(k)}, distinct {k.nunique()}")
    for col in ("DeletedGenes", "OverexpressedGenes"):
        rows, gs = 0, set()
        for v in cf[col]:
            g = [x for x in v.split(",") if x.startswith("G") and int(x[1:]) <= nt] if v != "NA" else []
            if g:
                rows += 1; gs |= set(g)
        print(f"  {col}: chips naming a TF-list gene {rows}, distinct TFs {len(gs)}")
    Xa = pd.read_csv(f"{P}net{n}/net{n}_expression_data_avg.tsv", sep="\t")
    ca = pd.read_csv(f"{P}net{n}/net{n}_chip_features_avg.tsv", sep="\t", keep_default_na=False, dtype=str)
    m = pd.DataFrame(A).groupby(grp).mean().values
    Xt = pd.read_csv(f"{P}net{n}/net{n}_expression_data_t.tsv", sep="\t", index_col=0)
    print(f"averaged file {Xa.shape}; Repeat sum {int(ca['Repeat'].astype(int).sum())}; max |avg - mean of replicate chips| "
          f"{np.abs(m - Xa.values).max():.2g}; transposed archive == Synapse matrix: {np.array_equal(Xt.values.T.astype(float), A)}")
