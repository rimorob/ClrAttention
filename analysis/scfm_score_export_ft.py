#!/usr/bin/env python
"""Scoring step 0 for the fine-tuning arm (E17/E18): verify the frozen fine-tuned extraction against its SHA-256 manifest
(results/scfm_ft/sha256_manifest.csv, hash recorded in the change log) and export the network matrices to flat float32 files
(column-major, as analysis/scfm_score_export.py) that analysis/scfm_score_ft.R reads. NO TRUTH FILE IS READ HERE.
Any mismatch: exit 2 and nothing is exported.

  python -I analysis/scfm_score_export_ft.py --manifest-sha <sha> [--out ~/scfm_models/step5_ft/scoring] [--check-only]
Sources: scprint_ft_A (<- scprint_attn: R0 R0p R2), scgpt_ft_A (<- scgpt_attn: attn_last_mean tok_cos tok_cos_ln hid_cos).
"""
import argparse, csv, hashlib, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT = os.path.expanduser("~/scfm_models/step5_ft")
SOURCES = {"scprint_ft_A": ("scprint_attn", ".nets.npz", ["R0", "R0p", "R2"]),
           "scgpt_ft_A": ("scgpt_attn", ".nets.npz", ["attn_last_mean", "tok_cos", "tok_cos_ln", "hid_cos"])}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest-sha", required=True)
    ap.add_argument("--out", default="~/scfm_models/step5_ft/scoring")
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()
    mf = f"{ROOT}/results/scfm_ft/sha256_manifest.csv"
    if sha256(mf) != a.manifest_sha:
        print("MANIFEST FILE HASH MISMATCH", sha256(mf), "expected", a.manifest_sha)
        sys.exit(2)
    rows = list(csv.DictReader(open(mf)))
    bad = [r["path"] for r in rows if not os.path.exists(f"{FT}/{r['path']}") or os.path.getsize(f"{FT}/{r['path']}") != int(r["bytes"]) or sha256(f"{FT}/{r['path']}") != r["sha256"]]
    print(f"manifest: {len(rows)} files verified, {len(bad)} mismatches", flush=True)
    if bad:
        print("MISMATCH:", *bad, sep="\n  ")
        sys.exit(2)
    if a.check_only:
        return
    man = {r["path"]: r["sha256"] for r in rows}
    out = os.path.expanduser(a.out)
    index = []
    for src, (d, suf, keys) in SOURCES.items():
        cases = sorted({os.path.basename(p)[:-len(suf)] for p in man if p.startswith(d + "/") and p.endswith(suf)})
        if not cases:
            continue
        os.makedirs(f"{out}/{src}", exist_ok=True)
        for case in cases:
            rel = f"{d}/{case}{suf}"
            z = np.load(f"{FT}/{rel}", allow_pickle=True)
            genes = [str(g) for g in z["genes"]]
            gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t")
            exp = list(gm[gm.in_scgpt.values & gm.scgpt_token.notna().values].symbol) if src == "scgpt_ft_A" else list(gm[gm.status == "ok"].symbol)
            assert genes == exp, f"{rel}: gene order differs from the gene map"
            G = len(genes)
            with open(f"{out}/{src}/{case}.genes.txt", "w") as f:
                f.write("\n".join(genes) + "\n")
            for k in keys:
                m = np.asarray(z[k], dtype=np.float32)
                assert m.shape == (G, G) and np.isfinite(m).all(), (rel, k, m.shape)
                p = f"{out}/{src}/{case}.{k}.f32"
                with open(p, "wb") as f:
                    f.write(np.asfortranarray(m).tobytes(order="F"))
                index.append(dict(source=src, case=case, readout=k, G=G, bytes=os.path.getsize(p), sha256=sha256(p), from_file=rel, from_sha256=man[rel]))
    pd.DataFrame(index).to_csv(f"{out}/export_index.csv", index=False)
    print(f"exported {len(index)} matrices to {out}")


if __name__ == "__main__":
    main()
