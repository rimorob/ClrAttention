#!/usr/bin/env python
"""Scoring step 0 (2026-11-02): verify the frozen extraction against the E13 manifest and export the network matrices to flat
float32 files that analysis/scfm_score.R reads. NO TRUTH FILE IS READ HERE.

Stop rule (E13): the SHA-256 of results/scfm_step3/sha256_manifest.csv must equal the value in the change log, and every one of the 239
network/attention files under ~/scfm_models/step3 must match its manifest hash. Any mismatch: exit 2 and nothing is exported.

Output (--out, default ~/scfm_models/step4_scoring), one directory per source, files per case:
  <source>/<case>.<readout>.f32    G x G float32, COLUMN-major (so R's matrix(readBin(...), G, G) is exactly the python array),
                                   rows = query (A convention: A[query, key]); the orientation S = t(A) is applied at scoring
  <source>/<case>.genes.txt        gene symbols in the same order (the gene-map order of results/scfm_step2/gene_map/<case>.tsv, status ok)
  export_index.csv                 source, case, readout, G, bytes, sha256 of every exported file
Readouts: scprint_medium_{A,B}: R0 R0p R2; r1_{A,B}: R1; ctrl_*: R0; scprint_{large,small}_A: R0; esm2_cos: cos;
scgpt_A: attn_last_mean tok_cos tok_cos_ln hid_cos (the seed-1 stability rerun is not exported).

  python -I analysis/scfm_score_export.py [--out DIR] [--check-only]
"""
import argparse, csv, hashlib, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEP3 = os.path.expanduser("~/scfm_models/step3")
MANIFEST_SHA = "e27a8422f5cf868459c694407a0cd0e900a18ce8b3aaf0c33ee179260c2e3f15"   # change log E13
CASES = [f"{d}_{n}" for d in ("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L") for n in (500, 1000)]

# source directory -> (file suffix, [npz keys to export], cases)
SOURCES = {}
for s in ("scprint_medium_A", "scprint_medium_B"):
    SOURCES[s] = (".nets.npz", ["R0", "R0p", "R2"])
for s in ("r1_A", "r1_B"):
    SOURCES[s] = (".r1.npz", ["R1"])
for k in ("ctrl_perm_s1001", "ctrl_perm_s1002", "ctrl_perm_s1003", "ctrl_rand_esm_s1", "ctrl_rand_esm_s2", "ctrl_rand_esm_s3",
          "ctrl_rand_emb_s1", "ctrl_rand_emb_s2", "ctrl_rand_emb_s3", "scprint_large_A", "scprint_small_A"):
    SOURCES[k] = (".nets.npz", ["R0"])
SOURCES["esm2_cos"] = (".esm2cos.npz", ["cos"])
SOURCES["scgpt_A"] = (".nets.npz", ["attn_last_mean", "tok_cos", "tok_cos_ln", "hid_cos"])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def verify(mf=None, expect=None):
    mf = mf or f"{ROOT}/results/scfm_step3/sha256_manifest.csv"
    expect = expect or MANIFEST_SHA
    if sha256(mf) != expect:
        print("MANIFEST FILE HASH MISMATCH", sha256(mf), "expected", expect)
        sys.exit(2)
    rows = list(csv.DictReader(open(mf)))
    bad = []
    for r in rows:
        p = os.path.join(STEP3, r["path"])
        if not os.path.exists(p) or os.path.getsize(p) != int(r["bytes"]) or sha256(p) != r["sha256"]:
            bad.append(r["path"])
    print(f"manifest: {len(rows)} files verified, {len(bad)} mismatches", flush=True)
    if bad:
        print("MISMATCH:", *bad, sep="\n  ")
        sys.exit(2)
    return {r["path"]: r["sha256"] for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="~/scfm_models/step4_scoring")
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--k562", action="store_true", help="export the K562 extraction (~/scfm_models/step3/k562) frozen in results/scfm_k562/sha256_manifest.csv")
    ap.add_argument("--manifest-sha", default="", help="with --k562: SHA-256 of that manifest, as recorded in the change log")
    a = ap.parse_args()
    global STEP3
    if a.k562:
        assert a.manifest_sha, "--manifest-sha required with --k562"
        STEP3 = os.path.expanduser("~/scfm_models/step3/k562")
        for k in ("scprint_medium_B", "r1_B", "esm2_cos"):
            SOURCES.pop(k, None)
        man = verify(f"{ROOT}/results/scfm_k562/sha256_manifest.csv", a.manifest_sha)
    else:
        man = verify()
    if a.check_only:
        return
    out = os.path.expanduser(a.out)
    index = []
    for src, (suf, keys) in SOURCES.items():
        os.makedirs(f"{out}/{src}", exist_ok=True)
        cases = sorted({os.path.basename(p).split(suf)[0] for p in man if p.startswith(src + "/") and p.endswith(suf)})
        if src == "scgpt_A":
            cases = [c for c in cases if "seed" not in c]
        for case in cases:
            rel = f"{src}/{case}{suf}"
            assert rel in man, rel
            z = np.load(f"{STEP3}/{rel}", allow_pickle=True)         # string object arrays (E13 note d)
            genes = [str(g) for g in z["genes"]]
            gm = pd.read_csv(f"{ROOT}/results/scfm_step2/gene_map/{case}.tsv", sep="\t")
            if src == "scgpt_A":
                exp = list(gm[gm.in_scgpt.values & gm.scgpt_token.notna().values].symbol)
            elif src.endswith("_B") or src.startswith("r1_B"):
                exp = None            # arm B drops genes without deposited counts: checked as a subsequence of the arm-A order
            else:
                exp = list(gm[gm.status == "ok"].symbol)
            if exp is not None:
                assert genes == exp, f"{rel}: gene order differs from the gene map"
            else:
                ok = list(gm[gm.status == "ok"].symbol)
                it = iter(ok)
                assert all(g in it for g in genes), f"{rel}: arm B genes are not an ordered subset of the arm A genes"
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
