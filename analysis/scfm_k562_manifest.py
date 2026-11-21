#!/usr/bin/env python
"""K562 freeze: SHA-256 manifest of every network/attention file (.npz, .npy) under ~/scfm_models/step3/k562 (paths relative to that directory),
written to results/scfm_k562/sha256_manifest.csv (columns path,bytes,sha256, sorted). NO TRUTH FILE IS READ.
  python -I analysis/scfm_k562_manifest.py [--dir ~/scfm_models/step3/k562] [--out results/scfm_k562/sha256_manifest.csv]
(--dir/--out select another dataset, e.g. --dir ~/scfm_models/step3/hepg2 --out results/scfm_hepg2/sha256_manifest.csv; defaults = the K562 freeze.)
"""
import argparse, csv, hashlib, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="~/scfm_models/step3/k562")
ap.add_argument("--out", default=f"{ROOT}/results/scfm_k562/sha256_manifest.csv")
_a = ap.parse_args()
D = os.path.expanduser(_a.dir)
OUTF = os.path.abspath(_a.out)
rows = []
for dp, _, fn in os.walk(D):
    for f in fn:
        if f.endswith((".npz", ".npy")):
            p = os.path.join(dp, f)
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for b in iter(lambda: fh.read(1 << 24), b""):
                    h.update(b)
            rows.append((os.path.relpath(p, D), os.path.getsize(p), h.hexdigest()))
rows.sort()
os.makedirs(os.path.dirname(OUTF), exist_ok=True)
with open(OUTF, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["path", "bytes", "sha256"])
    w.writerows(rows)
print(len(rows), "files")
h = hashlib.sha256(open(OUTF, "rb").read()).hexdigest()
print("manifest sha256", h)
