#!/usr/bin/env python
"""K562 freeze: SHA-256 manifest of every network/attention file (.npz, .npy) under ~/scfm_models/step3/k562 (paths relative to that directory),
written to results/scfm_k562/sha256_manifest.csv (columns path,bytes,sha256, sorted). NO TRUTH FILE IS READ.
  python -I analysis/scfm_k562_manifest.py
"""
import csv, hashlib, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.expanduser("~/scfm_models/step3/k562")
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
with open(f"{ROOT}/results/scfm_k562/sha256_manifest.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["path", "bytes", "sha256"])
    w.writerows(rows)
print(len(rows), "files")
h = hashlib.sha256(open(f"{ROOT}/results/scfm_k562/sha256_manifest.csv", "rb").read()).hexdigest()
print("manifest sha256", h)
