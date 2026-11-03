#!/usr/bin/env python
"""E17/E18 freeze: SHA-256 manifest of every network/attention file (.npz, .npy) of the fine-tuning arm under ~/scfm_models/step5_ft/<arm>
(paths relative to ~/scfm_models/step5_ft), written to results/scfm_ft/sha256_manifest.csv (columns path,bytes,sha256, sorted). NO TRUTH FILE IS READ.
  python -I analysis/scfm_ft_manifest.py [arm_dir ...]      default arms: scprint_attn scgpt_attn (those that exist)
"""
import csv, hashlib, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.expanduser("~/scfm_models/step5_ft")
arms = sys.argv[1:] or [a for a in ("scprint_attn", "scgpt_attn") if os.path.isdir(f"{D}/{a}")]
rows = []
for arm in arms:
    for dp, _, fn in os.walk(f"{D}/{arm}"):
        for f in fn:
            if f.endswith((".npz", ".npy")):
                p = os.path.join(dp, f)
                h = hashlib.sha256()
                with open(p, "rb") as fh:
                    for b in iter(lambda: fh.read(1 << 24), b""):
                        h.update(b)
                rows.append((os.path.relpath(p, D), os.path.getsize(p), h.hexdigest()))
rows.sort()
os.makedirs(f"{ROOT}/results/scfm_ft", exist_ok=True)
with open(f"{ROOT}/results/scfm_ft/sha256_manifest.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["path", "bytes", "sha256"])
    w.writerows(rows)
print(len(rows), "files in", arms)
print("manifest sha256", hashlib.sha256(open(f"{ROOT}/results/scfm_ft/sha256_manifest.csv", "rb").read()).hexdigest())
