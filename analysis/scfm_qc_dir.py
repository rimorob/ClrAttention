#!/usr/bin/env python
"""Step 3 truth-free QC gate for one output directory of scfm_extract_scprint.py / scfm_extract_scgpt.py / scfm_extract_r1.py.
Exits 1 (so that dependent scheduler jobs are held) if any case has NaN/Inf, a row-sum outside [0.9999, 1.0001] (R0 files), a missing
case, or an error. Prints one line per case. No truth file.
  python -I analysis/scfm_qc_dir.py DIR NCASES [r0|r1|scgpt]
"""
import glob, json, os, sys
d, ncases = sys.argv[1], int(sys.argv[2])
kind = sys.argv[3] if len(sys.argv) > 3 else "r0"
pat = {"r0": "*.json", "r1": "*.r1.json", "scgpt": "*.json"}[kind]
bad = []
files = [f for f in sorted(glob.glob(os.path.join(d, pat))) if not os.path.basename(f).startswith(("run", "reproduce"))]
for f in files:
    m = json.load(open(f))
    c = m.get("case")
    if "error" in m and m["error"]:
        bad.append((c, "error")); continue
    if kind == "r0":
        ok = m["attn_nan"] == 0 and m["attn_inf"] == 0 and 0.9999 < m["row_sum_min"] and m["row_sum_max"] < 1.0001
        extra = f"orgacc={m.get('organism_accuracy')} r0p={m['r0p_n_heads']} splitR0={m['split_half']['R0']['relL1']:.3f} infer={m['infer_s']:.0f}s"
    elif kind == "r1":
        ok = m["nan"] == 0 and m["inf"] == 0 and m["case_genes_missing_from_output"] == 0
        extra = f"rowsum[{m['row_sum_min']:.3f},{m['row_sum_max']:.3f}] t={m['gninfer_s']:.0f}s"
    else:
        ok = m["attn_nan"] == 0 and m["attn_inf"] == 0 and 0.9999 < m["row_sum_min"] and m["row_sum_max"] < 1.0001 and m["hid_cos_nan"] == 0
        extra = f"clsmass_last_max={max(m['cls_token_mass'][-1]):.4f} infer={m['infer_s']:.0f}s"
    print(f"{os.path.basename(d)}/{c}: {'OK' if ok else 'FAIL'} n={m['n_cells']} G={m['n_genes']} {extra}")
    if not ok:
        bad.append((c, "check"))
if len(files) != ncases:
    bad.append(("count", f"{len(files)} != {ncases}"))
print("QC", "PASS" if not bad else f"FAIL {bad}")
sys.exit(1 if bad else 0)
