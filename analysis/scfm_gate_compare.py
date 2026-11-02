#!/usr/bin/env python
"""Step 3.2(b) helper: bitwise comparison of two scPRINT production outputs of the same case (determinism gate). No truth.
  python -I analysis/scfm_gate_compare.py DIR1 DIR2 CASE OUT.json
"""
import json, sys
import numpy as np
d1, d2, case, out = sys.argv[1:5]
a, b = np.load(f"{d1}/{case}.attn.npy"), np.load(f"{d2}/{case}.attn.npy")
n1, n2 = np.load(f"{d1}/{case}.nets.npz"), np.load(f"{d2}/{case}.nets.npz")
R = dict(case=case, attn_bitwise_identical=bool(np.array_equal(a, b)), attn_max_abs_diff=float(np.abs(a - b).max()),
         nets_bitwise_identical=all(bool(np.array_equal(n1[k], n2[k])) for k in ("R0", "R0p", "R2")),
         heads_identical=bool(np.array_equal(n1["r0p_heads"], n2["r0p_heads"])))
R["passed"] = R["attn_bitwise_identical"] and R["nets_bitwise_identical"] and R["heads_identical"]
json.dump(R, open(out, "w"), indent=1)
print(json.dumps(R))
sys.exit(0 if R["passed"] else 1)
