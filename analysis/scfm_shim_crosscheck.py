#!/usr/bin/env python
"""Step 3.2(a) gate: the step-1 shim (scprint 1.0.13, ~/venvs/scfm) must reproduce the production scprint 2.2.0 attention on one
mouse case (mHSC-E_500, arm A, all cells) to max abs difference <= 1e-6 (E3). NO TRUTH FILE IS READ.

  python -I analysis/scfm_shim_crosscheck.py --case mHSC-E_500 --prod ~/scfm_models/step3/_gates/prod --out ~/scfm_models/step3/_gates/shim --device 1
Writes <out>/<case>.shim.attn.npy and shim_vs_prod.json; exits 1 if the gate fails.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd, torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scfm_smoke_scprint as sm            # the frozen step-1 shim: load_scprint(), run()
import scfm_extract_scprint as ex          # production input builder (no scprint import at module level)

ap = argparse.ArgumentParser()
ap.add_argument("--case", default="mHSC-E_500")
ap.add_argument("--prod", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--device", type=int, default=1)
ap.add_argument("--tol", type=float, default=1e-6)
a = ap.parse_args()
assert a.device != 0
os.makedirs(a.out, exist_ok=True)
dev = torch.device(f"cuda:{a.device}")
t0 = time.time()
model, info = sm.load_scprint(dev)
vocab = list(model.genes)
use, x, cells, _ = ex.case_inputs(a.case, "A")
gp = torch.tensor([[vocab.index(e) for e in use.ensembl]], dtype=torch.long)
A, hid, mass = sm.run(model, gp, torch.tensor(x, dtype=torch.float32), dev, 8)
np.save(f"{a.out}/{a.case}.shim.attn.npy", A)
P = np.load(f"{a.prod}/{a.case}.attn.npy")
d = np.abs(A - P)
R = dict(case=a.case, scprint_shim="1.0.13 (step-1 shim)", n_cells=int(x.shape[0]), shape=list(A.shape), max_abs_diff=float(d.max()),
         rel_L1=float(d.sum() / np.abs(P).sum()), bitwise_identical=bool(np.array_equal(A, P)),
         per_layer_max_abs_diff=[float(d[l].max()) for l in range(d.shape[0])], tol=a.tol, passed=bool(d.max() <= a.tol),
         mass_max_abs_diff=float(np.abs(mass - json.load(open(f"{a.prod}/{a.case}.json"))["special_token_mass"]).max()),
         shim_info={k: info[k] for k in ("dropped_hparams", "missing_keys", "unexpected_keys")}, total_s=time.time() - t0)
json.dump(R, open(f"{a.out}/shim_vs_prod.json", "w"), indent=1)
print(json.dumps({k: R[k] for k in ("case", "n_cells", "max_abs_diff", "rel_L1", "bitwise_identical", "passed")}))
sys.exit(0 if R["passed"] else 1)
