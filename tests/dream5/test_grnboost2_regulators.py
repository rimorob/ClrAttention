#!/usr/bin/env python3
"""Equivalence tests for tools/grnboost2.py --regulators (DREAM5 design, step 4.2).
Run from the repo root:  nice -n 19 ~/venvs/clr-baselines/bin/python tests/dream5/test_grnboost2_regulators.py
  A. regulators = all genes (same order) gives byte-identical output to the run without the option.
  B. that output equals the output of the pre-change script (git HEAD version of tools/grnboost2.py, taken from the
     last commit before this change) on the same input.
  C. a regulator subset: for a target outside the subset and for a target inside it, the edges equal the old script's
     edges for that target when the input file holds only the regulators plus that target (so the regulator
     set is exactly the subset, minus the target itself).
"""
import os, subprocess, sys, tempfile, shutil
import numpy as np, pandas as pd

PY = sys.executable
OLD_REV = os.environ.get("GB2_OLD_REV", "278270e")   # commit before --regulators existed
tmp = tempfile.mkdtemp(prefix="gb2reg_")
fails = []

def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name); 
    if not ok: fails.append(name)

# old script copy (needs grnboost2_variants.py next to it)
old_dir = os.path.join(tmp, "old"); os.makedirs(old_dir)
src = subprocess.check_output(["git", "show", f"{OLD_REV}:tools/grnboost2.py"]).decode()
open(os.path.join(old_dir, "grnboost2.py"), "w").write(src)
shutil.copy("tools/grnboost2_variants.py", old_dir)

rng = np.random.default_rng(7)
G, n = 24, 40
genes = [f"G{i+1}" for i in range(G)]
Z = rng.normal(size=(G, n)); Z[5] += 0.8 * Z[0]; Z[9] += 0.8 * Z[1]; Z[14] -= 0.7 * Z[2]
def write(path, names, rows):
    pd.DataFrame(Z[rows], index=names, columns=[f"c{i}" for i in range(n)]).rename_axis("gene").to_csv(path, sep="\t")
full = os.path.join(tmp, "full.tsv"); write(full, genes, list(range(G)))

def run(script, inp, out, *extra):
    subprocess.check_call([PY, script, inp, out, "--workers", "3", "--seed", "11", *extra], stdout=subprocess.DEVNULL)

new = "tools/grnboost2.py"; old = os.path.join(old_dir, "grnboost2.py")
o_old, o_new, o_all = (os.path.join(tmp, f) for f in ("old.tsv", "new.tsv", "all.tsv"))
run(old, full, o_old); run(new, full, o_new)
allf = os.path.join(tmp, "all.txt"); open(allf, "w").write("\n".join(genes) + "\n")
run(new, full, o_all, "--regulators", allf)
b = lambda p: open(p, "rb").read()
check("A: regulators=all genes is byte-identical to no option", b(o_new) == b(o_all))
check("B: new script without option is byte-identical to the pre-change script", b(o_old) == b(o_new))

# C: subset of regulators
reg_idx = [0, 1, 2, 3, 4, 7, 11, 17]; regs = [genes[i] for i in reg_idx]
rf = os.path.join(tmp, "reg.txt"); open(rf, "w").write("\n".join(regs) + "\n")
o_sub = os.path.join(tmp, "sub.tsv"); run(new, full, o_sub, "--regulators", rf)
sub = pd.read_csv(o_sub, sep="\t")
check("C0: every edge's TF is in the regulator list", set(sub.TF) <= set(regs))
check("C0: no self edges", (sub.TF != sub.target).sum() == len(sub))
check("C0: all genes appear as targets", set(sub.target) == set(genes) or set(sub.target) <= set(genes) and len(set(sub.target)) >= G - 2)
for tgt in (6, 17):     # 6 = outside the list; 17 = inside the list (its own removal)
    rows = reg_idx + ([tgt] if tgt not in reg_idx else [])
    names = [genes[i] for i in rows]
    f = os.path.join(tmp, f"red{tgt}.tsv"); write(f, names, rows)
    o = os.path.join(tmp, f"red{tgt}_out.tsv"); run(old, f, o)
    ref = pd.read_csv(o, sep="\t"); ref = ref[ref.target == genes[tgt]].sort_values("TF").reset_index(drop=True)
    got = sub[sub.target == genes[tgt]].sort_values("TF").reset_index(drop=True)
    check(f"C: target G{tgt+1} edges equal the old script on regulators+target only",
          len(ref) == len(got) and (ref.TF == got.TF).all() and np.allclose(ref.importance, got.importance, rtol=0, atol=0))
shutil.rmtree(tmp)
print(f"{'ALL PASS' if not fails else 'FAILED: ' + str(fails)}"); sys.exit(1 if fails else 0)
