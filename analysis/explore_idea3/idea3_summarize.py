#!/usr/bin/env python3
"""EXPLORATORY (docs/plans/explore_idea3_plan.md). Tables (markdown) and the criterion-curve figure for one dataset.
Usage: python3 analysis/explore_idea3/idea3_summarize.py m3d|k562
Writes results/explore_idea3/<ds>/tables.md and criterion_curves_per_decile.png. Statistics only; no interpretation."""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ds = sys.argv[1]
d = sys.argv[2] if len(sys.argv) > 2 else f"results/explore_idea3/{ds}"
GRID = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
out = []
def md(df, floatfmt=None):
    cols = list(df.columns)
    s = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append("" if np.isnan(v) else (floatfmt.get(c, "{:.4g}") if floatfmt else "{:.4g}").format(v))
            else:
                cells.append(str(v))
        s += "| " + " | ".join(cells) + " |\n"
    return s
def add(title, df, fmt=None):
    out.append(f"### {title}\n\n" + md(df, fmt) + "\n")

for f in ("equivalence_criterion.txt", "equivalence_scoring.txt"):
    p = f"{d}/{f}"
    if os.path.exists(p):
        out.append(f"### Equivalence: {f}\n\n```\n" + open(p).read().strip() + "\n```\n")

kn = f"{d}/knee_null_summary.csv"
if os.path.exists(kn):
    add("Knee null (step 0): between-decile statistics of log2 of the decile-median knee, real vs rewired operators (50) and vs label permutations (2,000)", pd.read_csv(kn))
    add("Per-decile knees (median of per-gene knee, log2; snapped grid depth) and the rewired-operator values", pd.read_csv(f"{d}/knee_decile.csv"))

sel = pd.read_csv(f"{d}/selected_depths.csv")
for ls in ["mi", "pi"]:
    s = sel[(sel.op == "real") & (sel.labelset == ls)].sort_values("decile")
    add(f"Selected depth per decile, stratifier {ls} (op real; mean criterion over {int(s.n_curves.iloc[0])} split-directions)",
        s[["decile", "t_d", "t_argmax", "t_lo", "t_hi", "width_octaves", "crit_max", "se_at_max", "splitdir_argmax_median", "splitdir_argmax_q25", "splitdir_argmax_q75"]])
s = sel[(sel.op == "rewired")].sort_values("decile")
if len(s):
    add("Selected depth per decile on the REWIRED operator (stratifier mi)", s[["decile", "t_d", "t_argmax", "t_lo", "t_hi", "width_octaves", "crit_max", "se_at_max"]])
add("Selection statistics (SD of log2 t_d over deciles, range, trend, deciles >= 1 octave from t = 9, deciles differing from the global-grid depth)", pd.read_csv(f"{d}/selection_stats.csv"))
add("Controls: parameters chosen by the same criterion (global grid depth, PPR alpha, heat-kernel tau)", pd.read_csv(f"{d}/selected_controls.csv"))

if os.path.exists(f"{d}/scores.csv"):
    sc = pd.read_csv(f"{d}/scores.csv"); con = pd.read_csv(f"{d}/contrasts.csv")
    order = ["t9", "global_grid", "form1_mi", "form1_pi", "form2_mi", "form2_pi", "form2raw", "form1_perm1", "ppr", "hk"]
    for e in sc.truth.unique():
        t = sc[sc.truth == e].pivot(index="readout", columns="mode", values="aupr")
        t = t.reindex([r for r in order if r in t.index])[["raw", "hub_mi", "hub_att"]].reset_index()
        t.columns = ["readout", "AUPR raw", "AUPR hub-strat (MI row mean)", "AUPR hub-strat (pi)"]
        add(f"Pooled AUPR, truth {e} (base rate {sc[sc.truth == e].base_rate.iloc[0]:.4f})", t, {c: "{:.5f}" for c in t.columns})
        c = con[con.truth == e].copy()
        c["diff (95% jackknife)"] = [f"{a:+.4f} [{b:+.4f}, {h:+.4f}]" for a, b, h in zip(c.diff_aupr, c.diff_lo, c.diff_hi)]
        cp = c.pivot(index="readout", columns="mode", values="diff (95% jackknife)")
        cp = cp.reindex([r for r in order if r in cp.index])[["raw", "hub_mi", "hub_att"]].reset_index()
        cp.columns = ["readout", "raw minus t=9", "hub MI minus t=9", "hub pi minus t=9"]
        add(f"AUPR difference to t = 9, truth {e} (paired delete-half gene jackknife, descriptive)", cp)
    add("Literal pass rule (hub-stratified diff >= +0.005 and raw diff >= -0.005, in AUPR)", pd.read_csv(f"{d}/pass_fail.csv"), {})
open(f"{d}/tables.md", "w").write("\n".join(out))

# ---- figure: criterion curves per decile on a log t axis ----
cur = pd.read_csv(f"{d}/heldout_curves.csv")
panels = [("real", "mi", "MI row mean deciles, real operator"), ("real", "pi", "attention pi deciles, real operator"),
          ("real", "perm1", "random decile labels (permutation 1), real operator")]
if (cur.op == "rewired").any():
    panels.append(("rewired", "mi", "MI row mean deciles, rewired operator"))
fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
cols = plt.cm.Blues(np.linspace(0.3, 1.0, 10))
ncurve = cur.groupby(["split", "direction"]).ngroups
for ax, (op, ls, ttl) in zip(axs.ravel(), panels):
    x = cur[(cur.op == op) & (cur.labelset == ls)]
    for k in range(1, 11):
        z = x[x.decile == k].groupby("t").criterion.mean().reindex(GRID)
        ax.plot(GRID, z / z.max(), color=cols[k - 1], lw=1.6, label=f"decile {k}")
        ax.plot([GRID[int(np.argmax(z.values))]], [1.0], marker="o", ms=4, color=cols[k - 1])
    ax.axvline(9, color="0.5", ls=":", lw=1)
    ax.set_xscale("log", base=2); ax.set_xticks(GRID); ax.set_xticklabels([str(g) for g in GRID], fontsize=8)
    ax.set_title(ttl, fontsize=10); ax.grid(alpha=0.2)
    ax.set_ylabel("held-out AP / its maximum"); ax.set_ylim(top=1.02)
for ax in axs[1]: ax.set_xlabel("depth t (log2 axis); dotted line t = 9")
axs[0, 0].legend(fontsize=7, ncol=2, frameon=False, loc="lower right", title="1 = least connected, 10 = most", title_fontsize=7)
fig.suptitle(f"Idea 3 ({ds}): label-free criterion curves per hub decile (mean of {ncurve} split-directions; exploratory)", fontsize=11)
fig.tight_layout(); fig.savefig(f"{d}/criterion_curves_per_decile.png", dpi=150)
print("wrote", f"{d}/tables.md")
