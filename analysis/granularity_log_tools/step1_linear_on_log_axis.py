#!/usr/bin/env python3
"""EXPLORATORY step 1: re-plot the saved LINEAR-grid curves (results/granularity/) on a log2 t axis.
Curves were already seen; nothing here is blind. No new compute.
Usage (repo root): python3 analysis/granularity_log_tools/step1_linear_on_log_axis.py"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
d = pd.read_csv("results/granularity/granularity_curves.csv")
OUT = "results/granularity_log/"
TMAX = int(d.t.max())
rows = []
fig, axs = plt.subplots(2, 4, figsize=(14, 7))
axs = axs.ravel()
BLUE, INK, MUTE = "#2a78d6", "#0b0b0b", "#52514e"
for ax, ((data, truth), g) in zip(axs, d.groupby(["data", "truth"], sort=False)):
    g = g.sort_values("t"); base = g[g.t == 0].aupr.iloc[0]
    a = g[g.t >= 1]; t = a.t.values; y = a.aupr.values
    ib = int(np.argmax(y)); tb = t[ib]; mx = y[ib]; ok = y >= 0.99 * mx
    lo = ib
    while lo > 0 and ok[lo - 1]: lo -= 1
    hi = ib
    while hi < len(t) - 1 and ok[hi + 1]: hi += 1
    clr_best = base >= mx
    rows.append(dict(data=data, truth=truth, t_best_ge1=int(tb), log2_t_best=np.log2(tb), t_lo=int(t[lo]), t_hi=int(t[hi]),
                     width_octaves=np.log2(t[hi] / t[lo]), edge_censored=bool(tb == TMAX), interval_touches_edge=bool(t[hi] == TMAX),
                     clr_beats_all_t=bool(clr_best), gain_at_best=mx - base))
    gain = y - base
    ax.axvspan(np.log2(t[lo]), np.log2(t[hi]), color=BLUE, alpha=0.15, lw=0)
    ax.axhline(0, color=MUTE, lw=0.8, ls=":")
    ax.plot(np.log2(t), gain, color=BLUE, lw=2)
    ax.plot([np.log2(tb)], [gain[ib]], "o", color=BLUE, ms=8, mec="#fcfcfb", mew=2)
    ax.annotate(f"t = {tb}" + (" (edge)" if tb == TMAX else ""), (np.log2(tb), gain[ib]), xytext=(0, 7), textcoords="offset points",
                ha="center", fontsize=8.5, color=INK)
    if tb == TMAX: ax.axvline(np.log2(TMAX), color="#e34948", lw=1.2, ls="--")
    ax.margins(y=0.18); ax.tick_params(labelbottom=True)
    ax.set_title(f"{'E. coli' if data.startswith('E.') else 'K562'}  {truth}", fontsize=10, color=INK, loc="left", pad=14)
    ax.set_xlim(-0.3, 5.3); ax.set_xticks(np.arange(0, 6)); ax.set_xticklabels([str(2 ** k) for k in range(6)], fontsize=8.5)
    ax.tick_params(colors=MUTE, labelsize=8.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#e5e4e0", lw=0.6)
for ax in axs[:8]: ax.set_xlabel("walk depth t (log2 axis)", fontsize=9, color=MUTE)
for ax in (axs[0], axs[4]): ax.set_ylabel("AUPR minus CLR (t = 0)", fontsize=9, color=MUTE)
fig.suptitle("EXPLORATORY (saved linear grid t = 1..32, already seen): AUPR gain over CLR against log2 t.  Dot = best t; shaded = within 1% of max AUPR; dashed red = grid edge",
             fontsize=10, color=INK, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(OUT + "step1_linear_grid_on_log_axis.png", dpi=170, facecolor="#fcfcfb")
r = pd.DataFrame(rows); r.to_csv(OUT + "step1_best_t_linear_on_log.csv", index=False)
pd.set_option("display.width", 200); print(r.to_string(index=False))
for org in ("E. coli", "K562"):
    s = r[r.data.str.startswith(org)]; n = len(s)
    rho, p = spearmanr(np.arange(1, n + 1), s.log2_t_best)
    print(org, "Spearman(rank, log2 t_best>=1) =", round(rho, 3), "n =", n)
    rho0, _ = spearmanr(np.arange(1, n + 1), np.where(s.clr_beats_all_t, -1, s.log2_t_best))
    print("   with CLR-best coded below t=1:", round(rho0, 3))
print("edge-censored at t=32:", int(r.edge_censored.sum()), "of", len(r))
