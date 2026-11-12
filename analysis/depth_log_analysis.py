#!/usr/bin/env python3
"""Pre-specified analysis of the log-grid depth re-test (docs/preregistration_depth_log.md, sections 5-10).
Reads results/granularity_log/{curves,boot,truth_sizes}_{ecoli,k562}.*; writes tables, two figures and the numbers
section of docs/depth_log_results.md. Usage (repo root): python3 analysis/depth_log_analysis.py [--dir results/granularity_log] [--doc docs/depth_log_results.md]"""
import sys, itertools, math, argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import rankdata

ap = argparse.ArgumentParser(); ap.add_argument("--dir", default="results/granularity_log"); ap.add_argument("--doc", default="docs/depth_log_results.md")
ap.add_argument("--no_doc", action="store_true"); a = ap.parse_args(); D = a.dir + "/"
EDGE = 256
ORGS = [("ecoli", "E. coli M3D"), ("k562", "K562 essential compendium")]

def rho_of(x, y):
    rx, ry = rankdata(x), rankdata(y)
    if np.ptp(rx) == 0 or np.ptp(ry) == 0: return 0.0     # undefined -> 0 (prereg section 6)
    return float(np.corrcoef(rx, ry)[0, 1])

def exact_perm(x, y):
    obs = rho_of(x, y); r = np.array([rho_of(x, np.array(p)) for p in itertools.permutations(y)])
    return obs, float(np.mean(r >= obs - 1e-12)), float(np.mean(np.abs(r) >= abs(obs) - 1e-12)), len(r)

def interval(t, y):
    ib = int(np.argmax(y)); ok = y >= 0.99 * y[ib]; lo = hi = ib
    while lo > 0 and ok[lo - 1]: lo -= 1
    while hi < len(t) - 1 and ok[hi + 1]: hi += 1
    return t[lo], t[hi]

def slope(x, y):
    if len(x) < 3 or np.ptp(x) == 0: return np.nan
    return float(np.polyfit(x, y, 1)[0])

curves, boots, sizes = {}, {}, {}
for key, _ in ORGS:
    curves[key] = pd.read_csv(D + f"curves_{key}.csv"); boots[key] = pd.read_csv(D + f"boot_{key}.csv.gz"); sizes[key] = pd.read_csv(D + f"truth_sizes_{key}.csv")

rows = []
for key, lab in ORGS:
    c = curves[key]
    for truth in sorted(c.truth.unique()):
        g = c[c.truth == truth].sort_values("t"); t = g.t.values; y = g.aupr.values
        base = y[t == 0][0]; tb = int(t[int(np.argmax(y))])       # primary: argmax over {0} + grid (first max)
        m = t >= 1; tb1 = int(t[m][int(np.argmax(y[m]))]); lo, hi = interval(t[m], y[m])
        s = sizes[key].set_index("truth").loc[truth]
        rows.append(dict(org=key, truth=truth, rank=int(truth.split()[0]), t_best=tb, t_best_ge1=tb1, aupr_clr=base, aupr_best=y.max(),
                         gain_best=y.max() - base, t_lo=int(lo), t_hi=int(hi), width_octaves=math.log2(hi / lo), interval_touches_edge=bool(hi == EDGE),
                         censored=bool(tb == EDGE), clr_best=bool(tb == 0), median_size=float(s.median_size), n_modules=int(s.n_modules), n_pairs=int(s.n_pairs)))
R = pd.DataFrame(rows)
bw = []
for key, _ in ORGS:
    b = boots[key]
    for truth, g in b.groupby("truth"):
        tb = g.t_best.values; pos = tb[tb >= 1]
        bw.append(dict(org=key, truth=truth, t_best_p2_5=np.percentile(tb, 2.5), t_best_median=np.median(tb), t_best_p97_5=np.percentile(tb, 97.5), frac_clr_best=float(np.mean(tb == 0)),
                       frac_at_edge=float(np.mean(tb == EDGE)), med_size_p2_5=np.percentile(g.median_size, 2.5), med_size_p97_5=np.percentile(g.median_size, 97.5)))
BW = pd.DataFrame(bw); R = R.merge(BW, on=["org", "truth"])

# ---- primary (section 6) and reading (section 10)
prim = []
for key, lab in ORGS:
    s = R[R.org == key].sort_values("rank"); x = s["rank"].values; y = s.t_best.values; y1 = s.t_best_ge1.values
    rho, p1, p2, nperm = exact_perm(x, y); rho1, p1_1, p2_1, _ = exact_perm(x, y1)
    prim.append(dict(org=key, n=len(s), rho=rho, p_one_sided=p1, p_two_sided=p2, n_perm=nperm, rho_t_best_ge1=rho1, p_one_sided_ge1=p1_1,
                     n_censored=int(s.censored.sum()), n_interval_touches_edge=int(s.interval_touches_edge.sum()), n_clr_best=int(s.clr_best.sum()),
                     min_attainable_p=1 / nperm))
P = pd.DataFrame(prim)
rho_e, rho_k = P.set_index("org").rho["ecoli"], P.set_index("org").rho["k562"]; any_cens = bool(R.censored.any())
if rho_e >= 0.9 and rho_k >= 0.9 and not any_cens: reading = "(a) resolution-like on a log axis"
elif rho_e <= 0.5 or rho_k <= 0.5: reading = "(b) not resolution-like"
else: reading = "(c) inconclusive"

# ---- secondary: slope vs size, bootstrap (section 7)
def fit_slopes(df):
    out = {}
    for name, sel in (("pooled", slice(None)), ("ecoli", "ecoli"), ("k562", "k562")):
        d = df if name == "pooled" else df[df.org == name]; d = d[d.t_best >= 1]
        need = 3
        out[name] = slope(np.log2(d.median_size.values), np.log2(d.t_best.values)) if len(d) >= need and (name != "ecoli" or len(d) == 3) else np.nan
    return out
obs = fit_slopes(R)
B = int(boots["ecoli"].draw.max()); bs = {"pooled": [], "ecoli": [], "k562": []}; rhoB = {"ecoli": [], "k562": []}
big = pd.concat([boots[k].assign(org=k, rank=boots[k].truth.str.split().str[0].astype(int)) for k, _ in ORGS])
for d, g in big.groupby("draw"):
    f = fit_slopes(g.rename(columns={"t_best": "t_best"}))
    for k in bs: bs[k].append(f[k])
    for k in rhoB:
        s = g[g.org == k].sort_values("rank"); rhoB[k].append(rho_of(s["rank"].values, s.t_best.values))
sl = []
for k in ("pooled", "ecoli", "k562"):
    v = np.array(bs[k]); ok = ~np.isnan(v)
    sl.append(dict(set=k, slope=obs[k], ci_lo=np.percentile(v[ok], 2.5) if ok.sum() > 10 else np.nan, ci_hi=np.percentile(v[ok], 97.5) if ok.sum() > 10 else np.nan,
                   n_draws_used=int(ok.sum()), n_draws_dropped=int((~ok).sum()), n_truths_in_fit=int(((R.t_best >= 1) & ((R.org == k) | (k == "pooled"))).sum())))
SL = pd.DataFrame(sl)
rb = pd.DataFrame([dict(org=k, frac_rho_ge_0_9=float(np.mean(np.array(v) >= 0.9)), frac_rho_le_0_5=float(np.mean(np.array(v) <= 0.5)),
                        median_rho=float(np.median(v))) for k, v in rhoB.items()])
cens = pd.DataFrame([dict(scope=k, n=int((R.org == k).sum() if k != "total" else len(R)), n_censored=int(R[(R.org == k) | (k == "total")].censored.sum()),
                          n_interval_touches_edge=int(R[(R.org == k) | (k == "total")].interval_touches_edge.sum())) for k in ("ecoli", "k562", "total")])
R.to_csv(D + "depth_log_per_truth.csv", index=False); P.to_csv(D + "depth_log_primary.csv", index=False); SL.to_csv(D + "depth_log_slope.csv", index=False)
rb.to_csv(D + "depth_log_boot_rho.csv", index=False); cens.to_csv(D + "depth_log_censoring.csv", index=False)

# ---- figures
INK, MUTE, BLUE, RED, SURF = "#0b0b0b", "#52514e", "#2a78d6", "#e34948", "#fcfcfb"
grid = sorted(t for t in curves["ecoli"].t.unique() if t >= 1); X0 = -1.3
xt = [X0] + [math.log2(t) for t in (1, 2, 4, 8, 16, 32, 64, 128, 256)]; xl = ["CLR"] + [str(t) for t in (1, 2, 4, 8, 16, 32, 64, 128, 256)]
fig, axs = plt.subplots(2, 4, figsize=(15, 7.2)); axs = axs.ravel(); n = 0
for key, lab in ORGS:
    for truth in sorted(curves[key].truth.unique(), key=lambda s: int(s.split()[0])):
        ax = axs[n]; n += 1; g = curves[key][curves[key].truth == truth].sort_values("t"); r = R[(R.org == key) & (R.truth == truth)].iloc[0]
        base = g.aupr.iloc[0]; gg = g[g.t >= 1]; lx = np.log2(gg.t.values); gain = gg.aupr.values - base
        ax.axvspan(math.log2(r.t_lo) - 0.08, math.log2(r.t_hi) + 0.08, color=BLUE, alpha=0.15, lw=0)
        ax.axhline(0, color=MUTE, lw=0.8, ls=":"); ax.plot(lx, gain, color=BLUE, lw=2)
        ax.plot([X0], [0], "o", mfc=SURF, mec=MUTE, ms=7, mew=1.5)
        if r.t_best == 0: ax.plot([X0], [0], "o", color=BLUE, ms=9, mec=SURF, mew=2); txt, xx, yy = "CLR best", X0, 0.0
        else: yy = float(g[g.t == r.t_best].aupr.iloc[0] - base); xx = math.log2(r.t_best); ax.plot([xx], [yy], "o", color=BLUE, ms=8, mec=SURF, mew=2); txt = f"t = {r.t_best}" + (" (edge)" if r.censored else "")
        ax.annotate(txt, (xx, yy), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8.5, color=INK)
        ax.axvline(math.log2(EDGE), color=RED, lw=1.0, ls="--", alpha=0.7)
        ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=7.5, color=MUTE); ax.set_xlim(-1.7, 8.4); ax.margins(y=0.18)
        ax.set_title(f"{'E. coli' if key == 'ecoli' else 'K562'}  {truth}", fontsize=10, color=INK, loc="left", pad=12)
        ax.set_xlabel("walk depth t (log2 axis)", fontsize=9, color=MUTE); ax.tick_params(axis="y", colors=MUTE, labelsize=8.5)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        ax.grid(axis="y", color="#e5e4e0", lw=0.6)
axs[0].set_ylabel("AUPR minus CLR (t = 0)", fontsize=9, color=MUTE); axs[4].set_ylabel("AUPR minus CLR (t = 0)", fontsize=9, color=MUTE)
fig.suptitle("Geometric depth grid t = 1..256 (+ CLR). AUPR gain over CLR. Dot = best t (incl. CLR); shaded = within 1% of max AUPR (t >= 1); dashed red = grid edge 256.  Truths from fine to coarse.\n"
             f"Pre-specified reading: {reading};  E. coli rho = {rho_e:.2f} (n = {int(P.n[0])}),  K562 rho = {rho_k:.2f} (n = {int(P.n[1])});  edge-censored truths: {int(R.censored.sum())} of {len(R)}",
             fontsize=10, color=INK, x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(D + "depth_log_curves.png", dpi=170, facecolor=SURF); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
cols = {"ecoli": BLUE, "k562": "#eb6834"}
for key, lab in ORGS:
    s = R[R.org == key].sort_values("rank")
    yy = np.where(s.t_best.values >= 1, np.log2(np.maximum(s.t_best.values, 1)), -1.3)
    axs[0].plot(s["rank"], yy, "-o", color=cols[key], lw=2, ms=8, mec=SURF, mew=2, label=f"{lab.split()[0]} {lab.split()[1] if key=='ecoli' else ''}".strip())
axs[0].set_yticks([-1.3, 0, 2, 4, 6, 8]); axs[0].set_yticklabels(["CLR", "1", "4", "16", "64", "256"]); axs[0].axhline(8, color=RED, ls="--", lw=1, alpha=.7)
axs[0].set_xlabel("truth coarseness rank (fine to coarse)", color=MUTE); axs[0].set_ylabel("t_best (log2 axis)", color=MUTE); axs[0].legend(frameon=False, fontsize=9)
axs[0].set_title("Best depth by coarseness", loc="left", fontsize=10)
for key, lab in ORGS:
    s = R[R.org == key]; m = s.t_best >= 1
    x = np.log2(s.median_size); yv = np.log2(np.where(m, s.t_best, np.nan))
    lo = np.log2(np.maximum(s.t_best_p2_5, 1)); hi = np.log2(s.t_best_p97_5)
    axs[1].errorbar(x[m], yv[m], yerr=[(yv - lo)[m], (hi - yv)[m]], fmt="o", color=cols[key], ms=8, mec=SURF, mew=2, lw=1.2, capsize=0, label=lab.split(" M3D")[0].split(" essential")[0])
    for xi, yi, nm in zip(x[m], yv[m], s.truth[m]): axs[1].annotate(nm.split(" ", 1)[1], (xi, yi), xytext=(4, 4), textcoords="offset points", fontsize=7, color=MUTE)
sp = SL.set_index("set").loc["pooled"]
axs[1].set_xlabel("median module size (log2 genes)", color=MUTE); axs[1].set_ylabel("t_best (log2 axis, bars: bootstrap 2.5-97.5%)", color=MUTE, fontsize=8.5)
axs[1].set_title(f"Pooled slope {sp.slope:.2f} [{sp.ci_lo:.2f}, {sp.ci_hi:.2f}]", loc="left", fontsize=10); axs[1].legend(frameon=False, fontsize=9)
for key, lab in ORGS:
    s = R[R.org == key].sort_values("rank"); axs[2].plot(s["rank"], s.width_octaves, "-o", color=cols[key], lw=2, ms=8, mec=SURF, mew=2)
axs[2].set_xlabel("truth coarseness rank", color=MUTE); axs[2].set_ylabel("width of within-1% interval (octaves)", color=MUTE); axs[2].set_title("Interval width (grid-limited)", loc="left", fontsize=10)
for ax in axs:
    for spn in ("top", "right"): ax.spines[spn].set_visible(False)
    ax.grid(color="#e5e4e0", lw=0.6); ax.tick_params(colors=MUTE)
fig.tight_layout(); fig.savefig(D + "depth_log_summary.png", dpi=170, facecolor=SURF); plt.close(fig)

# ---- numbers section of the results doc
def md(df, fmt):
    cols = list(df.columns); L = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows(): L.append("| " + " | ".join(fmt.get(c, lambda v: f"{v}")(r[c]) for c in cols) + " |")
    return "\n".join(L)
f3 = lambda v: f"{v:.3f}"; f2 = lambda v: f"{v:.2f}"; fi = lambda v: f"{int(v)}"
t1 = R[["org", "truth", "t_best", "t_best_ge1", "gain_best", "t_lo", "t_hi", "width_octaves", "censored", "interval_touches_edge", "median_size", "n_modules", "t_best_p2_5", "t_best_p97_5"]].copy()
body = [f"## Numbers (generated by `analysis/depth_log_analysis.py`, do not edit by hand)\n",
        "### Per truth (t_best over {CLR} + grid; t_best_ge1 over t >= 1; 1%-interval over t >= 1; bootstrap interval of t_best, 1,000 module draws)\n",
        md(t1, {"gain_best": f3, "width_octaves": f2, "median_size": lambda v: f"{v:g}", "t_best_p2_5": lambda v: f"{v:g}", "t_best_p97_5": lambda v: f"{v:g}"}),
        "\n### Primary: Spearman rho of coarseness rank with log2 t_best (exact permutation)\n",
        md(P, {"rho": f3, "p_one_sided": f3, "p_two_sided": f3, "rho_t_best_ge1": f3, "p_one_sided_ge1": f3, "min_attainable_p": f3}),
        f"\n**Pre-specified reading (a: rho >= 0.9 in both and no edge censoring; b: rho <= 0.5 in either; else c): {reading}.**  E. coli rho = {rho_e:.3f}, K562 rho = {rho_k:.3f}, edge-censored truths = {int(R.censored.sum())} of {len(R)}.\n",
        "### Secondary: slope of log2 t_best on log2 median module size (truths with t_best >= 1), bootstrap over modules\n", md(SL, {"slope": f3, "ci_lo": f3, "ci_hi": f3}),
        "\n### Secondary: censoring at t = 256\n", md(cens, {}),
        "\n### Secondary: bootstrap distribution of the per-organism rho\n", md(rb, {"frac_rho_ge_0_9": f3, "frac_rho_le_0_5": f3, "median_rho": f3}),
        "\n### Tertiary: interval widths are in the per-truth table (width_octaves; grid-limited to about 0.4-0.6 octave steps).\n"]
open(D + "depth_log_numbers.md", "w").write("\n".join(body))
print("\n".join(body)); print("READING:", reading)
if not a.no_doc and os.path.exists(a.doc):
    s = open(a.doc).read(); m0 = "<!-- NUMBERS-BEGIN -->"; m1 = "<!-- NUMBERS-END -->"
    if m0 in s and m1 in s:
        s = s[:s.index(m0) + len(m0)] + "\n" + "\n".join(body) + "\n" + s[s.index(m1):]; open(a.doc, "w").write(s)
