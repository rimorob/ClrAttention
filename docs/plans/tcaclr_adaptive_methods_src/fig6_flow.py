"""Section 6: how the pieces fit (flow diagram). Not data; a schematic of the plan. Drawn at its print size (about 7.2 in wide)."""
from common import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig = plt.figure(figsize=(7.2, 4.8))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(-4, 70); ax.axis("off")

def box(x, y, w, h, title, body, fc="#FFFFFF", ec=BLACK, fs=7.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0", fc=fc, ec=ec, lw=1.2))
    ax.text(x + w / 2, y + h - 1.0, title, ha="center", va="top", fontsize=8.4, fontweight="bold")
    ax.text(x + w / 2, y + h - (6.8 if "\n" in title else 4.4), body, ha="center", va="top", fontsize=fs, color="#222222", linespacing=1.35)

def arrow(x1, y1, x2, y2, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11, lw=1.2, color=BLACK, linestyle=ls))

box(10, 60, 80, 9.5, "Stage 0: simulation and negative controls",
    "planted modules (sizes, strengths, N, G varied) and per-gene permuted\ndata. Run BEFORE any real-data setting is trusted.", fc="#E8F1F8", ec=BLUE)
W = 24
xs = [0.5, 25.5, 50.5, 75.5]
Y2, H2 = 33, 22
box(xs[0], Y2, W, H2, "Spectral-edge\ndepth (PA)", "permute genes, rebuild\nthe operator, take the\nedge eigenvalue:\nt* = ln(eps) / ln(edge)\nsets t, given keff")
box(xs[1], Y2, W, H2, "Graph-signal\nlikelihood", "held-out, blocked by\nlaboratory; energy\nprofile of experiments;\nsets keff and t (a\nridge), 1-SE rule")
box(xs[2], Y2, W, H2, "Instability curve\n(StARS-style)", "stability of each gene's\ntop partners across\nsubsamples; a second\nreading for keff only")
box(xs[3], Y2, W, H2, "Island\ndiagnostics", "conductance, community\nprofile, leak window;\nsets a t window, checks\nkeff, does not set keff\nalone")
for x in xs: arrow(x + W / 2, 60, x + W / 2, 55.6)
box(0.5, 0, 31, 23.5, "Revised grid and\ncorrected criteria", "17 depths x 8 keff; cap =\n4 x keff; clip before\nStouffer (2 levels);\nC1 diagnostic only, C2\nnull-adjusted, C3 held-out\nlikelihood, no penalty;\nsmoothed surface + 1-SE", fs=7.2)
box(34.5, 0, 31, 23.5, "Stage 2: spent truths", "descriptive curves only:\nregret; criterion validity\n(Spearman across cells);\nAUPR by complex size.\nCannot override stage 1.", fs=7.2)
box(68.5, 0, 31, 23.5, "Pre-written decision rule", "success needs all four:\nregret <= 5% in >= 80% of\nsimulated settings; scale\ntracks module size\n(Spearman >= 0.7); permuted\ncontrols flagged >= 95%;\nsplit stability 7 of 8.", fc="#FBEFD9", ec=ORANGE, fs=7.2)
arrow(xs[0] + 8, Y2, 8, 24.8); arrow(xs[1] + 4, Y2, 20, 24.8); arrow(xs[2] + 3, Y2, 28, 24.8)
arrow(xs[3] + 6, Y2, 52, 24.8, ls="--")
arrow(31.8, 11.5, 34.2, 11.5); arrow(65.8, 11.5, 68.2, 11.5)
ax.text(99.5, 28.6, "dashed arrow:\ndiagnostic only,\nnot a selector", ha="right", va="center", fontsize=7, color="#555555", style="italic")
ax.text(38, -2.2, "Status of every box: proposed, not yet run on any data.", fontsize=8.6, color=VERM, fontweight="bold", ha="center", va="center")
fig.text(0.995, -0.005, "schematic, no data", ha="right", va="top", fontsize=7, color="#555555", style="italic")
save(fig, "fig6_flow.png")
