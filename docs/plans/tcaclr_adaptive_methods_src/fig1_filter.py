"""Section 1: the walk as a low-pass filter on a toy 60-node network (3 modules + noise edges)."""
from common import *

rng = np.random.default_rng(11)
sizes = [20, 20, 20]
A, lab = sbm_adj(sizes, [0.25, 0.25, 0.25], 0.03, rng)
n = len(A)
lam, Phi, d = sym_norm_eig(A)
pos = circle_layout(lab, np.random.default_rng(3), spread=0.9)
print("lambda (first 8):", np.round(lam[:8], 3), " last:", np.round(lam[-3:], 3))

# ---- (a) three eigenvectors ---------------------------------------------------------------------
# mode numbers are counted from 1 (mode 1 = constant, lambda = 1). Choose: smooth = mode 2, medium, high frequency
picks = [(1, "Smooth mode"), (9, "Medium mode"), (n - 3, "High-frequency mode")]
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.4))
for ax, (k, ttl) in zip(axs, picks):
    v = Phi[:, k] / np.sqrt(d)           # right eigenvector of the walk = a pattern of values over genes
    v = v / np.abs(v).max()
    draw_edges(ax, A, pos)
    ax.scatter(pos[:, 0], pos[:, 1], c=v, cmap=DIV, vmin=-1, vmax=1, s=120, edgecolor="#333333", linewidth=0.6, zorder=3)
    ax.set_title(f"{ttl} (number {k + 1})\neigenvalue = {lam[k]:.2f}", fontsize=13.5)
    ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.1); ax.autoscale_view()
sm = plt.cm.ScalarMappable(cmap=DIV, norm=plt.Normalize(-1, 1))
cb = fig.colorbar(sm, ax=axs, orientation="horizontal", fraction=0.05, pad=0.02, aspect=40)
cb.set_label("value of the pattern at each gene (blue = negative, orange = positive)", fontsize=12)
stamp(fig, "60 genes, 3 modules of 20")
save(fig, "fig1a_modes.png")

# ---- (b) gain curves ----------------------------------------------------------------------------
ks = np.arange(1, n + 1)
fig, ax = plt.subplots(figsize=(8.6, 4.6))
cols = [SKY, BLUE, ORANGE, VERM]
for t, c in zip([1, 3, 9, 27], cols):
    g = lam ** t
    ax.plot(ks, g, "-o", ms=3.5, lw=2.2, color=c, label=f"t = {t}")
    cut = np.argmax(g < 0.25)           # first mode whose gain is below 1/4
    ax.plot([cut + 1], [g[cut]], marker="v", ms=11, color=c, mec="black", zorder=5)
ax.axhline(0.25, color=GREY, lw=1, ls="--")
ax.text(n + 0.5, 0.27, "gain 1/4", ha="right", va="bottom", color="#555555", fontsize=11)
ax.set_xlabel("mode number (1 = smoothest)"); ax.set_ylabel("gain  =  (eigenvalue)$^t$")
ax.set_xlim(0, n + 1); ax.set_ylim(0, 1.04)
ax.legend(title="steps taken", loc="upper right", bbox_to_anchor=(1.0, 0.9), frameon=False)
ax.set_title("Each extra step shrinks the high modes more: a low-pass filter")
ax.annotate("triangles mark the 'cutoff mode':\nthe first mode with gain below 1/4", xy=(0.45, 0.33), xycoords="axes fraction", fontsize=11.5, color="#333333")
stamp(fig)
save(fig, "fig1b_gain.png")
print("cutoff modes:", {t: int(np.argmax(lam ** t < 0.25)) + 1 for t in [1, 3, 9, 27]})

# ---- (c) attention mass from one gene ------------------------------------------------------------
P = lazy_P(A)
src = 5   # a gene in module 0
fig, axs = plt.subplots(1, 4, figsize=(13, 3.6))
Pt = np.eye(n)
done = 0
for ax, t in zip(axs, [1, 3, 9, 27]):
    Pt = np.linalg.matrix_power(P, t)
    m = Pt[src]
    inside = m[lab == lab[src]].sum()
    draw_edges(ax, A, pos)
    sc = ax.scatter(pos[:, 0], pos[:, 1], c=m / m.max(), cmap=SEQ, vmin=0, vmax=1, s=95, edgecolor="#555555", linewidth=0.4, zorder=3)
    ax.scatter(*pos[src], s=260, facecolor="none", edgecolor=BLACK, linewidth=2.2, zorder=4)
    ax.set_title(f"t = {t}\n{inside * 100:.0f}% of mass in own module", fontsize=12.5)
    ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.1); ax.autoscale_view()
fig.text(0.5, 0.06, "start gene circled; darker = more of the walk's mass (each panel scaled to its own maximum)", ha="center", fontsize=12)
stamp(fig)
save(fig, "fig1c_spread.png")
