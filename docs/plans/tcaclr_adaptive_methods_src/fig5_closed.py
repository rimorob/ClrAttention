"""Section 5: a closed module breaks the reachable-only second CLR; the diffusion-PMI readout (toy illustration, not the planned test)."""
from common import *
from toy_lib import *

rng = np.random.default_rng(31)
nC, nO, nB = 20, 20, 80
G = nC + nO + nB
lab = np.array([0] * nC + [1] * nO + [2] * nB)       # 0 closed module, 1 open module, 2 background
A = np.zeros((G, G))
def block(i, j, p):
    u = rng.random((len(i), len(j))) < p
    A[np.ix_(i, j)] = np.maximum(A[np.ix_(i, j)], u)
iC, iO, iB = (np.where(lab == k)[0] for k in range(3))
block(iC, iC, 0.5); block(iO, iO, 0.5); block(iB, iB, 0.07)
block(iO, iB, 0.03); block(iC, iO, 0.0); block(iC, iB, 0.0)      # the closed module has no edge to the rest; the open one leaks into the background
A = np.triu(A, 1); A = A + A.T
for v in range(G):                                                # give every gene at least one partner (so the walk is defined)
    if A[v].sum() == 0:
        u = rng.choice(iB[iB != v]); A[v, u] = A[u, v] = 1
P = lazy_P(A, 0.5)
T = 9
Pt = np.linalg.matrix_power(P, T)
d = A.sum(1); pi = d / d.sum()

def a3r(Pt):
    """Port of sk_a3r: bilateral CLR on log P^t, row statistics over reachable pairs only, clipped at 0."""
    with np.errstate(divide="ignore"):
        L = np.log(Pt)
    L[~np.isfinite(L)] = np.nan; np.fill_diagonal(L, np.nan)
    mu = np.nanmean(L, 1, keepdims=True); n_ = np.sum(~np.isnan(L), 1, keepdims=True)
    Z = (L - mu) / np.sqrt(np.nansum((L - mu) ** 2, 1, keepdims=True) / (n_ - 1))
    Z[np.isnan(Z) | (Z < 0)] = 0
    S = (Z + Z.T) / np.sqrt(2); np.fill_diagonal(S, 0)
    return S, L

S_a, L = a3r(Pt)
reach = np.isfinite(L)
print("reachable entries per row: closed", reach[iC].sum(1).mean().round(1), " open", reach[iO].sum(1).mean().round(1), " background", reach[iB].sum(1).mean().round(1))

EPS = 1e-3
with np.errstate(divide="ignore"):
    W = np.log(np.maximum(Pt, EPS * pi[None, :]) / pi[None, :])
np.fill_diagonal(W, np.nan)
med = np.nanmedian(W); mad = 1.4826 * np.nanmedian(np.abs(W - med))
Zb = (W - med) / mad; Zb[np.isnan(Zb) | (Zb < 0)] = 0
S_b = (Zb + Zb.T) / np.sqrt(2); np.fill_diagonal(S_b, 0)

iu = np.triu_indices(G, 1)
def recov(S):
    sc = S[iu]; N = nC * (nC - 1) // 2 + nO * (nO - 1) // 2
    top = np.argsort(-sc, kind="stable")[:N]
    mk = np.zeros(len(sc), bool); mk[top] = True
    li, lj = lab[iu[0]], lab[iu[1]]
    out = []
    for m in (0, 1):
        inm = (li == m) & (lj == m); out.append(mk[inm].mean())
    return out
ra, rb = recov(S_a), recov(S_b)
print("share of module pairs in the global top-N: reachable-only (closed, open) =", np.round(ra, 2), " diffusion-PMI =", np.round(rb, 2))
# z of in-module pairs
def zin(S, m):
    ix = np.where(lab == m)[0]; sub = S[np.ix_(ix, ix)]; return sub[np.triu_indices(len(ix), 1)]
print("median in-module z: reachable-only closed/open", np.median(zin(S_a, 0)).round(2), np.median(zin(S_a, 1)).round(2),
      " PMI closed/open", np.median(zin(S_b, 0)).round(2), np.median(zin(S_b, 1)).round(2))

# ---- sensitivity (printed only): other standardisations of the same W --------------------------------------------------------
Wr = W.copy()
rm = np.nanmedian(Wr, 1, keepdims=True); rs = 1.4826 * np.nanmedian(np.abs(Wr - rm), 1, keepdims=True); rs[rs == 0] = np.nanstd(Wr)
Zr = (Wr - rm) / rs; Zr[np.isnan(Zr) | (Zr < 0)] = 0; Sr = (Zr + Zr.T) / np.sqrt(2); np.fill_diagonal(Sr, 0)
rr_ = recov(Sr)
fig = plt.figure(figsize=(13, 7.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], hspace=0.42, wspace=0.18)
def heat(ax, M, cmap, title, vmin=None, vmax=None, cbl=""):
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_yticks([]); ax.set_title(title, fontsize=13)
    for b in (nC, nC + nO):
        ax.axhline(b - 0.5, color="black", lw=0.8); ax.axvline(b - 0.5, color="black", lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.1, orientation="horizontal"); cb.set_label(cbl, fontsize=11)
    ax.set_xticks([nC / 2, nC + nO / 2, nC + nO + nB / 2]); ax.set_xticklabels(["C", "O", "bg"], fontsize=12)
    ax.tick_params(axis="x", length=0)
ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1]); ax3 = fig.add_subplot(gs[0, 2])
Lm = np.ma.masked_invalid(L)
cm = plt.get_cmap("viridis").copy(); cm.set_bad("#DDDDDD")
heat(ax1, Lm, cm, "log of the walk probability, t = 9\n(grey = unreachable)", cbl="log P$^t$")
heat(ax2, S_a, SEQ, "reachable-only row z (current)", vmin=0, vmax=3, cbl="score (clipped at 0)")
heat(ax3, S_b, SEQ, "stationary-null version A (proposal)", vmin=0, vmax=np.percentile(S_b, 99.5), cbl="score (clipped at 0)")
ax4 = fig.add_subplot(gs[1, :2])
x = np.arange(2); w = 0.27
ax4.bar(x - w, ra, w, color=ORANGE, edgecolor="black", label="reachable-only z (current)")
ax4.bar(x, rb, w, color=BLUE, edgecolor="black", label="stationary-null, one global scale (proposal A)")
ax4.bar(x + w, rr_, w, color=SKY, edgecolor="black", label="stationary-null, scale per row (proposal B)")
for i in range(2):
    for dx, v in zip((-w, 0, w), (ra[i], rb[i], rr_[i])):
        ax4.text(i + dx, v + 0.02, f"{v:.2f}", ha="center", fontsize=12)
ax4.set_xticks(x); ax4.set_xticklabels(["closed module (20 genes)", "open module (20 genes)"], fontsize=12.5)
ax4.set_ylabel("share of the module's gene pairs\namong the top-ranked pairs overall"); ax4.set_ylim(0, 1.45)
ax4.set_title("Do the module's pairs rank near the top? (all pairs ranked together)", fontsize=13)
ax4.legend(loc="upper right", frameon=False, fontsize=10.5, ncol=1, bbox_to_anchor=(1.0, 1.0))
ax5 = fig.add_subplot(gs[1, 2]); ax5.axis("off")
ax5.text(0, 1.0, "Reading the toy", fontsize=13, fontweight="bold", va="top")
ax5.text(0, 0.90, f"Current score: a closed gene reaches only\n{reach[iC].sum(1).mean():.0f} genes (its own module), so its z is\ncomputed against that module and the\nmodule has no background to stand out\nfrom. Its pairs rank low.\n\nProposal: the closed module is rescued.\nBut in this toy the open module slips,\nbecause pairs of low-degree genes get\nhigh scores too. Which scaling is right\nis open; the planted simulation has to\ndecide.", fontsize=11.5, va="top")
stamp(fig, "120 genes; walk depth t = 9; illustration of a proposal, not the planned test")
save(fig, "fig5_closed.png")

