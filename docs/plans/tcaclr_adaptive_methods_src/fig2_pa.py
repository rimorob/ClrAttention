"""Section 2: Horn's parallel analysis on toy expression data, and the same idea applied to the operator (spectral-edge depth)."""
from common import *
from toy_lib import *
import json

rng = np.random.default_rng(5)
G, E = 300, 120
sizes = [70, 50, 40, 30, 22, 18, 6, 5]
loads = [0.7, 0.8, 0.65, 0.7, 0.8, 0.75, 0.55, 0.5]
X, lab = simulate_modules(G, E, sizes, loads, rng)
ev, C = corr_eigs(X)
evals_full, vecs = np.linalg.eigh(C); o = np.argsort(-evals_full); vecs = vecs[:, o]
NPERM = 20
perms = [permute_genes(X, rng) for _ in range(NPERM)]
null = np.array([corr_eigs(Xp)[0] for Xp in perms])
edge = np.percentile(null[:, 0], 95)
k = int((ev > edge).sum())

# which planted modules are recovered: the module's gene-indicator pattern (centred) must have at least 50% of its energy
# inside the subspace spanned by the k retained eigenvectors (chance level is about k / G = 2%)
rec, share = {}, {}
Vk = vecs[:, :k]
for m in range(len(sizes)):
    ind = (lab == m).astype(float); ind -= ind.mean()
    share[m] = float(((Vk.T @ ind) ** 2).sum() / (ind ** 2).sum())
    rec[m] = share[m] >= 0.5
print("energy share of each module in the kept subspace:", {m: round(v, 2) for m, v in share.items()})
print("PA dimension k =", k, " edge =", round(edge, 2), " recovered:", rec)
print("null band top eigenvalue min/max:", null[:, 0].min().round(2), null[:, 0].max().round(2))

# ---- figure 2a ----------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]
nshow = 16
x = np.arange(1, nshow + 1)
ax.fill_between(x, null[:, :nshow].min(0), null[:, :nshow].max(0), color=GREY, alpha=0.45, label="permuted data (20 runs): min to max")
ax.axhline(edge, color=VERM, lw=2.2, ls="--", label=f"edge = {edge:.1f}\n(95th percentile of the\nlargest permuted eigenvalue)")
above = ev[:nshow] > edge
ax.plot(x, ev[:nshow], "-", color=BLACK, lw=1.2, zorder=3)
ax.scatter(x[above], ev[:nshow][above], s=90, color=BLUE, edgecolor="black", zorder=4, label=f"above the edge: {k} modes kept")
ax.scatter(x[~above], ev[:nshow][~above], s=90, color="white", edgecolor="black", zorder=4, label="below the edge: discarded")
ax.set_xlabel("eigenvalue number"); ax.set_ylabel("eigenvalue of the gene-gene\ncorrelation matrix")
ax.set_title(f"Scree plot: {k} eigenvalues beat chance")
ax.legend(loc="upper right", frameon=False, fontsize=11)

ax = axs[1]
order = np.arange(len(sizes))
strength = [s * b * b for s, b in zip(sizes, loads)]
cols = [BLUE if rec[m] else "white" for m in order]
ax.barh(order, strength, color=cols, edgecolor="black", height=0.65)
for m in order:
    ax.text(strength[m] + 0.6, m, f"{sizes[m]} genes, load {loads[m]}  -  " + ("recovered" if rec[m] else "MISSED"),
            va="center", fontsize=11.5, color="black" if rec[m] else VERM, fontweight="bold" if not rec[m] else "normal")
ax.set_yticks(order); ax.set_yticklabels([f"module {m + 1}" for m in order]); ax.invert_yaxis()
ax.set_xlabel("planted signal strength (genes $\\times$ load$^2$)")
ax.set_xlim(0, 70)
ax.set_title("Which planted modules were found")
stamp(fig, "300 genes, 120 experiments, 8 planted modules")
fig.tight_layout()
save(fig, "fig2a_pa.png")

# ---- operator-level PA (spectral edge) ----------------------------------------------------------
KEFF, CAP = 10, 40
A_obs = operator_from_corr(C, KEFF, CAP)
lam_obs, U_obs, _ = sym_lazy_eigs(A_obs)
lam_null = []
for Xp in perms:
    Cp = np.corrcoef(Xp)
    lp, _, _ = sym_lazy_eigs(operator_from_corr(Cp, KEFF, CAP))
    lam_null.append(lp)
lam_null = np.array(lam_null)
# edge = median over permutations of the largest NON-trivial eigenvalue (the first one, lambda = 1, is trivial)
lam_edge = float(np.median(lam_null[:, 1]))
eps = 0.05
tstar = np.log(eps) / np.log(lam_edge)
nk = int((lam_obs[1:] > lam_edge).sum())
print("operator: lam_edge =", round(lam_edge, 3), " null band lam2 min/max:", lam_null[:, 1].min().round(3), lam_null[:, 1].max().round(3),
      " t* =", round(tstar, 1), " modes above edge (nontrivial):", nk)
for e in (0.1, 0.05, 0.01):
    print("eps", e, "t* =", round(np.log(e) / np.log(lam_edge), 1))
for kk in (6, 10, 20):
    le = float(np.median(np.array([sym_lazy_eigs(operator_from_corr(np.corrcoef(Xp), kk, 4 * kk))[0][1] for Xp in perms[:5]])))
    print("keff", kk, "lam_edge (5 perms)", round(le, 3), "t* (eps .05) =", round(np.log(.05) / np.log(le), 1))
json.dump({"lam_edge": lam_edge, "tstar": tstar, "n_modes_edge": nk, "k_corr": k, "edge_corr": float(edge)},
          open(os.path.join(HERE, "fig2_numbers.json"), "w"), indent=1)

fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axs[0]
nshow = 30
x = np.arange(1, nshow + 1)
ax.fill_between(x, lam_null[:, :nshow].min(0), lam_null[:, :nshow].max(0), color=GREY, alpha=0.45, label="operator rebuilt on permuted data\n(20 runs): min to max")
ax.axhline(lam_edge, color=VERM, lw=2.2, ls="--", label=f"edge $\\lambda_{{edge}}$ = {lam_edge:.2f}")
ab = lam_obs[:nshow] > lam_edge
ax.plot(x, lam_obs[:nshow], "-", color=BLACK, lw=1.2)
ax.scatter(x[ab], lam_obs[:nshow][ab], s=70, color=BLUE, edgecolor="black", zorder=4, label=f"observed, above edge ({nk} besides the first)")
ax.scatter(x[~ab], lam_obs[:nshow][~ab], s=70, color="white", edgecolor="black", zorder=4, label="observed, below edge")
ax.set_xlabel("mode number"); ax.set_ylabel("eigenvalue of the walk operator")
ax.set_title("PA on the operator: observed vs permuted")
ax.legend(loc="upper right", frameon=False, fontsize=10.5)

ax = axs[1]
lg = np.linspace(0.3, 1, 400)
for t, c, ls in [(1, SKY, "-"), (3, BLUE, "-"), (9, ORANGE, "-")]:
    ax.plot(lg, lg ** t, color=c, lw=2.2, ls=ls, label=f"t = {t}")
ax.plot(lg, lg ** tstar, color=VERM, lw=3.2, label=f"t* = {tstar:.1f}  (chosen)")
ax.axvline(lam_edge, color=VERM, lw=1.6, ls="--")
ax.axhline(eps, color=GREY, lw=1.4, ls=":")
ax.scatter([lam_edge], [eps], s=140, color=VERM, edgecolor="black", zorder=6)
ax.text(0.32, eps + 0.02, f"$\\epsilon$ = {eps}", color="#444444", fontsize=12)
ax.annotate("at the edge, gain = $\\epsilon$", xy=(lam_edge, eps), xytext=(lam_edge - 0.40, 0.30), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="black"))
ax.axvspan(lam_edge, 1, color=BLUE, alpha=0.10, label="shaded: modes above the edge (kept)")
ax.text(0.5 * (0.3 + lam_edge), 0.98, "modes below the edge:\nfaded to gain < 0.05 at t*", ha="center", va="top", fontsize=12, color="#555555")
ax.set_xlabel("mode eigenvalue $\\lambda$"); ax.set_ylabel("gain after t steps  =  $\\lambda^t$")
ax.set_xlim(0.3, 1); ax.set_ylim(0, 1.0)
ax.set_title("Turning the edge into a depth")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=11)
stamp(fig, "same toy data as the previous figure; keff 10")
fig.tight_layout()
save(fig, "fig2b_edge.png")
