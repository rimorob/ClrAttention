"""Section 3: graph-signal likelihood on toy data. Graph built from half 1 (labs 1-5), likelihood scored on half 2 (labs 6-10, blocked by lab)."""
from common import *
from toy_lib import *
from scipy.optimize import minimize

rng = np.random.default_rng(21)
G, NLAB, PER = 240, 10, 40
N = NLAB * PER
sizes = [55, 45, 35, 30, 20, 15]; loads = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
X, lab = simulate_modules(G, N, sizes, loads, rng)
labid = np.repeat(np.arange(NLAB), PER)
# laboratory (batch) factors: each lab shifts a random 25% of genes by a lab-specific amount
for L in range(NLAB):
    genes = rng.random(G) < 0.25
    X[np.ix_(genes, np.where(labid == L)[0])] += 0.8 * rng.standard_normal()
h1 = labid < 5; h2 = ~h1
mu1, sd1 = X[:, h1].mean(1, keepdims=True), X[:, h1].std(1, keepdims=True)
Xs = (X - mu1) / sd1
C1 = np.corrcoef(Xs[:, h1])

KEFFS = [5, 7, 10, 14, 20, 28, 40]
TS = [1, 2, 3, 4, 6, 9, 13, 18, 27, 40, 64]

def fit_cs(e, n_tr, lam_t):
    """ML fit of f_k = c * lam_t + s2 to mean energies e (n_tr experiments); returns (c, s2)."""
    def nll(p):
        f = np.exp(p[0]) * lam_t + np.exp(p[1])
        return 0.5 * np.sum(np.log(f) + e / f)
    best = None
    for p0 in ([np.log(5.0), np.log(0.5)], [np.log(1.0), np.log(0.1)]):
        r = minimize(nll, p0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 400})
        if best is None or r.fun < best.fun: best = r
    return np.exp(best.x[0]), np.exp(best.x[1])

surf = np.zeros((len(KEFFS), len(TS), 5))      # per held-out lab block: mean log-likelihood per experiment (constant dropped)
store = {}
for a, kf in enumerate(KEFFS):
    A = operator_from_corr(C1, kf, min(4 * kf, G - 2))
    lam, U, d = sym_lazy_eigs(A)
    Y = U.T @ Xs[:, h2]                 # projections of held-out experiments on the modes
    blocks = [np.where(labid[h2] == L)[0] for L in range(5, 10)]
    E = [(Y[:, b] ** 2).mean(1) for b in blocks]; Nb = [len(b) for b in blocks]
    store[kf] = (lam, E, Nb)
    for bi, tt in enumerate(TS):
        lt = np.clip(lam, 1e-9, 1) ** tt
        for f in range(5):
            tr = [j for j in range(5) if j != f]
            e_tr = sum(E[j] * Nb[j] for j in tr) / sum(Nb[j] for j in tr)
            c, s2 = fit_cs(e_tr, sum(Nb[j] for j in tr), lt)
            fk = c * lt + s2
            surf[a, bi, f] = -0.5 * np.sum(np.log(fk) + E[f] / fk)
mean = surf.mean(2)
se = surf.std(2, ddof=1) / np.sqrt(5)
a_best, b_best = np.unravel_index(np.argmax(mean), mean.shape)
print("best keff", KEFFS[a_best], "best t", TS[b_best], "loglik/exp", mean[a_best, b_best].round(2))
for a, kf in enumerate(KEFFS):
    print("keff", kf, "best t", TS[int(np.argmax(mean[a]))], "max", mean[a].max().round(1))

# ---- (a) energy profile ------------------------------------------------------------------------
kf = KEFFS[a_best]; lam, E, Nb = store[kf]
e_all = sum(E[j] * Nb[j] for j in range(5)) / sum(Nb)
fig, ax = plt.subplots(figsize=(8.6, 4.9))
ax.scatter(lam, e_all, s=26, color=BLACK, alpha=0.75, label="measured energy of each mode", zorder=3)
xs = np.linspace(max(lam.min(), 0.01), 1, 300)
for tt, c_, ls_, nm in [(1, SKY, "--", "t = 1 (too little smoothing)"), (TS[b_best], VERM, "-", f"t = {TS[b_best]} (best held-out)"), (64, ORANGE, "--", "t = 64 (too much smoothing)")]:
    c, s2 = fit_cs(e_all, sum(Nb), np.clip(lam, 1e-9, 1) ** tt)
    ax.plot(xs, c * xs ** tt + s2, color=c_, lw=3 if ls_ == "-" else 2.2, ls=ls_, label=nm)
    if ls_ == "-": s2_best = s2; c_best = c
ax.axhline(s2_best, color=VERM, lw=1.2, ls=":")
ax.text(0.27, 0.36, f"dotted line: noise floor $\\sigma^2$ = {s2_best:.2f}", color=VERM, fontsize=12)
ax.set_ylim(0.3, 2000)
ax.set_yscale("log"); ax.set_xlabel("mode eigenvalue $\\lambda_k$  (1 = smooth)")
ax.set_ylabel("energy $e_k$ (mean squared\nprojection on mode k)")
ax.set_title(f"Energy profile of held-out experiments (keff = {kf})")
ax.legend(loc="upper left", frameon=False, fontsize=11)
stamp(fig, "240 genes, 400 experiments, 10 labs")
save(fig, "fig3a_energy.png")

# ---- (b),(c) held-out likelihood ----------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axs[0]
m = mean[a_best] - mean[a_best, b_best]
d = surf[a_best] - surf[a_best, b_best][None, :]          # paired differences to the best depth, per held-out lab
s = d.std(1, ddof=1) / np.sqrt(5)
ax.errorbar(TS, m, yerr=s, fmt="-o", color=BLUE, lw=2.2, ms=7, capsize=4)
ax.axhline(0, color=GREY, lw=1, ls=":")
ax.scatter([TS[b_best]], [0], s=230, marker="*", color=VERM, edgecolor="black", zorder=5, label="maximum")
ok = np.where(m + s >= 0)[0]
one_se = ok.min()
ax.set_xscale("log"); ax.set_xticks(TS); ax.set_xticklabels([str(t) for t in TS], fontsize=11)
ax.set_xlabel("depth t (log scale)"); ax.set_ylabel("held-out log-likelihood minus the\nmaximum (bars: paired SE over 5 labs)")
ax.set_title(f"At keff = {KEFFS[a_best]}: one peak in depth")
ax.text(1.0, 3.0, "left: under-smoothed (noisy)", fontsize=12, color="#444444")
ax.text(64, 3.0, "right: over-smoothed\n(blurred)", fontsize=12, color="#444444", ha="right", va="top")
ax.set_ylim(-11, 3.5)
print("1-SE choice at best keff: t =", TS[one_se], "; paired SE at that point", s[one_se].round(2), "diff", m[one_se].round(2))
# how far are the ridge cells from the maximum, in paired SEs
for a, kf_ in enumerate(KEFFS):
    bt = int(np.argmax(mean[a])); dd = surf[a, bt] - surf[a_best, b_best]
    print("ridge keff", kf_, "t", TS[bt], "diff to max", dd.mean().round(2), "paired SE", (dd.std(ddof=1) / np.sqrt(5)).round(2))
ax = axs[1]
rel = mean - mean.max()
im = ax.imshow(np.maximum(rel, -40), origin="lower", aspect="auto", cmap="viridis")
ax.set_xticks(range(len(TS))); ax.set_xticklabels(TS, fontsize=11); ax.set_yticks(range(len(KEFFS))); ax.set_yticklabels(KEFFS)
ax.set_xlabel("depth t"); ax.set_ylabel("breadth keff")
ridge = [int(np.argmax(mean[a])) for a in range(len(KEFFS))]
ax.plot(ridge, range(len(KEFFS)), "-o", color="white", lw=2.5, ms=8, mec="black", label="best t for each keff (the ridge)")
ax.scatter([b_best], [a_best], s=330, marker="*", color=VERM, edgecolor="black", zorder=6, label="overall maximum")
cb = fig.colorbar(im, ax=ax); cb.set_label("log-likelihood relative to the maximum\n(clipped at -40)", fontsize=11)
ax.set_title("Breadth and depth trade off along a ridge")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=11)
stamp(fig, "same data; 5 held-out labs")
fig.tight_layout()
save(fig, "fig3bc_heldout.png")
np.save(os.path.join(HERE, "fig3_surface.npy"), mean)
