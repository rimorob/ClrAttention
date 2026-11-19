"""Section 4: islands and conductance on a toy graph (islands of 5, 20, 80 genes + a hub, plus noise edges)."""
from common import *
from island_lib import *
import json, collections

def build(seed=4, pout=0.002, pin=(1.0, 0.5, 0.13), hubdeg=30):
    rng = np.random.default_rng(seed)
    A, lab = sbm_adj([5, 20, 80], list(pin), pout, rng)
    n0 = len(A); n = n0 + 1
    A2 = np.zeros((n, n)); A2[:n0, :n0] = A
    lab = np.concatenate([lab, [-2]])                       # -2 = hub
    nb = rng.choice(n0, hubdeg, replace=False); A2[n - 1, nb] = 1; A2[nb, n - 1] = 1
    return A2, lab

A, lab = build(); n = len(A)
cls = np.where(lab == 0, 0, np.where(lab == 1, 1, np.where(lab == 2, 2, 3)))
CNAME = ["5-gene island", "20-gene island", "80-gene island", "hub gene"]
CCOL = [VERM, GREEN, BLUE, BLACK]
ALPHA = 0.05
print("n genes", n, "edges", int(A.sum() / 2), "mean degree", A.sum(1).mean().round(1))

# ---- null: degree-preserving rewiring; per-seed global minimum conductance and all local minima ------------------------------
def local_minima(sz, phi, kmin=3):
    """Local minima of the sweep profile that are also the lowest value within a factor of 2 in size (k/2 to 2k)."""
    out = []
    for i in range(1, len(phi) - 1):
        k = int(sz[i])
        if k < kmin or not (phi[i] < phi[i - 1] and phi[i] <= phi[i + 1]): continue
        w = (sz >= k / 2) & (sz <= 2 * k)
        if phi[i] <= phi[w].min() + 1e-12: out.append((k, float(phi[i])))
    return out

null_ncp = collections.defaultdict(lambda: 9.0); null_min = []
for r in range(3):
    An = rewire(A, np.random.default_rng(10 + r), 5)
    for sd in range(0, n, 2):
        order, sz, phi = sweep(An, ppr(An, sd, ALPHA))
        for s_, p_ in zip(sz, phi):
            if s_ >= 3: null_ncp[s_] = min(null_ncp[s_], p_)
        ok = sz >= 3; null_min.append(phi[ok].min())
gate = float(np.quantile(null_min, 0.01))
print("null: median per-seed minimum conductance", np.median(null_min).round(3), " 1% quantile (gate)", round(gate, 3))

# ---- real graph: all genes are start genes -----------------------------------------------------------------------------------
ncp = collections.defaultdict(lambda: 9.0)
seed_isl = {}                       # seed -> (set, phi) smallest significant local minimum
allsets = []
profiles = {}
for sd in range(n):
    p = ppr(A, sd, ALPHA)
    order, sz, phi = sweep(A, p)
    profiles[sd] = (sz, phi)
    for s_, p_ in zip(sz, phi):
        if s_ >= 3: ncp[s_] = min(ncp[s_], p_)
    mins = [(k, f) for k, f in local_minima(sz, phi) if f < gate]
    if mins:
        k, f = min(mins)                                   # smallest significant scale
        seed_isl[sd] = (set(order[:k].tolist()), f)
        for k2, f2 in mins: allsets.append((set(order[:k2].tolist()), f2))
print("start genes with an island:", len(seed_isl), "of", n)
sizes_seed = np.array([len(s) for s, _ in seed_isl.values()])
kept = greedy_merge(allsets)
print("unique islands after greedy Jaccard merge:", sorted([(len(k[0]), round(k[1], 3), k[2]) for k in kept]))
for k in sorted(kept, key=lambda k: len(k[0])):
    print("  size", len(k[0]), "composition", {int(a): int(b) for a, b in collections.Counter(lab[list(k[0])]).items()})

# ---- depth window for each found island ------------------------------------------------------------------------------------
win = []
def conductance(A, idx):
    m = np.zeros(len(A), bool); m[list(idx)] = True
    cut = A[np.ix_(m, ~m)].sum(); vol = A[m].sum()
    return cut / min(vol, A.sum() - vol)
for c in range(3):
    idx = np.where(lab == c)[0]; sub = A[np.ix_(idx, idx)]
    lam, _, _ = sym_norm_eig(sub)
    gap = 1 - lam[1]; phi_ = conductance(A, idx)
    win.append(dict(size=len(idx), phi=round(float(phi_), 3), gap_within=round(float(gap), 3), t_low=round(1 / gap, 1),
                    t_high_nonlazy=round(1 / phi_, 1), t_high_lazy=round(2 / phi_, 1)))
print("depth windows (planted sets):", win)
json.dump(win, open(os.path.join(HERE, "fig4_windows.json"), "w"), indent=1)

# ---- (a) network + conductance profiles -------------------------------------------------------------------------------------
pos = np.zeros((n, 2)); rng = np.random.default_rng(2)
def disc(idx, c, r):
    ang = rng.uniform(0, 2 * np.pi, len(idx)); rad = r * np.sqrt(rng.random(len(idx)))
    return np.c_[c[0] + rad * np.cos(ang), c[1] + rad * np.sin(ang)]
pos[lab == 2] = disc(np.where(lab == 2)[0], (3.2, -0.6), 2.1)
ang = np.linspace(0, 2 * np.pi, 20, endpoint=False); pos[lab == 1] = np.c_[-2.6 + 0.9 * np.cos(ang), 2.2 + 0.9 * np.sin(ang)]
ang = np.linspace(0, 2 * np.pi, 5, endpoint=False); pos[lab == 0] = np.c_[-3.4 + 0.45 * np.cos(ang), -2.6 + 0.45 * np.sin(ang)]
pos[lab == -2] = [0.0, 0.0]
fig, axs = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1, 1.1]})
ax = axs[0]
draw_edges(ax, A, pos, color="#C8C8C8", lw=0.35, alpha=0.5)
hub = np.where(lab == -2)[0][0]
hn = np.where(A[hub] > 0)[0]
from matplotlib.collections import LineCollection
ax.add_collection(LineCollection([[pos[hub], pos[j]] for j in hn], colors=BLACK, linewidths=0.6, alpha=0.35, zorder=2))
for c in range(4):
    m = cls == c
    ax.scatter(pos[m, 0], pos[m, 1], s=[260 if c == 3 else 55][0], color=CCOL[c], edgecolor="white", linewidth=0.5, zorder=3, label=CNAME[c])
ax.set_aspect("equal"); ax.axis("off"); ax.margins(0.06)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=11.5)
ax.set_title("Toy network: three islands, one hub, noise edges")
ax = axs[1]
seeds = [(2, 0), (12, 1), (60, 2)]
found = []
for sd, c in seeds:
    sz, phi = profiles[sd]
    ax.plot(sz, phi, color=CCOL[c], lw=2.4)
    mins = [(k, f) for k, f in local_minima(sz, phi) if f < gate]
    k, f = min(mins)
    ax.scatter([k], [f], s=170, marker="v", color=CCOL[c], edgecolor="black", zorder=5)
    found.append((k, f, c))
    ax.plot([], [], color=CCOL[c], lw=2.4, label=f"start gene in the {CNAME[c]}: lowest at {k} genes, leak {f:.3f}")
ax.axhline(gate, color=GREY, lw=1.5, ls="--", label="gate: 1% quantile of the same statistic\non a degree-preserving rewired network")
ax.set_xscale("log"); ax.set_xlim(3, 150); ax.set_ylim(-0.02, 1.3); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlabel("size of the set (genes, sorted by walk mass per unit degree)")
ax.set_ylabel("conductance = share of edge\nweight leaving the set")
ax.set_title("Sweep: the island is where conductance is lowest")
ax.legend(loc="upper right", frameon=False, fontsize=10.5, bbox_to_anchor=(1.0, 1.02))
stamp(fig, "106 genes; personalised walk, restart 0.05")
fig.tight_layout()
save(fig, "fig4a_islands.png")

# ---- (b) NCP and size distribution -------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.05, 1]})
ax = axs[0]
ks = np.array(sorted(ncp)); ks = ks[ks <= 150]
ax.plot(ks, [ncp[k] for k in ks], color=BLUE, lw=2.4, label="this network: best set of each size")
kn = np.array(sorted(null_ncp)); kn = kn[kn <= 150]
ax.plot(kn, [null_ncp[k] for k in kn], color=GREY, lw=2.4, ls="--", label="degree-preserving rewired network")
for k in (5, 20, 80):
    ax.axvline(k, color=VERM, lw=1.3, ls=":")
    ax.text(k * 1.04, 0.93, f"planted\n{k}", color=VERM, fontsize=11, va="top")
ax.set_xscale("log"); ax.set_xlim(3, 150); ax.set_ylim(0, 1.0)
ax.set_xlabel("set size (genes)"); ax.set_ylabel("lowest conductance found\nfrom any start gene")
ax.set_title("Community profile: dips at the module scales")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, fontsize=10.5, ncol=1)
ax = axs[1]
grp = np.where(sizes_seed <= 10, 0, np.where(sizes_seed <= 40, 1, 2))
cnts = np.bincount(grp, minlength=3); labs_ = ["about 5 genes", "about 20 genes", "about 80 genes\n(81 to 89 found)"]
ax.bar(np.arange(3), cnts, color=BLUE, edgecolor="black", width=0.55, label="all 106 genes as start genes (one island each)")
usz = np.array(sorted(len(k[0]) for k in kept)); ug = np.where(usz <= 10, 0, np.where(usz <= 40, 1, 2))
ax.bar(np.arange(3), np.bincount(ug, minlength=3), color=ORANGE, edgecolor="black", width=0.25, label="distinct islands after the Jaccard merge")
ax.set_xticks(range(3)); ax.set_xticklabels(labs_, fontsize=12)
for i_, v_ in enumerate(cnts): ax.text(i_, v_ * 1.12, str(v_), ha="center", fontsize=12)
vals = np.unique(sizes_seed); 
ax.set_ylabel("number of start genes / islands")
ax.set_title("Island sizes over all start genes")
ax.legend(loc="upper left", frameon=False, fontsize=10.5)
ax.set_yscale("log"); ax.set_ylim(0.7, 600)
stamp(fig, "same network")
fig.tight_layout()
save(fig, "fig4b_ncp_sizes.png")
print("size counts", cnts.tolist())

# ---- (c) start-gene sampling -------------------------------------------------------------------------------------------------
rng = np.random.default_rng(7)
iu = np.triu_indices(n, 1)
S = rng.standard_normal(len(iu[0]))
li, lj = lab[iu[0]], lab[iu[1]]
same = (li == lj) & (li >= 0)
for m, b in [(0, 3.0), (1, 2.0), (2, 1.0)]:
    S[same & (li == m)] += b
o = np.argsort(-S); I = iu[0][o]; J = iu[1][o]
R, M = 1500, 25
r = np.arange(1, R + 1); P = (1 / r) / np.sum(1 / r)            # log-uniform in rank = density proportional to 1/rank
q = np.zeros(n); np.add.at(q, I[:R], P); np.add.at(q, J[:R], P)
pi = 1 - (1 - q) ** M                                            # inclusion probability of a gene in one sample of M pairs (both endpoints)
NR = 400; cnt = np.zeros(4); cntw = np.zeros(4)
for _ in range(NR):
    rr = rng.choice(R, size=M, p=P); sd = np.unique(np.concatenate([I[rr], J[rr]]))
    cnt += np.bincount(cls[sd], minlength=4) / len(sd)
    w = 1 / pi[sd]; cntw += np.bincount(cls[sd], weights=w, minlength=4) / w.sum()
truth = np.bincount(cls) / n; cnt /= NR; cntw /= NR
pi_cls = [pi[cls == c].mean() for c in range(4)]
print("inclusion prob by class", np.round(pi_cls, 3)); print("shares truth", truth.round(3), "unweighted", cnt.round(3), "HT", cntw.round(3))
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1, 1.3]})
ax = axs[0]
ax.bar(range(3), pi_cls[:3], color=CCOL[:3], edgecolor="black", width=0.6)
ax.set_xticks(range(3)); ax.set_xticklabels(["5-gene\nisland", "20-gene\nisland", "80-gene\nisland"], fontsize=12)
ax.set_ylabel("chance that a member gene becomes\na start gene in one sample")
ax.set_title("Not every gene is equally likely to be drawn")
for i, v in enumerate(pi_cls[:3]): ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=12)
ax.set_ylim(0, 0.62)
ax = axs[1]
w_ = 0.26; x = np.arange(3)
ax.bar(x - w_, truth[:3], w_, color="#BBBBBB", edgecolor="black", label="truth: every gene counted once")
ax.bar(x, cnt[:3], w_, color=ORANGE, edgecolor="black", label="log-uniform tail sample, unweighted")
ax.bar(x + w_, cntw[:3], w_, color=BLUE, edgecolor="black", label="same sample, inverse-probability weights")
ax.set_xticks(x); ax.set_xticklabels(["5-gene island", "20-gene island", "80-gene island"], fontsize=12)
ax.set_ylabel("share of start genes in each island")
ax.set_title("Unweighted shares are off; weights repair them")
ax.legend(loc="upper left", frameon=False, fontsize=10.5)
ax.set_ylim(0, 1.05)
stamp(fig, "25 pairs per sample, 400 repeats")
fig.tight_layout()
save(fig, "fig4c_seeds.png")
json.dump({"pi_by_class": [float(v) for v in pi_cls], "truth": truth.tolist(), "unweighted": cnt.tolist(), "HT": cntw.tolist(), "gate": gate},
          open(os.path.join(HERE, "fig4_numbers.json"), "w"), indent=1)
