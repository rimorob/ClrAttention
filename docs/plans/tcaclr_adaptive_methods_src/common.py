"""Shared helpers for the explainer figures. Everything is toy and simulated; every figure is seeded."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# Okabe-Ito colour-blind-safe palette
BLUE, ORANGE, GREEN, VERM, SKY, PURPLE, YELLOW, BLACK, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#F0E442", "#000000", "#999999")
MODCOL = [BLUE, ORANGE, GREEN, PURPLE, SKY, VERM]
DIV = LinearSegmentedColormap.from_list("bw_v", [BLUE, "#FFFFFF", VERM])
SEQ = LinearSegmentedColormap.from_list("seq", ["#FFFFFF", SKY, BLUE, "#002b45"])

plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 11.5, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 100, "savefig.dpi": 170, "savefig.bbox": "tight", "font.family": "DejaVu Sans"})

TOY = "toy, simulated"

def stamp(fig, extra="", label=TOY):
    fig.text(0.995, -0.01, label + (" | " + extra if extra else ""), ha="right", va="top", fontsize=10,
             color="#555555", style="italic")

def save(fig, name):
    stamp_done = getattr(fig, "_stamped", False)
    p = os.path.join(FIG, name)
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)

def sbm_adj(sizes, p_in, p_out, rng, extra_nodes=0, p_bg=None):
    """Undirected 0/1 stochastic-block adjacency. Returns A, labels (module index; -1 = background)."""
    labels = np.concatenate([np.full(s, i) for i, s in enumerate(sizes)] + [np.full(extra_nodes, -1)]).astype(int)
    n = len(labels)
    P = np.full((n, n), p_out)
    for i, s in enumerate(sizes):
        idx = np.where(labels == i)[0]
        P[np.ix_(idx, idx)] = p_in[i]
    if extra_nodes and p_bg is not None:
        idx = np.where(labels == -1)[0]
        P[np.ix_(idx, idx)] = p_bg
    U = rng.random((n, n)) < P
    A = np.triu(U, 1)
    A = (A | A.T).astype(float)
    return A, labels

def sym_norm_eig(A):
    """Eigen-decomposition of D^-1/2 A D^-1/2 (symmetric); lazy eigenvalues lam=(1+mu)/2, sorted descending."""
    d = A.sum(1)
    d[d == 0] = 1
    M = A / np.sqrt(np.outer(d, d))
    mu, Phi = np.linalg.eigh(M)
    o = np.argsort(-mu)
    mu, Phi = mu[o], Phi[:, o]
    return (1 + mu) / 2, Phi, d

def lazy_P(A, alpha=0.5):
    d = A.sum(1, keepdims=True)
    d[d == 0] = 1
    return (1 - alpha) * np.eye(len(A)) + alpha * A / d

def circle_layout(labels, rng, centres=None, spread=0.55, bg_spread=2.2):
    """Deterministic cluster layout: each module on a circle around its own centre."""
    mods = sorted(set(labels) - {-1})
    k = len(mods)
    if centres is None:
        centres = [(2.4 * np.cos(2 * np.pi * i / k + np.pi / 2), 2.4 * np.sin(2 * np.pi * i / k + np.pi / 2)) for i in range(k)]
    pos = np.zeros((len(labels), 2))
    for m, c in zip(mods, centres):
        idx = np.where(labels == m)[0]
        ang = np.linspace(0, 2 * np.pi, len(idx), endpoint=False) + rng.uniform(0, 1)
        r = spread * (0.6 + 0.4 * rng.random(len(idx))) * (1 + 0.04 * len(idx) ** 0.5)
        pos[idx] = np.c_[c[0] + r * np.cos(ang), c[1] + r * np.sin(ang)]
    bg = np.where(labels == -1)[0]
    if len(bg):
        ang = rng.uniform(0, 2 * np.pi, len(bg)); r = bg_spread * np.sqrt(rng.random(len(bg)))
        pos[bg] = np.c_[r * np.cos(ang), r * np.sin(ang)]
    return pos

def draw_edges(ax, A, pos, color="#BBBBBB", lw=0.5, alpha=0.6):
    i, j = np.where(np.triu(A, 1) > 0)
    segs = np.stack([pos[i], pos[j]], axis=1)
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, alpha=alpha, zorder=1))
