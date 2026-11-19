"""Toy island (local clustering) helpers: personalised PageRank, conductance sweep, island extraction. Toy only."""
import numpy as np

def ppr(A, seed, alpha):
    """Personalised PageRank vector p = alpha e_s + (1-alpha) p W, W = D^-1 A (solved exactly; toy graphs are small)."""
    n = len(A); d = A.sum(1); d[d == 0] = 1
    W = A / d[:, None]
    e = np.zeros(n); e[seed] = 1
    return np.linalg.solve(np.eye(n) - (1 - alpha) * W.T, alpha * e)

def sweep(A, p, max_frac=0.9, max_size=150):
    """Sort genes by p_j / d_j (walk mass per unit degree); conductance of every prefix set up to max_frac of the volume.
    Returns order, sizes, conductance."""
    d = A.sum(1); d_ = np.where(d == 0, 1, d)
    order = np.argsort(-(p / d_))
    vol_tot = d.sum()
    inset = np.zeros(len(A), bool)
    cut = 0.0; vol = 0.0
    sizes, phi = [], []
    for r, v in enumerate(order):
        # adding v: cut += (edges from v to outside) - (edges from v to inside)
        to_in = A[v, inset].sum(); to_out = d[v] - to_in - A[v, v]
        cut += to_out - to_in
        vol += d[v]; inset[v] = True
        if vol > max_frac * vol_tot or r + 1 > max_size: break
        sizes.append(r + 1); phi.append(cut / min(vol, vol_tot - vol))
    return order, np.array(sizes), np.array(phi)

def best_island(A, seed, alpha, min_size=3):
    p = ppr(A, seed, alpha)
    order, sizes, phi = sweep(A, p)
    ok = sizes >= min_size
    j = np.argmin(np.where(ok, phi, np.inf))
    return set(order[:sizes[j]].tolist()), float(phi[j]), sizes, phi

def rewire(A, rng, nswap_factor=10):
    """Degree-preserving rewiring (double edge swaps) of an undirected 0/1 graph without self loops."""
    A = A.copy(); n = len(A)
    ii, jj = np.where(np.triu(A, 1) > 0)
    edges = list(zip(ii.tolist(), jj.tolist()))
    m = len(edges); done = 0; tries = 0
    while done < nswap_factor * m and tries < 50 * nswap_factor * m:
        tries += 1
        a, b = rng.integers(m, size=2)
        (u, v), (x, y) = edges[a], edges[b]
        if rng.random() < 0.5: x, y = y, x
        if len({u, v, x, y}) < 4 or A[u, y] or A[x, v]: continue
        A[u, v] = A[v, u] = 0; A[x, y] = A[y, x] = 0
        A[u, y] = A[y, u] = 1; A[x, v] = A[v, x] = 1
        edges[a] = (min(u, y), max(u, y)); edges[b] = (min(x, v), max(x, v)); done += 1
    return A

def jaccard(a, b): return len(a & b) / len(a | b)

def greedy_merge(islands, thr=0.5):
    """islands: list of (set, phi). Greedy in order of increasing conductance: an island joins an earlier kept one if Jaccard >= thr."""
    kept = []
    for s, ph in sorted(islands, key=lambda x: x[1]):
        for k in kept:
            if jaccard(s, k[0]) >= thr: k[2] += 1; break
        else:
            kept.append([s, ph, 1])
    return kept
