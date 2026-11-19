"""Toy expression simulation, the attention operator (a numpy port of R/clr_attention.R .softmax_rows + lazy walk), and PA helpers.
All toy and simulated."""
import numpy as np

def simulate_modules(G, E, sizes, loads, rng):
    """x_g(e) = load_m * z_m(e) + noise (unit variance noise). Returns X (G x E), module labels (-1 = none)."""
    lab = -np.ones(G, int)
    start = 0
    for m, s in enumerate(sizes):
        lab[start:start + s] = m; start += s
    Z = rng.standard_normal((len(sizes), E))
    X = rng.standard_normal((G, E))
    for m, b in enumerate(loads):
        X[lab == m] += b * Z[m]
    return X, lab

def corr_eigs(X):
    C = np.corrcoef(X)
    return np.sort(np.linalg.eigvalsh(C))[::-1], C

def permute_genes(X, rng):
    Xp = X.copy()
    for g in range(X.shape[0]):
        Xp[g] = rng.permutation(Xp[g])
    return Xp

def clr_like(C):
    """Stand-in for the CLR score: row z-scores of the correlation (toy: correlation, not mutual information),
    negatives set to 0, Stouffer-combined over the two directions, zero diagonal (as in sk_a3r / rz_zmatrix)."""
    R = C.copy(); np.fill_diagonal(R, np.nan)
    Z = (R - np.nanmean(R, 1, keepdims=True)) / np.nanstd(R, 1, ddof=1, keepdims=True)
    Z = np.nan_to_num(Z); Z[Z < 0] = 0
    S = (Z + Z.T) / np.sqrt(2); np.fill_diagonal(S, 0)
    return S

def softmax_rows(S, keff=10, cap=50):
    """Port of clr:::.softmax_rows: softmax over each row's top-`cap` positive scores, temperature set (bisection on log T)
    so that exp(entropy) = keff. Rows with <= keff candidates get uniform weights."""
    G = S.shape[0]; A = np.zeros((G, G))
    for i in range(G):
        row = S[i].copy(); row[i] = -np.inf
        top = np.argsort(-row)[:min(cap, G - 1)]
        top = top[row[top] > 0]
        if len(top) == 0: continue
        s = row[top]
        if len(top) <= keff:
            A[i, top] = 1.0 / len(top); continue
        def neff(lt):
            w = np.exp((s - s.max()) / np.exp(lt)); w /= w.sum(); w = w[w > 0]
            return np.exp(-(w * np.log(w)).sum())
        lo, hi = np.log(1e-4), np.log(1e4)
        for _ in range(60):
            mid = (lo + hi) / 2
            if neff(mid) < keff: lo = mid
            else: hi = mid
        w = np.exp((s - s.max()) / np.exp((lo + hi) / 2)); A[i, top] = w / w.sum()
    return A

def operator_from_corr(C, keff=10, cap=50):
    return softmax_rows(clr_like(C), keff, cap)

def sym_lazy_eigs(A):
    """Symmetrise the attention matrix (A + A^T)/2, normalise by degree, lazy eigenvalues lam = (1 + mu)/2 (descending)
    and eigenvectors of the symmetric normalised matrix. The real operator is directed; this is the review's 'symmetrised' reading."""
    As = (A + A.T) / 2
    d = As.sum(1); d[d == 0] = 1
    M = As / np.sqrt(np.outer(d, d))
    mu, U = np.linalg.eigh(M)
    o = np.argsort(-mu)
    return (1 + mu[o]) / 2, U[:, o], d
