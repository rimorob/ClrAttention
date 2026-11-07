"""Faithful Python port of the DREAM5 Challenge 4 MATLAB scoring code (Evaluation scripts/matlab/*.m).
Run with python3 -I.  Usage: d5_official_port.py <scripts_dir> <team_dir_or_none>
Lines are cited as file:line of the MATLAB sources.  Used here only to test the port on the SAMPLE 'myteam' files
shipped inside the archive (compared with the shipped OUTPUT/myteam.txt)."""
import sys, numpy as np, scipy.io as sio


def load_dream_network(path):
    """load_dream_network.m: strips the leading 'G', returns int arrays and confidence."""
    I, J, C = [], [], []
    for l in open(path):
        f = l.rstrip("\r\n").split("\t")
        if len(f) < 2 or not f[0]:
            continue
        I.append(int(f[0][1:])); J.append(int(f[1][1:])); C.append(float(f[2]) if len(f) > 2 and f[2] != "" else np.nan)
    return np.array(I), np.array(J), np.array(C)


def evaluate(gold, pred, pdf_aupr, pdf_auroc):
    gi, gj, gc = gold
    keep = gc != 0                                              # DREAM5_Challenge4_Evaluation.m:4-5
    gi, gj = gi[keep], gj[keep]
    regs = np.unique(gi)                                        # edgelist_2_D5C4_gs.m:3
    tgts = np.unique(np.concatenate([gi, gj]))                  # edgelist_2_D5C4_gs.m:4  (unique over BOTH columns)
    posset = set(zip(gi.tolist(), gj.tolist()))                 # A = edgelist2sparse (:7)
    P = len(posset)
    T = len(regs) * len(tgts) - len(np.intersect1d(regs, tgts))  # all (reg, tgt) cells except i == j (:15-19)
    N = T - P
    regset, tgtset = set(regs.tolist()), set(tgts.tolist())
    pi, pj, _ = pred
    disc = []
    for a, b in zip(pi.tolist(), pj.tolist()):                  # remove_edges_not_in_gs.m:10-20: file order kept, no dedupe
        if a in regset and b in tgtset and a != b:
            disc.append(1.0 if (a, b) in posset else 0.0)
    disc = np.array(disc)
    L = len(disc)
    TPL = disc.sum()
    p = (P - TPL) / (T - L) if L < T else 0.0                   # :33-37
    pos = np.concatenate([disc, np.full(T - L, p)])             # :40-45
    neg = np.concatenate([1 - disc, np.full(T - L, 1 - p)])
    TPk = np.cumsum(pos); FPk = np.cumsum(neg)
    K = np.arange(1, T + 1)
    TPR = TPk / P; FPR = FPk / N; REC = TPR; PREC = TPk / K     # :51-57
    trapz = getattr(np, "trapezoid", None) or np.trapz
    AUROC = trapz(TPR, FPR)                                     # :76 trapz(FPR,TPR)
    AUPR = trapz(PREC, REC) / (1 - 1.0 / P)                     # :104
    def prob(pdf, x):                                           # probability.m: P(X>=x) = sum((X>=x).*Y*dx)
        X = pdf["X"].ravel(); Y = pdf["Y"].ravel(); dx = X[1] - X[0]
        return float(np.sum((X >= x) * Y * dx))
    return dict(P=P, N=N, T=T, L=L, TPL=int(TPL), AUROC=AUROC, AUPR=AUPR, P_AUROC=prob(pdf_auroc, AUROC), P_AUPR=prob(pdf_aupr, AUPR))


def main(root, team):
    res = {}
    for n in (1, 3, 4):
        gold = load_dream_network(f"{root}/INPUT/gold_standard_edges_only/DREAM5_NetworkInference_Edges_Network{n}.tsv")
        pred = load_dream_network(f"{root}/INPUT/predictions/{team}/DREAM5_NetworkInference_{team}_Network{n}.txt")
        pa = sio.loadmat(f"{root}/INPUT/probability_densities/Network{n}_AUPR.mat")
        po = sio.loadmat(f"{root}/INPUT/probability_densities/Network{n}_AUROC.mat")
        res[n] = evaluate(gold, pred, pa, po); print(n, res[n])
    aupr_score = -np.mean(np.log10([res[n]["P_AUPR"] for n in res]))      # go_all.m:66
    auroc_score = -np.mean(np.log10([res[n]["P_AUROC"] for n in res]))    # :67
    print("OVERALL %.3f AUPR_score %.3f AUROC_score %.3f" % ((aupr_score + auroc_score) / 2, aupr_score, auroc_score))  # :68


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
