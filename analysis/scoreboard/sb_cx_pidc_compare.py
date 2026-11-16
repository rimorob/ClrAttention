#!/usr/bin/env python3
"""Bit-for-bit comparison of tools/pidc_chunked.jl output (packed upper-triangle Float64) with an edge file written by the package's own
serial run (tools/pidc.jl: label1 label2 weight, both orientations). Every pair weight must be identical (tolerance 0).
Usage: sb_cx_pidc_compare.py expression.tsv orig_edges.tsv packed.f64 out.csv label   (exit 0 only if identical)"""
import sys, csv, numpy as np, pandas as pd
expr, orig, packed, outcsv, label = sys.argv[1:6]
genes = pd.read_csv(expr, sep='\t', usecols=[0]).iloc[:, 0].astype(str).tolist(); G = len(genes)
W = np.fromfile(packed, dtype='<f8'); assert W.size == G * (G - 1) // 2, (W.size, G)
e = pd.read_csv(orig, sep='\t', header=None, names=['a', 'b', 'w'], dtype={'a': str, 'b': str, 'w': np.float64}, float_precision='round_trip')
idx = pd.Series(np.arange(G), index=genes)
i = idx[e.a].to_numpy().astype(np.int64); j = idx[e.b].to_numpy().astype(np.int64); m = i < j
i = i[m]; j = j[m]; w = e.w.to_numpy()[m]
pos = i * G - (i + 1) * i // 2 + (j - i) - 1          # 0-based index of pair (i, j), i < j
ref = W[pos]; n = len(w)
ident = int((ref == w).sum()); mx = float(np.abs(ref - w).max())
ok = (ident == n) and (n == W.size) and len(np.unique(pos)) == n and not np.isnan(W).any()
print(f"{label}: genes {G}, pairs in the original file {n}, pairs in the new file {W.size}, identical {ident}, differing {n - ident}, max abs diff {mx}, PASS={ok}")
with open(outcsv, 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(['label', 'genes', 'pairs_original', 'pairs_new', 'identical', 'differing', 'max_abs_diff', 'pass'])
    wr.writerow([label, G, n, W.size, ident, n - ident, mx, ok])
sys.exit(0 if ok else 1)
