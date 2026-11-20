# tcaCLR2: definition and reading rules (exploratory; written before any tcaCLR2 score exists)

**Name and status.** tcaCLR2 is tcaCLR with its walk depth set by a data-derived rule instead of the constant t = 9. It is a second, separately named column next to the registered tcaCLR on every dataset. It is not part of the frozen preregistration (tag tcaclr-confirm-freeze): it enters none of H1, H2 or H3, takes no part in any multiplicity correction, and the registered tcaCLR results stay attributed to t = 9, keff = 10, cap 50, clipped.

**Rule (fixed here).**
1. Build the aCLR operator exactly as registered: keff 10, cap 50, lazy weight 0.5.
2. Edge: for 10 per-gene permutations of the expression matrix (each gene's values permuted across experiments independently, seeds fixed in the script), take the largest nontrivial eigenvalue modulus of the symmetrised lazy operator (W = (A + A^T)/2, P_sym = 0.5 I + 0.5 D^(-1/2) W D^(-1/2), as in `analysis/tcaclr_depth/adaptive/spectral_edge.R`) and average over the permutations: edge.
3. Depth: t2 = max(1, round(ln(0.05) / ln(edge))). The convention 0.05 is fixed now; 0.01 is not used.
4. Everything else is the registered tcaCLR: P^t2, the second CLR on log P^t2 over reachable pairs, clipped at 0, Stouffer combination, symmetric read-out.

No truth is read in computing t2. keff and cap are not adapted.

**Reading in advance.** Each dataset reports t2 and the edge. tcaCLR2 is read as "not worse than tcaCLR" on a truth if its pooled co-membership AUPR is at least 0.95 times the registered tcaCLR's AUPR on that truth, with the gene-jackknife interval of the existing scorer shown for the two cells. No superiority claim is made. Because the depth-24 and depth-32 curves of the grid (cap 40) were seen on every spent truth before this note, this reading is descriptive on those truths; a registered non-inferiority test needs a truth that has not been used.

**Wording allowed in documents.** Designed to adapt to the dataset; equivalent to tcaCLR within the stated tolerance on the already-used truths; not tested on fresh data. Not allowed: that it will work, or work better, on a new dataset.
