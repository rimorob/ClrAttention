# Translational symmetry (arXiv 2602.15029) and its relevance to association networks

## What the paper claims
The paper (Karkada, Korchinski, Nava, Wyart, Bahri; "Symmetry in language statistics shapes the geometry of model representations") is about language models, not networks or biology. Its claim is that word embeddings are governed by a normalised co-occurrence (PMI-like) matrix M*, with M*_ij = (P_ij - P_i P_j) / (0.5 (P_ij + P_i P_j)), roughly log P_ij/(P_i P_j). Embeddings are its top eigenmodes, W = Phi sqrt|Lambda|.

Assumption 3.1 is translation symmetry. Each word sits at a coordinate x_i on a lattice of semantic continuum (months on a circle, years on a line, places on a map), and M*_ij = C(dist(x_i, x_j)). Calendar months co-occur according to their time separation only.

On a periodic lattice M* is circulant, and so diagonalised by the discrete Fourier transform. The eigenvalues are the Fourier coefficients of the kernel C. Embeddings are then sinusoids along each principal component, with amplitudes and frequencies predicted analytically. Circles are the slowest nonconstant mode pair. "Ripples" are higher harmonics. Open boundaries (years) give a different, exponential-kernel solution, and a linear probe decodes the coordinate (Proposition 4).

Section 4 shows robustness: if a large vocabulary shares a latent periodic variable t, with P(i|t) = P(i)(1 + g(t - t_i)), then the seasonal structure appears as a few large eigenvalues of the full matrix. The embedding of months survives deleting all month-month co-occurrences, by Weyl/Davis-Kahan, because the eigenvalue gaps grow with vocabulary size N.

The mathematics is standard Toeplitz/circulant spectral theory plus a latent-variable model. The contribution is applying it to language. For us the relevant part is the recipe: a latent coordinate, a kernel that depends only on coordinate difference, a Fourier spectrum, and recovery through many indirectly related items.

## (a) Does this have a counterpart in our matrices?
Only conditionally, and mostly as metaphor.

- Lattice: there is none given. It would have to be a latent coordinate that genes occupy, for instance cell-cycle phase (G1/S, G2/M, M/G1 modules ordered around a circle), circadian phase, or developmental pseudotime (open line). Genes have a preferred phase t_i, the direct analogue of the seasonal centre.
- Translation: a shift of phase. Symmetry means association depends only on |t_i - t_j| (mod T). That is a circulant or Toeplitz structure after ordering genes by t_i. It is real for periodic transcription (yeast cell cycle, circadian); there the paper's latent-variable model transfers almost literally, with expression as P(i|t).
- Unit cell: one period T of the coordinate, or the repeated module. It is not a unit cell in the crystal sense; there is no spatial repetition of identical genes.
- Metaphor: a graph automorphism or exact symmetry of the gene network is not expected. Complexes, regulons and paralog families are not translation-related; they are cliques or stars, nearly the opposite of a kernel in distance. Gene-family repeats give approximate block structure, not translation. Our CLR, aCLR and tcaCLR matrices are driven by degree-heterogeneous, sparse, clustered structure, so exact circulance would be a surprising finding.
- A caution: the walk P^9 on a ring with a smooth kernel is itself circulant and has eigenvalues cos(k)-like and decaying. So the operator will show Fourier-like spectra if a ring exists, but the same smoothing makes any 1D manifold look "mode-like". Smoothing is not evidence.

## (b) Falsifiable tests
Undirected scores cannot encode directed cycles. CLR, aCLR (after symmetrisation) and tcaCLR are symmetric, so their eigenvalues are real and a feedback loop A->B->C->A is indistinguishable from the reverse loop or from a triangle of co-regulation. Symmetric matrices can reveal cyclic geometry (a ring of phases), not cyclic causality (a feedback loop). Say so to the PI.

Tests for cyclic geometry on the symmetric matrices:
1. Spectral: look for degenerate eigenvalue pairs (cos/sin) in the walk operator, with eigenvector pairs that trace a closed curve in the plane and whose angular coordinate is monotone in recovered phase. Compare with a degree-preserving or configuration-model null, and a rewired null that keeps the top-50 partner structure.
2. Seriation then circulant fit: order genes by the angle from the top eigenvector pair, then fit a symmetrised-Toeplitz/circulant model; report the fraction of off-diagonal variance explained against a shuffled-order null.
3. Topology: persistent homology H1 on the gene graph (1 - tcaCLR as distance, or diffusion-map coordinates) and on the cell-by-gene expression manifold; require a long-lived H1 bar exceeding the null bar distribution from permuted or degree-matched data.
4. Fourier along recovered phase: expression of known periodic genes against phase, with harmonics 1 to 3.
5. Directed information: the aCLR per-row operator before symmetrisation is a directed walk. Its complex eigenvalues (nonzero imaginary parts) and circulation (the antisymmetric part of the flux, for example the curl of P) can be tested against a null. But CLR z-scores are built from a symmetric MI, so directionality in the row operator arises only from row normalisation and top-50 truncation, i.e. hubs and degree. Expect that any imaginary part there is an artefact of degree asymmetry unless shown otherwise. This is my judgement, not tested.

Feedback loops proper need directional information: time series with lags or RNA velocity (Granger-type or lagged MI), Perturb-seq knockdown effects (A knocked down moves B gives a directed edge; cycles then found in the directed perturbation graph), or interventional data from yeast deletions.

## (c) Datasets
- Best fit for a ring: a yeast cell-cycle time series or synchronised dataset. I did not find one in data/ (data/yeast holds the deleteome, GSE42526/42527 and wild-type controls, none a time course; I did not inspect those matrices in detail, so this is uncertain). It would be a small download, and it would be the right test because periodic genes have known phases (Spellman-type lists).
- Single-cell K562/HepG2 Perturb-seq: cells carry a cell-cycle phase (S/G2M scores), so a ring exists in the cell manifold and is testable (H1 on cells, phase genes). It does not give gene-level translation symmetry in our gene graph unless we use the phase-ordered genes.
- Directed: K562 essential Perturb-seq (knockdown to response) for directed cycles.
- E. coli M3D and Kemmeren: steady-state compendia; no latent circular coordinate expected. Use them as negative controls, not as targets.

Pre-registered label-free first pass: fix the statistic (angular-coverage of the top eigenvector pair, circulant R^2 after seriation, longest H1 bar), the nulls and the thresholds before any run, using no truth sets. Negative control: M3D or Kemmeren, plus shuffled and degree-preserving nulls on every dataset. Positive control: a planted ring (simulate genes with P(i|t) = P(i)(1 + g(t - t_i)) with the same size and noise, run through the full CLR to aCLR to tcaCLR pipeline) to measure how much the pipeline and the P^9 smoothing preserve or destroy the signal. A second planted control is a directed 3-cycle for the pre-symmetrisation operator.

## (d) Recommendation
A short exploratory probe, not a programme, and it should be framed as cyclic geometry, not feedback loops. The paper's relevance is thin: it is a clean spectral-theory template (latent coordinate, kernel, Fourier modes, robustness through many shared items), not a method for networks or cycles. Our matrices lack a given lattice, and the undirected scores cannot hold cycle directions.

Effort: about 2 to 3 days. One day for the planted-ring and planted-3-cycle simulation through the pipeline (answers whether our operators could even see a ring), one day for Perturb-seq cell-cycle H1 and the label-free test on the existing matrices, and up to one day if a yeast cell-cycle time series must be fetched. If the planted ring is not recovered after tcaCLR, stop. Not worth pursuing now as a main line. Speculative points above (artefactual imaginary parts, expected failure on M3D) are my judgement and have not been run.
