# Figures

All figures are generated from cached results by the scripts below (run from the repo root) and share the theme and colour key in `analysis/plot_theme.R`: blue = dependence (MI), light blue = histogram bins, orange = CLR normalization, violet = softmax, magenta = depth/walk, aqua = the full method of this work only, grey = competitors/neutral, hatched or greyed yellow = post-hoc PC removal (outside the main pipeline). Platforms and batches use hues outside this key (red, green, brown, dark grey) plus point shapes.

| Figure | Script | Message |
|---|---|---|
| `graphical_abstract.png/.pdf` | `analysis/fig_graphical_abstract.R` | Published CLR is only the context normalization; adding softmax and a multi-step walk completes it as attention (schematic, with real pooled AUPR bars). |
| `fig1_method.png/.pdf` (main Fig. 1) | `analysis/fig_method.R` (`nice -n 19 Rscript analysis/fig_method.R`) | One traced cell through a real 14-gene toy computed with the package (D, CLR, softmax, lazy walk), each step paired with its transformer counterpart; no learned parameters. |
| `src/fig1_method.txt` | (drawing brief) | Fig. 1 brief (superseded by the script above). |
| `fig_shapley_SC.png/.pdf` (main Fig. 3) | `analysis/fig_shapley.R --evidence SC --out docs/whitepapers/figures` | CLR normalization and softmax drive the coherence gain; depth and CLR drive pooled AUPR (strong/confirmed regulons). |
| `fig_shapley_all.png/.pdf` (Supplementary Fig. 3) | `analysis/fig_shapley.R --evidence all --out docs/whitepapers/figures` | The same decomposition on all-evidence regulons; the caption lists where it disagrees with the strong/confirmed one. |
| `scatter_lexA_stats.csv`, `mi_vs_rho_summary.csv` | written by the two scripts below | n, r, MI per scatter panel; per-panel shares above/below the Gaussian curve. |
| `shapley_SC.csv`, `shapley_all.csv` | copies of `results/shapley/` (`analysis/shapley.R`) | Numbers behind Fig. 3 and Supplementary Fig. 3. |
| `scatter_lexA.png/.pdf` (Supplementary Fig. 1) | `analysis/ecomac_scatter.R --out docs/whitepapers/figures` | Mixing platforms creates the "Christmas tree" and dilutes a known LexA–target relationship. |
| `mi_vs_rho.png/.pdf` (Supplementary Fig. 2) | `analysis/mi_vs_rho.R --out docs/whitepapers/figures` (MI cached in `results/mi_vs_rho/`) | Excess MI over correlation comes from mixing platforms (8.4% of pairs vs 0.2–0.3%). |
