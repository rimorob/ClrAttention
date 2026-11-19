### Results since last night (independently QC-checked at 06:50: all numbers confirmed)

| Job | Headline |
|---|---|
| 087 BEELINE, all expressed genes as CLR context | No gain from a wider context. Median AUPRC ratio CLR (MI) 1.579 -> 1.532, CLR-NS 1.654 -> 1.641. The nominally significant drop for CLR (p = 0.02) is driven by one dataset (mHSC-L, two sparse STRING truths); without it p = 0.19. |
| 088 BEELINE, driven-gene context (multimodal or PC-driven, BIC) | Better than the variance-ranked 5,000-gene context, not better than BEELINE's own gene sets: CLR 1.561, CLR-NS 1.640. |
| All contexts | CLR-NS beats MI-based CLR in every context (+0.07 to +0.12; 53-57 of 72; p <= 0.002). Attention on CLR-NS ties CLR-NS. |
| Caveats (QC) | The 72 combinations come from 7 datasets (paired tests overstate the evidence; a blocked test is needed for the paper). "All genes" equals "5,000" for the three mHSC datasets (only about 4,750 genes in their files). |
| Reading | Single-cell datasets explore few network states (your point): a fuller CLR background adds no regulatory signal there. BEELINE stays a secondary result; the human perturbation compendia carry the human story. |

### Decisions for you

None blocking. One to consider: the paper's single-cell section could drop the context comparison to a supplementary sentence ("no gain from wider contexts") rather than a figure.

### Queue problem found and fixed this morning

The runner built its job list once per pass, so the three jobs added last night (0871, 0872 compendia, 0873 granularity) would have waited behind the 2-3-day FDR bootstrap. Fixed without interrupting the running job: 090/095/100 held until 0871 starts, then returned to the queue; the runner script now rescans after every job (applies at next runner restart).

### Overnight engineering (no results changed)

| Item | Status | Detail |
|---|---|---|
| Packed (tiled) MI -> CLR in the C++ core | Done, tested, installed | Built on your original `miSubMarix` band primitive (no scheduler). Double-precision path bit-identical to the dense path; all 57 package tests pass. |
| Benchmark, 8,563 genes x 707 profiles | 5x less memory | Peak 3.77 GB -> 0.76 GB, 53 s -> 47 s; identical attention operator. |
| BLAS threads, beeline.R memory, restricted attention walks | Done | Identical results; used by 088 onward. |

### Granularity (is depth a resolution parameter?)

Not answered yet: the sweep (0873) runs after the two compendia, early this afternoon.

### Expected today

089 (lab null without perturbation experiments) until about 10:30; then 0871 genome-wide K562 and 0872 RPE1 compendia (about 1 h each), 0873 granularity sweep (about 1 h, early afternoon), then 090 BEELINE bootstrap, 095 SSEM-Lasso, 100 FDR bootstrap.
