# Short pre-read for Jim Collins: internal history notes

**Internal notes, not for circulation.** This file holds what was removed from `docs/jim_short_preread.docx` so that the circulated document reads as a statement of the project as it stands, with no versioning. Written 16 Nov 2026.

## Where the document and its pieces live

- Deliverable: `docs/jim_short_preread.docx` (landscape; main body = pages 1 to 3, then appendix A to G).
- Source and build: `docs/jim_short_preread_src/jim_short_preread.md` (tokens `{{ACLR_NAME}}`, `{{LEAD_NAME}}`, `{{LEAD_DESC}}`, `{{PAGEBREAK}}`) and `docs/jim_short_preread_src/build_short_preread.py` (pandoc + post-processing; reads the landscape reference docx `export/ref/landscape.docx` of the session scratchpad; autocrops the three graphical-abstract pngs). Names live in one block of the build script and once in each figure script (`ACLR_NAME`, `LEAD_NAME`).
- Figure scripts: `analysis/fig_graphical_abstract_short.R` (Figure 1), `..._short_p2.R` (Figure A1), `..._short_p3.R` (Figure 2 raw scoreboard and Figure A3 hub-controlled). Outputs in `docs/whitepapers/figures/graphical_abstract_short_p1|p2|p3|p3b.{png,pdf}`. Sources of the originals they were copied from: `analysis/fig_graphical_abstract.R`, `_p2.R`, `_p3.R`.
- Depth plot: `results/granularity_log/depth_log_curves.png` (unchanged); source tables `docs/depth_log_results.md`, `docs/preregistration_depth_log.md`.
- The longer first version of the short pre-read (14 pages, with sections "The method in one paragraph", "The evidence", "Losses and caveats", "Possible extension", "Decisions for you and me", appendices A and B) is in git at commit 7873398.

## Method names (history)

Working names during drafting: "attended CLR" (PI, 15 Nov, briefly), then "a3r" as a neutral label, then the PI's names of 15 Nov: **aCLR** = CLR + attention (CLR, softmax, lazy walk); **tcaCLR** = "twice cooked aCLR" = aCLR, then the log of the walk P^9, then a second CLR. Internal name of tcaCLR in the repository, result folders and the scoreboard CSV: **a3r** (stacked variants a3r2, a3r3 are dropped from the short pre-read). Names are temporary; the circulated document says so in a footnote and on every figure.

## Dates and status wording removed from the document

| Removed | Where it was | Replaced by |
|---|---|---|
| Title block date "15 November 2026"; author line | front matter | title only |
| "26 October" (preregistered log-then-CLR test; complex superiority over aCLR at t = 9 needed one-sided p < 0.025; M3D +0.026 p = 0.079, EcoMAC +0.029 p = 0.051; regulon non-inferiority passed +0.004 and +0.003; verdict "not adopted"; K562 7.95 selected not fitted) | section 2.5 of the first version | "a preregistered test"; status: "lead candidate; did not reach the preregistered threshold" |
| "13 Nov" yeast estimate (plan fixed 13 Nov; PI decision 12 and 13 Nov that a3r is not adopted and that the study is an estimation, not a test) | section 2.6, scoreboard footnote | "exploratory estimate", "plan fixed in advance" |
| "11 Nov" DREAM5 yeast test (design frozen 4 Nov) | losses section | "preregistered DREAM5 yeast test" |
| "2 Nov" scPRINT comparison (preregistered, run once) | caveats | "preregistered, run once" |
| "not adopted under the preregistered rule", "main lead, not yet adopted", "earlier-adopted method" | text, figures | "lead candidate" (tcaCLR); "adopted" label kept for aCLR in figures only |
| "a3r in the repository" | names paragraph | removed (kept here) |

## Old "Decisions for you and me" list (as it stood in the first short version)

1. Venue: Nature Machine Intelligence (150-word abstract) or Nature Communications.
2. Scope: one paper on aCLR (with tcaCLR), with perturbation targets and differential CLR as a possible second paper.
3. The central claim: "better recovery of regulatory programs from expression alone", not "better network inference"; is that strong enough for the target venue? A statistical review suggests scoping the claim to *E. coli*, or presenting the work as a reframing of CLR with characterised trade-offs, since every positive result comes from the development organism and lineage.
4. The human single-cell story: CLR-NS rather than aCLR is the single-cell result; does it belong in this paper?
5. The yeast result: co-membership estimate positive but under the 10% margin, edge read-out inconclusive, E. coli network 3 pattern (module gain, edge cost) not independent of the headline; how should it enter the paper and abstract?
6. What to do with tcaCLR: headline method aCLR, tcaCLR, or both with trade-offs? Not decided; names temporary.
7. Whether the yeast estimate is the last evaluation of tcaCLR on a new complex truth. Not decided; the yeast Complex Portal truth is consumed.
8. AI disclosure: code, runs and drafts produced with Claude (Anthropic), independent AI statistical reviews, under the author's direction; contributions statement drafted for Nature Portfolio policy (kept as appendix G).

The circulated document keeps only three decision points (headline method; whether a further tcaCLR evaluation on a new complex truth is warranted; naming) in "Potential follow-ups and decision points". Items 1 to 5 (venue, scope, central claim, single-cell story, how the yeast result enters the paper) are not in the circulated document.

## Repo pointers and commit ids removed from the document

- Commit ids cited on the old pre-read page: depth log-scale re-test preregistration 8861ae4 (11 Nov); stacked-layer plan 98f3b0d; yeast complex design rev 2 099c26d; log-scale depth results commit 4460846.
- Files cited: `docs/whitepapers/sota_table.md` (full method table with every truth and footnotes), `docs/depth_log_results.md`, `docs/preregistration_depth_log.md`, `docs/plans/corum_analysis.md` (reference code for the tcaCLR read-out, ~line 127), `R/clr_attention.R` (`.softmax_rows`, `build_operator`), `docs/plans/complex_confirm_design.md`, `docs/plans/complex_confirm_review.md`, `docs/complex_confirm_code_review.md`, `docs/explore_stack_results.md`, `docs/plans/scoreboard_plan.md`, `docs/plans/stats_review_2026-11-12.md`, `docs/clr_logpt_verdict/`, `docs/keff_verdict/`, `docs/scfm_scoring.md`.
- Verbatim sentences the yeast study script wrote (not in the circulated document, which paraphrases the first three): "The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test)."; "The whole interval lies above the smallest effect of interest (+0.005)."; "The whole interval lies above the E. coli estimate (+0.027)."; "Exploratory estimation, not a confirmatory test, a replication or a validation. a3r is not adopted; this study cannot change that by itself (PI decisions 12 and 13 Nov)."
- Process facts dropped from appendix D: the code review of the code that ran was done after the run, "not before it as the protocol required" (the reviewer seeing the estimate by accident is kept).

## Which sections moved where (first short version -> circulated version)

| First version | Circulated version |
|---|---|
| 1. The method in one paragraph (+ Names paragraph) | Page 1 paragraph + footnote |
| Figure 1 graphical abstract (long caption) | Page 1, Figure 1 (short caption) |
| 2.1 Methods compared (Figure 2) | Appendix A, Figure A1 |
| 2.2 Scoreboard (Figure 3) | Page 2, Figure 2 |
| 2.3 aCLR against established methods (table, headline bullets, why it works, controls) | Appendix A.1 (table A1) and A.2 |
| 2.4 Walk depth (Figure 4, table) | Appendix B (Figure A2, table B1) |
| 2.5 Preregistered test of tcaCLR on E. coli | One sentence on page 2; full in appendix C |
| 2.6 Exploratory yeast estimate | One sentence on page 2; full in appendix D |
| 3. Losses and caveats (bulleted) | Page 3 paragraph; remaining caveats in appendix A.3 |
| 4. Possible extension | Page 3, "Potential follow-ups and decision points" |
| 5. Decisions for you and me | Page 3 (three neutral decision points); full old list above |
| Appendix A pseudocode (with "Reference code" line) | Appendix E, repo-path line removed |
| Appendix B hub-controlled scoreboard | Appendix F, Figure A3 |
| (AI disclosure sat in decision item 8) | Appendix G |

## Content of the older pre-read page not carried into the short document at all

Update notes by date, v14 to v20 change logs, follow-ups A to D framing, the ledger of every preregistered test and the "many looks" paragraph (one sentence of it survives in appendix A.3), the "what did not work" list, the in-progress list, robust row-standardisation (median/IQR) follow-up, stacking (a3r2, a3r3) results, the batch-effects supplement, networks 1 and 3 descriptive yeast results beyond one line, the official DREAM5 challenge score ranking, and the PIDC tie-block deviation. The Euclidean "ranks the same" claim from the old figure 1 panel 2 was dropped because it could not be verified from the sources; the figure now shows the two combinations without any claim about their origin or equivalence.
