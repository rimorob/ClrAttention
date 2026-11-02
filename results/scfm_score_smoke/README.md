Smoke test of analysis/scfm_score.R + scfm_score_stats.R on SYNTHETIC ground truth only (2026-11-02, job 1011); no real truth file was read.
Truths were planted from the networks by analysis/scfm_score_synth_truth.R:
- mixed: celltype_ChIP and STRING_regulatory from CLR-NS-like scores, nonspecific_ChIP and STRING_curated from scPRINT R0 in the orientation S = t(A), STRING from both, LOFGOF random.
  Expected and observed: d_i > 0 for the CLR-NS-planted truths (HL +1.85), d_i < 0 for the R0-planted truth (HL -0.54), S = t(A) beats S = A on the R0-planted truth.
- null: all random edges. Expected and observed: all d_i near 0 (primary HL +0.02), AUPRC ratios about 1.
The exact interval and the one-sided p were checked against wilcox.test(conf.int = TRUE, exact = TRUE) and against the Walsh averages (identical).
