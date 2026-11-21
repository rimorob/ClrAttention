#!/usr/bin/env Rscript
# scGPT extension (docs/plans/scgpt_extension_plan.md, "Statistics and multiplicity", with Revision 1): statistics over datasets. Reads no truth file: only the long metric tables of
# analysis/scgpt_ext/score_ext.R (results/scgpt_ext/<dataset>/metrics_ext.csv), or, with --legacy 1, the earlier tables (results/scfm_score/metrics_long.csv for hESC and hHep, results/scfm_score_k562/
# metrics_long.csv for K562) to reproduce docs/scfm_scoring.md sections 6 and 7.
#
# Definitions (all from the plan; the interval machinery is copied unchanged from analysis/scfm_score_stats.R, preregistration docs/preregistration_scfm.md):
#  d_i  = for dataset i, the mean over its units (subsamples or cases) and over the primary truths (DoRothEA_AB, CollecTRI, Complexes; equal weight: the mean of the unit x truth cells)
#         of log2(AUPR ratio of the reference method / AUPR ratio of the scGPT readout); positive favours the reference (CLR-family) method.
#  Interval: exact Wilcoxon / Hodges-Lehmann interval [L, U] = the C-th and (M+1-C)-th ordered Walsh averages, M = n(n+1)/2, C the largest integer with P(W <= C-1) <= 0.05
#         (n = 6: 3rd and 19th of 21, coverage 90.6%; n = 5: C = 1, L = min d_i, U = max Walsh average, 93.75%; n <= 4: C = 0, no interval, no decision).
#  Margin delta = 0.15 (log2). p = exact one-sided signed-rank p of H0: median <= -delta (on d_i + delta).
#  Decision table (scPRINT arm of the preregistration, "learned attention" = scGPT readout): NI if L > -delta; superior if L > 0; tie if L > -delta and U < delta; scGPT wins if U < 0
#         (materially if U < -delta); otherwise inconclusive. Attached only to P1 and the two key secondaries (and only where n >= 5 gives an interval).
#  P1 = CLR-NS vs scGPT attention. Key secondaries, Holm 0.05 over two, confirmatory only if P1 shows NI: S1 = CLR-NS vs the best of the three scGPT readouts (attention, token cosine,
#         hidden cosine) per dataset (smallest d_i, chosen on the truths), S2 = tcaCLR vs scGPT attention. Everything else is descriptive (intervals, no claim). tcaCLR2 is exploratory:
#         a column with d_i, t2 and edge, never in a test, a Holm family or a decision label. aCLR (Revision 3: the attention walk of the CLR-NS operator at its label-free held-out depth,
#         without the second CLR) is likewise a descriptive column with d_i, never in a test, a Holm family or a decision label.
# Usage (repo root):
#   Rscript analysis/scgpt_ext/stats_ext.R --in results/scgpt_ext --datasets hESC,hHep,K562,HepG2 --out results/scgpt_ext/stats
#   Rscript analysis/scgpt_ext/stats_ext.R --legacy 1 --out results/scgpt_ext/stats_legacy_check      (reproduces docs/scfm_scoring.md sections 6 and 7; exit 1 on a mismatch)

## ---- exact signed-rank machinery (copied unchanged from analysis/scfm_score_stats.R, "decide" renamed for scGPT) -------------------------------------------------
DELTA <- 0.15
walsh <- function(d) { n <- length(d); w <- outer(d, d, "+") / 2; sort(w[upper.tri(w, diag = TRUE)]) }
c_crit <- function(n) { C <- 0L; while (C + 1L <= n * (n + 1) / 2 && stats::psignrank(C, n) <= 0.05) C <- C + 1L; C }   # largest C with P(W <= C-1) <= 0.05
hl <- function(d) {
  d <- d[is.finite(d)]; n <- length(d); w <- walsh(d); Mw <- length(w); C <- c_crit(n)
  list(n = n, hl = stats::median(w), L = if (C >= 1L) w[C] else NA_real_, U = if (C >= 1L) w[Mw + 1L - C] else NA_real_,
       C = C, coverage = if (C >= 1L) 1 - 2 * stats::psignrank(C - 1L, n) else NA_real_, min = min(d), max = max(d))
}
p_ni <- function(d, delta = DELTA) {          # exact one-sided signed-rank p, H0: median <= -delta (on d + delta)
  x <- d + delta; x <- x[x != 0]; n <- length(x); if (n < 1L) return(NA_real_)
  r <- rank(abs(x)); W <- sum(r[x > 0])
  1 - stats::psignrank(ceiling(W - 1e-9) - 1, n)
}
decide <- function(L, U, delta = DELTA) {
  if (!is.finite(L) || !is.finite(U)) return("no interval (n too small): no decision")
  fl <- c(NI = L > -delta, superior = L > 0, tie = (L > -delta) && (U < delta), scGPT_wins = U < 0, scGPT_wins_materially = U < -delta)
  paste(c(names(fl)[fl], if (!any(fl[c("NI", "superior", "tie", "scGPT_wins")])) "inconclusive"), collapse = " + ")
}

## ---- long tables ----------------------------------------------------------------------------------------------------------------------
# normalised long table: dataset, unit, truth, gset, method, auprc_ratio. Methods: abs_pearson, clr, clr_ns, aCLR, tcaCLR, tcaCLR2, scgpt_attn, scgpt_tok, scgpt_hid, scgpt_tokln, scprint_R0.
SCGPT_RD <- c("scgpt_attn", "scgpt_tok", "scgpt_hid")
load_ext <- function(dir, datasets) {
  do.call(rbind, lapply(datasets, function(d) { x <- utils::read.csv(file.path(dir, d, "metrics_ext.csv"), stringsAsFactors = FALSE)
    data.frame(dataset = d, unit = as.character(x$unit), truth = x$truth, gset = x$gset, method = x$method, auprc_ratio = x$auprc_ratio, epr = x$epr, stringsAsFactors = FALSE) }))
}
legacy_map <- c(attn_last_mean = "scgpt_attn", tok_cos = "scgpt_tok", hid_cos = "scgpt_hid", tok_cos_ln = "scgpt_tokln")
norm_legacy <- function(x, unit, gset_key = NULL) {            # x: a legacy metrics_long table; returns the normalised rows of the sources used here
  m <- ifelse(x$source == "cmp", x$readout, ifelse(x$source == "scgpt_A", legacy_map[x$readout], ifelse(x$source == "scprint_medium_A" & x$readout == "R0" & x$orient == "T", "scprint_R0", NA)))
  keep <- !is.na(m) & (x$source != "scgpt_A" | x$orient %in% c("T", "sym")) & (x$source != "cmp" | x$readout %in% c("clr_ns", "clr", "abs_pearson"))
  data.frame(unit = as.character(unit[keep]), truth = x$truth[keep], gset = if (is.null(gset_key)) x$gset[keep] else gset_key[keep], method = unname(m[keep]),
             auprc_ratio = x$auprc_ratio[keep], epr = x$epr[keep], stringsAsFactors = FALSE)
}
load_legacy <- function() {
  b <- utils::read.csv("results/scfm_score/metrics_long.csv", stringsAsFactors = FALSE); cm <- utils::read.csv("results/scfm_score/cases_meta.csv", stringsAsFactors = FALSE)
  b <- b[b$dataset %in% c("hESC", "hHep") & b$truth %in% c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF"), ]
  # gene-set key -> "S" (scGPT set) / "AS" (three-way set) through cases_meta; "S" and "AS" are the same set in hESC and hHep
  kS <- cm$key_S[match(paste(b$dataset, b$N), paste(cm$dataset, cm$N))]; kAS <- cm$key_AS[match(paste(b$dataset, b$N), paste(cm$dataset, cm$N))]
  out <- list()
  for (g in c("S", "AS")) { sel <- which(b$gset == (if (g == "S") kS else kAS)); x <- b[sel, ]
    m <- ifelse(x$source == "cmp", x$readout, ifelse(x$source == "scgpt_A", legacy_map[x$readout], ifelse(x$source == "scprint_medium_A" & x$readout == "R0" & x$orient == "T", "scprint_R0", NA)))
    keep <- !is.na(m) & (x$source != "scgpt_A" | x$orient %in% c("T", "sym")) & (x$source != "cmp" | x$readout %in% c("clr_ns", "clr", "abs_pearson"))
    out[[g]] <- data.frame(dataset = x$dataset[keep], unit = as.character(x$N[keep]), truth = x$truth[keep], gset = g, method = unname(m[keep]), auprc_ratio = x$auprc_ratio[keep], epr = x$epr[keep], stringsAsFactors = FALSE) }
  bee <- do.call(rbind, out)
  k <- utils::read.csv("results/scfm_score_k562/metrics_long.csv", stringsAsFactors = FALSE); k <- k[k$gset %in% c("S", "AS"), ]
  kk <- norm_legacy(k, k$subsample); kk$dataset <- "K562"; kk <- kk[, names(bee)]
  # K562 gene sets: S (scGPT), AS (three-way); the cmp rows exist on S; scPRINT R0 on AS; scGPT readouts on both
  rbind(bee, kk)
}

## ---- contrasts ----------------------------------------------------------------------------------------------------------------------
# cells of a contrast: join reference and candidate rows on dataset, unit, truth; lr = log2(ratio_ref / ratio_cand)
cells <- function(L, ref, cand, truths, gset_ref = "S", gset_cand = gset_ref, datasets = unique(L$dataset)) {
  a <- L[L$method == ref & L$gset == gset_ref & L$truth %in% truths & L$dataset %in% datasets, c("dataset", "unit", "truth", "auprc_ratio")]
  b <- L[L$method == cand & L$gset == gset_cand & L$truth %in% truths & L$dataset %in% datasets, c("dataset", "unit", "truth", "auprc_ratio")]
  m <- merge(a, b, by = c("dataset", "unit", "truth"), suffixes = c("_ref", "_cand")); m$lr <- log2(m$auprc_ratio_ref / m$auprc_ratio_cand); m
}
dvec <- function(m, datasets) { v <- tapply(m$lr, m$dataset, mean); v[datasets[datasets %in% names(v)]] }
dtab <- function(m, label) {                       # per dataset: n cells, mean AUPR ratios beside d_i
  do.call(rbind, lapply(unique(m$dataset), function(ds) { x <- m[m$dataset == ds, ]
    data.frame(contrast = label, dataset = ds, truth = "mean over truths", cells = nrow(x), ratio_ref = mean(x$auprc_ratio_ref), ratio_cand = mean(x$auprc_ratio_cand), d = mean(x$lr), stringsAsFactors = FALSE) }))
}
dtab_truth <- function(m, label) {
  do.call(rbind, lapply(split(m, list(m$dataset, m$truth), drop = TRUE), function(x)
    data.frame(contrast = label, dataset = x$dataset[1], truth = x$truth[1], cells = nrow(x), ratio_ref = mean(x$auprc_ratio_ref), ratio_cand = mean(x$auprc_ratio_cand), d = mean(x$lr), stringsAsFactors = FALSE)))
}
summ_row <- function(label, d, kind = "descriptive", with_decision = FALSE) {
  h <- hl(d)
  data.frame(label = label, kind = kind, n = h$n, HL = h$hl, L = h$L, U = h$U, coverage = h$coverage, dmin = h$min, dmax = h$max, sd = if (h$n > 1L) stats::sd(d) else NA_real_,
             p_NI = if (kind %in% c("primary", "key secondary")) p_ni(d) else NA_real_,
             decision = if (with_decision) decide(h$L, h$U) else NA_character_,
             values = paste(sprintf("%s=%+.3f", names(d), d), collapse = "; "), stringsAsFactors = FALSE)
}
fmt <- function(x, k = 3) ifelse(is.finite(x), sprintf(paste0("%+.", k, "f"), x), "NA")

## ---- the statistics -------------------------------------------------------------------------------------------------------------------
run_stats <- function(L, datasets, truths, fresh = c("HepG2", "Jurkat", "RPE1"), out = NULL, edge = NULL, adepth = NULL, say = function(...) cat(..., "\n")) {
  n <- length(datasets); rows <- list(); dt <- list(); add <- function(...) rows[[length(rows) + 1L]] <<- summ_row(...)
  say(sprintf("== n = %d datasets (%s); primary truths %s; margin delta = %.2f", n, paste(datasets, collapse = ", "), paste(truths, collapse = ", "), DELTA))
  hh <- hl(rep(0, max(n, 1L))); say(sprintf("   interval at n = %d: C = %d, coverage %s, smallest attainable one-sided p = 1/2^n = %.4f%s", n, hh$C,
      if (is.finite(hh$coverage)) sprintf("%.4f", hh$coverage) else "NA", 1 / 2^n, if (hh$C < 1L) "  -> no interval, NO DECISION possible" else ""))
  ctr <- function(ref, cand, gset = "S", tr = truths) cells(L, ref, cand, tr, gset, gset, datasets)
  ## P1
  m1 <- ctr("clr_ns", "scgpt_attn"); d1 <- dvec(m1, datasets); dt$P1 <- dtab(m1, "P1 CLR-NS vs scGPT attention"); dt$P1t <- dtab_truth(m1, "P1 CLR-NS vs scGPT attention")
  h1 <- hl(d1); p1 <- p_ni(d1); dec1 <- decide(h1$L, h1$U); ni1 <- is.finite(h1$L) && h1$L > -DELTA
  say("== P1: CLR-NS vs scGPT attention"); say("   d_i: ", paste(sprintf("%s %s", names(d1), fmt(d1)), collapse = ", "))
  say(sprintf("   HL %s; interval [%s, %s] (coverage %s, C = %d); p_NI = %.5f; rows: %s", fmt(h1$hl), fmt(h1$L), fmt(h1$U), if (is.finite(h1$coverage)) sprintf("%.4f", h1$coverage) else "NA", h1$C, p1, dec1))
  add("P1 CLR-NS vs scGPT attention", d1, "primary", TRUE)
  ## S1: best of three scGPT readouts per dataset (smallest d_i)
  dS1 <- sapply(SCGPT_RD, function(r) dvec(ctr("clr_ns", r), datasets)); if (is.null(dim(dS1))) dS1 <- matrix(dS1, nrow = 1, dimnames = list(datasets, SCGPT_RD))
  best <- apply(dS1, 1, which.min); d_s1 <- setNames(apply(dS1, 1, min), rownames(dS1))
  h_s1 <- hl(d_s1); p_s1 <- p_ni(d_s1)
  say("== S1: CLR-NS vs the best of the three scGPT readouts per dataset (smallest d_i)"); say("   best readout: ", paste(sprintf("%s %s (%s)", rownames(dS1), SCGPT_RD[best], fmt(d_s1)), collapse = ", "))
  say(sprintf("   HL %s; interval [%s, %s]; p_NI = %.5f; rows: %s", fmt(h_s1$hl), fmt(h_s1$L), fmt(h_s1$U), p_s1, decide(h_s1$L, h_s1$U))); add("S1 CLR-NS vs best scGPT readout", d_s1, "key secondary", TRUE)
  ## S2: tcaCLR vs scGPT attention
  has_tca <- any(L$method == "tcaCLR")
  m2 <- ctr("tcaCLR", "scgpt_attn"); d_s2 <- if (has_tca) dvec(m2, datasets) else setNames(rep(NA_real_, n), datasets); if (has_tca) { dt$S2 <- dtab(m2, "S2 tcaCLR vs scGPT attention"); dt$S2t <- dtab_truth(m2, "S2 tcaCLR vs scGPT attention") }
  h_s2 <- if (has_tca) hl(d_s2) else list(hl = NA_real_, L = NA_real_, U = NA_real_); p_s2 <- if (has_tca) p_ni(d_s2) else NA_real_
  say("== S2: tcaCLR vs scGPT attention"); if (!has_tca) say("   (no tcaCLR rows in these tables)") else {
  say("   d_i: ", paste(sprintf("%s %s", names(d_s2), fmt(d_s2)), collapse = ", "))
  say(sprintf("   HL %s; interval [%s, %s]; p_NI = %.5f; rows: %s", fmt(h_s2$hl), fmt(h_s2$L), fmt(h_s2$U), p_s2, decide(h_s2$L, h_s2$U))); add("S2 tcaCLR vs scGPT attention", d_s2, "key secondary", TRUE) }
  ## Holm over S1, S2 (alpha 0.05): thresholds 0.025 then 0.05; confirmatory only if P1 shows NI
  ps <- c(S1 = p_s1, S2 = p_s2); ord <- order(ps); thr <- c(0.025, 0.05); holm <- setNames(logical(2), names(ps))
  for (k in seq_along(ord)) { if (is.finite(ps[ord[k]]) && ps[ord[k]] <= thr[k]) holm[ord[k]] <- TRUE else break }
  attain <- 1 / 2^n <= 0.025
  say(sprintf("== Holm (0.05 over S1, S2): P1 NI = %s -> secondaries %s; smallest attainable p at n = %d is %.4f (%s the 0.025 first-step threshold)", ni1,
              if (ni1 && attain) "CONFIRMATORY" else "DESCRIPTIVE", n, 1 / 2^n, if (attain) "below" else "above"))
  say(sprintf("   S1: p = %.5f, Holm-rejected %s; S2: p = %.5f, Holm-rejected %s; counted as confirmatory only if confirmatory = TRUE: S1 %s, S2 %s", p_s1, holm[["S1"]], p_s2, holm[["S2"]],
              ni1 && attain && holm[["S1"]], ni1 && attain && holm[["S2"]]))
  holm_tab <- data.frame(test = c("P1", "S1", "S2"), n = n, p_NI = c(p1, p_s1, p_s2), NI = c(ni1, is.finite(h_s1$L) && h_s1$L > -DELTA, is.finite(h_s2$L) && h_s2$L > -DELTA),
                         holm_threshold = c(NA, thr[rank(ps)]), holm_rejected = c(NA, holm[["S1"]], holm[["S2"]]), confirmatory = c(NA, ni1 && attain && holm[["S1"]], ni1 && attain && holm[["S2"]]),
                         p_min_attainable = 1 / 2^n)
  ## descriptive: every reference method x scGPT readout (primary-truth d_i), tcaCLR2 included (exploratory column), and per-truth versions
  say("== DESCRIPTIVE (intervals, no claims): method vs scGPT readout, d_i over the primary truths")
  for (ref in intersect(c("abs_pearson", "clr", "clr_ns", "tcaCLR", "tcaCLR2"), unique(L$method))) for (rd in c("scgpt_attn", "scgpt_tok", "scgpt_hid", "scgpt_tokln")) {
    m <- ctr(ref, rd); if (!nrow(m)) next; add(sprintf("%s vs %s", ref, rd), dvec(m, datasets), if (ref == "tcaCLR2") "exploratory (tcaCLR2)" else "descriptive") }
  if (any(L$method == "tcaCLR2")) { m <- ctr("tcaCLR2", "tcaCLR"); if (nrow(m)) { add("tcaCLR2 vs tcaCLR", dvec(m, datasets), "exploratory (tcaCLR2)"); dt$c2 <- dtab_truth(m, "tcaCLR2 vs tcaCLR") }
    m <- ctr("tcaCLR2", "scgpt_attn"); if (nrow(m)) dt$c2a <- dtab(m, "tcaCLR2 vs scGPT attention (exploratory)") }
  # aCLR (Revision 3): descriptive column, d_i and intervals only; no p, no decision label, no Holm row (kind "descriptive (aCLR)")
  if (any(L$method == "aCLR")) {
    for (rd in c("scgpt_attn", "scgpt_tok", "scgpt_hid", "scgpt_tokln")) { m <- ctr("aCLR", rd); if (nrow(m)) add(sprintf("aCLR vs %s", rd), dvec(m, datasets), "descriptive (aCLR)") }
    m <- ctr("aCLR", "scgpt_attn"); if (nrow(m)) { dt$a1 <- dtab(m, "aCLR vs scGPT attention (descriptive)"); dt$a1t <- dtab_truth(m, "aCLR vs scGPT attention (descriptive)") }
    for (pr in list(c("tcaCLR", "aCLR"), c("clr_ns", "aCLR"))) if (any(L$method == pr[1])) { m <- ctr(pr[1], pr[2]); if (nrow(m)) { add(sprintf("%s vs aCLR", pr[1]), dvec(m, datasets), "descriptive (aCLR)"); dt[[paste0("a_", pr[1])]] <- dtab_truth(m, sprintf("%s vs aCLR", pr[1])) } }
    for (tn in truths) { m <- ctr("aCLR", "scgpt_attn", tr = tn); if (nrow(m)) add(sprintf("aCLR vs scgpt_attn, truth %s only", tn), dvec(m, datasets), "descriptive (aCLR)") } }
  for (tn in truths) for (rd in c("scgpt_attn", "scgpt_tok", "scgpt_hid")) for (ref in intersect(c("clr_ns", "tcaCLR"), unique(L$method))) {
    m <- ctr(ref, rd, tr = tn); if (!nrow(m)) next; add(sprintf("%s vs %s, truth %s only", ref, rd, tn), dvec(m, datasets), "per-truth") }
  # leave-one-dataset-out for P1 and S2
  if (n >= 3L) for (omit in datasets) { ds <- setdiff(datasets, omit)
    add(sprintf("P1 CLR-NS vs scGPT attention, leaving out %s", omit), dvec(cells(L, "clr_ns", "scgpt_attn", truths, "S", "S", ds), ds), "leave-one-out")
    if (has_tca) add(sprintf("S2 tcaCLR vs scGPT attention, leaving out %s", omit), dvec(cells(L, "tcaCLR", "scgpt_attn", truths, "S", "S", ds), ds), "leave-one-out") }
  # fresh-only blocks (no interval at n <= 3)
  fr <- intersect(fresh, datasets); if (length(fr) >= 1L) for (pr in Filter(function(p) any(L$method == p[1]), list(c("clr_ns", "scgpt_attn"), c("tcaCLR", "scgpt_attn"))))
    add(sprintf("%s vs %s, fresh-only blocks (%s)", pr[1], pr[2], paste(fr, collapse = ", ")), dvec(cells(L, pr[1], pr[2], truths, "S", "S", fr), fr), "fresh-only")
  S <- do.call(rbind, rows)
  say("== table"); for (k in seq_len(nrow(S))) say(sprintf("  [%s] %-70s n=%d HL %s%s%s", S$kind[k], S$label[k], S$n[k], fmt(S$HL[k]),
      if (is.finite(S$L[k])) sprintf(" [%s, %s]", fmt(S$L[k]), fmt(S$U[k])) else sprintf(" range [%s, %s]", fmt(S$dmin[k]), fmt(S$dmax[k])), if (!is.na(S$decision[k])) paste0("  -> ", S$decision[k]) else ""))
  # tcaCLR2 depth per matrix (exploratory column)
  if (!is.null(edge) && nrow(edge)) { say("== tcaCLR2 (exploratory): t2 and edge per matrix"); for (ds in unique(edge$dataset)) { e <- edge[edge$dataset == ds, ]
    say(sprintf("   %-8s t2 = %s; edge %s", ds, paste(e$t2, collapse = "/"), paste(sprintf("%.4f", e$edge), collapse = "/"))) } }
  if (!is.null(adepth) && nrow(adepth)) { say("== aCLR (descriptive): label-free held-out depth t* per matrix"); for (ds in unique(adepth$dataset)) { e <- adepth[adepth$dataset == ds, ]
    say(sprintf("   %-8s t* = %s", ds, paste(e$t_heldout, collapse = "/"))) } }
  dtt <- do.call(rbind, dt); rownames(dtt) <- NULL
  if (!is.null(out)) { dir.create(out, recursive = TRUE, showWarnings = FALSE); utils::write.csv(S, file.path(out, "contrasts_ext.csv"), row.names = FALSE); utils::write.csv(dtt, file.path(out, "d_i_table.csv"), row.names = FALSE)
    utils::write.csv(holm_tab, file.path(out, "holm_P1_S1_S2.csv"), row.names = FALSE); utils::write.csv(dS1, file.path(out, "S1_d_by_readout.csv")); if (!is.null(edge)) utils::write.csv(edge, file.path(out, "tcaclr2_t2_edge.csv"), row.names = FALSE)
    if (!is.null(adepth) && nrow(adepth)) utils::write.csv(adepth, file.path(out, "aclr_heldout_depth.csv"), row.names = FALSE) }
  invisible(list(S = S, dtab = dtt, holm = holm_tab, d1 = d1, d_s1 = d_s1, d_s2 = d_s2))
}

## ---- reproduction of docs/scfm_scoring.md sections 6 and 7 (--legacy 1) ----------------------------------------------------------------
legacy_check <- function(L, say) {
  fails <- 0L
  chk <- function(label, got, ref, digits) { ok <- round(got, digits) == round(ref, digits); if (!ok) fails <<- fails + 1L
    say(sprintf("  %-4s %-62s computed %s  docs %s", if (ok) "ok" else "FAIL", label, formatC(got, format = "f", digits = digits + 1, flag = "+"), formatC(ref, format = "f", digits = digits, flag = "+"))) }
  say("== Section 6 (n = 2: hESC, hHep): CLR-NS (or scPRINT R0) minus the scGPT readout; HL with the range of the two d_i, 2 decimals")
  six <- list(list("clr_ns", "scgpt_attn", "S", c(0.37, 0.18, 0.57)), list("clr_ns", "scgpt_tok", "S", c(-0.78, -0.91, -0.64)), list("clr_ns", "scgpt_hid", "S", c(-0.29, -0.43, -0.16)),
              list("scprint_R0", "scgpt_attn", "AS", c(0.08, 0.07, 0.09)), list("scprint_R0", "scgpt_tok", "AS", c(-1.07, -1.12, -1.01)), list("scprint_R0", "scgpt_hid", "AS", c(-0.58, -0.63, -0.53)))
  tr4 <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
  for (s in six) { d <- dvec(cells(L, s[[1]], s[[2]], tr4, s[[3]], s[[3]], c("hESC", "hHep")), c("hESC", "hHep")); h <- hl(d); r <- s[[4]]
    chk(sprintf("%s minus %s: HL", s[[1]], s[[2]]), h$hl, r[1], 2); chk(sprintf("%s minus %s: min d_i", s[[1]], s[[2]]), h$min, r[2], 2); chk(sprintf("%s minus %s: max d_i", s[[1]], s[[2]]), h$max, r[3], 2) }
  say("== Section 6, LayerNorm'd token cosine (exploratory): CLR-NS minus it -0.754, range [-0.910, -0.598]")
  d <- dvec(cells(L, "clr_ns", "scgpt_tokln", tr4, "S", "S", c("hESC", "hHep")), c("hESC", "hHep")); h <- hl(d); chk("clr_ns minus scgpt_tokln: HL", h$hl, -0.754, 3); chk("  min", h$min, -0.910, 3); chk("  max", h$max, -0.598, 3)
  say("== Section 7 (K562, 5 subsamples; per truth and mean of the three truths; 3 decimals)")
  k7 <- list(list("clr_ns", "scgpt_attn", "S", c(1.002, 0.712, 0.160, 0.625)), list("clr_ns", "scgpt_tok", "S", c(-0.714, -0.325, 0.096, -0.314)), list("clr_ns", "scgpt_hid", "S", c(0.233, 0.214, 0.298, 0.248)),
             list("scprint_R0", "scgpt_attn", "AS", c(0.722, 0.371, 0.012, 0.369)), list("scprint_R0", "scgpt_tok", "AS", c(-1.003, -0.664, -0.052, -0.573)), list("scprint_R0", "scgpt_hid", "AS", c(-0.045, -0.120, 0.201, 0.012)))
  k3 <- c("DoRothEA_AB", "CollecTRI", "KD_response")
  for (s in k7) { m <- cells(L, s[[1]], s[[2]], k3, s[[3]], s[[3]], "K562"); r <- s[[4]]
    for (j in 1:3) chk(sprintf("%s minus %s, %s", s[[1]], s[[2]], k3[j]), mean(m$lr[m$truth == k3[j]]), r[j], 3)
    chk(sprintf("%s minus %s, mean of the three truths (d_i)", s[[1]], s[[2]]), dvec(m, "K562")[["K562"]], r[4], 3) }
  say(sprintf("== reproduction: %s", if (fails == 0L) "ALL MATCH to printed precision" else sprintf("%d MISMATCHES", fails))); fails
}

if (sys.nframe() == 0L) {
  args0 <- commandArgs(trailingOnly = TRUE)
  o <- list(`in` = "results/scgpt_ext", datasets = "hESC,hHep,K562,HepG2", out = "results/scgpt_ext/stats", legacy = 0L, truths = "DoRothEA_AB,CollecTRI,Complexes", delta = 0.15)
  i <- 1L; while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(o)); o[[k]] <- if (k %in% c("legacy")) as.integer(args0[i + 1L]) else if (k == "delta") as.numeric(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
  DELTA <- o$delta; dir.create(o$out, recursive = TRUE, showWarnings = FALSE); logf <- file.path(o$out, "stats_ext.log"); if (file.exists(logf)) file.remove(logf)
  say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
  ds <- strsplit(o$datasets, ",")[[1]]; truths <- strsplit(o$truths, ",")[[1]]
  if (o$legacy == 1L) {
    L <- load_legacy(); fails <- legacy_check(L, say)
    say("== the same statistics machinery on the old numbers (hESC, hHep, K562; truths as scored then; n = 3: no interval)")
    run_stats(L, c("hESC", "hHep", "K562"), c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF", "DoRothEA_AB", "CollecTRI", "KD_response"), out = o$out, say = say)
    if (fails > 0L) quit(status = 1L)
  } else {
    L <- load_ext(o$`in`, ds)
    edge <- do.call(rbind, lapply(ds, function(d) { f <- file.path(o$`in`, d, "tcaclr2_edge_t2.csv"); if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE)[, c("dataset", "unit", "G", "N", "edge", "t2")] }))
    adepth <- do.call(rbind, lapply(ds, function(d) { f <- file.path(o$`in`, d, "aclr_depth.csv"); if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) }))
    run_stats(L, ds, truths, out = o$out, edge = edge, adepth = adepth, say = say)
  }
  say("done")
}
