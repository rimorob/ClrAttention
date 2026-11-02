#!/usr/bin/env Rscript
# Stage 2 of the preregistered scoring (docs/preregistration_scfm.md, change log E1-E13): from the long metric table of
# analysis/scfm_score.R to the preregistered primary test, key secondaries, robust NI, and the descriptive tables. It reads no truth
# file; it only reads <dir>/metrics_long.csv and <dir>/cases_meta.csv.
#
# Definitions (all from the preregistration):
#  d_i  = for dataset i, the mean over gene sets and primary truths (celltype_ChIP, nonspecific_ChIP, STRING, LOFGOF [mESC only]) of
#         log2( metric(reference) / metric(candidate) ); higher values favour the reference (CLR-NS in the primary).
#         Only cells (dataset x gene set x truth) present for both arms are used. EPR is floored at 1e-9 (as beeline_baselines.R).
#  Interval: exact Wilcoxon / Hodges-Lehmann interval [L, U] = the C-th and (M+1-C)-th ordered Walsh averages, M = n(n+1)/2, C the largest
#         integer with P(W <= C-1) <= 0.05 (n = 7: C = 4, the 4th and 25th of 28; coverage 1 - 2 * 5/128 = 92.2%).
#  Margin delta = 0.15 (log2 units). Primary p = exact one-sided signed-rank p of H0: median <= -delta, on d_i + delta.
#  Decision table: L > -delta NI; L > 0 superior; -delta < L and U < delta tie; U < 0 learned wins; U < -delta learned wins materially;
#         otherwise inconclusive.
#  Robust NI: scPRINT takes, per dataset, the best of {R0, R0', R1, R2} x {S = t(A), A, max(A, t(A))} x {arm A, arm B where available}
#         (chosen on the truths = the configuration with the smallest d_i); NI is robust iff L > -delta on those d_i (E11).
# Usage:  Rscript analysis/scfm_score_stats.R --dir results/scfm_score [--tag ""]
args0 <- commandArgs(trailingOnly = TRUE)
o <- list(dir = "results/scfm_score", delta = 0.15)
i <- 1L; while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(o))
  o[[k]] <- if (k == "delta") as.numeric(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
DELTA <- o$delta
M <- utils::read.csv(file.path(o$dir, "metrics_long.csv"), stringsAsFactors = FALSE)
CM <- utils::read.csv(file.path(o$dir, "cases_meta.csv"), stringsAsFactors = FALSE)
PRIMARY <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
DS <- c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")
logf <- file.path(o$dir, "stats.log"); if (file.exists(logf)) file.remove(logf)
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

## ---- exact signed-rank machinery ------------------------------------------------------------------------------------------
walsh <- function(d) { n <- length(d); w <- outer(d, d, "+") / 2; sort(w[upper.tri(w, diag = TRUE)]) }
c_crit <- function(n) { C <- 0L; while (C + 1L <= n * (n + 1) / 2 && stats::psignrank(C, n) <= 0.05) C <- C + 1L; C }  # largest C with P(W <= C-1) <= 0.05
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
  fl <- c(NI = L > -delta, superior = L > 0, tie = (L > -delta) && (U < delta), learned_wins = U < 0, learned_wins_materially = U < -delta)
  fl <- fl[!is.na(fl)]
  paste(c(names(fl)[fl], if (!any(fl[c("NI", "superior", "tie", "learned_wins")])) "inconclusive"), collapse = " + ")
}

## ---- cell-level contrasts -----------------------------------------------------------------------------------------------
key3 <- function(x) paste(x$dataset, x$N, x$truth, x$gset)
pick <- function(src, rd, ori = NULL, gset = NULL, cases = NULL) {
  x <- M[M$source == src & M$readout == rd, ]
  if (!is.null(ori)) x <- x[x$orient == ori, ]
  if (!is.null(gset)) x <- x[x$gset == gset, ]
  x
}
# cells: join reference and candidate rows on dataset, N, truth and gene set
cells <- function(ref, cand, metric = "auprc_ratio", truths = PRIMARY) {
  r <- ref[ref$truth %in% truths, ]; cd_ <- cand[cand$truth %in% truths, ]
  m <- merge(r[, c("dataset", "N", "truth", "gset", metric)], cd_[, c("dataset", "N", "truth", "gset", metric)],
             by = c("dataset", "N", "truth", "gset"), suffixes = c("_ref", "_cand"))
  f <- if (metric == "epr") function(v) pmax(v, 1e-9) else identity
  m$lr <- log2(f(m[[paste0(metric, "_ref")]]) / f(m[[paste0(metric, "_cand")]]))
  m
}
dvec <- function(m) { v <- tapply(m$lr, m$dataset, mean); v[order(match(names(v), DS))] }
# scGPT rows live on two gene sets per case; select the scGPT-own set (key_S) or the three-way set (key_AS)
restrict_gset <- function(x, which) {
  if (!nrow(x)) return(x)
  k <- CM[, c("dataset", "N", paste0("key_", which))]; names(k)[3] <- "gset"
  merge(x, k, by = c("dataset", "N", "gset"))
}
ref_cmp <- function(rd, gset = NULL) pick("cmp", rd, gset = gset)

summ_row <- function(label, d, kind = "descriptive", extra = list()) {
  h <- hl(d)
  data.frame(label = label, kind = kind, n = h$n, HL = h$hl, L = h$L, U = h$U, coverage = h$coverage, dmin = h$min, dmax = h$max,
             p_NI = if (h$n == 7L) p_ni(d) else NA_real_,
             decision = if (h$n == 7L && kind %in% c("primary", "key secondary", "robust")) decide(h$L, h$U) else NA_character_,   # decision rule only where preregistered
             values = paste(sprintf("%s=%+.3f", names(d), d), collapse = "; "), stringsAsFactors = FALSE)
}
out_rows <- list()
add_row <- function(...) { r <- summ_row(...); out_rows[[length(out_rows) + 1L]] <<- r; invisible(r) }

## ---- 1. PRIMARY: CLR-NS vs scPRINT medium-v1.5 R0 (arm A, S = t(A)), AUPRC ratio -------------------------------------------
say("== PRIMARY: CLR-NS vs scPRINT medium-v1.5 R0, arm A, S = t(A); d_i = mean over gene sets and primary truths of log2(AUPRC ratio CLR-NS / AUPRC ratio scPRINT)")
cand0 <- pick("scprint_medium_A", "R0", "T", "A")
prim_cells <- cells(ref_cmp("clr_ns", "A"), cand0)
d0 <- dvec(prim_cells)
stopifnot(length(d0) == 7L)
h0 <- hl(d0); pni <- p_ni(d0); dec0 <- decide(h0$L, h0$U)
say(sprintf("  cells used: %d (14 gene-set x dataset combinations, primary truths)", nrow(prim_cells)))
say("  d_i: ", paste(sprintf("%s %+.4f", names(d0), d0), collapse = ", "))
say(sprintf("  HL estimate %+.4f; interval [%+.4f, %+.4f] (coverage %.4f, C = %d of M = %d); margin delta = %.2f", h0$hl, h0$L, h0$U, h0$coverage, h0$C, 28L, DELTA))
say(sprintf("  primary one-sided exact p (H0: median <= -delta) = %.5f   (reject at 0.05 iff W+ >= 25; smallest attainable %.5f)", pni, 1 / 128))
say("  decision-rule rows that apply: ", dec0)
say(sprintf("  NI (L > -delta): %s; superiority (L > 0), tested next in the fixed sequence%s: %s", h0$L > -DELTA,
            if (h0$L > -DELTA) "" else " (NOT reached: NI failed)", h0$L > 0))
write.csv(prim_cells, file.path(o$dir, "primary_cells.csv"), row.names = FALSE)
add_row("PRIMARY CLR-NS vs R0 (arm A, S=t(A), AUPRC)", d0, "primary")
primary <- list(d = d0, hl = h0, p = pni, decision = dec0)

## ---- 2. Key secondaries (Holm over two tests, alpha 0.05, each n = 7) --------------------------------------------------------
say("== KEY SECONDARIES (confirmatory only if primary NI holds; Holm at 0.05 over S1, S2)")
S1_cells <- cells(ref_cmp("clr", "A"), cand0); d_s1 <- dvec(S1_cells)
S2_cells <- cells(ref_cmp("clr_ns", "A"), cand0, metric = "epr"); d_s2 <- dvec(S2_cells)
h1 <- hl(d_s1); h2 <- hl(d_s2); p1 <- p_ni(d_s1); p2 <- p_ni(d_s2)
ps <- c(S1 = p1, S2 = p2); ord <- order(ps); holm_rej <- setNames(logical(2), names(ps))
thr <- c(0.025, 0.05); for (k in seq_along(ord)) { if (ps[ord[k]] <= thr[k]) holm_rej[ord[k]] <- TRUE else break }
for (nm in c("S1", "S2")) {
  d <- if (nm == "S1") d_s1 else d_s2; h <- if (nm == "S1") h1 else h2
  say(sprintf("  %s (%s): d_i %s", nm, if (nm == "S1") "CLR vs R0, AUPRC ratio" else "CLR-NS vs R0, EPR",
              paste(sprintf("%s %+.3f", names(d), d), collapse = ", ")))
  say(sprintf("     HL %+.4f, interval [%+.4f, %+.4f], p_NI = %.5f, Holm-rejected: %s, rows: %s", h$hl, h$L, h$U, ps[[nm]], holm_rej[[nm]], decide(h$L, h$U)))
  add_row(paste0(nm, if (nm == "S1") " CLR vs R0 (AUPRC)" else " CLR-NS vs R0 (EPR)"), d, "key secondary")
}
say(sprintf("  confirmatory status: primary NI %s -> the secondaries are %s", if (h0$L > -DELTA) "holds" else "does NOT hold",
            if (h0$L > -DELTA) "interpreted confirmatorily" else "descriptive only"))

## ---- 3. Robust NI -----------------------------------------------------------------------------------------------------------
say("== ROBUST NI: scPRINT takes, per dataset, the best of {R0,R0',R1,R2} x {T,A,max} x {arm A, arm B where available}")
cfg <- rbind(
  expand.grid(src = "scprint_medium_A", rd = c("R0", "R0p", "R2"), ori = c("T", "A", "max"), stringsAsFactors = FALSE),
  expand.grid(src = "r1_A", rd = "R1", ori = c("T", "A", "max"), stringsAsFactors = FALSE),
  expand.grid(src = "scprint_medium_B", rd = c("R0", "R0p", "R2"), ori = c("T", "A", "max"), stringsAsFactors = FALSE),
  expand.grid(src = "r1_B", rd = "R1", ori = c("T", "A", "max"), stringsAsFactors = FALSE))
dm <- NULL
for (k in seq_len(nrow(cfg))) {
  cd <- pick(cfg$src[k], cfg$rd[k], cfg$ori[k], gset = if (cfg$src[k] %in% c("scprint_medium_A", "r1_A")) "A" else NULL); if (!nrow(cd)) next
  d <- dvec(cells(ref_cmp("clr_ns"), cd))
  dm <- rbind(dm, data.frame(config = sprintf("%s|%s|%s", cfg$src[k], cfg$rd[k], cfg$ori[k]), dataset = names(d), d = as.numeric(d)))
}
write.csv(dm, file.path(o$dir, "robust_configs_d.csv"), row.names = FALSE)
best <- do.call(rbind, lapply(DS, function(ds) { x <- dm[dm$dataset == ds, ]; x[which.min(x$d), ] }))
d_rob <- setNames(best$d, best$dataset)
hr <- hl(d_rob)
say("  best-for-scPRINT configuration per dataset (smallest d_i):"); for (k in seq_len(nrow(best))) say(sprintf("    %-8s %-40s d = %+.4f", best$dataset[k], best$config[k], best$d[k]))
say(sprintf("  robust: HL %+.4f; interval [%+.4f, %+.4f]; p_NI = %.5f; rows: %s", hr$hl, hr$L, hr$U, p_ni(d_rob), decide(hr$L, hr$U)))
say(sprintf("  ROBUST NI (L > -delta): %s", hr$L > -DELTA))
add_row("ROBUST NI: best of {R0,R0',R1,R2}x3 orient x arms", d_rob, "robust"); write.csv(best, file.path(o$dir, "robust_best_config.csv"), row.names = FALSE)
# arm-A-only version of the maximum (the preregistered maximum before E11 extended it with R0' and arm B), for reference
dmA <- dm[grepl("_A\\|", dm$config), ]; bestA <- do.call(rbind, lapply(DS, function(ds) { x <- dmA[dmA$dataset == ds, ]; x[which.min(x$d), ] }))
add_row("robust max, arm A only, {R0,R0',R1,R2}x3", setNames(bestA$d, bestA$dataset), "robust (reference)")
# the original preregistered maximum {R0,R1,R2} x 3 orientations (arm A)
dmP <- dmA[!grepl("R0p", dmA$config), ]; bestP <- do.call(rbind, lapply(DS, function(ds) { x <- dmP[dmP$dataset == ds, ]; x[which.min(x$d), ] }))
add_row("robust max, original prereg set {R0,R1,R2}x3, arm A", setNames(bestP$d, bestP$dataset), "robust (reference)")

## ---- 4. DESCRIPTIVE: other readouts / orientations / arms of scPRINT vs CLR-NS (AUPRC) ----------------------------------------
say("== DESCRIPTIVE (no claims): CLR-NS vs each scPRINT readout / orientation / arm (AUPRC), n = 7 unless stated")
for (k in seq_len(nrow(cfg))) {
  cd <- pick(cfg$src[k], cfg$rd[k], cfg$ori[k], gset = if (cfg$src[k] %in% c("scprint_medium_A", "r1_A")) "A" else NULL); if (!nrow(cd)) next
  d <- dvec(cells(ref_cmp("clr_ns"), cd))
  add_row(sprintf("CLR-NS vs %s %s S=%s", cfg$src[k], cfg$rd[k], cfg$ori[k]), d, if (grepl("_B$", cfg$src[k])) "descriptive arm B" else "descriptive")
}
# arm B: per-dataset d_i under B, paired difference of scPRINT R0 (B minus A) in log2 AUPRC ratio; no p-value
say("  Arm B (n = 4: hESC and mHSC-E/GM/L): paired delta_i = log2 AUPRC ratio of scPRINT R0 under B minus under A")
cB <- pick("scprint_medium_B", "R0", "T"); cA <- pick("scprint_medium_A", "R0", "T", gset = "A")
mm <- merge(cA[cA$truth %in% PRIMARY, c("dataset", "N", "truth", "auprc_ratio")], cB[cB$truth %in% PRIMARY, c("dataset", "N", "truth", "auprc_ratio")],
            by = c("dataset", "N", "truth"), suffixes = c("_A", "_B"))
mm$delta <- log2(mm$auprc_ratio_B) - log2(mm$auprc_ratio_A)
dl <- tapply(mm$delta, mm$dataset, mean); dl <- dl[order(match(names(dl), DS))]
hb <- hl(dl)
say("    delta_i: ", paste(sprintf("%s %+.4f", names(dl), dl), collapse = ", "))
say(sprintf("    HL %+.4f, range [%+.4f, %+.4f], median %+.4f; reading rule: median delta_i %s delta = %.2f -> %s", hb$hl, hb$min, hb$max,
            stats::median(dl), if (stats::median(dl) > DELTA) ">" else "<=", DELTA,
            if (stats::median(dl) > DELTA) "the primary result depends on the input scale on the count datasets" else "no input-scale dependence beyond the margin on the count datasets"))
add_row("ARM B: delta_i = R0(B) - R0(A) log2 AUPRC ratio (n=4)", dl, "descriptive arm B")
write.csv(mm, file.path(o$dir, "armB_paired_cells.csv"), row.names = FALSE)
# secondary scPRINT checkpoints
for (s in c("scprint_large_A", "scprint_small_A")) add_row(sprintf("CLR-NS vs %s R0 S=t(A)", s), dvec(cells(ref_cmp("clr_ns", "A"), pick(s, "R0", "T", "A"))), "descriptive checkpoint")
# controls
say("  Controls (scPRINT, R0, S = t(A)): d_i vs CLR-NS, and scPRINT R0 minus control in log2 AUPRC ratio")
for (grp in list(perm = paste0("ctrl_perm_s", 1001:1003), rand_esm = paste0("ctrl_rand_esm_s", 1:3), rand_emb = paste0("ctrl_rand_emb_s", 1:3))) {
  nm <- sub("_s[0-9]+$", "", grp[1])
  dd <- sapply(grp, function(s) dvec(cells(ref_cmp("clr_ns", "A"), pick(s, "R0", "T", "A"))))
  add_row(sprintf("CLR-NS vs control %s (mean over 3 seeds)", nm), rowMeans(dd), "descriptive control")
  for (s in grp) add_row(sprintf("CLR-NS vs control %s", s), dd[, s], "descriptive control seed")
  da <- sapply(grp, function(s) dvec(cells(pick("scprint_medium_A", "R0", "T", "A"), pick(s, "R0", "T", "A"))))
  add_row(sprintf("scPRINT R0 minus control %s (mean over 3 seeds; log2 AUPRC ratio)", nm), rowMeans(da), "descriptive control")
}
add_row("scPRINT R0 minus ESM2 cosine (data-free prior; log2 AUPRC ratio)", dvec(cells(pick("scprint_medium_A", "R0", "T", "A"), pick("esm2_cos", "cos", gset = "A"))), "descriptive control")
add_row("CLR-NS vs ESM2 cosine", dvec(cells(ref_cmp("clr_ns", "A"), pick("esm2_cos", "cos", gset = "A"))), "descriptive control")
add_row("CLR-NS vs CLR applied to A (hybrid; symmetrised R0)", dvec(cells(ref_cmp("clr_ns", "A"), pick("scprint_medium_A", "R0_clr_hybrid", gset = "A"))), "descriptive control")
add_row("CLR-NS vs |r| (G-intersection)", dvec(cells(ref_cmp("clr_ns", "A"), ref_cmp("abs_pearson", "A"))), "descriptive")
add_row("CLR vs |r| (G-intersection)", dvec(cells(ref_cmp("clr", "A"), ref_cmp("abs_pearson", "A"))), "descriptive")
add_row("CLR-NS vs CLR (G-intersection)", dvec(cells(ref_cmp("clr_ns", "A"), ref_cmp("clr", "A"))), "descriptive")
# CLR-NS soft10 attention at its held-out depth vs scPRINT R0 (preregistered item 5)
if (nrow(ref_cmp("clr_ns_att", "A"))) {
  add_row("CLR-NS soft10 attention (held-out depth) vs scPRINT R0 [item 5]", dvec(cells(ref_cmp("clr_ns_att", "A"), cand0)), "descriptive")
  add_row("CLR-NS soft10 attention vs CLR-NS", dvec(cells(ref_cmp("clr_ns_att", "A"), ref_cmp("clr_ns", "A"))), "descriptive")
}
# GENIE3 / GRNBoost2 / PIDC (fitted on the full gene set, scored on G-A) vs scPRINT R0
for (b in c("genie3", "grnboost2", "pidc")) add_row(sprintf("%s (full-set fit, scored on G-A) vs scPRINT R0", b), dvec(cells(ref_cmp(b, "A"), cand0)), "descriptive baseline")
for (b in c("genie3", "grnboost2", "pidc")) add_row(sprintf("CLR-NS vs %s (full-set fit, scored on G-A)", b), dvec(cells(ref_cmp("clr_ns", "A"), ref_cmp(b, "A"))), "descriptive baseline")
# derived truths (secondary truths): the primary contrast on the primary + derived truths, and on the derived truths alone
ALLT <- c(PRIMARY, "STRING_regulatory", "STRING_curated")
add_row("CLR-NS vs R0, primary + derived STRING truths", dvec(cells(ref_cmp("clr_ns", "A"), cand0, truths = ALLT)), "secondary truths")
add_row("CLR-NS vs R0, derived STRING truths only", dvec(cells(ref_cmp("clr_ns", "A"), cand0, truths = c("STRING_regulatory", "STRING_curated"))), "secondary truths")
for (tn in c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")) {
  v <- tryCatch(dvec(cells(ref_cmp("clr_ns", "A"), cand0, truths = tn)), error = function(e) NULL)
  if (!is.null(v) && length(v)) add_row(sprintf("CLR-NS vs R0, truth %s only", tn), v, "per-truth")
}
for (nn in c(500, 1000)) {
  cc <- prim_cells[prim_cells$N == nn, ]; add_row(sprintf("CLR-NS vs R0, TFs+%d only", nn), dvec(cc), "per-gene-set")
}

## ---- 5. scGPT (human: hESC, hHep) and scPRINT vs scGPT (descriptive; no p-value, no NI claim) -------------------------------
say("== scGPT (descriptive; n = 2 BEELINE datasets plus K562 handled separately; no p-value)")
for (rd in c("attn_last_mean", "tok_cos", "hid_cos", "tok_cos_ln")) {
  ori <- if (rd == "attn_last_mean") "T" else "sym"
  cs <- restrict_gset(pick("scgpt_A", rd, ori), "S"); if (!nrow(cs)) next
  rf <- restrict_gset(ref_cmp("clr_ns"), "S")
  add_row(sprintf("CLR-NS vs scGPT %s (G-scGPT)%s", rd, if (rd == "tok_cos_ln") " [exploratory readout]" else ""), dvec(cells(rf, cs)), "descriptive scGPT")
  if (rd != "tok_cos_ln") {
    rp <- restrict_gset(pick("scprint_medium_A", "R0", "T"), "AS"); cs3 <- restrict_gset(pick("scgpt_A", rd, ori), "AS")
    add_row(sprintf("scPRINT R0 vs scGPT %s (three-way intersection)", rd), dvec(cells(rp, cs3)), "descriptive scPRINT vs scGPT")
  }
}

## ---- write the summary table ------------------------------------------------------------------------------------------------
SUM <- do.call(rbind, out_rows)
write.csv(SUM, file.path(o$dir, "contrasts_summary.csv"), row.names = FALSE)
for (k in seq_len(nrow(SUM))) say(sprintf("  [%s] %-80s n=%d HL %+.3f%s%s", SUM$kind[k], SUM$label[k], SUM$n[k], SUM$HL[k],
  if (is.finite(SUM$L[k])) sprintf(" [%+.3f, %+.3f]", SUM$L[k], SUM$U[k]) else sprintf(" range [%+.3f, %+.3f]", SUM$dmin[k], SUM$dmax[k]),
  if (!is.na(SUM$decision[k])) paste0("  -> ", SUM$decision[k]) else ""))

## ---- per-case tables: key metrics, dropped genes and truth edges ------------------------------------------------------------
key <- rbind(
  cbind(method = "CLR-NS (G-A)", ref_cmp("clr_ns", "A")[, c("dataset", "N", "truth", "auprc_ratio", "epr")]),
  cbind(method = "CLR (G-A)", ref_cmp("clr", "A")[, c("dataset", "N", "truth", "auprc_ratio", "epr")]),
  cbind(method = "|r| (G-A)", ref_cmp("abs_pearson", "A")[, c("dataset", "N", "truth", "auprc_ratio", "epr")]),
  cbind(method = "scPRINT R0 S=t(A) arm A", cand0[, c("dataset", "N", "truth", "auprc_ratio", "epr")]))
write.csv(key, file.path(o$dir, "key_metrics_per_case.csv"), row.names = FALSE)
full <- ref_cmp("clr_ns", "full")[, c("dataset", "N", "truth", "n_edges", "n_regulators", "auprc_ratio", "epr")]
names(full)[4:7] <- paste0(c("edges", "regs", "auprc_ratio", "epr"), "_full")
ga <- ref_cmp("clr_ns", "A")[, c("dataset", "N", "truth", "n_edges", "n_regulators")]; names(ga)[4:5] <- c("edges_A", "regs_A")
drop <- merge(merge(CM, full, by = c("dataset", "N")), ga, by = c("dataset", "N", "truth"), all.x = TRUE)
drop$edges_dropped <- drop$edges_full - drop$edges_A
write.csv(drop, file.path(o$dir, "dropped_genes_and_truth_edges.csv"), row.names = FALSE)
say("== dropped truth edges (full case gene set -> scPRINT-covered set), primary truths:")
dp <- drop[drop$truth %in% PRIMARY, ]
say(sprintf("   truth edges in the full gene sets: %d; retained on G-A: %d; dropped: %d (%.1f%%)", sum(dp$edges_full), sum(dp$edges_A, na.rm = TRUE),
            sum(dp$edges_dropped, na.rm = TRUE), 100 * sum(dp$edges_dropped, na.rm = TRUE) / sum(dp$edges_full)))

## ---- sensitivity: mixed model (nlme in place of lme4, which is not installed) ----------------------------------------------
say("== SENSITIVITY (descriptive): mixed model, log2 AUPRC ratio ~ method + N + (1 | dataset) + (1 | dataset:truth); method = CLR-NS vs scPRINT R0")
mx <- NULL
try({
  library(nlme)
  a <- prim_cells; long <- rbind(data.frame(dataset = a$dataset, N = factor(a$N), truth = a$truth, method = "CLR-NS", y = log2(a$auprc_ratio_ref)),
                                 data.frame(dataset = a$dataset, N = factor(a$N), truth = a$truth, method = "scPRINT", y = log2(a$auprc_ratio_cand)))
  long$method <- factor(long$method, levels = c("scPRINT", "CLR-NS")); long$dataset <- factor(long$dataset); long$truth <- factor(long$truth)
  fit <- nlme::lme(y ~ method + N, random = ~ 1 | dataset/truth, data = long, method = "REML")
  ci <- nlme::intervals(fit, which = "fixed")$fixed
  say(sprintf("   method effect (CLR-NS minus scPRINT) = %+.4f, 95%% CI [%+.4f, %+.4f] (log2 units; nlme::lme with random intercepts dataset and truth within dataset)",
              ci["methodCLR-NS", "est."], ci["methodCLR-NS", "lower"], ci["methodCLR-NS", "upper"]))
  mx <- as.data.frame(ci); write.csv(mx, file.path(o$dir, "mixed_model_fixed.csv"))
}, silent = FALSE)

say("== leaked-dataset sensitivity: empty (E9: no scPRINT v1 or scGPT pair is classified 'in')")
saveRDS(list(primary = primary, secondaries = list(S1 = list(d = d_s1, hl = h1, p = p1), S2 = list(d = d_s2, hl = h2, p = p2), holm = holm_rej),
             robust = list(d = d_rob, hl = hr, p = p_ni(d_rob)), armB = list(delta = dl, hl = hb)), file.path(o$dir, "stats_results.rds"))
say("done")
