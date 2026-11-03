#!/usr/bin/env Rscript
# Descriptive fine-tuning arm (E10/E17/E18). Stage 2 of 2: from the metric tables to the descriptive contrasts. Reads no truth file.
#   fine-tuned R0 vs zero-shot R0:      d_FT-ZS = per dataset, mean over gene sets and primary truths of log2(AUPRC ratio fine-tuned / zero-shot)
#                                       (positive = fine-tuning raised the score)
#   CLR-NS minus fine-tuned R0:         the SAME d_i definition as the primary (analysis/scfm_score_stats.R): mean over gene sets and primary
#                                       truths of log2(AUPRC ratio CLR-NS / AUPRC ratio candidate); higher favours CLR-NS
#   reference rows from the E16 table:  CLR-NS minus zero-shot R0 (= the primary d_i; HL +1.032, interval [+0.343, +1.897])
# The exact Wilcoxon / Hodges-Lehmann machinery, the cell join (cells()) and dvec() are EVALUATED from analysis/scfm_score_stats.R. Descriptive only:
# no p-value, no non-inferiority claim, no decision-table label (E10). n = 7 datasets (scPRINT), n = 2 (scGPT: per-dataset values, no interval).
# Usage: Rscript analysis/scfm_score_ft_stats.R --zs results/scfm_score --ft results/scfm_score_ft [--out results/scfm_score_ft]
args0 <- commandArgs(trailingOnly = TRUE)
o <- list(zs = "results/scfm_score", ft = "results/scfm_score_ft", out = NA_character_, delta = 0.15)
i <- 1L; while (i <= length(args0)) { k <- sub("^--", "", args0[i]); stopifnot(k %in% names(o))
  o[[k]] <- if (k == "delta") as.numeric(args0[i + 1L]) else args0[i + 1L]; i <- i + 2L }
if (is.na(o$out)) o$out <- o$ft
DELTA <- o$delta; PRIMARY <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF"); DS <- c("hESC", "hHep", "mDC", "mESC", "mHSC-E", "mHSC-GM", "mHSC-L")
Z <- utils::read.csv(file.path(o$zs, "metrics_long.csv"), stringsAsFactors = FALSE)
Fz <- utils::read.csv(file.path(o$ft, "metrics_ft.csv"), stringsAsFactors = FALSE)
CM <- utils::read.csv(file.path(o$zs, "cases_meta.csv"), stringsAsFactors = FALSE)
common <- intersect(names(Z), names(Fz))
M <- rbind(Z[, common], Fz[, common])
logf <- file.path(o$out, "stats_ft.log"); if (file.exists(logf)) file.remove(logf)
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
# machinery evaluated from the primary statistics script
st <- parse("analysis/scfm_score_stats.R")
for (e in st) if (is.call(e) && identical(as.character(e[[1]]), "<-") && as.character(e[[2]])[1] %in% c("walsh", "c_crit", "hl", "p_ni", "cells", "dvec", "pick", "restrict_gset", "ref_cmp")) eval(e, globalenv())
stopifnot(all(vapply(c("walsh", "c_crit", "hl", "cells", "dvec", "pick", "restrict_gset", "ref_cmp"), exists, TRUE)))
out_rows <- list()
add_row <- function(label, d, kind = "descriptive", metric = "AUPRC ratio") {
  d <- d[is.finite(d)]; n <- length(d)
  if (n >= 5L) h <- hl(d) else h <- list(hl = stats::median(d), L = NA_real_, U = NA_real_, coverage = NA_real_, min = min(d), max = max(d))
  out_rows[[length(out_rows) + 1L]] <<- data.frame(label = label, kind = kind, metric = metric, n = n, HL = h$hl, L = h$L, U = h$U, coverage = h$coverage,
    dmin = h$min, dmax = h$max, n_positive = sum(d > 0), values = paste(sprintf("%s=%+.3f", names(d), d), collapse = "; "), stringsAsFactors = FALSE)
  r <- out_rows[[length(out_rows)]]
  say(sprintf("  [%s] %-72s n=%d HL %+.3f%s  (%d of %d positive)  %s", kind, label, n, r$HL, if (is.finite(r$L)) sprintf(" [%+.3f, %+.3f]", r$L, r$U) else sprintf(" range [%+.3f, %+.3f]", r$dmin, r$dmax),
              r$n_positive, n, r$values))
  invisible(r)
}
say("== FINE-TUNING ARM (descriptive; no p-value, no NI claim). d_i definitions: see the header of analysis/scfm_score_ft_stats.R")
cmp <- function(rd, gset = "A") ref_cmp(rd, gset)
FTR0 <- pick("scprint_ft_A", "R0", "T", "A"); ZSR0 <- pick("scprint_medium_A", "R0", "T", "A")
stopifnot(nrow(FTR0) > 0, nrow(ZSR0) > 0)
## reference: the primary contrast reproduced from the E16 table
prim <- dvec(cells(cmp("clr_ns"), ZSR0))
add_row("REFERENCE (E16) CLR-NS minus zero-shot R0, S=t(A), AUPRC ratio", prim, "reference")
## A. the two preregistered-descriptive contrasts
dft <- dvec(cells(FTR0, ZSR0))
add_row("A1  fine-tuned R0 minus zero-shot R0, S=t(A), AUPRC ratio (log2)", dft, "fine-tuning effect")
dcf <- dvec(cells(cmp("clr_ns"), FTR0))
add_row("A2  CLR-NS minus fine-tuned R0, S=t(A), AUPRC ratio (same d_i as the primary)", dcf, "CLR-NS vs fine-tuned")
add_row("A2b CLR-NS minus fine-tuned R0, S=t(A), EPR (log2, floor 1e-9)", dvec(cells(cmp("clr_ns"), FTR0, metric = "epr")), "CLR-NS vs fine-tuned", "EPR")
add_row("A1b fine-tuned R0 minus zero-shot R0, S=t(A), EPR (log2, floor 1e-9)", dvec(cells(FTR0, ZSR0, metric = "epr")), "fine-tuning effect", "EPR")
add_row("A3  CLR minus fine-tuned R0, S=t(A), AUPRC ratio", dvec(cells(cmp("clr"), FTR0)), "CLR-NS vs fine-tuned")
## B. per truth
for (tn in PRIMARY) {
  v <- tryCatch(dvec(cells(FTR0, ZSR0, truths = tn)), error = function(e) NULL); if (!is.null(v) && length(v)) add_row(sprintf("B1  fine-tuned minus zero-shot R0, truth %s only", tn), v, "fine-tuning effect per truth")
  v <- tryCatch(dvec(cells(cmp("clr_ns"), FTR0, truths = tn)), error = function(e) NULL); if (!is.null(v) && length(v)) add_row(sprintf("B2  CLR-NS minus fine-tuned R0, truth %s only", tn), v, "CLR-NS vs fine-tuned per truth")
}
for (nn in c(500, 1000)) {
  cc <- cells(FTR0, ZSR0); cc <- cc[cc$N == nn, ]; add_row(sprintf("B3  fine-tuned minus zero-shot R0, TFs+%d only", nn), dvec(cc), "fine-tuning effect per gene set")
  cc <- cells(cmp("clr_ns"), FTR0); cc <- cc[cc$N == nn, ]; add_row(sprintf("B4  CLR-NS minus fine-tuned R0, TFs+%d only", nn), dvec(cc), "CLR-NS vs fine-tuned per gene set")
}
## C. other readouts and orientations
for (cfgi in list(c("R0", "A"), c("R0", "max"), c("R0p", "T"), c("R2", "T"))) {
  rd <- cfgi[1]; ori <- cfgi[2]
  f <- pick("scprint_ft_A", rd, ori, "A"); zz <- pick("scprint_medium_A", rd, ori, "A")
  add_row(sprintf("C1  fine-tuned minus zero-shot %s S=%s", rd, ori), dvec(cells(f, zz)), "fine-tuning effect other readout")
  add_row(sprintf("C2  CLR-NS minus fine-tuned %s S=%s", rd, ori), dvec(cells(cmp("clr_ns"), f)), "CLR-NS vs fine-tuned other readout")
}
## D. scGPT (hESC, hHep; n = 2; per-dataset values, no interval)
say("== scGPT fine-tuned (n = 2 datasets: per-dataset values only)")
for (rd in c("attn_last_mean", "tok_cos", "hid_cos", "tok_cos_ln")) {
  ori <- if (rd == "attn_last_mean") "T" else "sym"
  f <- restrict_gset(pick("scgpt_ft_A", rd, ori), "S"); zz <- restrict_gset(pick("scgpt_A", rd, ori), "S"); if (!nrow(f) || !nrow(zz)) next
  rf <- restrict_gset(ref_cmp("clr_ns"), "S")
  add_row(sprintf("D1  scGPT fine-tuned minus zero-shot %s (G-scGPT)", rd), dvec(cells(f, zz)), "scGPT fine-tuning effect")
  add_row(sprintf("D2  CLR-NS minus scGPT fine-tuned %s (G-scGPT)", rd), dvec(cells(rf, f)), "CLR-NS vs scGPT fine-tuned")
  add_row(sprintf("D0  (E16) CLR-NS minus scGPT zero-shot %s (G-scGPT)", rd), dvec(cells(rf, zz)), "reference")
}
SUM <- do.call(rbind, out_rows)
write.csv(SUM, file.path(o$out, "contrasts_ft_summary.csv"), row.names = FALSE)
write.csv(cells(FTR0, ZSR0), file.path(o$out, "cells_ft_vs_zs_R0.csv"), row.names = FALSE)
write.csv(cells(cmp("clr_ns"), FTR0), file.path(o$out, "cells_clrns_vs_ft_R0.csv"), row.names = FALSE)
say("done")
