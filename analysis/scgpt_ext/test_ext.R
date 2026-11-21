#!/usr/bin/env Rscript
# Unit and smoke tests of the scGPT extension scripts. Reads no real truth score: the statistics tests use fabricated d_i, the complex-AUPR test uses the REAL complex structure
# (pair universe built on a gene list, structure only, no method score) with an oracle-plus-noise and a random score, never a method.
#   Rscript analysis/scgpt_ext/test_ext.R [--cx DIR_WITH_cx_members.tsv --genes genes.txt]
args0 <- commandArgs(TRUE); cxdir <- if ("--cx" %in% args0) args0[which(args0 == "--cx") + 1L] else ""; gfile <- if ("--genes" %in% args0) args0[which(args0 == "--genes") + 1L] else ""
suppressPackageStartupMessages(source("analysis/tcaclr_confirm/tc_lib.R")); source("analysis/scgpt_ext/ext_lib.R"); source("analysis/scgpt_ext/stats_ext.R")
fails <- 0L; ok <- function(name, cond) { cat(sprintf("%-5s %s\n", if (isTRUE(cond)) "ok" else "FAIL", name)); if (!isTRUE(cond)) fails <<- fails + 1L }

## 1. exact interval = wilcox.test(conf.int, exact) ------------------------------------------------------------------------------------
set.seed(3)
for (n in c(5L, 6L, 7L, 8L)) for (rep in 1:3) {
  d <- stats::rnorm(n, 0.3, 0.4); h <- hl(d)
  w <- suppressWarnings(stats::wilcox.test(d, conf.int = TRUE, conf.level = h$coverage - 1e-9, exact = TRUE))
  ok(sprintf("n = %d: [L, U] and HL equal wilcox.test(conf.int, exact = TRUE) at coverage %.4f", n, h$coverage),
     abs(h$L - w$conf.int[1]) < 1e-9 && abs(h$U - w$conf.int[2]) < 1e-9 && abs(h$hl - unname(w$estimate)) < 1e-9)
}
ok("n = 6: C = 3, coverage 0.9063 (3rd and 19th of 21 Walsh averages)", hl(1:6)$C == 3L && abs(hl(1:6)$coverage - 0.90625) < 1e-9 && length(walsh(1:6)) == 21L)
d6 <- c(0.1, 0.5, 0.2, 0.9, 0.4, 0.3); w6 <- walsh(d6); ok("n = 6: L and U are Walsh averages 3 and 19", hl(d6)$L == w6[3] && hl(d6)$U == w6[19])
d5 <- c(0.1, 0.5, 0.2, 0.9, 0.4); ok("n = 5: C = 1, coverage 0.9375, L = min d_i, U = max Walsh average", hl(d5)$C == 1L && abs(hl(d5)$coverage - 0.9375) < 1e-9 && hl(d5)$L == min(d5) && hl(d5)$U == max(walsh(d5)))
ok("n = 4: no interval (C = 0), decision refused", hl(c(1, 2, 3, 4))$C == 0L && grepl("no decision", decide(hl(c(1, 2, 3, 4))$L, hl(c(1, 2, 3, 4))$U)))
ok("smallest attainable p: n = 6 is 1/64, n = 5 is 1/32, n = 4 is 1/16 (all d_i + delta > 0)", abs(p_ni(rep(1, 6)) - 1 / 64) < 1e-12 && abs(p_ni(rep(1, 5)) - 1 / 32) < 1e-12 && abs(p_ni(rep(1, 4)) - 1 / 16) < 1e-12)
ok("decision table rows", decide(0.2, 0.9) == "NI + superior" && decide(-0.1, 0.1) == "NI + tie" && decide(-0.5, -0.2) == "scGPT_wins + scGPT_wins_materially" && decide(-0.3, 0.3) == "inconclusive" &&
   decide(-0.5, -0.05) == "scGPT_wins" && decide(-0.1, 0.4) == "NI")

## 2. complex AUPR on the real structure with oracle and random scores ------------------------------------------------------------------
if (nzchar(cxdir) && nzchar(gfile)) {
  genes <- readLines(gfile); tr <- tc_read_truth(file.path(cxdir, "cx_members.tsv"), file.path(cxdir, "cx_excluded_pairs.tsv")); G <- length(genes)
  p <- tc_pairs(tc_sets(tr$mem, "main", genes), genes, tc_excl_keys(tr$excl, "all", genes)); cat(sprintf("complex structure: %d pairs, %d positive, base %.4f\n", length(p$label), sum(p$label), mean(p$label)))
  set.seed(5); S0 <- matrix(stats::runif(G * G), G, G); S0 <- (S0 + t(S0)) / 2; Sor <- S0 * 0.2; Sor[cbind(p$i, p$j)] <- Sor[cbind(p$i, p$j)] + p$label; Sor <- pmax(Sor, t(Sor))
  ro <- ext_score_complex(Sor, p); rr <- ext_score_complex(S0, p)
  ok(sprintf("oracle score: AUPR ratio %.1f >> 1; random score: %.2f near 1; base rate equals mean label", ro[["auprc_ratio"]], rr[["auprc_ratio"]]), ro[["auprc_ratio"]] > 10 && abs(rr[["auprc_ratio"]] - 1) < 0.35 && abs(ro[["density"]] - mean(p$label)) < 1e-12)
}

## 3. fabricated multi-dataset tables: n = 6, 5, 4 through run_stats --------------------------------------------------------------------
fab <- function(ds, truths = EXT_PRIMARY) {          # CLR-NS beats attention by ~0.3 log2, tcaCLR by ~0.2, token cosine beats both, per-dataset shift
  do.call(rbind, lapply(seq_along(ds), function(i) do.call(rbind, lapply(1:2, function(u) do.call(rbind, lapply(truths, function(t) {
    base <- 1.4; sh <- 0.05 * i; data.frame(dataset = ds[i], unit = as.character(u), truth = t, gset = "S",
      method = c("clr_ns", "tcaCLR", "tcaCLR2", "abs_pearson", "clr", "scgpt_attn", "scgpt_tok", "scgpt_hid", "scgpt_tokln", "aCLR"),
      auprc_ratio = base * c(1.25 + sh, 1.15 + sh, 1.1, 1.0, 1.05, 1.0, 1.5, 1.2, 1.5, 1.12 + sh) + 0.01 * u, epr = NA_real_, stringsAsFactors = FALSE) }))))) )
}
for (n in c(6L, 5L, 4L)) {
  ds <- c("hESC", "hHep", "K562", "HepG2", "RPE1", "Jurkat")[seq_len(n)]; L <- fab(ds)
  cat(sprintf("\n-- fabricated n = %d\n", n)); r <- run_stats(L, ds, EXT_PRIMARY, out = NULL, say = function(...) invisible(NULL))
  h <- r$holm; ok(sprintf("n = %d: P1 NI computed from d_i (HL %.3f)", n, hl(r$d1)$hl), if (n >= 5L) h$NI[1] else is.na(hl(r$d1)$L))
  ok(sprintf("n = %d: Holm confirmatory flags %s", n, if (n == 6L) "TRUE for S2 (p = 1/64 below 0.025), FALSE for S1 (token cosine better)" else "all FALSE (smallest p above 0.025)"),
     if (n == 6L) isTRUE(h$confirmatory[3]) && !isTRUE(h$confirmatory[2]) else !any(h$confirmatory[2:3], na.rm = TRUE))
  ok(sprintf("n = %d: S1 d_i is the smallest over the three readouts (never above the attention d_i)", n), all(r$d_s1 <= r$d1 + 1e-12))
  ok(sprintf("n = %d: tcaCLR2 enters no Holm row and no decision", n), !any(grepl("tcaCLR2", r$holm$test)) && all(is.na(r$S$decision[grepl("tcaCLR2", r$S$label)])))
  ka <- r$S$kind == "descriptive (aCLR)"          # (the label regex "aCLR" would also match "tcaCLR")
  ok(sprintf("n = %d: aCLR is present as descriptive rows, enters no Holm row, no p and no decision", n), any(ka) && !any(grepl("(^| )aCLR", r$holm$test)) &&
     all(is.na(r$S$decision[ka])) && all(is.na(r$S$p_NI[ka])) && identical(r$holm$test, c("P1", "S1", "S2")))
}
cat(sprintf("\n%s\n", if (fails == 0L) "ALL TESTS PASSED" else sprintf("%d TESTS FAILED", fails))); if (fails > 0L) quit(status = 1L)
