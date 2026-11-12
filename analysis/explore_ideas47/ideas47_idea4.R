#!/usr/bin/env Rscript
# EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED. Idea 4 of docs/plans/fable_ideation.md: topology-gated hybrid readout.
# Readouts, gates and scoring fixed in docs/plans/explore_ideas47_plan.md section 1 (committed before compute).
# Usage (repo root): Rscript analysis/explore_ideas47/ideas47_idea4.R [--out results/explore_ideas47/idea4] [--threads 8] [--workers 8] [--jack 100] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(out = "results/explore_ideas47/idea4", threads = 8L, workers = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "ideas47_lib.R"))

## ---- references and equivalence gate (nothing new is computed before it passes) -----------------------------------
Ta3r <- a3r(Pt9)
res_clr <- score_readout("CLR", C); res_k10 <- score_readout("K10", S9); res_a3r <- score_readout("a3r", Ta3r)
if (!QUICK) e47_gate(res_clr, res_k10, res_a3r) else say("QUICK run: gate skipped (gene subset)")

## ---- gate c_i: Onnela weighted clustering coefficient of the top-10 attention neighbourhood ----------------------
Aoff <- A; diag(Aoff) <- 0
W <- (Aoff + t(Aoff)) / 2; wmax <- max(W)
nn <- t(apply(Aoff, 1, function(v) order(v, decreasing = TRUE)[1:10]))
cvec <- vapply(seq_len(G), function(i) { n <- nn[i, ]; w <- W[n, n] / wmax; diag(w) <- 0; wi <- W[i, n] / wmax
  sum((outer(wi, wi) * w)^(1 / 3)) / 90 }, 0)
dens <- stats::density(cvec); c_mode <- dens$x[which.max(dens$y)]; c_sd <- stats::sd(cvec)
s_soft <- stats::pnorm((cvec - c_mode) / c_sd); r_pct <- rank(cvec) / G
set.seed(20261112L); s_perm <- s_soft[sample.int(G)]
qs <- c(0, .1, .25, .5, .75, .9, 1)
diag4 <- rbind(data.frame(quantity = "c (Onnela)", t(stats::setNames(stats::quantile(cvec, qs), paste0("q", qs * 100)))),
               data.frame(quantity = "s (soft threshold)", t(stats::setNames(stats::quantile(s_soft, qs), paste0("q", qs * 100)))),
               data.frame(quantity = "r (percentile)", t(stats::setNames(stats::quantile(r_pct, qs), paste0("q", qs * 100)))))
diag4$mode_c <- c(c_mode, NA, NA); diag4$sd_c <- c(c_sd, NA, NA)
diag4$share_c_above_mode <- c(mean(cvec > c_mode), NA, NA)
diag4$spearman_with_mi_rowmean <- c(stats::cor(cvec, mi_rm, method = "spearman"), NA, NA)
diag4$spearman_with_pi <- c(stats::cor(cvec, q_att, method = "spearman"), NA, NA)
gq <- function(g) stats::quantile(g[UT[seq(1, length(UT), by = 97L)]], c(.1, .5, .9))
gd <- function(nm, v, f) { o <- outer(v, v, f); data.frame(g = nm, t(stats::setNames(gq(o), c("q10", "q50", "q90")))) }
gdist <- rbind(gd("B1 c_i c_j", cvec, `*`), gd("B2 max(c_i,c_j)", cvec, pmax), gd("B3 s_i s_j", s_soft, `*`), gd("B4 max(s)", s_soft, pmax), gd("B5 r_i r_j", r_pct, `*`))
utils::write.csv(diag4, file.path(opt$out, "gate_c_diagnostics.csv"), row.names = FALSE); utils::write.csv(gdist, file.path(opt$out, "gate_g_distribution.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene = genes, c = cvec, s_soft = s_soft, r_pct = r_pct, mi_rowmean = mi_rm, pi = q_att), file.path(opt$out, "gate_c_genes.csv"), row.names = FALSE)
print(diag4, row.names = FALSE); print(gdist, row.names = FALSE)
rm(W, nn, Aoff, dens); invisible(gc(verbose = FALSE))

## ---- common-scale terms and the eight blends -----------------------------------------------------------------------------
m_C <- mean(C[UT]); s_C <- stats::sd(C[UT]); m_a <- mean(Ta3r[UT]); s_a <- stats::sd(Ta3r[UT])
utils::write.csv(data.frame(term = c("CLR", "a3r"), mean = c(m_C, m_a), sd = c(s_C, s_a)), file.path(opt$out, "common_scale.csv"), row.names = FALSE)
say(sprintf("common scale: CLR mean %.4f sd %.4f; a3r mean %.4f sd %.4f", m_C, s_C, m_a, s_a))
Zc <- (C - m_C) / s_C; Za <- (Ta3r - m_a) / s_a
blend <- function(g, Zc, Za) { B <- g * Zc + (1 - g) * Za; diag(B) <- 0; B }
specs <- list(B1 = function() blend(outer(cvec, cvec), Zc, Za), B2 = function() blend(outer(cvec, cvec, pmax), Zc, Za),
              B3 = function() blend(outer(s_soft, s_soft), Zc, Za), B4 = function() blend(outer(s_soft, s_soft, pmax), Zc, Za),
              B5 = function() blend(outer(r_pct, r_pct), Zc, Za), B6 = function() blend(0.5, Zc, Za),
              B7 = function() blend(outer(s_perm, s_perm), Zc, Za),
              B8 = function() blend(outer(cvec, cvec), C / s_C, Ta3r / s_a))
for (nm in names(specs)) { score_readout(nm, specs[[nm]]()); invisible(gc(verbose = FALSE)) }
saveRDS(list(base_rate = base_rate, truths = TRUTHS, opt = opt, c_mode = c_mode, c_sd = c_sd), file.path(opt$out, "cache", "meta.rds"))
say("done")
