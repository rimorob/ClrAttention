#!/usr/bin/env Rscript
# Timing and memory of the scoring stages and of one refit at the sizes of the two arms, on SYNTHETIC data only (random scores with a weak planted
# signal, planted modules in the expression matrix). No real data. Usage (repo root):
#   Rscript analysis/tcaclr_confirm/scale_test.R --mode score|refit --G 9000 --n_truth 3500 --n_prof 1500 --workers 16 --out DIR
suppressPackageStartupMessages(source("analysis/tcaclr_confirm/tc_lib.R"))
opt <- tc_options(commandArgs(TRUE), list(mode = "score", G = 9000L, n_truth = 3500L, n_prof = 1500L, workers = 16L, out = "", n_readouts = 4L))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE); log <- file.path(opt$out, paste0("scale_", opt$mode, ".log")); say <- function(...) tc_say(log, ...)
G <- opt$G; gn <- sprintf("S%05d", seq_len(G)); set.seed(99)
if (opt$mode == "refit") {
  n_h <- floor(opt$n_prof / 2); X <- matrix(rnorm(G * n_h), G, n_h, dimnames = list(gn, NULL)); for (k in 1:40) { m <- ((k - 1) * 30 + 1):(k * 30); X[m, ] <- X[m, ] + matrix(rep(rnorm(n_h, sd = 1), each = 30), 30) }
  say(sprintf("refit test: %d genes x %d profiles (one half of %d)", G, n_h, opt$n_prof)); t0 <- Sys.time(); R <- tc_fit_three(X, opt$workers); t1 <- Sys.time(); R$absr <- tc_absr(X)
  ix <- cbind(sample(G, 5e6, TRUE), sample(G, 5e6, TRUE)); for (m in names(R)) invisible(R[[m]][ix])
  say(sprintf("tc_fit_three %.0f s, absr + gather %.0f s; peak RSS %.2f GB", as.numeric(difftime(t1, t0, units = "secs")), as.numeric(difftime(Sys.time(), t1, units = "secs")), tc_vmhwm_gb())); quit(save = "no")
}
# score mode: synthetic truth of n_truth genes in net4-like complex sizes, seven method files plus read-outs, then the real driver
sizes <- c(80, 70, 60, 45, 30, 25, 22, pmax(3L, round(exp(rnorm(600, 1.5, 0.75))))); sizes <- sizes[cumsum(sizes) <= opt$n_truth]
tg <- gn[sample(G, opt$n_truth)]; at <- 1; cx <- list(); for (k in seq_along(sizes)) { cx[[k]] <- tg[at:(at + sizes[k] - 1)]; at <- at + sizes[k] }
mm <- do.call(rbind, lapply(seq_along(cx), function(k) data.frame(variant = "main", dataset = "syn", complex = sprintf("CX%04d", k), orf = cx[[k]])))
td <- file.path(opt$out, "truth"); dir.create(td, showWarnings = FALSE); utils::write.table(mm, file.path(td, "members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(scope = "all", orf_a = cx[[1]][1], orf_b = cx[[1]][2], reason = "paralog_set"), file.path(td, "excluded.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste(tc_sha256(file.path(td, "members.tsv")), "members.tsv"), paste(tc_sha256(file.path(td, "excluded.tsv")), "excluded.tsv")), file.path(td, "t.sha256"))
sd_ <- file.path(opt$out, "scores"); dir.create(sd_, showWarnings = FALSE); writeLines(gn, file.path(sd_, "genes.txt")); np <- G * (G - 1) / 2
idx_all <- unlist(lapply(cx, function(m) { mi <- sort(match(m, gn)); a <- rep(mi, times = length(mi)); b <- rep(mi, each = length(mi)); s <- a < b; tc_pair_index(a[s], b[s]) }))   # within-complex pairs, once
for (m in c(TC$methods, paste0("rewire_", sprintf("%02d", seq_len(opt$n_readouts))))) { v <- rnorm(np); v[idx_all] <- v[idx_all] + 0.5; tc_write_bin(v, file.path(sd_, paste0(m, ".scores"))); rm(v) }
say(sprintf("synthetic score directory: %d genes, %d methods/read-outs, truth %d genes in %d complexes", G, length(TC$methods) + opt$n_readouts, length(unique(unlist(cx))), length(cx)))
for (stage in c("gates", "score", "report")) { t0 <- Sys.time()
  rc <- system2("Rscript", c("analysis/tcaclr_confirm/score_confirm.R", "--stage", stage, "--arm", "A", "--scores", sd_, "--truth_dir", td, "--members", file.path(td, "members.tsv"), "--excluded", file.path(td, "excluded.tsv"),
    "--sha", file.path(td, "t.sha256"), "--out", file.path(opt$out, "run"), "--workers", opt$workers, "--validation", "1"), stdout = FALSE, stderr = FALSE)
  say(sprintf("stage %-6s rc %d  %.0f s", stage, rc, as.numeric(difftime(Sys.time(), t0, units = "secs")))) }
l <- readLines(file.path(opt$out, "run", "run_score.log")); say("score log peak line: ", tail(grep("peak RSS", l, value = TRUE), 1)); say("per-method lines: ", paste(sub("^\\[[^]]*\\] ", "", grep("^\\[.*method ", l, value = TRUE)), collapse = " | "))
