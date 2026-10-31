#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Follow-up to sim_width.R: the PI's "background" idea applied
# to CLR's own standardization instead of to the width. CLR z-scores each MI row
# by its mean and SD; a gene with many strong partners (the ribosome) inflates
# its own row SD and compresses its z-scores. rCLR replaces the Gaussian row
# background by a robust one estimated from the bulk only: cube-root MI,
# symmetric two-way median polish, per-gene MAD, Stouffer (EB_additive_score in
# sim_width.R). Here the method's attention (softmax K10/K25, cap 50, lazy walk,
# t = 9) is built on rCLR instead of CLR, plus rCLR + reciprocity.
# Same generator, seeds and metrics as sim_width.R (functions sourced from it).
# Usage: Rscript sim_rclr.R [seeds=1:4] [scenario=messy]
args <- commandArgs(TRUE)
SEEDS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:4
SCEN <- if (length(args) >= 2) args[2] else "messy"
src <- readLines("analysis/synthetic_width/sim_width.R")
src <- src[seq_len(grep("^## ---- run", src) - 1L)]
src <- src[!grepl("^args <-|^SEEDS <-|^SCEN <-|^NSCALE <-|^OUT <-|^dir.create", src)]
eval(parse(text = src))
OUT <- file.path("analysis/synthetic_width/out", SCEN)
res <- list()
for (seed in SEEDS) {
  t0 <- Sys.time(); sd0 <- 20261030L + seed
  sim <- simulate(sd0); X <- sim$X; tr <- make_truth(sim)
  M <- mi_of(X); C <- clr_of(M); U <- eb_additive(M)$U; rm(M)
  rows <- list(); add <- function(name, S) rows[[length(rows) + 1L]] <<- cbind(op = name, score(S, tr))
  ops <- list(K10 = selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L)),
              rCLR_K10 = selfloop(clr:::.softmax_rows(U, keff = 10, cap = 50L)),
              rCLR_K25 = selfloop(clr:::.softmax_rows(U, keff = 25, cap = 50L)))
  ops$rCLR_K10_RECIP <- recip(ops$rCLR_K10)
  for (nm in names(ops)) { P9 <- walk9(lazy(ops[[nm]]))
    add(paste0(nm, "|att_t9"), att(P9)); add(paste0(nm, "|a3r_t9"), a3r(P9)) }
  r <- do.call(rbind, rows); r$seed <- seed; res[[length(res) + 1L]] <- r
  say(sprintf("%s seed %d: %.0f s", SCEN, seed, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- do.call(rbind, res)
utils::write.csv(res, file.path(OUT, "metrics_rclr_by_seed.csv"), row.names = FALSE)
mets <- c("cpx_x", "s5_x", "s20_x", "s80_x", "small_x", "mid_x", "reg_x")
ref <- res[res$op == "K10|att_t9", ]
pr <- do.call(rbind, lapply(split(res, factor(res$op, unique(res$op))), function(d) {
  d <- d[match(ref$seed, d$seed), ]
  data.frame(op = d$op[1], t(sapply(mets, function(v) { q <- d[[v]] / ref[[v]]
    sprintf("%.1f | %.2f (%d/%d)", mean(d[[v]]), exp(mean(log(q))), sum(q > 1), length(q)) }))) }))
utils::write.csv(pr, file.path(OUT, "rclr_vs_K10.csv"), row.names = FALSE)
print(pr, row.names = FALSE)
say("done")
