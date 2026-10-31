#!/usr/bin/env Rscript
# SYNTHETIC ONLY. Follow-up to sim_width.R: a realizable detector for the
# "symmetric sink" (partner-less genes neither send nor receive attention).
# Row gate from cross-half reproducibility: repro_n_i = number of row i's top-50
# positive CLR candidates shared by the two lab halves. Under no shared signal
# the overlap is about Poisson(lambda0 = 50 * 50 / (G - 1)); the gate is open
# (rho_i = 1) when repro_n_i >= k0, the smallest k with P(Poisson(lambda0) >= k)
# < 0.05 (k0 = 4 at G = 2,000), else closed (rho_i = 0). Unlike REPRO10/25 in
# sim_width.R, the gate does not filter candidates: an open row attends exactly
# as K10 does on the full data.
# Same generator, seeds, metrics as sim_width.R (functions sourced from it).
# Usage: Rscript sim_gate.R [seeds=1:4] [scenario=messy]
args <- commandArgs(TRUE)
SEEDS <- if (length(args) >= 1) eval(parse(text = args[1])) else 1:4
SCEN <- if (length(args) >= 2) args[2] else "messy"
src <- readLines("analysis/synthetic_width/sim_width.R")
src <- src[seq_len(grep("^## ---- run", src) - 1L)]
src <- src[!grepl("^args <-|^SEEDS <-|^SCEN <-|^NSCALE <-|^OUT <-|^dir.create", src)]
eval(parse(text = src))
OUT <- file.path("analysis/synthetic_width/out", SCEN)
lambda0 <- 50 * 50 / (G - 1); k0 <- min(which(stats::ppois(0:50 - 1, lambda0, lower.tail = FALSE) < 0.05)) - 1L
res <- list(); gates <- list()
for (seed in SEEDS) {
  t0 <- Sys.time(); sd0 <- 20261030L + seed
  sim <- simulate(sd0); X <- sim$X; tr <- make_truth(sim)
  C <- clr_of(mi_of(X)); hA <- sim$lab <= 5L
  TA <- top50(clr_of(mi_of(X[, hA]))); TB <- top50(clr_of(mi_of(X[, !hA])))
  rep_n <- rowSums(TA & TB); rm(TA, TB); rho <- as.numeric(rep_n >= k0)
  K10 <- selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L))
  ops <- list(K10 = K10, GATE_sink2 = sink2(K10, rho), GATE_sink2_RECIP = sink2(recip(K10), rho),
              ORACLE_sink2 = sink2(K10, as.numeric(sim$role != "background")))
  rows <- list(); add <- function(name, S) rows[[length(rows) + 1L]] <<- cbind(op = name, score(S, tr))
  for (nm in names(ops)) { P9 <- walk9(lazy(ops[[nm]]))
    add(paste0(nm, "|att_t9"), att(P9)); add(paste0(nm, "|a3r_t9"), a3r(P9)) }
  r <- do.call(rbind, rows); r$seed <- seed; res[[length(res) + 1L]] <- r
  roleg <- ifelse(sim$role %in% c("complex", "ribosome"), ifelse(sim$csize >= 20, "complex_>=20", "complex_<20"), sim$role)
  gates[[length(gates) + 1L]] <- data.frame(seed = seed, role = roleg, open = rho)
  say(sprintf("%s seed %d: %.0f s; k0 = %d; gate open for %.1f%% of genes", SCEN, seed,
              as.numeric(difftime(Sys.time(), t0, units = "secs")), k0, 100 * mean(rho)))
}
res <- do.call(rbind, res); gates <- do.call(rbind, gates)
utils::write.csv(res, file.path(OUT, "metrics_gate_by_seed.csv"), row.names = FALSE)
go <- stats::aggregate(open ~ role, data = gates, FUN = mean)
utils::write.csv(go, file.path(OUT, "gate_open_by_role.csv"), row.names = FALSE)
mets <- c("cpx_x", "s5_x", "s20_x", "s80_x", "small_x", "mid_x", "reg_x")
ref <- res[res$op == "K10|att_t9", ]
pr <- do.call(rbind, lapply(split(res, factor(res$op, unique(res$op))), function(d) {
  d <- d[match(ref$seed, d$seed), ]
  data.frame(op = d$op[1], t(sapply(mets, function(v) { q <- d[[v]] / ref[[v]]
    sprintf("%.1f | %.2f (%d/%d)", mean(d[[v]]), exp(mean(log(q))), sum(q > 1), length(q)) }))) }))
utils::write.csv(pr, file.path(OUT, "gate_vs_K10.csv"), row.names = FALSE)
print(go, row.names = FALSE); print(pr, row.names = FALSE)
say("done")
