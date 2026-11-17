#!/usr/bin/env Rscript
# Per-method time and memory estimates for the confirmatory plan, recomputed from the REAL matrix dimensions of both arms
# (docs/plans/tcaclr_confirm_freeze.md). Pure arithmetic on measured anchors; runs nothing.
# Anchors (32 workers unless stated): RPE1 8,745 genes x 1,023 profiles: GENIE3 435.6 min peak 34.6 GB; GRNBoost2 221.7 min peak 4.3 GB;
# PIDC 1,160 s peak 10.48 GB (read 178 s, node pairs to 311 s, puc to 845 s). net4 5,950 x 536 at 16 threads: GENIE3 12,435 s, GRNBoost2 11,764 s
# (results/dream5/SCORE_MANIFEST.tsv; halved here for 32 threads, a linear-scaling assumption that overstates the speed-up).
# Scaling laws (docs/plans/tcaclr_prereg_candidates.md section D): GENIE3 per target ~ n log n sqrt(G), G targets; GRNBoost2 ~ G^2 n x tree-count noise;
# PIDC: read ~ G n, node pairs ~ G^2, puc ~ G^3, tail ~ G^2.
args <- commandArgs(trailingOnly = TRUE); GA <- as.numeric(args[1]); nA <- as.numeric(args[2]); GB <- as.numeric(args[3]); nB <- as.numeric(args[4])
gates <- function(f) { x <- utils::read.csv(f); c(sec = x$seconds[1], gb = x$peak_gb[1], G = x$genes[1]) }
gA <- gates(args[5]); gB <- gates(args[6]); out <- dir <- args[7]
anch <- list(rpe1 = c(G = 8745, n = 1023), net4 = c(G = 5950, n = 536))
g3 <- function(G, n, a) { f <- (G / a[["G"]]) * ((n * log(n)) / (a[["n"]] * log(a[["n"]]))) * sqrt(G / a[["G"]]); f }
gb <- function(G, n, a) (G / a[["G"]])^2 * (n / a[["n"]])
pidc <- function(G, n) { r <- c(read = 178 * (G * n) / (8745 * 1023), nodes = 133 * (G / 8745)^2 * (n / 1023), puc = 534 * (G / 8745)^3, tail = 315 * (G / 8745)^2); sum(r) }
rows <- list()
for (arm in c("A", "B")) { G <- if (arm == "A") GA else GB; n <- if (arm == "A") nA else nB; g <- if (arm == "A") gA else gB
  t3r <- 435.6 * g3(G, n, anch$rpe1); t3n <- 12435 / 60 / 2 * g3(G, n, anch$net4) * 0 + (12435 / 60 / 2) * g3(G, n, anch$net4)
  tbr <- 221.7 * gb(G, n, anch$rpe1); tbn <- (11764 / 60 / 2) * gb(G, n, anch$net4); tp <- pidc(G, n) / 60
  m3 <- 34.6 * (G * n) / (8745 * 1023); mb <- 4.3 * (G / 8745); mp <- 10.48 * (G / 8745)^2
  mk <- function(method, lo, hi, mem, note) data.frame(arm = arm, genes = G, profiles = n, method = method, minutes_low = round(lo, 1), minutes_high = round(hi, 1),
    hours_central = round(sqrt(lo * hi) / 60, 2), peak_gb_scaled = round(mem, 1), declare_gb = round(1.4 * mem + 0.49), note = note)
  rows[[length(rows) + 1]] <- mk("GENIE3 (32 workers)", min(t3r, t3n), max(t3r, t3n), m3, "RPE1 and net4 anchors")
  rows[[length(rows) + 1]] <- mk("GRNBoost2 (32 workers)", min(tbr, tbn), max(tbr, tbn), mb, "RPE1 and net4 anchors; tree-count noise not included (net4 ran 1.8x the formula)")
  rows[[length(rows) + 1]] <- mk("PIDC (tools/pidc_chunked.jl, 32 threads)", tp * 0.9, tp * 1.2, mp, "RPE1 component scaling")
  # cheap block: measured gate job (MI + CLR + operator + P^9, 8 threads) scaled by G^2 n (MI) from Arm B's measured run when arm A
  rows[[length(rows) + 1]] <- mk("CLR+aCLR+tcaCLR+|r| fit incl. P^9 (8 threads), measured/scaled", g[["sec"]] / 60 * (if (arm == "A") 1 else 1), g[["sec"]] / 60 * 1.5, g[["gb"]], "gate job wall and VmHWM (arm's own measurement)")
  rows[[length(rows) + 1]] <- mk("one experiment refit of the cheap block (half the profiles)", g[["sec"]] / 60 * 0.4, g[["sec"]] / 60 * 0.7, g[["gb"]], "x100 replicates for SE_exp")
  rows[[length(rows) + 1]] <- mk("one experiment refit of PIDC (half the profiles)", tp * 0.5, tp * 0.8, mp, "x100 replicates: this is about 15 to 35 h per arm, not 'seconds to minutes'") }
res <- do.call(rbind, rows); utils::write.csv(res, file.path(out, "schedule_estimates.csv"), row.names = FALSE); print(res, row.names = FALSE)
