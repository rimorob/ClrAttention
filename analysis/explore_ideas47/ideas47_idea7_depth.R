#!/usr/bin/env Rscript
# EXPLORATORY, IN-SAMPLE, NOT PREREGISTERED. Idea 7 (soft MCL inflation between walk steps), depth stage:
# label-free held-out criterion (analysis/lab_splits.R as run by analysis/robust_z/robust_z.R on K562) over the grid
# r in {1, 1.1, 1.2, 1.5} x t = 1..20; selection of t*_r and of the joint (r*, t*). No truth is read.
# docs/plans/explore_ideas47_plan.md section 2. Usage: Rscript ideas47_idea7_depth.R [--out results/explore_ideas47/idea7] [--workers 6] [--nsplit 6] [--tmax 20] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(out = "results/explore_ideas47/idea7", workers = 6L, nsplit = 6L, tmax = 20L, quick = 0L, depth_seed = 20261014L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source("analysis/robust_z/robust_z_lib.R"); source("analysis/regulondb.R")
source(file.path(script_dir, "ideas47_inflate.R"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = file.path(opt$out, "depth.log"), append = TRUE) }
raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (opt$quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
G <- nrow(X); lab <- colnames(X); UT <- which(upper.tri(matrix(FALSE, G, G)))
RS <- c(1, 1.1, 1.2, 1.5)
lab_split <- function() {                      # analysis/robust_z/robust_z.R, verbatim
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab]) }
set.seed(opt$depth_seed); splits <- lapply(seq_len(opt$nsplit), function(s) lab_split())
depth_split <- function(si) {
  h <- splits[[si]]; set.seed(opt$depth_seed + si); t0 <- Sys.time()
  Ms <- lapply(1:2, function(d) { Xm <- X[, h == d]; rng <- apply(Xm, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9)
    bspline_mi(Xm, bins = "hg", spline_order = 3, threads = 1, transform = "none") })
  Cs <- lapply(1:2, function(d) rz_scores(Ms[[d]], "V0")); rm(Ms)
  rows <- list()
  for (dir in 1:2) { te <- 3L - dir
    sv <- Cs[[te]][UT]; target <- as.integer(sv >= stats::quantile(sv, 0.99)); rm(sv)
    P <- rz_lazy_P(rz_operator(Cs[[dir]]))
    for (r in RS) { crit <- numeric(opt$tmax)
      inflate_walk(P, r, opt$tmax, function(t, Q) { S <- (Q + t(Q)) / 2; crit[t] <<- ap_fast(S[UT], target); rm(S) })
      rows[[length(rows) + 1L]] <- data.frame(split = si, direction = dir, r = r, t = seq_len(opt$tmax), criterion = crit)
      invisible(gc(verbose = FALSE)) }
    rm(P, target) }
  out <- do.call(rbind, rows); out$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); out }
say(sprintf("depth stage: %d splits, r = %s, t = 1..%d, %d workers", opt$nsplit, paste(RS, collapse = ","), opt$tmax, min(opt$workers, opt$nsplit)))
invisible(gc(verbose = FALSE))
res <- parallel::mclapply(seq_len(opt$nsplit), depth_split, mc.cores = min(opt$workers, opt$nsplit), mc.preschedule = FALSE)
bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), NA)
if (any(bad)) stop("depth stage failed for splits: ", paste(which(bad), collapse = ","), " ", paste(unlist(res[bad]), collapse = " | "))
cur <- do.call(rbind, res); say(sprintf("depth stage done; per-split seconds: %s", paste(sprintf("%.0f", vapply(res, function(x) x$secs[1], 0)), collapse = " ")))
cur <- cur[, setdiff(names(cur), "secs")]; utils::write.csv(cur, file.path(opt$out, "heldout_curves.csv"), row.names = FALSE)

## ---- gate (d), part 2: r = 1 curves equal the robust_z V0 curves ------------------------------------------------------------
if (opt$quick == 0L) {
  ref <- utils::read.csv("results/robust_z/k562/heldout_curves.csv", stringsAsFactors = FALSE); ref <- ref[ref$variant == "V0", ]
  mm <- merge(cur[cur$r == 1, ], ref, by = c("split", "direction", "t"), suffixes = c("", ".ref"))
  stopifnot(nrow(mm) == nrow(ref)); dmax <- max(abs(mm$criterion - mm$criterion.ref))
  utils::write.csv(data.frame(check = "(d) r = 1 held-out criterion curves vs robust_z V0 heldout_curves.csv (6 splits x 2 directions x 20 t)", max_abs_diff = dmax, pass = dmax < 1e-9),
                   file.path(opt$out, "gate_d_curves.csv"), row.names = FALSE)
  say(sprintf("gate (d) curves: max abs diff %.3g", dmax)); if (dmax >= 1e-9) stop("EQUIVALENCE GATE (d) FAILED for the r = 1 depth curves")
}
## ---- selection ----------------------------------------------------------------------------------------------------------------
ts <- do.call(rbind, lapply(split(cur, list(cur$split, cur$direction, cur$r), drop = TRUE), function(x) {
  k <- which.max(x$criterion); data.frame(split = x$split[1], direction = x$direction[1], r = x$r[1], t_best = x$t[k], crit_best = x$criterion[k]) }))
utils::write.csv(ts, file.path(opt$out, "heldout_tstar.csv"), row.names = FALSE)
tsplit <- do.call(rbind, lapply(split(ts, list(ts$split, ts$r), drop = TRUE), function(x) data.frame(split = x$split[1], r = x$r[1], t_star = as.integer(round(exp(mean(log(x$t_best))))))))
tr <- do.call(rbind, lapply(split(tsplit, tsplit$r), function(x) data.frame(r = x$r[1], n_splits = nrow(x), t_star_r = as.integer(floor(stats::median(x$t_star) + 0.5)),
  t_q25 = stats::quantile(x$t_star, 0.25, names = FALSE), t_q75 = stats::quantile(x$t_star, 0.75, names = FALSE), t_min = min(x$t_star), t_max = max(x$t_star))))
tr <- tr[order(tr$r), ]; utils::write.csv(tr, file.path(opt$out, "tstar_by_r.csv"), row.names = FALSE)
mc <- stats::aggregate(criterion ~ r + t, cur, mean); mc <- mc[order(mc$r, mc$t), ]
utils::write.csv(mc, file.path(opt$out, "mean_criterion_grid.csv"), row.names = FALSE)
best <- mc[mc$criterion == max(mc$criterion), ]; best <- best[which.min(best$r), ]
joint <- data.frame(r_star = best$r, t_star = best$t, mean_criterion = best$criterion, tie_break = "smaller r")
utils::write.csv(joint, file.path(opt$out, "joint_selection.csv"), row.names = FALSE)
print(tr, row.names = FALSE); print(joint, row.names = FALSE)
say(sprintf("t*_r: %s; joint (r*, t*) = (%.1f, %d)", paste(sprintf("r%.1f=%d", tr$r, tr$t_star_r), collapse = " "), joint$r_star, joint$t_star))
say("done")
