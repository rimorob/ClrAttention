#!/usr/bin/env Rscript
# 907-chip M3D compendium: point estimates plus replicate-aware bootstrap.
#
# The 907 chips are technical/biological replicates of 466 experiments (207
# experiments have 1 chip, 107 have 2, 135 have 3, 17 have 4-5). Treating chips
# as independent samples is pseudo-replication. Two resampling schemes, both
# recomputing MI/CLR/attention from scratch on every draw:
#   replicate  one randomly chosen chip per experiment (all 466 experiments):
#              technical-noise variability only; no pseudo-replication and
#              no duplicated columns.
#   cluster    466 experiments drawn WITH replacement (conditions resampled),
#              one chip per occurrence, taking a different replicate chip on
#              repeated draws while any remain (duplicate columns only when
#              an experiment is drawn more often than it has chips): condition
#              + technical variability -- the headline error bars.
# Point estimates are reported on (i) all 907 chips (as in the 2007 paper) and
# (ii) the 466 replicate-averaged experiments (analysis/followup.R).
#
# Methods (fixed in advance; depths are the held-out-lab t* from
# results/avg_followup): |Pearson|, raw MI, CLR parity 2007 (10 bins,
# Euclidean), CLR HG-Stouffer (default), soft10 attention at t = 4, 8*, 16,
# pearson_top10 attention at t = 6*.
#
# Usage (repo root):
#   Rscript analysis/chips_bootstrap.R [--Bcluster 40] [--Brep 20]
#       [--out results/chips_bootstrap] [--threads N] [--quick 0] [--workers W]
# Draws run in parallel (foreach over local cores minus 2, capped by RAM; see
# analysis/parallel.R). Each finished draw is written to draws/<draw>.csv at
# once, so a partial run is usable and a restart skips finished draws.
# Creating the file jobs/control/stop makes workers skip remaining draws.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            Bcluster = 40L, Brep = 20L, out = "results/chips_bootstrap",
            threads = NULL, quick = 0L, seed = 20260924L, alpha = 0.5,
            workers = NULL, mem_gb = 3.3)
ints <- c("Bcluster", "Brep", "threads", "quick", "seed", "workers")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% ints) as.integer(args[i + 1L]) else
    if (key == "mem_gb") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "bootstrap.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))
set.seed(opt$seed)

## ---- data ------------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
rd <- function(f) {
  raw <- utils::read.delim(file.path(opt$m3d, f), check.names = FALSE,
                           stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
  rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
  X
}
Xc <- rd(paste0(base, "_chips907probes4297.tab"))
Xa <- rd(paste0("avg_", base, "_exps466probes4297.tab"))
stopifnot(identical(rownames(Xc), rownames(Xa)))
keep <- apply(Xc, 1, function(x) diff(range(x))) > 0 &
  apply(Xa, 1, function(x) diff(range(x))) > 0
if (opt$quick > 0L) {
  v <- apply(Xa, 1, stats::var); v[!keep] <- -Inf
  keep <- seq_len(nrow(Xa)) %in% order(-v)[seq_len(opt$quick)]
}
Xc <- Xc[keep, ]; Xa <- Xa[keep, ]
genes <- rownames(Xc); G <- length(genes)
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                          stringsAsFactors = FALSE)
exp_of <- stats::setNames(desc$experiment_name, desc$chip_name)[colnames(Xc)]
stopifnot(!anyNA(exp_of))
chips_by_exp <- split(colnames(Xc), exp_of)
say(sprintf("%d genes; %d chips of %d experiments (replicates per experiment: %s)",
            G, ncol(Xc), length(chips_by_exp),
            paste(names(table(lengths(chips_by_exp))), table(lengths(chips_by_exp)),
                  sep = "x", collapse = " ")))

to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck)
  p <- comembership_pairs(reg, operon_of, genes)
  list(reg = reg, pairs = p, co = comember_matrix(reg, G),
       idx = cbind(p$i, p$j))
})

## ---- one full evaluation of a sample matrix -----------------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
# Attention mass (P^t + P^t')/2 of the lazy operator P = (1-a)I + aA at each
# depth in ts, handed to fn(S, t) as soon as it exists (only P^t and one
# output matrix are alive at a time). P is built sparse without forming the
# dense (1-a)I + aA; values are identical to the dense construction.
attention <- function(A, ts, fn) {
  P <- Matrix::Diagonal(G, 1 - opt$alpha) +
    opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- NULL
  for (t in seq_len(max(ts))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- Pt + t(Pt); S <- S / 2; diag(S) <- 0; fn(S, t); rm(S) }
    par_release()
  }
  invisible(NULL)
}
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
evaluate <- function(X, label) {
  t0 <- Sys.time()
  rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0)) {  # a gene constant in this draw: tiny jitter, logged
    X[rng <= 0, ] <- X[rng <= 0, ] + matrix(stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9),
                                            sum(rng <= 0))
  }
  # Each score matrix (G x G doubles, ~150 MB at G = 4297) is scored as soon
  # as it exists and then dropped, keeping peak memory per worker low.
  rows <- list()
  score <- function(S, m) {
    diag(S) <- 0
    for (ev in names(bench)) {
      b <- bench[[ev]]
      ap <- pr_summary(S[b$idx], b$pairs$label)[["aupr"]]
      coh <- stats::median(regulon_coherence(S, b$reg, b$co, operon_of, genes),
                           na.rm = TRUE)
      rows[[length(rows) + 1L]] <<- data.frame(draw = label, evidence = ev,
                                               method = m, aupr = ap, coherence = coh)
    }
    rm(S); par_release()   # collect now rather than when R's GC gets to it
  }
  Z <- Zs(X); R <- tcrossprod(Z) / (ncol(Z) - 1); rm(Z)
  score(abs(R), "abs_pearson")
  Ar <- abs(R); diag(Ar) <- 0; rm(R)
  Apear <- clr:::.topk_rows(Ar, 10L); rm(Ar)
  par <- ClrAttention$new(X)$estimate_mi(bins = 10, transform = "none",
                                         threads = nthr())
  score(par$calibrate(method = "normal", combine = "euclidean")$clr_scores,
        "clr_parity2007")
  par$release(); rm(par)
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none",
                                       threads = nthr())
  score(f$mi, "mi")
  f$calibrate(method = "normal", combine = "stouffer")
  score(f$clr_scores, "clr_hg_stouffer")
  f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
  attention(f$operator, c(4, 8, 16), function(S, t) score(S, paste0("soft10_t", t)))
  f$set_operator(Apear, alpha = opt$alpha, label = "pearson_top10")
  attention(f$operator, 6, function(S, t) score(S, "pearson_top10_t6"))
  f$release(); rm(f, Apear); par_release()
  out <- do.call(rbind, rows)
  say(sprintf("%-14s N=%d  %.0f s  | SC AUPR: CLR %.4f  soft10_t8 %.4f  pearson_t6 %.4f",
              label, ncol(X), as.numeric(difftime(Sys.time(), t0, units = "secs")),
              out$aupr[out$evidence == "SC" & out$method == "clr_hg_stouffer"],
              out$aupr[out$evidence == "SC" & out$method == "soft10_t8"],
              out$aupr[out$evidence == "SC" & out$method == "pearson_top10_t6"]))
  out
}
draws_f <- file.path(opt$out, "draws.csv")
draws_dir <- file.path(opt$out, "draws")
dir.create(draws_dir, showWarnings = FALSE)
done <- unique(c(if (file.exists(draws_f)) unique(utils::read.csv(draws_f)$draw),
          sub("\\.csv$", "", list.files(draws_dir, "\\.csv$"))))

## ---- task list: point estimates, then resampling draws ------------------------
draw_cols <- function(mode) {
  ex <- names(chips_by_exp)
  pick <- if (mode == "replicate") ex else sample(ex, length(ex), replace = TRUE)
  used <- list()
  vapply(pick, function(e) {
    ch <- chips_by_exp[[e]]
    avail <- setdiff(ch, used[[e]])
    if (!length(avail)) avail <- ch
    c1 <- if (length(avail) == 1L) avail else sample(avail, 1L)
    used[[e]] <<- c(used[[e]], c1)
    c1
  }, "")
}
# Column sets are drawn in the master with the same per-draw seeds as the
# earlier sequential version, so draws already on disk stay valid.
tasks <- list(list(label = "point_chips907", set = "chips", cols = colnames(Xc), seed = opt$seed),
              list(label = "point_avg466", set = "avg", cols = colnames(Xa), seed = opt$seed))
for (mode in c("cluster", "replicate")) {
  B <- if (mode == "cluster") opt$Bcluster else opt$Brep
  for (b in seq_len(B)) {
    sd_b <- opt$seed + 1000L * (mode == "replicate") + b
    set.seed(sd_b)
    tasks[[length(tasks) + 1L]] <- list(label = sprintf("%s_%03d", mode, b), set = "chips",
                                        cols = draw_cols(mode), seed = sd_b)
  }
}
tasks <- Filter(function(t) !t$label %in% done, tasks)
say(sprintf("%d draws already done; %d to run", length(done), length(tasks)))

source(file.path(script_dir, "parallel.R"))
if (length(tasks)) {
  pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
  st <- foreach(task = tasks, .inorder = FALSE, .errorhandling = "pass") %dopar% {
    if (par_should_stop()) return("stopped")
    set.seed(task$seed)
    Xm <- if (task$set == "avg") Xa else Xc[, task$cols, drop = FALSE]
    d <- evaluate(Xm, task$label); rm(Xm)
    tmp <- file.path(draws_dir, paste0(task$label, ".csv.part"))
    utils::write.csv(d, tmp, row.names = FALSE)
    file.rename(tmp, file.path(draws_dir, paste0(task$label, ".csv")))  # atomic
    par_release()
    "ok"
  }
  par_stop(pc)
  bad <- vapply(st, function(x) inherits(x, "error"), NA)
  if (any(bad)) say("draws failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
  if (any(unlist(st[!bad]) == "stopped")) say("stop file seen: remaining draws skipped")
}
# Merge per-draw files into draws.csv (the single file the summary reads).
parts <- list.files(draws_dir, "\\.csv$", full.names = TRUE)
old <- if (file.exists(draws_f)) utils::read.csv(draws_f) else NULL
new <- do.call(rbind, lapply(parts, utils::read.csv))
all_d <- rbind(old, new[!new$draw %in% old$draw, ])
utils::write.csv(all_d, draws_f, row.names = FALSE)

## ---- summary ---------------------------------------------------------------
d <- utils::read.csv(draws_f)
d$mode <- sub("_.*$", "", d$draw)
ref <- d[d$method == "clr_hg_stouffer", c("draw", "evidence", "aupr", "coherence")]
names(ref)[3:4] <- c("aupr_ref", "coh_ref")
d <- merge(d, ref, by = c("draw", "evidence"))
d$d_aupr <- d$aupr - d$aupr_ref; d$d_coh <- d$coherence - d$coh_ref
q <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)
summ <- do.call(rbind, lapply(split(d, list(d$mode, d$evidence, d$method), drop = TRUE),
  function(x) {
    boot <- x$mode[1] %in% c("cluster", "replicate")
    data.frame(mode = x$mode[1], evidence = x$evidence[1], method = x$method[1],
               n = nrow(x), aupr = mean(x$aupr),
               aupr_lo = if (boot) q(x$aupr, .025) else NA,
               aupr_hi = if (boot) q(x$aupr, .975) else NA,
               d_aupr = mean(x$d_aupr),
               d_aupr_lo = if (boot) q(x$d_aupr, .025) else NA,
               d_aupr_hi = if (boot) q(x$d_aupr, .975) else NA,
               coherence = mean(x$coherence),
               d_coh = mean(x$d_coh),
               d_coh_lo = if (boot) q(x$d_coh, .025) else NA,
               d_coh_hi = if (boot) q(x$d_coh, .975) else NA)
  }))
summ <- summ[order(summ$mode, summ$evidence, -summ$aupr), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
for (m in unique(summ$mode)) for (ev in c("SC", "all")) {
  s <- summ[summ$mode == m & summ$evidence == ev, ]
  say(sprintf("== %s [%s], n = %d", m, ev, s$n[1]))
  for (k in seq_len(nrow(s)))
    say(sprintf("   %-18s AUPR %.4f  dAUPR vs CLR %+.4f [%+.4f, %+.4f]  coherence %.3f  dCoh %+.3f [%+.3f, %+.3f]",
                s$method[k], s$aupr[k], s$d_aupr[k], s$d_aupr_lo[k], s$d_aupr_hi[k],
                s$coherence[k], s$d_coh[k], s$d_coh_lo[k], s$d_coh_hi[k]))
}
say("done")
