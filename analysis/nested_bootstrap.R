#!/usr/bin/env Rscript
# Nested cluster bootstrap on the 907-chip M3D compendium: the depth is
# re-selected inside every draw (docs/COMPUTE_CONSTRAINTS.md item 2), and the
# TF-node benchmark is scored in every draw (item 8, E. coli part).
#
# chips_bootstrap.R scores attention at depths fixed in advance from the
# full-data held-out search, so its intervals leave out the uncertainty of
# choosing t. Here every draw repeats the whole procedure:
#   1. resample the 466 experiments with replacement, one replicate chip per
#      occurrence -- the same columns as chips_bootstrap.R's cluster draws
#      (same per-draw seeds), so the two runs pair draw by draw;
#   2. split the draw's labs at random into two balanced halves (an
#      experiment drawn twice stays in one half);
#   3. in both directions, build the operator on one half and evaluate the
#      held-out criterion (AUPR of attention mass for the other half's top-1%
#      CLR pairs; RegulonDB never used) at t = 1..32, 40, 48, 56, 64;
#      t*_draw = geometric mean of the two argmax depths;
#   4. fit the whole draw and score CLR, |Pearson|, raw MI, soft10 attention
#      at t*_draw and at the fixed t = 8, and pearson_top10 at t*_draw and at
#      t = 6, on (a) regulon co-membership AUPR, (b) median regulon coherence
#      and (c) the TF-node benchmark (pooled AUPR ratio, median per-TF AUROC,
#      share of TFs where the method beats CLR), at strong/confirmed and at
#      all evidence.
# Differences against CLR are paired within each draw; intervals are the
# 2.5-97.5% percentiles over draws.
#
# Usage (repo root):
#   Rscript analysis/nested_bootstrap.R [--B 500] [--out results/nested_bootstrap]
#       [--quick 0] [--workers W] [--mem_gb 4.5]
# Each finished draw is written to draws/<draw>.csv (and its depth search to
# tstar/<draw>.csv) at once; a restart skips finished draws.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            B = 500L, out = "results/nested_bootstrap", threads = NULL, quick = 0L,
            seed = 20260924L, alpha = 0.5, target_frac = 0.01,
            tgrid = "1:32,40,48,56,64", workers = NULL, mem_gb = 4.5)
ints <- c("B", "threads", "quick", "seed", "workers")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else
    if (key %in% c("mem_gb", "alpha", "target_frac")) as.numeric(v) else v
  i <- i + 2L
}
for (d in c("draws", "tstar")) dir.create(file.path(opt$out, d), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "nested_bootstrap.log")
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
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
tgrid <- sort(unique(unlist(lapply(strsplit(opt$tgrid, ",")[[1]],
                                   function(s) eval(parse(text = s))))))

## ---- data (as chips_bootstrap.R) ----------------------------------------------
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
keep <- apply(Xc, 1, function(x) diff(range(x))) > 0 &
  apply(Xa, 1, function(x) diff(range(x))) > 0
if (opt$quick > 0L) {
  v <- apply(Xa, 1, stats::var); v[!keep] <- -Inf
  keep <- seq_len(nrow(Xa)) %in% order(-v)[seq_len(opt$quick)]
}
Xc <- Xc[keep, ]; rm(Xa)
genes <- rownames(Xc); G <- length(genes)
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_descriptions")),
                          stringsAsFactors = FALSE)
exp_of <- stats::setNames(desc$experiment_name, desc$chip_name)[colnames(Xc)]
stopifnot(!anyNA(exp_of))
chips_by_exp <- split(colnames(Xc), exp_of)
exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)
exper[!grepl("experimenter:", desc$description)] <- NA
lab_of_exp <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])
lab_of_exp[is.na(lab_of_exp)] <- paste0("unknown_", which(is.na(lab_of_exp)))
lab_of_chip <- stats::setNames(unname(lab_of_exp[exp_of]), names(exp_of))
rm(desc, exper)

to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck)
  p <- comembership_pairs(reg, operon_of, genes)
  list(reg = reg, pairs = p, co = comember_matrix(reg, G), idx = cbind(p$i, p$j))
})
# TF-node universe (as tfnode_attention.R)
tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
tf_parts <- strsplit(.classic_tf_to_genes(tfr$regulator), ";", fixed = TRUE)
tfnet <- data.frame(tf_gene = to_bn(unlist(tf_parts)),
                    target = rep(tfr$bnumber, lengths(tf_parts)),
                    confidence = rep(tfr$confidence, lengths(tf_parts)),
                    stringsAsFactors = FALSE)
tfnet <- tfnet[!is.na(tfnet$tf_gene), ]
U <- list(SC = edge_universe(tfnet, genes, conf_keep = c("S", "C")),
          all = edge_universe(tfnet, genes))
rm(rtab, tfr, tf_parts, tfnet, to_bn)
say(sprintf("%d genes; %d chips of %d experiments; %d labs; TF-node universe SC %d TFs / %d positives",
            G, ncol(Xc), length(chips_by_exp), length(unique(lab_of_chip)),
            length(U$SC$tfs), U$SC$n_pos))

## ---- constructions ------------------------------------------------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
OMP_THREADS <- NULL
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
jitter_const <- function(X) {
  rng <- apply(X, 1, function(x) diff(range(x)))
  if (any(rng <= 0))
    X[rng <= 0, ] <- X[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(X), sd = 1e-9)
  X
}
build_fit <- function(X) {
  f <- ClrAttention$new(jitter_const(X))$estimate_mi(bins = "hg", transform = "none",
                                                     threads = nthr())
  f$calibrate(method = "normal", combine = "stouffer")
  f
}
pearson_top10 <- function(X) {
  Ar <- abs(tcrossprod(Zs(X)) / (ncol(X) - 1)); diag(Ar) <- 0
  clr:::.topk_rows(Ar, 10L)
}
make_A <- function(f, X, op) {
  if (op == "soft10") { f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
                        return(f$operator) }
  f$set_operator(pearson_top10(X), alpha = opt$alpha, label = op); f$operator
}
lazy_P <- function(A) Matrix::Diagonal(nrow(A), 1 - opt$alpha) +
  opt$alpha * Matrix::Matrix(A, sparse = TRUE)
# fn(S, t) for attention mass at each t in ts (ascending)
walk <- function(P, ts, fn) {
  Pt <- NULL
  for (t in seq_len(max(ts))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- Pt + t(Pt); S <- S / 2; diag(S) <- 0; fn(S, t); rm(S) }
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
}
per_tf_auroc <- function(S, u) vapply(names(u$targets), function(tn) {
  t <- as.integer(tn); y <- integer(G); y[unique(u$targets[[tn]])] <- 1L; y[t] <- NA
  if (sum(y, na.rm = TRUE) < 3) return(NA_real_)
  s <- S[t, ]; ok <- !is.na(y)
  r <- rank(s[ok]); n1 <- sum(y[ok]); n0 <- sum(ok) - n1
  (sum(r[y[ok] == 1]) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0)
}, numeric(1))

draw_cols <- function() {     # identical to chips_bootstrap.R, mode "cluster"
  ex <- names(chips_by_exp)
  pick <- sample(ex, length(ex), replace = TRUE)
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
lab_halves <- function(cols) {
  lab <- lab_of_chip[cols]
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab])
}

## ---- one nested draw ----------------------------------------------------------
run_draw <- function(task) {
  X <- Xc[, task$cols, drop = FALSE]
  set.seed(task$seed + 500000L)
  h <- lab_halves(task$cols)
  # 3. held-out depth search on lab halves
  fits <- lapply(1:2, function(d) build_fit(X[, h == d, drop = FALSE]))
  ut <- which(upper.tri(matrix(0, G, G)))
  curves <- list()
  for (dir in 1:2) {
    sv <- fits[[3L - dir]]$clr_scores[ut]
    target <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
    for (op in c("soft10", "pearson_top10")) {
      P <- lazy_P(make_A(fits[[dir]], X[, h == dir, drop = FALSE], op))
      crit <- stats::setNames(numeric(length(tgrid)), tgrid)
      walk(P, tgrid, function(S, t) crit[[as.character(t)]] <<- ap_fast(S[ut], target))
      curves[[length(curves) + 1L]] <- data.frame(draw = task$label, direction = dir,
        operator = op, t = tgrid, criterion = unname(crit),
        n_train = sum(h == dir), n_test = sum(h != dir))
      rm(P)
    }
  }
  for (f in fits) f$release(); rm(fits, ut); par_release()
  curves <- do.call(rbind, curves)
  tbest <- tapply(seq_len(nrow(curves)), list(curves$operator, curves$direction),
                  function(ix) curves$t[ix][which.max(curves$criterion[ix])])
  tsel <- apply(tbest, 1, function(v) as.integer(round(exp(mean(log(v))))))
  # 4. full-draw fit and scoring
  rows <- list(); clr_tf <- list()
  score <- function(S, m, t) {
    diag(S) <- 0
    for (ev in names(bench)) {
      b <- bench[[ev]]; u <- U[[ev]]
      au <- per_tf_auroc(S, u)
      if (m == "clr_hg_stouffer") clr_tf[[ev]] <<- au
      ref <- clr_tf[[ev]]
      okk <- !is.na(au) & !is.na(ref)
      rows[[length(rows) + 1L]] <<- data.frame(draw = task$label, evidence = ev, method = m, t = t,
        aupr = ap_fast(S[b$idx], b$pairs$label),
        coherence = stats::median(regulon_coherence(S, b$reg, b$co, operon_of, genes), na.rm = TRUE),
        tf_aupr_ratio = ap_fast(S[cbind(u$i, u$j)], u$label) / mean(u$label),
        tf_auroc_median = stats::median(au, na.rm = TRUE),
        tf_frac_beats_clr = if (m == "clr_hg_stouffer") NA_real_ else mean(au[okk] > ref[okk]))
    }
    rm(S); par_release()
  }
  f <- build_fit(X)
  score(f$clr_scores, "clr_hg_stouffer", NA)     # first: per-TF reference
  score(abs(tcrossprod(Zs(X)) / (ncol(X) - 1)), "abs_pearson", NA)
  score(f$mi, "mi", NA)
  P <- lazy_P(make_A(f, X, "soft10"))
  ts <- sort(unique(c(8L, tsel[["soft10"]])))
  walk(P, ts, function(S, t) {
    if (t == tsel[["soft10"]]) score(S, "soft10_tsel", t)
    if (t == 8L) score(S, "soft10_t8", t)
  })
  P <- lazy_P(make_A(f, X, "pearson_top10"))
  ts <- sort(unique(c(6L, tsel[["pearson_top10"]])))
  walk(P, ts, function(S, t) {
    if (t == tsel[["pearson_top10"]]) score(S, "pearson_top10_tsel", t)
    if (t == 6L) score(S, "pearson_top10_t6", t)
  })
  f$release(); rm(f, P, X); par_release()
  list(rows = do.call(rbind, rows), curves = curves,
       tstar = data.frame(draw = task$label, operator = rownames(tbest),
                          t_dir1 = tbest[, 1], t_dir2 = tbest[, 2], t_sel = tsel))
}

## ---- tasks ---------------------------------------------------------------------
done <- sub("\\.csv$", "", list.files(file.path(opt$out, "draws"), "\\.csv$"))
tasks <- list()
for (b in seq_len(opt$B)) {
  sd_b <- opt$seed + b              # chips_bootstrap.R cluster seed
  set.seed(sd_b)
  tasks[[b]] <- list(label = sprintf("cluster_%03d", b), cols = draw_cols(), seed = sd_b)
}
tasks <- Filter(function(t) !t$label %in% done, tasks)
say(sprintf("%d draws done; %d to run; depth grid %s", length(done), length(tasks),
            paste(range(tgrid), collapse = "..")))
source(file.path(script_dir, "parallel.R"))
if (length(tasks)) {
  pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
  st <- foreach(task = tasks, .inorder = FALSE, .errorhandling = "pass") %dopar% {
    if (par_should_stop()) return("stopped")
    t0 <- Sys.time()
    set.seed(task$seed)
    r <- run_draw(task)
    for (k in c("tstar", "draws")) {
      d <- if (k == "draws") r$rows else merge(r$tstar, r$curves, by = c("draw", "operator"))
      tmp <- file.path(opt$out, k, paste0(task$label, ".csv.part"))
      utils::write.csv(d, tmp, row.names = FALSE)
      file.rename(tmp, file.path(opt$out, k, paste0(task$label, ".csv")))
    }
    sprintf("%s %.0f s (t* soft10 %d, pearson %d)", task$label,
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            r$tstar$t_sel[r$tstar$operator == "soft10"],
            r$tstar$t_sel[r$tstar$operator == "pearson_top10"])
  }
  par_stop(pc)
  bad <- vapply(st, function(x) inherits(x, "error"), NA)
  if (any(bad)) say("draws failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
  say("finished: ", paste(unlist(st[!bad]), collapse = "; "))
}

## ---- summary -------------------------------------------------------------------
d <- do.call(rbind, lapply(list.files(file.path(opt$out, "draws"), "\\.csv$", full.names = TRUE),
                           utils::read.csv))
tz <- do.call(rbind, lapply(list.files(file.path(opt$out, "tstar"), "\\.csv$", full.names = TRUE),
                            utils::read.csv))
utils::write.csv(d, file.path(opt$out, "draws.csv"), row.names = FALSE)
utils::write.csv(tz, file.path(opt$out, "tstar_curves.csv"), row.names = FALSE)
q <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)
ts1 <- unique(tz[, c("draw", "operator", "t_dir1", "t_dir2", "t_sel")])
for (op in unique(ts1$operator)) {
  x <- ts1$t_sel[ts1$operator == op]
  say(sprintf("t*_draw %-14s n = %d  median %g  IQR [%g, %g]  2.5-97.5%% [%g, %g]",
              op, length(x), stats::median(x), q(x, .25), q(x, .75), q(x, .025), q(x, .975)))
}
ref <- d[d$method == "clr_hg_stouffer", c("draw", "evidence", "aupr", "coherence",
                                           "tf_aupr_ratio", "tf_auroc_median")]
names(ref)[-(1:2)] <- paste0(names(ref)[-(1:2)], "_ref")
m <- merge(d, ref, by = c("draw", "evidence"))
mets <- c("aupr", "coherence", "tf_aupr_ratio", "tf_auroc_median")
summ <- do.call(rbind, lapply(split(m, list(m$evidence, m$method), drop = TRUE), function(x) {
  out <- data.frame(evidence = x$evidence[1], method = x$method[1], n = nrow(x),
                    t_median = stats::median(x$t))
  for (k in mets) {
    dd <- x[[k]] - x[[paste0(k, "_ref")]]
    out[[k]] <- mean(x[[k]]); out[[paste0("d_", k)]] <- mean(dd)
    out[[paste0("d_", k, "_lo")]] <- q(dd, .025); out[[paste0("d_", k, "_hi")]] <- q(dd, .975)
    out[[paste0("p_", k, "_le0")]] <- mean(dd <= 0)
  }
  out$tf_frac_beats_clr <- mean(x$tf_frac_beats_clr)
  out
}))
summ <- summ[order(summ$evidence, -summ$aupr), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
for (ev in c("SC", "all")) {
  s <- summ[summ$evidence == ev, ]
  say(sprintf("== nested cluster bootstrap [%s], n = %d draws (differences vs CLR within draw, 2.5-97.5%%)", ev, s$n[1]))
  for (k in seq_len(nrow(s)))
    say(sprintf("   %-20s comember AUPR %.4f %+.4f [%+.4f, %+.4f] | coherence %.3f %+.3f [%+.3f, %+.3f] | TF-node AUPR ratio %.2f %+.2f [%+.2f, %+.2f] | per-TF AUROC %.3f %+.3f [%+.3f, %+.3f], beats CLR in %.0f%% of TFs",
                s$method[k], s$aupr[k], s$d_aupr[k], s$d_aupr_lo[k], s$d_aupr_hi[k],
                s$coherence[k], s$d_coherence[k], s$d_coherence_lo[k], s$d_coherence_hi[k],
                s$tf_aupr_ratio[k], s$d_tf_aupr_ratio[k], s$d_tf_aupr_ratio_lo[k], s$d_tf_aupr_ratio_hi[k],
                s$tf_auroc_median[k], s$d_tf_auroc_median[k], s$d_tf_auroc_median_lo[k],
                s$d_tf_auroc_median_hi[k], 100 * s$tf_frac_beats_clr[k]))
}
say("done")
