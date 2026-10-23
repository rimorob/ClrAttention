#!/usr/bin/env Rscript
# CLR(log P^t) test, as preregistered in docs/preregistration_clr_logpt.md (committed 90b0e3d).
# This file implements: the gates, the a3r operator (and a2r, a1) and the K562 reproduction check.
# Also implemented: --mode selections (label-free held-out t*, see below).
# NOT implemented yet: EcoCyc truth build, E. coli / RPE1 / K562 truth phases (gate-then-stop stubs).
#
# Modes (--mode):
#   repro   K562 implementation check. Reads ONLY the K562 raw-CORUM truth rows and must reproduce the
#           exploratory AUPR 0.1054351 (a_clr_bilateral_logP_t9_robust) to 6 significant digits.
#           Writes <out>/k562_repro/repro_check.csv (passed = TRUE/FALSE). --quick N uses the first N
#           genes (seeded sample) and reads NO truth (operator diagnostics only; never counts as the check).
#   selections  Label-free held-out t* (preregistration, "Secondary: held-out t*"). Never reads a truth.
#           --dataset m3d|ecomac|rpe1|k562ess. Criterion A (AP of S for the other half's top-1% CLR pairs,
#           upper triangle, ap_fast) on the keff splits (same seeds and generators as keff_test.R: E. coli seed
#           20261014, 20 lab / author splits; human 20261015, 10 random halves), both directions; operator = keff-10
#           softmax (cap 50, self-loop for empty rows), lazy walk alpha 0.5; readouts: a3r (new) and attention
#           (comparator) over the grid; per split-direction argmax (first max), per split geometric mean of the
#           two directions (rounded), median over splits (floor(x + 0.5)). Attention t* must equal the committed
#           K10 selections (m3d 9, ecomac 8) or the run stops. Writes <out>/<dataset>/{selections.csv, ...};
#           per-split curves are cached (resumable). Then label-free checks at t = 9 and at t* on the full data.
#   truth   gate-checked stub (see gates below).
#
# Gates (functions below), each failure is a stop():
#   gate_prereg()      preregistration tracked in git and identical to HEAD; hash of last commit touching
#                      it is written into every output.
#   gate_repro()       K562 reproduction output exists, non-quick, passed, and carries the same prereg hash.
#   gate_selections()  <out>/<dataset>/selections.csv byte-identical to docs/clr_logpt_selections/<dataset>_selections.csv.
#   gate_ecoli_truth() SHA-256 of the built E. coli complex member list equals the committed one
#                      (docs/clr_logpt_selections/ecoli_complexes.sha256).
# Truth phases (E. coli, RPE1, K562 scoring) call all four gates; the repro mode calls only gate_prereg.
#
# Usage (repo root):
#   Rscript analysis/clr_logpt_test.R --mode repro [--quick 0|N] [--threads 4] [--out results/clr_logpt_test]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(mode = "repro", dataset = "k562ess", out = "results/clr_logpt_test", threads = 4L, quick = 0L,
            t = 9L, alpha = 0.5, keff = 10, cap = 50L, seed = 20261022L,
            workers = 4L, nsplit = NA_integer_, mem_gb = 4, target_frac = 0.01,
            tgrid = "1,2,3,4,5,6,8,10,12,16,20,24,32", h5 = "",
            expr = "results/baselines_k562ess_truthgenes_genie3/expression.tsv",
            members = "results/baselines_k562ess_truthgenes_genie3/regulon_members.csv")
ints <- c("threads", "quick", "t", "cap", "seed", "workers", "nsplit")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("alpha", "keff", "mem_gb", "target_frac")) as.numeric(v) else v
  i <- i + 2L
}
PREREG <- "docs/preregistration_clr_logpt.md"
SELDIR <- "docs/clr_logpt_selections"
REPRO_AUPR <- 0.1054351          # exploratory a_clr_bilateral_logP_t9_robust, K562 raw CORUM, t = 9
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n") }

## ---- gates --------------------------------------------------------------------------------------
git <- function(...) system2("git", c(...), stdout = TRUE, stderr = TRUE)
gate_prereg <- function() {
  if (system2("git", c("ls-files", "--error-unmatch", PREREG), stdout = FALSE, stderr = FALSE) != 0L)
    stop("gate refused: ", PREREG, " is not tracked in git")
  if (system2("git", c("diff", "--quiet", "HEAD", "--", PREREG)) != 0L)
    stop("gate refused: ", PREREG, " differs from HEAD (commit the change and log it in the change log first)")
  h <- git("log", "-1", "--format=%H", "--", PREREG)
  if (!length(h) || !nzchar(h[1])) stop("gate refused: no commit touches ", PREREG)
  h[1]
}
rd <- function(f) readBin(f, "raw", file.info(f)$size)
gate_repro <- function(hash) {
  f <- file.path(opt$out, "k562_repro", "repro_check.csv")
  if (!file.exists(f)) stop("truth phase refused: K562 reproduction check output ", f, " does not exist")
  r <- utils::read.csv(f)
  if (nrow(r) != 1L || !isTRUE(as.logical(r$passed)) || isTRUE(as.logical(r$quick)))
    stop("truth phase refused: the K562 reproduction check did not pass (or was a --quick run)")
  if (r$prereg_commit != hash) stop("truth phase refused: reproduction check was run under another preregistration commit")
  invisible(TRUE)
}
gate_selections <- function(dataset) {
  ref <- file.path(SELDIR, paste0(dataset, "_selections.csv")); cur <- file.path(opt$out, dataset, "selections.csv")
  if (!file.exists(ref) || !file.exists(cur)) stop("truth phase refused: missing ", ref, " or ", cur)
  if (system2("git", c("ls-files", "--error-unmatch", ref), stdout = FALSE, stderr = FALSE) != 0L)
    stop("truth phase refused: ", ref, " is not committed")
  if (!identical(rd(cur), rd(ref))) stop("truth phase refused: ", cur, " is not byte-identical to ", ref)
  invisible(TRUE)
}
gate_ecoli_truth <- function(member_file) {
  ref <- file.path(SELDIR, "ecoli_complexes.sha256")
  if (!file.exists(ref) || !file.exists(member_file)) stop("truth phase refused: missing ", ref, " or ", member_file)
  want <- strsplit(readLines(ref, 1L), "[[:space:]]+")[[1]][1]
  got <- digest::digest(file = member_file, algo = "sha256")
  if (!identical(tolower(want), got)) stop("truth phase refused: E. coli complex member list SHA-256 ", got, " != committed ", want)
  invisible(TRUE)
}

## ---- operator: a3r (primary), a2r, a1 -----------------------------------------------------------
# P^t, dense double, no floor (preregistration: "Pt <- P; for (s in 2:t) Pt <- Pt %*% P").
lazy_walk_P <- function(A, alpha = 0.5) { P <- alpha * as.matrix(A); diag(P) <- diag(P) + 1 - alpha; P }
power_t <- function(P, t) { Pt <- P; if (t >= 2L) for (s in 2:t) Pt <- Pt %*% P; Pt }
# Steps 2-3: per-row mean/SD (n - 1) over non-missing off-diagonal entries; negative, missing, non-finite z -> 0.
rowz_na <- function(L) {
  diag(L) <- NA
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  sdv <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1))
  Z <- (L - mu) / sdv
  Z[!is.finite(Z) | Z < 0] <- 0     # includes rows with < 2 reachable entries or sd = 0
  Z
}
stouffer_sym <- function(Z) { S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }
# Step 1: L = log P^t, zeros (unreachable or underflowed) and diagonal missing.
log_na <- function(Pt) { L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; L }
a3r <- function(Pt) stouffer_sym(rowz_na(log_na(Pt)))                      # primary readout R
a2r <- function(Pt) { L <- log_na(Pt); stouffer_sym(rowz_na((L + t(L)) / 2)) }  # secondary (NA if either NA)
a1  <- function(Pt) stouffer_sym(rowz_na(Pt))                              # mechanism control (raw P^t)
attention <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }       # comparator
# Label-free checks: unreachable share, rows with < 2 reachable entries, sigma_i distribution.
operator_checks <- function(Pt) {
  L <- log_na(Pt); G <- nrow(L); n <- rowSums(!is.na(L))
  mu <- rowMeans(L, na.rm = TRUE); sdv <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1))
  q <- stats::quantile(sdv[is.finite(sdv)], c(0, .05, .5, .95, 1))
  data.frame(G = G, unreachable_share = sum(is.na(L) & row(L) != col(L)) / (as.numeric(G) * (G - 1)),
             rows_lt2_reachable = sum(n < 2L), sigma_min = q[[1]], sigma_p05 = q[[2]], sigma_med = q[[3]],
             sigma_p95 = q[[4]], sigma_max = q[[5]])
}
sig6 <- function(x) signif(x, 6)

## ---- common start: gate, output dir ------------------------------------------------------------
HASH <- gate_prereg()
say("preregistration gate passed; last commit touching it: ", HASH)
suppressPackageStartupMessages({ library(clr); library(Matrix) })
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
source(file.path(script_dir, "regulondb.R")); source(file.path(script_dir, "regulons.R"))

## ---- mode: selections (label-free held-out t*; reads NO truth) ------------------------------------
if (opt$mode == "selections") {
  ds <- opt$dataset
  if (!ds %in% c("m3d", "ecomac", "rpe1", "k562ess")) stop("selections: dataset must be m3d, ecomac, rpe1 or k562ess")
  ecoli <- ds %in% c("m3d", "ecomac")
  if (is.na(opt$nsplit)) opt$nsplit <- if (ecoli) 20L else 10L
  seed <- if (ecoli) 20261014L else 20261015L                  # keff_test.R seeds
  tgrid <- as.integer(strsplit(opt$tgrid, ",")[[1]])
  expected_att <- c(m3d = 9L, ecomac = 8L)                     # committed K10 selections (docs/keff_selections)
  dout <- file.path(opt$out, ds); dir.create(file.path(dout, "splits"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dout, "cache"), showWarnings = FALSE)
  logf <- file.path(dout, "selections.log")
  say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
  wcsv <- function(x, n) utils::write.csv(x, file.path(dout, n), row.names = FALSE)
  source(file.path(script_dir, "parallel.R"))
  say("selections: ", ds, "; options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))
  import_defs <- function(file, nms) {
    got <- character(0)
    for (e in parse(file.path(script_dir, file), keep.source = FALSE))
      if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
          as.character(e[[2]]) %in% nms && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
        eval(e, globalenv()); got <- c(got, as.character(e[[2]])) }
    if (length(setdiff(nms, got))) stop("not found in ", file, ": ", paste(setdiff(nms, got), collapse = ", "))
  }
  import_defs("lab_splits.R", c("lab_split", "build_fit"))
  import_defs("ecomac_batch.R", "author_split")
  THREADS <- opt$threads; nthr <- function() THREADS

  ## data: as keff_test.R (m3d, ecomac, k562ess) and perturbseq_bulk.R (rpe1); expression only
  if (ds == "m3d") {
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
    raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
    desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
    exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description)
    exper[!grepl("experimenter:", desc$description)] <- NA
    lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])
    lab <- lab[colnames(X)]; lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab))); rm(desc, exper)
  } else if (ds == "ecomac") {
    X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
    rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
    plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
    author <- utils::read.delim("data/ecomac/samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)[["Author"]]
    stopifnot(length(author) == ncol(X), length(plat) == ncol(X))
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
  } else {
    h5 <- if (nzchar(opt$h5)) opt$h5 else if (ds == "rpe1") "data/perturbseq/rpe1_normalized_bulk_01.h5ad" else
      "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"
    set.seed(20261015L)                                  # perturbseq_bulk.R split of significant knockdowns
    rdh <- function(p) rhdf5::h5read(h5, p)
    cat_of <- function(grp, field) { x <- rdh(paste0(grp, "/", field))
      if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
      cats <- tryCatch(rdh(paste0(grp, "/__categories/", field)), error = function(e) NULL)
      if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
    Z <- rdh("X"); rownames(Z) <- cat_of("/var", "gene_name"); pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
    pval <- as.numeric(rdh("/obs/energy_test_p_value"))
    K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
    Ahalf <- sort(K[seq_len(floor(length(K) / 2))])
    X <- Z[, Ahalf]; rm(Z); invisible(gc())
    nf <- !is.finite(X); bad_gene <- rowMeans(nf) > 0.01; X <- X[!bad_gene, , drop = FALSE]
    n_imp <- sum(!is.finite(X)); X[!is.finite(X)] <- 0; rm(nf)
    say(sprintf("non-finite: %d genes dropped (> 1%% missing), %d values set to 0", sum(bad_gene), n_imp))
    if (ds == "k562ess") {
      X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
      dup <- rownames(X) %in% rownames(X)[duplicated(rownames(X))]
      if (any(dup)) { say(sprintf("dropped %d rows with duplicated gene names", sum(dup))); X <- X[!dup, ] }
    }
    rm(rdh, cat_of, pert, pval, K, Ahalf)
  }
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  if (opt$quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
  storage.mode(X) <- "double"
  genes <- rownames(X); G <- length(genes)
  say(sprintf("%s: %d genes x %d samples%s", ds, G, ncol(X), if (opt$quick > 0L) " (--quick)" else ""))

  ## splits (keff_test.R generators and seeds)
  set.seed(seed)
  splits <- lapply(seq_len(opt$nsplit), function(s) switch(ds, m3d = lab_split(), ecomac = author_split(),
                                                           sample(rep(1:2, length.out = ncol(X)))))
  fp <- list(dataset = ds, quick = opt$quick, seed = seed, G = G, genes_md5 = digest::digest(genes, algo = "md5"),
             splits_md5 = digest::digest(splits, algo = "md5"), tgrid = tgrid, cap = opt$cap, keff = opt$keff,
             alpha = opt$alpha, target_frac = opt$target_frac, prereg_commit = HASH)
  fp_f <- file.path(dout, "splits", "fingerprint.rds")
  if (file.exists(fp_f)) { if (!isTRUE(all.equal(readRDS(fp_f), fp))) stop("splits/ fingerprint mismatch (", fp_f, "): use a new --out")
  } else { if (length(list.files(file.path(dout, "splits"), "\\.csv$"))) stop("splits/ has results without a fingerprint")
    saveRDS(fp, fp_f) }
  wcsv(data.frame(sample = colnames(X), do.call(cbind, stats::setNames(splits, sprintf("split_%02d", seq_along(splits))))), "split_assignments.csv")
  say(sprintf("%d splits; half sizes (first split): %d / %d", opt$nsplit, sum(splits[[1]] == 1L), sum(splits[[1]] == 2L)))

  selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
  run_split <- function(si) {
    h <- splits[[si]]
    Cs <- lapply(1:2, function(d) { set.seed(seed + 100L * si + d); f <- build_fit(X[, h == d, drop = FALSE])
      C <- f$clr_scores; f$release(); C })
    ut <- which(upper.tri(matrix(0, G, G))); rows <- list()
    for (dir in 1:2) {
      te <- 3L - dir
      sv <- Cs[[te]][ut]; tA <- as.integer(sv >= stats::quantile(sv, 1 - opt$target_frac)); rm(sv)
      P <- lazy_walk_P(selfloop(clr:::.softmax_rows(Cs[[dir]], keff = opt$keff, cap = opt$cap)), opt$alpha)
      cr <- matrix(NA_real_, length(tgrid), 2, dimnames = list(NULL, c("attention", "a3r"))); Pt <- P
      for (t in seq_len(max(tgrid))) {
        if (t > 1L) Pt <- as.matrix(Pt %*% P)
        if (t %in% tgrid) { k <- match(t, tgrid)
          S <- attention(Pt); cr[k, 1] <- ap_fast(S[ut], tA); rm(S)
          S <- a3r(Pt); cr[k, 2] <- ap_fast(S[ut], tA); rm(S); invisible(gc(verbose = FALSE)) }
      }
      for (r in colnames(cr)) rows[[length(rows) + 1L]] <- data.frame(split = si, direction = dir, readout = r, t = tgrid, crit_A = cr[, r])
    }
    do.call(rbind, rows)
  }
  done <- sub("\\.csv$", "", list.files(file.path(dout, "splits"), "\\.csv$"))
  todo <- setdiff(sprintf("split_%02d", seq_len(opt$nsplit)), done)
  say(sprintf("%d splits done, %d to run; t grid %s", length(done), length(todo), opt$tgrid))
  if (length(todo)) {
    W <- min(opt$workers, length(todo)); THREADS <- opt$threads
    pc <- par_start(mem_gb = opt$mem_gb, workers = W, reserve = max(0L, .par_cores() - W * opt$threads), say = say,
                    export = setdiff(ls(globalenv()), c("pc")))
    st <- foreach(nm = todo, .inorder = FALSE, .errorhandling = "pass") %dopar% {
      if (par_should_stop()) return("stopped")
      si <- as.integer(sub("split_", "", nm)); t0 <- Sys.time()
      d <- run_split(si)
      tmp <- file.path(dout, "splits", paste0(nm, ".csv.part"))
      utils::write.csv(d, tmp, row.names = FALSE); file.rename(tmp, file.path(dout, "splits", paste0(nm, ".csv")))
      par_release()
      sprintf("%s %.1f min", nm, as.numeric(difftime(Sys.time(), t0, units = "mins")))
    }
    par_stop(pc)
    bad <- vapply(st, function(x) inherits(x, "error"), NA)
    if (any(bad)) say("splits failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
    say("split timings: ", paste(unlist(st[!bad]), collapse = "; "))
    if (any(bad) || any(unlist(st[!bad]) == "stopped")) stop("not all splits finished; rerun to resume")
  }
  THREADS <- opt$threads
  cur <- do.call(rbind, lapply(sprintf("%s/splits/split_%02d.csv", dout, seq_len(opt$nsplit)), utils::read.csv))
  wcsv(cur, "heldout_curves.csv")

  ## aggregation as keff_test.R select(): argmax per split-direction (first max), geometric mean of the two
  ## directions per split (rounded), median over splits (floor(x + 0.5))
  sels <- list(); per <- list()
  for (r in c("a3r", "attention")) {
    x <- cur[cur$readout == r, ]
    per_sd <- do.call(rbind, lapply(split(x, list(x$split, x$direction), drop = TRUE), function(y) {
      y <- y[order(y$t), ]; k <- which.max(y$crit_A)
      data.frame(split = y$split[1], direction = y$direction[1], t_best = y$t[k], crit_best = y$crit_A[k]) }))
    per_s <- do.call(rbind, lapply(split(per_sd, per_sd$split), function(y)
      data.frame(split = y$split[1], t_split = as.integer(round(exp(mean(log(y$t_best))))))))
    ts <- as.integer(floor(stats::median(per_s$t_split) + 0.5))
    sels[[r]] <- data.frame(dataset = ds, readout = r, criterion = "crit_A", t_star = ts,
      t_q25 = stats::quantile(per_s$t_split, 0.25, names = FALSE), t_q75 = stats::quantile(per_s$t_split, 0.75, names = FALSE),
      t_min = min(per_s$t_split), t_max = max(per_s$t_split), n_splits = nrow(per_s), t_secondary = opt$t,
      flag_t_lt4 = ts < 4L, quick = opt$quick, prereg_commit = HASH, stringsAsFactors = FALSE)
    per[[r]] <- cbind(readout = r, per_sd)
  }
  sel <- do.call(rbind, sels); wcsv(do.call(rbind, per), "selection_per_split_direction.csv")
  say(sprintf("t*: a3r %d (IQR %.2f-%.2f, range %d-%d)%s; attention %d (IQR %.2f-%.2f, range %d-%d)",
              sel$t_star[1], sel$t_q25[1], sel$t_q75[1], sel$t_min[1], sel$t_max[1], if (sel$flag_t_lt4[1]) " FLAG t* < 4" else "",
              sel$t_star[2], sel$t_q25[2], sel$t_q75[2], sel$t_min[2], sel$t_max[2]))
  if (opt$quick == 0L && ds %in% names(expected_att) && sel$t_star[sel$readout == "attention"] != expected_att[[ds]])
    stop("STOP: attention t* ", sel$t_star[sel$readout == "attention"], " != committed K10 selection ", expected_att[[ds]],
         " for ", ds, "; selections.csv NOT written")
  wcsv(sel, "selections.csv")

  ## label-free checks on the full data at t = 9 and at the new readout's t*
  Cf <- file.path(dout, "cache", "clr_full.rds")
  if (file.exists(Cf)) { cf <- readRDS(Cf); stopifnot(identical(cf$genes_md5, fp$genes_md5)); C <- cf$C; rm(cf)
  } else { set.seed(seed); f <- build_fit(X); C <- f$clr_scores; f$release(); rm(f)
    saveRDS(list(C = C, genes_md5 = fp$genes_md5), Cf, compress = FALSE) }
  P <- lazy_walk_P(selfloop(clr:::.softmax_rows(C, keff = opt$keff, cap = opt$cap)), opt$alpha); rm(C)
  tt <- sort(unique(c(opt$t, sel$t_star[1]))); Pt <- P; chk <- list()
  for (t in seq_len(max(tt))) { if (t > 1L) Pt <- as.matrix(Pt %*% P)
    if (t %in% tt) chk[[as.character(t)]] <- cbind(dataset = ds, depth = if (t == opt$t && t == sel$t_star[1]) "t9=t*" else if (t == opt$t) "t9" else "t_star_a3r",
                                                  t = t, operator_checks(Pt)) }
  chk <- do.call(rbind, chk); chk$prereg_commit <- HASH; wcsv(chk, "label_free_checks.csv")
  say("label-free checks:\n", paste(utils::capture.output(print(chk[, 1:8], row.names = FALSE)), collapse = "\n"))
  say("done")
  quit(save = "no")
}

## ---- mode: selections / truth stubs -------------------------------------------------------------
if (opt$mode == "truth") {
  {
    gate_repro(HASH)
    if (opt$dataset %in% c("m3d", "ecomac")) {
      gate_selections(opt$dataset)
      gate_ecoli_truth(file.path(SELDIR, "ecoli_complexes_members.tsv"))
    } else if (opt$dataset == "rpe1") gate_selections("rpe1")
  }
  stop("mode '", opt$mode, "' is not implemented yet (preregistered steps 3-5 come after the reproduction check)")
}
if (opt$mode != "repro") stop("unknown mode: ", opt$mode)

## ---- mode: repro ---------------------------------------------------------------------------------
if (opt$dataset != "k562ess") stop("repro mode is K562 only")
dout <- file.path(opt$out, "k562_repro"); dir.create(dout, recursive = TRUE, showWarnings = FALSE)
T0 <- Sys.time()
ex <- utils::read.delim(opt$expr, check.names = FALSE)
genes <- ex[[1]]; X <- as.matrix(ex[, -1]); rownames(X) <- genes; rm(ex)
if (opt$quick > 0L) { set.seed(opt$seed); keep <- sort(sample(nrow(X), min(opt$quick, nrow(X)))); X <- X[keep, , drop = FALSE]; genes <- genes[keep] }
G <- length(genes); say(sprintf("K562 half A: %d genes x %d profiles%s", G, ncol(X), if (opt$quick > 0L) " (--quick)" else ""))
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$
  calibrate(method = "normal", combine = "stouffer")
f$build_operator(alpha = opt$alpha, softmax_keff = opt$keff, softmax_cap = opt$cap)
A <- f$operator; dimnames(A) <- NULL; f$release(); rm(f, X); invisible(gc())
say(sprintf("MI, CLR, softmax operator in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
Pt <- power_t(lazy_walk_P(A, opt$alpha), opt$t); rm(A)
chk <- operator_checks(Pt); chk$prereg_commit <- HASH; chk$quick <- opt$quick > 0L; chk$t <- opt$t
utils::write.csv(chk, file.path(dout, if (opt$quick > 0L) "label_free_quick.csv" else "label_free_checks.csv"), row.names = FALSE)
say(sprintf("t = %d: unreachable share %.4f, rows with <2 reachable %d, sigma median %.3f",
            opt$t, chk$unreachable_share, chk$rows_lt2_reachable, chk$sigma_med))
S <- a3r(Pt); rm(Pt)
stopifnot(isSymmetric(S), all(diag(S) == 0), all(S >= 0))
if (opt$quick > 0L) { say("--quick: operator ran; no truth read; no check output written"); quit(save = "no") }

## the ONLY truth read: K562 raw CORUM rows (other truth rows in the file are dropped immediately)
mem <- utils::read.csv(opt$members, stringsAsFactors = FALSE); mem <- mem[mem$truth == "CORUM", ]
m <- lapply(strsplit(mem$members, ";", fixed = TRUE), function(g) sort(match(g, genes))); names(m) <- mem$regulon
reg <- list(members = m, class = stats::setNames(rep("complex", length(m)), names(m)))
pr <- comembership_pairs(reg, stats::setNames(character(0), character(0)), genes)
say(sprintf("K562 raw CORUM: %d complexes, %d pairs, %d positives", length(m), length(pr$label), sum(pr$label)))
sc <- S[cbind(pr$i, pr$j)]
aupr <- ap_fast(sc, pr$label)
ok <- identical(sig6(aupr), sig6(REPRO_AUPR))
out <- data.frame(dataset = "k562ess", truth = "CORUM_raw_pooled", readout = "a3r", t = opt$t, aupr = aupr,
                  aupr_6sig = sig6(aupr), expected_6sig = sig6(REPRO_AUPR), abs_diff = abs(aupr - REPRO_AUPR),
                  n_pairs = length(pr$label), n_pos = sum(pr$label), passed = ok, quick = FALSE,
                  prereg_commit = HASH, git_head = git("rev-parse", "HEAD")[1],
                  time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
utils::write.csv(out, file.path(dout, "repro_check.csv"), row.names = FALSE)
say(sprintf("AUPR %.9f vs expected %.7f: %s", aupr, REPRO_AUPR, if (ok) "PASSED" else "FAILED"))
say(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
if (!ok) quit(save = "no", status = 2L)
