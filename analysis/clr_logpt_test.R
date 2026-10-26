#!/usr/bin/env Rscript
# CLR(log P^t) test, as preregistered in docs/preregistration_clr_logpt.md (committed 90b0e3d).
# This file implements: the gates, the a3r operator (and a2r, a1), the K562 reproduction check, --mode selections
# (label-free held-out t*), and the confirmatory scoring: --mode truth (per dataset), --mode boot (M3D cluster
# bootstrap) and --mode verdict (decision rules). Implemented, NOT yet reviewed or run on any E. coli / RPE1 truth.
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
#   truth   --dataset m3d|ecomac|rpe1|k562ess. Full-data scores (a3r, a2r, a1, attention at t = 9; a3r and attention at
#           the committed t*; CLR; |r|; hub-stratified versions), delete-half gene jackknife (100 replicates, seed
#           20261022, drawn per truth in table order), one-sided z p-values, endpoint rows (P1-P4, S1-S5), per-complex
#           AP, r-protein module coherence. Cached per stage (fit, pair scores, each truth's jackknife): resumable.
#           Writes <out>/truth/<dataset>/ and a DONE marker. Refuses a second run once DONE.
#   boot    --dataset m3d. 200-draw paired cluster bootstrap (chips_bootstrap.R cluster mode, seeds 20260924 + b), MI, CLR,
#           attention and a3r refitted per draw at t = 9; P1 and P3 differences; per-draw files (resumable).
#   verdict Reads the endpoint tables and bootstrap intervals and applies the preregistered rules (intersection-union
#           adoption, Holm over P1-P4 for single claims, Holm over S1-S5 only if adopted); writes <out>/verdict/verdict.csv.
#   --synthetic 1  smoke test of truth/boot/verdict on generated data and a generated truth (never confirmatory;
#           bypasses the repro/selection/E. coli-truth/script-clean gates, keeps the preregistration gate; use a scratch --out).
#
# Gates (functions below), each failure is a stop():
#   gate_prereg()      preregistration tracked in git and identical to HEAD; hash of last commit touching
#                      it is written into every output.
#   gate_repro()       K562 reproduction output exists, non-quick, passed, and carries the same prereg hash.
#   gate_selections()  <out>/<dataset>/selections.csv byte-identical to docs/clr_logpt_selections/<dataset>_selections.csv.
#   gate_ecoli_truth() SHA-256 of the built E. coli complex member list equals the committed one
#                      (docs/clr_logpt_selections/ecoli_complexes.sha256).
# Truth phases (E. coli, RPE1, K562 scoring) call all four gates plus gate_script_clean (the script is committed) and
# gate_order (E. coli, then RPE1, then K562); the repro mode calls only gate_prereg.
#
# Usage (repo root):
#   Rscript analysis/clr_logpt_test.R --mode repro [--quick 0|N] [--threads 4] [--out results/clr_logpt_test]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(mode = "repro", dataset = "k562ess", out = "results/clr_logpt_test", threads = 4L, quick = 0L,
            t = 9L, alpha = 0.5, keff = 10, cap = 50L, seed = 20261022L,
            workers = 4L, nsplit = NA_integer_, mem_gb = 4, target_frac = 0.01,
            synthetic = 0L, draws = 200L, jack = 100L, hubstrat = 1L,
            tgrid = "1,2,3,4,5,6,8,10,12,16,20,24,32", h5 = "",
            expr = "results/baselines_k562ess_truthgenes_genie3/expression.tsv",
            members = "results/baselines_k562ess_truthgenes_genie3/regulon_members.csv")
ints <- c("threads", "quick", "t", "cap", "seed", "workers", "nsplit", "synthetic", "draws", "jack", "hubstrat")
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

## ======================================================================================================
## ---- truth phase, bootstrap and verdict (modes truth / boot / verdict) --------------------------------
## ======================================================================================================
## Preregistered constants (docs/preregistration_clr_logpt.md, "Statistics" and "Decision rules"). They are
## fixed here, not options: the truth modes refuse other values (outside --synthetic smoke tests).
SYN <- opt$synthetic == 1L                  # smoke test on synthetic data and a synthetic truth; never confirmatory
NI_MARGIN <- 0.005; ALPHA1 <- 0.025; EFFECT_MIN <- 0.005
JACK_SEED <- 20261022L; JACK_N <- 100L; BOOT_N <- 200L; BOOT_SEED0 <- 20260924L; T_PRIMARY <- 9L
RP_RE <- "^RP[LS][0-9]+[A-Z]?[0-9]*$|^RPLP[0-2]$|^RPSA$"
DATASETS <- c("m3d", "ecomac", "rpe1", "k562ess")
GIT_HEAD <- git("rev-parse", "HEAD")[1]
SCRIPT_COMMIT <- { h <- git("log", "-1", "--format=%H", "--", "analysis/clr_logpt_test.R"); if (length(h)) h[1] else NA_character_ }

fix_check <- function() {
  if (SYN) return(invisible(TRUE))
  if (opt$quick != 0L || opt$jack != JACK_N || opt$draws != BOOT_N || opt$t != T_PRIMARY || opt$seed != JACK_SEED ||
      opt$alpha != 0.5 || opt$keff != 10 || opt$cap != 50L)
    stop("truth/boot/verdict modes refuse non-preregistered settings (--quick, --jack, --draws, --t, --seed, --alpha, --keff, --cap)")
}
# Extra gate (not in the preregistration, protects its audit trail): the script that produced an output is the committed one.
gate_script_clean <- function() {
  f <- "analysis/clr_logpt_test.R"
  if (system2("git", c("ls-files", "--error-unmatch", f), stdout = FALSE, stderr = FALSE) != 0L ||
      system2("git", c("diff", "--quiet", "HEAD", "--", f)) != 0L)
    stop("gate refused: ", f, " is not committed (commit the code before any truth-phase run)")
}
# E. coli first, then RPE1, then K562 (preregistration, protocol step 5).
gate_order <- function(ds) {
  need <- switch(ds, rpe1 = c("m3d", "ecomac"), k562ess = c("m3d", "ecomac", "rpe1"), character(0))
  for (n in need) if (!file.exists(file.path(opt$out, "truth", n, "DONE")))
    stop("truth phase refused: run order is E. coli, then RPE1, then K562; ", n, " is not finished")
}
gate_all <- function(ds) {
  fix_check()
  if (SYN) return(invisible(TRUE))
  gate_script_clean(); gate_repro(HASH)
  if (ds %in% c("m3d", "ecomac")) { gate_selections(ds); gate_ecoli_truth(file.path(SELDIR, "ecoli_complexes_members.tsv")) }
  else if (ds == "rpe1") gate_selections("rpe1")
  else if (file.exists(file.path(SELDIR, "k562ess_selections.csv"))) gate_selections("k562ess")
  gate_order(ds)
  invisible(TRUE)
}
stamp <- function(df) { df$prereg_commit <- HASH; df$git_head <- GIT_HEAD; df$script_commit <- SCRIPT_COMMIT; df$synthetic <- SYN; df }
wout <- function(x, f) {                      # every output carries the preregistration commit hash (protocol step 1)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".part"); utils::write.csv(stamp(x), tmp, row.names = FALSE); file.rename(tmp, f)
}
logger <- function(f) function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n")
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); cat(m, "\n", file = f, append = TRUE) }
none_op <- stats::setNames(character(0), character(0))
p_sup <- function(d, se) if (is.finite(se) && se > 0) 1 - stats::pnorm(d / se) else NA_real_
p_ni  <- function(d, se, m = NI_MARGIN) if (is.finite(se) && se > 0) 1 - stats::pnorm((d + m) / se) else NA_real_
p_loss <- function(d, se) if (is.finite(se) && se > 0) stats::pnorm(d / se) else NA_real_

# Overlap merge (preregistration "Overlap merge"; the rule of analysis/build_ecoli_complex_truth.R, verbatim):
# (1) drop exact duplicates and proper subsets; (2) complete-linkage hclust on 1 - Jaccard, cutree(h = 0.5), union.
merge_sets <- function(s) {
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]
  key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
    length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]
  if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
  out
}
# Hub-stratified score: within-stratum percentile in 10 strata of pi_i + pi_j (baselines.R, verbatim).
hub_stratify <- function(S, H, k = 10L) {
  S[!is.finite(S)] <- 0; ut <- which(upper.tri(S)); h <- H[ut]
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  r <- stats::ave(S[ut], g, FUN = function(x) rank(x) / length(x))
  R <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S)); R[ut] <- r; R + t(R)
}

## ---- synthetic data (smoke tests only; needs no truth) --------------------------------------------------
synth_data <- function(ds, G = 160L, n = 48L) {
  set.seed(match(ds, DATASETS) * 101L)
  ecoli <- ds %in% c("m3d", "ecomac")
  gn <- if (ecoli) sprintf("b%04d", seq_len(G)) else c(sprintf("RPL%d", 1:12), sprintf("GENE%03d", seq_len(G - 12L)))
  mod <- c(rep(1:20, each = 6L), rep(0L, G - 120L))
  Fm <- matrix(stats::rnorm(20 * n), 20, n)
  X <- matrix(stats::rnorm(G * n), G, n); X[mod > 0, ] <- X[mod > 0, ] + 1.2 * Fm[mod[mod > 0], ]
  rownames(X) <- gn; colnames(X) <- sprintf("s%03d", seq_len(n))
  cx <- lapply(1:20, function(k) gn[mod == k])
  cx[[21]] <- gn[13:17]                       # proper subset of module 3 (absorption)
  cx[[22]] <- cx[[4]]                         # exact duplicate of module 4
  cx[[23]] <- c(gn[19:23], gn[121])           # Jaccard 5/7 with module 4: merged with it
  names(cx) <- sprintf("CX%02d", seq_along(cx))
  regs <- lapply(1:10, function(k) gn[mod %in% c(2 * k - 1, 2 * k)]); names(regs) <- sprintf("REG%02d", 1:10)
  op <- if (ecoli) stats::setNames(sprintf("op%03d", rep(seq_len(G / 2), each = 2L)), gn) else none_op
  list(X = X, genes = gn, complexes = cx, regs = regs, operon_of = op)
}
synth_chips <- function(sd) {                 # chips of 3 per experiment for the cluster bootstrap smoke test
  set.seed(7L); X <- sd$X; ch <- lapply(seq_len(ncol(X)), function(j) X[, j] + matrix(stats::rnorm(nrow(X) * 3L, sd = 0.3), nrow(X)))
  Xc <- do.call(cbind, ch); rownames(Xc) <- rownames(X); colnames(Xc) <- sprintf("chip%03d", seq_len(ncol(Xc)))
  list(Xc = Xc, Xa = X, exp_of = stats::setNames(rep(colnames(X), each = 3L), colnames(Xc)))
}

## ---- expression input (the same construction as the committed selections) -------------------------------
h5_cat <- function(h5, grp, field) {          # anndata categorical, old or new layout (as keff_test.R / baselines.R)
  x <- rhdf5::h5read(h5, paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rhdf5::h5read(h5, paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L]
}
half_split <- function(h5) {                  # perturbseq_bulk.R / baselines.R split of the significant knockdowns
  set.seed(20261015L)
  pert <- h5_cat(h5, "/obs", "gene_transcript"); pval <- as.numeric(rhdf5::h5read(h5, "/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A <- sort(K[seq_len(floor(length(K) / 2))]); list(A = A, B = setdiff(K, A), pert = pert)
}
load_expr <- function(ds) {
  if (SYN) { s <- synth_data(ds); return(list(X = s$X, h5 = NULL, B = NULL, synth = s)) }
  if (ds == "m3d") {
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
    raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
  } else if (ds == "ecomac") {
    X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
    rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
    plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
  } else if (ds == "rpe1") {
    h5 <- "data/perturbseq/rpe1_normalized_bulk_01.h5ad"; hs <- half_split(h5)
    Z <- rhdf5::h5read(h5, "X"); rownames(Z) <- h5_cat(h5, "/var", "gene_name"); colnames(Z) <- hs$pert; X <- Z[, hs$A]; rm(Z); invisible(gc())
    nf <- !is.finite(X); bad <- rowMeans(nf) > 0.01; X <- X[!bad, , drop = FALSE]; X[!is.finite(X)] <- 0; rm(nf)
    return(list(X = {X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; storage.mode(X) <- "double"; X}, h5 = h5, B = hs$B, synth = NULL))
  } else {                                    # k562ess: network half A, truth-gene subset (baselines.R --gene_subset truths), as repro
    ex <- utils::read.delim(opt$expr, check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex)
    h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"
    return(list(X = X, h5 = h5, B = half_split(h5)$B, synth = NULL))
  }
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; storage.mode(X) <- "double"
  list(X = X, h5 = NULL, B = NULL, synth = NULL)
}

## ---- truths ---------------------------------------------------------------------------------------------
# Returns the ordered list of truths (endpoint-table order first: complexes, regulons SC; then descriptive ones).
# Each: reg (members = gene row indices), pairs (comembership_pairs), operon, role.
mk_truth <- function(r, o, genes, role) {
  if (!length(r$members)) return(NULL)
  p <- comembership_pairs(r, o, genes); if (sum(p$label) == 0) return(NULL)
  list(reg = r, pairs = p, operon = o, role = role)
}
as_reg <- function(sets, cls) list(members = sets, class = stats::setNames(rep(cls, length(sets)), names(sets)))
ecoli_regtab <- local({ cache <- NULL; function() {
  if (is.null(cache)) { to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
    cache <<- list(rtab = build_regulon_table("data/RegulonDBExtract", to_bn), operon_of = build_operon_map("data/RegulonDBExtract", to_bn)) }
  cache } })
ecoli_members <- function(ds, genes, synth) {  # long data frame: variant, dataset, complex, bnumber (the SHA-256-gated file)
  if (!is.null(synth)) {
    cx <- lapply(synth$complexes, function(g) g)
    mk <- function(s, v) do.call(rbind, lapply(names(s), function(k) data.frame(variant = v, dataset = ds, complex = k, bnumber = s[[k]])))
    return(rbind(mk(merge_sets(cx), "main"), mk(cx, "nomerge")))
  }
  utils::read.delim(file.path(SELDIR, "ecoli_complexes_members.tsv"), stringsAsFactors = FALSE)
}
complexes_with_pos <- function(sets, operon_of, genes, with_operon) {   # as the builder's evaluability count
  sum(vapply(sets, function(m) { o <- operon_of[genes[m]]
    if (!with_operon) length(m) >= 2 && (any(is.na(o)) || length(unique(o)) > 1) else length(m) >= 2 }, TRUE))
}
build_truths_ecoli <- function(ds, genes, synth = NULL, only = NULL) {
  if (is.null(synth)) { rt <- ecoli_regtab(); rtab <- rt$rtab; operon_of <- rt$operon_of } else operon_of <- synth$operon_of
  mem <- ecoli_members(ds, genes, synth); mem <- mem[mem$dataset == ds, ]
  miss <- setdiff(mem$bnumber, genes)
  if (length(miss) && is.null(synth)) stop("E. coli member list has ", length(miss), " genes outside the ", ds, " gene universe (e.g. ", paste(utils::head(miss, 3), collapse = ", "), ")")
  sets_of <- function(v) { x <- mem[mem$variant == v, ]
    lapply(split(x$bnumber, x$complex), function(g) sort(stats::na.omit(match(g, genes)))) }
  cx <- function(v, o, role) { s <- sets_of(v); if (!length(s)) NULL else mk_truth(as_reg(s, "complex"), o, genes, role) }
  main <- sets_of("main")
  p_main <- if (length(main)) comembership_pairs(as_reg(main, "complex"), operon_of, genes) else NULL
  n_pos <- if (is.null(p_main)) 0L else sum(p_main$label)
  n_cx <- complexes_with_pos(main, operon_of, genes, FALSE)
  fallback <- !(n_pos >= 200 && n_cx >= 10)    # evaluability, structure only (preregistration "Evaluability")
  if (!is.null(synth)) fallback <- FALSE
  reg_sc <- if (is.null(synth)) make_regulons(rtab, genes, conf_keep = c("C", "S")) else
    { s <- lapply(synth$regs, function(g) match(g, genes)); list(members = s, class = stats::setNames(rep("TF", length(s)), names(s))) }
  reg_all <- if (is.null(synth)) make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?")) else
    { s <- c(lapply(synth$regs, function(g) match(g, genes)), list(EXTRA = match(genes[c(121:132)], genes))); list(members = s, class = stats::setNames(rep("TF", length(s)), names(s))) }
  B <- list()
  B$complexes <- cx("main", if (fallback) none_op else operon_of, "primary")
  B$regulons_SC <- mk_truth(reg_sc, operon_of, genes, "primary")
  B$complexes_operon_incl <- cx("main", none_op, "sensitivity")
  B$complexes_size2 <- cx("size2", operon_of, "sensitivity")
  B$complexes_nomerge <- cx("nomerge", operon_of, "sensitivity")
  B$complexes_noribo <- cx("noribo", operon_of, "sensitivity")
  B$regulons_all <- mk_truth(reg_all, operon_of, genes, "descriptive")
  if (is.null(synth)) {                       # definitions of analysis/granularity.R, as keff_test.R load_truths
    ra <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
    pick <- function(keep) list(members = ra$members[keep], class = ra$class[keep])
    B$TF_regulons <- mk_truth(pick(ra$class == "TF" & lengths(ra$members) <= 500 & lengths(ra$members) >= 5), operon_of, genes, "descriptive")
    B$sigmulons <- mk_truth(pick(ra$class == "sigma" & lengths(ra$members) >= 5 & !names(ra$members) %in% c("sigma:sigma70", "sigma:sigma38")), operon_of, genes, "descriptive")
    op <- operon_of[names(operon_of) %in% genes]
    op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
    B$operons <- mk_truth(as_reg(op_sets, "operon"), none_op, genes, "descriptive")
  }
  B <- Filter(Negate(is.null), B)
  if (!is.null(only)) B <- B[intersect(only, names(B))]
  rp <- sets_of("rp_module")
  list(bench = B, operon_of = operon_of, rp_module = if (length(rp)) as_reg(rp, "rp") else NULL,
       evaluability = data.frame(dataset = ds, n_positive_pairs_main_operon_excluded = n_pos,
                                 n_complexes_with_positives = n_cx, fallback_operon_included = fallback),
       complex_truth = "complexes", sc_truth = "regulons_SC")
}
build_truths_human <- function(ds, genes, h5, B_ids, synth = NULL) {
  raw <- if (is.null(synth)) corum_sets(genes)$members else
    lapply(synth$complexes, function(g) sort(stats::na.omit(match(g, genes))))
  raw <- raw[lengths(raw) >= 3L & lengths(raw) <= 200L]
  mg <- merge_sets(lapply(raw, function(ix) genes[ix]))
  mg <- lapply(mg, function(g) sort(match(g, genes)))
  rpfrac <- vapply(mg, function(ix) mean(grepl(RP_RE, genes[ix])), 0)
  B <- list(CORUM_merged = mk_truth(as_reg(mg, "complex"), none_op, genes, "primary"),
            CORUM_raw = mk_truth(as_reg(raw, "complex"), none_op, genes, "sensitivity"),
            CORUM_noRP = mk_truth(as_reg(mg[rpfrac < 0.3], "complex"), none_op, genes, "sensitivity"))
  if (is.null(synth)) {
    B$CollecTRI <- mk_truth(collectri_sets(genes), none_op, genes, "descriptive")
    B$DoRothEA_AB <- mk_truth(dorothea_sets(genes), none_op, genes, "descriptive")
    B$DoRothEA_ABC <- mk_truth(dorothea_sets(genes, c("A", "B", "C")), none_op, genes, "descriptive")
    if (!is.null(h5)) B$KD_response <- mk_truth(response_sets(genes, path = h5, use = B_ids), none_op, genes, "descriptive")
  }
  list(bench = Filter(Negate(is.null), B), operon_of = none_op, rp_module = NULL, evaluability = NULL,
       complex_truth = "CORUM_merged", sc_truth = NULL)
}
build_truths <- function(ds, genes, ld, only = NULL) {
  if (ds %in% c("m3d", "ecomac")) build_truths_ecoli(ds, genes, ld$synth, only)
  else build_truths_human(ds, genes, ld$h5, ld$B, ld$synth)
}
truth_counts <- function(TR) do.call(rbind, lapply(names(TR$bench), function(nm) { b <- TR$bench[[nm]]
  data.frame(truth = nm, role = b$role, n_sets = length(b$reg$members), n_genes = length(b$pairs$genes_in), n_pairs = length(b$pairs$label),
             n_pos = sum(b$pairs$label), base_rate = mean(b$pairs$label), n_same_operon_removed = b$pairs$n_same_operon_removed) }))

## ---- full-data scores for one dataset (resumable at two points: fit, pair scores) ---------------------
score_dataset <- function(ds, X, genes, TR, sel, dout, say) {
  G <- length(genes); md5 <- digest::digest(genes, algo = "md5"); bench <- TR$bench
  fp <- list(ds = ds, genes_md5 = md5, truths = names(bench), npairs = vapply(bench, function(b) length(b$pairs$label), 0L),
             prereg = HASH, tstar = sel, alpha = opt$alpha, keff = opt$keff, cap = opt$cap, hubstrat = opt$hubstrat)
  cf <- file.path(dout, "cache", "scores.rds"); dir.create(file.path(dout, "cache"), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(cf)) { s <- readRDS(cf); if (!isTRUE(all.equal(s$fp, fp))) stop("cache/scores.rds fingerprint mismatch: use a new --out")
    say("pair scores: reusing ", cf); return(s) }
  fitf <- file.path(dout, "cache", "fit.rds")
  if (file.exists(fitf)) { ft <- readRDS(fitf); stopifnot(identical(ft$genes_md5, md5), identical(ft$prereg, HASH)); say("fit: reusing ", fitf)
  } else {
    t0 <- Sys.time()
    f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$
      calibrate(method = "normal", combine = "stouffer")
    f$build_operator(alpha = opt$alpha, softmax_keff = opt$keff, softmax_cap = opt$cap)   # exactly as repro mode (build_fit / build_operator)
    ft <- list(C = as.matrix(f$clr_scores), A = as.matrix(f$operator), genes_md5 = md5, prereg = HASH)
    dimnames(ft$C) <- NULL; dimnames(ft$A) <- NULL; f$release(); rm(f); invisible(gc())
    saveRDS(ft, fitf, compress = FALSE); say(sprintf("MI, CLR, operator in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  A <- ft$A; C <- ft$C; rm(ft)
  # stationary distribution of the lazy walk (teleport 1e-3), as baselines.R (sparse P)
  Ps <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  q0 <- rep(1 / G, G)
  for (it in 1:20000) { q <- (1 - 1e-3) * as.numeric(q0 %*% Ps) + 1e-3 / G; if (sum(abs(q - q0)) < 1e-13) break; q0 <- q }
  rm(Ps); Hm <- if (opt$hubstrat == 1L) outer(q, q, "+") else NULL
  cname <- TR$complex_truth; co <- NULL
  if (!is.null(cname) && cname %in% names(bench)) co <- comember_matrix(bench[[cname]]$reg, G)
  rp <- TR$rp_module; co_rp <- if (!is.null(rp)) comember_matrix(rp, G) else NULL
  HUB <- c("attention_t9", "a3r_t9", "clr", "abs_r"); CORE <- HUB; bench_role <- vapply(bench, function(b) b$role, ""); PCM <- c("attention_t9", "a3r_t9", "clr")
  SC <- list(); PC <- list(); RPC <- list()
  consume <- function(S, nm) {
    # pair scores: non-finite scores rank last (ap_fast treats NA as last); other uses below see 0
    Sna <- S; Sna[!is.finite(Sna)] <- NA; S[!is.finite(S)] <- 0; diag(S) <- 0
    # descriptive truths (huge pair sets on human data) keep only the core methods, to bound memory
    use <- names(bench)[bench_role != "descriptive" | nm %in% CORE]
    SC[[nm]] <<- lapply(bench[use], function(b) Sna[cbind(b$pairs$i, b$pairs$j)]); rm(Sna)
    if (!is.null(Hm) && nm %in% HUB) { hs <- hub_stratify(S, Hm)
      SC[[paste0(nm, "_hubstrat")]] <<- lapply(bench[use], function(b) hs[cbind(b$pairs$i, b$pairs$j)]); rm(hs) }
    if (nm %in% PCM && !is.null(co)) { b <- bench[[cname]]       # per-complex AP, complexes with >= 10 positive pairs (macro_aupr definitions)
      x <- regulon_ap_table(S, b$reg, co, b$operon, genes, min_pos = 10L); PC[[nm]] <<- cbind(method = nm, x) }
    if (nm %in% PCM && !is.null(rp)) RPC[[nm]] <<- data.frame(method = nm,
      coherence_auroc = regulon_coherence(S, rp, co_rp, TR$operon_of, genes))
    invisible(NULL)
  }
  ts <- sort(unique(c(T_PRIMARY, unlist(sel)))); Pd <- lazy_walk_P(A, opt$alpha); rm(A); Pt <- Pd; lf <- NULL
  for (t in seq_len(max(ts))) {
    if (t > 1L) Pt <- Pt %*% Pd
    if (!(t %in% ts)) next
    say("depth t = ", t)
    if (t == T_PRIMARY) { consume(attention(Pt), "attention_t9"); consume(a3r(Pt), "a3r_t9"); consume(a2r(Pt), "a2r_t9"); consume(a1(Pt), "a1_t9")
      lf <- cbind(dataset = ds, t = t, operator_checks(Pt)) }
    if (!is.null(sel$a3r) && t == sel$a3r) consume(a3r(Pt), "a3r_tstar")
    if (!is.null(sel$attention) && t == sel$attention) consume(attention(Pt), "attention_tstar")
    invisible(gc(verbose = FALSE))
  }
  rm(Pt, Pd); consume(C, "clr"); rm(C)
  r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; consume(r, "abs_r"); rm(r, Hm)
  s <- list(SC = SC, PC = PC, RPC = RPC, lf = lf, fp = fp); saveRDS(s, cf); s
}

## ---- jackknife, metrics, endpoints ---------------------------------------------------------------------
# Paired delete-half gene jackknife (baselines.R): per truth, 100 replicates; a random half of the truth's genes
# (genes_in) is kept, AUPR over pairs with both genes kept; the same halves serve every method.
# Halves are drawn for all truths of this dataset, in table order, right after set.seed(20261022).
jack_dataset <- function(bench, SC, G, dout, say) {
  set.seed(JACK_SEED)
  halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r)
    sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))))
  JR <- list(); dir.create(file.path(dout, "cache"), showWarnings = FALSE)
  for (nm in names(bench)) {
    jf <- file.path(dout, "cache", paste0("jack_", nm, ".rds"))
    methods <- names(SC)[vapply(SC, function(x) nm %in% names(x), NA)]
    if (file.exists(jf)) { J <- readRDS(jf); if (identical(colnames(J), methods)) { JR[[nm]] <- J; next } }
    p <- bench[[nm]]$pairs; J <- matrix(NA_real_, opt$jack, length(methods), dimnames = list(NULL, methods))
    for (r in seq_len(opt$jack)) {
      inh <- logical(G); inh[halves[[nm]][[r]]] <- TRUE; ix <- which(inh[p$i] & inh[p$j]); lab <- p$label[ix]
      for (m in methods) J[r, m] <- ap_fast(SC[[m]][[nm]][ix], lab)
    }
    saveRDS(J, jf); JR[[nm]] <- J; say("jackknife done: ", nm)
  }
  JR
}
comparators <- function(m) {                  # same-depth, same-adjustment comparators
  hs <- endsWith(m, "_hubstrat"); b <- sub("_hubstrat$", "", m)
  c(att = if (b == "a3r_tstar") paste0("attention_tstar", if (hs) "_hubstrat" else "") else paste0("attention_t9", if (hs) "_hubstrat" else ""),
    clr = paste0("clr", if (hs) "_hubstrat" else ""))
}
metrics_dataset <- function(ds, bench, SC, JR) {
  rows <- list()
  for (nm in names(bench)) { p <- bench[[nm]]$pairs; J <- JR[[nm]]
    ms <- colnames(J); au <- vapply(ms, function(m) ap_fast(SC[[m]][[nm]], p$label), 0)
    for (m in ms) { cp <- comparators(m); row <- list(dataset = ds, truth = nm, role = bench[[nm]]$role, method = m,
        aupr = au[[m]], base_rate = mean(p$label), n_pos = sum(p$label))
      for (k in names(cp)) { c0 <- cp[[k]]
        if (!c0 %in% ms) { d0 <- se <- NA_real_ } else { d0 <- au[[m]] - au[[c0]]; se <- stats::sd(J[, m] - J[, c0]) }
        row[[paste0("d_vs_", k)]] <- d0; row[[paste0("se_vs_", k)]] <- se
        row[[paste0("lo_vs_", k)]] <- d0 - 1.96 * se; row[[paste0("hi_vs_", k)]] <- d0 + 1.96 * se
        row[[paste0("p_sup_vs_", k)]] <- p_sup(d0, se); row[[paste0("p_ni_vs_", k)]] <- p_ni(d0, se) }
      rows[[length(rows) + 1L]] <- as.data.frame(row) } }
  do.call(rbind, rows)
}
# endpoint table (preregistration "Endpoints"): new readout a3r at t = 9 minus the comparator
ENDPOINTS <- data.frame(
  code = c("P1", "P2", "P3", "P4", "S1", "S2", "S3", "S4", "S5"),
  dataset = c("m3d", "ecomac", "m3d", "ecomac", "m3d", "ecomac", "m3d", "ecomac", "rpe1"),
  truth_kind = c("complex", "complex", "sc", "sc", "sc", "sc", "complex", "complex", "complex"),
  comparator = c("att", "att", "att", "att", "clr", "clr", "clr", "clr", "att"),
  hypothesis = c("sup", "sup", "ni", "ni", "sup", "sup", "ni", "ni", "sup"), stringsAsFactors = FALSE)
endpoints_dataset <- function(ds, met, TR) {
  e <- ENDPOINTS[ENDPOINTS$dataset == ds, ]; if (!nrow(e)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(e)), function(k) {
    tn <- if (e$truth_kind[k] == "complex") TR$complex_truth else TR$sc_truth
    x <- met[met$truth == tn & met$method == "a3r_t9", ]
    cmp <- e$comparator[k]; d <- x[[paste0("d_vs_", cmp)]]; se <- x[[paste0("se_vs_", cmp)]]
    data.frame(code = e$code[k], dataset = ds, truth = tn, comparator = if (cmp == "att") "attention_t9" else "clr", hypothesis = e$hypothesis[k],
               margin = if (e$hypothesis[k] == "ni") NI_MARGIN else 0, aupr_new = x$aupr, d = d, se = se, lo = d - 1.96 * se, hi = d + 1.96 * se,
               p_raw = if (e$hypothesis[k] == "sup") p_sup(d, se) else p_ni(d, se), p_loss = p_loss(d, se), stringsAsFactors = FALSE) }))
}
per_complex_summary <- function(PC) {
  if (!all(c("attention_t9", "a3r_t9") %in% names(PC))) return(NULL)
  a <- PC$attention_t9; b <- PC$a3r_t9; ok <- !is.na(a$ap) & !is.na(b$ap); d <- b$ap[ok] - a$ap[ok]
  r <- data.frame(comparison = "a3r_t9 minus attention_t9", n_complexes_ge10pos = sum(ok), mean_d_ap = mean(d), median_d_ap = stats::median(d),
                  n_improved = sum(d > 0), n_worse = sum(d < 0))
  if ("clr" %in% names(PC)) { cc <- PC$clr; ok2 <- !is.na(a$ap) & !is.na(cc$ap); d2 <- a$ap[ok2] - cc$ap[ok2]
    r <- rbind(r, data.frame(comparison = "attention_t9 minus clr", n_complexes_ge10pos = sum(ok2), mean_d_ap = mean(d2), median_d_ap = stats::median(d2),
                             n_improved = sum(d2 > 0), n_worse = sum(d2 < 0))) }
  r
}

## ---- mode: truth ---------------------------------------------------------------------------------------
if (opt$mode == "truth") {
  ds <- opt$dataset; if (!ds %in% DATASETS) stop("truth: --dataset must be one of ", paste(DATASETS, collapse = ", "))
  gate_all(ds)
  dout <- file.path(opt$out, "truth", ds); dir.create(dout, recursive = TRUE, showWarnings = FALSE)
  say <- logger(file.path(dout, "truth.log")); T0 <- Sys.time()
  say("truth phase: ", ds, if (SYN) " (SYNTHETIC smoke test)" else "", "; prereg commit ", HASH, "; HEAD ", GIT_HEAD)
  if (file.exists(file.path(dout, "DONE")) ) stop("truth phase for ", ds, " is already DONE; a rerun needs a logged reason (preregistration step 5)")
  ld <- load_expr(ds); X <- ld$X; genes <- rownames(X); G <- length(genes)
  fpf <- file.path(opt$out, ds, "splits", "fingerprint.rds")
  if (!SYN && ds %in% c("m3d", "ecomac", "rpe1")) {          # same input as the committed selections
    if (!file.exists(fpf) || !identical(readRDS(fpf)$genes_md5, digest::digest(genes, algo = "md5")))
      stop("gene set differs from the one the selections were computed on (", fpf, ")")
  }
  say(sprintf("%s: %d genes x %d samples", ds, G, ncol(X)))
  tsel <- list()
  if (SYN) tsel <- list(a3r = 4L, attention = 9L) else {
    sf <- file.path(SELDIR, paste0(ds, "_selections.csv"))
    if (file.exists(sf)) { s <- utils::read.csv(sf); tsel <- list(a3r = s$t_star[s$readout == "a3r"], attention = s$t_star[s$readout == "attention"]) }
    else say("no committed selections for ", ds, ": t* readouts skipped")
  }
  TR <- build_truths(ds, genes, ld)
  tc <- truth_counts(TR); wout(tc, file.path(dout, "truth_counts.csv"))
  if (!is.null(TR$evaluability)) { wout(TR$evaluability, file.path(dout, "evaluability.csv"))
    say("evaluability: ", paste(names(TR$evaluability), unlist(TR$evaluability[1, ]), sep = "=", collapse = " ")) }
  say("truths: ", paste(sprintf("%s (%d pos / %d pairs)", tc$truth, tc$n_pos, tc$n_pairs), collapse = "; "))
  need <- c(TR$complex_truth, TR$sc_truth); if (!all(need %in% names(TR$bench))) stop("primary truth missing: ", paste(setdiff(need, names(TR$bench)), collapse = ", "))
  S <- score_dataset(ds, X, genes, TR, tsel, dout, say); rm(X)
  wout(S$lf, file.path(dout, "label_free_t9.csv"))
  if (ds == "k562ess" && "CORUM_raw" %in% names(TR$bench)) {   # consistency with the implementation check (descriptive)
    a <- ap_fast(S$SC$a3r_t9$CORUM_raw, TR$bench$CORUM_raw$pairs$label)
    say(sprintf("K562 raw CORUM a3r t = 9 AUPR %.9f vs repro %.7f: %s", a, REPRO_AUPR, if (identical(sig6(a), sig6(REPRO_AUPR))) "matches at 6 digits" else "DIFFERS"))
  }
  JR <- jack_dataset(TR$bench, S$SC, G, dout, say)
  jr <- do.call(rbind, lapply(names(JR), function(nm) { J <- JR[[nm]]
    data.frame(truth = nm, method = rep(colnames(J), each = nrow(J)), rep = rep(seq_len(nrow(J)), ncol(J)), aupr = as.vector(J)) }))
  wout(jr, file.path(dout, "jack_replicates.csv"))
  met <- metrics_dataset(ds, TR$bench, S$SC, JR); wout(met, file.path(dout, "metrics.csv"))
  ep <- endpoints_dataset(ds, met, TR); if (!is.null(ep)) wout(ep, file.path(dout, paste0("endpoints_", ds, ".csv")))
  if (length(S$PC)) { wout(do.call(rbind, S$PC), file.path(dout, "per_complex_ap.csv"))
    pcs <- per_complex_summary(S$PC); if (!is.null(pcs)) wout(pcs, file.path(dout, "per_complex_summary.csv")) }
  if (length(S$RPC)) wout(do.call(rbind, S$RPC), file.path(dout, "rp_module_coherence.csv"))
  m9 <- met[met$method %in% c("a3r_t9", "attention_t9", "clr") & met$role != "descriptive", c("truth", "method", "aupr", "d_vs_att", "se_vs_att", "d_vs_clr", "se_vs_clr")]
  say("AUPR (primary and sensitivity truths):\n", paste(utils::capture.output(print(m9, row.names = FALSE, digits = 4)), collapse = "\n"))
  if (!is.null(ep)) say("endpoints:\n", paste(utils::capture.output(print(ep[, c("code", "truth", "comparator", "hypothesis", "d", "se", "p_raw")], row.names = FALSE, digits = 4)), collapse = "\n"))
  writeLines(sprintf("%s %s %s", HASH, GIT_HEAD, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file.path(dout, "DONE"))
  say(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
  quit(save = "no")
}

## ---- mode: boot (M3D experiment-level cluster bootstrap, chips_bootstrap.R cluster mode) ---------------
if (opt$mode == "boot") {
  if (opt$dataset != "m3d") stop("boot: M3D only (EcoMAC has no such bootstrap)")
  gate_all("m3d"); source(file.path(script_dir, "parallel.R"))
  dout <- file.path(opt$out, "boot"); dir.create(file.path(dout, "draws"), recursive = TRUE, showWarnings = FALSE)
  say <- logger(file.path(dout, "boot.log")); T0 <- Sys.time()
  say("cluster bootstrap: ", opt$draws, " draws, seeds ", BOOT_SEED0, " + b; prereg commit ", HASH, if (SYN) " (SYNTHETIC)" else "")
  if (SYN) { sd0 <- synth_data("m3d"); ch <- synth_chips(sd0); Xc <- ch$Xc; Xa <- ch$Xa; exp_of <- ch$exp_of; synth <- sd0
  } else {                                    # construction of chips_bootstrap.R
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d); synth <- NULL
    rdm <- function(f) { raw <- utils::read.delim(file.path(m3d, f), check.names = FALSE, stringsAsFactors = FALSE)
      X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"; rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); X }
    Xc <- rdm(paste0(base, "_chips907probes4297.tab")); Xa <- rdm(paste0("avg_", base, "_exps466probes4297.tab"))
    stopifnot(identical(rownames(Xc), rownames(Xa)))
    keep <- apply(Xc, 1, function(x) diff(range(x))) > 0 & apply(Xa, 1, function(x) diff(range(x))) > 0
    Xc <- Xc[keep, ]; Xa <- Xa[keep, ]
    desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
    exp_of <- stats::setNames(desc$experiment_name, desc$chip_name)[colnames(Xc)]; stopifnot(!anyNA(exp_of))
  }
  genes <- rownames(Xc); G <- length(genes); chips_by_exp <- split(colnames(Xc), exp_of)
  say(sprintf("%d genes; %d chips of %d experiments", G, ncol(Xc), length(chips_by_exp)))
  if (!SYN && !identical(genes, rownames(Xa)[apply(Xa, 1, function(x) diff(range(x))) > 0])) say("note: bootstrap gene set differs from the 466-experiment truth gene set")
  ld <- list(synth = synth, h5 = NULL, B = NULL)
  TR <- build_truths("m3d", genes, ld, only = c("complexes", "regulons_SC")); bench <- TR$bench
  stopifnot(all(c("complexes", "regulons_SC") %in% names(bench)))
  say("truths: ", paste(sprintf("%s %d pos / %d pairs", names(bench), vapply(bench, function(b) sum(b$pairs$label), 0), vapply(bench, function(b) length(b$pairs$label), 0L)), collapse = "; "))
  draw_cols <- function() {                   # chips_bootstrap.R draw_cols("cluster"), verbatim
    ex <- names(chips_by_exp); pick <- sample(ex, length(ex), replace = TRUE); used <- list()
    vapply(pick, function(e) { ch <- chips_by_exp[[e]]; avail <- setdiff(ch, used[[e]]); if (!length(avail)) avail <- ch
      c1 <- if (length(avail) == 1L) avail else sample(avail, 1L); used[[e]] <<- c(used[[e]], c1); c1 }, "")
  }
  boot_eval <- function(Xm, label) {
    rng <- apply(Xm, 1, function(x) diff(range(x)))
    if (any(rng <= 0)) Xm[rng <= 0, ] <- Xm[rng <= 0, ] + matrix(stats::rnorm(sum(rng <= 0) * ncol(Xm), sd = 1e-9), sum(rng <= 0))
    f <- ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$calibrate(method = "normal", combine = "stouffer")
    f$build_operator(alpha = opt$alpha, softmax_keff = opt$keff, softmax_cap = opt$cap)
    A <- as.matrix(f$operator); C <- as.matrix(f$clr_scores); f$release(); rm(f)
    Pt <- power_t(lazy_walk_P(A, opt$alpha), T_PRIMARY); rm(A)
    S <- list(attention = attention(Pt), a3r = a3r(Pt)); rm(Pt); S$clr <- C; rm(C)
    out <- list()
    for (m in names(S)) { s <- S[[m]]; s[!is.finite(s)] <- NA
      for (nm in names(bench)) { p <- bench[[nm]]$pairs
        out[[length(out) + 1L]] <- data.frame(draw = label, truth = nm, method = m, aupr = ap_fast(s[cbind(p$i, p$j)], p$label)) } }
    do.call(rbind, out)
  }
  tasks <- list()
  for (b in seq_len(opt$draws)) { sd_b <- BOOT_SEED0 + b; set.seed(sd_b)
    tasks[[b]] <- list(label = sprintf("cluster_%03d", b), seed = sd_b, cols = draw_cols()) }
  done <- sub("\\.csv$", "", list.files(file.path(dout, "draws"), "\\.csv$"))
  tasks <- Filter(function(t) !t$label %in% done, tasks)
  say(sprintf("%d draws done, %d to run", length(done), length(tasks)))
  if (length(tasks)) {
    W <- min(opt$workers, length(tasks))
    pc <- par_start(mem_gb = opt$mem_gb, workers = W, reserve = max(0L, .par_cores() - W * opt$threads), say = say,
                    export = setdiff(ls(globalenv()), c("pc", "tasks")))
    st <- foreach::foreach(task = tasks, .inorder = FALSE, .errorhandling = "pass") %dopar% {
      if (par_should_stop()) return("stopped")
      set.seed(task$seed); t0 <- Sys.time()
      d <- boot_eval(Xc[, task$cols, drop = FALSE], task$label)
      tmp <- file.path(dout, "draws", paste0(task$label, ".csv.part")); utils::write.csv(d, tmp, row.names = FALSE)
      file.rename(tmp, file.path(dout, "draws", paste0(task$label, ".csv"))); par_release()
      sprintf("%s %.1f min", task$label, as.numeric(difftime(Sys.time(), t0, units = "mins")))
    }
    par_stop(pc)
    bad <- vapply(st, function(x) inherits(x, "error"), NA)
    if (any(bad)) say("draws failed: ", paste(vapply(st[bad], conditionMessage, ""), collapse = " | "))
    if (any(bad) || any(unlist(st[!bad]) == "stopped")) stop("not all draws finished; rerun to resume")
  }
  parts <- list.files(file.path(dout, "draws"), "\\.csv$", full.names = TRUE)
  d <- do.call(rbind, lapply(parts, utils::read.csv)); stopifnot(length(unique(d$draw)) == opt$draws)
  wout(d, file.path(dout, "boot_draws.csv"))
  w <- stats::reshape(d, idvar = c("draw", "truth"), timevar = "method", direction = "wide")
  q <- function(x, p) stats::quantile(x, p, names = FALSE)
  sm <- do.call(rbind, lapply(c(P1 = "complexes", P3 = "regulons_SC"), function(tn) { x <- w[w$truth == tn, ]
    data.frame(code = c("P1", "P3")[match(tn, c("complexes", "regulons_SC"))], truth = tn, n_draws = nrow(x),
               d_vs_attention_mean = mean(x$aupr.a3r - x$aupr.attention), d_vs_attention_lo = q(x$aupr.a3r - x$aupr.attention, .025),
               d_vs_attention_hi = q(x$aupr.a3r - x$aupr.attention, .975), d_vs_clr_mean = mean(x$aupr.a3r - x$aupr.clr),
               d_vs_clr_lo = q(x$aupr.a3r - x$aupr.clr, .025), d_vs_clr_hi = q(x$aupr.a3r - x$aupr.clr, .975)) }))
  wout(sm, file.path(dout, "boot_summary.csv"))
  say("bootstrap intervals (a3r minus attention, 2.5-97.5%):\n", paste(utils::capture.output(print(sm, row.names = FALSE, digits = 4)), collapse = "\n"))
  say(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
  quit(save = "no")
}

## ---- mode: verdict (computed from the preregistered rules; no numbers are interpreted here) -----------
# decide(): pure function of the endpoint table, bootstrap intervals and the repro flag; unit-testable.
decide <- function(ep, boot, repro_ok) {
  g <- function(code, col) ep[[col]][ep$code == code]
  P <- c("P1", "P2", "P3", "P4"); have_p <- all(P %in% ep$code)
  if (!have_p) stop("verdict needs P1-P4 (run M3D and EcoMAC first)")
  pP <- vapply(P, function(k) g(k, "p_raw"), 0); dP <- vapply(P, function(k) g(k, "d"), 0)
  holm_p <- stats::p.adjust(pP, "holm"); loss_holm <- stats::p.adjust(vapply(P, function(k) g(k, "p_loss"), 0), "holm")
  rej <- pP < ALPHA1; rej[is.na(rej)] <- FALSE
  r1 <- isTRUE(rej[1] && dP[1] >= EFFECT_MIN && rej[2] && dP[2] >= EFFECT_MIN)
  r2 <- isTRUE(rej[3] && rej[4])
  r3 <- if (is.null(boot)) NA else isTRUE(boot$d_vs_attention_lo[boot$code == "P1"] > 0 && boot$d_vs_attention_lo[boot$code == "P3"] > -NI_MARGIN)
  r4 <- isTRUE(repro_ok)
  adopt <- if (is.na(r3)) NA else (r1 && r2 && r3 && r4)
  prim <- data.frame(code = P, d = dP, se = vapply(P, function(k) g(k, "se"), 0), p_raw = pP, reject_alpha = rej,
                     p_holm_primary = holm_p, p_loss = vapply(P, function(k) g(k, "p_loss"), 0), p_loss_holm = loss_holm,
                     effect_ge_min = ifelse(P %in% c("P1", "P2"), dP >= EFFECT_MIN, NA), row.names = NULL)
  loss <- any(dP[3:4] <= -NI_MARGIN & loss_holm[3:4] < ALPHA1, na.rm = TRUE)
  outc <- data.frame(outcome = c("Complexes improve, regulon NI not shown on M3D or EcoMAC", "Regulon loss on M3D or EcoMAC (d <= -0.005, Holm-adjusted p < 0.025)",
                                 "Complex superiority not shown on M3D or EcoMAC", "Gene jackknife passes, cluster bootstrap fails (rule 3)", "Implementation check fails"),
                     applies = c(r1 && !r2, loss, !r1, if (is.na(r3)) NA else (r1 && r2 && !r3), !r4),
                     conclusion = c("Not adopted. Reported as a complex-specific readout with an unresolved regulon cost.",
                                    "Not adopted. Reported as a trade-off: the readout harms regulons.",
                                    "Not adopted. K562 rescue reported as exploratory and not replicated in E. coli. (Attention vs CLR on E. coli complexes is descriptive; see metrics.csv.)",
                                    "Not adopted. Supported over genes, not robust to resampling experiments.",
                                    "No E. coli scoring until fixed; fix and cause logged in the preregistration first."), stringsAsFactors = FALSE)
  S <- ep[ep$code %in% paste0("S", 1:5), ]
  sec <- NULL
  if (nrow(S)) { sec <- S[, c("code", "d", "se", "p_raw")]; sec$interpretable <- isTRUE(adopt)
    sec$p_holm <- NA_real_; sec$reject <- NA
    if (isTRUE(adopt) && nrow(S) == 5L) { sec$p_holm <- stats::p.adjust(sec$p_raw, "holm"); sec$reject <- sec$p_holm < ALPHA1 } }
  s_pending <- !(nrow(S) == 5L)
  notes <- character(0)
  if (isTRUE(adopt) && !s_pending) {
    rj <- stats::setNames(sec$reject, sec$code)
    if (!isTRUE(rj[["S1"]]) || !isTRUE(rj[["S2"]])) notes <- c(notes, "S1 or S2 failed: the paper does not claim the new method beats CLR on regulons; the headline over CLR stays attributed to attention.")
    if (!isTRUE(rj[["S5"]])) notes <- c(notes, "S5 failed: the human complex result did not replicate in RPE1; the complex claim rests on E. coli plus in-sample K562.")
    sm <- which(dP[3:4] < 0 & rej[3:4] & !is.na(loss_holm[3:4]) & pP[3:4] < ALPHA1 & vapply(P[3:4], function(k) g(k, "p_loss"), 0) < ALPHA1)
    if (length(sm)) notes <- c(notes, "A small but significant regulon loss within the margin (NI shown, d < 0 significant) is reported as such.")
  }
  list(prim = prim, rules = data.frame(rule1_complexes = r1, rule2_regulon_ni = r2, rule3_bootstrap = r3, rule4_repro = r4, adopt = adopt,
                                       secondary_pending = s_pending), outcomes = outc, secondary = sec, notes = notes)
}
if (opt$mode == "verdict") {
  fix_check(); if (!SYN) { gate_script_clean(); gate_repro(HASH) }
  rd1 <- function(f) if (file.exists(f)) utils::read.csv(f) else NULL
  dd <- function(ds) file.path(opt$out, "truth", ds)
  for (n in c("m3d", "ecomac")) if (!file.exists(file.path(dd(n), "DONE"))) stop("verdict needs the finished truth phase of ", n)
  ep <- do.call(rbind, lapply(c("m3d", "ecomac", "rpe1"), function(n) rd1(file.path(dd(n), paste0("endpoints_", n, ".csv")))))
  hashes <- unique(ep$prereg_commit); if (length(hashes) != 1L || hashes != HASH) stop("an endpoint file was produced under another preregistration commit")
  boot <- rd1(file.path(opt$out, "boot", "boot_summary.csv"))
  repro_ok <- if (SYN) TRUE else isTRUE(tryCatch(gate_repro(HASH), error = function(e) FALSE))
  dec <- decide(ep, boot, repro_ok)
  vdir <- file.path(opt$out, "verdict"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
  wout(ep, file.path(vdir, "endpoints_all.csv")); wout(dec$prim, file.path(vdir, "primary.csv"))
  if (!is.null(dec$secondary)) wout(dec$secondary, file.path(vdir, "secondary.csv"))
  wout(dec$outcomes, file.path(vdir, "outcome_table.csv"))
  if (!is.null(boot)) wout(boot, file.path(vdir, "bootstrap_intervals.csv"))
  m3 <- rd1(file.path(dd("m3d"), "metrics.csv")); em <- rd1(file.path(dd("ecomac"), "metrics.csv"))
  vc <- do.call(rbind, lapply(list(m3d = m3, ecomac = em), function(m) m[m$truth == "complexes" & m$method == "attention_t9", c("dataset", "truth", "method", "d_vs_clr", "se_vs_clr", "lo_vs_clr", "hi_vs_clr")]))
  wout(vc, file.path(vdir, "descriptive_attention_vs_clr_complexes.csv"))
  k <- rd1(file.path(dd("k562ess"), "metrics.csv")); kn <- NULL
  if (!is.null(k)) { x <- k[k$truth == "CORUM_merged" & k$method == "a3r_t9", ]
    if (nrow(x) && x$d_vs_att <= 0) kn <- "K562 merged-CORUM difference vs attention is <= 0 (descriptive, in-sample): report prominently." }
  verdict <- if (is.na(dec$rules$adopt)) "PENDING (cluster bootstrap missing)" else if (dec$rules$adopt) "ADOPT CLR(log P^t) as main method; attention is the ablation" else "NOT ADOPTED: failure to confirm; main method stays attention"
  v <- cbind(data.frame(verdict = verdict), dec$rules, data.frame(notes = paste(c(dec$notes, kn), collapse = " | "),
        outcomes_applying = paste(dec$outcomes$outcome[which(dec$outcomes$applies %in% TRUE)], collapse = " | ")))
  wout(v, file.path(vdir, if (SYN) "verdict_SYNTHETIC.csv" else "verdict.csv"))
  cat(verdict, "\n"); print(dec$prim, digits = 4); print(dec$rules); quit(save = "no")
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
