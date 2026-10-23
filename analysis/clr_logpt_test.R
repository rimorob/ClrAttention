#!/usr/bin/env Rscript
# CLR(log P^t) test, as preregistered in docs/preregistration_clr_logpt.md (committed 90b0e3d).
# This file implements: the gates, the a3r operator (and a2r, a1) and the K562 reproduction check.
# NOT implemented yet (later commits, after the check passes): label-free t* phase, EcoCyc truth build,
# E. coli / RPE1 / K562 truth phases. Their entry points exist only as gate-then-stop stubs.
#
# Modes (--mode):
#   repro   K562 implementation check. Reads ONLY the K562 raw-CORUM truth rows and must reproduce the
#           exploratory AUPR 0.1054351 (a_clr_bilateral_logP_t9_robust) to 6 significant digits.
#           Writes <out>/k562_repro/repro_check.csv (passed = TRUE/FALSE). --quick N uses the first N
#           genes (seeded sample) and reads NO truth (operator diagnostics only; never counts as the check).
#   selections, truth   gate-checked stubs (see gates below).
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
            expr = "results/baselines_k562ess_truthgenes_genie3/expression.tsv",
            members = "results/baselines_k562ess_truthgenes_genie3/regulon_members.csv")
ints <- c("threads", "quick", "t", "cap", "seed")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% ints) as.integer(v) else if (key %in% c("alpha", "keff")) as.numeric(v) else v
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

## ---- mode: selections / truth stubs -------------------------------------------------------------
if (opt$mode %in% c("selections", "truth")) {
  if (opt$mode == "truth") {
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
