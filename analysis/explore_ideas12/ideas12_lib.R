# EXPLORATORY, NOT PREREGISTERED, NOTHING CONFIRMATORY. Shared code for the exploratory runs of ideation-memo ideas 1 and 2
# (docs/plans/explore_ideas12_plan.md). K562 truths are spent; nothing here goes into a paper or pre-read.
# Truth sets, jackknife stream, scoring and hub stratification are copied from analysis/robust_z/robust_z.R (read-only reuse of
# analysis/robust_z/robust_z_lib.R for the operator helpers).
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/robust_z/robust_z_lib.R")
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
EI_OUT <- "results/explore_ideas12"
EI_SEED <- 20261016L; EI_DEPTH_SEED <- 20261014L; EI_T <- 9L; EI_ALPHA <- 0.5; EI_JACK <- 100L
ei_say <- function(f) function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = f, append = TRUE) }
none <- stats::setNames(character(0), character(0))

## ---- K562 essential half A: data and truths (robust_z.R, k562 branch) ---------------------------------------------------
ei_setup_k562 <- function(say = function(...) cat(..., "\n"), workers = 6L) {
  raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; storage.mode(X) <- "double"
  genes <- rownames(X); G <- length(genes); lab <- colnames(X)
  say(sprintf("k562: %d genes x %d samples", G, ncol(X)))
  h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"
  set.seed(20261015L)
  rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  pert <- cat_of("/obs", "gene_transcript"); pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A_half <- sort(K[seq_len(floor(length(K) / 2))]); B_half <- setdiff(K, A_half); rm(rd, cat_of, pert, pval, K)
  tr <- list(DoRothEA_AB = dorothea_sets(genes), CORUM = corum_sets(genes), CollecTRI = collectri_sets(genes), KD_response = response_sets(genes, path = h5, use = B_half))
  cr <- tr$CORUM; cx <- utils::read.csv("analysis/corum_exploratory/corum_complexes.csv", stringsAsFactors = FALSE)
  sz <- lengths(cr$members); nm <- names(cr$members); rp <- nm %in% cx$complex[cx$group == "RP-dominated"]
  stopifnot(sum(rp) == 8L, all(sz == cx$size[match(nm, cx$complex)]))
  sub_of <- function(keep) list(members = cr$members[keep], class = cr$class[keep])
  mk <- function(reg) list(reg = reg, pairs = comembership_pairs(reg, none, genes))
  bench <- list(DoRothEA_AB = mk(tr$DoRothEA_AB), CORUM_all = mk(cr), KD_response = mk(tr$KD_response), CORUM_noRP = mk(sub_of(!rp)),
                CORUM_RP = mk(sub_of(rp)), CORUM_small = mk(sub_of(sz >= 3 & sz <= 6)), CORUM_large = mk(sub_of(sz >= 20)), CollecTRI = mk(tr$CollecTRI))
  for (e in names(bench)) say(sprintf("truth %-12s %5d regulons, %8d pairs, base rate %.4f", e, length(bench[[e]]$reg$members), length(bench[[e]]$pairs$label), mean(bench[[e]]$pairs$label)))
  base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
  ## jackknife streams: seed 20261016, order of baselines.R / robust_z.R (DoRothEA A-B, A-C, CollecTRI, CORUM, KD; then 100 d20 gene subsets; then the extra CORUM subsets)
  set.seed(EI_SEED)
  gin_of <- function(reg) sort(unique(unlist(reg$members)))
  draw_half <- function(gin) lapply(seq_len(EI_JACK), function(r) sort(sample(gin, floor(length(gin) / 2))))
  halves <- list()
  halves$DoRothEA_AB <- draw_half(bench$DoRothEA_AB$pairs$genes_in)
  invisible(draw_half(gin_of(dorothea_sets(genes, c("A", "B", "C")))))
  halves$CollecTRI <- draw_half(gin_of(collectri_sets(genes)))
  halves$CORUM_all <- draw_half(bench$CORUM_all$pairs$genes_in); halves$KD_response <- draw_half(bench$KD_response$pairs$genes_in)
  for (r in seq_len(100L)) invisible(sample(G, G - round(0.2 * G)))     # robust_z's d20 gene subsets: advance the stream only
  for (e in c("CORUM_noRP", "CORUM_RP", "CORUM_small", "CORUM_large")) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
  list(X = X, genes = genes, G = G, lab = lab, bench = bench, halves = halves, base_rate = base_rate, UT = which(upper.tri(matrix(FALSE, G, G))),
       nc = max(1L, min(workers, parallel::detectCores() - 1L)))
}

## ---- scoring (robust_z.R: pooled / jack_one / evaluate, macro and drift dropped) ----------------------------------------------
ei_make_eval <- function(env, H_att, H_mi, say = function(...) cat(..., "\n")) {
  force(H_att); force(H_mi); G <- env$G; bench <- env$bench; halves <- env$halves; nc <- env$nc
  jack_one <- function(sc, b, keeps) {
    i <- b$pairs$i; j <- b$pairs$j; lab_ <- b$pairs$label
    unlist(parallel::mclapply(keeps, function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[i] & inh[j]); ap_fast(sc[ix], lab_[ix]) }, mc.cores = nc)) }
  pooled <- function(S, jack = TRUE) lapply(stats::setNames(names(bench), names(bench)), function(e) { b <- bench[[e]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]
    list(aupr = ap_fast(sc, b$pairs$label), jack = if (jack) jack_one(sc, b, halves[[e]]) else NULL) })
  function(S, hub = TRUE, jack = TRUE) {
    S[!is.finite(S)] <- 0; diag(S) <- 0; invisible(gc(verbose = FALSE))
    out <- list(pooled = pooled(S, jack))
    if (hub) { out$hub_att <- pooled(rz_hub_stratify(S, H_att), jack); invisible(gc(verbose = FALSE))
               out$hub_mi <- pooled(rz_hub_stratify(S, H_mi), jack); invisible(gc(verbose = FALSE)) }
    out }
}

## ---- K10 walk helpers ------------------------------------------------------------------------------------------------------
ei_attention <- function(A, t = EI_T, alpha = EI_ALPHA) {   # A: row-stochastic operator (dense or sparse); returns S_t = (P^t + P^t')/2
  P <- Matrix::Diagonal(nrow(A), 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- as.matrix(P); for (s in seq_len(t - 1L)) Pt <- as.matrix(Pt %*% P)
  (Pt + t(Pt)) / 2 }

## top-k positive-CLR partners per gene among the columns in `cols` (diagonal excluded); returns a G x k integer matrix (0 = none)
ei_topk <- function(C, cols, k = 20L) {
  G <- nrow(C); out <- matrix(0L, G, k); inU <- logical(G); inU[cols] <- TRUE
  for (i in seq_len(G)) { v <- C[i, ]; v[!inU] <- 0; v[i] <- 0; pos <- which(v > 0); if (!length(pos)) next
    o <- pos[order(v[pos], decreasing = TRUE, method = "radix")]; o <- o[seq_len(min(k, length(o)))]; out[i, seq_along(o)] <- o }
  out }
## hypergeometric chance rule -> soft weight rho; Tx, Ty: G x k matrices from ei_topk, n_univ = |U|
ei_rho <- function(Tx, Ty, n_univ, p_ref = 0.05) {
  G <- nrow(Tx); ov <- integer(G); p <- rep(1, G)
  kx <- rowSums(Tx > 0); ky <- rowSums(Ty > 0)
  for (i in seq_len(G)) { a <- Tx[i, Tx[i, ] > 0]; b <- Ty[i, Ty[i, ] > 0]; ov[i] <- length(intersect(a, b))
    if (kx[i] > 0 && ky[i] > 0 && ov[i] > 0) p[i] <- stats::phyper(ov[i] - 1L, kx[i], n_univ - 1L - kx[i], ky[i], lower.tail = FALSE) }
  list(overlap = ov, p = p, rho = pmin(1, log(p) / log(p_ref))) }
