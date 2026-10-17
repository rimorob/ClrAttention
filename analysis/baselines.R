#!/usr/bin/env Rscript
# Competitor baselines on the same data and truths as CLR attention:
#   GENIE3     random-forest regression of each gene on all others (Huynh-Thu
#              et al. 2010; GENIE3 R package), DREAM5 winner class
#   GRNBoost2  gradient-boosting variant (Moerman et al. 2019; arboreto)
#   PIDC       partial information decomposition (Chan, Stumpf & Babtie 2017;
#              NetworkInference.jl), the information-theoretic state of the art
# Every gene is a candidate regulator, so each method yields a gene x gene
# score matrix; directed weights are symmetrized with max(w_ij, w_ji). The
# same matrices feed every benchmark: regulon co-membership AUPR and median
# coherence (E. coli: strong/confirmed and all evidence; K562: DoRothEA A-B,
# A-C, CollecTRI, CORUM, knockdown-response modules), and for E. coli the
# TF-node benchmark (pooled AUPR ratio, median per-TF AUROC).
# References computed in the same run: |Pearson r|, CLR (MI), and soft10
# attention on CLR at --t. Differences against attention and against CLR have
# paired delete-half gene jackknife intervals (AUPR).
#
# Datasets (same construction as the main analyses):
#   m3d     M3D, 466 replicate-averaged experiments
#   ecomac  EcoMAC, 2,198 arrays, ComBat by platform
#   k562ess Perturb-seq K562 essential pseudobulk, network half A of the
#           significant knockdowns (seed as perturbseq_bulk.R)
#
# Usage (repo root):
# --gene_subset truths (k562ess): every method is run on the genes covered by
#   at least one literature truth (DoRothEA A-C, CollecTRI, CORUM; 5,672 of
#   8,563), the same input for all methods. PIDC's cost grows with the cube of
#   the gene count (about 5.5 days single-threaded on all genes, about 40 h on
#   this subset), so the like-for-like comparison with PIDC uses this subset.
#
#   Rscript analysis/baselines.R --dataset m3d [--methods genie3,grnboost2,pidc]
#       [--t 9] [--workers 24] [--out results/baselines_m3d]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(dataset = "m3d", methods = "genie3,grnboost2,pidc", t = 9L, workers = 24L, out = NULL,
            jack = 100L, seed = 20261016L, alpha = 0.5, ntrees = 1000L, quick = 0L,
            python = "~/venvs/clr-baselines/bin/python", julia = "~/.juliaup/bin/julia",
            gene_subset = "none", export_only = 0L)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i]); if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  v <- args[i + 1L]
  opt[[key]] <- if (key %in% c("t", "workers", "jack", "seed", "ntrees", "quick", "export_only")) as.integer(v) else v
  i <- i + 2L
}
if (is.null(opt$out)) opt$out <- file.path("results", paste0("baselines_", opt$dataset))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "baselines.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=", collapse = " "))

## ---- data and truths -------------------------------------------------------------------
none <- stats::setNames(character(0), character(0))
if (opt$dataset %in% c("m3d", "ecomac")) {
  if (opt$dataset == "m3d") {
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
    raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
  } else {
    X <- as.matrix(utils::read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
    rownames(X) <- utils::read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)[[2]]
    plat <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)$platform
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    X <- suppressMessages(sva::ComBat(X, batch = factor(plat)))
  }
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  if (opt$quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ] }
  genes <- rownames(X); G <- length(genes)
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
  operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes), co = comember_matrix(reg, G), operon = operon_of)
  })
  tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
  tf_parts <- strsplit(.classic_tf_to_genes(tfr$regulator), ";", fixed = TRUE)
  tfnet <- data.frame(tf_gene = to_bn(unlist(tf_parts)), target = rep(tfr$bnumber, lengths(tf_parts)),
                      confidence = rep(tfr$confidence, lengths(tf_parts)), stringsAsFactors = FALSE)
  tfnet <- tfnet[!is.na(tfnet$tf_gene), ]
  U <- list(SC = edge_universe(tfnet, genes, conf_keep = c("S", "C")), all = edge_universe(tfnet, genes))
} else if (opt$dataset == "k562ess") {
  h5 <- "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"
  set.seed(20261015L)                                   # perturbseq_bulk.R split
  rd <- function(p) rhdf5::h5read(h5, p)
  cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
  Z <- rd("X"); rownames(Z) <- cat_of("/var", "gene_name"); pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
  pval <- as.numeric(rd("/obs/energy_test_p_value"))
  K <- pert[pval < 1e-3 & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
  A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
  X <- Z[, A]; rm(Z); X[!is.finite(X)] <- 0
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  dup <- rownames(X) %in% rownames(X)[duplicated(rownames(X))]   # ambiguous symbols (e.g. HSPA14, TBCE): drop all copies
  if (any(dup)) { say(sprintf("dropped %d rows with duplicated gene names: %s", sum(dup),
                              paste(unique(rownames(X)[dup]), collapse = ", "))); X <- X[!dup, ] }
  if (opt$gene_subset == "truths") {
    g0 <- rownames(X)
    lit <- Reduce(union, lapply(list(dorothea_sets(g0, c("A", "B", "C")), collectri_sets(g0), corum_sets(g0)),
                                function(s) unique(unlist(s$members))))
    X <- X[sort(lit), ]
  }
  if (opt$quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ] }
  genes <- rownames(X); G <- length(genes)
  tr <- list(DoRothEA_AB = dorothea_sets(genes), DoRothEA_ABC = dorothea_sets(genes, c("A", "B", "C")),
             CollecTRI = collectri_sets(genes), CORUM = corum_sets(genes),
             KD_response = response_sets(genes, path = h5, use = B))
  bench <- lapply(tr, function(r) list(reg = r, pairs = comembership_pairs(r, none, genes), co = comember_matrix(r, G), operon = none))
  U <- NULL
} else stop("unknown dataset: ", opt$dataset)
say(sprintf("%s: %d genes x %d samples", opt$dataset, G, ncol(X)))
# --export_only 1: write the input matrix and stop (the PIDC lane runs Julia
# directly on it, outside R, and the main lane scores the edges later).
tsv <- file.path(opt$out, "expression.tsv")
if (!file.exists(tsv)) utils::write.table(data.frame(gene = genes, X, check.names = FALSE), tsv, sep = "\t",
                                          quote = FALSE, row.names = FALSE)
if (opt$export_only == 1L) { say("exported ", tsv); quit(save = "no") }
set.seed(opt$seed)
halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r) {
  keep <- sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))
  inh <- logical(G); inh[keep] <- TRUE; which(inh[b$pairs$i] & inh[b$pairs$j]) }))

## ---- scoring ---------------------------------------------------------------------------
per_tf_median_auroc <- function(S, u) stats::median(vapply(names(u$targets), function(tn) {
  t <- as.integer(tn); y <- integer(G); y[unique(u$targets[[tn]])] <- 1L; y[t] <- NA
  if (sum(y, na.rm = TRUE) < 3) return(NA_real_)
  s <- S[t, ]; ok <- !is.na(y); r <- rank(s[ok]); n1 <- sum(y[ok]); n0 <- sum(ok) - n1
  (sum(r[y[ok] == 1]) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0) }, 0), na.rm = TRUE)
evaluate <- function(S) {
  S[!is.finite(S)] <- 0; diag(S) <- 0
  lapply(stats::setNames(names(bench), names(bench)), function(nm) {
    b <- bench[[nm]]; sc <- S[cbind(b$pairs$i, b$pairs$j)]
    out <- list(aupr = ap_fast(sc, b$pairs$label),
                jack = vapply(halves[[nm]], function(ix) ap_fast(sc[ix], b$pairs$label[ix]), 0),
                coherence = stats::median(regulon_coherence(S, b$reg, b$co, b$operon, genes), na.rm = TRUE))
    if (!is.null(U) && nm %in% names(U)) { u <- U[[nm]]
      out$tf_aupr_ratio <- ap_fast(S[cbind(u$i, u$j)], u$label) / mean(u$label)
      out$tf_auroc_median <- per_tf_median_auroc(S, u) }
    out })
}
sym_max <- function(W) { W[is.na(W)] <- 0; pmax(W, t(W)) }
edges_to_matrix <- function(a, b, w) { M <- matrix(0, G, G, dimnames = list(genes, genes))
  ia <- match(a, genes); ib <- match(b, genes); ok <- !is.na(ia) & !is.na(ib)
  M[cbind(ia[ok], ib[ok])] <- w[ok]; M }

## ---- references: |r|, CLR, attention ---------------------------------------------------
E <- list()
r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; E$abs_pearson <- evaluate(r); rm(r)
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$workers)$
  calibrate(method = "normal", combine = "stouffer")
E$clr <- evaluate(f$clr_scores)
f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(f$operator, sparse = TRUE)
Pt <- as.matrix(P); for (s in seq_len(opt$t - 1L)) Pt <- as.matrix(Pt %*% P)
E$attention <- evaluate((Pt + t(Pt)) / 2); rm(P, Pt); f$release(); invisible(gc())
say("references scored")

## ---- baselines ---------------------------------------------------------------------------
run_method <- function(m) {
  rdsf <- file.path(opt$out, paste0(m, "_scores.rds"))
  if (file.exists(rdsf)) { say(m, ": reusing ", rdsf); return(readRDS(rdsf)) }
  t0 <- Sys.time()
  S <- switch(m,
    genie3 = { set.seed(opt$seed)
      W <- GENIE3::GENIE3(X, regulators = NULL, targets = NULL, treeMethod = "RF", K = "sqrt",
                          nTrees = opt$ntrees, nCores = opt$workers, verbose = FALSE)
      W <- W[genes, genes]; sym_max(W) },
    grnboost2 = { outf <- file.path(opt$out, "grnboost2_edges.tsv")
      st <- system2(path.expand(opt$python), c("tools/grnboost2.py", tsv, outf, "--workers", opt$workers,
                                                "--seed", opt$seed), stdout = TRUE, stderr = TRUE)
      say("grnboost2: ", utils::tail(st, 1))
      e <- utils::read.delim(outf, stringsAsFactors = FALSE)
      sym_max(edges_to_matrix(e$TF, e$target, e$importance)) },
    pidc = { outf <- file.path(opt$out, "pidc_edges.tsv")
      if (!file.exists(outf)) {
        # Normally produced by the PIDC lane (Julia launched directly, not from
        # R: NetworkInference.jl on Julia 1.6 segfaults when started from R).
        stop("PIDC edges not found: ", outf, " (run tools/pidc.jl on ", tsv, " first)")
      }
      e <- utils::read.delim(outf, header = FALSE, stringsAsFactors = FALSE)
      sym_max(edges_to_matrix(e[[1]], e[[2]], e[[3]])) },
    stop("unknown method ", m))
  say(sprintf("%s computed in %.1f min", m, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(S, rdsf); S
}
for (m in strsplit(opt$methods, ",")[[1]]) {
  S <- tryCatch(run_method(m), error = function(e) { say(m, " FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(S)) E[[m]] <- evaluate(S)
  rm(S); invisible(gc())
}

## ---- summary -------------------------------------------------------------------------------
rows <- list()
for (m in names(E)) for (nm in names(bench)) {
  e <- E[[m]][[nm]]; da <- e$jack - E$attention[[nm]]$jack; dc <- e$jack - E$clr[[nm]]$jack
  d0a <- e$aupr - E$attention[[nm]]$aupr; d0c <- e$aupr - E$clr[[nm]]$aupr
  rows[[length(rows) + 1L]] <- data.frame(dataset = opt$dataset, method = m, truth = nm, aupr = e$aupr,
    base_rate = mean(bench[[nm]]$pairs$label), coherence = e$coherence,
    d_vs_attention = d0a, d_att_lo = d0a - 1.96 * stats::sd(da), d_att_hi = d0a + 1.96 * stats::sd(da),
    d_vs_clr = d0c, d_clr_lo = d0c - 1.96 * stats::sd(dc), d_clr_hi = d0c + 1.96 * stats::sd(dc),
    tf_aupr_ratio = if (is.null(e$tf_aupr_ratio)) NA else e$tf_aupr_ratio,
    tf_auroc_median = if (is.null(e$tf_auroc_median)) NA else e$tf_auroc_median)
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(opt$out, "metrics.csv"), row.names = FALSE)
for (nm in names(bench)) {
  x <- res[res$truth == nm, ]; x <- x[order(-x$aupr), ]
  say(sprintf("== %s %s (base rate %.4f)", opt$dataset, nm, x$base_rate[1]))
  for (k in seq_len(nrow(x))) say(sprintf("   %-12s AUPR %.4f  coherence %.3f  vs attention %+.4f [%+.4f, %+.4f]%s",
    x$method[k], x$aupr[k], x$coherence[k], x$d_vs_attention[k], x$d_att_lo[k], x$d_att_hi[k],
    if (is.na(x$tf_aupr_ratio[k])) "" else sprintf("  | TF-node AUPR x%.2f, per-TF AUROC %.3f", x$tf_aupr_ratio[k], x$tf_auroc_median[k])))
}
say("done")
