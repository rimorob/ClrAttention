# DREAM5 pipeline, label-free part (docs/plans/dream5_design.md, FIXED 2026-11-04).
# Data loading, replicate groups, the seven methods, score caches with SHA-256 manifests, and the
# experiment-jackknife refit. Nothing here reads a gold standard or any scoring file.
#
# Conventions (all of them from the design):
#   * Input: official replicate-level net{n}_expression_data.tsv (chips x genes, header G1...),
#     no further normalisation; stored here as X = genes x chips (the package's orientation).
#   * Candidate edges: TF list x all genes, no self-edges; TF-TF pairs in both directions. Symmetric
#     methods give (a,b) and (b,a) the same score. A "TF-row matrix" M is |TF| x G, NA at self pairs.
#   * Replicate groups are runs of consecutive chips starting at Repeat == 1 (never #Experiment).
#   * Symmetric scores are cached as the column-major upper triangle (doubles, no rounding).
#
# Source from the repo root:  source("analysis/dream5_lib.R")

suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R")   # ap_fast
source("analysis/regulons.R")    # make_regulons, comembership_pairs

D5 <- list(
  synapse = "data/dream5/synapse", archive = "data/dream5/public_gnw/DREAM5_NetworkInferenceChallenge_AlternativeDataFormats",
  out = "results/dream5", seed_trees = 20261016L, ntrees = 1000L,
  alpha = 0.5, keff = 10, cap = 50L, t = 9L,
  refit_seed_base = 20263000L, n_refit = 100L,
  series_seed_base = 20264000L, n_series_refit = 100L,   # supplementary series-level jackknife (net4; design change log 2026-11-05, item 4)
  python = "~/venvs/clr-baselines/bin/python", julia = "~/.juliaup/bin/julia")
D5$nets <- c(1L, 3L, 4L)
# method codes; "sym" methods have a gene-gene cache, "tfrow" methods a TF-row (or directed all-gene) cache
D5$methods <- c(abs_pearson = "sym", clr = "sym", clrns = "sym", attention = "sym", pidc = "sym",
                genie3_tf = "tfrow", grnboost2_tf = "tfrow", genie3_all = "dir", grnboost2_all = "dir")
D5$primary <- c("abs_pearson", "clr", "clrns", "attention", "genie3_tf", "grnboost2_tf", "pidc")  # required by the truth guard
D5$refit_methods <- c("abs_pearson", "clr", "clrns", "attention")

## ---- data -------------------------------------------------------------------------------------
# Replicate group index from the Repeat column: a new group starts at every Repeat == 1.
d5_replicate_groups <- function(repeat_col) {
  r <- as.integer(repeat_col)
  stopifnot(!anyNA(r), length(r) > 0, r[1] == 1L)
  grp <- cumsum(r == 1L)
  # Repeat runs 1..k within every group
  ok <- unlist(lapply(split(r, grp), function(x) identical(x, seq_along(x))))
  if (!all(ok)) stop("Repeat is not 1..k within every replicate group")
  grp
}

d5_read_net <- function(net, input = c("replicate", "avg"), root = D5$synapse, archive = D5$archive) {
  input <- match.arg(input)
  d <- file.path(root, paste0("net", net), paste0("net", net, "_"))
  Xc <- as.matrix(utils::read.delim(paste0(d, "expression_data.tsv"), check.names = FALSE, colClasses = "numeric"))
  genes <- colnames(Xc)
  cf <- utils::read.delim(paste0(d, "chip_features.tsv"), check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  stopifnot(nrow(cf) == nrow(Xc), "Repeat" %in% names(cf))
  tfs <- trimws(readLines(paste0(d, "transcription_factors.tsv"))); tfs <- tfs[nzchar(tfs)]
  gid <- utils::read.delim(paste0(d, "gene_ids.tsv"), check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
  stopifnot(identical(gid[[1]], genes), all(tfs %in% genes), !anyDuplicated(genes), !anyDuplicated(tfs))
  grp <- d5_replicate_groups(cf$Repeat)
  X <- t(Xc); colnames(X) <- paste0("c", seq_len(ncol(X)))
  if (input == "avg") {   # the organisers' replicate-averaged file (sensitivity); layout groups x genes
    Xa <- as.matrix(utils::read.delim(file.path(archive, paste0("net", net), paste0("net", net, "_expression_data_avg.tsv")),
                                      check.names = FALSE, colClasses = "numeric"))
    stopifnot(identical(colnames(Xa), genes), nrow(Xa) == max(grp))
    X <- t(Xa); colnames(X) <- paste0("g", seq_len(ncol(X))); grp <- seq_len(ncol(X))
  }
  rownames(X) <- genes
  zv <- apply(X, 1, function(x) diff(range(x)) == 0)
  if (any(zv)) stop(sum(zv), " zero-variance genes (the design specifies 0 in all three networks)")
  # series: the coarse `#Experiment` index of each chip (107, 107, 59 values); used only by the supplementary series-level
  # jackknife, never by a method. Not defined for the replicate-averaged input.
  series <- if (input == "replicate" && "#Experiment" %in% names(cf)) as.integer(cf[["#Experiment"]]) else NULL
  list(net = net, input = input, X = X, genes = genes, tf_ids = tfs, tf_idx = match(tfs, genes),
       group = grp, n_groups = max(grp), series = series, gene_names = gid[[2]], chip_features = cf)
}

# TF-row matrix (|TF| x G, NA on self pairs) template
d5_selfmask <- function(nd) cbind(seq_along(nd$tf_idx), nd$tf_idx)

## ---- random halves ---------------------------------------------------------------------------
.with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()))
  set.seed(seed); force(expr)
}
# SE_B replicate b: keep a random half of the replicate groups; seed 20263000 + b (design, Statistics)
d5_group_half <- function(n_groups, b) .with_seed(D5$refit_seed_base + as.integer(b), sort(sample.int(n_groups, floor(n_groups / 2))))
d5_subset_groups <- function(nd, keep_groups) {
  chips <- which(nd$group %in% keep_groups)
  nd$X <- nd$X[, chips, drop = FALSE]; nd$group <- nd$group[chips]; nd
}
# Supplementary series-level jackknife replicate b: keep a random half (floor(n/2)) of the `#Experiment` series, all chips of a
# kept series; seed 20264000 + b. Returns indices into sort(unique(series)).
d5_series_half <- function(n_series, b) .with_seed(D5$series_seed_base + as.integer(b), sort(sample.int(n_series, floor(n_series / 2))))
d5_subset_series <- function(nd, b) {
  stopifnot(!is.null(nd$series)); ids <- sort(unique(nd$series)); keep <- ids[d5_series_half(length(ids), b)]
  chips <- which(nd$series %in% keep)
  nd$X <- nd$X[, chips, drop = FALSE]; nd$group <- nd$group[chips]; nd$series <- nd$series[chips]; nd$n_series_kept <- length(keep); nd
}

## ---- symmetric-matrix packing ------------------------------------------------------------------
d5_clean_sym <- function(S) { S[!is.finite(S)] <- 0; diag(S) <- 0; S }
d5_pack <- function(S) S[upper.tri(S)]
d5_unpack <- function(v, G) { S <- matrix(0, G, G); S[upper.tri(S)] <- v; S + t(S) }
# index of the pair (i, j), i != j, in the column-major upper triangle (vectorised)
d5_pair_index <- function(i, j) { lo <- pmin(i, j); hi <- pmax(i, j); lo + (as.numeric(hi) - 1) * (as.numeric(hi) - 2) / 2 }

## ---- methods (each takes X = genes x samples) --------------------------------------------------
d5_m_abs_pearson <- function(X) { r <- abs(stats::cor(t(X))); r[is.na(r)] <- 0; d5_clean_sym(r) }

# CLR and attention from one fit: exactly analysis/baselines.R (bins "hg", transform "none", normal + Stouffer;
# alpha 0.5 softmax keff 10 cap 50; lazy walk P = (1-alpha) I + alpha A; score (P^t + t(P^t)) / 2, t = 9)
d5_m_clr_attention <- function(X, workers = 16L, alpha = D5$alpha, t = D5$t, keff = D5$keff, cap = D5$cap, want = c("clr", "attention")) {
  G <- nrow(X)
  f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = workers)$calibrate(method = "normal", combine = "stouffer")
  out <- list()
  if ("clr" %in% want) out$clr <- d5_clean_sym(f$clr_scores)
  if ("attention" %in% want) {
    f$build_operator(alpha = alpha, softmax_keff = keff, softmax_cap = as.integer(cap))
    P <- Matrix::Diagonal(G, 1 - alpha) + alpha * Matrix::Matrix(f$operator, sparse = TRUE)
    Pt <- as.matrix(P); for (s in seq_len(t - 1L)) Pt <- as.matrix(Pt %*% P)
    out$attention <- d5_clean_sym((Pt + t(Pt)) / 2); rm(Pt, P)
  }
  f$release(); invisible(gc(verbose = FALSE)); out
}

# Blom normal scores with average ranks per gene (identical to blom_avg in analysis/beeline.R; a test checks the text)
d5_blom_avg <- function(X) {
  n <- ncol(X)
  t(apply(X, 1, function(x) stats::qnorm((rank(x, ties.method = "average") - 3 / 8) / (n + 1 / 4))))
}
# CLR-NS: |r| of the normal scores, then CLR (normal + Stouffer), as ns_clr in analysis/beeline.R
d5_m_clrns <- function(X) {
  r <- abs(stats::cor(t(d5_blom_avg(X)))); r[is.na(r)] <- 0; diag(r) <- 0
  d5_clean_sym(clr_calibrate(r, method = "normal", combine = "stouffer"))
}

# GENIE3: R package, RF, K = sqrt, 1,000 trees, regulators = TF list (|TF| x G, NA on self) or all genes (G x G)
d5_m_genie3 <- function(X, tf_ids = NULL, workers = 16L, seed = D5$seed_trees, ntrees = D5$ntrees) {
  # GENIE3 seeds its trees through doRNG only when nCores >= 2 (then the result does not depend on the number of
  # cores); with nCores = 1 it takes another code path and gives different trees. All runs here use >= 2.
  if (workers < 2L) stop("GENIE3 needs workers >= 2 for the seeded (doRNG) code path")
  set.seed(seed)
  W <- GENIE3::GENIE3(X, regulators = tf_ids, targets = NULL, treeMethod = "RF", K = "sqrt",
                      nTrees = ntrees, nCores = workers, verbose = FALSE)
  genes <- rownames(X); W <- W[if (is.null(tf_ids)) genes else tf_ids, genes, drop = FALSE]
  ri <- match(rownames(W), genes); W[cbind(seq_len(nrow(W)), ri)] <- NA
  W
}

# GRNBoost2: tools/grnboost2.py (legacy early stopping, seed 20261016); writes the expression file, runs python,
# returns TF-row matrix (regulators = TF list) or G x G directed matrix (all genes).
d5_expr_tsv <- function(X, path) {
  if (!file.exists(path)) {
    tmp <- paste0(path, ".tmp", Sys.getpid())
    utils::write.table(data.frame(gene = rownames(X), X, check.names = FALSE), tmp, sep = "\t", quote = FALSE, row.names = FALSE)
    file.rename(tmp, path)
  }
  invisible(path)
}
d5_edges_to_matrix <- function(a, b, w, rows, cols) {
  M <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  ia <- match(a, rows); ib <- match(b, cols); ok <- !is.na(ia) & !is.na(ib)
  M[cbind(ia[ok], ib[ok])] <- w[ok]; M
}
d5_m_grnboost2 <- function(X, tsv, edges_out, tf_ids = NULL, workers = 16L, seed = D5$seed_trees, python = D5$python) {
  d5_expr_tsv(X, tsv)
  args <- c("tools/grnboost2.py", tsv, edges_out, "--workers", workers, "--seed", seed)
  if (!is.null(tf_ids)) { rf <- paste0(edges_out, ".regulators.txt"); writeLines(tf_ids, rf); args <- c(args, "--regulators", rf) }
  st <- system2(path.expand(python), args, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) stop("grnboost2.py failed: ", paste(utils::tail(st, 5), collapse = " | "))
  message("grnboost2: ", utils::tail(st, 1))
  e <- utils::read.delim(edges_out, stringsAsFactors = FALSE)
  genes <- rownames(X)
  M <- d5_edges_to_matrix(e$TF, e$target, e$importance, if (is.null(tf_ids)) genes else tf_ids, genes)
  ri <- match(rownames(M), genes); M[cbind(seq_len(nrow(M)), ri)] <- NA
  M
}

# PIDC: NetworkInference.jl, Bayesian blocks, all genes. The Julia process is launched by the job script (not from R:
# NetworkInference.jl on Julia 1.6 segfaults when started from R); this reads its edge file (gene1 gene2 weight, no header).
d5_m_pidc_read <- function(edges_file, genes) {
  e <- utils::read.delim(edges_file, header = FALSE, stringsAsFactors = FALSE)
  M <- matrix(0, length(genes), length(genes))
  ia <- match(e[[1]], genes); ib <- match(e[[2]], genes); ok <- !is.na(ia) & !is.na(ib)
  M[cbind(ia[ok], ib[ok])] <- e[[3]][ok]
  M <- pmax(M, t(M))     # undirected: one weight per pair
  d5_clean_sym(M)
}

## ---- caches and SHA-256 manifests -------------------------------------------------------------
d5_sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
d5_sha256_raw <- function(x) digest::digest(x, algo = "sha256", serialize = FALSE)
# deterministic binary: int32 header (ndim, dims...) then raw doubles
d5_write_scores <- function(x, path) {
  dims <- if (is.null(dim(x))) length(x) else dim(x)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  con <- file(tmp, "wb"); writeBin(c(length(dims), dims), con, size = 4L); writeBin(as.numeric(x), con, size = 8L); close(con)
  file.rename(tmp, path); invisible(path)
}
d5_read_scores <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  nd <- readBin(con, "integer", 1L, size = 4L); dims <- readBin(con, "integer", nd, size = 4L)
  x <- readBin(con, "double", prod(as.numeric(dims)), size = 8L)
  if (nd > 1L) dim(x) <- dims
  x
}
d5_dir <- function(net, out = D5$out) file.path(out, paste0("net", net))
d5_cache_path <- function(net, method, input = "replicate", out = D5$out)
  file.path(d5_dir(net, out), sprintf("%s__%s.scores", method, input))
# sidecar manifest (one per score file; concurrent jobs never share a file); atomic write
d5_sidecar_write <- function(path, meta) {
  row <- data.frame(net = meta$net, method = meta$method, input = meta$input, file = path, bytes = file.size(path),
                    sha256 = d5_sha256(path), dims = paste(meta$dims, collapse = "x"), seconds = round(meta$seconds, 1),
                    threads = meta$threads, finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                    stringsAsFactors = FALSE)
  tmp <- paste0(path, ".sha256.tmp", Sys.getpid())
  utils::write.table(row, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, paste0(path, ".sha256"))
  invisible(row)
}
# consolidated manifest: sort by net/method/input, atomic write. Includes only score files (not refit hash sidecars).
d5_collect_manifest <- function(out = D5$out, file = file.path(out, "SCORE_MANIFEST.tsv")) {
  sc <- Sys.glob(file.path(out, "net*", "*.scores.sha256"))
  rows <- do.call(rbind, lapply(sc, utils::read.delim, stringsAsFactors = FALSE, colClasses = c(sha256 = "character", dims = "character")))
  if (is.null(rows)) stop("no score sidecars under ", out)
  rows <- rows[order(rows$net, rows$method, rows$input), ]
  tmp <- paste0(file, ".tmp", Sys.getpid())
  utils::write.table(rows, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, file)
  invisible(rows)
}

# compute (or reuse) one method's primary cache for a network; returns the file path
d5_run_method <- function(nd, method, input = "replicate", workers = 16L, out = D5$out, force = FALSE) {
  stopifnot(method %in% names(D5$methods))
  net <- nd$net; dir.create(d5_dir(net, out), recursive = TRUE, showWarnings = FALSE)
  path <- d5_cache_path(net, method, input, out)
  if (file.exists(path) && file.exists(paste0(path, ".sha256")) && !force) { message(method, ": reusing ", path); return(path) }
  X <- nd$X; G <- nrow(X); t0 <- Sys.time()
  edges_dir <- d5_dir(net, out)
  res <- switch(method,
    abs_pearson = d5_pack(d5_m_abs_pearson(X)),
    clr = d5_pack(d5_m_clr_attention(X, workers, want = "clr")$clr),
    attention = d5_pack(d5_m_clr_attention(X, workers, want = "attention")$attention),
    clrns = d5_pack(d5_m_clrns(X)),
    genie3_tf = d5_m_genie3(X, nd$tf_ids, workers),
    genie3_all = d5_m_genie3(X, NULL, workers),
    grnboost2_tf = d5_m_grnboost2(X, file.path(edges_dir, sprintf("expression__%s.tsv", input)),
                                  file.path(edges_dir, sprintf("grnboost2_tf__%s_edges.tsv", input)), nd$tf_ids, workers),
    grnboost2_all = d5_m_grnboost2(X, file.path(edges_dir, sprintf("expression__%s.tsv", input)),
                                   file.path(edges_dir, sprintf("grnboost2_all__%s_edges.tsv", input)), NULL, workers),
    pidc = d5_pack(d5_m_pidc_read(file.path(edges_dir, sprintf("pidc__%s_edges.tsv", input)), nd$genes)))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d5_write_scores(res, path)
  d5_sidecar_write(path, list(net = net, method = method, input = input, dims = if (is.null(dim(res))) length(res) else dim(res),
                              seconds = secs, threads = workers))
  message(sprintf("%s net%d: %.1f s -> %s", method, net, secs, path)); path
}

## ---- experiment-jackknife refits ----------------------------------------------------------------
# One refit replicate: CLR, attention, |r| and CLR-NS refitted on a random half of the replicate groups (seed
# 20263000 + b). Returns the upper-triangle vectors. The truth and gene universe are fixed (not touched here).
# unit = "group": SE_B (random half of the replicate groups, seed 20263000 + b); unit = "series": the supplementary
# series-level jackknife (random half of the `#Experiment` series, seed 20264000 + b).
d5_refit_scores <- function(nd, b, workers = 16L, methods = D5$refit_methods, unit = c("group", "series")) {
  unit <- match.arg(unit)
  nb <- if (unit == "group") d5_subset_groups(nd, d5_group_half(nd$n_groups, b)) else d5_subset_series(nd, b); X <- nb$X
  zv <- apply(X, 1, function(x) diff(range(x)) == 0)
  if (any(zv)) stop("replicate ", b, ": ", sum(zv), " zero-variance genes in the kept half")
  out <- list()
  ca <- intersect(c("clr", "attention"), methods)
  if (length(ca)) { f <- d5_m_clr_attention(X, workers, want = ca); for (m in ca) out[[m]] <- d5_pack(f[[m]]); rm(f) }
  if ("abs_pearson" %in% methods) out$abs_pearson <- d5_pack(d5_m_abs_pearson(X))
  if ("clrns" %in% methods) out$clrns <- d5_pack(d5_m_clrns(X))
  out[methods]
}
d5_refit_dir <- function(net, out = D5$out, unit = c("group", "series")) file.path(d5_dir(net, out), if (match.arg(unit) == "group") "refit" else "refit_series")
d5_refit_path <- function(net, b, out = D5$out, unit = "group") file.path(d5_refit_dir(net, out, unit), sprintf("rep%03d.sha256.tsv", b))
d5_refit_manifest_file <- function(net, out = D5$out, unit = "group") file.path(d5_refit_dir(net, out, unit), if (unit == "group") "REFIT_MANIFEST.tsv" else "REFIT_SERIES_MANIFEST.tsv")
d5_hash_scores <- function(v) { con <- rawConnection(raw(0), "wb"); on.exit(close(con)); writeBin(v, con, size = 8L); d5_sha256_raw(rawConnectionValue(con)) }
# one replicate: compute, write the hash sidecar (one row per method). Resumable: an existing sidecar is kept.
# keep = "none" stores only the hashes (the design's storage ceiling does not allow 100 x 4 gene-gene matrices; the
# truth stage recomputes the refit and verifies the hashes before scoring). keep = "scores" also stores the vectors.
d5_refit_run <- function(nd, b, workers = 16L, out = D5$out, methods = D5$refit_methods, keep = c("none", "scores"), unit = c("group", "series")) {
  keep <- match.arg(keep); unit <- match.arg(unit); net <- nd$net
  dir.create(d5_refit_dir(net, out, unit), recursive = TRUE, showWarnings = FALSE)
  p <- d5_refit_path(net, b, out, unit)
  if (file.exists(p)) return(invisible(utils::read.delim(p, stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))))
  t0 <- Sys.time(); sc <- d5_refit_scores(nd, b, workers, methods, unit)
  kept <- if (unit == "group") length(d5_group_half(nd$n_groups, b)) else length(d5_series_half(length(unique(nd$series)), b))
  rows <- data.frame(net = net, rep = b, method = names(sc), kept = kept,
                     sha256 = vapply(sc, d5_hash_scores, ""), mean = vapply(sc, mean, 0), sd = vapply(sc, stats::sd, 0),
                     threads = workers, seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), stringsAsFactors = FALSE)
  names(rows)[names(rows) == "kept"] <- if (unit == "group") "n_groups_kept" else "n_series_kept"
  if (keep == "scores") for (m in names(sc)) d5_write_scores(sc[[m]], file.path(d5_refit_dir(net, out, unit), sprintf("rep%03d_%s.scores", b, m)))
  tmp <- paste0(p, ".tmp", Sys.getpid()); utils::write.table(rows, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, p)
  invisible(rows)
}
d5_refit_manifest <- function(net, out = D5$out, unit = c("group", "series"), file = d5_refit_manifest_file(net, out, match.arg(unit))) {
  fs <- sort(Sys.glob(file.path(d5_refit_dir(net, out, match.arg(unit)), "rep*.sha256.tsv")))
  rows <- do.call(rbind, lapply(fs, utils::read.delim, stringsAsFactors = FALSE, colClasses = c(sha256 = "character")))
  tmp <- paste0(file, ".tmp", Sys.getpid()); utils::write.table(rows, tmp, sep = "\t", quote = FALSE, row.names = FALSE); file.rename(tmp, file)
  invisible(rows)
}

## ---- synthetic data for tests and smoke runs ----------------------------------------------------
# A small network-shaped dataset with planted TF regulation and replicate groups; never real data.
d5_synth_net <- function(G = 60L, nTF = 8L, n_groups = 40L, seed = 1L, net = 0L) {
  .with_seed(seed, {
    sizes <- sample(1:3, n_groups, replace = TRUE, prob = c(.5, .3, .2)); n <- sum(sizes)
    grp <- rep(seq_len(n_groups), sizes); rep_col <- unlist(lapply(sizes, seq_len))
    genes <- paste0("G", seq_len(G)); tf_idx <- seq_len(nTF)
    Z <- matrix(stats::rnorm(G * n_groups), G, n_groups)       # group-level latent state (replicates share it)
    X <- Z[, grp] + matrix(stats::rnorm(G * n, sd = 0.5), G, n)
    for (t in tf_idx) { tg <- sample(setdiff(seq_len(G), t), 6); for (g in tg) X[g, ] <- X[g, ] + 0.9 * sample(c(-1, 1), 1) * X[t, ] }
    rownames(X) <- genes; colnames(X) <- paste0("c", seq_len(n))
    list(net = net, input = "replicate", X = X, genes = genes, tf_ids = genes[tf_idx], tf_idx = tf_idx,
         group = grp, n_groups = n_groups, series = (grp - 1L) %/% 3L + 1L, gene_names = genes, chip_features = data.frame(Repeat = rep_col))
  })
}
