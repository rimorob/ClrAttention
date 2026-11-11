#!/usr/bin/env Rscript
# Depth on a GEOMETRIC grid (log-scale re-test of analysis/granularity.R).
# Preregistration: docs/preregistration_depth_log.md (commit 8861ae4). Same truths, seeds and operators as
# analysis/granularity.R (which and whose outputs are untouched); only the depth grid and the module bootstrap differ.
#   grid t in {1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256} plus t = 0 (CLR). P^t is built by repeated
#   multiplication by P (t = 1..256); only grid points are scored.
#   Equivalence gate: at t in {0,1,2,4,8,16,32} the AUPR, base rate and coherence must equal
#   results/granularity/granularity_curves.csv to 1e-9, or the script stops.
#   Module bootstrap (secondary statistic): per truth B resamples of the modules; positive pair weight =
#   (sum of multiplicities of the resampled modules containing the pair) / (number of modules containing it);
#   positives with weight 0 dropped; negatives weight 1; weighted AP with the ap_fast tie convention.
#
# Usage (repo root):  Rscript analysis/granularity_log.R --org ecoli|k562 [--out results/granularity_log]
#   [--grid 1,2,...] [--truths "1 operons,3 sigmulons"] [--boot 1000] [--cores 8] [--check results/granularity/granularity_curves.csv]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(org = "ecoli", out = "results/granularity_log", grid = "1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256",
            truths = "", boot = 1000L, cores = 8L, alpha = 0.5, seed = 20261015L, p_max = 1e-3,
            check = "results/granularity/granularity_curves.csv", check_depths = "0,1,2,4,8,16,32",
            h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad")
i <- 1L
while (i <= length(args)) { k <- sub("^--", "", args[i]); opt[[k]] <- if (k %in% c("boot", "cores", "seed")) as.integer(args[i + 1L]) else args[i + 1L]; i <- i + 2L }
grid <- as.integer(strsplit(opt$grid, ",")[[1]]); tmax <- max(grid)
want <- if (nzchar(opt$truths)) strsplit(opt$truths, ",")[[1]] else NULL
stopifnot(opt$org %in% c("ecoli", "k562"))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, sprintf("granularity_log_%s.log", opt$org))
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
say("start: org ", opt$org, ", grid ", opt$grid, ", boot ", opt$boot, ", cores ", opt$cores, if (!is.null(want)) paste0(", truths ", paste(want, collapse = " | ")))

soft <- function(S) { A <- clr:::.softmax_rows(S, keff = 10, cap = 50L); e <- which(rowSums(A) == 0)
                      A[cbind(e, e)] <- 1; A }
sel <- function(bench) if (is.null(want)) bench else bench[names(bench) %in% want]

# score every truth at t = 0 and at each grid depth; keep the tie-aware ordering for the bootstrap
sweep <- function(C, bench, genes, label) {
  G <- nrow(C); rows <- list(); ordl <- lapply(bench, function(b) list()); lastl <- ordl
  scoreS <- function(S, t) for (nm in names(bench)) { b <- bench[[nm]]; p <- b$pairs; sc <- S[cbind(p$i, p$j)]
    rows[[length(rows) + 1L]] <<- data.frame(data = label, truth = nm, t = t, aupr = ap_fast(sc, p$label), base_rate = mean(p$label),
      coherence = stats::median(regulon_coherence(S, b$reg, b$co, b$operon, genes), na.rm = TRUE))
    if (opt$boot > 0L) { sc[is.na(sc)] <- -Inf; o <- order(sc, decreasing = TRUE, method = "radix"); s <- sc[o]
      ordl[[nm]][[as.character(t)]] <<- o; lastl[[nm]][[as.character(t)]] <<- c(s[-1L] != s[-length(s)], TRUE) } }
  scoreS(C, 0L)
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(soft(C), sparse = TRUE); Pt <- NULL
  for (t in seq_len(tmax)) { Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% grid) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; scoreS(S, t); say(sprintf("%s t=%d scored", label, t)) }
    else if (t %% 16L == 0L) say(sprintf("%s t=%d (walk only)", label, t)) }
  list(curves = do.call(rbind, rows), ord = ordl, last = lastl)
}

# equivalence gate against the saved linear-grid curves
check_equiv <- function(curves) {
  old <- utils::read.csv(opt$check, check.names = FALSE); td <- as.integer(strsplit(opt$check_depths, ",")[[1]])
  td <- intersect(td, c(0L, grid)); m <- merge(curves[curves$t %in% td, ], old, by = c("data", "truth", "t"), suffixes = c(".new", ".old"))
  stopifnot(nrow(m) == sum(curves$t %in% td))
  d <- c(aupr = max(abs(m$aupr.new - m$aupr.old)), base_rate = max(abs(m$base_rate.new - m$base_rate.old)),
         coherence = max(abs(m$coherence.new - m$coherence.old)))
  say(sprintf("EQUIVALENCE vs %s: %d rows (depths %s), max |diff| AUPR %.3g, base rate %.3g, coherence %.3g",
              opt$check, nrow(m), paste(td, collapse = ","), d[1], d[2], d[3]))
  writeLines(sprintf("rows %d depths %s max_abs_diff aupr %.3e base_rate %.3e coherence %.3e tolerance 1e-9 %s", nrow(m), paste(td, collapse = ","),
                     d[1], d[2], d[3], if (all(d <= 1e-9)) "PASS" else "FAIL"), file.path(opt$out, sprintf("equivalence_%s.txt", opt$org)))
  if (!all(d <= 1e-9)) stop("equivalence check FAILED (tolerance 1e-9)")
}

# weighted average precision on a stored ordering (same tie convention as ap_fast)
wap <- function(o, last, wpos, wall) {
  x <- cumsum(wpos[o]); tp <- x[last]; k <- cumsum(wall[o])[last]
  sum(diff(c(0, tp)) * tp / pmax(k, 1e-300)) / sum(wpos)   # k = 0 only before the first undropped pair (then tp = 0)
}

# module bootstrap of t_best for every truth
bootstrap <- function(bench, genes, st) {
  G <- length(genes); tv <- c(0L, grid); set.seed(opt$seed)
  sizes <- list(); out <- list()
  for (nm in names(bench)) {
    b <- bench[[nm]]; p <- b$pairs; npr <- length(p$i); mem <- lapply(b$reg$members, function(m) sort(unique(m))); nm_mod <- length(mem)
    msz <- lengths(mem)
    idx <- matrix(0L, G, G); idx[cbind(p$i, p$j)] <- seq_len(npr)
    qi <- vector("list", nm_mod); qj <- qi
    for (k in seq_len(nm_mod)) { m <- mem[[k]]; if (length(m) < 2L) next
      u <- which(upper.tri(diag(length(m))), arr.ind = TRUE); pid <- idx[cbind(m[u[, 1]], m[u[, 2]])]; pid <- pid[pid > 0L]
      qi[[k]] <- pid; qj[[k]] <- rep.int(k, length(pid)) }
    rm(idx); Q <- Matrix::sparseMatrix(i = unlist(qi), j = unlist(qj), x = 1, dims = c(npr, nm_mod)); rm(qi, qj)
    cnt <- as.numeric(Matrix::rowSums(Q)); lab <- p$label
    stopifnot(all(cnt[lab == 1L] >= 1), all(cnt[lab == 0L] == 0))
    wfun <- function(mk) { w <- rep(1, npr); pos <- lab == 1L; w[pos] <- as.numeric(Q %*% mk)[pos] / cnt[pos]; w }
    # assertion: m_k = 1 for all modules reproduces ap_fast exactly
    w1 <- wfun(rep(1, nm_mod)); wp1 <- w1 * lab
    for (tt in c(0L, tail(grid, 1L))) { a <- wap(st$ord[[nm]][[as.character(tt)]], st$last[[nm]][[as.character(tt)]], wp1, w1)
      ref <- st$curves$aupr[st$curves$truth == nm & st$curves$t == tt]
      if (abs(a - ref) > 1e-12) stop(sprintf("weighted AP at m=1 differs from ap_fast: %s t=%d %.15g vs %.15g", nm, tt, a, ref)) }
    say(sprintf("bootstrap %s: %d modules, %d pairs, weighted-AP identity at m=1 OK", nm, nm_mod, npr))
    dr <- matrix(sample.int(nm_mod, nm_mod * opt$boot, replace = TRUE), nm_mod, opt$boot)
    ordn <- st$ord[[nm]]; lastn <- st$last[[nm]]
    res <- parallel::mclapply(seq_len(opt$boot), function(d) { mk <- tabulate(dr[, d], nm_mod); w <- wfun(mk); wp <- w * lab
      a <- vapply(as.character(tv), function(tt) wap(ordn[[tt]], lastn[[tt]], wp, w), numeric(1))
      c(tv[which.max(a)], stats::median(msz[dr[, d]])) }, mc.cores = opt$cores)
    bad <- which(vapply(res, function(z) length(z) != 2L, TRUE)); if (length(bad)) stop("a bootstrap worker failed for ", nm, ": ", paste(capture.output(str(res[[bad[1]]])), collapse = " "), " (", length(bad), " bad of ", length(res), ")")
    r <- do.call(rbind, res)
    out[[nm]] <- data.frame(truth = nm, draw = seq_len(opt$boot), t_best = r[, 1], median_size = r[, 2])
    sizes[[nm]] <- data.frame(truth = nm, n_modules = nm_mod, median_size = stats::median(msz), mean_size = mean(msz), n_pairs = npr,
                              base_rate = mean(lab), n_genes_universe = length(p$genes_in))
    say(sprintf("bootstrap %s done: t_best median %g, 2.5%% %g, 97.5%% %g", nm, stats::median(r[, 1]), stats::quantile(r[, 1], .025), stats::quantile(r[, 1], .975)))
    rm(Q, dr, res, r); invisible(gc())
  }
  list(draws = do.call(rbind, out), sizes = do.call(rbind, sizes))
}

finish <- function(st, bench, genes) {
  utils::write.csv(st$curves, file.path(opt$out, sprintf("curves_%s.csv", opt$org)), row.names = FALSE)
  check_equiv(st$curves)
  if (opt$boot > 0L) { bs <- bootstrap(bench, genes, st)
    gz <- gzfile(file.path(opt$out, sprintf("boot_%s.csv.gz", opt$org)), "w"); utils::write.csv(bs$draws, gz, row.names = FALSE); close(gz)
    utils::write.csv(bs$sizes, file.path(opt$out, sprintf("truth_sizes_%s.csv", opt$org)), row.names = FALSE) }
  say("done")
}

if (opt$org == "ecoli") {
## ---- E. coli ------------------------------------------------------------------------
m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X); G <- length(genes)
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
reg_all <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
pick <- function(keep) list(members = reg_all$members[keep], class = reg_all$class[keep])
tf_reg <- pick(reg_all$class == "TF" & lengths(reg_all$members) <= 500 & lengths(reg_all$members) >= 5)
# Sigmulons without the two global ones: sigma70 (housekeeping, 1,596 genes, 37%) and
# sigmaS/sigma38 (general stress, 331 genes): sigma32, sigma54, sigma24, sigma28, sigma19.
sig_reg <- pick(reg_all$class == "sigma" & lengths(reg_all$members) >= 5 &
                  !names(reg_all$members) %in% c("sigma:sigma70", "sigma:sigma38"))
op <- operon_of[names(operon_of) %in% genes]
op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
op_reg <- list(members = op_sets, class = stats::setNames(rep("operon", length(op_sets)), names(op_sets)))
none <- stats::setNames(character(0), character(0))
mk <- function(r, o) list(reg = r, pairs = comembership_pairs(r, o, genes), co = comember_matrix(r, G), operon = o)
bench_e <- list(`1 operons` = mk(op_reg, none), `2 TF regulons` = mk(tf_reg, operon_of),
                `3 sigmulons` = mk(sig_reg, operon_of))
for (nm in names(bench_e)) say(sprintf("E. coli %s: %d sets, %d pairs, base rate %.4f", nm,
  length(bench_e[[nm]]$reg$members), length(bench_e[[nm]]$pairs$label), mean(bench_e[[nm]]$pairs$label)))
bench_e <- sel(bench_e)
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
st <- sweep(f$clr_scores, bench_e, genes, "E. coli M3D"); f$release(); rm(f, X); invisible(gc())
finish(st, bench_e, genes)

} else {
## ---- K562 compendium -------------------------------------------------------------------
none <- stats::setNames(character(0), character(0))
set.seed(opt$seed)                       # same split as perturbseq_bulk.R
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
Z <- rd("X"); hg <- cat_of("/var", "gene_name"); rownames(Z) <- hg
pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
pval <- as.numeric(rd("/obs/energy_test_p_value"))
K <- pert[pval < opt$p_max & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
Xh <- Z[, A]; rm(Z); Xh[!is.finite(Xh)] <- 0
Xh <- Xh[apply(Xh, 1, function(x) diff(range(x))) > 0, ]; hgenes <- rownames(Xh); H <- length(hgenes)
mkh <- function(r) list(reg = r, pairs = comembership_pairs(r, none, hgenes), co = comember_matrix(r, H), operon = none)
mkl <- list(`1 CORUM complexes` = function() mkh(corum_sets(hgenes)),
            `2 KD-response modules` = function() mkh(response_sets(hgenes, path = opt$h5, use = B)),
            `3 DoRothEA A-B` = function() mkh(dorothea_sets(hgenes)),
            `4 DoRothEA A-C` = function() mkh(dorothea_sets(hgenes, c("A", "B", "C"))),
            `5 CollecTRI` = function() mkh(collectri_sets(hgenes)))
bench_h <- lapply(if (is.null(want)) mkl else mkl[names(mkl) %in% want], function(f) f())
for (nm in names(bench_h)) say(sprintf("K562 %s: %d sets, base rate %.4f", nm,
  length(bench_h[[nm]]$reg$members), mean(bench_h[[nm]]$pairs$label)))
fh <- ClrAttention$new(Xh)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
st <- sweep(fh$clr_scores, bench_h, hgenes, "K562 essential compendium")
finish(st, bench_h, hgenes)
}
