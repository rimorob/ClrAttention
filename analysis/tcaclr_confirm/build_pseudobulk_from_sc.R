#!/usr/bin/env Rscript
# Label-free reconstruction of the Replogle-style "normalized bulk" Perturb-seq pseudobulk from a deposited single-cell h5ad
# (docs/plans/tcaclr_confirm_plan.md, Arm A; PROPOSED deviation D1/D2, see docs/plans/tcaclr_confirm_freeze.md).
# The K562/RPE1 recipe (analysis/perturbseq_bulk.R) starts from the authors' normalized-bulk h5ad (genes x perturbations of mean
# z-normalised expression, plus obs/energy_test_p_value). GSE264667 HepG2 ships only the single-cell h5ad, so this script builds the same
# object. Definition (validated on K562 against the shipped file, see freeze document):
#   1. size normalisation: counts / UMI_count (obs) x median(UMI_count);
#   2. z-normalisation of every gene against the non-targeting (control) cells of the SAME gem group (mean and SD over those cells,
#      SD with --ddof); gem groups with fewer than --min_ctrl control cells use the pooled control mean and SD; SD 0 or non-finite -> z 0;
#   3. perturbation profile = mean z over the cells carrying that gene_transcript label (all labels, non-targeting groups included, as shipped);
#   4. energy-distance permutation test of each non-targeting-free perturbation against a fixed random sample of --nctrl control cells
#      in the first --npc principal components of the z-normalised cells (PCA on a --pca_cells random cell sample), --nperm permutations,
#      p = (1 + #{perm >= observed}) / (1 + nperm).
# No truth is read; nothing is scored. Output: an h5ad-like file readable by the unchanged loader of perturbseq_bulk.R
# (/X genes x perturbations, /obs/gene_transcript, /obs/energy_test_p_value, /obs/num_cells, /var/gene_name).
# Usage: Rscript build_pseudobulk_from_sc.R --h5 <raw sc h5ad> --out <bulk h5ad> [--threads 8] [--limit_cells 0]
args <- commandArgs(trailingOnly = TRUE)
opt <- list(h5 = NULL, out = NULL, threads = 8L, ddof = 1L, min_ctrl = 30L, npc = 50L, pca_cells = 25000L, nctrl = 1000L, nperm = 10000L,
            seed = 20261120L, chunk = 4000L, limit_cells = 0L, limit_perts = 0L, ref = "")
i <- 1L; while (i <= length(args)) { k <- sub("^--", "", args[i]); stopifnot(k %in% names(opt)); v <- args[i + 1L]
  opt[[k]] <- if (is.character(opt[[k]]) || is.null(opt[[k]])) v else if (is.integer(opt[[k]])) as.integer(v) else as.numeric(v); i <- i + 2L }
stopifnot(!is.null(opt$h5), !is.null(opt$out))
Sys.setenv(OMP_NUM_THREADS = as.character(opt$threads), OPENBLAS_NUM_THREADS = as.character(opt$threads))
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n") }
vmhwm <- function() { l <- readLines("/proc/self/status"); as.numeric(sub("^VmHWM:[[:space:]]*([0-9]+).*", "\\1", grep("^VmHWM:", l, value = TRUE))) / 1048576 }
suppressPackageStartupMessages({ library(rhdf5); library(Rcpp); library(parallel) })
cppsrc <- '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
// [[Rcpp::export]]
List energy_perm(NumericMatrix Y, int n, int nperm, int seed) {   // Y: d x (n + m) coordinates; first n columns = perturbed cells
  int N = Y.ncol(), d = Y.nrow(), m = N - n;
  std::vector<double> D((size_t)N * N, 0.0), s(N, 0.0);
  for (int a = 0; a < N; a++) for (int b = a + 1; b < N; b++) { double q = 0; for (int k = 0; k < d; k++) { double z = Y(k, a) - Y(k, b); q += z * z; }
    double v = std::sqrt(q); D[(size_t)a * N + b] = v; D[(size_t)b * N + a] = v; s[a] += v; s[b] += v; }
  double T = 0; for (int a = 0; a < N; a++) T += s[a];
  auto stat = [&](const std::vector<int>& I) { double W = 0, R = 0;
    for (int a = 0; a < n; a++) { R += s[I[a]]; const double* row = &D[(size_t)I[a] * N]; for (int b = 0; b < n; b++) W += row[I[b]]; }
    double cross = R - W, WY = T - 2 * R + W; return 2 * cross / ((double)n * m) - W / ((double)n * n) - WY / ((double)m * m); };
  std::vector<int> I(n); for (int a = 0; a < n; a++) I[a] = a; double obs = stat(I);
  std::vector<int> perm(N); for (int a = 0; a < N; a++) perm[a] = a;
  unsigned long long st = 88172645463325252ULL ^ (unsigned long long)seed * 2654435761ULL; auto rnd = [&]() { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; };
  int ge = 0;
  for (int it = 0; it < nperm; it++) { for (int a = 0; a < n; a++) { size_t j = a + rnd() % (N - a); std::swap(perm[a], perm[j]); }
    for (int a = 0; a < n; a++) I[a] = perm[a]; if (stat(I) >= obs - 1e-12) ge++; }
  return List::create(_["E"] = obs, _["p"] = (1.0 + ge) / (1.0 + nperm));
}'
sourceCpp(code = cppsrc, verbose = FALSE)
f <- opt$h5; t0 <- Sys.time()
cat_of <- function(grp, field) { x <- h5read(f, paste0(grp, "/", field)); if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(h5read(f, paste0(grp, "/__categories/", field)), error = function(e) NULL); if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
gt <- cat_of("/obs", "gene_transcript"); um <- as.numeric(h5read(f, "/obs/UMI_count")); gem <- as.integer(h5read(f, "/obs/gem_group")); genes <- cat_of("/var", "gene_name")
Ncell <- length(gt); G <- length(genes); keep_cells <- seq_len(Ncell)
if (opt$limit_cells > 0L) { set.seed(1); keep_cells <- sort(sample(Ncell, opt$limit_cells)) }   # smoke tests only
nt <- grepl("non-targeting", gt); scale0 <- stats::median(um)
say(sprintf("%d cells x %d genes; %d non-targeting cells; %d gem groups; UMI median %.0f", Ncell, G, sum(nt), length(unique(gem)), scale0))
readX <- function(idx) { x <- h5read(f, "X", index = list(NULL, idx)); sweep(x, 2, um[idx], "/") * scale0 }   # genes x cells, size-normalised
## ---- control statistics per gem group --------------------------------------------------------------------------------------
ntc <- intersect(which(nt), keep_cells); gem_ids <- sort(unique(gem)); S1 <- S2 <- matrix(0, G, length(gem_ids)); NC <- integer(length(gem_ids))
for (st in seq(1L, length(ntc), by = opt$chunk)) { idx <- ntc[st:min(st + opt$chunk - 1L, length(ntc))]; x <- readX(idx); gi <- match(gem[idx], gem_ids)
  for (g in unique(gi)) { k <- which(gi == g); S1[, g] <- S1[, g] + rowSums(x[, k, drop = FALSE]); S2[, g] <- S2[, g] + rowSums(x[, k, drop = FALSE]^2); NC[g] <- NC[g] + length(k) } }
sdev <- function(s1, s2, n) { v <- (s2 - s1^2 / n) / (n - opt$ddof); sqrt(pmax(v, 0)) }
Mu <- S1 / rep(pmax(NC, 1), each = G); Sd <- sdev(S1, S2, rep(NC, each = G)); pooled_mu <- rowSums(S1) / sum(NC); pooled_sd <- sdev(rowSums(S1), rowSums(S2), sum(NC))
small <- NC < opt$min_ctrl; Mu[, small] <- pooled_mu; Sd[, small] <- pooled_sd
say(sprintf("control cells per gem group: min %d, median %g, max %d; groups using pooled statistics (< %d controls): %d", min(NC), stats::median(NC), max(NC), opt$min_ctrl, sum(small)))
znorm <- function(x, idx) { gi <- match(gem[idx], gem_ids); z <- (x - Mu[, gi, drop = FALSE]) / Sd[, gi, drop = FALSE]; z[!is.finite(z)] <- 0; z }
## ---- PCA on a random cell sample (randomized SVD) ----------------------------------------------------------------------------
set.seed(opt$seed); ps <- sort(sample(keep_cells, min(opt$pca_cells, length(keep_cells))))
Zs <- matrix(0, G, length(ps)); for (st in seq(1L, length(ps), by = opt$chunk)) { k <- st:min(st + opt$chunk - 1L, length(ps)); Zs[, k] <- znorm(readX(ps[k]), ps[k]) }
cen <- rowMeans(Zs); Zs <- Zs - cen
Om <- matrix(stats::rnorm(length(ps) * 0 + length(ps) * (opt$npc + 20L)), length(ps), opt$npc + 20L)
Q <- qr.Q(qr(Zs %*% Om)); for (q in 1:4) { Q <- qr.Q(qr(Zs %*% crossprod(Zs, Q))) }
sv <- svd(crossprod(Q, Zs), nu = opt$npc, nv = 0); V <- Q %*% sv$u; rm(Zs, Om, Q)
say(sprintf("PCA on %d cells, %d components; singular values %.1f ... %.1f; %.0f s", length(ps), opt$npc, sv$d[1], sv$d[opt$npc], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
## ---- one pass over all cells: perturbation means and PC coordinates -----------------------------------------------------------
labs <- sort(unique(gt[keep_cells]), method = "radix"); li <- match(gt, labs)
SUM <- matrix(0, G, length(labs)); CNT <- integer(length(labs)); Ycoord <- matrix(0, opt$npc, Ncell)
for (st in seq(1L, length(keep_cells), by = opt$chunk)) { idx <- keep_cells[st:min(st + opt$chunk - 1L, length(keep_cells))]
  z <- znorm(readX(idx), idx); Ycoord[, idx] <- crossprod(V, z - cen)
  ind <- Matrix::sparseMatrix(i = seq_along(idx), j = li[idx], x = 1, dims = c(length(idx), length(labs)))
  SUM <- SUM + as.matrix(z %*% ind); CNT <- CNT + as.integer(Matrix::colSums(ind)) }
Xbulk <- SUM / rep(CNT, each = G); dimnames(Xbulk) <- list(genes, labs)
say(sprintf("%d labels (%d non-targeting), cells per label: median %g; %.0f s; VmHWM %.1f GB", length(labs), sum(grepl("non-targeting", labs)), stats::median(CNT), as.numeric(difftime(Sys.time(), t0, units = "secs")), vmhwm()))
## ---- energy-distance test ------------------------------------------------------------------------------------------------------
pv <- rep(NA_real_, length(labs)); En <- rep(NA_real_, length(labs)); tst <- which(!grepl("non-targeting", labs) & CNT >= 2L)
if (opt$limit_perts > 0L) tst <- head(tst, opt$limit_perts)
one <- function(k) { idx <- which(li == k & seq_len(Ncell) %in% keep_cells); set.seed(opt$seed + k); cs <- sample(ntc, min(opt$nctrl, length(ntc)))
  r <- energy_perm(Ycoord[, c(idx, cs), drop = FALSE], length(idx), opt$nperm, opt$seed + k); c(r$E, r$p) }
ord <- tst[order(-CNT[tst])]; res <- mclapply(ord, one, mc.cores = opt$threads, mc.preschedule = TRUE)
for (j in seq_along(ord)) { pv[ord[j]] <- res[[j]][2]; En[ord[j]] <- res[[j]][1] }
say(sprintf("energy test: %d perturbations tested, %d with p < 1e-3; %.0f s; VmHWM %.1f GB", length(tst), sum(pv < 1e-3, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "secs")), vmhwm()))
## ---- write h5ad-like output ----------------------------------------------------------------------------------------------------
if (file.exists(opt$out)) file.remove(opt$out)
h5createFile(opt$out); h5write(Xbulk, opt$out, "X"); h5createGroup(opt$out, "obs"); h5createGroup(opt$out, "var")
h5write(labs, opt$out, "obs/gene_transcript"); h5write(pv, opt$out, "obs/energy_test_p_value"); h5write(En, opt$out, "obs/energy_statistic"); h5write(as.integer(CNT), opt$out, "obs/num_cells")
h5write(genes, opt$out, "var/gene_name"); H5close()
say(sprintf("wrote %s; wall %.0f s; VmHWM %.1f GB", opt$out, as.numeric(difftime(Sys.time(), t0, units = "secs")), vmhwm()))
