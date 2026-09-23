#!/usr/bin/env Rscript
# M3D E. coli (v4 build 6) x RegulonDB evaluation of CLR-as-attention.
#
# Usage (from the repo root):
#   Rscript analysis/run_m3d_regulondb.R \
#       --m3d data/E_coli_v4_Build_6 --regulondb data/NetworkRegulatorGene.tsv \
#       [--set chips|avg] [--B 100] [--threads N] [--tmax 20] [--alpha 0.5]
#       [--quick 600] [--mi-null] [--out results/<set>]
#
# Stages (all outputs under --out):
#   1. Load compendium (genes = b-numbers) and RegulonDB; build the
#      Faith-2007-style evaluation universe (every known-TF x gene pair).
#   2. Edge-level PR for four MI/CLR configurations:
#        parity2007  transform none, 10 bins, Euclidean   (historical CLR)
#        none_fd     transform none, FD bins,  Euclidean   (pre-review default)
#        rank_fd     transform rank, FD bins,  Euclidean   (new default)
#        rank_fd_st  transform rank, FD bins,  Stouffer
#      scored by raw MI and by CLR, against all and strong/confirmed-only
#      RegulonDB interactions.
#   3. Primary config (rank_fd): permutation edge selection (BH on the CLR
#      null, B replicates), precision/recall of the selected set; optional
#      MI-null comparison (--mi-null, costs another B MI builds).
#   4. Design-A diffusion depth sweep with the fixed operator: per depth t,
#      per-TF regulon average precision when genes are ranked by
#        (a) multi-hop attention mass (P^t)[TF, ], and
#        (b) correlation of diffused profiles cor(E_t[TF, ], E_t[g, ]),
#      plus edge-level AUPR of |cor(E_t)| on the universe.
#   5. Artifacts for reproducibility (see docs/diffusion_depth_design_review.md
#      section 9): MI/CLR/edges/operator of the primary config, parameters,
#      seed, sessionInfo, timings.

## ---- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", regulondb = NULL, set = "chips",
            B = 100L, threads = NULL, tmax = 20L, alpha = 0.5, quick = 0L,
            mi_null = FALSE, out = NULL, seed = 20260922L)
i <- 1L
while (i <= length(args)) {
  a <- args[i]
  key <- gsub("-", "_", sub("^--", "", a))
  if (key == "mi_null") { opt$mi_null <- TRUE; i <- i + 1L; next }
  if (!key %in% names(opt)) stop("unknown argument: ", a)
  val <- args[i + 1L]
  opt[[key]] <- if (key %in% c("B", "threads", "tmax", "quick", "seed"))
    as.integer(val) else if (key == "alpha") as.numeric(val) else val
  i <- i + 2L
}
if (is.null(opt$regulondb)) stop("--regulondb <path> is required")
if (is.null(opt$out)) opt$out <- file.path("results", opt$set)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "run.log")
say <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(msg, "\n"); cat(msg, "\n", file = logf, append = TRUE)
}

## ---- package ---------------------------------------------------------------
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
repo_root <- normalizePath(file.path(script_dir, ".."))
if (requireNamespace("clr", quietly = TRUE)) {
  library(clr)
} else {
  pkgload::load_all(repo_root, quiet = TRUE)
}
suppressPackageStartupMessages(library(Matrix))
source(file.path(script_dir, "regulondb.R"))
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=",
                       collapse = " "))
timings <- list()
tic <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  timings[[label]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  say(sprintf("%s: %.1f s", label, timings[[label]]))
  invisible(val)
}

## ---- 1. data -----------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
f <- if (opt$set == "chips")
  file.path(opt$m3d, paste0(base, "_chips907probes4297.tab")) else
  file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab"))
raw <- tic("read_m3d", utils::read.delim(f, check.names = FALSE,
                                         stringsAsFactors = FALSE))
probe <- raw[[1]]
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
bnum <- vapply(strsplit(probe, "_"), function(p) p[length(p) - 1L], "")
desc <- utils::read.delim(file.path(opt$m3d, paste0(base, ".probe_set_descriptions")),
                          stringsAsFactors = FALSE)
if (!all(c("probe_set_name", "gene_symbol") %in% names(desc)))
  stop("probe_set_descriptions lacks probe_set_name/gene_symbol columns: ",
       paste(names(desc), collapse = ", "))
sym <- desc$gene_symbol[match(probe, desc$probe_set_name)]
from_probe <- vapply(strsplit(probe, "_"), function(p) p[1], "")
miss <- is.na(sym) | !nzchar(sym)
sym[miss] <- from_probe[miss]
say(sprintf("symbols: %d from descriptions, %d from probe names", sum(!miss), sum(miss)))
rownames(X) <- bnum
stopifnot(!anyDuplicated(bnum))
rng <- apply(X, 1, function(x) diff(range(x)))
if (any(rng <= 0)) {
  say("dropping ", sum(rng <= 0), " constant genes")
  X <- X[rng > 0, , drop = FALSE]; sym <- sym[rng > 0]
}
if (anyNA(X)) stop("NA values in the compendium")
if (opt$quick > 0L) {
  v <- apply(X, 1, stats::var)
  keep <- sort(order(-v)[seq_len(min(opt$quick, nrow(X)))])
  X <- X[keep, ]; sym <- sym[keep]
  say("QUICK MODE: top-", nrow(X), " variance genes only")
}
say("compendium: ", nrow(X), " genes x ", ncol(X), " arrays (", opt$set, ")")

net <- read_regulondb(opt$regulondb)
u_all <- edge_universe(net, sym)
u_str <- edge_universe(net, sym, conf_keep = c("S", "C"))
say(sprintf("universe (all evidence): %d TFs, %d pairs, %d positives (%d mapped rows)",
            length(u_all$tfs), length(u_all$label), u_all$n_pos, u_all$n_known_rows))
say(sprintf("universe (strong/confirmed): %d TFs, %d pairs, %d positives",
            length(u_str$tfs), length(u_str$label), u_str$n_pos))

## ---- 2. edge-level PR for four configurations -----------------------------------
configs <- list(
  parity2007 = list(transform = "none", bins = 10, combine = "euclidean"),
  none_fd    = list(transform = "none", bins = "fd", combine = "euclidean"),
  rank_fd    = list(transform = "rank", bins = "fd", combine = "euclidean"),
  rank_fd_st = list(transform = "rank", bins = "fd", combine = "stouffer")
)
edge_rows <- list()
fits <- list()
for (nm in names(configs)) {
  cf <- configs[[nm]]
  fit <- ClrAttention$new(X)
  tic(paste0("mi_", nm), fit$estimate_mi(bins = cf$bins, transform = cf$transform,
                                         threads = opt$threads))
  fit$calibrate(method = "normal", combine = cf$combine)
  b <- fit$params$mi$bins_used
  say(sprintf("%s: bins median %g (range %d-%d)", nm, stats::median(b), min(b), max(b)))
  for (ref in c("all", "strong")) {
    u <- if (ref == "all") u_all else u_str
    for (st in c("mi", "clr")) {
      V <- if (st == "mi") fit$mi else fit$clr_scores
      s <- pr_summary(V[cbind(u$i, u$j)], u$label)
      edge_rows[[length(edge_rows) + 1L]] <- data.frame(
        config = nm, reference = ref, score = st, t(s), check.names = FALSE)
    }
  }
  if (nm == "rank_fd") fits[[nm]] <- fit
  rm(fit); invisible(gc())
}
edge_tab <- do.call(rbind, edge_rows)
utils::write.csv(edge_tab, file.path(opt$out, "edge_pr_by_config.csv"), row.names = FALSE)
print(edge_tab, digits = 3, row.names = FALSE)

## ---- 3. permutation edge selection (primary config) ----------------------------
fit <- fits$rank_fd
set.seed(opt$seed)
tic("select_threshold_clr", fit$select_threshold(B = opt$B, method = "fdr", q = 0.05,
                                                 threads = opt$threads,
                                                 statistic = "clr"))
E_sel <- fit$edges
sel_row <- function(nm, f, E) data.frame(
  null = nm, tau = f$threshold, n_edges = sum(E) / 2,
  mean_degree = mean(rowSums(E)), t(selected_pr(E, u_all)),
  strong_recall = selected_pr(E, u_str)[["recall"]])
sel_rows <- list(sel_row("clr", fit, E_sel))
if (opt$mi_null) {
  fit_m <- fit$clone()
  set.seed(opt$seed)
  tic("select_threshold_mi", fit_m$select_threshold(B = opt$B, method = "fdr",
                                                    threads = opt$threads,
                                                    statistic = "mi"))
  sel_rows[[2]] <- sel_row("mi", fit_m, fit_m$edges)
  rm(fit_m)
}
sel_tab <- do.call(rbind, sel_rows)
utils::write.csv(sel_tab, file.path(opt$out, "selected_edges.csv"), row.names = FALSE)
print(sel_tab, digits = 3, row.names = FALSE)

## ---- 4. diffusion depth sweep (Design A: fixed operator) ------------------------
fit$build_operator(alpha = opt$alpha)
A <- Matrix::Matrix(fit$operator, sparse = TRUE)
G <- nrow(A)
P <- (1 - opt$alpha) * Matrix::Diagonal(G) + opt$alpha * A
say(sprintf("operator: %d nonzeros off-diagonal, %d isolated genes (self-loop)",
            sum(fit$operator > 0 & row(fit$operator) != col(fit$operator)),
            sum(diag(fit$operator) == 1)))
Z <- X - rowMeans(X); Z <- Z / sqrt(rowSums(Z^2) / (ncol(Z) - 1))
tfs <- u_all$tfs
Rt <- Matrix::sparseMatrix(i = seq_along(tfs), j = tfs, x = 1, dims = c(length(tfs), G))
Rt <- as.matrix(Rt)                    # (P^t)[tfs, ], starts at identity rows
Et <- Z
depth <- list()
for (t in 0:opt$tmax) {
  if (t > 0) {
    Rt <- as.matrix(Rt %*% P)
    Et <- as.matrix(P %*% Et)
  }
  Ctf <- stats::cor(t(Et[tfs, , drop = FALSE]), t(Et))       # |tfs| x G
  ap_att <- regulon_ap(Rt, u_all)
  ap_cor <- regulon_ap(abs(Ctf), u_all)
  # edge-level |cor| AUPR on the universe
  tf_row <- match(u_all$i, tfs); other <- u_all$j
  swap <- is.na(tf_row); tf_row[swap] <- match(u_all$j[swap], tfs); other[swap] <- u_all$i[swap]
  ci <- abs(Ctf[cbind(tf_row, other)])
  e <- pr_summary(ci, u_all$label)
  depth[[t + 1L]] <- data.frame(
    t = t,
    regulon_ap_attention_median = if (t == 0) NA_real_ else stats::median(ap_att, na.rm = TRUE),
    regulon_ap_cor_median = stats::median(ap_cor, na.rm = TRUE),
    edge_aupr_abscor = e[["aupr"]],
    n_tfs_scored = sum(!is.na(ap_cor)),
    median_abs_cor_all = stats::median(abs(Ctf)),
    row_sd_min = min(apply(Et, 1, stats::sd)))
  if (t %% 5 == 0 || t == opt$tmax)
    say(sprintf("t=%2d  regulonAP(att)=%.4f  regulonAP(|cor|)=%.4f  edgeAUPR(|cor|)=%.4f  median|cor|=%.3f",
                t, depth[[t + 1L]]$regulon_ap_attention_median,
                depth[[t + 1L]]$regulon_ap_cor_median, e[["aupr"]],
                depth[[t + 1L]]$median_abs_cor_all))
}
depth_tab <- do.call(rbind, depth)
utils::write.csv(depth_tab, file.path(opt$out, "depth_sweep.csv"), row.names = FALSE)

png(file.path(opt$out, "depth_sweep.png"), width = 900, height = 560, res = 110)
par(mar = c(5, 4.5, 3, 1))
yl <- range(unlist(depth_tab[, c("regulon_ap_attention_median",
                                 "regulon_ap_cor_median", "edge_aupr_abscor")]),
            na.rm = TRUE)
plot(depth_tab$t, depth_tab$regulon_ap_cor_median, type = "b", pch = 16,
     col = "#1b6ca8", ylim = yl, xlab = "diffusion depth t",
     ylab = "average precision", main = paste0("M3D (", opt$set,
                                              ") x RegulonDB: depth sweep"))
lines(depth_tab$t, depth_tab$regulon_ap_attention_median, type = "b", pch = 17,
      col = "#e8833a")
lines(depth_tab$t, depth_tab$edge_aupr_abscor, type = "b", pch = 15, col = "#2ca02c")
legend("topright", bty = "n", pch = c(16, 17, 15),
       col = c("#1b6ca8", "#e8833a", "#2ca02c"),
       legend = c("per-TF regulon AP, |cor| of diffused profiles (median)",
                  "per-TF regulon AP, attention mass (P^t) (median)",
                  "edge AUPR, |cor| of diffused profiles"))
dev.off()

## ---- 5. artifacts ------------------------------------------------------------
saveRDS(list(mi = fit$mi, clr = fit$clr_scores, edges = E_sel,
             operator = A, params = fit$params, genes = rownames(X),
             symbols = sym, seed = opt$seed, options = opt),
        file.path(opt$out, "primary_fit.rds"), compress = "gzip")
utils::write.csv(data.frame(stage = names(timings), seconds = unlist(timings)),
                 file.path(opt$out, "timings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(opt$out, "sessionInfo.txt"))
say("done; outputs in ", normalizePath(opt$out))
