#!/usr/bin/env Rscript
# M3D E. coli (v4 build 6) x RegulonDB evaluation of CLR-as-attention.
#
# Usage (from the repo root):
#   Rscript analysis/run_m3d_regulondb.R \
#       [--m3d data/E_coli_v4_Build_6] [--rdb data/RegulonDBExtract] \
#       [--set chips|avg] [--B 100] [--threads N] [--alpha 0.5] \
#       [--depths 0,1,2,3,5,8,12,20] [--min-size 5] [--max-size 500] \
#       [--quick 600] [--mi-null] [--out results/<set>]
#       [--primary none_hg (default: raw values, Hacine-Gharbi joint-histogram
#        bin rule) | none_scott2d | none_median_scott | parity2007 | rank_fd | ...]
#       [--bin-sweep 6,8,12,16]    (raw-value fixed-bin sensitivity; "" to skip)
#       [--reuse results/<set>/primary_fit.rds]  skip the configuration sweep
#        and the permutation null: recompute the primary MI/CLR (seconds) and
#        rebuild the operator from the saved CLR threshold (deterministic).
#
# PRIMARY BENCHMARK (regulator-agnostic; see analysis/regulons.R): regulons of
# every regulator type in RegulonDB -- TFs, sRNAs, small molecules (ppGpp),
# other proteins, and sigma factors (sigmulons) -- expanded from promoter/TU
# targets to genes, size window [min-size, max-size]. A gene pair is positive
# if the two genes share a regulon; same-operon pairs are excluded. Scored by
#   * co-membership AUPR/AUROC over all annotated gene pairs, overall and per
#     regulator class, for strong/confirmed ("SC") and all ("all") evidence;
#   * per-regulon coherence AUROC (within-regulon pairs vs member-to-non-
#     co-regulated pairs), summarized as the median over regulons.
# SECONDARY (continuity with Faith et al. 2007): TF-node edge PR, where the TF's
# own mRNA stands in for its activity.
#
# Stages:
#   1. data + regulons
#   2. four MI/CLR configurations (parity2007, none_fd, rank_fd, rank_fd_st),
#      each scored by raw MI and by CLR; plus plain |Pearson| as a baseline
#   3. permutation edge selection (BH, CLR null) on rank_fd; optional MI null
#   4. Design-A depth sweep with the fixed operator: similarity of diffused
#      profiles |cor(E_t)| and symmetrized attention mass (P^t + P^t') / 2;
#      diffusion with raw, signed and conditional-expectation values
#   5. artifacts (MI/CLR/edges/operator, parameters, seed, sessionInfo, timings)

## ---- arguments ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            set = "chips", B = 100L, threads = NULL, alpha = 0.5,
            depths = "0,1,2,3,5,8,12,20", min_size = 5L, max_size = 500L,
            quick = 0L, mi_null = FALSE, out = NULL, seed = 20260922L,
            primary = "none_hg", bin_sweep = "6,8,12,16", reuse = "")
int_opts <- c("B", "threads", "quick", "seed", "min_size", "max_size")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (key == "mi_null") { opt$mi_null <- TRUE; i <- i + 1L; next }
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  val <- args[i + 1L]
  opt[[key]] <- if (key %in% int_opts) as.integer(val) else
    if (key == "alpha") as.numeric(val) else val
  i <- i + 2L
}
depths <- sort(unique(as.integer(strsplit(opt$depths, ",")[[1]])))
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
if (requireNamespace("clr", quietly = TRUE)) library(clr) else
  pkgload::load_all(repo_root, quiet = TRUE)
suppressPackageStartupMessages(library(Matrix))
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))
say("options: ", paste(names(opt), unlist(lapply(opt, format)), sep = "=",
                       collapse = " "))
say("clr OpenMP: ", paste(names(clr_openmp_info()), clr_openmp_info(),
                          sep = "=", collapse = " "))
timings <- list()
tic <- function(label, expr) {
  t0 <- Sys.time(); val <- force(expr)
  timings[[label]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  say(sprintf("%s: %.1f s", label, timings[[label]]))
  invisible(val)
}
wcsv <- function(x, name) utils::write.csv(x, file.path(opt$out, name),
                                           row.names = FALSE)

## ---- 1. data and regulons -------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
f <- if (opt$set == "chips")
  file.path(opt$m3d, paste0(base, "_chips907probes4297.tab")) else
  file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab"))
raw <- tic("read_m3d", utils::read.delim(f, check.names = FALSE,
                                         stringsAsFactors = FALSE))
probe <- raw[[1]]
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
bnum <- vapply(strsplit(probe, "_"), function(p) p[length(p) - 1L], "")
stopifnot(all(grepl("^b[0-9]+$", bnum)), !anyDuplicated(bnum))
rownames(X) <- bnum
rng <- apply(X, 1, function(x) diff(range(x)))
if (any(rng <= 0)) {
  say("dropping ", sum(rng <= 0), " constant genes")
  X <- X[rng > 0, , drop = FALSE]
}
if (anyNA(X)) stop("NA values in the compendium")
if (opt$quick > 0L) {
  v <- apply(X, 1, stats::var)
  X <- X[sort(order(-v)[seq_len(min(opt$quick, nrow(X)))]), ]
  say("QUICK MODE: top-", nrow(X), " variance genes only")
}
genes <- rownames(X); G <- length(genes)
say("compendium: ", G, " genes x ", ncol(X), " arrays (", opt$set, ")")

to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)
operon_of <- build_operon_map(opt$rdb, to_bn)
say(sprintf("RegulonDB: %d regulator-target rows, %d target names unmapped to b-numbers; %d regulated genes, %d on the array",
            nrow(rtab), sum(is.na(rtab$bnumber)), length(unique(na.omit(rtab$bnumber))),
            length(intersect(unique(rtab$bnumber), genes))))
bench <- list()
for (ev in c("SC", "all")) {
  ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
  reg <- make_regulons(rtab, genes, conf_keep = ck, min_size = opt$min_size,
                       max_size = opt$max_size)
  pairs <- comembership_pairs(reg, operon_of, genes)
  co <- comember_matrix(reg, G)
  say(sprintf("regulons [%s]: %d (%s); dropped by size: %s; annotated genes %d, pairs %d, positive %d (%.1f%%), same-operon removed %d",
              ev, length(reg$members),
              paste(names(table(reg$class)), table(reg$class), sep = "=", collapse = " "),
              paste(sprintf("%s(%d)", reg$dropped$regulon[reg$dropped$size > opt$max_size],
                            reg$dropped$size[reg$dropped$size > opt$max_size]), collapse = " "),
              length(pairs$genes_in), length(pairs$label), sum(pairs$label),
              100 * mean(pairs$label), pairs$n_same_operon_removed))
  bench[[ev]] <- list(reg = reg, pairs = pairs, co = co)
}
wcsv(data.frame(regulon = names(bench$all$reg$members), class = bench$all$reg$class,
                size = bench$all$reg$size), "regulons_all.csv")

# Score a similarity matrix on every benchmark; returns long data.frames.
score_all <- function(S, label) {
  cm <- lapply(names(bench), function(ev) {
    r <- comembership_eval(S, bench[[ev]]$pairs)
    cbind(method = label, evidence = ev, r)
  })
  coh <- lapply(names(bench), function(ev) {
    b <- bench[[ev]]
    a <- regulon_coherence(S, b$reg, b$co, operon_of, genes)
    rbind(data.frame(method = label, evidence = ev, class = "all",
                     median_auroc = stats::median(a, na.rm = TRUE),
                     frac_gt_0.6 = mean(a > 0.6, na.rm = TRUE),
                     n = sum(!is.na(a))),
          do.call(rbind, lapply(sort(unique(b$reg$class)), function(c) {
            ac <- a[b$reg$class == c]
            data.frame(method = label, evidence = ev, class = c,
                       median_auroc = stats::median(ac, na.rm = TRUE),
                       frac_gt_0.6 = mean(ac > 0.6, na.rm = TRUE),
                       n = sum(!is.na(ac)))
          })))
  })
  list(comembership = do.call(rbind, cm), coherence = do.call(rbind, coh))
}
cm_rows <- list(); coh_rows <- list()
add_scores <- function(sc) {
  cm_rows[[length(cm_rows) + 1L]] <<- sc$comembership
  coh_rows[[length(coh_rows) + 1L]] <<- sc$coherence
}
headline <- function(label) {
  cm <- do.call(rbind, cm_rows); coh <- do.call(rbind, coh_rows)
  a <- cm[cm$method == label & cm$stratum == "all", ]
  h <- coh[coh$method == label & coh$class == "all", ]
  say(sprintf("  %-22s co-membership AUPR SC %.4f (base %.4f) all %.4f (base %.4f) | coherence median AUROC SC %.3f all %.3f",
              label, a$aupr[a$evidence == "SC"], a$base_rate[a$evidence == "SC"],
              a$aupr[a$evidence == "all"], a$base_rate[a$evidence == "all"],
              h$median_auroc[h$evidence == "SC"], h$median_auroc[h$evidence == "all"]))
}

# Secondary: TF-node network (TF class only), TF gene resolved by name.
tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
tf_genes <- .classic_tf_to_genes(tfr$regulator)   # CRP -> cRP, IHF -> ihfA;ihfB ...
tf_parts <- strsplit(tf_genes, ";", fixed = TRUE)
tfnet <- data.frame(tf_gene = to_bn(unlist(tf_parts)),
                    target = rep(tfr$bnumber, lengths(tf_parts)),
                    confidence = rep(tfr$confidence, lengths(tf_parts)),
                    stringsAsFactors = FALSE)
say(sprintf("TF-node network: %d/%d TF names resolved to a b-number",
            length(unique(tfr$regulator[!is.na(to_bn(vapply(tf_parts, `[`, "", 1L)))])),
            length(unique(tfr$regulator))))
tfnet <- tfnet[!is.na(tfnet$tf_gene), ]
u_all <- edge_universe(tfnet, genes)
u_str <- edge_universe(tfnet, genes, conf_keep = c("S", "C"))
tfnode <- function(S, label) {
  do.call(rbind, lapply(c("SC", "all"), function(ev) {
    u <- if (ev == "SC") u_str else u_all
    data.frame(method = label, evidence = ev,
               t(pr_summary(S[cbind(u$i, u$j)], u$label)), check.names = FALSE)
  }))
}
tf_rows <- list()

## ---- 2. configurations ---------------------------------------------------------
Z <- X - rowMeans(X); Z <- Z / sqrt(rowSums(Z^2) / (ncol(Z) - 1))
C0 <- abs(tcrossprod(Z) / (ncol(Z) - 1))
add_scores(score_all(C0, "abs_pearson")); headline("abs_pearson")
tf_rows[[length(tf_rows) + 1L]] <- tfnode(C0, "abs_pearson")
rm(C0); invisible(gc())

configs <- list(
  parity2007        = list(transform = "none", bins = 10, combine = "euclidean"),
  none_fd           = list(transform = "none", bins = "fd", combine = "euclidean"),
  # historical practice: per-gene optimal count, median used for all genes
  none_median_fd    = list(transform = "none", bins = "median_fd", combine = "euclidean"),
  none_median_scott = list(transform = "none", bins = "median_scott", combine = "euclidean"),
  # 2-D-aware counts: the budget is the joint histogram, not the marginal
  none_scott2d      = list(transform = "none", bins = "median_scott2d", combine = "euclidean"),
  none_hg           = list(transform = "none", bins = "hg", combine = "euclidean"),
  rank_fd           = list(transform = "rank", bins = "fd", combine = "euclidean"),
  rank_fd_st        = list(transform = "rank", bins = "fd", combine = "stouffer")
)
sweep <- if (nzchar(opt$bin_sweep)) as.integer(strsplit(opt$bin_sweep, ",")[[1]]) else integer()
for (nb in sweep)
  configs[[sprintf("none_b%02d", nb)]] <- list(transform = "none", bins = nb,
                                                combine = "euclidean")
if (!opt$primary %in% names(configs))
  stop("--primary must be one of: ", paste(names(configs), collapse = ", "))

# Parametric (Gaussian) MI: I = -1/2 log2(1 - r^2) on the raw values, then CLR.
{
  r2 <- pmin((tcrossprod(Z) / (ncol(Z) - 1))^2, 1 - 1e-12)
  Mg <- -0.5 * log2(1 - r2); diag(Mg) <- 0
  add_scores(score_all(Mg, "gaussian:mi")); headline("gaussian:mi")
  Sg <- clr_calibrate(Mg, method = "normal", combine = "euclidean")
  add_scores(score_all(Sg, "gaussian:clr")); headline("gaussian:clr")
  tf_rows[[length(tf_rows) + 1L]] <- tfnode(Sg, "gaussian:clr")
  rm(r2, Mg, Sg); invisible(gc())
}
reuse <- NULL
if (nzchar(opt$reuse)) {
  reuse <- readRDS(opt$reuse)
  stopifnot(identical(reuse$genes, genes),
            identical(reuse$params$threshold$statistic, "clr"))
  opt$primary <- reuse$options$primary
  configs <- configs[opt$primary]
  say("REUSE: ", opt$reuse, " (primary ", opt$primary, ", tau ",
      format(reuse$params$threshold$tau, digits = 6), ")")
}
fit <- NULL
for (nm in names(configs)) {
  cf <- configs[[nm]]
  f_ <- ClrAttention$new(X)
  tic(paste0("mi_", nm), f_$estimate_mi(bins = cf$bins, transform = cf$transform,
                                        threads = opt$threads))
  f_$calibrate(method = "normal", combine = cf$combine)
  b <- f_$params$mi$bins_used
  say(sprintf("%s: bins median %g (range %d-%d)", nm, stats::median(b), min(b), max(b)))
  M <- f_$mi; diag(M) <- 0
  if (nm != "rank_fd_st") {   # MI is identical to rank_fd's (only CLR differs)
    add_scores(score_all(M, paste0(nm, ":mi"))); headline(paste0(nm, ":mi"))
    tf_rows[[length(tf_rows) + 1L]] <- tfnode(M, paste0(nm, ":mi"))
  }
  add_scores(score_all(f_$clr_scores, paste0(nm, ":clr"))); headline(paste0(nm, ":clr"))
  tf_rows[[length(tf_rows) + 1L]] <- tfnode(f_$clr_scores, paste0(nm, ":clr"))
  if (nm == opt$primary) fit <- f_
  rm(f_, M); invisible(gc())
}
wcsv(do.call(rbind, cm_rows), "comembership_by_config.csv")
wcsv(do.call(rbind, coh_rows), "coherence_by_config.csv")
wcsv(do.call(rbind, tf_rows), "tfnode_edge_pr_by_config.csv")
n_cfg_cm <- length(cm_rows)

## ---- 3. permutation edge selection (primary config) ----------------------------
say("primary configuration for selection and diffusion: ", opt$primary)
sel_eval <- function(E, label, tau) {
  do.call(rbind, lapply(names(bench), function(ev) {
    p <- bench[[ev]]$pairs
    s <- E[cbind(p$i, p$j)]
    data.frame(null = label, evidence = ev, tau = tau, n_edges = sum(E) / 2,
               mean_degree = mean(rowSums(E)),
               selected_annotated_pairs = sum(s),
               comember_precision = if (sum(s)) mean(p$label[s]) else NA_real_,
               comember_base_rate = mean(p$label),
               comember_recall = sum(s & p$label == 1) / sum(p$label),
               tfnode_precision = selected_pr(E, if (ev == "SC") u_str else u_all)[["precision"]])
  }))
}
if (!is.null(reuse)) {
  E_sel <- reuse$edges
} else {
set.seed(opt$seed)
tic("select_threshold_clr", fit$select_threshold(B = opt$B, method = "fdr", q = 0.05,
                                                 threads = opt$threads,
                                                 statistic = "clr"))
E_sel <- fit$edges
sel_rows <- list(sel_eval(E_sel, "clr", fit$threshold))
if (opt$mi_null) {
  fit_m <- fit$clone()
  set.seed(opt$seed)
  tic("select_threshold_mi", fit_m$select_threshold(B = opt$B, method = "fdr",
                                                    threads = opt$threads,
                                                    statistic = "mi"))
  sel_rows[[2]] <- sel_eval(fit_m$edges, "mi", fit_m$threshold)
  rm(fit_m)
}
sel_tab <- do.call(rbind, sel_rows)
wcsv(sel_tab, "selected_edges.csv")
print(sel_tab, digits = 3, row.names = FALSE)
}

## ---- 4. diffusion depth sweep (Design A: fixed operator) ------------------------
if (is.null(reuse)) fit$build_operator(alpha = opt$alpha) else {
  fit$build_operator(tau = reuse$params$threshold$tau, alpha = opt$alpha)
  E_now <- fit$clr_scores >= reuse$params$threshold$tau
  diag(E_now) <- FALSE
  if (!identical(unname(E_now), unname(reuse$edges)))
    stop("reused threshold does not reproduce the saved edge set")
}
A <- Matrix::Matrix(fit$operator, sparse = TRUE)
P <- (1 - opt$alpha) * Matrix::Diagonal(G) + opt$alpha * A
say(sprintf("operator: %d off-diagonal nonzeros, %d isolated genes (self-loop), alpha %.2f",
            Matrix::nnzero(A) - sum(Matrix::diag(A) != 0), sum(Matrix::diag(A) == 1),
            opt$alpha))
Pt <- diag(G)                        # dense P^t
Et <- Z                              # E^(0): row-standardized data
depth_rows <- list()
for (t in 0:max(depths)) {
  if (t > 0) {
    Pt <- as.matrix(Pt %*% P)
    Et <- as.matrix(P %*% Et)
  }
  if (!t %in% depths) next
  Ez <- Et - rowMeans(Et)
  sdv <- sqrt(rowSums(Ez^2) / (ncol(Ez) - 1)); sdv[sdv <= 0] <- 1
  Ez <- Ez / sdv
  Ct <- abs(tcrossprod(Ez) / (ncol(Ez) - 1))
  sc <- score_all(Ct, sprintf("depth%02d:abscor", t)); add_scores(sc)
  headline(sprintf("depth%02d:abscor", t))
  if (t > 0) {
    Sa <- (Pt + t(Pt)) / 2; diag(Sa) <- 0
    add_scores(score_all(Sa, sprintf("depth%02d:attention", t)))
    headline(sprintf("depth%02d:attention", t))
  }
  depth_rows[[length(depth_rows) + 1L]] <- data.frame(
    t = t, median_abs_cor_offdiag = stats::median(Ct[upper.tri(Ct)]),
    effective_rank = { ev <- svd(Ez, nu = 0, nv = 0)$d^2; sum(ev)^2 / sum(ev^2) })
  rm(Ct); invisible(gc())
}
# Value transforms (the continuous W_V): same operator and readout, but the
# message from gene j to gene i is sign(cor) * E_j ("signed") or the B-spline
# conditional expectation E[z_i | z_j] applied to E_j ("conditional").
for (vt in c("signed", "conditional")) {
  tic(paste0("diffuse_", vt), fit$diffuse(steps = max(depths), values = vt))
  tr <- fit$trajectory
  for (t in depths[depths > 0]) {
    Ez <- tr[[t + 1L]] - rowMeans(tr[[t + 1L]])
    sdv <- sqrt(rowSums(Ez^2) / (ncol(Ez) - 1)); sdv[sdv <= 0] <- 1
    Ez <- Ez / sdv
    Ct <- abs(tcrossprod(Ez) / (ncol(Ez) - 1))
    lab <- sprintf("depth%02d:abscor_%s", t, vt)
    add_scores(score_all(Ct, lab)); headline(lab)
    rm(Ct); invisible(gc())
  }
  fit$reset_diffusion(); rm(tr); invisible(gc())
}

cm_all <- do.call(rbind, cm_rows); coh_all <- do.call(rbind, coh_rows)
dep_cm <- cm_all[grepl("^depth", cm_all$method), ]
dep_coh <- coh_all[grepl("^depth", coh_all$method), ]
wcsv(dep_cm, "depth_comembership.csv")
wcsv(dep_coh, "depth_coherence.csv")
wcsv(do.call(rbind, depth_rows), "depth_collapse.csv")

png(file.path(opt$out, "depth_sweep.png"), width = 1000, height = 560, res = 110)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
for (what in c("aupr", "coh")) {
  for (kind in c("abscor", "attention")) {
    if (what == "aupr") {
      d <- dep_cm[dep_cm$stratum == "all" & dep_cm$evidence == "SC" &
                    grepl(paste0(":", kind, "$"), dep_cm$method), ]
      y <- d$aupr
    } else {
      d <- dep_coh[dep_coh$class == "all" & dep_coh$evidence == "SC" &
                     grepl(paste0(":", kind, "$"), dep_coh$method), ]
      y <- d$median_auroc
    }
    tt <- as.integer(sub("^depth([0-9]+):.*$", "\\1", d$method))
    if (kind == "abscor") {
      plot(tt, y, type = "b", pch = 16, col = "#1b6ca8",
           ylim = range(c(y, if (what == "aupr")
             dep_cm$aupr[dep_cm$stratum == "all" & dep_cm$evidence == "SC"] else
               dep_coh$median_auroc[dep_coh$class == "all" & dep_coh$evidence == "SC"]),
             na.rm = TRUE),
           xlab = "diffusion depth t",
           ylab = if (what == "aupr") "co-membership AUPR (SC)" else
             "median regulon coherence AUROC (SC)",
           main = if (what == "aupr") "Regulon co-membership" else "Regulon coherence")
    } else lines(tt, y, type = "b", pch = 17, col = "#e8833a")
  }
  legend("bottomright", bty = "n", pch = c(16, 17), col = c("#1b6ca8", "#e8833a"),
         legend = c("|cor| of diffused profiles", "attention mass (P^t sym.)"))
}
dev.off()

## ---- 5. artifacts ------------------------------------------------------------
saveRDS(list(mi = fit$mi, clr = fit$clr_scores, edges = E_sel,
             operator = A, params = fit$params, genes = genes,
             seed = opt$seed, options = opt),
        file.path(opt$out, "primary_fit.rds"), compress = "gzip")
wcsv(data.frame(stage = names(timings), seconds = unlist(timings)), "timings.csv")
writeLines(capture.output(sessionInfo()), file.path(opt$out, "sessionInfo.txt"))
say("done; outputs in ", normalizePath(opt$out))
