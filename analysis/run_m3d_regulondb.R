#!/usr/bin/env Rscript
# M3D E. coli (v4 build 6) x RegulonDB evaluation of CLR-as-attention.
#
# Usage (from the repo root):
#   Rscript analysis/run_m3d_regulondb.R \
#       [--m3d data/E_coli_v4_Build_6] [--rdb data/RegulonDBExtract] \
#       [--set avg|chips]  (default avg: 466 replicate-averaged experiments;
#        chips = 907 arrays with technical replicates counted as samples) [--B 100] [--threads N] [--alpha 0.5] \
#       [--depths 0,1,2,3,5,8,12,20] [--min-size 5] [--max-size 500] \
#       [--quick 600] [--mi-null] [--out results/<set>]
#       [--primary none_hg (default: raw values, Hacine-Gharbi joint-histogram
#        bin rule, Stouffer CLR) | none_hg_eu | none_scott2d | parity2007 | ...]
#       [--bin-sweep 6,8,12,16]    (raw-value fixed-bin sensitivity; "" to skip)
#       [--operators fdr05,fdr20,fdr05_top5,fdr05_top10,soft10,pearson_top5,pearson_top10]
#        attention operators for the depth sweep (pre-declared 2026-09-23):
#          fdrQQ        CLR edges passing BH at q = QQ/100 (permutation null)
#          fdr05_topK   fdr05 plus each gene's own top-K CLR neighbours
#                       (directed: row i attends to its K highest S[i, ])
#          softK        softmax over each row's top-50 CLR scores, per-row
#                       temperature giving K effective neighbours
#          pearson_topK control: directed top-K by |Pearson|, weights |r|
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
# SECONDARY (continuity with Faith, Hayete et al. 2007): TF-node edge PR, where the TF's
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
            set = "avg", B = 100L, threads = NULL, alpha = 0.5,
            depths = "0,1,2,3,5,8,12,20", min_size = 5L, max_size = 500L,
            quick = 0L, mi_null = FALSE, out = NULL, seed = 20260922L,
            primary = "none_hg", bin_sweep = "6,8,12,16", reuse = "",
            operators = "fdr05,fdr20,fdr05_top5,fdr05_top10,soft10,pearson_top5,pearson_top10")
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
  none_hg           = list(transform = "none", bins = "hg", combine = "stouffer"),
  none_hg_eu        = list(transform = "none", bins = "hg", combine = "euclidean"),
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
  if (!is.null(reuse$null)) {
    fit$restore_threshold(reuse$null, reuse$params$threshold$tau,
                          reuse$params$threshold$q %||% 0.05)
  } else fit$restore_threshold(list(suf = 0, w = 1, n_null = 0, nbins = 1,
                                    M = 1, statistic = "clr"),
                               reuse$params$threshold$tau)
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

## ---- 4. diffusion depth sweep (Design A: fixed operators) ---------------------
# Every pre-declared operator gets the same sweep: depths t, value modes
# raw / signed / conditional, readout |cor| of diffused profiles; plus the
# value-independent attention mass (P^t + P^t')/2. Nothing here is chosen by
# RegulonDB; all cells are reported.
R0 <- tcrossprod(Z) / (ncol(Z) - 1)           # Pearson, for signs and control
q_base <- fit$params$threshold$q %||% 0.05
build_op <- function(name) {
  if (grepl("^fdr[0-9]+$", name)) {
    fit$reselect(as.numeric(sub("^fdr", "", name)) / 100)
    fit$build_operator(alpha = opt$alpha)
  } else if (grepl("^fdr[0-9]+_top[0-9]+$", name)) {
    fit$reselect(as.numeric(sub("^fdr([0-9]+)_top.*$", "\\1", name)) / 100)
    fit$build_operator(alpha = opt$alpha,
                       topk_union = as.integer(sub("^.*_top", "", name)))
  } else if (grepl("^soft[0-9]+$", name)) {
    fit$build_operator(alpha = opt$alpha,
                       softmax_keff = as.numeric(sub("^soft", "", name)),
                       softmax_cap = 50L)
  } else if (grepl("^pearson_top[0-9]+$", name)) {
    k <- as.integer(sub("^pearson_top", "", name))
    Ar <- abs(R0); diag(Ar) <- 0
    fit$set_operator(clr:::.topk_rows(Ar, k), alpha = opt$alpha, label = name)
  } else stop("unknown operator: ", name)
  invisible(fit$operator)
}
readout <- function(E) {
  Ez <- E - rowMeans(E)
  sdv <- sqrt(rowSums(Ez^2) / (ncol(Ez) - 1)); sdv[sdv <= 0] <- 1
  Ez <- Ez / sdv
  abs(tcrossprod(Ez) / (ncol(Ez) - 1))
}
op_rows <- list(); depth_rows <- list()
Ct0 <- readout(Z)
add_scores(score_all(Ct0, "op=none;values=none;t=00;readout=abscor"))
headline("op=none;values=none;t=00;readout=abscor")
depth_rows[[1]] <- data.frame(operator = "none", values = "none", t = 0L,
                              median_abs_cor_offdiag = stats::median(Ct0[upper.tri(Ct0)]))
rm(Ct0); invisible(gc())
ops <- strsplit(opt$operators, ",")[[1]]
if (!is.null(reuse) && is.null(reuse$null)) {
  drop <- grepl("^fdr", ops) & !grepl("^fdr05", ops)
  if (any(drop)) say("reused fit has no stored null: skipping ",
                     paste(ops[drop], collapse = ", "))
  ops <- ops[!drop]
}
for (op in ops) {
  A <- build_op(op)
  off <- A > 0 & row(A) != col(A)
  op_rows[[length(op_rows) + 1L]] <- data.frame(
    operator = op, edges = sum(off), mean_out_degree = mean(rowSums(off)),
    isolated_genes = sum(rowSums(off) == 0),
    neg_cor_edge_share = mean(R0[off] < 0), alpha = opt$alpha)
  say(sprintf("operator %-14s edges %6d, mean out-degree %5.1f, isolated %4d, negative-cor share %.4f",
              op, sum(off), mean(rowSums(off)), sum(rowSums(off) == 0),
              mean(R0[off] < 0)))
  # attention mass (value-independent)
  P <- Matrix::Matrix((1 - opt$alpha) * diag(G) + opt$alpha * A, sparse = TRUE)
  Pt <- diag(G)
  for (t in seq_len(max(depths))) {
    Pt <- as.matrix(Pt %*% P)
    if (t %in% depths) {
      Sa <- (Pt + t(Pt)) / 2; diag(Sa) <- 0
      lab <- sprintf("op=%s;values=none;t=%02d;readout=attention", op, t)
      add_scores(score_all(Sa, lab)); headline(lab)
    }
  }
  rm(Pt, Sa, P); invisible(gc())
  for (vt in c("raw", "signed", "conditional")) {
    tic(sprintf("diffuse_%s_%s", op, vt),
        fit$diffuse(steps = max(depths), values = vt))
    tr <- fit$trajectory
    for (t in depths[depths > 0]) {
      Ct <- readout(tr[[t + 1L]])
      lab <- sprintf("op=%s;values=%s;t=%02d;readout=abscor", op, vt, t)
      add_scores(score_all(Ct, lab)); headline(lab)
      depth_rows[[length(depth_rows) + 1L]] <- data.frame(
        operator = op, values = vt, t = t,
        median_abs_cor_offdiag = stats::median(Ct[upper.tri(Ct)]))
      rm(Ct); invisible(gc())
    }
    fit$reset_diffusion(); rm(tr); invisible(gc())
  }
}
wcsv(do.call(rbind, op_rows), "operators.csv")
# restore the primary operator for the saved artifacts
build_op(ops[1]); A <- Matrix::Matrix(fit$operator, sparse = TRUE)

cm_all <- do.call(rbind, cm_rows); coh_all <- do.call(rbind, coh_rows)
parse_lab <- function(d) {
  m <- grepl("^op=", d$method)
  d <- d[m, ]
  kv <- strsplit(d$method, ";", fixed = TRUE)
  get <- function(key) vapply(kv, function(x) sub(paste0("^", key, "="), "",
                                                  grep(paste0("^", key, "="), x, value = TRUE)), "")
  cbind(operator = get("op"), values = get("values"), t = as.integer(get("t")),
        readout = get("readout"), d)
}
dep_cm <- parse_lab(cm_all); dep_coh <- parse_lab(coh_all)
wcsv(dep_cm, "depth_comembership.csv")
wcsv(dep_coh, "depth_coherence.csv")
wcsv(do.call(rbind, depth_rows), "depth_collapse.csv")

# Figure: co-membership AUPR (SC, all classes) vs depth, one panel per value
# mode, one line per operator (|cor| readout); horizontal reference = CLR of
# the primary configuration.
clr_ref <- cm_all$aupr[cm_all$method == paste0(opt$primary, ":clr") &
                         cm_all$stratum == "all" & cm_all$evidence == "SC"]
d0 <- dep_cm[dep_cm$stratum == "all" & dep_cm$evidence == "SC" &
               dep_cm$readout == "abscor", ]
png(file.path(opt$out, "depth_sweep.png"), width = 1500, height = 520, res = 110)
par(mfrow = c(1, 3), mar = c(5, 4.5, 3, 1))
cols <- stats::setNames(grDevices::hcl.colors(length(ops), "Dark 3"), ops)
base0 <- d0$aupr[d0$operator == "none"]
yl <- range(c(d0$aupr, clr_ref), na.rm = TRUE)
for (vt in c("raw", "signed", "conditional")) {
  plot(NA, xlim = c(0, max(depths)), ylim = yl, xlab = "diffusion depth t",
       ylab = "co-membership AUPR (SC)", main = paste("values =", vt))
  abline(h = clr_ref, lty = 2, col = "grey40")
  for (op in ops) {
    d <- d0[d0$operator == op & d0$values == vt, ]
    d <- d[order(d$t), ]
    lines(c(0, d$t), c(base0, d$aupr), type = "b", pch = 16, cex = 0.7,
          col = cols[op])
  }
  if (vt == "raw") legend("bottomright", bty = "n", cex = 0.8, lty = 1,
                          col = cols, legend = ops)
}
dev.off()

## ---- 5. artifacts ------------------------------------------------------------
saveRDS(list(mi = fit$mi, clr = fit$clr_scores, edges = E_sel,
             null = tryCatch(fit$null_distribution, error = function(e) NULL),
             operator = A, params = fit$params, genes = genes,
             seed = opt$seed, options = opt),
        file.path(opt$out, "primary_fit.rds"), compress = "gzip")
wcsv(data.frame(stage = names(timings), seconds = unlist(timings)), "timings.csv")
writeLines(capture.output(sessionInfo()), file.path(opt$out, "sessionInfo.txt"))
say("done; outputs in ", normalizePath(opt$out))
