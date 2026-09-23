#!/usr/bin/env Rscript
# Perturbation-target identification on M3D (466 replicate-averaged
# experiments) by the delete-k jackknife influence of the perturbation's
# experiments on the network (user-approved design, 2026-09-23).
#
# For a perturbation P (a regulator knocked out / over-expressed / mutated in
# k experiments), each gene's influence is how much its network row changes
# when P's k experiments are removed:
#     I_g(P) = RMS_j ( S_all[g, j] - S_{-P}[g, j] )
# for a network S (|Pearson|, CLR, or CLR-attention mass at the held-out depth
# t = 8). Influence is calibrated per gene against B random k-subsets of
# experiments (the delete-k jackknife null, shared across perturbations of
# the same k): z_g = (I_g - mean_null_g) / sd_null_g. Genes with large z are
# predicted targets. No parametric (Gaussian / Hotelling) assumption is made;
# k enters through the matched null.
#
# Baselines: differential expression (mean over P's experiments of |z| of the
# gene's value against all other experiments -- the "trivial gene-level
# feature" of Kendiukhov 2026), and gene variance (a P-independent prior).
# Pre-declared combination: Stouffer of rank-normalized DE and CLR-attention
# influence. Reference method for this task on M3D: SSEM-Lasso network
# filtering (Cosgrove, Zhou, Gardner & Kolaczyk, Bioinformatics 2008,
# doi:10.1093/bioinformatics/btn476), which reports sensitivity among the top
# 100 ranked genes; that metric is reported here too (not yet re-run).
#
# Truth: the RegulonDB regulon of the perturbed regulator, mapped by a table
# fixed in advance (gene -> regulator, including relA -> ppGpp,
# rpoS -> sigma38 and recA -> LexA/SOS). The perturbed genes themselves are
# excluded from scoring.
#
# Usage (repo root):
#   Rscript analysis/perturbation.R [--B 30] [--out results/perturbation]
#       [--threads N] [--quick 0] [--workers W] [--mem_gb 4]
# The null draws and the per-perturbation fits run in parallel (foreach over
# local cores minus 2, capped by RAM; analysis/parallel.R). The delete-k
# subsets are drawn in the master in the original order, so they match the
# sequential version; each fit then runs under its own seed (used only for
# the 1e-9 jitter of a gene that is constant in a subset).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            B = 30L, out = "results/perturbation", threads = NULL,
            quick = 0L, seed = 20260926L, alpha = 0.5, t_att = 8L,
            workers = NULL, mem_gb = 4)
ints <- c("B", "threads", "quick", "seed", "t_att", "workers")
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key %in% ints) as.integer(args[i + 1L]) else
    if (key == "mem_gb") as.numeric(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "perturbation.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))
set.seed(opt$seed)

## ---- data ------------------------------------------------------------------
base <- basename(normalizePath(opt$m3d))
raw <- utils::read.delim(file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (opt$quick > 0L) {
  v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]
}
genes <- rownames(X); G <- length(genes); E <- colnames(X)
to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)

feat <- utils::read.delim(file.path(opt$m3d, paste0(base, ".experiment_feature_descriptions")),
                          stringsAsFactors = FALSE, quote = "")
pg <- feat[feat$feature_name == "perturbation_gene", c("experiment_name", "value")]
pt <- feat[feat$feature_name == "perturbation", c("experiment_name", "value")]
names(pt)[2] <- "type"
pg <- merge(pg, pt, all.x = TRUE)
pg <- pg[pg$experiment_name %in% E, ]

# Pre-declared gene -> regulator mapping (regulator names as in RISet /
# NetworkSigmaGene). Multi-gene entries are split on ",", "and".
reg_map <- c(rpoS = "sigma:sigma38", relA = "compound:ppGpp", hupA = "TF:HU",
             hupB = "TF:HU", fis = "TF:Fis", fnr = "TF:FNR", soxS = "TF:SoxS",
             oxyR = "TF:OxyR", crp = "TF:CRP", arcA = "TF:ArcA", appY = "TF:AppY",
             ryhB = "sRNA:RyhB", recA = "TF:LexA", lexA = "TF:LexA",
             marA = "TF:MarA", rob = "TF:Rob", lrp = "TF:Lrp", hns = "TF:H-NS",
             cpxR = "TF:CpxR", gadX = "TF:GadX", lacI = "TF:LacI", dksA = "TF:DksA",
             rpoH = "sigma:sigma32", rpoE = "sigma:sigma24", rpoN = "sigma:sigma54",
             fliA = "sigma:sigma28", ompR = "TF:OmpR", phoB = "TF:PhoB")
pg$genes <- lapply(strsplit(pg$value, ",|\\band\\b"), trimws)
pg$keys <- lapply(pg$genes, function(g) unique(unname(reg_map[g[g %in% names(reg_map)]])))
pg <- pg[lengths(pg$keys) > 0, ]
grp_key <- vapply(pg$keys, function(k) paste(sort(k), collapse = "+"), "")
groups <- split(pg$experiment_name, grp_key)
pgenes <- lapply(split(pg$genes, grp_key), function(l) unique(unlist(l)))
rk <- paste(rtab$reg_class, rtab$regulator, sep = ":")
truth <- lapply(names(groups), function(k) {
  ks <- strsplit(k, "+", fixed = TRUE)[[1]]
  tg <- unique(rtab$bnumber[rk %in% ks & !is.na(rtab$bnumber)])
  lab <- as.integer(genes %in% tg)
  self <- genes %in% to_bn(pgenes[[k]])
  list(label = lab, exclude = self)
})
names(truth) <- names(groups)
ok <- vapply(truth, function(t) sum(t$label[!t$exclude]) >= 5, NA)
groups <- groups[ok]; truth <- truth[ok]
say(sprintf("%d genes x %d experiments; %d perturbation groups with >= 5 mapped targets:",
            G, length(E), length(groups)))
for (k in names(groups))
  say(sprintf("   %-26s k = %d experiments (%s); %d targets on array", k,
              length(groups[[k]]),
              paste(unique(pg$type[pg$experiment_name %in% groups[[k]]]), collapse = "/"),
              sum(truth[[k]]$label)))

## ---- networks ------------------------------------------------------------
Zs <- function(M) { Z <- M - rowMeans(M); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1))
                    s[s <= 0] <- 1; Z / s }
OMP_THREADS <- NULL                  # master: bspline_mi default (cores - 2)
nthr <- function() if (is.null(opt$threads)) OMP_THREADS else opt$threads
networks <- function(cols) {
  Y <- X[, cols, drop = FALSE]
  rng <- apply(Y, 1, function(x) diff(range(x)))
  if (any(rng <= 0))
    Y[rng <= 0, ] <- Y[rng <= 0, ] + stats::rnorm(sum(rng <= 0) * ncol(Y), sd = 1e-9)
  Z <- Zs(Y)
  pear <- abs(tcrossprod(Z) / (ncol(Z) - 1)); rm(Z)
  f <- ClrAttention$new(Y)$estimate_mi(bins = "hg", transform = "none",
                                       threads = nthr())
  rm(Y)
  f$calibrate(method = "normal", combine = "stouffer")
  f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
  # Sparse lazy operator built directly (identical values to the dense
  # construction), so no dense (1-a)I + aA temporary.
  P <- Matrix::Diagonal(G, 1 - opt$alpha) +
    opt$alpha * Matrix::Matrix(f$operator, sparse = TRUE)
  Pt <- as.matrix(P)
  if (opt$t_att > 1L) for (t in 2:opt$t_att) Pt <- as.matrix(Pt %*% P)
  rm(P)
  A <- Pt + t(Pt); rm(Pt); A <- A / 2
  out <- list(pearson = pear, clr = f$clr_scores, clr_attention = A)
  f$release(); rm(f, pear, A)
  lapply(out, function(s) { diag(s) <- 0; s })
}
influence <- function(full, part) vapply(names(full), function(m)
  sqrt(rowMeans((full[[m]] - part[[m]])^2)), numeric(G))

t0 <- Sys.time()
full <- networks(E)
say(sprintf("full networks built in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Task list: B random delete-k draws per k (drawn here, in the original
# order), then one delete-P fit per perturbation group.
ks <- sort(unique(lengths(groups)))
tasks <- list()
for (k in ks) for (b in seq_len(opt$B))
  tasks[[length(tasks) + 1L]] <- list(kind = "null", k = k, b = b,
                                      keep = setdiff(E, sample(E, k)))
for (g in names(groups))
  tasks[[length(tasks) + 1L]] <- list(kind = "pert", g = g,
                                      keep = setdiff(E, groups[[g]]))
say(sprintf("%d network fits (%d null draws for k in {%s}, %d perturbations)",
            length(tasks), opt$B * length(ks), paste(ks, collapse = ","),
            length(groups)))
source(file.path(script_dir, "parallel.R"))
pc <- par_start(mem_gb = opt$mem_gb, workers = opt$workers, say = say)
t0 <- Sys.time()
fits <- foreach(task = tasks, ti = seq_along(tasks), .errorhandling = "stop") %dopar% {
  if (par_should_stop()) stop("stop file jobs/control/stop present")
  set.seed(opt$seed + ti)
  I <- influence(full, part <- networks(task$keep)); rm(part); par_release()
  I
}
par_stop(pc)
say(sprintf("fits done in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
kind <- vapply(tasks, `[[`, "", "kind")
null <- list()
for (k in ks) {
  sel <- which(kind == "null" & vapply(tasks, function(t) if (is.null(t$k)) NA_integer_ else t$k, 0L) == k)
  arr <- array(NA_real_, c(G, length(full), opt$B), dimnames = list(genes, names(full), NULL))
  for (j in seq_along(sel)) arr[, , j] <- fits[[sel[j]]]
  null[[as.character(k)]] <- arr
}
pert_I <- stats::setNames(fits[kind == "pert"],
                          vapply(tasks[kind == "pert"], `[[`, "", "g"))
rm(fits); invisible(gc())

## ---- scores and evaluation ---------------------------------------------------
qn <- function(x) stats::qnorm((rank(x) - 0.5) / length(x))
res <- list(); scores_out <- list()
mu_ref <- rowMeans(X)
for (k in names(groups)) {
  cols <- groups[[k]]
  ref <- setdiff(E, cols)
  mu <- rowMeans(X[, ref]); sdv <- apply(X[, ref], 1, stats::sd); sdv[sdv <= 0] <- 1
  de <- rowMeans(abs((X[, cols, drop = FALSE] - mu) / sdv))
  I <- pert_I[[k]]
  nl <- null[[as.character(length(cols))]]
  zI <- (I - apply(nl, c(1, 2), mean)) / pmax(apply(nl, c(1, 2), stats::sd), 1e-12)
  sc <- list(differential_expression = de, gene_variance = apply(X, 1, stats::var),
             influence_pearson = zI[, "pearson"], influence_clr = zI[, "clr"],
             influence_clr_attention = zI[, "clr_attention"])
  sc$combined_de_clr_attention <- (qn(de) + qn(zI[, "clr_attention"])) / sqrt(2)
  tr <- truth[[k]]; use <- !tr$exclude
  for (m in names(sc)) {
    ps <- pr_summary(sc[[m]][use], tr$label[use])
    top100 <- order(-sc[[m]][use])[1:100]
    res[[length(res) + 1L]] <- data.frame(perturbation = k, k = length(cols),
                                          n_targets = sum(tr$label[use]),
                                          method = m, auroc = ps[["auroc"]],
                                          aupr_ratio = ps[["aupr"]] / mean(tr$label[use]),
                                          # Cosgrove et al. 2008 (SSEM-Lasso on M3D)
                                          # report sensitivity among the top 100 genes
                                          sens_top100 = sum(tr$label[use][top100]) /
                                            sum(tr$label[use]))
  }
  scores_out[[k]] <- as.data.frame(sc)
  cur <- do.call(rbind, res); cur <- cur[cur$perturbation == k, ]
  say(sprintf("%-26s %s", k, paste(sprintf("%s %.3f", sub("influence_", "I_", cur$method),
                                           cur$auroc), collapse = " | ")))
}
res <- do.call(rbind, res)
utils::write.csv(res, file.path(opt$out, "perturbation_metrics.csv"), row.names = FALSE)
saveRDS(list(scores = scores_out, genes = genes, groups = groups),
        file.path(opt$out, "perturbation_scores.rds"))
de_ref <- res[res$method == "differential_expression", c("perturbation", "auroc", "aupr_ratio")]
names(de_ref)[2:3] <- c("auroc_de", "apr_de")
m <- merge(res, de_ref)
summ <- do.call(rbind, lapply(split(m, m$method), function(x) data.frame(
  method = x$method[1], n = nrow(x), median_auroc = stats::median(x$auroc),
  median_aupr_ratio = stats::median(x$aupr_ratio),
  mean_sens_top100 = mean(x$sens_top100),
  wins_vs_de_auroc = sum(x$auroc > x$auroc_de),
  wilcoxon_p_vs_de = if (all(x$auroc == x$auroc_de)) NA_real_ else
    suppressWarnings(stats::wilcox.test(x$auroc, x$auroc_de, paired = TRUE)$p.value))))
summ <- summ[order(-summ$median_auroc), ]
utils::write.csv(summ, file.path(opt$out, "summary.csv"), row.names = FALSE)
say("summary over perturbations (AUROC; paired vs differential expression):")
for (j in seq_len(nrow(summ)))
  say(sprintf("  %-28s median AUROC %.3f  median AUPR ratio %.2f  mean sens@100 %.3f  wins vs DE %d/%d  p = %.3g",
              summ$method[j], summ$median_auroc[j], summ$median_aupr_ratio[j],
              summ$mean_sens_top100[j],
              summ$wins_vs_de_auroc[j], summ$n[j], summ$wilcoxon_p_vs_de[j]))
say("done")
