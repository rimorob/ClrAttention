#!/usr/bin/env Rscript
# E. coli TF-node benchmark for the attention operators (466 replicate-averaged
# M3D experiments). This is the Faith, Hayete et al. 2007 convention: a
# candidate edge is (TF gene, other gene), scored by the TF's own mRNA, and it
# is positive if RegulonDB lists that TF as regulating that gene. It asks
# "is X a direct regulator of Y?", the same question BEELINE asks.
#
# Depths are a diagnostic scan (t = 1 ... 32), reported next to the held-out
# t* = 8; nothing here is used to choose t.
#
# Outputs (results/tfnode_attention/): tfnode_pooled.csv (pooled AUPR/AUROC
# over all candidate edges, SC and all evidence) and tfnode_per_tf.csv
# (per-TF AUROC of the TF's row: known targets vs other genes).
#
# Usage (repo root):  Rscript analysis/tfnode_attention.R [--out results/tfnode_attention]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(m3d = "data/E_coli_v4_Build_6", rdb = "data/RegulonDBExtract",
            out = "results/tfnode_attention", alpha = 0.5, threads = NULL)
i <- 1L
while (i <= length(args)) {
  key <- gsub("-", "_", sub("^--", "", args[i]))
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- if (key == "threads") as.integer(args[i + 1L]) else args[i + 1L]
  i <- i + 2L
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "tfnode.log")
say <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
script_dir <- (function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else "analysis"
})()
suppressPackageStartupMessages({ library(clr); library(Matrix) })
invisible(compiler::enableJIT(3L))
source(file.path(script_dir, "regulondb.R"))
source(file.path(script_dir, "regulons.R"))

base <- basename(normalizePath(opt$m3d))
raw <- utils::read.delim(file.path(opt$m3d, paste0("avg_", base, "_exps466probes4297.tab")),
                         check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); storage.mode(X) <- "double"
rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
genes <- rownames(X); G <- length(genes)
to_bn <- make_bnumber_mapper(file.path(script_dir, "ecoli_k12_genes.tsv"))
rtab <- build_regulon_table(opt$rdb, to_bn)

tfr <- rtab[rtab$reg_class == "TF" & !is.na(rtab$bnumber), ]
tf_parts <- strsplit(.classic_tf_to_genes(tfr$regulator), ";", fixed = TRUE)
tfnet <- data.frame(tf_gene = to_bn(unlist(tf_parts)),
                    target = rep(tfr$bnumber, lengths(tf_parts)),
                    confidence = rep(tfr$confidence, lengths(tf_parts)),
                    stringsAsFactors = FALSE)
tfnet <- tfnet[!is.na(tfnet$tf_gene), ]
U <- list(SC = edge_universe(tfnet, genes, conf_keep = c("S", "C")),
          all = edge_universe(tfnet, genes))
say(sprintf("%d genes x %d experiments; TF-node universe: SC %d TFs, %d positives of %d pairs; all %d TFs, %d positives",
            G, ncol(X), length(U$SC$tfs), U$SC$n_pos, length(U$SC$label),
            length(U$all$tfs), U$all$n_pos))

pooled <- list(); per_tf <- list()
score <- function(S, m) {
  diag(S) <- 0
  for (ev in names(U)) {
    u <- U[[ev]]
    ps <- pr_summary(S[cbind(u$i, u$j)], u$label)
    pooled[[length(pooled) + 1L]] <<- data.frame(method = m, evidence = ev,
      aupr = ps[["aupr"]], aupr_ratio = ps[["aupr"]] / mean(u$label),
      auroc = ps[["auroc"]])
    # per TF: its row, known targets (either direction) vs all other genes
    au <- vapply(names(u$targets), function(tn) {
      t <- as.integer(tn); y <- integer(G); y[unique(u$targets[[tn]])] <- 1L; y[t] <- NA
      if (sum(y, na.rm = TRUE) < 3) return(NA_real_)
      s <- S[t, ]; ok <- !is.na(y)
      r <- rank(s[ok]); n1 <- sum(y[ok]); n0 <- sum(ok) - n1
      (sum(r[y[ok] == 1]) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0)
    }, numeric(1))
    per_tf[[length(per_tf) + 1L]] <<- data.frame(method = m, evidence = ev,
      tf = genes[as.integer(names(u$targets))], n_targets = lengths(u$targets),
      auroc = au)
  }
  invisible(gc(verbose = FALSE))
}
attention <- function(A, ts, fn) {
  P <- Matrix::Diagonal(G, 1 - opt$alpha) + opt$alpha * Matrix::Matrix(A, sparse = TRUE)
  Pt <- NULL
  for (t in seq_len(max(ts))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts) { S <- Pt + t(Pt); S <- S / 2; fn(S, t); rm(S) }
  }
}

t0 <- Sys.time()
Z <- X - rowMeans(X); Z <- Z / sqrt(rowSums(Z^2) / (ncol(Z) - 1))
R <- tcrossprod(Z) / (ncol(Z) - 1); rm(Z)
score(abs(R), "abs_pearson")
Ar <- abs(R); diag(Ar) <- 0; Apear <- clr:::.topk_rows(Ar, 10L); rm(Ar, R)
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
score(f$mi, "mi")
f$calibrate(method = "normal", combine = "stouffer")
score(f$clr_scores, "clr")
f$build_operator(alpha = opt$alpha, softmax_keff = 10, softmax_cap = 50L)
attention(f$operator, c(1, 2, 4, 8, 16, 32), function(S, t) score(S, sprintf("soft10_t%02d", t)))
f$set_operator(Apear, alpha = opt$alpha, label = "pearson_top10")
attention(f$operator, c(2, 6), function(S, t) score(S, sprintf("pearson_top10_t%02d", t)))
f$release()

pooled <- do.call(rbind, pooled); per_tf <- do.call(rbind, per_tf)
utils::write.csv(pooled, file.path(opt$out, "tfnode_pooled.csv"), row.names = FALSE)
utils::write.csv(per_tf, file.path(opt$out, "tfnode_per_tf.csv"), row.names = FALSE)
clr_tf <- per_tf[per_tf$method == "clr", c("evidence", "tf", "auroc")]
names(clr_tf)[3] <- "auroc_clr"
m <- merge(per_tf, clr_tf)
for (ev in names(U)) {
  say(sprintf("== %s evidence: pooled AUPR ratio (x random), pooled AUROC, median per-TF AUROC, TFs beating CLR", ev))
  for (k in which(pooled$evidence == ev)) {
    x <- m[m$evidence == ev & m$method == pooled$method[k] & !is.na(m$auroc), ]
    say(sprintf("   %-18s AUPR ratio %5.2f  AUROC %.3f  per-TF median AUROC %.3f  beats CLR in %d/%d TFs",
                pooled$method[k], pooled$aupr_ratio[k], pooled$auroc[k],
                stats::median(x$auroc), sum(x$auroc > x$auroc_clr), nrow(x)))
  }
}
say(sprintf("done in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
