#!/usr/bin/env Rscript
# GRNBoost2 diagnostic, scored with the main comparison's evaluation (no new
# truth or metric). Inputs: the variant edge files of analysis/grnboost2_check.R
# (results/grnboost2_check/<case>/<variant>_s<seed>.tsv).
#
# M3D: analysis/baselines.R's own definitions, evaluated from its source (data,
#   truths, delete-half gene jackknife splits with its seed, evaluate(),
#   sym_max(), hub_stratify(), and the walk's stationary distribution for the
#   hub strata): pooled co-membership AUPR / base rate (RegulonDB SC and all
#   evidence), raw and hub-stratified, as in Table 2 / sota_table.md. Each
#   variant replaces GRNBoost2 in its row; paired jackknife differences against
#   the stored replicates of the other methods (results/baselines_m3d/
#   jack_replicates.csv, same splits) and the Bonferroni tie rule of
#   analysis/sota_table.R. Per regulon: regulon_ap_table() (min_pos 10, as
#   macro_aupr.R's regulon set), mean AP difference vs attention, descriptive.
#   Reproduction gate: the reported edges (results/baselines_m3d/
#   grnboost2_edges.tsv) must give the reported AUPR, jackknife replicates and
#   hub-stratified AUPR, or the script stops.
# BEELINE: per case, mean over BEELINE's primary truths of log2(AUPRC ratio /
#   CLR's) for each method, directed GRNBoost2 scores as beeline_baselines.R
#   (variant scores from results/grnboost2_check/metrics.csv, which used
#   beeline_baselines.R's score_truth; references from
#   results/beeline_baselines/metrics.csv).
#
# Usage (repo root; light: one thread, about 15 min, < 3 GB):
#   nice -n 19 Rscript analysis/grnboost2_check_score.R [--out results/grnboost2_check]

a0 <- commandArgs(trailingOnly = TRUE)
gout <- if (length(a0) >= 2 && a0[1] == "--out") a0[2] else "results/grnboost2_check"
work <- file.path(gout, "score_work"); dir.create(work, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

## ---- BEELINE (cheap; first) ---------------------------------------------------------------
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
bm <- utils::read.csv("results/beeline_baselines/metrics.csv", stringsAsFactors = FALSE)
bm$case <- paste0(bm$dataset, "_", bm$N)
gm <- utils::read.csv(file.path(gout, "metrics.csv"), stringsAsFactors = FALSE)
gm <- gm[gm$case != "m3d" & gm$orientation == "directed" & gm$truth %in% primary, ]
brows <- list()
for (cs in unique(gm$case)) {
  ref <- bm[bm$case == cs & bm$truth %in% primary, ]
  clr <- stats::setNames(ref$auprc_ratio[ref$method == "clr"], ref$truth[ref$method == "clr"])
  clr_e <- stats::setNames(ref$epr[ref$method == "clr"], ref$truth[ref$method == "clr"])
  # gate: the diagnostic's grnboost2_reported must equal beeline_baselines' grnboost2
  rep_d <- gm[gm$case == cs & gm$variant == "grnboost2_reported", ]
  rep_b <- ref[ref$method == "grnboost2", ]
  if (max(abs(rep_d$auprc_ratio - rep_b$auprc_ratio[match(rep_d$truth, rep_b$truth)])) > 1e-9)
    stop("BEELINE reproduction gate failed for ", cs)
  add <- function(method, seed, x) brows[[length(brows) + 1L]] <<- data.frame(case = cs, method = method, seed = seed,
    n_truths = nrow(x), mean_auprc_ratio = mean(x$auprc_ratio),
    log2_auprc_vs_clr = mean(log2(x$auprc_ratio / clr[x$truth])),
    log2_epr_vs_clr = mean(log2(pmax(x$epr, 1e-9) / pmax(clr_e[x$truth], 1e-9))))
  for (m in unique(ref$method)) add(m, NA, ref[ref$method == m, ])
  v <- gm[gm$case == cs & !gm$variant %in% c("abs_pearson", "grnboost2_reported"), ]
  vs <- unique(v[, c("variant", "seed")])
  for (k in seq_len(nrow(vs))) {
    x <- v[v$variant == vs$variant[k] & v$seed == vs$seed[k], ]; add(paste0("grnboost2_", vs$variant[k]), vs$seed[k], x) }
}
bres <- do.call(rbind, brows)
utils::write.csv(bres, file.path(gout, "beeline_vs_clr.csv"), row.names = FALSE)
cat("BEELINE written\n")

## ---- M3D: baselines.R's definitions ---------------------------------------------------------
src <- parse("analysis/baselines.R", keep.source = TRUE)
dp <- vapply(src, function(e) deparse(e)[1], "")
stop_at <- which(dp == "E <- list()")[1]
args <- c("--dataset", "m3d", "--workers", "1", "--out", work)
for (k in seq_len(stop_at - 1L)) {
  if (dp[k] == "args <- commandArgs(trailingOnly = TRUE)") next
  if (startsWith(dp[k], "if (!file.exists(tsv))")) next          # no expression export
  eval(src[[k]], globalenv())
}
# hub strata: the walk's stationary distribution, as in baselines.R (cached)
hf <- file.path(work, "H_m3d.rds")
if (file.exists(hf)) H <- readRDS(hf) else {
  want <- c("^f <- ClrAttention\\$new", "^f\\$build_operator", "^P <- Matrix::Diagonal", "^pi <- rep",
            "^eps <- ", "^for \\(it in 1:20000\\)", "^H <- outer\\(q, q")
  idx <- sort(unlist(lapply(want, function(w) which(grepl(w, dp) & seq_along(dp) > stop_at))))
  if (length(idx) != length(want)) stop("hub-strata expressions not found in baselines.R: ", paste(dp[idx], collapse = " / "))
  for (k in idx) eval(src[[k]], globalenv())
  saveRDS(H, hf); rm(P); f$release()
}
say("H ready")

jr <- utils::read.csv("results/baselines_m3d/jack_replicates.csv", stringsAsFactors = FALSE)
mr <- utils::read.csv("results/baselines_m3d/metrics.csv", stringsAsFactors = FALSE)
pr <- utils::read.csv("results/baselines_m3d/per_regulon.csv", stringsAsFactors = FALSE)
others <- c("attention", "clr", "genie3", "pidc", "abs_pearson", "hub")

read_edges <- function(f) utils::read.delim(f, colClasses = c("character", "character", "numeric"),
                                            na.strings = character(0))
score_one <- function(f, variant, seed) {
  e <- read_edges(f); S <- sym_max(edges_to_matrix(e$TF, e$target, e$importance)); rm(e)
  nz <- mean(S[upper.tri(S)] != 0)
  ev <- list(raw = evaluate(S, macro = FALSE), hubstrat = evaluate(hub_stratify(S, H), macro = FALSE))
  preg <- lapply(names(bench), function(nm) { b <- bench[[nm]]
    cbind(truth = nm, regulon_ap_table(S, b$reg, b$co, b$operon, genes)) })
  rm(S); invisible(gc())
  rows <- list()
  for (sc in names(ev)) for (nm in names(bench)) {
    e1 <- ev[[sc]][[nm]]; suf <- if (sc == "hubstrat") "_hubstrat" else ""
    ms <- if (sc == "hubstrat") setdiff(others, "hub") else others
    br <- mean(bench[[nm]]$pairs$label)
    a <- c(grnboost2 = e1$aupr, stats::setNames(mr$aupr[match(paste0(ms, suf, nm), paste0(mr$method, mr$truth))], ms))
    best <- names(which.max(a)); zc <- stats::qnorm(1 - 0.025 / (length(a) - 1L))
    jk <- function(m) { y <- jr[jr$method == paste0(m, suf) & jr$truth == nm, ]; y$aupr[order(y$rep)] }
    sdd <- function(m) stats::sd(e1$jack - jk(m))
    z_best <- if (best == "grnboost2") NA_real_ else (a[[best]] - a[["grnboost2"]]) / stats::sd(jk(best) - e1$jack)
    r <- data.frame(case = "m3d", variant = variant, seed = seed, scoring = sc, truth = nm, nonzero_pairs = nz,
                    aupr = e1$aupr, base_rate = br, ratio = e1$aupr / br, coherence = e1$coherence,
                    rank_in_row = sum(a > a[["grnboost2"]]) + 1L, n_methods = length(a), best = best,
                    tie_with_best = if (is.na(z_best)) "best" else if (z_best < 1.96) "tie" else if (z_best < zc) "tie_bonferroni" else "no",
                    tf_aupr_ratio = if (is.null(e1$tf_aupr_ratio)) NA else e1$tf_aupr_ratio)
    for (m in setdiff(ms, "hub")) { d <- e1$aupr - a[[m]]; s <- sdd(m)
      r[[paste0("d_vs_", m)]] <- d; r[[paste0("d_", m, "_lo")]] <- d - 1.96 * s; r[[paste0("d_", m, "_hi")]] <- d + 1.96 * s }
    rows[[length(rows) + 1L]] <- r
  }
  out <- do.call(rbind, rows)
  pg <- do.call(rbind, preg)
  mac <- do.call(rbind, lapply(names(bench), function(nm) {
    att <- pr[pr$method == "attention" & pr$truth == nm, ]; x <- pg[pg$truth == nm, ]
    k <- match(x$regulon, att$regulon); ok <- !is.na(x$ap) & !is.na(k) & !is.na(att$ap[k])
    d <- x$ap[ok] - att$ap[k[ok]]
    data.frame(truth = nm, macro_n_regulons = sum(ok), macro_mean_ap = mean(x$ap[ok]),
               macro_mean_dap_vs_attention = mean(d), macro_median_dap_vs_attention = stats::median(d),
               macro_frac_regulons_ahead = mean(d > 0)) }))
  list(pooled = merge(out, mac, by = "truth", sort = FALSE), per_regulon = cbind(variant = variant, seed = seed, pg))
}

# reproduction gate on the reported GRNBoost2 network
rep <- score_one("results/baselines_m3d/grnboost2_edges.tsv", "reported", NA)
for (nm in names(bench)) {
  x <- rep$pooled[rep$pooled$truth == nm, ]
  for (sc in c("raw", "hubstrat")) {
    want <- mr$aupr[mr$method == paste0("grnboost2", if (sc == "hubstrat") "_hubstrat" else "") & mr$truth == nm]
    got <- x$aupr[x$scoring == sc]
    say(sprintf("gate %s %s: reported %.9f, rescored %.9f", nm, sc, want, got))
    if (abs(want - got) > 1e-9) stop("M3D reproduction gate failed (", nm, " ", sc, ")")
  }
}
pg0 <- rep$per_regulon; p0 <- pr[pr$method == "grnboost2", ]
k <- match(paste(pg0$truth, pg0$regulon), paste(p0$truth, p0$regulon))
if (max(abs(pg0$ap - p0$ap[k]), na.rm = TRUE) > 1e-9 || any(is.na(pg0$ap) != is.na(p0$ap[k])))
  stop("per-regulon reproduction gate failed")
say("reproduction gates passed")

res <- list(rep$pooled); preg <- list(rep$per_regulon)
for (f in sort(list.files(file.path(gout, "m3d"), pattern = "_s[0-9]+\\.tsv$", full.names = TRUE))) {
  b <- sub("\\.tsv$", "", basename(f)); v <- sub("_s[0-9]+$", "", b); sd0 <- as.integer(sub(".*_s", "", b))
  t0 <- Sys.time(); r <- score_one(f, v, sd0); res[[length(res) + 1L]] <- r$pooled; preg[[length(preg) + 1L]] <- r$per_regulon
  x <- r$pooled[r$pooled$scoring == "raw", ]
  say(sprintf("m3d %s seed %d: ratio SC %.3f, all %.3f (rank %s); %.1f min", v, sd0, x$ratio[x$truth == "SC"],
              x$ratio[x$truth == "all"], paste(x$rank_in_row, collapse = "/"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  utils::write.csv(do.call(rbind, res), file.path(gout, "m3d_main_eval.csv"), row.names = FALSE)
}
utils::write.csv(do.call(rbind, preg), file.path(gout, "m3d_per_regulon.csv"), row.names = FALSE)
say("done")
