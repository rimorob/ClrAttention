#!/usr/bin/env Rscript
# Scoreboard (docs/plans/scoreboard_plan.md), E. coli jobs. DESCRIPTIVE: point estimates, no intervals, no tests.
# Cells: pooled AUPR on RegulonDB SC / all (raw, hub pi, hub MI), direct TF -> target edge metrics (tf_aupr_ratio, tf_auroc_median),
# per-regulon mean AP, for CLR, attention (t = 9), a3r, a3r2, a3r3 (t = 9 in every layer), |r|. Equivalence gate first: attention, CLR and |r|
# reproduce the stored baselines outputs (results/baselines_<dataset>/metrics.csv, per_regulon.csv), a3r/a3r2/a3r3 the stored explore_stack
# pooled values (m3d) and endpoints_all.csv (ecomac a3r); the job aborts before writing results if any check fails.
# Usage (repo root): Rscript analysis/scoreboard/sb_ecoli.R --dataset m3d|ecomac [--out results/scoreboard/<dataset>] [--threads 8]
args <- commandArgs(TRUE); opt <- list(dataset = "m3d", out = NULL, threads = 8L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "threads") as.integer(args[k + 1]) else args[k + 1] }
if (is.null(opt$out)) opt$out <- file.path("results/scoreboard", opt$dataset)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
sb_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(sb_dir, "sb_lib.R")); sb_source()
say <- function(...) sb_say(file.path(opt$out, "run.log"), ...)
say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
DS <- opt$dataset
D <- sb_load_ecoli(DS); genes <- D$genes; G <- length(genes); bench <- D$bench; U <- D$U
say(sprintf("%s: %d genes x %d samples; truths %s", DS, G, ncol(D$X), paste(sprintf("%s (%d pairs, %d pos)", names(bench), vapply(bench, function(b) length(b$pairs$label), 0L), vapply(bench, function(b) sum(b$pairs$label), 0)), collapse = "; ")))
R <- abs(stats::cor(t(D$X))); R[is.na(R)] <- 0
Pr <- sb_prep(D$X, opt$threads, say); D$X <- NULL
SS <- c(Pr$S, list(abs_pearson = R)); rm(R); invisible(gc(verbose = FALSE))

pooled <- list(); tfm <- list(); perreg <- list(); prtab <- list()
for (m in names(SS)) {
  S <- SS[[m]]; S[!is.finite(S)] <- 0; diag(S) <- 0
  pooled[[m]] <- sb_pooled(m, S, bench, Pr$H_att, Pr$H_mi)
  for (nm in names(U)) { u <- U[[nm]]
    tfm[[length(tfm) + 1L]] <- data.frame(readout = m, truth = nm, tf_aupr_ratio = ap_fast(S[cbind(u$i, u$j)], u$label) / mean(u$label), tf_auroc_median = sb_per_tf_median_auroc(S, u, G), stringsAsFactors = FALSE) }
  for (nm in names(bench)) { b <- bench[[nm]]; x <- regulon_ap_table(S, b$reg, b$co, b$operon, genes); x <- cbind(readout = m, truth = nm, x); prtab[[length(prtab) + 1L]] <- x
    y <- x[!is.na(x$ap), ]
    perreg[[length(perreg) + 1L]] <- data.frame(readout = m, truth = nm, n_regulons = nrow(y), mean_ap = mean(y$ap), mean_base_rate = mean(y$base_rate), enrichment = mean(y$ap) / mean(y$base_rate),
      mean_ap_norm = mean((y$ap - y$base_rate) / (1 - y$base_rate)), median_fold = stats::median(y$ap / y$base_rate), stringsAsFactors = FALSE) }
  say("scored ", m); rm(S); invisible(gc(verbose = FALSE))
}
pooled <- do.call(rbind, pooled); tfm <- do.call(rbind, tfm); perreg <- do.call(rbind, perreg); prtab <- do.call(rbind, prtab)

## ---- equivalence gate -----------------------------------------------------------------------------------------------------
rows <- list()
bm <- utils::read.csv(sprintf("results/baselines_%s/metrics.csv", DS), stringsAsFactors = FALSE)
mp <- c(attention = "attention", CLR = "clr", abs_pearson = "abs_pearson")
for (m in names(mp)) {
  for (tn in c("SC", "all")) {
    mine <- pooled[pooled$readout == m & pooled$truth == tn, ]; bs <- bm[bm$truth == tn, ]
    rows <- sb_gate_add(rows, sprintf("%s %s raw AUPR vs baselines metrics.csv (%s)", m, tn, mp[[m]]), abs(mine$aupr[mine$scope == "raw"] - bs$aupr[bs$method == mp[[m]]]))
    rows <- sb_gate_add(rows, sprintf("%s %s hub_att AUPR vs baselines metrics.csv (%s_hubstrat)", m, tn, mp[[m]]), abs(mine$aupr[mine$scope == "hub_att"] - bs$aupr[bs$method == paste0(mp[[m]], "_hubstrat")]))
    t <- tfm[tfm$readout == m & tfm$truth == tn, ]
    rows <- sb_gate_add(rows, sprintf("%s %s tf_aupr_ratio / tf_auroc_median vs baselines metrics.csv", m, tn), max(abs(t$tf_aupr_ratio - bs$tf_aupr_ratio[bs$method == mp[[m]]]), abs(t$tf_auroc_median - bs$tf_auroc_median[bs$method == mp[[m]]])))
  }
}
pr0 <- utils::read.csv(sprintf("results/baselines_%s/per_regulon.csv", DS), stringsAsFactors = FALSE)
for (m in names(mp)) for (tn in c("SC", "all")) {
  a <- prtab[prtab$readout == m & prtab$truth == tn, ]; b <- pr0[pr0$method == mp[[m]] & pr0$truth == tn, ]
  j <- match(a$regulon, b$regulon); stopifnot(!anyNA(j), nrow(a) == nrow(b))
  d <- max(abs(a$ap - b$ap[j]), na.rm = TRUE); stopifnot(identical(is.na(a$ap), is.na(b$ap[j])))
  rows <- sb_gate_add(rows, sprintf("%s %s per-regulon AP (%d regulons scored) vs per_regulon.csv", m, tn, sum(!is.na(a$ap))), max(d, abs(a$n_pos - b$n_pos[j]), abs(a$base_rate - b$base_rate[j])))
}
if (DS == "m3d") {
  pl <- utils::read.csv("results/explore_stack/m3d/pooled_levels.csv", stringsAsFactors = FALSE)
  mp2 <- c(attention = "K10", CLR = "CLR", a3r = "a3r", a3r2 = "a3r2", a3r3 = "a3r3"); d <- 0
  for (m in names(mp2)) for (tn in c("regulons_SC", "regulons_all")) for (sc in c("raw", "hub_att", "hub_mi")) {
    mine <- pooled$aupr[pooled$readout == m & pooled$truth == sub("regulons_", "", tn) & pooled$scope == sc]
    d <- max(d, abs(mine - pl$aupr[pl$readout == mp2[[m]] & pl$scope == sc & pl$truth == tn])) }
  rows <- sb_gate_add(rows, "attention, CLR, a3r, a3r2, a3r3 pooled AUPR (SC, all; raw, hub_att, hub_mi) vs explore_stack/m3d/pooled_levels.csv", d)
} else {
  ep <- utils::read.csv("docs/clr_logpt_verdict/endpoints_all.csv", stringsAsFactors = FALSE); ep <- ep[ep$dataset == "ecomac" & ep$truth == "regulons_SC", ]
  rows <- sb_gate_add(rows, "a3r SC raw AUPR vs docs/clr_logpt_verdict/endpoints_all.csv (ecomac regulons_SC, S2 and P4)", max(abs(pooled$aupr[pooled$readout == "a3r" & pooled$truth == "SC" & pooled$scope == "raw"] - ep$aupr_new)))
}
sb_gate_finish(rows, file.path(opt$out, "equivalence_gate.csv"), say)

## ---- outputs ----------------------------------------------------------------------------------------------------------------
utils::write.csv(pooled, file.path(opt$out, "pooled.csv"), row.names = FALSE)
utils::write.csv(tfm, file.path(opt$out, "tf_edges.csv"), row.names = FALSE)
utils::write.csv(perreg, file.path(opt$out, "per_regulon_summary.csv"), row.names = FALSE)
utils::write.csv(prtab[prtab$readout %in% c("a3r", "a3r2", "a3r3"), ], file.path(opt$out, "per_regulon_a3r_family.csv"), row.names = FALSE)
say("pooled enrichment (AUPR / base rate):"); print(pooled[pooled$scope != "hub_mi", c("readout", "scope", "truth", "enrichment")], digits = 3, row.names = FALSE)
say("TF edges:"); print(tfm, digits = 4, row.names = FALSE)
say("per-regulon:"); print(perreg, digits = 4, row.names = FALSE)
say("done")
