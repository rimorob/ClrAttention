#!/usr/bin/env Rscript
# Spectral-embedding exploratory pass: descriptive comparison tables (docs/plans/spectral_embedding_pass_plan.md section 5). No test, no interval, no pass/fail.
# Reads <in>/<dataset>/metrics_spectral.csv (this pass) and results/scgpt_ext/<dataset>/metrics_ext.csv (comparison columns: tcaCLR2, scgpt_tok = scGPT token cosine, scgpt_attn = scGPT attention,
# and for context aCLR, tcaCLR, clr_ns, scgpt_hid). Jurkat is read only with --datasets Jurkat (after the lead froze the variants); it is not in the default list.
# Criterion (fixed in the plan): for each variant v, dataset d, truth T: ratio_vs_X = mean_u AUPR/base(v; d, T, unit u) / mean_u AUPR/base(X; d, T, u), units u as in score_ext.R.
# Also reported: the per-unit log2 ratio averaged over units and over the three primary truths (the d_i convention of stats_ext.R), as a convenience column.
# Output: <out>/ratio_by_dataset_truth.csv (long: variant, dataset, truth, auprc_variant, auprc_<comparator> and ratio_vs_<comparator>), <out>/ratio_primary_mean_log2.csv, <out>/comparison_means.csv.
a0 <- commandArgs(trailingOnly = TRUE); so <- list(datasets = "hESC,hHep,K562,HepG2,RPE1", `in` = "results/spectral_pass", out = "results/spectral_pass/summary", ext = "results/scgpt_ext")
i <- 1L; while (i <= length(a0)) { k <- sub("^--", "", a0[i]); stopifnot(k %in% names(so)); so[[k]] <- a0[i + 1L]; i <- i + 2L }
DS <- strsplit(so$datasets, ",")[[1]]; PRIM <- c("DoRothEA_AB", "CollecTRI", "Complexes"); COMP <- c(tcaCLR2 = "tcaCLR2", scgpt_tok = "scgpt_tok", scgpt_attn = "scgpt_attn", aCLR = "aCLR", tcaCLR = "tcaCLR", clr_ns = "clr_ns", scgpt_hid = "scgpt_hid")
dir.create(so$out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) utils::read.csv(f, stringsAsFactors = FALSE)
V <- do.call(rbind, lapply(DS, function(d) { x <- rd(file.path(so$`in`, d, "metrics_spectral.csv")); x$unit <- as.character(x$unit); x[, c("dataset", "unit", "truth", "method", "auprc_ratio")] }))
Cm <- do.call(rbind, lapply(DS, function(d) { x <- rd(file.path(so$ext, d, "metrics_ext.csv")); x$unit <- as.character(x$unit); x <- x[x$method %in% COMP, c("dataset", "unit", "truth", "method", "auprc_ratio")]; x }))
cm <- stats::aggregate(auprc_ratio ~ dataset + truth + method, Cm, mean); utils::write.csv(cm, file.path(so$out, "comparison_means.csv"), row.names = FALSE)
vm <- stats::aggregate(auprc_ratio ~ dataset + truth + method, V, mean); names(vm)[names(vm) == "method"] <- "variant"
out <- vm; names(out)[names(out) == "auprc_ratio"] <- "auprc_variant"
nu <- stats::aggregate(unit ~ dataset + truth + method, V, function(u) length(unique(u))); names(nu) <- c("dataset", "truth", "variant", "n_units"); out <- merge(out, nu, by = c("dataset", "truth", "variant"))
for (cn in names(COMP)) { x <- cm[cm$method == COMP[[cn]], c("dataset", "truth", "auprc_ratio")]; names(x)[3] <- paste0("auprc_", cn); out <- merge(out, x, by = c("dataset", "truth"), all.x = TRUE); out[[paste0("ratio_vs_", cn)]] <- out$auprc_variant / out[[paste0("auprc_", cn)]] }
out <- out[order(out$variant, out$dataset, out$truth), ]; utils::write.csv(out, file.path(so$out, "ratio_by_dataset_truth.csv"), row.names = FALSE)
# convenience: mean over units and the three primary truths of log2(variant / comparator), per dataset
res <- list()
for (cn in c("scgpt_tok", "tcaCLR2", "scgpt_attn")) { m <- merge(V, Cm[Cm$method == COMP[[cn]], c("dataset", "unit", "truth", "auprc_ratio")], by = c("dataset", "unit", "truth"), suffixes = c("", "_cmp"))
  m <- m[m$truth %in% PRIM, ]; m$l2 <- log2(m$auprc_ratio / m$auprc_ratio_cmp); a <- stats::aggregate(l2 ~ dataset + truth + method, m, mean); b <- stats::aggregate(l2 ~ dataset + method, a, mean)
  b$comparator <- cn; b$n_truths <- stats::aggregate(truth ~ dataset + method, a, length)$truth; res[[cn]] <- b }
r <- do.call(rbind, res); names(r)[names(r) == "method"] <- "variant"; names(r)[names(r) == "l2"] <- "mean_log2_ratio_over_units_and_primary_truths"; utils::write.csv(r, file.path(so$out, "ratio_primary_mean_log2.csv"), row.names = FALSE)
cat(sprintf("wrote %s: %d rows ratio table, %d rows primary log2 table\n", so$out, nrow(out), nrow(r)))
