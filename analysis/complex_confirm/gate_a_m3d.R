#!/usr/bin/env Rscript
# Gate Ga, M3D part (docs/plans/complex_confirm_design.md, Gates, Ga part 3): the code path that the new driver reuses must reproduce, from scratch
# (no explore_stack cache is read), a3r complexes 0.19209, SC 0.16222 and the P1 difference and SE (0.025659, 0.018169) against
# docs/clr_logpt_verdict/endpoints_all.csv, plus the stored 26 Oct metrics and 100 jackknife replicates. Spent E. coli data only.
# Set-up and rows are those of analysis/explore_stack/stack_m3d.R up to its equivalence gate (copied; the stacked layers are not run).
# Usage (repo root): Rscript analysis/complex_confirm/gate_a_m3d.R [--out results/complex_confirm/gate_a/m3d] [--threads 8] [--workers 8]
args <- commandArgs(TRUE)
opt <- list(out = "results/complex_confirm/gate_a/m3d", threads = 8L, workers = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
source("analysis/explore_idea3/idea3_lib.R"); i3_source_all()      # i3_load, rz_*, regulondb.R, regulons.R, human_truths.R, ap_fast
source("analysis/explore_stack/stack_lib.R")
dir.create(file.path(opt$out, "cache"), recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "run.log"); say <- function(...) i3_say(logf, ...)
QUICK <- opt$quick > 0L; none_op <- stats::setNames(character(0), character(0))
SELDIR <- "docs/clr_logpt_selections"; SCORED <- c("complexes", "regulons_SC", "regulons_all")
say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))

D <- i3_load("m3d", opt$quick); X <- D$X; genes <- D$genes; G <- length(genes); rm(D)
UT <- which(upper.tri(matrix(FALSE, G, G))); nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))
say(sprintf("m3d: %d genes x %d experiments (quick = %d)", G, ncol(X), opt$quick))

## ---- CLR, operator, K10 walk, stratifiers (package code path, as ideas47_lib.R) ---------------------------------------------
fit <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
fit$calibrate(method = "normal", combine = "stouffer")
M <- fit$mi; C <- fit$clr_scores
fit$build_operator(alpha = 0.5, softmax_keff = 10, softmax_cap = 50L); A_pkg <- fit$operator; fit$release(); rm(fit, X)
A <- rz_operator(C); eq_op <- max(abs(A - A_pkg)); rm(A_pkg); say(sprintf("max|A(rz_operator) - package operator| = %.3g", eq_op)); stopifnot(eq_op < 1e-12)
P <- rz_lazy_P(A); q_att <- rz_stationary(P); H_att <- outer(q_att, q_att, "+")
mi_rm <- (rowSums(M) - diag(M)) / (G - 1); H_mi <- outer(mi_rm, mi_rm, "+"); rm(M)
Pt9 <- as.matrix(P); for (s in 2:9) Pt9 <- as.matrix(Pt9 %*% P)
S9 <- (Pt9 + t(Pt9)) / 2; diag(S9) <- 0; S1 <- sk_a3r(Pt9); rm(Pt9, P); invisible(gc(verbose = FALSE))

## ---- truths (analysis/clr_logpt_test.R, verbatim logic for m3d) ------------------------------------------------------------
mk_truth <- function(r, o, genes, role) { if (!length(r$members)) return(NULL)
  p <- comembership_pairs(r, o, genes); if (sum(p$label) == 0) return(NULL); list(reg = r, pairs = p, operon = o, role = role) }
as_reg <- function(sets, cls) list(members = sets, class = stats::setNames(rep(cls, length(sets)), names(sets)))
complexes_with_pos <- function(sets, operon_of, genes, with_operon) sum(vapply(sets, function(m) { o <- operon_of[genes[m]]
  if (!with_operon) length(m) >= 2 && (any(is.na(o)) || length(unique(o)) > 1) else length(m) >= 2 }, TRUE))
want <- strsplit(readLines(file.path(SELDIR, "ecoli_complexes.sha256"), 1L), "[[:space:]]+")[[1]][1]
got <- digest::digest(file = file.path(SELDIR, "ecoli_complexes_members.tsv"), algo = "sha256")
if (!identical(tolower(want), got)) stop("E. coli complex member list SHA-256 mismatch: ", got, " != ", want)
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv"); rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
mem <- utils::read.delim(file.path(SELDIR, "ecoli_complexes_members.tsv"), stringsAsFactors = FALSE); mem <- mem[mem$dataset == "m3d", ]
if (!QUICK) { miss <- setdiff(mem$bnumber, genes); if (length(miss)) stop("E. coli member list has ", length(miss), " genes outside the m3d gene universe") }
sets_of <- function(v) { x <- mem[mem$variant == v, ]; lapply(split(x$bnumber, x$complex), function(g) sort(stats::na.omit(match(g, genes)))) }
cx <- function(v, o, role) { s <- sets_of(v); if (!length(s)) NULL else mk_truth(as_reg(s, "complex"), o, genes, role) }
main <- sets_of("main"); p_main <- if (length(main)) comembership_pairs(as_reg(main, "complex"), operon_of, genes) else NULL
n_pos <- if (is.null(p_main)) 0L else sum(p_main$label); n_cx <- complexes_with_pos(main, operon_of, genes, FALSE); fallback <- !(n_pos >= 200 && n_cx >= 10)
if (!QUICK) stopifnot(!fallback, n_pos == 1150L, n_cx == 54L)
reg_sc <- make_regulons(rtab, genes, conf_keep = c("C", "S")); reg_all <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"))
B <- list()
B$complexes <- cx("main", if (fallback) none_op else operon_of, "primary")
B$regulons_SC <- mk_truth(reg_sc, operon_of, genes, "primary")
B$complexes_operon_incl <- cx("main", none_op, "sensitivity"); B$complexes_size2 <- cx("size2", operon_of, "sensitivity")
B$complexes_nomerge <- cx("nomerge", operon_of, "sensitivity"); B$complexes_noribo <- cx("noribo", operon_of, "sensitivity")
B$regulons_all <- mk_truth(reg_all, operon_of, genes, "descriptive")
ra <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
pick <- function(keep) list(members = ra$members[keep], class = ra$class[keep])
B$TF_regulons <- mk_truth(pick(ra$class == "TF" & lengths(ra$members) <= 500 & lengths(ra$members) >= 5), operon_of, genes, "descriptive")
B$sigmulons <- mk_truth(pick(ra$class == "sigma" & lengths(ra$members) >= 5 & !names(ra$members) %in% c("sigma:sigma70", "sigma:sigma38")), operon_of, genes, "descriptive")
op <- operon_of[names(operon_of) %in% genes]; op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
B$operons <- mk_truth(as_reg(op_sets, "operon"), none_op, genes, "descriptive")
B <- Filter(Negate(is.null), B); rm(ra, op, op_sets)
set.seed(20261022L)                      # clr_logpt_test.R::jack_dataset: halves for all truths in table order
halves_all <- lapply(B, function(b) lapply(seq_len(opt$jack), function(r) sample(b$pairs$genes_in, floor(length(b$pairs$genes_in) / 2))))
bench <- B[SCORED]; halves <- halves_all[SCORED]; rm(B, halves_all)
base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
for (e in names(bench)) say(sprintf("truth %-14s %8d pairs, %7d positive, base rate %.4f", e, length(bench[[e]]$pairs$label), sum(bench[[e]]$pairs$label), base_rate[[e]]))
score_readout <- sk_make_scorer(bench, halves, G, H_att, H_mi, nc, opt$out, say)

## ---- references scored first, then the equivalence gate ---------------------------------------------------------------------
res_clr <- score_readout("CLR", C); res_k10 <- score_readout("K10", S9); res_a3r <- score_readout("a3r", S1)
rows <- list()
L1 <- sk_layer(C, ts_out = 9L); d_e <- max(abs(L1$outs[["9"]] - S1)); d_A <- max(abs(L1$A - A)); rm(L1)
rows <- sk_gate_add(rows, "(e) generic stack N = 1 vs a3r of the package path: max |S1' - S1|", d_e, 1e-12)
rows <- sk_gate_add(rows, "(e) generic stack N = 1 operator vs package operator: max |A' - A|", d_A, 1e-12)
if (!QUICK) {
  mine <- list(clr = res_clr, attention_t9 = res_k10, a3r_t9 = res_a3r)
  cm <- utils::read.csv("results/clr_logpt_test/truth/m3d/metrics.csv", stringsAsFactors = FALSE)
  jr <- utils::read.csv("results/clr_logpt_test/truth/m3d/jack_replicates.csv", stringsAsFactors = FALSE)
  d_au <- 0; d_jk <- 0; d_hu <- 0; d_hj <- 0
  for (m in names(mine)) for (e in SCORED) {
    d_au <- max(d_au, abs(mine[[m]]$raw[[e]]$aupr - cm$aupr[cm$truth == e & cm$method == m]))
    j <- jr[jr$truth == e & jr$method == m, ]; d_jk <- max(d_jk, abs(mine[[m]]$raw[[e]]$jack - j$aupr[order(j$rep)]))
    mh <- paste0(m, "_hubstrat"); d_hu <- max(d_hu, abs(mine[[m]]$hub_att[[e]]$aupr - cm$aupr[cm$truth == e & cm$method == mh]))
    j <- jr[jr$truth == e & jr$method == mh, ]; d_hj <- max(d_hj, abs(mine[[m]]$hub_att[[e]]$jack - j$aupr[order(j$rep)])) }
  rows <- sk_gate_add(rows, "(b) CLR, K10, a3r pooled AUPR (3 truths) vs results/clr_logpt_test/truth/m3d/metrics.csv", d_au)
  rows <- sk_gate_add(rows, "(c) CLR, K10, a3r: 100 stored jackknife replicates x 3 truths vs jack_replicates.csv", d_jk)
  rows <- sk_gate_add(rows, "(b) hub-pi stratified AUPR of CLR, K10, a3r (3 truths) vs metrics.csv (_hubstrat)", d_hu)
  rows <- sk_gate_add(rows, "(c) hub-pi stratified replicates of CLR, K10, a3r vs jack_replicates.csv (_hubstrat)", d_hj)
  ep <- utils::read.csv("docs/clr_logpt_verdict/endpoints_all.csv", stringsAsFactors = FALSE); ep <- ep[ep$dataset == "m3d", ]
  d_ep <- 0
  for (k in seq_len(nrow(ep))) { e <- ep$truth[k]; cmp <- ep$comparator[k]; ref <- if (cmp == "clr") res_clr else res_k10
    d <- res_a3r$raw[[e]]$aupr - ref$raw[[e]]$aupr; se <- stats::sd(res_a3r$raw[[e]]$jack - ref$raw[[e]]$jack)
    d_ep <- max(d_ep, abs(res_a3r$raw[[e]]$aupr - ep$aupr_new[k]), abs(d - ep$d[k]), abs(se - ep$se[k])) }
  rows <- sk_gate_add(rows, "(d) a3r AUPR, a3r minus comparator d and SE vs docs/clr_logpt_verdict/endpoints_all.csv (P1, P3, S1, S3)", d_ep)
  bm <- utils::read.csv("results/baselines_m3d/metrics.csv", stringsAsFactors = FALSE); bm <- bm[bm$method == "attention", ]
  d_a <- max(abs(res_k10$raw$regulons_SC$aupr - bm$aupr[bm$truth == "SC"]), abs(res_k10$raw$regulons_all$aupr - bm$aupr[bm$truth == "all"]))
  rows <- sk_gate_add(rows, "(a) K10 pooled AUPR (SC, all) vs results/baselines_m3d/metrics.csv", d_a)
  rzp <- utils::read.csv("results/robust_z/m3d/pooled_metrics.csv", stringsAsFactors = FALSE); rzp <- rzp[rzp$variant == "V0" & rzp$kind == "att9", ]
  d_r <- 0; for (pr in list(c("regulons_SC", "SC"), c("regulons_all", "all"))) { r <- rzp[rzp$truth == pr[2], ]; br <- base_rate[[pr[1]]]
    d_r <- max(d_r, abs(res_k10$hub_mi[[pr[1]]]$aupr / br - r$enr_hub_mi), abs(res_k10$hub_att[[pr[1]]]$aupr / br - r$enr_hub_att), abs(res_k10$raw[[pr[1]]]$aupr / br - r$enrichment)) }
  rows <- sk_gate_add(rows, "(a) K10 enrichment raw, hub MI, hub pi (SC, all) vs results/robust_z/m3d/pooled_metrics.csv (V0 att9)", d_r)
} else say("QUICK run: gates (a)-(d) skipped (gene subset)")
if (!QUICK) {      # the values the design names
  ep <- utils::read.csv("docs/clr_logpt_verdict/endpoints_all.csv", stringsAsFactors = FALSE); p1 <- ep[ep$code == "P1", ]
  dP1 <- res_a3r$raw$complexes$aupr - res_k10$raw$complexes$aupr; sP1 <- stats::sd(res_a3r$raw$complexes$jack - res_k10$raw$complexes$jack)
  rows <- sk_gate_add(rows, sprintf("(design) a3r M3D complexes AUPR %.6f vs 0.19209", res_a3r$raw$complexes$aupr), abs(round(res_a3r$raw$complexes$aupr, 5) - 0.19209), 1e-9)
  rows <- sk_gate_add(rows, sprintf("(design) a3r M3D SC AUPR %.6f vs 0.16222", res_a3r$raw$regulons_SC$aupr), abs(round(res_a3r$raw$regulons_SC$aupr, 5) - 0.16222), 1e-9)
  rows <- sk_gate_add(rows, sprintf("(design) P1 d %.6f vs 0.025659", dP1), abs(round(dP1, 6) - 0.025659), 1e-9)
  rows <- sk_gate_add(rows, sprintf("(design) P1 SE %.6f vs 0.018169", sP1), abs(round(sP1, 6) - 0.018169), 1e-9)
  rows <- sk_gate_add(rows, "(design) P1 d and SE vs endpoints_all.csv row P1 (full precision)", max(abs(dP1 - p1$d), abs(sP1 - p1$se)))
}
g <- do.call(rbind, rows); g$part <- "M3D"
utils::write.csv(g, file.path(opt$out, "equivalence_gate_m3d.csv"), row.names = FALSE); print(g, row.names = FALSE)
if (!all(g$pass)) { say("M3D GATE FAILED"); quit(save = "no", status = 1L) }
say("M3D GATE PASSED")
