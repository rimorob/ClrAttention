# scGPT extension (docs/plans/scgpt_extension_plan.md, Revision 1): shared pieces of analysis/scgpt_ext/score_ext.R. Source from the repo root AFTER the set-up of
# score_ext.R (beeline.R set-up and helpers, score_truth, stack_lib.R, robust_z_lib.R, analysis/tcaclr2/edge_depth.R, tsb_lib.R). Nothing here defines a metric, an operator or
# an eigen-solver: the CLR-NS score is beeline.R::ns_clr, the operator gl_make_P (grid_lib.R), the walk read-out sk_a3r (stack_lib.R), the registered chain tsb_walk
# (tsb_lib.R), the edge machinery ad_permute_genes / ad_Psym / ad_eigs / ed_t2 (adaptive_lib.R, edge_depth.R), the directed scorer score_truth (beeline_baselines.R) and the
# complex truth sets / pairs tc_sets, tc_pairs (tc_lib.R).

## ---- dataset registry ------------------------------------------------------------------------------------------------------------
# sc: single-cell control-cell datasets built by analysis/scfm_k562_prep.py (prep directory, case prefix, exports, knockdown-response bulk file)
# beeline: BEELINE cases (case_X of beeline_baselines.R, scGPT extraction of step 3, exports in step4_scoring)
EXT_REG <- list(
  hESC   = list(type = "beeline", ds = "hESC", cases = c("hESC_500", "hESC_1000"), exports = "~/scfm_models/step4_scoring"),
  hHep   = list(type = "beeline", ds = "hHep", cases = c("hHep_500", "hHep_1000"), exports = "~/scfm_models/step4_scoring"),
  K562   = list(type = "sc", h5 = "data/perturbseq/K562_essential_raw_singlecell_01.h5ad", dir = "results/scfm_k562", prefix = "hK562_s", exports = "~/scfm_models/step4_k562",
                bulk = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad"),
  HepG2  = list(type = "sc", h5 = "data/hepg2_perturbseq/GSE264667_hepg2_raw_singlecell_01.h5ad", dir = "results/scfm_hepg2", prefix = "hHepG2_s", exports = "~/scfm_models/step4_hepg2",
                bulk = "results/tcaclr_confirm/pseudobulk/HepG2_rebuilt_bulk.h5ad"),
  RPE1   = list(type = "sc", h5 = "data/perturbseq/rpe1_raw_singlecell_01.h5ad", dir = "results/scfm_rpe1", prefix = "hRPE1_s", exports = "~/scfm_models/step4_rpe1",
                bulk = "data/perturbseq/rpe1_normalized_bulk_01.h5ad"),
  Jurkat = list(type = "sc", h5 = "data/jurkat_perturbseq/GSE264667_jurkat_raw_singlecell_01.h5ad", dir = "results/scfm_jurkat", prefix = "hJurkat_s", exports = "~/scfm_models/step4_jurkat",
                bulk = NA_character_))
EXT_PRIMARY <- c("DoRothEA_AB", "CollecTRI", "Complexes")     # Revision 1: three primary truths, equal weight
EXT_SCGPT <- c(scgpt_attn = "attn_last_mean", scgpt_tok = "tok_cos", scgpt_hid = "hid_cos", scgpt_tokln = "tok_cos_ln")

## ---- small helpers ---------------------------------------------------------------------------------------------------------------
ext_clean <- function(S) { S[!is.finite(S)] <- 0; diag(S) <- 0; S }
ext_rd_mat <- function(dir, src, case, readout) {
  f <- file.path(dir, src, sprintf("%s.%s.f32", case, readout)); g <- readLines(file.path(dir, src, sprintf("%s.genes.txt", case))); G <- length(g)
  m <- matrix(readBin(f, "numeric", n = G * G, size = 4L), G, G); dimnames(m) <- list(g, g); stopifnot(all(is.finite(m))); m
}
ext_rd_genes <- function(dir, src, case) readLines(file.path(dir, src, sprintf("%s.genes.txt", case)))

## ---- tcaCLR and tcaCLR2 on one expression matrix --------------------------------------------------------------------------------------
# Edge rule of docs/plans/tcaclr2_definition.md, CLR-NS extension of analysis/tcaclr2/edge_scoreboard.R (--id beeline), recomputed on THIS matrix: 10 per-gene permutations
# (seeds ED$seed + 7000 + p), operator = gl_make_P(ns_clr(permuted), keff 10, cap 50), edge = mean over permutations of the second eigenvalue of ad_Psym(P),
# t2 = ed_t2(edge). op_diff = max |gl_make_P - rz_lazy_P(rz_operator)| on permutation 1 (the check of edge_scoreboard.R; must be < 1e-12). No truth is read.
ext_edge <- function(X, perms = ED$perms) {
  op_diff <- 0
  m <- do.call(rbind, lapply(seq_len(perms), function(p) {
    Xp <- ad_permute_genes(X, ED$seed + 7000L + p); C <- ns_clr(Xp); P <- gl_make_P(C, ED$keff, ED$cap)
    if (p == 1L) { op_diff <<- max(abs(P - rz_lazy_P(rz_operator(C)))); stopifnot(op_diff < 1e-12) }
    v <- ad_eigs(ad_Psym(P), ED$kvals)$values; c(top = v[1L], edge = v[2L]) }))
  edge <- mean(m[, "edge"])
  data.frame(G = nrow(X), N = ncol(X), edge = edge, edge_min = min(m[, "edge"]), edge_max = max(m[, "edge"]), t2 = ed_t2(edge), tstar = log(ED$eps) / log(edge),
             top_eig_min = min(m[, "top"]), op_diff = op_diff, perms = perms)
}
# tcaCLR (t = 9) and tcaCLR2 (t = t2) from one registered sequential chain on the CLR-NS operator. gate: max abs difference of the t = 9 read-out to the independent
# rz_operator path (sk_layer of stack_lib.R, the path of analysis/scoreboard/sb_beeline.R that reproduced the stored scoreboard); must be < 1e-12.
ext_tcaclr <- function(Cn, t2, gate = TRUE) {
  P <- gl_make_P(Cn, ED$keff, ED$cap)
  o <- tsb_walk(P, c(9L, t2), function(t, Pt) sk_a3r(Pt))
  g <- NA_real_
  if (gate) { ref <- sk_layer(Cn, ts_out = 9L)$outs[["9"]]; g <- max(abs(o[["9"]] - ref)) }
  list(tcaCLR = o[["9"]], tcaCLR2 = o[[as.character(t2)]], gate_rz = g)
}

## ---- complex truth ---------------------------------------------------------------------------------------------------------------------
# Complex Portal human on the unit's gene vector (the G-intersection), by analysis/scgpt_ext/build_complex_truth.R (the builder of analysis/tcaclr_confirm with the one
# documented change for small universes); cached per gene vector. Returns the pair universe of the primary variant `main` (paralog pairs excluded, as arm A).
ext_complex_truth <- function(genes, dir, id, say = function(...) invisible(NULL)) {
  d <- file.path(dir, id); mf <- file.path(d, "cx_members.tsv"); xf <- file.path(d, "cx_excluded_pairs.tsv")
  if (!file.exists(mf)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE); uf <- file.path(d, "universe.txt"); writeLines(genes, uf); pf <- file.path(d, "perturbed_none.txt"); file.create(pf)
    r <- system2("Rscript", c("analysis/scgpt_ext/build_complex_truth.R", uf, pf, d, "cx"), stdout = file.path(d, "build.log"), stderr = file.path(d, "build.log"))
    if (r != 0L || !file.exists(mf)) stop("complex truth build failed, see ", file.path(d, "build.log")) }
  tr <- tc_read_truth(mf, xf); k_all <- tc_excl_keys(tr$excl, "all", genes)
  p <- tc_pairs(tc_sets(tr$mem, "main", genes), genes, k_all); p$sets <- tc_sets(tr$mem, "main", genes)
  say(sprintf("complex truth %s: %d complexes, %d pairs, %d positive (%.4f), %d paralog pairs removed", id, length(p$sets), length(p$label), sum(p$label), mean(p$label), p$n_removed)); p
}
# pooled co-membership AUPR (ap_fast, the registered tie convention) of a symmetric score on the pair universe; ratio = AUPR / base rate
ext_score_complex <- function(Ssym, p) {
  sc <- Ssym[cbind(p$i, p$j)]; dens <- mean(p$label); ap <- ap_fast(sc, p$label)
  c(n_regulators = NA_real_, n_edges = sum(p$label), density = dens, auprc_ratio = ap / dens, epr = NA_real_, auroc = NA_real_)
}
# synthetic complexes for the smoke test (never from a truth file): 40 seed genes, each with its 4 nearest CLR-NS partners plus noise
ext_plant_complexes <- function(genes, sig, seed = 11L, n = 40L) {
  set.seed(seed); G <- length(genes); seeds <- sort(sample.int(G, n)); sets <- lapply(seeds, function(r) { s <- rank(sig[r, ]) / G + stats::rnorm(G, sd = 0.2); s[r] <- -Inf; sort(c(r, order(s, decreasing = TRUE)[1:4])) })
  names(sets) <- paste0("syn", seq_along(sets)); p <- tc_pairs(sets, genes, NULL); p$sets <- sets; p
}

## ---- gates against stored values --------------------------------------------------------------------------------------------------------
# tcaCLR (t = 9) and tcaCLR2 on the full BEELINE case matrix (CLR-NS operator), cell-type ChIP truth: AUPR ratio and EPR vs results/scoreboard/beeline/metrics.csv (method a3r) and
# results/tcaclr2_scoreboard/beeline_metrics.csv (tcaCLR_t9, tcaCLR2), t2 vs the stored edge file. These truths were seen before (BEELINE primary truths); this is a gate, not a new score.
ext_gate_stored <- function(case, say = function(...) invisible(NULL)) {
  ds <- sub("_[0-9]+$", "", case); N <- as.integer(sub("^.*_", "", case)); X <- case_X(ds, N); genes <- rownames(X)
  e <- ext_edge(X); Cn <- ns_clr(X); r <- ext_tcaclr(Cn, e$t2, gate = TRUE)
  stored_edge <- utils::read.delim(sprintf("results/tcaclr2_scoreboard/edge_beeline_%s_%d.tsv", ds, N), stringsAsFactors = FALSE)
  sa <- utils::read.csv("results/scoreboard/beeline/metrics.csv", stringsAsFactors = FALSE); sa <- sa[sa$method == "a3r" & sa$dataset == ds & sa$N == N, ]
  sb <- utils::read.csv("results/tcaclr2_scoreboard/beeline_metrics.csv", stringsAsFactors = FALSE); sb <- sb[sb$dataset == ds & sb$N == N, ]
  gts <- lapply(truths(ds), function(f) utils::read.csv(f, stringsAsFactors = FALSE)); tn <- intersect(c("celltype_ChIP", "nonspecific_ChIP"), names(gts))
  out <- list(data.frame(gate = sprintf("%s: edge vs stored edge file", case), value = abs(e$edge - stored_edge$edge), tol = 1e-9),
              data.frame(gate = sprintf("%s: t2 equal to stored (%d)", case, stored_edge$t2), value = as.numeric(e$t2 != stored_edge$t2), tol = 0),
              data.frame(gate = sprintf("%s: tcaCLR t=9 vs rz_operator path", case), value = r$gate_rz, tol = 1e-12))
  for (m in c("tcaCLR", "tcaCLR2")) for (t in tn) {
    S <- ext_clean(r[[m]]); x <- score_truth(S, genes, NULL, gts[[t]])
    ref <- if (m == "tcaCLR") sa[sa$truth == t, ] else sb[sb$truth == t & sb$method == "tcaCLR2", ]
    if (nrow(ref) == 1L) out[[length(out) + 1L]] <- data.frame(gate = sprintf("%s %s %s: AUPR ratio and EPR vs stored scoreboard", case, m, t),
      value = max(abs(x[["auprc_ratio"]] - ref$auprc_ratio), abs(x[["epr"]] - ref$epr)), tol = 1e-9) }
  g <- do.call(rbind, out); g$pass <- g$value <= g$tol; g
}
