# Truth loaders and matrix loaders for stage 2 (descriptive curves on truths that are already used up). They reuse the registered builders:
#   Arm A / Arm B / net4 : tc_read_truth + tc_build_truths(...)$main with the frozen members / excluded files (SHA-256 guard as score_confirm.R),
#                          as in analysis/tcaclr_confirm/jobs/7380*.sh, 7381*.sh and 7101*.sh (pooled AUPR over all pairs among the truth genes)
#   K562 CORUM           : corum_sets(genes) (analysis/human_truths.R) -> comembership_pairs, as analysis/baselines.R
#   M3D TF regulons      : RegulonDB regulons, strong/confirmed (SC) and all evidence, same-operon pairs removed, as analysis/baselines.R
# Source after grid_lib.R, from the repo root.
suppressPackageStartupMessages({ source("analysis/tcaclr_confirm/tc_lib.R"); source("analysis/human_truths.R") })
GT <- list(
  k562_corum     = list(matrix = "results/baselines_k562ess_truthgenes_genie3/expression.tsv"),
  net4_cp_yeast  = list(matrix = "net4"),
  armA_cp_human  = list(matrix = "results/tcaclr_confirm/armA/expression.tsv", dir = "docs/tcaclr_confirm_selections", sha = "armA.sha256",
                        members = "armA_complex_portal_members.tsv", excluded = "armA_complex_portal_excluded_pairs.tsv"),
  armB_cyc2008   = list(matrix = "results/tcaclr_confirm/armB/expression.tsv", dir = "docs/tcaclr_confirm_selections", sha = "armB.sha256",
                        members = "armB_cyc2008_members.tsv", excluded = "armB_cyc2008_excluded_pairs.tsv"),
  m3d_tf_regulons = list(matrix = "results/baselines_m3d/expression.tsv"))
# the expression matrix of a truth (genes x profiles, rownames = the gene names of the truth); "net4" is read from the DREAM5 files as validate_net4_prep.R
gl_truth_matrix <- function(name, matrix = NULL) {
  m <- if (is.null(matrix) || !nzchar(matrix)) GT[[name]]$matrix else matrix
  if (identical(m, "net4")) { source("analysis/dream5_lib.R"); nd <- d5_read_net(4L); X <- nd$X; rownames(X) <- nd$gene_names; return(X) }
  gl_read_matrix(m)
}
# list of truths for the genes of X: each list(scope, i, j, label, ix, ixT, n_pairs, n_pos, genes_in); ix = i + (j - 1) G, ixT = j + (i - 1) G
gl_truths <- function(name, genes) {
  G <- length(genes); spec <- GT[[name]]; if (is.null(spec)) stop("unknown truth: ", name)
  none <- stats::setNames(character(0), character(0)); out <- list()
  add <- function(scope, p) out[[scope]] <<- list(scope = scope, i = p$i, j = p$j, label = as.integer(p$label), ix = as.integer(p$i + (as.numeric(p$j) - 1) * G),
    ixT = as.integer(p$j + (as.numeric(p$i) - 1) * G), n_pairs = length(p$label), n_pos = sum(p$label), genes_in = p$genes_in)
  if (name == "k562_corum") add("CORUM", comembership_pairs(corum_sets(genes), none, genes))
  else if (name == "net4_cp_yeast") {
    d <- "docs/complex_confirm_selections"; tc_hash_guard(file.path(d, "yeast_complexes.sha256"), c("yeast_complexes_members.tsv", "yeast_complexes_excluded_pairs.tsv"))
    tr <- tc_read_truth(file.path(d, "yeast_complexes_members.tsv"), file.path(d, "yeast_complexes_excluded_pairs.tsv"))
    add("main", tc_build_truths(tr, genes, "net4")$main)
  } else if (name %in% c("armA_cp_human", "armB_cyc2008")) {
    tc_hash_guard(file.path(spec$dir, spec$sha), c(spec$members, spec$excluded))
    tr <- tc_read_truth(file.path(spec$dir, spec$members), file.path(spec$dir, spec$excluded))
    add("main", tc_build_truths(tr, genes, "plan")$main)
  } else if (name == "m3d_tf_regulons") {
    to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv"); rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
    for (ev in c("SC", "all")) { ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
      add(ev, comembership_pairs(make_regulons(rtab, genes, conf_keep = ck), operon_of, genes)) }
  } else stop("unknown truth: ", name)
  out
}
