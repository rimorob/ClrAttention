# Human ground truths as gene sets for the regulator-agnostic co-membership and
# coherence benchmarks (sourced by analysis/perturbseq.R). Each builder returns
# list(members = list(<set> = integer indices into `genes`), class = named
# character), the structure analysis/regulons.R expects.
#
#   dorothea_sets   DoRothEA TF regulons (OmniPath), by confidence level
#                   (A-B is the analogue of RegulonDB strong/confirmed)
#   collectri_sets  CollecTRI TF regulons (OmniPath)
#   corum_sets      protein complexes with CORUM among their sources
#   response_sets   Perturb-seq knockdown-response modules: for each knockdown
#                   with a significant transcriptional phenotype (energy test),
#                   the genes with |z| > z_min (default 0.5; the 99th percentile of
#                   |z| over all knockdowns and genes is 0.59) (top max_size,
#                   at least min_size), the knocked-down gene excluded. Built
#                   only from perturbed cells, so independent of networks built
#                   on control cells.
# Data: data/human_truth/{dorothea,collectri,complexes}.tsv (OmniPath web
# service), data/perturbseq/K562_essential_normalized_bulk_01.h5ad.

.as_sets <- function(src, tgt, genes, min_size, max_size, class) {
  ok <- !is.na(src) & !is.na(tgt)
  idx <- match(toupper(tgt[ok]), toupper(genes))
  m <- lapply(split(idx, src[ok]), function(v) sort(unique(v[!is.na(v)])))
  m <- m[lengths(m) >= min_size & lengths(m) <= max_size]
  list(members = m, class = stats::setNames(rep(class, length(m)), names(m)))
}
dorothea_sets <- function(genes, levels = c("A", "B"), min_size = 5L, max_size = 500L,
                          path = "data/human_truth/dorothea.tsv") {
  d <- utils::read.delim(path, stringsAsFactors = FALSE)
  lv <- vapply(strsplit(d$dorothea_level, ";"), function(x) any(x %in% levels), NA)
  d <- d[lv, ]
  .as_sets(d$source_genesymbol, d$target_genesymbol, genes, min_size, max_size, "TF")
}
collectri_sets <- function(genes, min_size = 5L, max_size = 500L, path = "data/human_truth/collectri.tsv") {
  d <- utils::read.delim(path, stringsAsFactors = FALSE)
  .as_sets(d$source_genesymbol, d$target_genesymbol, genes, min_size, max_size, "TF")
}
corum_sets <- function(genes, min_size = 3L, max_size = 200L, path = "data/human_truth/complexes.tsv") {
  d <- utils::read.delim(path, stringsAsFactors = FALSE, quote = "")
  d <- d[grepl("CORUM", d$sources), ]
  parts <- strsplit(d$components_genesymbols, "_", fixed = TRUE)
  nm <- make.unique(d$name)
  .as_sets(rep(nm, lengths(parts)), unlist(parts), genes, min_size, max_size, "complex")
}
response_sets <- function(genes, z_min = 0.5, min_size = 10L, max_size = 100L, p_max = 1e-3,
                          path = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad",
                          use = NULL) {   # use: knockdown ids (gene_transcript) allowed as modules
  rd <- function(p) rhdf5::h5read(path, p)
  cat_of <- function(grp, field) {                    # anndata categorical, old or new layout
    x <- rd(paste0(grp, "/", field))
    if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
    cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
    if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L]
  }
  Z <- rd("X")                                        # genes x knockdowns
  vg <- cat_of("/var", "gene_name")
  pert <- cat_of("/obs", "gene_transcript")
  target <- vapply(strsplit(pert, "_"), `[`, "", 2)
  pval <- as.numeric(rd("/obs/energy_test_p_value"))
  keep <- which(pval < p_max & (is.null(use) | pert %in% use))
  gi <- match(toupper(vg), toupper(genes))
  m <- list()
  for (j in keep) {
    z <- abs(Z[, j]); z[toupper(vg) == toupper(target[j])] <- 0
    o <- order(-z); o <- o[z[o] > z_min]; o <- utils::head(o, max_size)
    s <- sort(unique(gi[o][!is.na(gi[o])]))
    if (length(s) >= min_size) m[[pert[j]]] <- s
  }
  list(members = m, class = stats::setNames(rep("response", length(m)), names(m)))
}
