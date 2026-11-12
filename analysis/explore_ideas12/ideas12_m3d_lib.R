# EXPLORATORY. M3D (E. coli) data and truths for the conditional lab-split stage of idea 1 (copied from analysis/robust_z/robust_z.R, m3d branch).
ei_setup_m3d <- function(say = function(...) cat(..., "\n"), workers = 6L) {
  m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
  raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
  X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
  X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
  desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
  exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description); exper[!grepl("experimenter:", desc$description)] <- NA
  lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])[colnames(X)]
  lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab))); rm(desc, exper)
  storage.mode(X) <- "double"; genes <- rownames(X); G <- length(genes)
  say(sprintf("m3d: %d genes x %d samples (%d labs)", G, ncol(X), length(unique(lab))))
  to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
  rtab <- build_regulon_table("data/RegulonDBExtract", to_bn); operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
  bench <- lapply(c(SC = "SC", all = "all"), function(ev) {
    ck <- if (ev == "SC") c("C", "S") else c("C", "S", "W", "?")
    reg <- make_regulons(rtab, genes, conf_keep = ck)
    list(reg = reg, pairs = comembership_pairs(reg, operon_of, genes)) })
  for (e in names(bench)) say(sprintf("truth %-4s %5d regulons, %8d pairs, base rate %.4f", e, length(bench[[e]]$reg$members), length(bench[[e]]$pairs$label), mean(bench[[e]]$pairs$label)))
  base_rate <- vapply(bench, function(b) mean(b$pairs$label), 0)
  set.seed(EI_SEED)
  draw_half <- function(gin) lapply(seq_len(EI_JACK), function(r) sort(sample(gin, floor(length(gin) / 2))))
  halves <- list(); for (e in names(bench)) halves[[e]] <- draw_half(bench[[e]]$pairs$genes_in)
  list(X = X, genes = genes, G = G, lab = lab, bench = bench, halves = halves, base_rate = base_rate,
       nc = max(1L, min(workers, parallel::detectCores() - 1L)))
}
