# Regulator-agnostic regulons from a RegulonDB (v12+, release 14.x) extract,
# and co-membership evaluation of gene-gene similarity matrices.
#
# Rationale (user, 2026-09-23): regulators are almost always latent -- a TF's
# mRNA need not track its activity, ppGpp is not a gene, sRNAs are not on the
# genes-only array -- so the primary benchmark asks whether genes that share a
# regulator (of ANY kind: TF, sRNA, small molecule, other protein, sigma
# factor) are recovered as related, without ever using the regulator as a
# node. Same-operon pairs are excluded because co-transcription makes them
# trivially co-expressed.
#
# Inputs (data/RegulonDBExtract/, as downloaded from RegulonDB "Datasets"):
#   RISet.tsv            all regulatory interactions; `type` is one of
#                        TF-/sRNA-/compound-/regulator- x promoter/TU/gene;
#                        targetTuOrGene is "<RegulonDB id>:<name>" where the id
#                        is a TU (RDBECOLITUC...) or a gene (RDBECOLIGNC...)
#   TUSet.tsv            TU id -> tuGenes (";"-separated)
#   OperonSet.tsv        operonGenes ("|"-separated)
#   NetworkSigmaGene.tsv sigma factor -> gene (sigmulons)
# Gene names are mapped to Blattner numbers with analysis/ecoli_k12_genes.tsv,
# derived from NCBI gene_info for E. coli K-12 MG1655 (current symbol +
# synonyms), because RegulonDB 2026 uses current names while the 2008 M3D
# compendium uses the names of its time (e.g. acrZ was ybhT, cra was fruR).

# Read a RegulonDB TSV: '#' comment lines, header "1)id\t2)type...".
read_rdb_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[!startsWith(lines, "#") & nzchar(trimws(lines))]
  rows <- strsplit(sub("[ \t]+$", "", lines), "\t", fixed = TRUE)
  hdr <- trimws(sub("^[0-9]+\\)", "", rows[[1]]))
  n <- length(hdr)
  body <- lapply(rows[-1], function(r) { length(r) <- n; r })
  df <- as.data.frame(do.call(rbind, body), stringsAsFactors = FALSE)
  names(df) <- hdr
  df[] <- lapply(df, function(x) { x <- trimws(x); x[is.na(x)] <- ""; x })
  df
}

# Name -> b-number mapper. Resolution order: exact current symbol, exact
# b-number, unique synonym, then the same three case-insensitively.
# Ambiguous synonyms (claimed by several genes) are not used.
make_bnumber_mapper <- function(gene_tab_path) {
  gt <- utils::read.delim(gene_tab_path, stringsAsFactors = FALSE,
                          quote = "", comment.char = "")
  syn <- strsplit(gt$synonyms, "|", fixed = TRUE)
  syn_long <- data.frame(name = unlist(syn),
                         bnumber = rep(gt$bnumber, lengths(syn)),
                         stringsAsFactors = FALSE)
  syn_long <- syn_long[syn_long$name != "-" & !grepl("^ECK", syn_long$name), ]
  amb <- names(which(tapply(syn_long$bnumber, syn_long$name,
                            function(b) length(unique(b))) > 1))
  syn_long <- syn_long[!syn_long$name %in% amb, ]
  lut_sym <- stats::setNames(gt$bnumber, gt$symbol)
  lut_bn <- stats::setNames(gt$bnumber, gt$bnumber)
  lut_syn <- stats::setNames(syn_long$bnumber, syn_long$name)
  lower <- function(l) stats::setNames(l, tolower(names(l)))
  luts <- list(lut_sym, lut_bn, lut_syn, lower(lut_sym), lower(lut_bn),
               lower(lut_syn))
  function(x) {
    out <- rep(NA_character_, length(x))
    for (k in seq_along(luts)) {
      key <- if (k <= 3) x else tolower(x)
      hit <- is.na(out) & key %in% names(luts[[k]])
      out[hit] <- luts[[k]][key[hit]]
    }
    unname(out)
  }
}

# Build the long regulon table: one row per (regulator, target b-number).
# Columns: regulator, reg_class (TF, sRNA, compound, protein, sigma),
# bnumber, target_name, confidence (C/S/W/?), effect, via (gene/TU/promoter).
build_regulon_table <- function(rdb_dir, to_bnumber) {
  ri <- read_rdb_tsv(file.path(rdb_dir, "RISet.tsv"))
  tu <- read_rdb_tsv(file.path(rdb_dir, "TUSet.tsv"))
  tu_genes <- stats::setNames(strsplit(tu$tuGenes, ";", fixed = TRUE), tu$id)
  tid <- sub(":.*$", "", ri$targetTuOrGene)
  tnm <- sub("^[^:]*:", "", ri$targetTuOrGene)
  is_tu <- startsWith(tid, "RDBECOLITUC")
  if (any(is_tu & !tid %in% names(tu_genes)))
    warning(sum(is_tu & !tid %in% names(tu_genes)),
            " RISet rows point at TU ids absent from TUSet", call. = FALSE)
  targets <- as.list(tnm)
  targets[is_tu] <- lapply(tid[is_tu], function(i) tu_genes[[i]] %||% character())
  cls_raw <- sub("-.*$", "", ri$type)
  cls <- c(TF = "TF", sRNA = "sRNA", compound = "compound",
           regulator = "protein")[cls_raw]
  cls[is.na(cls)] <- cls_raw[is.na(cls)]
  via <- sub("^.*-", "", ri$type)
  n_t <- lengths(targets)
  tab <- data.frame(regulator = rep(ri$regulatorName, n_t),
                    reg_class = rep(unname(cls), n_t),
                    target_name = trimws(unlist(targets)),
                    confidence = rep(ri$confidenceLevel, n_t),
                    effect = rep(ri$riFunction, n_t),
                    via = rep(via, n_t), stringsAsFactors = FALSE)
  sg <- read_rdb_tsv(file.path(rdb_dir, "NetworkSigmaGene.tsv"))
  tab <- rbind(tab, data.frame(regulator = sg$sigmaName, reg_class = "sigma",
                               target_name = sg$regulatedGeneName,
                               confidence = sg$confidenceLevel,
                               effect = sg$`function`, via = "gene",
                               stringsAsFactors = FALSE))
  tab <- tab[nzchar(tab$target_name), ]
  tab$bnumber <- to_bnumber(tab$target_name)
  tab
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Operon id per b-number (genes absent from OperonSet get their own id).
build_operon_map <- function(rdb_dir, to_bnumber) {
  op <- read_rdb_tsv(file.path(rdb_dir, "OperonSet.tsv"))
  g <- strsplit(op$operonGenes, "|", fixed = TRUE)
  df <- data.frame(operon = rep(op$operonId, lengths(g)),
                   bnumber = to_bnumber(trimws(unlist(g))),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$bnumber) & !duplicated(df$bnumber), ]
  stats::setNames(df$operon, df$bnumber)
}

# Regulons restricted to the compendium's genes, filtered by confidence and
# size. genes: b-numbers of the expression matrix rows (in order).
# Returns list(members = named list of integer row indices,
#              class = named character, size = named integer).
make_regulons <- function(tab, genes, conf_keep = c("C", "S", "W", "?"),
                          min_size = 5L, max_size = 500L) {
  t <- tab[tab$confidence %in% conf_keep & !is.na(tab$bnumber), ]
  t$row <- match(t$bnumber, genes)
  t <- t[!is.na(t$row), ]
  key <- paste(t$reg_class, t$regulator, sep = ":")
  mem <- lapply(split(t$row, key), function(v) sort(unique(v)))
  cls <- vapply(split(t$reg_class, key), `[`, "", 1L)
  sz <- lengths(mem)
  keep <- sz >= min_size & sz <= max_size
  list(members = mem[keep], class = cls[keep], size = sz[keep],
       dropped = data.frame(regulon = names(sz)[!keep], size = sz[!keep],
                            stringsAsFactors = FALSE))
}

# Pair universe for co-membership: all unordered pairs of genes that belong
# to at least one retained regulon, minus same-operon pairs. label = 1 if the
# pair shares any regulon; per-class labels mark sharing a regulon of that
# class. Returns list(i, j, label, class_labels (matrix), genes_in).
comembership_pairs <- function(reg, operon_of, genes) {
  gin <- sort(unique(unlist(reg$members)))
  n <- length(gin)
  pos_in <- integer(length(genes)); pos_in[gin] <- seq_len(n)
  shared <- matrix(FALSE, n, n)
  classes <- sort(unique(reg$class))
  by_class <- lapply(classes, function(c) matrix(FALSE, n, n))
  names(by_class) <- classes
  for (k in seq_along(reg$members)) {
    m <- pos_in[reg$members[[k]]]
    shared[m, m] <- TRUE
    by_class[[reg$class[k]]][m, m] <- TRUE
  }
  ut <- which(upper.tri(shared), arr.ind = TRUE)
  i <- gin[ut[, 1]]; j <- gin[ut[, 2]]
  op_i <- operon_of[genes[i]]; op_j <- operon_of[genes[j]]
  same_op <- !is.na(op_i) & !is.na(op_j) & op_i == op_j
  keep <- !same_op
  cl <- vapply(by_class, function(M) M[ut][keep], logical(sum(keep)))
  list(i = i[keep], j = j[keep], label = as.integer(shared[ut][keep]),
       class_labels = cl, n_same_operon_removed = sum(same_op),
       genes_in = gin)
}

# Co-membership AUPR/AUROC of a symmetric similarity matrix S (G x G):
# overall, and per regulator class (positives: pairs sharing a regulon of
# that class; negatives: pairs sharing no regulon at all).
comembership_eval <- function(S, pairs) {
  s <- S[cbind(pairs$i, pairs$j)]
  out <- list(data.frame(stratum = "all", t(pr_summary(s, pairs$label)),
                         n_pos = sum(pairs$label), check.names = FALSE))
  neg <- pairs$label == 0L
  for (c in colnames(pairs$class_labels)) {
    pc <- pairs$class_labels[, c]
    sel <- pc | neg
    if (sum(pc) < 20) next
    out[[length(out) + 1L]] <- data.frame(
      stratum = c, t(pr_summary(s[sel], as.integer(pc[sel]))),
      n_pos = sum(pc), check.names = FALSE)
  }
  do.call(rbind, out)
}

# Gene x gene logical matrix: TRUE if the two genes share a retained regulon.
comember_matrix <- function(reg, G) {
  inc <- matrix(0, G, length(reg$members))
  for (k in seq_along(reg$members)) inc[reg$members[[k]], k] <- 1
  co <- tcrossprod(inc) > 0
  diag(co) <- TRUE
  co
}

# Per-regulon coherence: AUROC of within-regulon pairs versus pairs linking a
# member to a gene that shares NO retained regulon with that member
# (same-operon pairs excluded on both sides). 0.5 = no coherence. Average
# ranks, so tied scores are treated fairly.
regulon_coherence <- function(S, reg, co, operon_of, genes) {
  op <- operon_of[genes]
  vapply(seq_along(reg$members), function(k) {
    m <- reg$members[[k]]
    if (length(m) < 3) return(NA_real_)
    same_m <- outer(op[m], op[m], "==")
    w <- S[m, m][upper.tri(same_m) & (is.na(same_m) | !same_m)]
    ok_b <- !co[m, , drop = FALSE]
    same_b <- outer(op[m], op, "==")
    ok_b <- ok_b & (is.na(same_b) | !same_b)
    b <- S[m, , drop = FALSE][ok_b]
    if (length(w) < 3 || length(b) < 10) return(NA_real_)
    r <- rank(c(w, b))
    n1 <- length(w); n0 <- length(b)
    (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (as.numeric(n1) * n0)
  }, numeric(1))
}
