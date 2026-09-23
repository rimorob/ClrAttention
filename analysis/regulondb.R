# RegulonDB gold-standard loading and edge/regulon evaluation helpers for the
# M3D x RegulonDB run. Sourced by run_m3d_regulondb.R; kept free of side
# effects so it can be unit-checked on its own.

# Read a RegulonDB TF -> gene network file into a data.frame with columns
# tf_gene (regulator's *gene* symbol, lower-case first letter as in E. coli
# gene names), target (target gene symbol), confidence ("S"/"W"/"C"/NA).
#
# Supports
#   * RegulonDB >= 12 "NetworkRegulatorGene" TSV: comment lines start with
#     '#'; columns regulatorId, regulatorName, RegulatorGeneName,
#     regulatedId, regulatedName, function, confidenceLevel. Heteromeric
#     regulators list several genes in RegulatorGeneName (e.g. "ihfA;ihfB");
#     each subunit gene is expanded to its own row.
#   * classic network_tf_gene.txt (TF name, gene name, effect, evidence,
#     evidence type); the TF protein name is converted to a gene symbol by
#     lower-casing its first letter (ArcA -> arcA), and heteromers written
#     "IhfA-IhfB" are split.
# The format is detected from the column count/header, and a short summary
# is printed so a wrong guess is visible immediately.
read_regulondb <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  rows <- strsplit(lines, "\t", fixed = TRUE)
  ncol <- stats::median(lengths(rows))
  rows <- rows[lengths(rows) >= min(ncol, 5)]
  hdr <- tolower(rows[[1]])
  has_header <- any(grepl("regulator|regulated|gene|tf", hdr)) &&
    !any(grepl("^[a-z]{3}[A-Z]?$", rows[[1]]))
  if (has_header) rows <- rows[-1]
  tab <- do.call(rbind, lapply(rows, function(r) r[seq_len(max(lengths(rows)))]))
  if (ncol >= 7) {
    # NetworkRegulatorGene: 3 = RegulatorGeneName, 5 = regulatedName, 7 = conf
    if (has_header) {
      ci <- function(pat, default) {
        j <- grep(pat, hdr)
        if (length(j)) j[1] else default
      }
      c_reg <- ci("regulatorgenename", 3L)
      c_tgt <- ci("regulatedname", 5L)
      c_conf <- ci("confidence", 7L)
    } else {
      c_reg <- 3L; c_tgt <- 5L; c_conf <- 7L
    }
    reg <- tab[, c_reg]; tgt <- tab[, c_tgt]; conf <- tab[, c_conf]
    fmt <- "NetworkRegulatorGene"
  } else {
    reg <- tab[, 1]; tgt <- tab[, 2]
    conf <- if (ncol(tab) >= 5) tab[, 5] else NA_character_
    reg <- gsub("-", ";", reg)
    fmt <- "network_tf_gene"
  }
  conf <- toupper(substr(trimws(conf), 1, 1))  # Strong/Weak/Confirmed -> S/W/C
  parts <- strsplit(reg, "[;,/ ]+")
  out <- data.frame(tf_gene = unlist(parts),
                    target = rep(trimws(tgt), lengths(parts)),
                    confidence = rep(conf, lengths(parts)),
                    stringsAsFactors = FALSE)
  if (fmt == "network_tf_gene")
    out$tf_gene <- paste0(tolower(substr(out$tf_gene, 1, 1)),
                          substring(out$tf_gene, 2))
  out$tf_gene <- trimws(out$tf_gene)
  out <- out[nzchar(out$tf_gene) & nzchar(out$target), ]
  out <- unique(out)
  message(sprintf("RegulonDB (%s format): %d TF-gene rows, %d regulator genes, %d targets; confidence: %s",
                  fmt, nrow(out), length(unique(out$tf_gene)),
                  length(unique(out$target)),
                  paste(names(table(out$confidence, useNA = "ifany")),
                        table(out$confidence, useNA = "ifany"),
                        sep = "=", collapse = " ")))
  out
}

# Map gene symbols to row indices of the expression matrix (case-insensitive
# on symbol). sym: character vector of the matrix's gene symbols.
map_symbols <- function(x, sym) {
  match(tolower(x), tolower(sym))
}

# Build the undirected evaluation universe in the style of Faith et al. 2007:
# all pairs (tf, g), tf a regulator gene with >= 1 known target present in the
# data, g any other gene. A pair is positive if either direction is a known
# interaction. Returns list(i, j, label) with i < j (matrix indices), plus the
# TF indices and a per-TF target list.
edge_universe <- function(net, sym, conf_keep = NULL) {
  if (!is.null(conf_keep)) net <- net[net$confidence %in% conf_keep, ]
  ti <- map_symbols(net$tf_gene, sym)
  gi <- map_symbols(net$target, sym)
  ok <- !is.na(ti) & !is.na(gi) & ti != gi
  ti <- ti[ok]; gi <- gi[ok]
  G <- length(sym)
  tfs <- sort(unique(ti))
  pos <- unique(cbind(pmin(ti, gi), pmax(ti, gi)))
  # candidate pairs: every TF with every other gene, undirected, deduped
  cand <- do.call(rbind, lapply(tfs, function(t) {
    g <- setdiff(seq_len(G), t)
    cbind(pmin(t, g), pmax(t, g))
  }))
  cand <- unique(cand)
  key <- function(m) m[, 1] * (G + 1) + m[, 2]
  label <- as.integer(key(cand) %in% key(pos))
  targets <- split(gi, ti)
  list(i = cand[, 1], j = cand[, 2], label = label, tfs = tfs,
       targets = targets, n_pos = sum(label), n_known_rows = sum(ok))
}

# Precision-recall summary of a score vector against 0/1 labels.
# Ties are broken pessimistically (negatives first within a tie) so that a
# constant score cannot look informative.
pr_summary <- function(score, label, prec_levels = c(0.8, 0.6, 0.4)) {
  o <- order(-score, label)
  l <- label[o]
  tp <- cumsum(l); k <- seq_along(l)
  prec <- tp / k; rec <- tp / sum(l)
  aupr <- sum(prec[l == 1]) / sum(l)          # average precision
  # AUROC via Mann-Whitney
  r <- rank(score)
  n1 <- sum(label); n0 <- length(label) - n1
  auroc <- (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  at <- vapply(prec_levels, function(p) {
    ok <- which(prec >= p & k >= 10)
    if (!length(ok)) 0 else tp[max(ok)]
  }, numeric(1))
  names(at) <- paste0("tp_at_prec", prec_levels)
  c(aupr = aupr, auroc = auroc, base_rate = n1 / length(label), at)
}

# Precision/recall of a selected (logical) edge set restricted to the universe.
selected_pr <- function(sel, u) {
  s <- sel[cbind(u$i, u$j)]
  tp <- sum(s & u$label == 1)
  c(n_selected_universe = sum(s), tp = tp,
    precision = if (sum(s)) tp / sum(s) else NA_real_,
    recall = tp / sum(u$label))
}

# Per-TF regulon ranking: for each TF with >= min_targets targets, rank all
# other genes by a TF-row score vector and compute average precision.
# score_rows: |tfs| x G matrix (rows aligned to u$tfs).
regulon_ap <- function(score_rows, u, min_targets = 5L) {
  G <- ncol(score_rows)
  res <- vapply(seq_along(u$tfs), function(k) {
    t <- u$tfs[k]
    tg <- unique(u$targets[[as.character(t)]])
    tg <- setdiff(tg, t)
    if (length(tg) < min_targets) return(NA_real_)
    lab <- integer(G); lab[tg] <- 1L
    keep <- setdiff(seq_len(G), t)
    s <- score_rows[k, keep]; l <- lab[keep]
    o <- order(-s, l)
    l <- l[o]
    sum((cumsum(l) / seq_along(l))[l == 1]) / sum(l)
  }, numeric(1))
  res
}
