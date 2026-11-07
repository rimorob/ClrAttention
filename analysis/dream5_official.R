# DREAM5 Challenge 4: file-level reimplementation of the organisers' scoring scripts in R, and the scorer-gate helpers.
# (analysis/dream5_scorer_gate.R is the command-line gate; the matrix-based Y-E metrics used by the pipeline are in
# analysis/dream5_score_lib.R: d5_official_expected, d5_official_literal, d5_official_p, d5_official_overall.)
# Nothing here touches any score of ours: inputs are the organisers' gold-standard edge files, prediction files (a team directory:
# DREAM5_NetworkInference_<team>_Network{1,3,4}.txt, ranked best first, file order = rank) and the density tables.
# This file is NOT sourced by the truth driver (the driver's code-file guard does not list it).
#
# Lines of the MATLAB sources are cited as file:line (Evaluation scripts/matlab/).

if (!exists("D5")) source("analysis/dream5_score_lib.R")

# load_dream_network.m: tab-separated G<i> G<j> [confidence]; the leading letter is stripped, ids are integers. The confidence is
# never used by the scorer (the file order is the ranking), a missing third column gives NA.
d5_off_read_edges <- function(path) {
  ln <- readLines(path, warn = FALSE); ln <- ln[nzchar(trimws(ln))]
  f <- strsplit(ln, "\t", fixed = TRUE)
  id <- function(k) as.integer(substring(vapply(f, function(x) if (length(x) >= k) x[k] else NA_character_, ""), 2L))
  cc <- suppressWarnings(as.numeric(vapply(f, function(x) if (length(x) >= 3L) x[3] else NA_character_, "")))
  list(i = id(1L), j = id(2L), c = cc)
}

# DREAM5_Challenge4_Evaluation.m (+ edgelist_2_D5C4_gs.m, remove_edges_not_in_gs.m), literally.
#   gold, pred: d5_off_read_edges() lists. Returns P, N, T, L, TPL, AUROC, AUPR (no p-values; see d5_off_score_team).
d5_off_evaluate <- function(gold, pred) {
  keep <- is.na(gold$c) | gold$c != 0                                    # :4-5 (NaN counts as non-zero, as in MATLAB find)
  gi <- gold$i[keep]; gj <- gold$j[keep]
  regs <- sort(unique(gi)); tgts <- sort(unique(c(gi, gj)))              # edgelist_2_D5C4_gs.m:3-4
  key <- function(a, b) as.numeric(a) * 2^24 + b
  posset <- unique(key(gi, gj)); P <- length(posset)
  T <- length(regs) * length(tgts) - length(intersect(regs, tgts)); N <- T - P   # all (regulator, target) cells but i == j (:15-19)
  inbox <- pred$i %in% regs & pred$j %in% tgts & pred$i != pred$j        # remove_edges_not_in_gs.m:10-20; file order kept, no dedupe
  disc <- as.numeric(key(pred$i[inbox], pred$j[inbox]) %in% posset)
  L <- length(disc)
  m <- d5_official_literal(disc, P, T)
  c(P = P, N = N, T = T, L = L, TPL = sum(disc), AUROC = m[["auroc"]], AUPR = m[["aupr"]])
}

# one team: the three networks, p-values by density lookup, go_all.m scores. `files(net)` returns c(gold =, pred =) paths.
d5_off_score_team <- function(files, dens_dir, nets = D5$nets) {
  rows <- lapply(nets, function(n) {
    f <- files(n); r <- d5_off_evaluate(d5_off_read_edges(f[["gold"]]), d5_off_read_edges(f[["pred"]]))
    p <- d5_official_p(n, r[["AUPR"]], r[["AUROC"]], dens_dir)
    data.frame(net = n, t(r), P_AUPR = p[["P_AUPR"]], P_AUROC = p[["P_AUROC"]], row.names = NULL)
  })
  tab <- do.call(rbind, rows)
  list(table = tab, overall = d5_official_overall(tab$P_AUPR, tab$P_AUROC))
}

# ---- scorer gate ---------------------------------------------------------------------------------------------------------------
# Published results table (CSV exported from DREAM5_NetworkInference_MethodScoringResults.xls). The layout of the xls cannot be seen
# here (it is on Synapse only), so columns are recognised by name after normalisation (lower case, non-alphanumerics to "_"):
# a team column (team / method / name / submission) and, per network n in 1,3,4, columns containing "aupr" or "auroc" and the
# network number, with or without a leading "p" / "p_" for the p-value (e.g. AUPR_Net1, P_AUPR_Net1, net1_aupr, "AUPR Network 1").
# Anything else is given explicitly with colmap = c(team = "col", aupr1 = "col", ..., p_aupr1 = "col", ...).
d5_gate_norm <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
d5_gate_read_results <- function(csv, colmap = NULL) {
  tb <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  nm <- d5_gate_norm(names(tb)); find <- function(key) {
    if (!is.null(colmap) && !is.na(colmap[key])) return(match(colmap[[key]], names(tb)))
    metric <- if (grepl("auroc", key)) "auroc" else "aupr"; isp <- grepl("^p_", key); net <- substring(key, nchar(key))
    hit <- which(grepl(metric, nm) & grepl(paste0("(^|[^0-9])", net, "($|[^0-9])"), nm) & (grepl("(^|_)p(_|$)|pval|p_value", nm) == isp) & !grepl("rank|score$", nm))
    if (length(hit) == 1L) hit else NA_integer_
  }
  keys <- c(paste0("aupr", D5$nets), paste0("auroc", D5$nets), paste0("p_aupr", D5$nets), paste0("p_auroc", D5$nets))
  idx <- vapply(keys, find, 0L); names(idx) <- keys
  tcol <- if (!is.null(colmap) && !is.na(colmap["team"])) match(colmap[["team"]], names(tb)) else which(grepl("team|method|name|submission", nm))[1]
  miss <- c(if (is.na(tcol)) "team", names(idx)[is.na(idx)])
  if (length(miss)) stop("cannot find these columns in the results CSV: ", paste(miss, collapse = ", "), "; its (normalised) headers are: ",
                         paste(nm, collapse = ", "), ". Pass colmap = c(team = '<header>', aupr1 = '<header>', ...).")
  num <- function(x) suppressWarnings(as.numeric(x))
  out <- data.frame(team = tb[[tcol]], stringsAsFactors = FALSE)
  for (k in keys) out[[k]] <- num(tb[[idx[[k]]]])
  out
}
d5_gate_team_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
# tolerance rule: AUPR and AUROC absolute `tol`; p-values absolute `tol` and, where the published p is below `tol` (the absolute test
# is vacuous there), also within 5% relative; a published p of exactly 0 needs ours to be 0.
d5_gate_close <- function(ours, pub, tol, kind = c("metric", "p")) {
  kind <- match.arg(kind); ok <- abs(ours - pub) <= tol
  if (kind == "p") ok <- ok & ifelse(pub == 0, ours == 0, ifelse(pub < tol, abs(ours - pub) <= 0.05 * pub, TRUE))
  ok
}
# the gate comparison for a set of teams. `ours`, `port` data.frames: team, net, AUPR, AUROC, P_AUPR, P_AUROC (R path = file-level R
# scorer; port = Python port; `matrix` optional: the pipeline's matrix path). `pub`: d5_gate_read_results() output.
d5_gate_compare <- function(impls, pub, tol = 1e-3) {
  pk <- d5_gate_team_key(pub$team); rows <- list()
  for (nm in names(impls)) { d <- impls[[nm]]
    for (i in seq_len(nrow(d))) {
      k <- match(d5_gate_team_key(d$team[i]), pk); if (is.na(k)) next
      n <- d$net[i]
      for (m in c("AUPR", "AUROC", "P_AUPR", "P_AUROC")) {
        key <- paste0(if (grepl("^P_", m)) "p_" else "", tolower(sub("^P_", "", m)), n); pubv <- pub[[key]][k]
        rows[[length(rows) + 1L]] <- data.frame(impl = nm, team = d$team[i], net = n, quantity = m, ours = d[[m]][i], published = pubv,
          abs_diff = abs(d[[m]][i] - pubv), pass = d5_gate_close(d[[m]][i], pubv, tol, if (grepl("^P_", m)) "p" else "metric"), stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

# The pipeline's matrix path on a public prediction file: the file order becomes the score (position 1 = highest, no ties), the cap is
# the file length, and d5_edge_metrics() is run with the design's Y-E primary switches (gold-derived universe, cap then universe,
# keep_both, drop_tied, expected ties, official readout). Returns NULL if the file repeats a pair (the matrix has one cell per pair;
# the official scorer does not de-duplicate) else c(AUPR, AUROC, N, n_pos, K). `gold_tab`: the gold edge file as read.delim(header = FALSE).
d5_off_matrix_path <- function(pred, gold_tab, nd) {
  edges <- d5_edges_from_ids(gold_tab, nd); u <- d5_gold_universe(edges, nd)
  tfrow <- match(paste0("G", pred$i), nd$tf_ids); gcol <- match(paste0("G", pred$j), nd$genes)
  ok <- !is.na(tfrow) & !is.na(gcol) & nd$tf_idx[tfrow] != gcol
  cells <- cbind(tfrow[ok], gcol[ok]); if (anyDuplicated(cells)) return(NULL)
  nTF <- length(nd$tf_idx); G <- length(nd$genes)
  M <- matrix(-Inf, nTF, G); M[d5_selfmask(nd)] <- NA; M[cells] <- -seq_len(nrow(cells))
  sw <- do.call(d5_edge_switches, D5$ye_primary)
  r <- d5_edge_metrics(M, d5_truth_matrix(edges, nTF, G), nd$tf_idx, sw, cap = nrow(cells), universe = u)
  c(AUPR = r[["aupr_official"]], AUROC = r[["auroc_official"]], N = r[["N"]], n_pos = r[["n_pos"]], K = r[["K"]])
}
