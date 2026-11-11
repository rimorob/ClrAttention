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
# Published results table: the tidy CSV written by analysis/dream5_published_table.py from the organisers' sheet survey4.xls (syn2787254).
# Columns (exact names): team, score_overall, score_aupr, score_auroc, aupr1, aupr3, aupr4, auroc1, auroc3, auroc4. The sheet has NO
# per-network p-values: what can be compared is the per-network AUPR and AUROC and the three aggregate scores, which are the
# organisers' go_all.m scores (aupr_score = -mean(log10 P_AUPR) over networks 1, 3, 4, auroc_score likewise, overall = their mean).
# The sheet stores every value rounded to 3 decimals (a few cells 4), so the rounding error of a stored number is at most 5e-4.
d5_gate_cols <- c("team", "score_overall", "score_aupr", "score_auroc", paste0("aupr", D5$nets), paste0("auroc", D5$nets))
d5_gate_read_results <- function(csv) {
  tb <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  miss <- setdiff(d5_gate_cols, names(tb))
  if (length(miss)) stop("the results CSV lacks the columns ", paste(miss, collapse = ", "), " (it must be the tidy table written by analysis/dream5_published_table.py; its columns are: ", paste(names(tb), collapse = ", "), ")")
  if (anyDuplicated(tb$team)) stop("duplicated team names in the results CSV: ", paste(tb$team[duplicated(tb$team)], collapse = ", "))
  out <- data.frame(team = tb$team, stringsAsFactors = FALSE)
  for (k in setdiff(d5_gate_cols, "team")) {
    v <- suppressWarnings(as.numeric(tb[[k]])); if (anyNA(v)) stop("non-numeric or empty values in column ", k, " of the results CSV")
    out[[k]] <- v
  }
  out
}
d5_gate_team_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
# PASS RULE for one quantity (fixed before the full gate was run): |ours - published| <= tol (default 1e-3) in absolute value, for the
# six per-network metrics and for the three aggregate scores alike. Justification: the table stores values rounded to 3 decimals
# (rounding error <= 5e-4; truncation error < 1e-3, which two cells of the sheet turned out to need), so 1e-3 covers rounding or truncation of the last stored decimal; the "within the printed
# rounding" alternative is not used because trailing zeros are not recoverable from a stored number (8.8 may be 8.800), so the printed
# precision is never coarser than 3 decimals. Two exactly equal values (both Inf) match; a missing value never does.
d5_gate_close <- function(ours, pub, tol) { ok <- ours == pub | abs(ours - pub) <= tol; ok[is.na(ok)] <- FALSE; ok }
# Aggregate scores exactly as the official go_all.m:66-68 computes them: no floor and no cap on p or on -log10 p anywhere in the official
# code (probability.m:2-5 is a plain sum; DREAM5_Challenge4_Evaluation.m:110-116 and go_all.m:66-68 take log10 of it), so p = 0 gives Inf
# and the gate compares it as Inf; a published finite value then fails for that team and the failure is reported, never capped.
d5_gate_aggregates <- function(p_aupr, p_auroc) { o <- d5_official_overall(p_aupr, p_auroc); c(score_overall = o[["overall"]], score_aupr = o[["aupr_score"]], score_auroc = o[["auroc_score"]]) }
# The gate comparison. `impls`: named list of data.frames (team, net, AUPR, AUROC, P_AUPR, P_AUROC): R = file-level R scorer, port =
# Python port, matrix = the pipeline's matrix path (rows of skipped team-networks absent). `pub`: d5_gate_read_results() output. Teams
# are matched by name only. Returns impl, team, net (NA for an aggregate), quantity, ours, published, abs_diff, pass. An implementation
# that lacks a network for a team (matrix path skipped it) gets no aggregate rows for that team.
d5_gate_compare <- function(impls, pub, tol = 1e-3) {
  pk <- d5_gate_team_key(pub$team); rows <- list()
  add <- function(impl, team, net, q, ours, pubv) rows[[length(rows) + 1L]] <<- data.frame(impl = impl, team = team, net = net, quantity = q, ours = ours, published = pubv,
    abs_diff = abs(ours - pubv), pass = d5_gate_close(ours, pubv, tol), stringsAsFactors = FALSE)
  for (nm in names(impls)) { d <- impls[[nm]]
    for (t in unique(d$team)) {
      k <- match(d5_gate_team_key(t), pk); if (is.na(k)) next
      dt <- d[d$team == t, , drop = FALSE]
      for (i in seq_len(nrow(dt))) { n <- dt$net[i]
        add(nm, t, n, "AUPR", dt$AUPR[i], pub[[paste0("aupr", n)]][k]); add(nm, t, n, "AUROC", dt$AUROC[i], pub[[paste0("auroc", n)]][k]) }
      if (setequal(dt$net, D5$nets)) { dt <- dt[match(D5$nets, dt$net), ]; a <- d5_gate_aggregates(dt$P_AUPR, dt$P_AUROC)
        for (q in names(a)) add(nm, t, NA_integer_, q, a[[q]], pub[[q]][k]) }
    }
  }
  if (!length(rows)) return(data.frame(impl = character(), team = character(), net = integer(), quantity = character(), ours = numeric(), published = numeric(), abs_diff = numeric(), pass = logical()))
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
