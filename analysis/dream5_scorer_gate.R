#!/usr/bin/env Rscript
# DREAM5 SCORER GATE: does our implementation reproduce the organisers' published scores of public teams?
# It runs only when BOTH data inputs exist (they are on Synapse only and are supplied by the PI; nothing is fetched here):
#   --pred_dir     the organisers' public predictions as a FLAT directory: <pred_dir>/<team>/DREAM5_NetworkInference_<team>_Network{1,3,4}.txt
#                  (Network_predictions.zip, syn2787251, extracts into three folders; analysis/dream5_published_table.py flat-pred
#                  builds this layout as symlinks without touching the extracted files)
#   --results_csv  the published results table as the tidy CSV written by `python3 -I analysis/dream5_published_table.py tidy` from the
#                  sheet survey4.xls of DREAM5_NetworkInference_MethodScoringResults.xls (syn2787254): team, score_overall, score_aupr,
#                  score_auroc, aupr1/3/4, auroc1/3/4. The sheet has NO per-network p-values, so p-values are not compared directly;
#                  they enter through the three aggregate scores (go_all.m: -mean log10 p over networks 1, 3, 4; overall = mean of the two).
# plus the local extraction of the scoring archive (gold edges, .mat densities; --scripts_root) and the converted density TSVs.
# For every team present in both (matched by NAME only, never by value), AUPR and AUROC per network and the three aggregate scores are
# computed by (a) our R file-level scorer, (b) the Python port of the official MATLAB scripts and (c) the pipeline's matrix path
# (d5_edge_metrics with the design's Y-E switches; skipped with --matrix 0 and for files that repeat a pair), compared with the published
# values at --tol (default 1e-3 absolute; the rule and its justification are at d5_gate_close), and R vs port vs matrix to 1e-9 (AUPR,
# AUROC, P_AUPR, P_AUROC and -log10 of both p-values). No floor or cap on p is applied anywhere (the official code has none; p = 0 gives Inf).
# PASS needs >= --min_teams (default 3) teams whose every quantity matches the published table for R, the port and, where it ran, the
# matrix path, no disagreement between R, port and matrix beyond 1e-9, and no team with a failing matrix path. Writes <out_dir>/
# gate_report.md and gate_comparison.csv; exit status 0 = PASS, 1 = FAIL, 2 = refused (missing inputs).
#
#   Rscript analysis/dream5_scorer_gate.R --pred_dir <dir> --results_csv <tidy.csv> [--out_dir results/dream5/scorer_gate]
#       [--teams a,b,c] [--min_teams 3] [--tol 1e-3] [--matrix 1] [--python python3]
#       [--scripts_root "data/dream5/truth_extracted/Evaluation scripts"] [--densities data/dream5/truth_extracted_tools/densities]
source("analysis/dream5_official.R")

d5_gate_args <- function(args) {
  o <- list(pred_dir = "", results_csv = "", out_dir = "results/dream5/scorer_gate", teams = "", min_teams = "3", tol = "1e-3",
            matrix = "1", python = "python3", scripts_root = "data/dream5/truth_extracted/Evaluation scripts", densities = D5$densities_dir)
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); if (!grepl("^--", args[i]) || !key %in% names(o)) stop("unknown argument: ", args[i])
    if (i + 1L > length(args)) stop("argument ", args[i], " has no value"); o[[key]] <- args[i + 1L]; i <- i + 2L
  }
  o$min_teams <- as.integer(o$min_teams); o$tol <- as.numeric(o$tol); o$matrix <- o$matrix == "1"
  if (is.na(o$min_teams) || o$min_teams < 1L || is.na(o$tol) || o$tol <= 0) stop("--min_teams must be a positive integer and --tol a positive number")
  o
}
# Returns list(pass, status, report lines, comparison table). `refuse` conditions return status 2 without computing anything.
d5_gate_run <- function(o, say = message) {
  gold_dir <- file.path(o$scripts_root, "INPUT", "gold_standard_edges_only"); mat_dir <- file.path(o$scripts_root, "INPUT", "probability_densities")
  need <- c(`--pred_dir (the organisers' public predictions)` = o$pred_dir, `--results_csv (the published results table)` = o$results_csv)
  miss <- names(need)[!nzchar(need) | !file.exists(need)]
  if (length(miss)) return(list(status = 2L, msg = paste0("REFUSED: required data input missing: ", paste(miss, collapse = "; "),
    ". The scorer gate runs only when the PI has supplied both Synapse files (Network_predictions.zip, syn2787251, extracted as INPUT/predictions; ",
    "DREAM5_NetworkInference_MethodScoringResults.xls, syn2787254, converted to the tidy CSV by analysis/dream5_published_table.py).")))
  if (!dir.exists(o$pred_dir)) return(list(status = 2L, msg = paste0("REFUSED: --pred_dir is not a directory: ", o$pred_dir)))
  loc <- c(gold_dir, mat_dir, file.path(o$densities, sprintf("Network%d_%s.tsv", rep(D5$nets, each = 2), c("AUPR", "AUROC"))))
  if (!all(file.exists(loc))) return(list(status = 2L, msg = paste0("REFUSED: local scoring-archive extraction or converted densities missing: ", paste(loc[!file.exists(loc)], collapse = "; "))))
  pub <- d5_gate_read_results(o$results_csv)
  pf <- function(t, n) file.path(o$pred_dir, t, sprintf("DREAM5_NetworkInference_%s_Network%d.txt", t, n))
  teams <- sort(list.dirs(o$pred_dir, full.names = FALSE, recursive = FALSE)); teams <- teams[nzchar(teams)]
  full <- vapply(teams, function(t) all(file.exists(vapply(D5$nets, function(n) pf(t, n), ""))), NA); incomplete <- teams[!full]; teams <- teams[full]
  if (nzchar(o$teams)) teams <- intersect(teams, strsplit(o$teams, ",", fixed = TRUE)[[1]])
  inpub <- !is.na(match(d5_gate_team_key(teams), d5_gate_team_key(pub$team))); notpub <- teams[!inpub]; teams <- teams[inpub]
  if (!length(teams)) return(list(status = 2L, msg = "REFUSED: no team directory with all three network files matches a team in the results table"))
  say(sprintf("gate: %d teams (R file scorer, Python port, matrix path)", length(teams)))
  # (a) R file-level scorer
  rr <- do.call(rbind, lapply(teams, function(t) { r <- d5_off_score_team(function(n) c(gold = file.path(gold_dir, sprintf("DREAM5_NetworkInference_Edges_Network%d.tsv", n)), pred = pf(t, n)), o$densities)$table
    data.frame(team = t, r[, c("net", "AUPR", "AUROC", "P_AUPR", "P_AUROC", "T", "L")], stringsAsFactors = FALSE) }))
  # (b) Python port
  ptsv <- tempfile(fileext = ".tsv")
  out <- suppressWarnings(system2(o$python, c("-I", "analysis/dream5_official_tools.py", "port-run", shQuote("data/dream5/truth_extracted_tools/d5_official_port.py"), shQuote(gold_dir),
                                              shQuote(mat_dir), shQuote(o$pred_dir), ptsv, teams), stdout = TRUE, stderr = TRUE))
  if (!file.exists(ptsv)) return(list(status = 1L, msg = paste0("Python port failed: ", paste(utils::tail(out, 5), collapse = " | "))))
  pp <- utils::read.delim(ptsv, stringsAsFactors = FALSE)
  # (c) matrix path
  mm <- NULL; mnote <- "not run (--matrix 0)"
  if (o$matrix) {
    nds <- lapply(stats::setNames(D5$nets, D5$nets), d5_read_net)
    mm <- do.call(rbind, lapply(teams, function(t) do.call(rbind, lapply(D5$nets, function(n) {
      gold_tab <- utils::read.delim(file.path(gold_dir, sprintf("DREAM5_NetworkInference_Edges_Network%d.tsv", n)), header = FALSE, stringsAsFactors = FALSE)
      r <- d5_off_matrix_path(d5_off_read_edges(pf(t, n)), gold_tab, nds[[as.character(n)]])
      if (is.null(r)) return(data.frame(team = t, net = n, AUPR = NA_real_, AUROC = NA_real_, P_AUPR = NA_real_, P_AUROC = NA_real_, N = NA_real_))
      p <- d5_official_p(n, r[["AUPR"]], r[["AUROC"]], o$densities)
      data.frame(team = t, net = n, AUPR = r[["AUPR"]], AUROC = r[["AUROC"]], P_AUPR = p[["P_AUPR"]], P_AUROC = p[["P_AUROC"]], N = r[["N"]]) }))))
    mnote <- sprintf("run; %d of %d team-networks skipped (file repeats a pair)", sum(is.na(mm$AUPR)), nrow(mm))
  }
  impls <- list(R = rr, port = pp); if (!is.null(mm)) impls$matrix <- mm[!is.na(mm$AUPR), ]
  cmp <- d5_gate_compare(impls, pub, o$tol)
  # R vs port vs matrix agreement to 1e-9 (AUPR, AUROC, both p-values, and -log10 of both p-values; equal Inf counts as agreeing)
  key <- function(d) paste(d$team, d$net)
  diffs <- function(a, b) { k <- intersect(key(a), key(b)); ia <- match(k, key(a)); ib <- match(k, key(b)); dd <- function(x, y) { z <- abs(x - y); z[x == y] <- 0; z }
    unlist(lapply(c("AUPR", "AUROC", "P_AUPR", "P_AUROC"), function(m) c(dd(a[[m]][ia], b[[m]][ib]), if (grepl("^P_", m)) dd(-log10(a[[m]][ia]), -log10(b[[m]][ib]))))) }
  agree <- function(a, b) max(diffs(a, b), 0, na.rm = TRUE)
  mm_ok <- if (!is.null(mm)) mm[!is.na(mm$AUPR), ] else NULL
  d_rp <- agree(rr, pp); d_rm <- if (!is.null(mm_ok)) agree(rr, mm_ok) else 0
  if (anyNA(unlist(rr[, -1])) || anyNA(unlist(pp[, -1]))) d_rp <- Inf   # a missing value in R or port is a failure, not silence
  ok_t <- vapply(teams, function(t) all(cmp$pass[cmp$team == t]), NA); names(ok_t) <- teams
  pass <- sum(ok_t) >= o$min_teams && d_rp <= 1e-9 && d_rm <= 1e-9
  mx <- function(t, im, q) { v <- cmp$abs_diff[cmp$team == t & cmp$impl == im & cmp$quantity %in% q]; if (all(is.na(v)) || !length(v)) NA_real_ else max(v, na.rm = TRUE) }
  fmt <- function(x) ifelse(is.na(x), "n/a", sprintf("%.3g", x))
  failq <- function(t) { f <- cmp[cmp$team == t & !cmp$pass, ]; if (!nrow(f)) "" else paste(unique(sprintf("%s:%s%s", f$impl, f$quantity, ifelse(is.na(f$net), "", f$net))), collapse = " ") }
  skipped <- if (!is.null(mm)) unique(sprintf("%s net%d", mm$team[is.na(mm$AUPR)], mm$net[is.na(mm$AUPR)])) else character(0)
  mets <- c("AUPR", "AUROC"); scs <- c("score_overall", "score_aupr", "score_auroc")
  rpt <- c(sprintf("# DREAM5 scorer gate: %s", if (pass) "PASS" else "FAIL"), "",
    sprintf("Run: %s. Pass rule: |ours - published| <= %g for the six per-network AUPR/AUROC values and the three aggregate scores of a team (R, Python port, and matrix path where it ran); implementations vs each other 1e-9.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), o$tol),
    sprintf("Inputs: pred_dir `%s`, results `%s`. Matrix path: %s.", o$pred_dir, o$results_csv, mnote),
    sprintf("Teams compared: %d (%s). Teams fully matching the published table: %d; required %d.", length(teams), paste(teams, collapse = ", "), sum(ok_t), o$min_teams),
    sprintf("Max |R - port| = %.3g; max |R - matrix path| = %.3g (AUPR, AUROC, P_AUPR, P_AUROC and their -log10).", d_rp, d_rm),
    if (length(skipped)) sprintf("Matrix path skipped (file repeats a pair; the official scorer does not de-duplicate and the R and port scorers follow it): %s.", paste(skipped, collapse = ", ")),
    if (length(notpub)) sprintf("In pred_dir but not in the results table (not compared): %s.", paste(notpub, collapse = ", ")),
    if (length(incomplete)) sprintf("Directories without all three network files (ignored): %s.", paste(incomplete, collapse = ", ")),
    if (!nzchar(o$teams)) sprintf("Teams in the results table without a prediction directory: %s.", { m <- setdiff(d5_gate_team_key(pub$team), d5_gate_team_key(c(teams, notpub, incomplete))); if (length(m)) paste(pub$team[match(m, d5_gate_team_key(pub$team))], collapse = ", ") else "none" }), "",
    "Columns: max abs difference to the published table over the six per-network AUPR/AUROC values and over the three aggregate scores, for R, the Python port and the matrix path; the last column lists failing quantities (impl:quantity[net]).", "",
    "| team | matches | R: AUPR/AUROC | R: scores | port: AUPR/AUROC | port: scores | matrix: AUPR/AUROC | matrix: scores | failing |", "|---|---|---|---|---|---|---|---|---|",
    vapply(teams, function(t) sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s |", t, if (ok_t[[t]]) "yes" else "NO", fmt(mx(t, "R", mets)), fmt(mx(t, "R", scs)),
           fmt(mx(t, "port", mets)), fmt(mx(t, "port", scs)), fmt(mx(t, "matrix", mets)), fmt(mx(t, "matrix", scs)), failq(t)), ""))
  dir.create(o$out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(rpt, file.path(o$out_dir, "gate_report.md")); utils::write.csv(cmp, file.path(o$out_dir, "gate_comparison.csv"), row.names = FALSE)
  list(status = if (pass) 0L else 1L, msg = rpt[1], pass = pass, teams = teams, ok_teams = teams[ok_t], comparison = cmp, d_rp = d_rp, d_rm = d_rm)
}

if (sys.nframe() == 0L) {
  r <- tryCatch(d5_gate_run(d5_gate_args(commandArgs(trailingOnly = TRUE))), error = function(e) { message("ERROR: ", conditionMessage(e)); list(status = 2L, msg = conditionMessage(e)) })
  message(r$msg); quit(save = "no", status = r$status)
}
