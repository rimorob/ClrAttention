#!/usr/bin/env Rscript
# DREAM5 SCORER GATE: does our implementation reproduce the organisers' published scores of public teams?
# It runs only when BOTH data inputs exist (they are on Synapse only and are supplied by the PI; nothing is fetched here):
#   --pred_dir     the organisers' public predictions, laid out as INPUT/predictions: <pred_dir>/<team>/DREAM5_NetworkInference_<team>_Network{1,3,4}.txt
#                  (from Network_predictions.zip, syn2787251)
#   --results_csv  the published results table, exported to CSV from DREAM5_NetworkInference_MethodScoringResults.xls (syn2787254)
# plus the local extraction of the scoring archive (gold edges, .mat densities; --scripts_root) and the converted density TSVs.
# For every team present in both, AUPR, AUROC, P_AUPR and P_AUROC are computed by (a) our R file-level scorer, (b) the Python port of the
# official MATLAB scripts and (c) the pipeline's matrix path (d5_edge_metrics with the design's Y-E switches; skipped with --matrix 0),
# compared with the published values at --tol (default 1e-3; see d5_gate_close for the p-value rule), and R vs port vs matrix to 1e-9.
# PASS needs >= --min_teams (default 3) teams whose every quantity matches the published table for both R and the port, no
# disagreement between R, port and matrix beyond 1e-9 in any team, and no team with a matrix path that fails. Writes <out_dir>/
# gate_report.md and gate_comparison.csv; exit status 0 = PASS, 1 = FAIL, 2 = refused (missing inputs).
#
#   Rscript analysis/dream5_scorer_gate.R --pred_dir <dir> --results_csv <file.csv> [--out_dir results/dream5/scorer_gate]
#       [--teams a,b,c] [--min_teams 3] [--tol 1e-3] [--colmap "team=Team;aupr1=AUPR_Net1;..."] [--matrix 1] [--python python3]
#       [--scripts_root "data/dream5/truth_extracted/Evaluation scripts"] [--densities data/dream5/truth_extracted_tools/densities]
source("analysis/dream5_official.R")

d5_gate_args <- function(args) {
  o <- list(pred_dir = "", results_csv = "", out_dir = "results/dream5/scorer_gate", teams = "", min_teams = "3", tol = "1e-3", colmap = "",
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
    "DREAM5_NetworkInference_MethodScoringResults.xls, syn2787254, exported to CSV).")))
  if (!dir.exists(o$pred_dir)) return(list(status = 2L, msg = paste0("REFUSED: --pred_dir is not a directory: ", o$pred_dir)))
  loc <- c(gold_dir, mat_dir, file.path(o$densities, sprintf("Network%d_%s.tsv", rep(D5$nets, each = 2), c("AUPR", "AUROC"))))
  if (!all(file.exists(loc))) return(list(status = 2L, msg = paste0("REFUSED: local scoring-archive extraction or converted densities missing: ", paste(loc[!file.exists(loc)], collapse = "; "))))
  colmap <- if (nzchar(o$colmap)) { kv <- strsplit(strsplit(o$colmap, ";", fixed = TRUE)[[1]], "=", fixed = TRUE); stats::setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1)) } else NULL
  pub <- d5_gate_read_results(o$results_csv, colmap)
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
  # R vs port vs matrix agreement to 1e-9
  q <- c("AUPR", "AUROC", "P_AUPR", "P_AUROC"); key <- function(d) paste(d$team, d$net)
  agree <- function(a, b) { k <- intersect(key(a), key(b)); max(abs(unlist(lapply(q, function(m) a[[m]][match(k, key(a))] - b[[m]][match(k, key(b))]))), 0, na.rm = TRUE) }
  d_rp <- agree(rr, pp); d_rm <- if (!is.null(mm)) agree(rr, mm[!is.na(mm$AUPR), ]) else 0
  ok_t <- vapply(teams, function(t) all(cmp$pass[cmp$team == t & cmp$impl %in% c("R", "port")]) && all(cmp$pass[cmp$team == t & cmp$impl == "matrix"]), NA)
  pass <- sum(ok_t) >= o$min_teams && d_rp <= 1e-9 && d_rm <= 1e-9
  rpt <- c(sprintf("# DREAM5 scorer gate: %s", if (pass) "PASS" else "FAIL"), "",
    sprintf("Run: %s. Tolerance vs published %g (p-values: see d5_gate_close); implementations vs each other 1e-9.", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), o$tol),
    sprintf("Inputs: pred_dir `%s`, results `%s`. Matrix path: %s.", o$pred_dir, o$results_csv, mnote),
    sprintf("Teams compared: %d (%s). Teams fully matching the published table (R and Python port, and matrix path where run): %d; required %d.", length(teams), paste(teams, collapse = ", "), sum(ok_t), o$min_teams),
    sprintf("Max |R - port| = %.3g; max |R - matrix path| = %.3g.", d_rp, d_rm),
    if (length(notpub)) sprintf("In pred_dir but not in the results table (not compared): %s.", paste(notpub, collapse = ", ")),
    if (length(incomplete)) sprintf("Directories without all three network files (ignored): %s.", paste(incomplete, collapse = ", ")), "",
    "| team | matches published | max abs diff (R) | max abs diff (port) |", "|---|---|---|---|",
    vapply(teams, function(t) sprintf("| %s | %s | %.3g | %.3g |", t, if (ok_t[[t]]) "yes" else "NO",
           max(cmp$abs_diff[cmp$team == t & cmp$impl == "R"], 0, na.rm = TRUE), max(cmp$abs_diff[cmp$team == t & cmp$impl == "port"], 0, na.rm = TRUE)), ""))
  dir.create(o$out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(rpt, file.path(o$out_dir, "gate_report.md")); utils::write.csv(cmp, file.path(o$out_dir, "gate_comparison.csv"), row.names = FALSE)
  list(status = if (pass) 0L else 1L, msg = rpt[1], pass = pass, teams = teams, ok_teams = teams[ok_t], comparison = cmp, d_rp = d_rp, d_rm = d_rm)
}

if (sys.nframe() == 0L) {
  r <- tryCatch(d5_gate_run(d5_gate_args(commandArgs(trailingOnly = TRUE))), error = function(e) { message("ERROR: ", conditionMessage(e)); list(status = 2L, msg = conditionMessage(e)) })
  message(r$msg); quit(save = "no", status = r$status)
}
