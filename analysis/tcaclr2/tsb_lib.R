# tcaCLR2 on the scoreboard rows (exploratory; rule fixed in docs/plans/tcaclr2_definition.md, not changed here): shared helpers of
# edge_scoreboard.R and score_scoreboard.R. Source from the repo root. Nothing here defines a metric or an operator: the operator, the walk
# read-out (sk_a3r), the pooled / TF-edge / per-regulon scoring and the BEELINE scorer are those of analysis/scoreboard/ and analysis/explore_stack/.
TSB <- list(out = "results/tcaclr2_scoreboard", sb_csv = "docs/whitepapers/figures/scoreboard_data.csv")

tsb_say <- function(logf = NULL) function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); if (!is.null(logf)) cat(m, "\n", file = logf, append = TRUE) }
tsb_args <- function(defaults, args = commandArgs(TRUE)) {
  o <- defaults; i <- 1L
  while (i <= length(args)) { k <- sub("^--", "", args[i]); if (!k %in% names(o)) stop("unknown option --", k)
    v <- args[i + 1L]; o[[k]] <- if (is.character(defaults[[k]])) v else if (is.integer(defaults[[k]])) as.integer(v) else as.numeric(v); i <- i + 2L }
  o
}
tsb_refuse_frozen <- function(path) if (grepl("tcaclr_confirm|docs/whitepapers|scoreboard_data|results/scoreboard/|results/tcaclr2/", normalizePath(path, mustWork = FALSE)))
  stop("refusing to write outside results/tcaclr2_scoreboard: ", path)

# t2 and edge of one edge TSV (edge_depth.R output format)
tsb_edge_read <- function(f) { e <- utils::read.delim(f, stringsAsFactors = FALSE); stopifnot(nrow(e) == 1L, e$t2 >= 1L, is.finite(e$edge), e$keff == 10, e$cap == 50); list(t2 = as.integer(e$t2), edge = e$edge, row = e) }

# registered sequential chain: Pt = P, P %*% P, ... one sparse product per step (as sb_prep, ideas47_lib walk_to, make_scores.R); fn(t, Pt) at each depth of `depths`
tsb_walk <- function(P, depths, fn) {
  depths <- sort(unique(as.integer(depths))); Pt <- as.matrix(P); out <- list()
  if (1L %in% depths) out[["1"]] <- fn(1L, Pt)
  if (max(depths) > 1L) for (s in 2:max(depths)) { Pt <- as.matrix(Pt %*% P); if (s %in% depths) out[[as.character(s)]] <- fn(s, Pt); if (s %% 8L == 0L) invisible(gc(verbose = FALSE)) }
  rm(Pt); invisible(gc(verbose = FALSE)); out
}

# row of the scoreboard_data.csv schema
tsb_row <- function(dataset, truth, metric, scope, method, value, base, source, note = "")
  data.frame(dataset = dataset, truth = truth, metric = metric, scope = scope, method = method, enrichment_over_random = value, base_rate = base,
             source_file = source, status = "computed_now", note = note, stringsAsFactors = FALSE)

# stored tcaCLR (key a3r) value of a cell of scoreboard_data.csv (NA if absent)
tsb_stored <- function(D, dataset, truth, metric, scope, key = "a3r") { x <- D[D$dataset == dataset & D$truth == truth & D$metric == metric & D$scope == scope & D$method == key & D$status != "not_run", ]
  if (nrow(x) == 1L) x$enrichment_over_random else NA_real_ }

# gate: recomputed tcaCLR (t = 9) cells against the stored scoreboard values; returns the table, stops if any compared cell differs by more than tol
tsb_gate <- function(cells, D, tol = 1e-9, say = cat) {
  x <- cells[cells$method == "tcaCLR_t9", ]
  x$stored <- mapply(function(d, t, m, s) tsb_stored(D, d, t, m, s), x$dataset, x$truth, x$metric, x$scope)
  x$abs_diff <- abs(x$enrichment_over_random - x$stored); x$compared <- is.finite(x$stored); x$pass <- !x$compared | x$abs_diff <= tol
  say(sprintf("gate: %d tcaCLR (t = 9) cells recomputed, %d with a stored scoreboard value, max abs diff %.3g (tol %g), %d failing\n", nrow(x), sum(x$compared), suppressWarnings(max(x$abs_diff, na.rm = TRUE)), tol, sum(!x$pass)))
  g <- x[, c("dataset", "truth", "metric", "scope", "enrichment_over_random", "stored", "abs_diff", "compared", "pass")]
  names(g)[5] <- "recomputed"; g
}
tsb_hwm <- function() { l <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0)); x <- grep("^VmHWM:", l, value = TRUE); if (length(x)) as.numeric(gsub("[^0-9]", "", x)) / 1024^2 else NA_real_ }
