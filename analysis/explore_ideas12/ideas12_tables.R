# EXPLORATORY. Tables from scored readouts (point estimates and intervals only; no p-values, no winner selection).
# R: named list of scored readouts (each list(pooled, hub_att, hub_mi) of per-truth list(aupr, jack)); must contain "CLR" and "K10".
ei_tables <- function(R, br, outdir, tag, readouts = setdiff(names(R), c("CLR", "K10")), say = function(...) cat(..., "\n")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE); truths <- names(br)
  scopes <- c(pooled = "raw", hub_att = "hub-stratified (pi of V0 attention)", hub_mi = "hub-stratified (MI row mean)")
  pm <- do.call(rbind, lapply(names(R), function(k) do.call(rbind, lapply(names(scopes), function(sc) do.call(rbind, lapply(truths, function(e) {
    x <- R[[k]][[sc]]; if (is.null(x)) return(NULL); data.frame(readout = k, scope = sc, truth = e, base_rate = br[[e]], aupr = x[[e]]$aupr, enrichment = x[[e]]$aupr / br[[e]]) }))))))
  utils::write.csv(pm, file.path(outdir, paste0(tag, "_pooled_metrics.csv")), row.names = FALSE)
  ct <- do.call(rbind, lapply(readouts, function(k) do.call(rbind, lapply(names(scopes), function(sc) do.call(rbind, lapply(c("CLR", "K10"), function(v) do.call(rbind, lapply(truths, function(e) {
    a <- R[[k]][[sc]][[e]]; b <- R[[v]][[sc]][[e]]; if (is.null(a) || is.null(a$jack)) return(NULL)
    d <- a$aupr - b$aupr; dj <- a$jack - b$jack; s <- stats::sd(dj)
    data.frame(readout = k, scope = sc, truth = e, versus = v, d_aupr = d, lo_aupr = d - 1.96 * s, hi_aupr = d + 1.96 * s,
               d_enr = d / br[[e]], lo_enr = (d - 1.96 * s) / br[[e]], hi_enr = (d + 1.96 * s) / br[[e]], n_rep = length(dj)) }))))))))
  utils::write.csv(ct, file.path(outdir, paste0(tag, "_contrasts.csv")), row.names = FALSE)
  mt <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"), apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  md <- character(0); allr <- c("CLR", "K10", readouts)
  for (sc in names(scopes)) { x <- pm[pm$scope == sc, ]
    w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = sprintf("%s (base %.4f)", e, br[[e]]), t(stats::setNames(sprintf("%.2f", y$enrichment), y$readout))[, allr, drop = FALSE], check.names = FALSE) }))
    md <- c(md, sprintf("### Pooled enrichment (AUPR / base rate), %s", scopes[[sc]]), "", mt(w), "")
    w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = e, t(stats::setNames(sprintf("%.5f", y$aupr), y$readout))[, allr, drop = FALSE], check.names = FALSE) }))
    md <- c(md, sprintf("### Pooled AUPR, %s", scopes[[sc]]), "", mt(w), "") }
  for (sc in names(scopes)) for (v in c("K10", "CLR")) { x <- ct[ct$scope == sc & ct$versus == v, ]; if (!nrow(x)) next
    x$cell <- sprintf("%+.4f [%+.4f, %+.4f]", x$d_aupr, x$lo_aupr, x$hi_aupr)
    w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = e, t(stats::setNames(y$cell, y$readout))[, readouts, drop = FALSE], check.names = FALSE) }))
    md <- c(md, sprintf("### Difference in AUPR vs %s [95%% paired delete-half gene jackknife], %s", if (v == "K10") "K10 attention (t = 9)" else "CLR", scopes[[sc]]), "", mt(w), "")
    x$cell <- sprintf("%+.2f [%+.2f, %+.2f]", x$d_enr, x$lo_enr, x$hi_enr)
    w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = e, t(stats::setNames(y$cell, y$readout))[, readouts, drop = FALSE], check.names = FALSE) }))
    md <- c(md, sprintf("### Difference in enrichment vs %s [95%% interval], %s", if (v == "K10") "K10 attention (t = 9)" else "CLR", scopes[[sc]]), "", mt(w), "") }
  writeLines(md, file.path(outdir, paste0(tag, "_tables.md"))); invisible(list(pm = pm, ct = ct))
}
