#!/usr/bin/env Rscript
# EXPLORATORY summary of analysis/robust_z/robust_z.R caches: point estimates and intervals only,
# no p-values, no winner selection. Usage: Rscript analysis/robust_z/robust_z_summarize.R <dir> [baselines_dir]
args <- commandArgs(TRUE); dir <- args[1]
bdir <- if (length(args) >= 2) args[2] else NA
cd <- file.path(dir, "cache"); meta <- readRDS(file.path(cd, "meta.rds"))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
VARS <- c("V0", "V1", "V2", "V3", "V4"); KINDS <- c(pair = "pair", att9 = sprintf("att%d", meta$opt$t), attT = "attT")
ld <- function(v, k) { f <- file.path(cd, sprintf("%s__%s.rds", v, k)); if (file.exists(f)) readRDS(f) else NULL }
R <- list(); for (v in c(VARS, "V0x")) for (k in names(KINDS)) { x <- ld(v, KINDS[[k]]); if (!is.null(x)) R[[paste(v, k)]] <- x }
truths <- names(meta$base_rate); br <- meta$base_rate
set.seed(20261111L)

## 1. pooled metrics ----------------------------------------------------------------------------------
pm <- do.call(rbind, lapply(names(R), function(key) { v <- sub(" .*", "", key); k <- sub(".* ", "", key); r <- R[[key]]
  do.call(rbind, lapply(truths, function(e) data.frame(truth = e, variant = v, kind = k, base_rate = br[[e]], aupr = r$pooled[[e]]$aupr,
    enrichment = r$pooled[[e]]$aupr / br[[e]],
    enr_hub_att = if (!is.null(r$hub_att)) r$hub_att[[e]]$aupr / br[[e]] else NA, enr_hub_mi = if (!is.null(r$hub_mi)) r$hub_mi[[e]]$aupr / br[[e]] else NA))) }))
utils::write.csv(pm, file.path(dir, "pooled_metrics.csv"), row.names = FALSE)

## 2. paired delete-half jackknife contrasts (difference in enrichment; 95% = d +- 1.96 sd of paired replicate differences)
getp <- function(key, scope, e) { r <- R[[key]]; if (is.null(r) || is.null(r[[scope]])) NULL else r[[scope]][[e]] }
contr <- list()
add <- function(scope, e, type, a, b) { pa <- getp(a, scope, e); pb <- getp(b, scope, e); if (is.null(pa) || is.null(pb)) return()
  d <- (pa$aupr - pb$aupr) / br[[e]]; dj <- (pa$jack - pb$jack) / br[[e]]; s <- stats::sd(dj)
  contr[[length(contr) + 1L]] <<- data.frame(scope = switch(scope, pooled = "raw", hub_att = "hub-strat (pi of V0 attention)", hub_mi = "hub-strat (MI row mean)"),
    truth = e, contrast = type, a = a, b = b, d_enrichment = d, lo = d - 1.96 * s, hi = d + 1.96 * s, n_rep = length(dj)) }
for (scope in c("pooled", "hub_att", "hub_mi")) for (e in truths) {
  for (v in c("V0x", "V1", "V2", "V3", "V4")) { add(scope, e, "pairwise: variant - V0", paste(v, "pair"), "V0 pair")
    add(scope, e, "attention t=9: variant - V0", paste(v, "att9"), "V0 att9"); if (v != "V0x") add(scope, e, "attention t*: variant - V0", paste(v, "attT"), "V0 attT") }
  for (v in VARS) add(scope, e, "attention t=9 minus CLR(V0 pairwise)", paste(v, "att9"), "V0 pair")
}
ct <- do.call(rbind, contr); utils::write.csv(ct, file.path(dir, "contrasts_pooled.csv"), row.names = FALSE)

## 3. macro (per-regulon / per-complex AP) ----------------------------------------------------------------
msrc <- unique(meta$macro_of); mrows <- list(); mlev <- list()
boot_ci <- function(d) { b <- replicate(2000, mean(d[sample.int(length(d), replace = TRUE)])); stats::quantile(b, c(.025, .975), names = FALSE) }
for (e in truths) { src <- meta$macro_of[[e]]; sel <- meta$macro_sel[[e]]
  get_tab <- function(key) { r <- R[[key]]; if (is.null(r) || is.null(r$macro)) return(NULL); t <- r$macro[[src]]$table; t <- t[!is.na(t$ap), ]; if (!is.null(sel)) t <- t[t$regulon %in% sel, ]; t }
  get_jk <- function(key) { r <- R[[key]]; j <- r$macro[[src]]$jack; if (!is.null(sel)) j <- j[j$regulon %in% sel, ]; j }
  for (key in names(R)) { t <- get_tab(key); if (is.null(t)) next
    mlev[[length(mlev) + 1L]] <- data.frame(truth = e, variant = sub(" .*", "", key), kind = sub(".* ", "", key), n_regulons = nrow(t), mean_ap = mean(t$ap), median_ap = stats::median(t$ap),
      mean_ap_norm = mean((t$ap - t$base_rate) / (1 - t$base_rate))) }
  mc <- function(type, a, b) { ta <- get_tab(a); tb <- get_tab(b); if (is.null(ta) || is.null(tb)) return()
    m <- merge(ta[, c("regulon", "ap")], tb[, c("regulon", "ap")], by = "regulon", suffixes = c(".a", ".b")); d <- m$ap.a - m$ap.b
    ja <- get_jk(a); jb <- get_jk(b); jj <- merge(ja[, c("rep", "regulon", "ap")], jb[, c("rep", "regulon", "ap")], by = c("rep", "regulon"), suffixes = c(".a", ".b"))
    jj <- jj[jj$regulon %in% m$regulon, ]; R_ <- length(unique(jj$rep)); okc <- tapply(!is.na(jj$ap.a) & !is.na(jj$ap.b), jj$regulon, sum); fixed <- names(okc)[okc == R_]
    jj <- jj[jj$regulon %in% fixed, ]; th <- vapply(split(jj$ap.a - jj$ap.b, jj$rep), mean, 0); se <- sqrt(4 * mean((th - mean(th))^2))   # delete-20%: (n-d)/d = 4
    est <- mean(d[m$regulon %in% fixed]); bc <- boot_ci(d)
    mrows[[length(mrows) + 1L]] <<- data.frame(truth = e, contrast = type, a = a, b = b, n_regulons = length(d), mean_d_ap = mean(d), boot_lo_regulons = bc[1], boot_hi_regulons = bc[2],
      n_regulons_genejack = length(fixed), mean_d_genejack_set = est, genejack_lo = est - 1.96 * se, genejack_hi = est + 1.96 * se, n_rep = R_, frac_a_better = mean(d > 0), frac_tied = mean(d == 0)) }
  for (v in c("V1", "V2", "V3", "V4")) { mc("pairwise: variant - V0", paste(v, "pair"), "V0 pair"); mc("attention t=9: variant - V0", paste(v, "att9"), "V0 att9") }
  for (v in VARS) mc("attention t=9 minus CLR(V0 pairwise)", paste(v, "att9"), "V0 pair")
}
mt_ <- do.call(rbind, mrows); ml_ <- do.call(rbind, mlev)
utils::write.csv(mt_, file.path(dir, "contrasts_macro.csv"), row.names = FALSE); utils::write.csv(ml_, file.path(dir, "macro_levels.csv"), row.names = FALSE)

## 4. hub drift ----------------------------------------------------------------------------------------------
hd <- do.call(rbind, lapply(names(R), function(key) { r <- R[[key]]; if (is.null(r$drift)) return(NULL)
  cbind(variant = sub(" .*", "", key), kind = sub(".* ", "", key), r$drift) }))
utils::write.csv(hd, file.path(dir, "hub_drift.csv"), row.names = FALSE)

## 5. checks against the project's own baselines run -----------------------------------------------------------
chk <- NULL
if (!is.na(bdir) && file.exists(file.path(bdir, "metrics.csv"))) { bm <- utils::read.csv(file.path(bdir, "metrics.csv"), stringsAsFactors = FALSE)
  mp <- c(CORUM_all = "CORUM"); tn <- function(e) if (e %in% names(mp)) mp[[e]] else e
  chk <- do.call(rbind, lapply(truths, function(e) { b <- bm[bm$truth == tn(e), ]; if (!nrow(b)) return(NULL)
    data.frame(truth = e, V0_pair_aupr = pm$aupr[pm$truth == e & pm$variant == "V0" & pm$kind == "pair"], baselines_clr_aupr = b$aupr[b$method == "clr"],
               V0_att_aupr = pm$aupr[pm$truth == e & pm$variant == "V0" & pm$kind == "att9"], baselines_attention_aupr = b$aupr[b$method == "attention"]) }))
  chk$max_abs_diff <- pmax(abs(chk$V0_pair_aupr - chk$baselines_clr_aupr), abs(chk$V0_att_aupr - chk$baselines_attention_aupr))
  utils::write.csv(chk, file.path(dir, "check_vs_baselines.csv"), row.names = FALSE) }

## 6. markdown fragments -----------------------------------------------------------------------------------------
mt <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
fmt <- function(d, lo, hi, dig = 3) sprintf("%+.*f [%+.*f, %+.*f]", dig, d, dig, lo, dig, hi)
md <- c(sprintf("# Tables for %s (generated by analysis/robust_z/robust_z_summarize.R; EXPLORATORY)", meta$dataset), "")
lv <- function(kind, col = "enrichment") { x <- pm[pm$kind == kind & pm$variant %in% VARS, ]
  do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = sprintf("%s (base %.4f)", e, br[[e]]), t(stats::setNames(sprintf("%.2f", y[[col]]), y$variant))[, VARS, drop = FALSE], check.names = FALSE) })) }
md <- c(md, "## Pooled enrichment (AUPR / base rate), pairwise score", "", mt(lv("pair")), "", "## Pooled enrichment, attention t = 9", "", mt(lv("att9")), "",
        "## Pooled enrichment, attention at each variant's own held-out depth", "", mt(lv("attT")), "")
for (scope in unique(ct$scope)) for (typ in unique(ct$contrast)) { x <- ct[ct$scope == scope & ct$contrast == typ, ]; if (!nrow(x)) next
  x$cell <- fmt(x$d_enrichment, x$lo, x$hi); x$a <- sub(" (pair|att9|attT)$", "", x$a)
  w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = e, t(stats::setNames(y$cell, y$a))[, unique(x$a), drop = FALSE], check.names = FALSE) }))
  md <- c(md, sprintf("## Paired contrast, %s, %s (difference in enrichment [95%% delete-half gene jackknife])", typ, scope), "", mt(w), "") }
ml_$mean_ap <- sprintf("%.3f", ml_$mean_ap)
for (kd in c("pair", "att9")) { x <- ml_[ml_$kind == kd, ]
  w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = sprintf("%s (n=%d)", e, y$n_regulons[1]), t(stats::setNames(y$mean_ap, y$variant))[, VARS, drop = FALSE], check.names = FALSE) }))
  md <- c(md, sprintf("## Macro: mean per-regulon AP, %s", if (kd == "pair") "pairwise" else "attention t = 9"), "", mt(w), "") }
for (typ in unique(mt_$contrast)) { x <- mt_[mt_$contrast == typ, ]; x$cell <- sprintf("%+.3f; reg [%+.3f, %+.3f]; gene [%+.3f, %+.3f]", x$mean_d_ap, x$boot_lo_regulons, x$boot_hi_regulons, x$genejack_lo, x$genejack_hi)
  x$a <- sub(" (pair|att9)$", "", x$a)
  w <- do.call(rbind, lapply(truths, function(e) { y <- x[x$truth == e, ]; data.frame(truth = e, t(stats::setNames(y$cell, y$a))[, unique(x$a), drop = FALSE], check.names = FALSE) }))
  md <- c(md, sprintf("## Macro contrast, %s (mean difference in per-regulon AP; regulon bootstrap and delete-20%% gene jackknife intervals)", typ), "", mt(w), "") }
hd2 <- hd; for (cc in c("share_touching_top5pct", "chance_share", "mean_endpoint_rank", "frac_positive_score")) hd2[[cc]] <- sprintf("%.3f", hd2[[cc]])
md <- c(md, "## Hub drift (top pairs of each score against two connectedness rankings; chance share = 0.0975, chance mean rank = 0.5)", "",
        mt(hd2[, c("variant", "kind", "stratifier", "share_touching_top5pct", "mean_endpoint_rank", "frac_positive_score")]), "")
md <- c(md, "## Zero-scale guard (rows with scale 0)", "", mt(data.frame(variant = names(meta$guard), rows = meta$guard)), "")
if (!is.null(chk)) md <- c(md, "## V0 against the project's baselines run", "", mt(chk), "")
rs <- utils::read.csv(file.path(dir, "row_scale_spearman.csv")); rs$spearman <- sprintf("%+.3f", rs$spearman)
md <- c(md, "## Row scale versus number of strong partners (Spearman)", "", mt(rs[, c("row_quantity", "versus", "spearman", "n_genes")]), "")
rd <- utils::read.csv(file.path(dir, "row_scale_deciles.csv")); for (cc in names(rd)[-(1:2)]) rd[[cc]] <- signif(rd[[cc]], 3)
md <- c(md, "## Row scale by decile of the number of globally strong partners", "", mt(rd), "")
if (file.exists(file.path(dir, "tstar_summary.csv"))) { ts <- utils::read.csv(file.path(dir, "tstar_summary.csv")); md <- c(md, "## Label-free held-out depth", "", mt(ts), "") }
writeLines(md, file.path(dir, "tables.md")); say("done")
