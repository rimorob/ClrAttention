#!/usr/bin/env Rscript
# Scoreboard data matrix (docs/plans/scoreboard_plan.md). One row per dataset x truth x metric x scope x method, descriptive only.
# Reads existing outputs (status from_existing), the four scoreboard job folders results/scoreboard/{m3d,ecomac,k562,beeline} (computed_now),
# and marks everything else not_run. Also writes the K10 / sota_table cross-check (docs/whitepapers/figures/scoreboard_crosscheck.csv).
# Output: docs/whitepapers/figures/scoreboard_data.csv
# Run:    Rscript analysis/scoreboard_matrix.R
out_dir <- "docs/whitepapers/figures"
M9 <- c("Attention", "CLR", "a3r", "a3r2", "a3r3", "GENIE3", "PIDC", "GRNBoost2", "abs_r")
old6 <- c(attention = "Attention", clr = "CLR", genie3 = "GENIE3", pidc = "PIDC", grnboost2 = "GRNBoost2", abs_pearson = "abs_r")
rd <- function(f) utils::read.csv(f, stringsAsFactors = FALSE)
L <- list()
add <- function(dataset, truth, metric, scope, method, value, base, source, status, note = "")
  L[[length(L) + 1L]] <<- data.frame(dataset = dataset, truth = truth, metric = metric, scope = scope, method = method,
    enrichment_over_random = value, base_rate = base, source_file = source, status = status, note = note, stringsAsFactors = FALSE)
cells <- list()   # every (dataset, truth, metric, scope) that must exist for all 9 methods
cell <- function(dataset, truth, metric, scope) cells[[length(cells) + 1L]] <<- data.frame(dataset = dataset, truth = truth, metric = metric, scope = scope, stringsAsFactors = FALSE)
PM <- "pooled_aupr_over_base"

## ---- pooled co-membership AUPR: the six methods of sota_table (generator input metrics.csv) --------------------------------------------
runs <- c(m3d = "results/baselines_m3d", ecomac = "results/baselines_ecomac", k562 = "results/baselines_k562ess_truthgenes_genie3")
bm <- lapply(runs, function(r) rd(file.path(r, "metrics.csv")))
for (ds in names(runs)) { x <- bm[[ds]]; src <- file.path(runs[[ds]], "metrics.csv")
  for (tn in unique(x$truth)) for (sc in c("raw", "hub_att")) {
    cell(ds, tn, PM, sc)
    for (m in names(old6)) { key <- if (sc == "raw") m else paste0(m, "_hubstrat"); y <- x[x$truth == tn & x$method == key, ]
      if (nrow(y) == 1L) add(ds, tn, PM, sc, old6[[m]], y$aupr / y$base_rate, y$base_rate, src, "from_existing",
        if (sc == "hub_att") "within-stratum percentile in 10 strata of pi_i + pi_j (sota_table hub-controlled)" else "") } } }
# ---- the explore_stack K10 / CLR / a3r / a3r2 / a3r3 readouts (pooled_levels.csv) -------------------------------------------------------
pl <- list(m3d = rd("results/explore_stack/m3d/pooled_levels.csv"), k562 = rd("results/explore_stack/k562/pooled_levels.csv"))
plname <- c(K10 = "Attention", CLR = "CLR", a3r = "a3r", a3r2 = "a3r2", a3r3 = "a3r3")
tmap <- list(m3d = c(regulons_SC = "SC", regulons_all = "all", complexes = "complexes"),
             k562 = c(DoRothEA_AB = "DoRothEA_AB", CollecTRI = "CollecTRI", CORUM_all = "CORUM", KD_response = "KD_response", CORUM_noRP = "CORUM_noRP", CORUM_RP = "CORUM_RP",
                      CORUM_small = "CORUM_small", CORUM_large = "CORUM_large", CORUM_merged = "CORUM_merged", CORUM_merged_noRP = "CORUM_merged_noRP"))
have <- function(ds, tn, sc, m) any(vapply(L, function(d) d$dataset == ds && d$truth == tn && d$metric == PM && d$scope == sc && d$method == m, NA))
for (ds in names(pl)) { x <- pl[[ds]]; src <- sprintf("results/explore_stack/%s/pooled_levels.csv", ds)
  for (tn in names(tmap[[ds]])) for (sc in c("raw", "hub_att", "hub_mi")) {
    ct <- tmap[[ds]][[tn]]; cell(ds, ct, PM, sc)
    for (rn in names(plname)) { y <- x[x$readout == rn & x$scope == sc & x$truth == tn, ]; stopifnot(nrow(y) == 1L)
      m <- plname[[rn]]
      if (sc != "hub_mi" && ct %in% bm[[ds]]$truth && m %in% c("Attention", "CLR")) next                # sota_table input is the source of record (cross-checked below)
      add(ds, ct, PM, sc, m, y$enrichment, y$base_rate, src, "from_existing",
          if (ct == "CORUM" && ds == "k562") "CORUM_all of the stack run (= CORUM raw)" else if (sc == "hub_mi") "within-stratum percentile in 10 strata of the raw-MI row mean (hub_mi)" else "") } } }
# ---- computed now: M3D has all five readouts already; EcoMAC SC / all raw + hub_att; K562 DoRothEA ABC -----------------------------------
pe <- rd("results/scoreboard/ecomac/pooled.csv"); pe <- pe[pe$readout %in% c("a3r", "a3r2", "a3r3"), ]
for (k in seq_len(nrow(pe))) if (pe$scope[k] != "hub_mi") add("ecomac", pe$truth[k], PM, pe$scope[k], pe$readout[k], pe$enrichment[k], pe$base_rate[k], "results/scoreboard/ecomac/pooled.csv", "computed_now",
  if (pe$scope[k] == "hub_att") "within-stratum percentile in 10 strata of pi_i + pi_j" else "")
pk <- rd("results/scoreboard/k562/pooled.csv"); pk <- pk[pk$truth == "DoRothEA_ABC" & pk$readout %in% c("a3r", "a3r2", "a3r3"), ]
for (k in seq_len(nrow(pk))) { if (pk$scope[k] == "hub_mi") next
  add("k562", "DoRothEA_ABC", PM, pk$scope[k], pk$readout[k], pk$enrichment[k], pk$base_rate[k], "results/scoreboard/k562/pooled.csv", "computed_now", "") }
# ---- complexes of E. coli (M3D, EcoMAC) and RPE1 CORUM: clr_logpt_test (|r|, and attention / CLR / a3r for EcoMAC and RPE1) -------------------
lp <- function(ds) rd(sprintf("results/clr_logpt_test/truth/%s/metrics.csv", ds))
lmap <- c(attention_t9 = "Attention", clr = "CLR", a3r_t9 = "a3r", abs_r = "abs_r")
for (spec in list(c("m3d", "complexes", "complexes"), c("ecomac", "complexes", "complexes"), c("rpe1", "CORUM_merged", "CORUM_merged"))) {
  ds <- spec[1]; x <- lp(ds); x <- x[x$truth == spec[2], ]; src <- sprintf("results/clr_logpt_test/truth/%s/metrics.csv", ds)
  for (sc in c("raw", "hub_att")) { cell(ds, spec[3], PM, sc)
    for (mm in names(lmap)) { m <- lmap[[mm]]; key <- if (sc == "raw") mm else paste0(mm, "_hubstrat"); y <- x[x$method == key, ]
      if (ds == "m3d" && m != "abs_r") next                                                    # M3D: the five readouts come from pooled_levels
      if (nrow(y) == 1L) add(ds, spec[3], PM, sc, m, y$aupr / y$base_rate, y$base_rate, src, "from_existing",
        if (ds == "rpe1") "RPE1 Perturb-seq, merged CORUM (descriptive; script exists, clr_logpt_test.R)" else if (sc == "hub_att") "hub_att = pi stratifier (clr_logpt_test _hubstrat)" else "") } } }
# EcoMAC / RPE1 a3r2, a3r3 are not_run (filled below)

## ---- direct TF -> target edges (RegulonDB TF-node benchmark) and per-regulon mean AP -----------------------------------------------------
for (ds in c("m3d", "ecomac")) { x <- bm[[ds]]; tf <- rd(sprintf("results/scoreboard/%s/tf_edges.csv", ds)); pr0 <- rd(file.path(runs[[ds]], "per_regulon.csv"))
  pe_ <- rd(sprintf("results/scoreboard/%s/per_regulon_summary.csv", ds))
  for (tn in c("SC", "all")) {
    for (mt in c("tf_aupr_ratio", "tf_auroc_median")) { cell(ds, tn, mt, "raw")
      for (m in names(old6)) { y <- x[x$truth == tn & x$method == m, ]; add(ds, tn, mt, "raw", old6[[m]], y[[mt]], NA, file.path(runs[[ds]], "metrics.csv"), "from_existing",
        if (mt == "tf_auroc_median") "median per-TF AUROC; chance = 0.50, not an enrichment" else "directed TF -> target edges, 199 TFs, pooled AUPR / positive rate") }
      for (m in c("a3r", "a3r2", "a3r3")) { y <- tf[tf$readout == m & tf$truth == tn, ]; add(ds, tn, mt, "raw", m, y[[mt]], NA, sprintf("results/scoreboard/%s/tf_edges.csv", ds), "computed_now",
        if (mt == "tf_auroc_median") "median per-TF AUROC; chance = 0.50, not an enrichment" else "directed TF -> target edges, 199 TFs, pooled AUPR / positive rate") } }
    for (mt in c("per_regulon_mean_ap_over_base", "per_regulon_mean_ap", "per_regulon_mean_ap_norm")) { cell(ds, tn, mt, "raw")
      nt <- switch(mt, per_regulon_mean_ap_over_base = "mean over regulons of AP / mean of the regulons' base rates (1.00 = random); regulons with >= 10 positive pairs",
                   per_regulon_mean_ap = "mean over regulons of AP (regulons with >= 10 positive pairs)", per_regulon_mean_ap_norm = "mean of (AP - b) / (1 - b) per regulon; 0 = random, not an enrichment")
      for (m in names(old6)) { y <- pr0[pr0$method == m & pr0$truth == tn & !is.na(pr0$ap), ]
        v <- switch(mt, per_regulon_mean_ap_over_base = mean(y$ap) / mean(y$base_rate), per_regulon_mean_ap = mean(y$ap), per_regulon_mean_ap_norm = mean((y$ap - y$base_rate) / (1 - y$base_rate)))
        add(ds, tn, mt, "raw", old6[[m]], v, mean(y$base_rate), file.path(runs[[ds]], "per_regulon.csv"), "from_existing", paste0(nt, "; n = ", nrow(y))) }
      for (m in c("a3r", "a3r2", "a3r3")) { y <- pe_[pe_$readout == m & pe_$truth == tn, ]
        v <- switch(mt, per_regulon_mean_ap_over_base = y$enrichment, per_regulon_mean_ap = y$mean_ap, per_regulon_mean_ap_norm = y$mean_ap_norm)
        add(ds, tn, mt, "raw", m, v, y$mean_base_rate, sprintf("results/scoreboard/%s/per_regulon_summary.csv", ds), "computed_now", paste0(nt, "; n = ", y$n_regulons)) } } } }

## ---- BEELINE: median AUPRC ratio over the complete primary case x truth cells (14 cases x BEELINE primary truths) and over the 14 cell-type ChIP cases ---
bb <- rd("results/beeline_baselines/metrics.csv"); sb <- rd("results/scoreboard/beeline/summary.csv")
primary <- c("celltype_ChIP", "nonspecific_ChIP", "STRING", "LOFGOF")
cellsb <- unique(bb[, c("dataset", "N", "truth")]); comp <- merge(cellsb, stats::aggregate(method ~ dataset + N + truth, bb, length)); comp <- comp[comp$method == length(unique(bb$method)), c("dataset", "N", "truth")]
rcb <- merge(bb, comp)
bmap <- c(abs_pearson = "abs_r", clr = "CLR", clr_ns = "CLR-NS", genie3 = "GENIE3", grnboost2 = "GRNBoost2", pidc = "PIDC")
for (ts in c("primary", "celltype_ChIP")) { mt <- if (ts == "primary") "beeline_median_auprc_ratio_primary" else "beeline_median_auprc_ratio_celltypeChIP"; cell("beeline", "all_14_cases", mt, "raw")
  x <- rcb[if (ts == "primary") rcb$truth %in% primary else rcb$truth == ts, ]
  nt <- if (ts == "primary") "median over the complete case x truth cells (14 cases x BEELINE primary truths: cell-type ChIP, non-specific ChIP, STRING, mESC LOF/GOF = 44 cells); the '14 cases' of the figures" else "median over the 14 cases, cell-type-specific ChIP truth only"
  for (m in names(bmap)) { y <- x[x$method == m, ]; add("beeline", "all_14_cases", mt, "raw", bmap[[m]], stats::median(y$auprc_ratio), NA, "results/beeline_baselines/metrics.csv", "from_existing", paste0(nt, "; n = ", nrow(y), if (m == "clr_ns") "; extra method, not a scoreboard column" else "")) }
  for (m in c("attention_tstar", "attention_t9", "a3r", "a3r2", "a3r3")) { y <- sb[sb$truths == ts & sb$method == m, ]
    nm <- c(attention_tstar = "Attention", attention_t9 = "Attention_t9", a3r = "a3r", a3r2 = "a3r2", a3r3 = "a3r3")[[m]]
    add("beeline", "all_14_cases", mt, "raw", nm, y$median_auprc_ratio, NA, "results/scoreboard/beeline/summary.csv", "computed_now",
        paste0(nt, "; n = ", y$n_cells, "; ", switch(m, attention_tstar = "soft10 attention on CLR-NS at the held-out t* of each case (BEELINE protocol; equals att_clr_absr_blomavg of beeline_blom_test, gated)",
          attention_t9 = "soft10 attention on CLR-NS at t = 9; extra method, not a scoreboard column", "a3r family on CLR-NS at t = 9 in every layer, not tuned"))) } }
cells[[length(cells) + 1L]] <- data.frame(dataset = "beeline", truth = "all_14_cases", metric = "beeline_median_auprc_ratio_primary", scope = "raw")   # keep Attention_t9 / CLR-NS rows

## ---- DREAM5 yeast network 4 (held-out test, run once): not_run for every method ---------------------------------------------------------
cell("dream5_yeast", "network4_truth", PM, "raw")

D <- do.call(rbind, L); C <- unique(do.call(rbind, cells))
# complete: every cell x 9 methods; missing -> not_run
grid <- merge(C, data.frame(method = M9), by = NULL)
key <- function(x) paste(x$dataset, x$truth, x$metric, x$scope, x$method, sep = "|")
miss <- grid[!key(grid) %in% key(D), ]
if (nrow(miss)) {
  nr <- function(r) if (r$dataset == "dream5_yeast") "held-out test truth (run once); not scored for any method here" else if (r$dataset == "beeline") "not run" else
    if (r$scope %in% c("hub_att", "hub_mi") && r$metric != PM) "no hub-stratified run exists for this metric" else
    if (r$metric %in% c("tf_aupr_ratio", "tf_auroc_median", "per_regulon_mean_ap_over_base", "per_regulon_mean_ap", "per_regulon_mean_ap_norm")) "not run" else
    if (r$truth %in% c("complexes", "CORUM_merged") || r$dataset == "k562" && grepl("^CORUM_", r$truth)) "no competitor run on this truth (GENIE3, PIDC, GRNBoost2 were not run on complex truths)" else
    if (r$scope == "hub_mi") "no competitor hub_mi run (the sota_table stratifier is hub_att)" else "not run"
  D <- rbind(D, data.frame(miss, enrichment_over_random = NA_real_, base_rate = NA_real_, source_file = "", status = "not_run",
    note = vapply(seq_len(nrow(miss)), function(i) nr(miss[i, ]), ""), stringsAsFactors = FALSE)) }
# base_rate for not_run cells: take it from another method of the same cell when known
for (i in which(is.na(D$base_rate) & D$status == "not_run")) { j <- which(D$dataset == D$dataset[i] & D$truth == D$truth[i] & D$metric == D$metric[i] & D$scope == D$scope[i] & !is.na(D$base_rate)); if (length(j)) D$base_rate[i] <- D$base_rate[j[1]] }
D <- D[order(D$dataset, D$truth, D$metric, D$scope, match(D$method, c(M9, "Attention_t9", "CLR-NS"))), ]
stopifnot(!anyDuplicated(key(D)))
utils::write.csv(D, file.path(out_dir, "scoreboard_data.csv"), row.names = FALSE)

## ---- cross-check: K10 / CLR of pooled_levels against the sota_table (generator input metrics.csv and the printed table) ---------------------------------
md <- readLines("docs/whitepapers/sota_table.md"); sect <- cumsum(grepl("^## ", md)); parse_tab <- function(s) { x <- md[sect == s & grepl("^\\| (M3D|EcoMAC|K562)", md)]
  do.call(rbind, lapply(strsplit(x, "\\|"), function(p) { p <- trimws(p[-1]); data.frame(dataset = p[1], truth = gsub(" ", "_", p[2]), printed = I(list(p[-(1:3)])), stringsAsFactors = FALSE) })) }
raw_tab <- parse_tab(1L); hub_tab <- parse_tab(2L)
pr_cols <- list(raw = c("Attention", "CLR", "GENIE3", "PIDC", "GRNBoost2", "abs_r", "Hub"), hub_att = c("Attention", "CLR", "GENIE3", "PIDC", "GRNBoost2", "abs_r"))
printed <- function(ds, tn, sc, m) { tb <- if (sc == "raw") raw_tab else hub_tab; r <- tb[tb$dataset == ds & tb$truth == tn, ]; v <- r$printed[[1]][match(m, pr_cols[[sc]])]
  as.numeric(gsub("[^0-9.]", "", v)) }
cc <- list(); dsn <- c(M3D = "m3d", EcoMAC = "ecomac", K562 = "k562")
for (ds in c("m3d", "k562")) for (tn in names(tmap[[ds]])) for (sc in c("raw", "hub_att")) {
  ct <- tmap[[ds]][[tn]]; if (!ct %in% bm[[ds]]$truth) next
  for (rn in c("K10", "CLR")) { m <- plname[[rn]]
    y <- pl[[ds]][pl[[ds]]$readout == rn & pl[[ds]]$scope == sc & pl[[ds]]$truth == tn, ]
    b <- bm[[ds]][bm[[ds]]$truth == ct & bm[[ds]]$method == (if (sc == "raw") c(Attention = "attention", CLR = "clr")[[m]] else paste0(c(Attention = "attention", CLR = "clr")[[m]], "_hubstrat")), ]
    pv <- printed(c(m3d = "M3D", k562 = "K562")[[ds]], ct, sc, m)
    cc[[length(cc) + 1L]] <- data.frame(dataset = ds, truth = ct, scope = sc, method = m, stack_enrichment = y$enrichment, sota_metrics_enrichment = b$aupr / b$base_rate,
      abs_diff_aupr = abs(y$aupr - b$aupr), sota_table_printed = pv, stack_rounded_2dp = round(y$enrichment, 2), matches_printed_2dp = isTRUE(all.equal(round(y$enrichment, 2), pv)),
      matches_full_precision_1e9 = abs(y$aupr - b$aupr) < 1e-9, stringsAsFactors = FALSE) } }
cc <- do.call(rbind, cc); utils::write.csv(cc, file.path(out_dir, "scoreboard_crosscheck.csv"), row.names = FALSE)
cat("cross-check rows:", nrow(cc), "| match printed 2dp:", sum(cc$matches_printed_2dp), "| match full precision 1e-9:", sum(cc$matches_full_precision_1e9), "\n")
print(cc[!cc$matches_printed_2dp | !cc$matches_full_precision_1e9, ], row.names = FALSE)
cat("status counts:\n"); print(table(D$status)); print(table(D$dataset, D$status))
