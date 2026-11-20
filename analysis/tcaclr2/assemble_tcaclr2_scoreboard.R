#!/usr/bin/env Rscript
# tcaCLR2 on the scoreboard rows, step 3: assemble every tcaCLR2 cell of the scoreboard (docs/whitepapers/figures/scoreboard_data.csv schema, method key "tcaCLR2") from
#   results/tcaclr2_scoreboard/cells_{m3d,ecomac,k562,rpe1,beeline}.csv   (score_scoreboard.R; its tcaCLR t = 9 gate has passed)
#   results/tcaclr2/score_{net4_cp_yeast,armA_cp_human,armB_cyc2008}/results.csv   (score_tcaclr2.R, already run: pooled AUPR of tcaCLR2 and tcaCLR on the same truths)
# Rows of the DREAM5 yeast, HepG2 and Kemmeren cells: enrichment = pooled AUPR / base rate (n_pos / n_pairs), as the scoreboard; the registered tcaCLR value of the same file is checked
# against the stored scoreboard cell (scoreboard_data.csv for net4; docs/tcaclr_confirm_results/{armA,armB}/metrics.csv for HepG2 / Kemmeren) to 1e-9.
# Also cross-checks the family runs against the earlier results/tcaclr2 AUPR values (K562 CORUM, M3D SC, all) and writes t2_by_dataset.csv for the figure footnote.
#   Rscript analysis/tcaclr2/assemble_tcaclr2_scoreboard.R [--out results/tcaclr2_scoreboard]
# Output: cells_tcaclr2_all.csv (tcaCLR2 cells), cells_reference_tcaclr_all.csv (tcaCLR t = 9 recomputed, for the audit), crosscheck.csv, t2_by_dataset.csv.
source("analysis/tcaclr2/tsb_lib.R")
topt <- tsb_args(list(out = TSB$out)); tsb_refuse_frozen(topt$out); say <- tsb_say()
SBD <- utils::read.csv(TSB$sb_csv, stringsAsFactors = FALSE); PM <- "pooled_aupr_over_base"
fam <- c("m3d", "ecomac", "k562", "rpe1", "beeline"); fs <- file.path(topt$out, sprintf("cells_%s.csv", fam)); miss <- fam[!file.exists(fs)]
if (length(miss)) stop("missing family cells: ", paste(miss, collapse = ", "))
cells <- do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
chk <- list(); ck <- function(what, a, b, tol = 1e-9) { chk[[length(chk) + 1L]] <<- data.frame(check = what, value = a, reference = b, abs_diff = abs(a - b), pass = abs(a - b) <= tol, stringsAsFactors = FALSE) }

# ---- the three rows already scored by score_tcaclr2.R ---------------------------------------------------------------------------------------
extra <- list(); ref <- list()
one <- function(truth_key, scope_key, ds, truth_name, stored_fn) {
  r <- utils::read.csv(sprintf("results/tcaclr2/score_%s/results.csv", truth_key), stringsAsFactors = FALSE); r <- r[r$scope == scope_key, ]; stopifnot(nrow(r) == 1L)
  base <- r$n_pos / r$n_pairs
  ck(sprintf("%s: tcaCLR (t = 9) AUPR / base of score_tcaclr2.R vs the stored scoreboard cell", ds), r$aupr_tcaCLR / base, stored_fn(base))
  note <- sprintf("tcaCLR2: depth t2 = %d (edge %.5f; rule of docs/plans/tcaclr2_definition.md); exploratory, outside the preregistration", r$t_tcaCLR2, r$edge)
  extra[[length(extra) + 1L]] <<- tsb_row(ds, truth_name, PM, "raw", "tcaCLR2", r$aupr_tcaCLR2 / base, base, sprintf("results/tcaclr2/score_%s/results.csv", truth_key), note)
  ref[[length(ref) + 1L]] <<- tsb_row(ds, truth_name, PM, "raw", "tcaCLR_t9", r$aupr_tcaCLR / base, base, sprintf("results/tcaclr2/score_%s/results.csv", truth_key), "tcaCLR at t = 9 of score_tcaclr2.R")
}
one("net4_cp_yeast", "main", "dream5_yeast", "complexes_main", function(b) tsb_stored(SBD, "dream5_yeast", "complexes_main", PM, "raw"))
conf <- function(arm) { m <- utils::read.csv(sprintf("docs/tcaclr_confirm_results/%s/metrics.csv", arm), stringsAsFactors = FALSE); m <- m[m$truth == "main" & m$scope == "raw" & m$method == "tcaCLR", ]; stopifnot(nrow(m) == 1L); m$aupr / m$base_rate }
one("armA_cp_human", "main", "hepg2", "complexes_main", function(b) conf("armA")); one("armB_cyc2008", "main", "kemmeren", "complexes_main", function(b) conf("armB"))
cells <- rbind(cells, do.call(rbind, extra), do.call(rbind, ref))

# ---- cross-check of the family runs against the earlier AUPR of results/tcaclr2 (K562 CORUM, M3D SC / all) ----------------------------------------------------
for (spec in list(c("k562_corum", "CORUM", "k562", "CORUM"), c("m3d_tf_regulons", "SC", "m3d", "SC"), c("m3d_tf_regulons", "all", "m3d", "all"))) {
  r <- utils::read.csv(sprintf("results/tcaclr2/score_%s/results.csv", spec[1]), stringsAsFactors = FALSE); r <- r[r$scope == spec[2], ]; base <- r$n_pos / r$n_pairs
  for (m in c("tcaCLR2", "tcaCLR_t9")) { v <- cells$enrichment_over_random[cells$dataset == spec[3] & cells$truth == spec[4] & cells$metric == PM & cells$scope == "raw" & cells$method == m]
    stopifnot(length(v) == 1L); ck(sprintf("%s %s: %s AUPR (family run, scoreboard path) vs results/tcaclr2 (grid path)", spec[3], spec[4], m), v * base, if (m == "tcaCLR2") r$aupr_tcaCLR2 else r$aupr_tcaCLR) } }
chk <- do.call(rbind, chk); utils::write.csv(chk, file.path(topt$out, "crosscheck.csv"), row.names = FALSE); print(chk, row.names = FALSE, digits = 6)
if (!all(chk$pass)) stop("cross-check failed (crosscheck.csv)")

utils::write.csv(cells[cells$method == "tcaCLR2", ], file.path(topt$out, "cells_tcaclr2_all.csv"), row.names = FALSE)
utils::write.csv(cells[cells$method == "tcaCLR_t9", ], file.path(topt$out, "cells_reference_tcaclr_all.csv"), row.names = FALSE)

# ---- t2 by data set (figure footnote) ---------------------------------------------------------------------------------------------------------------------
rd <- function(f) tsb_edge_read(f)
bl <- utils::read.csv(file.path(topt$out, "beeline_t2.csv"), stringsAsFactors = FALSE)
t2tab <- data.frame(label = c("M3D", "EcoMAC", "K562", "RPE1", "yeast net4", "HepG2", "Kemmeren", "BEELINE (14 cases, CLR-NS)"),
  edge = c(rd("results/tcaclr2/edge_m3d_tf_regulons.tsv")$edge, rd(file.path(topt$out, "edge_ecomac.tsv"))$edge, rd("results/tcaclr2/edge_k562_corum.tsv")$edge, rd(file.path(topt$out, "edge_rpe1.tsv"))$edge,
           rd("results/tcaclr2/edge_net4_cp_yeast.tsv")$edge, rd("results/tcaclr2/edge_armA_cp_human.tsv")$edge, rd("results/tcaclr2/edge_armB_cyc2008.tsv")$edge, NA),
  t2 = c(rd("results/tcaclr2/edge_m3d_tf_regulons.tsv")$t2, rd(file.path(topt$out, "edge_ecomac.tsv"))$t2, rd("results/tcaclr2/edge_k562_corum.tsv")$t2, rd(file.path(topt$out, "edge_rpe1.tsv"))$t2,
         rd("results/tcaclr2/edge_net4_cp_yeast.tsv")$t2, rd("results/tcaclr2/edge_armA_cp_human.tsv")$t2, rd("results/tcaclr2/edge_armB_cyc2008.tsv")$t2, NA), stringsAsFactors = FALSE)
t2tab$t2 <- as.character(t2tab$t2); t2tab$t2[8] <- if (min(bl$t2) == max(bl$t2)) as.character(min(bl$t2)) else sprintf("%d to %d", min(bl$t2), max(bl$t2))
utils::write.csv(t2tab, file.path(topt$out, "t2_by_dataset.csv"), row.names = FALSE); print(t2tab, row.names = FALSE)
say("wrote cells_tcaclr2_all.csv (", sum(cells$method == "tcaCLR2"), " tcaCLR2 cells), cells_reference_tcaclr_all.csv, crosscheck.csv, t2_by_dataset.csv")
