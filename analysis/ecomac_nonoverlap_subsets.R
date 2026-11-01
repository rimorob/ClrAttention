#!/usr/bin/env Rscript
# Sensitivity analysis (NOT preregistered): array subsets of EcoMAC defined by similarity to M3D
# (results/ecomac_m3d_overlap/per_array_best_match.csv, headline call = centred rule AND mutual
# nearest neighbour; no truth involved). Writes 1-based array indices (one per line) for
# analysis/baselines.R --array_subset, and the platform composition.
#   A        all 2,198 arrays (gate: must reproduce the existing EcoMAC numbers)
#   B        not called (1,370)
#   C        called (828)
#   Bg3154   B restricted to GPL3154 (735)
#   Bmix     B restricted to GPL199 and no-GEO-record arrays (635)
#   R1370    (reference) random 1,370 arrays from A (seed 20261101): effect of simply dropping arrays
out <- "results/ecomac_nonoverlap/subsets"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
pm <- utils::read.csv("results/ecomac_m3d_overlap/per_array_best_match.csv", stringsAsFactors = FALSE)
ap <- utils::read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)
smp <- utils::read.delim("data/ecomac/samples.tsv", stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(pm) == 2198, nrow(ap) == 2198, identical(as.integer(pm$array), seq_len(2198L)),
          identical(pm$id, ap$sample), identical(smp[["GEO accession"]], ap$geo),
          identical(sub(" .*", "", ap$platform) == "Antisense", pm$platform == "GPL199"),
          sum(pm$call) == 828L)
called <- pm$call
subs <- list(A = seq_len(2198L), B = which(!called), C = which(called),
             Bg3154 = which(!called & pm$platform == "GPL3154"),
             Bmix = which(!called & pm$platform != "GPL3154"))
set.seed(20261101L); subs$R1370 <- sort(sample(2198L, 1370L))
for (n in names(subs)) writeLines(as.character(subs[[n]]), file.path(out, paste0(n, ".txt")))
comp <- do.call(rbind, lapply(names(subs), function(n) { ix <- subs[[n]]; p <- ap$platform[ix]
  data.frame(subset = n, n = length(ix), antisense_GPL199 = sum(grepl("GPL199", p)), genome2_GPL3154 = sum(grepl("GPL3154", p)),
             no_geo_record = sum(!grepl("GPL", p)), post_build_GSE20000plus = sum(pm$post_build[ix]),
             n_series = length(unique(pm$series[ix])), n_authors = length(unique(ap$author[ix]))) }))
utils::write.csv(comp, "results/ecomac_nonoverlap/subset_composition.csv", row.names = FALSE)
print(comp); print(table(ap$platform, pm$call))
