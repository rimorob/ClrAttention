# EcoMAC / M3D overlap audit.  Label-free: no truth file is opened.
#
# Question: do the EcoMAC (Carrera et al. 2014; 2,198 arrays) and E. coli M3D v4 Build 6
# (907 arrays, 466 replicate-averaged experiments) compendia share arrays?
#
# Thresholds fixed BEFORE any similarity was computed (see docs/plans/ecomac_m3d_overlap.md):
#   primary near-duplicate rule: Pearson r > 0.99 between an EcoMAC array and some M3D array,
#     over the shared genes, on the log2 scale after per-array standardisation;
#     sensitivity 0.98 and 0.995.  Spearman is reported alongside.
#   secondary (compendium-centred) rule: each gene is centred on its own compendium mean over
#     arrays before correlating (removes the shared mean-expression profile that inflates r);
#     a near-duplicate is flagged when the centred r exceeds the largest centred max-r seen in
#     the post-build null (EcoMAC arrays from series released after the M3D build; below).
#   null A (time-matched): EcoMAC arrays whose GEO series number is >= GSE20000 (released 2010
#     or later, after the Sept 2009 M3D build); they cannot be in M3D Build 6.  Reported by
#     platform (GPL199 Antisense vs GPL3154 Genome 2.0) so the null is platform-matched.
#   null B: within-compendium nearest neighbour excluding same series / same experiment.
# Identifier check: M3D carries no GEO accessions (chip and experiment names only), so GSE/GSM
#   cannot be compared directly; experimenter / author surnames are compared as weak evidence.
# Run: nice -n 19 systemd-run --user --scope -p MemoryMax=16G -p CPUQuota=400% \
#        Rscript analysis/ecomac_m3d_overlap.R
suppressMessages(library(stats))
out <- "results/ecomac_m3d_overlap"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
THR <- c(0.98, 0.99, 0.995)
say <- function(...) cat(sprintf(...), "\n", sep = "")

## ---- load EcoMAC ---------------------------------------------------------------------
E <- as.matrix(read.delim("data/ecomac/EcoMAC.tsv", header = FALSE))
gn <- read.delim("data/ecomac/gene_names.tsv", header = FALSE, stringsAsFactors = FALSE)
rownames(E) <- gn[[2]]
smp <- read.delim("data/ecomac/samples.tsv", check.names = FALSE, stringsAsFactors = FALSE)
plt <- read.csv("data/ecomac/array_platform.csv", stringsAsFactors = FALSE)
stopifnot(ncol(E) == nrow(smp), ncol(E) == nrow(plt), all(smp[[2]] == plt$sample))
colnames(E) <- smp[[2]]                       # GSM (where known) or other CEL id
meta <- data.frame(array = seq_len(ncol(E)), id = smp[[2]], series = smp[[3]], author = smp[[4]],
                   date = smp[[5]], strain = smp[[6]], platform = sub(".*\\((GPL[0-9]+)\\)", "\\1", plt$platform),
                   stringsAsFactors = FALSE)
meta$platform[!grepl("^GPL", meta$platform)] <- "unknown (no GEO record / Isalan)"
gse_num <- suppressWarnings(as.integer(sub("^GSE", "", ifelse(grepl("^GSE", meta$series), meta$series, NA))))
meta$post_build <- !is.na(gse_num) & gse_num >= 20000L

## ---- load M3D -------------------------------------------------------------------------
m3dir <- "data/E_coli_v4_Build_6"; base <- basename(m3dir)
rd <- function(f) { r <- read.delim(file.path(m3dir, f), check.names = FALSE); X <- as.matrix(r[, -1])
  rownames(X) <- vapply(strsplit(r[[1]], "_"), function(p) p[length(p) - 1L], ""); X }
M  <- rd(paste0(base, "_chips907probes4297.tab"))
Ma <- rd(paste0("avg_", base, "_exps466probes4297.tab"))
chipmeta <- read.delim(file.path(m3dir, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE, quote = "")
chipmeta <- chipmeta[match(colnames(M), chipmeta$chip_name), ]
stopifnot(!anyNA(chipmeta$chip_name))
chipmeta$experimenter <- sub(",.*", "", sub(".*experimenter:", "", ifelse(grepl("experimenter:", chipmeta$description), chipmeta$description, "")))

## ---- (c) genes and normalisation ------------------------------------------------------
sh <- intersect(rownames(E), rownames(M))
say("EcoMAC genes %d; M3D genes %d; shared b-numbers %d (%.1f%% of EcoMAC, %.1f%% of M3D)",
    nrow(E), nrow(M), length(sh), 100 * length(sh) / nrow(E), 100 * length(sh) / nrow(M))
say("Duplicated b-numbers: EcoMAC %d, M3D %d", sum(duplicated(rownames(E))), sum(duplicated(rownames(M))))
qs <- function(X) c(min = min(X), q25 = quantile(X, .25, names = FALSE), med = median(X), mean = mean(X),
                    q75 = quantile(X, .75, names = FALSE), max = max(X))
colmean_range <- function(X) range(colMeans(X)); colsd_range <- function(X) range(apply(X, 2, sd))
norm_tab <- rbind(EcoMAC = c(qs(E), cm_lo = colmean_range(E)[1], cm_hi = colmean_range(E)[2], csd_lo = colsd_range(E)[1], csd_hi = colsd_range(E)[2]),
                  M3D_907 = c(qs(M), cm_lo = colmean_range(M)[1], cm_hi = colmean_range(M)[2], csd_lo = colsd_range(M)[1], csd_hi = colsd_range(M)[2]))
write.csv(round(norm_tab, 4), file.path(out, "scale_summary.csv"))
print(round(norm_tab, 3))

## ---- (a) identifier / metadata comparison ---------------------------------------------
eco_sur <- tolower(sub("^.*[ .]([A-Za-z\\-]+)\\s*$", "\\1", sub(",.*", "", meta$author)))   # "K. Aggarwal" -> aggarwal
m3d_sur <- tolower(sub("\\s+.*$", "", chipmeta$experimenter))                              # "Wierzbowski J" -> wierzbowski
eco_authors <- sort(table(eco_sur[nzchar(eco_sur)]), decreasing = TRUE)
m3d_authors <- sort(table(m3d_sur[nzchar(m3d_sur)]), decreasing = TRUE)
common <- intersect(names(eco_authors), names(m3d_authors))
auth_tab <- data.frame(surname = common, ecomac_arrays = as.integer(eco_authors[common]), m3d_chips = as.integer(m3d_authors[common]))
auth_tab <- auth_tab[order(-auth_tab$ecomac_arrays), ]
write.csv(auth_tab, file.path(out, "shared_author_surnames.csv"), row.names = FALSE)
say("Surnames shared between EcoMAC 'Author' and M3D 'experimenter': %d (EcoMAC arrays by those authors %d of %d; M3D chips %d of %d)",
    length(common), sum(eco_sur %in% common), ncol(E), sum(m3d_sur %in% common), ncol(M))
say("M3D chip/experiment names containing GSE/GSM: %d; EcoMAC ids that are GSM: %d of %d",
    sum(grepl("GS[EM][0-9]", c(colnames(M), chipmeta$experiment_name, chipmeta$description))),
    sum(grepl("^GSM", meta$id)), nrow(meta))

## ---- (b) label-free similarity --------------------------------------------------------
stdz <- function(X) { X <- scale(X); X[is.na(X)] <- 0; X }           # per-array (column) standardisation
pear <- function(A, B) crossprod(stdz(A), stdz(B)) / (nrow(A) - 1)
rk   <- function(X) apply(X, 2, rank)
Es <- E[sh, ]; Ms <- M[sh, ]; Mas <- Ma[sh, ]
cen <- function(X) X - rowMeans(X)                                    # gene-centre within compendium
C_P  <- pear(Es, Ms)                                                  # 2198 x 907 Pearson
C_S  <- pear(rk(Es), rk(Ms))                                          # Spearman
C_Pc <- pear(cen(Es), cen(Ms))                                        # compendium-centred Pearson
C_Pa <- pear(Es, Mas)                                                 # vs the 466 averaged experiments
C_Pac<- pear(cen(Es), cen(Mas))
best <- function(Cm) { j <- max.col(Cm, ties.method = "first"); v <- Cm[cbind(seq_len(nrow(Cm)), j)]
  s <- apply(Cm, 1, function(x) sort(x, decreasing = TRUE)[2]); list(j = j, v = v, second = s) }
bP <- best(C_P); bS <- best(C_S); bPc <- best(C_Pc); bPa <- best(C_Pa); bPac <- best(C_Pac)

## within-compendium nearest neighbours (excluding self and same series / experiment)
nn_within <- function(X, grp) { C <- pear(X, X); diag(C) <- -Inf
  same <- outer(grp, grp, "=="); same[is.na(same)] <- FALSE; C[same] <- -Inf
  apply(C, 1, max) }
grp_e <- ifelse(is.na(meta$series) | meta$series %in% c("Not found", ""), paste0("single", seq_len(nrow(meta))), meta$series)
nnE <- nn_within(Es, grp_e)
nnM <- nn_within(Ms, chipmeta$experiment_name)
nnE_any <- { C <- pear(Es, Es); diag(C) <- -Inf; apply(C, 1, max) }   # incl. same series (replicates)

## flag rates
flag <- function(v, t) mean(v > t)
res <- data.frame(array = meta$array, id = meta$id, series = meta$series, platform = meta$platform, post_build = meta$post_build,
                  best_pearson = bP$v, second_pearson = bP$second, best_m3d_chip = colnames(M)[bP$j], best_m3d_experiment = chipmeta$experiment_name[bP$j],
                  best_m3d_experimenter = chipmeta$experimenter[bP$j],
                  best_spearman = bS$v, best_centred = bPc$v, best_centred_m3d_chip = colnames(M)[bPc$j],
                  best_pearson_avg466 = bPa$v, best_centred_avg466 = bPac$v,
                  nn_within_ecomac_other_series = nnE, nn_within_ecomac_any = nnE_any, stringsAsFactors = FALSE)
write.csv(res, file.path(out, "per_array_best_match.csv"), row.names = FALSE)

grpdefs <- list(all = rep(TRUE, nrow(meta)), `GPL199` = meta$platform == "GPL199", `GPL3154` = meta$platform == "GPL3154",
                `post_build (GSE>=20000)` = meta$post_build, `post_build GPL199` = meta$post_build & meta$platform == "GPL199",
                `post_build GPL3154` = meta$post_build & meta$platform == "GPL3154",
                `pre_build-eligible (GSE<20000)` = !is.na(gse_num) & gse_num < 20000L,
                `non-GEO / unknown series` = is.na(gse_num))
mk <- function(v, name) do.call(rbind, lapply(names(grpdefs), function(g) { i <- grpdefs[[g]]
  data.frame(measure = name, group = g, n = sum(i), median = median(v[i]), q05 = quantile(v[i], .05, names = FALSE), q95 = quantile(v[i], .95, names = FALSE),
             max = max(v[i]), n_gt_0.98 = sum(v[i] > .98), n_gt_0.99 = sum(v[i] > .99), n_gt_0.995 = sum(v[i] > .995),
             frac_gt_0.99 = mean(v[i] > .99)) }))
tab <- rbind(mk(bP$v, "EcoMAC->M3D907 max Pearson"), mk(bS$v, "EcoMAC->M3D907 max Spearman"),
             mk(bPa$v, "EcoMAC->M3D466avg max Pearson"), mk(bPc$v, "EcoMAC->M3D907 max centred Pearson"),
             mk(nnE, "EcoMAC->EcoMAC other series max Pearson"), mk(nnE_any, "EcoMAC->EcoMAC any (incl. same series)"))
tab <- cbind(tab[1:2], round(tab[-(1:2)], 4)); write.csv(tab, file.path(out, "similarity_summary.csv"), row.names = FALSE)
nm <- data.frame(measure = "M3D->M3D other experiment max Pearson", n = ncol(M), median = median(nnM), q05 = quantile(nnM, .05, names=FALSE), q95 = quantile(nnM, .95, names=FALSE), max = max(nnM),
                 n_gt_0.98 = sum(nnM > .98), n_gt_0.99 = sum(nnM > .99), n_gt_0.995 = sum(nnM > .995))
write.csv(round(nm[-1], 4), file.path(out, "m3d_within_nn_summary.csv"), row.names = FALSE)

## secondary rule: centred r > largest centred max-r in the post-build null
null_max_c <- max(bPc$v[meta$post_build]); null_q999 <- quantile(bPc$v[meta$post_build], .999, names = FALSE)
res$flag_centred_rule <- bPc$v > null_max_c
say("Centred-r null: post-build arrays n=%d, max=%.4f, 99.9th pct=%.4f; EcoMAC arrays exceeding the null max: %d (post-build %d)",
    sum(meta$post_build), null_max_c, null_q999, sum(res$flag_centred_rule), sum(res$flag_centred_rule & meta$post_build))
write.csv(res, file.path(out, "per_array_best_match.csv"), row.names = FALSE)

## per-series table: which GEO series have flagged arrays, and which M3D experiments they match
res$f99 <- res$best_pearson > .99
ser <- do.call(rbind, lapply(split(res, ifelse(is.na(res$series), "NA", res$series)), function(d)
  data.frame(series = d$series[1], n = nrow(d), platform = paste(unique(d$platform), collapse = ","), n_gt_0.98 = sum(d$best_pearson > .98), n_gt_0.99 = sum(d$f99),
             n_gt_0.995 = sum(d$best_pearson > .995), median_best = round(median(d$best_pearson), 4), n_centred_rule = sum(d$flag_centred_rule),
             median_best_centred = round(median(d$best_centred), 4),
             top_m3d_experiment = names(sort(table(d$best_m3d_experiment), decreasing = TRUE))[1])))
ser <- ser[order(-ser$n_gt_0.99, -ser$median_best), ]; write.csv(ser, file.path(out, "per_series_overlap.csv"), row.names = FALSE)

## reciprocal-best: mutual nearest neighbour across compendia
rb <- max.col(t(C_P), ties.method = "first")                          # best EcoMAC for each M3D chip
res$mutual <- rb[bP$j] == seq_len(nrow(res))
say("EcoMAC arrays that are mutual nearest neighbours of their best M3D chip: %d; with r>0.99: %d", sum(res$mutual), sum(res$mutual & res$f99))
write.csv(res, file.path(out, "per_array_best_match.csv"), row.names = FALSE)

## ---- POST-HOC additions (made after inspecting the distributions above; flagged as such) ---------
## The pre-stated Pearson thresholds cannot be reached across the two pipelines (max r < 0.995),
## so the overlap call below combines the pre-stated centred-r rule with a mutual-nearest-neighbour
## requirement.  Both are reported separately as well.
res$call <- res$flag_centred_rule & res$mutual
res$chip_for_call <- ifelse(res$call, res$best_m3d_chip, NA)
ncall <- sum(res$call)
chips_cov <- unique(na.omit(res$chip_for_call))
exps_cov <- unique(chipmeta$experiment_name[match(chips_cov, chipmeta$chip_name)])
rules <- data.frame(
  rule = c("Pearson > 0.98 (pre-stated sensitivity)", "Pearson > 0.99 (pre-stated primary)", "Pearson > 0.995 (pre-stated sensitivity)",
           "Spearman > 0.99", "centred Pearson > post-build null max (pre-stated secondary)", "mutual nearest neighbour (post hoc)",
           "centred rule AND mutual NN (post hoc overlap call)"),
  n_ecomac = c(sum(bP$v > .98), sum(bP$v > .99), sum(bP$v > .995), sum(bS$v > .99), sum(res$flag_centred_rule), sum(res$mutual), ncall))
rules$frac_of_2198 <- round(rules$n_ecomac / nrow(res), 4)
rules$false_positives_in_post_build_null <- c(sum(bP$v[meta$post_build] > .98), sum(bP$v[meta$post_build] > .99), sum(bP$v[meta$post_build] > .995),
  sum(bS$v[meta$post_build] > .99), sum(res$flag_centred_rule[meta$post_build]), sum(res$mutual[meta$post_build]), sum(res$call[meta$post_build]))
write.csv(rules, file.path(out, "overlap_rules.csv"), row.names = FALSE)
say("Overlap call (centred rule & mutual NN): %d EcoMAC arrays (%.1f%%); covers %d of %d M3D chips and %d of %d M3D experiments",
    ncall, 100 * ncall / nrow(res), length(chips_cov), ncol(M), length(exps_cov), length(unique(chipmeta$experiment_name)))
byplat <- table(platform = meta$platform, called = res$call); print(byplat); write.csv(as.data.frame.matrix(byplat), file.path(out, "call_by_platform.csv"))
by_era <- table(era = ifelse(meta$post_build, "GSE>=20000", ifelse(is.na(gse_num), "non-GEO/unknown", "GSE<20000")), called = res$call)
print(by_era); write.csv(as.data.frame.matrix(by_era), file.path(out, "call_by_era.csv"))
## bimodality of the centred correlation on GPL199: pre-build-eligible versus post-build
brk <- c(-1, .5, .6, .7, .75, .8, .85, .9, 1)
h <- rbind(`GPL199 GSE<20000` = table(cut(bPc$v[meta$platform == "GPL199" & !is.na(gse_num) & gse_num < 20000], brk)),
           `GPL199 GSE>=20000` = table(cut(bPc$v[meta$platform == "GPL199" & meta$post_build], brk)),
           `GPL3154` = table(cut(bPc$v[meta$platform == "GPL3154"], brk)))
write.csv(h, file.path(out, "centred_r_histogram.csv")); print(h)
## per-series one-to-one structure
res$call_series <- ifelse(is.na(res$series), "NA", res$series)
sc <- do.call(rbind, lapply(split(res, res$call_series), function(d) data.frame(series = d$call_series[1], n = nrow(d), n_called = sum(d$call),
      distinct_m3d_chips_among_called = length(unique(na.omit(d$chip_for_call))), frac_called = round(mean(d$call), 3))))
sc <- sc[order(-sc$n_called), ]; write.csv(sc, file.path(out, "call_by_series.csv"), row.names = FALSE)
cat("series with >=1 called array: ", sum(sc$n_called > 0), " of ", nrow(sc), "; fully (>=90%) called: ", sum(sc$frac_called >= .9), "\n", sep = "")
## normalisation check: do all arrays share one distribution (joint quantile normalisation)?
qn <- function(X) { S <- apply(X, 2, sort); max(apply(S, 1, function(r) diff(range(r)))) }
say("Max over rank positions of the range across arrays of sorted values: EcoMAC %.2e (0 = identical distributions, quantile-normalised); M3D 907 chips %.3f",
    qn(E), qn(M))
say("Per-platform EcoMAC array mean (range of means): %s", paste(sapply(split(colMeans(E), meta$platform), function(x) sprintf("%.4f-%.4f", min(x), max(x))), collapse = "; "))
## context for the centred correlation: how similar are true replicates inside one compendium?
Cm <- pear(cen(Ms), cen(Ms)); same_exp <- outer(chipmeta$experiment_name, chipmeta$experiment_name, "==") & upper.tri(Cm)
Pm <- pear(Ms, Ms)
say("M3D replicate chips (same experiment, %d pairs): median centred r %.3f (IQR %.3f-%.3f), median raw r %.4f",
    sum(same_exp), median(Cm[same_exp]), quantile(Cm[same_exp], .25), quantile(Cm[same_exp], .75), median(Pm[same_exp]))
rm(Cm, Pm)
write.csv(res, file.path(out, "per_array_best_match.csv"), row.names = FALSE)

## figure-free text summary
sink(file.path(out, "summary.txt"))
say("shared genes %d / EcoMAC %d / M3D %d", length(sh), nrow(E), nrow(M))
print(tab[tab$measure %in% c("EcoMAC->M3D907 max Pearson", "EcoMAC->M3D907 max centred Pearson"), ], row.names = FALSE)
print(nm, row.names = FALSE)
print(head(ser, 25), row.names = FALSE)
sink()
cat(readLines(file.path(out, "summary.txt")), sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
