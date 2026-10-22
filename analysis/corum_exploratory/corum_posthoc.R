#!/usr/bin/env Rscript
# EXPLORATORY. Post-hoc summaries from corum_diag.R outputs (no truths beyond those already scored):
# dilution slopes, per-complex loss vs size and member degree, local vs pooled.
setwd("/home/boris/Projects/ClrAttention")
SP <- "/tmp/claude-1000/-home-boris-Projects-ClrAttention/b33914c9-0cfd-485d-8dd4-cf8a2a22b632/scratchpad"
OUT <- file.path(SP, "corum"); logf <- file.path(OUT, "posthoc.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
source("analysis/regulondb.R")
RUN <- "results/baselines_k562ess_truthgenes_genie3"
genes <- utils::read.delim(file.path(RUN, "expression.tsv"), colClasses = c("character", rep("NULL", 707)))[[1]]
G <- length(genes)
core <- readRDS(file.path(OUT, "core.rds")); C <- core$C; A <- core$A; rm(core)
cc <- utils::read.csv(file.path(OUT, "corum_complexes.csv"), stringsAsFactors = FALSE)
rmt <- utils::read.csv(file.path(RUN, "regulon_members.csv"), stringsAsFactors = FALSE)
x <- rmt[rmt$truth == "CORUM", ]; cm <- lapply(strsplit(x$members, ";", fixed = TRUE), function(g) match(g, genes)); names(cm) <- x$regulon
# member "degree": CLR partners within 1 unit of the row maximum (near-equal plateau) and number with CLR >= 3
plateau <- vapply(seq_len(G), function(i) { r <- C[i, ]; r[i] <- -Inf; sum(r >= max(r) - 1) }, 0)
strong <- rowSums(C >= 3)
perp <- apply(A, 1, function(a) { a <- a[a > 0]; exp(-sum(a * log(a))) })
say(sprintf("softmax row perplexity: median %.2f (should be 10)", stats::median(perp)))
lo <- readRDS(file.path(OUT, "local.rds")); MA <- lo$MA; RA <- lo$RA
mac <- data.frame(complex = MA$clr$regulon, size = MA$clr$size, clr = MA$clr$ap, t1 = MA$att_t1$ap, t9 = MA$att_t9$ap)
for (v in setdiff(names(MA), c("clr", "att_t1", "att_t9", "att_t2"))) mac[[v]] <- MA[[v]]$ap
mac$group <- cc$group[match(mac$complex, cc$complex)]
mac$plateau <- vapply(cm[mac$complex], function(m) stats::median(plateau[m]), 0)
mac$strong <- vapply(cm[mac$complex], function(m) stats::median(strong[m]), 0)
ok <- !is.na(mac$clr) & !is.na(mac$t9)
say(sprintf("macro (per-complex AP, %d complexes with >= 10 pairs): mean CLR %.3f, t1 %.3f, t9 %.3f; CLR > t9 in %.0f%%",
            sum(ok), mean(mac$clr[ok]), mean(mac$t1[ok]), mean(mac$t9[ok]), 100 * mean(mac$clr[ok] > mac$t9[ok])))
for (v in c("a_clr_on_logatt_t9", "a_clr_on_att_t9", "d_density_t9", "b_coarse_r06", "b0_boost_on_att_t9"))
  if (!is.null(mac[[v]])) say(sprintf("   %-20s mean %.3f; beats t9 in %.0f%%, beats CLR in %.0f%%", v, mean(mac[[v]][ok]),
                                      100 * mean(mac[[v]][ok] > mac$t9[ok]), 100 * mean(mac[[v]][ok] > mac$clr[ok])))
d <- mac$clr - mac$t9
say(sprintf("per-complex loss (CLR - t9) vs size: Spearman %.2f (p %.2g); vs member plateau %.2f (p %.2g); vs members with CLR >= 3: %.2f (p %.2g)",
            stats::cor(d[ok], mac$size[ok], method = "spearman"), stats::cor.test(d[ok], mac$size[ok], method = "spearman", exact = FALSE)$p.value,
            stats::cor(d[ok], mac$plateau[ok], method = "spearman"), stats::cor.test(d[ok], mac$plateau[ok], method = "spearman", exact = FALSE)$p.value,
            stats::cor(d[ok], mac$strong[ok], method = "spearman"), stats::cor.test(d[ok], mac$strong[ok], method = "spearman", exact = FALSE)$p.value))
nr <- ok & mac$group == "other"
say(sprintf("   non-RP complexes only (%d): Spearman with size %.2f, with plateau %.2f; mean loss %.4f, median %.4f, CLR > t9 in %.0f%%",
            sum(nr), stats::cor(d[nr], mac$size[nr], method = "spearman"), stats::cor(d[nr], mac$plateau[nr], method = "spearman"),
            mean(d[nr]), stats::median(d[nr]), 100 * mean(d[nr] > 0)))
for (b in list(c(5, 9), c(10, 24), c(25, 49), c(50, 200))) { s <- ok & mac$size >= b[1] & mac$size <= b[2]
  say(sprintf("   size %3d-%3d: %3d complexes, mean AP CLR %.3f t1 %.3f t9 %.3f logCLR %.3f", b[1], b[2], sum(s),
              mean(mac$clr[s]), mean(mac$t1[s]), mean(mac$t9[s]),
              if (is.null(mac$a_clr_on_logatt_t9)) NA else mean(mac$a_clr_on_logatt_t9[s]))) }
utils::write.csv(mac, file.path(OUT, "corum_macro_ap.csv"), row.names = FALSE)
# row-local AP (no cross-row comparison) vs pooled
ra <- do.call(rbind, lapply(names(RA), function(v) data.frame(variant = v, RA[[v]])))
agg <- stats::aggregate(ap ~ variant + group, ra, mean)
say("row-local AP (member rows; complex partners vs genes sharing no complex), mean over member rows:")
for (v in unique(agg$variant)) say(sprintf("   %-22s RP-dominated %.3f   other %.3f", v, agg$ap[agg$variant == v & agg$group == "RP-dominated"],
                                           agg$ap[agg$variant == v & agg$group == "other"]))
# dilution slope: log within-complex pair score vs log complex size (median per complex)
wpc <- do.call(rbind, lapply(seq_along(cm), function(k) { m <- cm[[k]]; if (length(m) < 3) return(NULL); ut <- upper.tri(diag(length(m)))
  a <- (A[m, m] + t(A[m, m]))[ut] / 2
  data.frame(size = length(m), clr = stats::median(C[m, m][ut]), a1 = mean(a), a1_nz = mean(a > 0)) }))
f1 <- stats::lm(log(a1 + 1e-12) ~ log(size), wpc[wpc$a1 > 0, ]); f0 <- stats::lm(clr ~ log(size), wpc)
say(sprintf("within-complex pairs: log(mean softmax weight) ~ log(size) slope %.2f (pure 1/n dilution = -1); CLR median ~ log(size) slope %.2f z-units per e-fold",
            stats::coef(f1)[2], stats::coef(f0)[2]))
say(sprintf("fraction of within-complex pairs with any softmax weight (j in i's top 50 or i in j's): size <= 9 %.2f, 10-24 %.2f, 25-49 %.2f, >= 50 %.2f",
            mean(wpc$a1_nz[wpc$size <= 9]), mean(wpc$a1_nz[wpc$size >= 10 & wpc$size <= 24]),
            mean(wpc$a1_nz[wpc$size >= 25 & wpc$size <= 49]), mean(wpc$a1_nz[wpc$size >= 50])))
