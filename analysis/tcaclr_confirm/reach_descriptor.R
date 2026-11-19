#!/usr/bin/env Rscript
# Registered reachability descriptor (plan section 3, mandatory descriptors) for one arm: per complex size class, the reachable-entry count and the
# in-module share of reachable entries for member rows of log P^9 (P = 0.5 I + 0.5 A, t = 9; A as in make_scores.R). STRUCTURE ONLY: the expression
# matrix and the complex membership are read; no score is compared with any truth and no AUPR is computed. A class is flagged when its median
# reachable count is under 1.5 times its median size.
#   Rscript analysis/tcaclr_confirm/reach_descriptor.R --expr FILE --members FILE --out PREFIX [--threads 8]
suppressPackageStartupMessages(source("analysis/tcaclr_confirm/tc_lib.R"))
opt <- tc_options(commandArgs(TRUE), list(expr = "", members = "", out = "", threads = 8L))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
ex <- utils::read.delim(opt$expr, check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; rm(ex); G <- nrow(X)
mem <- utils::read.delim(opt$members, stringsAsFactors = FALSE); main <- mem[mem$variant == "main", ]; ribo <- unique(mem$complex[mem$variant == "ribo_list"])
say(sprintf("%d genes x %d profiles; %d main complexes", G, ncol(X), length(unique(main$complex))))
a <- TC$alpha
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)$calibrate(method = "normal", combine = "stouffer")
f$build_operator(alpha = a, softmax_keff = TC$keff, softmax_cap = as.integer(TC$cap))
P <- Matrix::Diagonal(G, 1 - a) + a * Matrix::Matrix(f$operator, sparse = TRUE); f$release()
Pt <- as.matrix(P); for (s in seq_len(TC$t - 1L)) Pt <- as.matrix(Pt %*% P); rm(P)
reach <- Pt > 0 & is.finite(log(Pt)); diag(reach) <- FALSE; rm(Pt); nreach <- rowSums(reach)
say(sprintf("unreachable share of all pairs %.4f; rows with fewer than 2 reachable: %d; median reachable per row %d of %d",
            1 - sum(nreach) / (G * (G - 1)), sum(nreach < 2), as.integer(median(nreach)), G - 1))
rows <- lapply(split(main$orf, main$complex), function(g) { m <- match(g, rownames(X)); m <- m[!is.na(m)]; n <- length(m)
  inm <- rowSums(reach[m, m, drop = FALSE]); data.frame(size = n, reachable_offdiag = mean(nreach[m]), inmodule_share = mean(inm / pmax(nreach[m], 1)), module_size_minus1 = n - 1) })
d <- do.call(rbind, Map(function(x, nm) cbind(complex = nm, x), rows, names(rows))); rownames(d) <- NULL
d$class <- ifelse(d$complex %in% ribo, "ribosomal", ifelse(d$size <= 6, "3-6", ifelse(d$size <= 19, "7-19", "20+")))
s <- do.call(rbind, lapply(split(d, d$class), function(x) data.frame(class = x$class[1], complexes = nrow(x), median_size = median(x$size),
  median_reachable = median(x$reachable_offdiag), median_inmodule_share = median(x$inmodule_share), min_inmodule_share = min(x$inmodule_share),
  flagged_under_1.5x_size = median(x$reachable_offdiag) < 1.5 * median(x$size), complexes_closed = sum(x$reachable_offdiag < 1.5 * x$size))))
print(format(s, digits = 3), row.names = FALSE)
write.csv(d, paste0(opt$out, "_by_complex.csv"), row.names = FALSE); write.csv(s, paste0(opt$out, "_by_class.csv"), row.names = FALSE)
