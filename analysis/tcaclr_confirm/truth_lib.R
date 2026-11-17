# Shared truth-building code for docs/plans/tcaclr_confirm_plan.md (label-free build stage).
# Functions filt(), merge_sets() and structure_of() are verbatim copies of
# analysis/complex_confirm/build_yeast_complex_truth.R (net4 truth); structure_of() reads the global `paralog_sets`.
# make_variants() reproduces the variant list of that builder plus the plan's size classes (mid split 7-10 / 11-19);
# the only new arguments are the dataset-specific ribosome predicate and the ribosome-subunit complex finders.
# Builds and counts structure only; no expression value or score is read.

filt <- function(sets, u, lo = 3L, hi = 200L) { s <- lapply(sets, function(g) intersect(g, u)); s[lengths(s) >= lo & lengths(s) <= hi] }
merge_sets <- function(s) {
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]
  key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
    length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]
  if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
  out
}
structure_of <- function(s, xgroups = list()) {
  if (!length(s)) return(list(w = numeric(0), n_pos = 0, n_pairs = 0, n_genes = 0, comp = 0, n_excl = 0))
  g <- sort(unique(unlist(s))); n <- length(g)
  cnt <- matrix(0L, n, n)
  for (k in seq_along(s)) { m <- match(s[[k]], g); cnt[m, m] <- cnt[m, m] + 1L }
  excl <- matrix(FALSE, n, n)
  for (ps in paralog_sets) { m <- stats::na.omit(match(ps, g)); if (length(m) >= 2) excl[m, m] <- TRUE }
  for (xg in xgroups) { a <- match(xg[[1]], g); b <- match(xg[[2]], g); excl[a, b] <- TRUE; excl[b, a] <- TRUE }  # cross-subunit
  diag(excl) <- FALSE; cnt[excl] <- 0L
  ut <- upper.tri(cnt); pos <- cnt > 0 & ut
  w <- vapply(seq_along(s), function(k) { m <- match(s[[k]], g); sub <- cnt[m, m]; x <- sub[upper.tri(sub)]; sum(1 / x[x > 0]) }, 0)
  adj <- crossprod(vapply(s, function(x) g %in% x, logical(n))) > 0
  cl <- seq_along(s); repeat { new <- vapply(seq_along(s), function(k) min(cl[adj[k, ]]), 1L); if (all(new == cl)) break; cl <- new }
  list(w = stats::setNames(w, names(s)), n_pos = sum(pos), n_pairs = n * (n - 1) / 2 - sum(excl & ut), n_genes = n,
       comp = length(unique(cl)), n_excl = sum(excl & ut))
}

# gen: named list of complexes -> genes (before universe filter); u: universe; ribo_frac(g): r-protein fraction;
# sub_pairs: named list of c(complex-id-regex-a, complex-id-regex-b) identifying the cytoplasmic and mitochondrial
# large/small ribosomal subunit complexes among the MERGED complexes (by source id or name; label-free);
# perturbed: character vector of truth-universe genes that are perturbed in the panel (counts only).
make_variants <- function(gen, u, ribo_frac, sub_pairs, perturbed, dataset, ribo_cut = 0.3) {
  s3 <- filt(gen, u); main <- merge_sets(s3)
  isrp <- vapply(main, function(g) ribo_frac(g) >= ribo_cut, TRUE)
  find1 <- function(rx) { h <- grep(rx, names(main), perl = TRUE); stopifnot(length(h) == 1L); names(main)[h] }
  sp <- lapply(sub_pairs, function(p) c(find1(p[1]), find1(p[2])))
  xgroups <- lapply(sp, function(p) list(main[[p[1]]], main[[p[2]]]))
  joined <- main[setdiff(names(main), unlist(sp))]
  for (p in sp) joined[[paste(p, collapse = "+")]] <- sort(unique(c(main[[p[1]]], main[[p[2]]])))
  joined <- joined[order(names(joined))]
  variants <- list(main = main, main_xsub_excl = main, main_ribo_joined = joined, noribo = main[!isrp], ribo_list = main[isrp], nomerge = s3,
                   small_3to6 = main[lengths(main) <= 6],
                   mid_7to19 = main[lengths(main) >= 7 & lengths(main) <= 19],
                   mid_7to10 = main[lengths(main) >= 7 & lengths(main) <= 10],
                   mid_11to19 = main[lengths(main) >= 11 & lengths(main) <= 19],
                   large_ge20 = main[lengths(main) >= 20])
  rows <- list(); counts <- list(); sizes <- NULL
  for (v in names(variants)) {
    s <- variants[[v]]
    for (k in names(s)) rows[[length(rows) + 1]] <- data.frame(variant = v, dataset = dataset, complex = k, gene = s[[k]])
    st <- structure_of(s, if (v == "main_xsub_excl") xgroups else list()); w <- st$w / max(st$n_pos, 1)
    tg <- sort(unique(unlist(s)))
    counts[[length(counts) + 1]] <- data.frame(dataset = dataset, variant = v, n_universe = length(u), n_complexes = length(s),
      n_genes = st$n_genes, n_genes_perturbed = sum(tg %in% perturbed), n_pairs = st$n_pairs, n_paralog_or_xsub_pairs_excluded = st$n_excl,
      n_pos = st$n_pos, base_rate = st$n_pos / max(st$n_pairs, 1), n_complexes_pos = sum(lengths(s) >= 2), n_components = st$comp,
      n_eff = if (length(w)) 1 / sum(w^2) else 0, largest_share = if (length(w)) max(w) else 0,
      ribo_dominated_share = if (v == "main") sum(w[isrp]) else NA_real_)
    if (v == "main") sizes <- data.frame(dataset = dataset, complex = names(s), n_genes = lengths(s), pos_weight = st$w, share = w,
      rp_frac = vapply(s, ribo_frac, 0), ribo_dominated = isrp, row.names = NULL)
  }
  list(main = main, variants = variants, members = do.call(rbind, rows), counts = do.call(rbind, counts),
       sizes = sizes[order(-sizes$share), ], sub_pairs = sp, xgroups = xgroups, isrp = isrp, s3 = s3)
}
excluded_pairs <- function(paralog_sets, sp, main) {
  pp <- do.call(rbind, lapply(paralog_sets, function(ps) if (length(ps) >= 2) { cb <- utils::combn(sort(ps), 2)
    data.frame(scope = "all", gene_a = cb[1, ], gene_b = cb[2, ], reason = "paralog_set") }))
  xs <- do.call(rbind, lapply(names(sp), function(nm) { a <- main[[sp[[nm]][1]]]; b <- main[[sp[[nm]][2]]]
    e <- expand.grid(gene_a = a, gene_b = b, stringsAsFactors = FALSE); e <- e[e$gene_a != e$gene_b, ]
    data.frame(scope = "main_xsub_excl", gene_a = pmin(e$gene_a, e$gene_b), gene_b = pmax(e$gene_a, e$gene_b), reason = paste0("cross_subunit_", nm)) }))
  ex <- rbind(pp, xs); ex[order(ex$scope, ex$reason, ex$gene_a, ex$gene_b), ]
}
