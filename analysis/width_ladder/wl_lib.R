# Width ladder (docs/plans/width_ladder_plan.md): shared helpers. EXPLORATORY, DESCRIPTIVE; spent truths. Source from the repo root.
# Nothing here defines an operator, a walk, a read-out or a metric: gl_make_P / gl_walk / gl_zmat (analysis/tcaclr_depth/grid_lib.R), the edge machinery of
# analysis/tcaclr2/edge_depth.R and analysis/tcaclr_depth/adaptive/adaptive_lib.R, ap_fast, sb_per_tf_median_auroc, rz_hub_stratify, the truth builders of
# sb_load_ecoli / build_truths_ecoli / ideas47_lib.R / grid_truths.R are the existing ones. New: loops over widths, multi-scale accumulation, tables.
suppressPackageStartupMessages({ source("analysis/tcaclr2/edge_depth.R"); source("analysis/tcaclr2/tsb_lib.R") })   # grid_lib, grid_truths (tc_lib), adaptive_lib, ED, ed_t2
source("analysis/robust_z/robust_z_lib.R")

WL <- list(widths = c(2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128), ms3 = c(4, 6, 8, 12, 16, 24, 32, 48), cap_mult = 4L, tmax = 256L, t_fixed = 9L,
           ref = list(keff = 10, cap = 50L), out = "results/width_ladder")
WLD <- list(   # dataset -> kind, edge file of the registered tcaCLR2 (t2 of the comparator), names in the reference tables
  m3d      = list(kind = "ecoli", edge = "results/tcaclr2/edge_m3d_tf_regulons.tsv", ref_ds = "m3d"),
  ecomac   = list(kind = "ecoli", edge = file.path(TSB$out, "edge_ecomac.tsv"), ref_ds = "ecomac"),
  k562     = list(kind = "k562", edge = "results/tcaclr2/edge_k562_corum.tsv", ref_ds = "k562"),
  yeast    = list(kind = "plain", gt = "net4_cp_yeast", edge = "results/tcaclr2/edge_net4_cp_yeast.tsv", ref_ds = "dream5_yeast", truth = "complexes_main"),
  hepg2    = list(kind = "plain", gt = "armA_cp_human", edge = "results/tcaclr2/edge_armA_cp_human.tsv", ref_ds = "hepg2", truth = "complexes_main"),
  kemmeren = list(kind = "plain", gt = "armB_cyc2008", edge = "results/tcaclr2/edge_armB_cyc2008.tsv", ref_ds = "kemmeren", truth = "complexes_main"))
K562_TRUTHS <- c(DoRothEA_AB = "DoRothEA_AB", DoRothEA_ABC = "DoRothEA_ABC", CollecTRI = "CollecTRI", CORUM = "CORUM_all", KD_response = "KD_response")   # name in tables = name in ideas47 bench

wl_cap <- function(keff, G) as.integer(min(WL$cap_mult * keff, G - 1L))
wl_ops <- function(G, widths = WL$widths, ref = TRUE) {
  o <- lapply(widths, function(k) list(id = sprintf("W%g", k), keff = k, cap = wl_cap(k, G), ref = FALSE))
  if (ref) o <- c(list(list(id = "ref", keff = WL$ref$keff, cap = WL$ref$cap, ref = TRUE)), o); o }
# depth from the edge: t = max(1, round(ln 0.05 / ln edge)); no decay (edge >= 1 - 1e-6) or t* > tmax -> tmax, flagged (plan, fixed choice 1)
wl_t_from_edge <- function(edge) {
  if (!is.finite(edge) || edge >= 1 - 1e-6) return(c(t = WL$tmax, capped = 1L, tstar = NA_real_))
  if (edge <= 0) return(c(t = 1L, capped = 0L, tstar = 0))
  ts <- log(ED$eps) / log(edge); if (ts > WL$tmax) c(t = WL$tmax, capped = 1L, tstar = ts) else c(t = ed_t2(edge), capped = 0L, tstar = ts) }

## ---- matrices (X only; the edge job needs no truth) ---------------------------------------------------------------------------
wl_matrix <- function(ds) {
  s <- WLD[[ds]]
  if (s$kind == "ecoli") { source("analysis/scoreboard/sb_lib.R"); sb_source(); return(sb_load_ecoli(ds)$X) }
  if (s$kind == "k562") return(gl_truth_matrix("k562_corum"))
  gl_truth_matrix(s$gt)
}

## ---- edge (permutation) stage --------------------------------------------------------------------------------------------------
# One fit per permuted matrix (seeds ED seed + 7000 + p, as the registered edge); every operator of `ops` is read from that fit: edge = 2nd eigenvalue of the
# symmetrised lazy operator (ad_Psym, ad_eigs with ED$kvals). All clr work inside the tasks (forked children of a process that has fitted hang).
wl_edges <- function(X, ops, perms, workers, threads, say = function(...) invisible(NULL)) {
  t0 <- Sys.time()
  one <- function(p) {
    f <- ad_fit(ad_permute_genes(X, ED$seed + 7000L + p), threads); on.exit(f$release())
    m <- t(vapply(ops, function(o) { v <- ad_eigs(ad_Psym(ad_P_from_fit(f, o$keff, cap = o$cap, clipS = TRUE)), ED$kvals)$values; c(top = v[1L], edge = v[2L]) }, c(top = 0, edge = 0)))
    list(m = m, rss = gl_vmhwm_gb()) }
  all <- ad_lapply(seq_len(perms), one, workers)
  bad <- vapply(all, function(x) !is.list(x), NA); if (any(bad)) stop("permutation task failed: ", all[bad][[1]])
  say(sprintf("%d permutations x %d operators: %.0f s (worker peak RSS up to %.1f GB)", perms, length(ops), as.numeric(difftime(Sys.time(), t0, units = "secs")), max(vapply(all, `[[`, 0, "rss"))))
  per <- do.call(rbind, lapply(seq_along(all), function(p) data.frame(op = vapply(ops, `[[`, "", "id"), perm = p, top = all[[p]]$m[, "top"], edge = all[[p]]$m[, "edge"])))
  tab <- do.call(rbind, lapply(seq_along(ops), function(k) { o <- ops[[k]]; e <- per$edge[per$op == o$id]; tt <- wl_t_from_edge(mean(e))
    data.frame(op = o$id, keff = o$keff, cap = o$cap, edge = mean(e), edge_min = min(e), edge_max = max(e), tstar = tt[["tstar"]], t = as.integer(tt[["t"]]), capped = as.integer(tt[["capped"]]),
               top_eig_min = min(per$top[per$op == o$id]), perms = perms, stringsAsFactors = FALSE) }))
  list(tab = tab, per = per)
}

## ---- datasets and truths -------------------------------------------------------------------------------------------------------
wl_pidx <- function(tr) { sw <- tr$i > tr$j; if (any(sw)) { x <- tr$i[sw]; tr$i[sw] <- tr$j[sw]; tr$j[sw] <- x }; tr$pidx <- tc_pair_index(tr$i, tr$j); tr }
wl_pooled_from <- function(name, i, j, label) wl_pidx(list(name = name, i = as.integer(i), j = as.integer(j), label = as.integer(label)))
# synthetic truths (smoke tests only): planted random modules / random TF target sets, no real label read
wl_synth_pooled <- function(name, G, seed, nmod = 40L) {
  set.seed(seed); mods <- lapply(seq_len(nmod), function(k) sort(sample(G, sample(4:25, 1L)))); gin <- sort(unique(unlist(mods)))
  M <- vapply(mods, function(m) gin %in% m, logical(length(gin))); co <- (M %*% t(M)) > 0; ut <- which(upper.tri(co), arr.ind = TRUE)
  wl_pooled_from(name, gin[ut[, 1]], gin[ut[, 2]], co[ut]) }
wl_synth_edge <- function(G, seed, ntf = 40L) {
  set.seed(seed); tfs <- sort(sample(G, ntf)); tg <- lapply(tfs, function(t) sort(sample(setdiff(seq_len(G), t), sample(8:60, 1L)))); names(tg) <- as.character(tfs)
  cand <- unique(do.call(rbind, lapply(tfs, function(t) { g <- setdiff(seq_len(G), t); cbind(pmin(t, g), pmax(t, g)) })))
  pos <- unique(do.call(rbind, lapply(seq_along(tfs), function(k) cbind(pmin(tfs[k], tg[[k]]), pmax(tfs[k], tg[[k]])))))
  lab <- as.integer(paste(cand[, 1], cand[, 2]) %in% paste(pos[, 1], pos[, 2])); list(i = cand[, 1], j = cand[, 2], label = lab, targets = tg, tfs = tfs) }

wl_load <- function(ds, synthetic = FALSE) {
  s <- WLD[[ds]]; D <- list(ds = ds, kind = s$kind, hub = s$kind != "plain", synthetic = synthetic, pooled = list(), edge = list())
  if (s$kind == "ecoli") {
    source("analysis/scoreboard/sb_lib.R"); sb_source(); source("analysis/scoreboard/sb_cx_lib.R"); sb_import_truth_defs()
    E <- sb_load_ecoli(ds); D$X <- E$X; D$genes <- E$genes; G <- length(D$genes)
    if (synthetic) { for (k in c("SC", "all", "complexes")) D$pooled[[k]] <- wl_synth_pooled(k, G, 1L + match(k, c("SC", "all", "complexes"))); D$edge <- list(SC = wl_synth_edge(G, 11L), all = wl_synth_edge(G, 12L)) }
    else { TR <- build_truths_ecoli(ds, D$genes, NULL, only = c("complexes", "regulons_SC", "regulons_all")); E$bench$complexes <- TR$bench$complexes
      for (k in c("SC", "all", "complexes")) { b <- E$bench[[k]]; D$pooled[[k]] <- wl_pooled_from(k, b$pairs$i, b$pairs$j, b$pairs$label) }; D$edge <- E$U }
  } else if (s$kind == "k562") {
    ex <- parse("analysis/explore_ideas47/ideas47_lib.R", keep.source = FALSE)
    key <- function(pat) which(vapply(ex, function(e) grepl(pat, paste(deparse(e)[1], collapse = "")), NA))[1]
    last <- if (synthetic) key("^lib_truths <- function") - 1L else key("^rm\\(TT\\)")
    opt <<- list(out = file.path(tempdir(), "wl_k562"), quick = 0L, jack = 1L, threads = 1L, workers = 1L)
    for (e in ex[seq_len(last)]) eval(e, globalenv())
    D$X <- X; D$genes <- genes; G <- length(genes); stopifnot(G == 5671L)
    if (synthetic) { for (k in names(K562_TRUTHS)) D$pooled[[k]] <- wl_synth_pooled(k, G, 20L + match(k, names(K562_TRUTHS))) }
    else { abc <- dorothea_sets(genes, c("A", "B", "C")); bench$DoRothEA_ABC <- list(reg = abc, pairs = comembership_pairs(abc, none, genes))
      for (k in names(K562_TRUTHS)) { b <- bench[[K562_TRUTHS[[k]]]]; D$pooled[[k]] <- wl_pooled_from(k, b$pairs$i, b$pairs$j, b$pairs$label) } }
  } else {
    X <- gl_truth_matrix(s$gt); D$X <- X; D$genes <- rownames(X); G <- nrow(X)
    if (synthetic) D$pooled[[s$truth]] <- wl_synth_pooled(s$truth, G, 31L)
    else { tr <- gl_truths(s$gt, D$genes)$main; D$pooled[[s$truth]] <- wl_pooled_from(s$truth, tr$i, tr$j, tr$label) }
  }
  D$G <- G; D$U <- gl_ut(G); D }

## ---- metrics ------------------------------------------------------------------------------------------------------------------
wl_unpack <- function(v, D) { S <- matrix(0, D$G, D$G); S[D$U$ut] <- v; S[D$U$utT] <- v; S }
# metric rows of one packed score vector v (upper triangle, gl_ut order). Same arithmetic as score_scoreboard.R / direct_edge_diag/03_diag.R.
wl_eval <- function(v, D) {
  v[!is.finite(v)] <- 0; rows <- list(); add <- function(truth, scope, metric, val) rows[[length(rows) + 1L]] <<- data.frame(truth = truth, scope = scope, metric = metric, value = val, stringsAsFactors = FALSE)
  for (tr in D$pooled) add(tr$name, "raw", "pooled_aupr_over_base", ap_fast(v[tr$pidx], tr$label) / mean(tr$label))
  if (D$hub || length(D$edge)) S <- wl_unpack(v, D)
  if (D$hub) { SH <- rz_hub_stratify(S, D$H)
    for (tr in D$pooled) add(tr$name, "hub_att", "pooled_aupr_over_base", ap_fast(SH[cbind(tr$i, tr$j)], tr$label) / mean(tr$label)) }
  for (nm in names(D$edge)) { u <- D$edge[[nm]]
    add(nm, "raw", "tf_aupr_ratio", ap_fast(S[cbind(u$i, u$j)], u$label) / mean(u$label))
    add(nm, "raw", "tf_auroc_median", sb_per_tf_median_auroc(S, u, D$G))
    add(nm, "hub_att", "tf_aupr_ratio", ap_fast(SH[cbind(u$i, u$j)], u$label) / mean(u$label)) }
  do.call(rbind, rows) }

## ---- scores of one operator ----------------------------------------------------------------------------------------------------
wl_tca_vec <- function(Pt, U) { Z <- gl_zmat(Pt); (pmax(Z[U$ut], 0) + pmax(Z[U$utT], 0)) / sqrt(2) }   # clipped Stouffer read-out, packed (= sk_a3r upper triangle)
wl_read_t2 <- function(ds) { e <- utils::read.delim(WLD[[ds]]$edge, stringsAsFactors = FALSE); stopifnot(nrow(e) == 1L, e$keff == 10, e$cap == 50); e }
wl_save <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); gl_wcsv(x, path) }
