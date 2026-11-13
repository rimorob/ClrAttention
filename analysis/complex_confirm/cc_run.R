#!/usr/bin/env Rscript
# Driver of the exploratory estimation of docs/plans/complex_confirm_design.md (rev 2, FIXED 2026-11-13): a3r minus attention (and CLR) on
# yeast protein complexes (EBI Complex Portal 559292) in the DREAM5 network 4 compendium. One run; every decision is in the design.
#
# Modes (repo root):
#   gate_a        Ga part 1 (library identity) + collection of the K562 and M3D gate tables (gate_a_k562.R, gate_a_m3d.R run first)
#   refit_timing  one label-free group refit (rep 1): wall time, peak memory, SHA-256 against the DREAM5 manifest (fresh process)
#   labelfree     guard (no prior a3r score on net4), operator, P^9, caches with SHA-256, Ga parts 4 and 4b, Gb, decoy diagnostic
#   truth         primary and all secondary scores, 500 gene and 500 complex delete-half jackknife replicates
#   expjack       series-level (59) and group-level (321) refit jackknives, 100 replicates each, hash-checked (--lane/--nlanes, --collect 1)
#   summarize     writes estimate.csv ONCE and estimate_final.csv, tables.md
# --synthetic 1 --quick N: smoke test on a gene subset with synthetic truths only (never production).
suppressPackageStartupMessages(source("analysis/complex_confirm/cc_lib.R"))
opt <- cc_options(commandArgs(TRUE), list(mode = "", root = CC$out, quick = 0L, synthetic = 0L, jack = CC$jack_n, t = CC$t, seed = CC$jack_seed,
  cjack = CC$cjack_n, cseed = CC$cjack_seed, exp_n = CC$exp_n, workers = 8L, threads = 8L, lane = 0L, nlanes = 1L, collect = 0L))
SYN <- opt$synthetic == 1L; ROOT <- opt$root; G_HEAD <- cc_git_rc("rev-parse", "--short", "HEAD")$out[1]
LF <- file.path(ROOT, "net4/labelfree"); TR <- file.path(ROOT, "net4/truth"); EJ <- file.path(ROOT, "net4/expjack"); GA <- file.path(ROOT, "gate_a")
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))
common_guards <- function() { if (SYN) return(invisible(TRUE)); cc_design_guard(); cc_code_guard(); invisible(TRUE) }

## ---- the guard of the design (m9): no a3r score on net4 anywhere under results/ -------------------------------------------------
cc_guard_prior_a3r <- function(file, say) {
  fs <- list.files("results", recursive = TRUE, all.files = TRUE, full.names = TRUE, include.dirs = FALSE)
  fs <- fs[!startsWith(fs, paste0(CC$out, "/")) & !startsWith(fs, paste0(ROOT, "/"))]
  nm <- grepl("net4|dream5", fs, ignore.case = TRUE); in_dir <- grepl("net4|dream5", dirname(fs), ignore.case = TRUE)
  path_hit <- fs[nm & grepl("a3r", fs, ignore.case = TRUE)]
  textlike <- grepl("\\.(csv|tsv|txt|log|md|json|sha256|yaml|yml|R|py|sh|rc)$", fs, ignore.case = TRUE) | !grepl("\\.[A-Za-z0-9]+$", basename(fs))
  cand <- fs[in_dir & textlike & file.size(fs) > 0]
  header_hit <- cand[vapply(cand, function(f) { l <- tryCatch(readLines(f, n = 1L, warn = FALSE), error = function(e) ""); length(l) > 0 && grepl("a3r", l[1], ignore.case = TRUE) }, NA)]
  bin <- fs[in_dir & !textlike]
  pass <- !length(path_hit) && !length(header_hit)
  d <- data.frame(check = c("files searched under results/ (outside results/complex_confirm)", "files in a net4/dream5 path or directory", "path contains net4|dream5 together with a3r",
                            "header (first line) of a text file in a net4/dream5 directory contains a3r", "binary files in net4/dream5 directories (name checked only)", "GUARD PASSED"),
                  n = c(length(fs), sum(nm | in_dir), length(path_hit), length(header_hit), length(bin), as.integer(pass)),
                  detail = c("results/ recursive, hidden files included", "", paste(path_hit, collapse = ";"), paste(header_hit, collapse = ";"), "", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  cc_wcsv(d, file); print(d, row.names = FALSE)
  if (!pass) stop("GUARD FAILED: an a3r score on net4 exists under results/: ", paste(c(path_hit, header_hit), collapse = "; "))
  say("guard passed: no a3r score on net4 under results/ (", length(fs), " files searched)"); invisible(d)
}

get_truth_files <- function(nd) { tf <- cc_read_truth_files(); if (SYN) cc_synth_files(nd, tf) |> (\(x) c(x, list(counts = tf$counts, src = tf$src)))() else tf }

## ======================================================================================================================================
if (opt$mode == "gate_a") {
  dir.create(GA, recursive = TRUE, showWarnings = FALSE); log <- file.path(GA, "run.log"); say <- function(...) cc_say(log, ...)
  common_guards(); say("gate a, collection; git ", G_HEAD)
  li <- cc_lib_identity(); rows <- list(data.frame(part = "libraries", check = paste0("byte-identical to commit ", CC$libs_commit, ": ", li$file, " (sha256 ", substr(li$sha256, 1, 12), ")"),
                                                   max_abs_diff = as.numeric(!li$identical), pass = li$identical, stringsAsFactors = FALSE))
  for (p in c("k562", "m3d")) { f <- file.path(GA, p, paste0("equivalence_gate_", p, ".csv"))
    if (file.exists(f)) { g <- utils::read.csv(f, stringsAsFactors = FALSE); rows[[length(rows) + 1L]] <- g[, c("part", "check", "max_abs_diff", "pass")] }
    else rows[[length(rows) + 1L]] <- data.frame(part = toupper(p), check = paste("gate table missing:", f), max_abs_diff = Inf, pass = FALSE) }
  g <- do.call(rbind, rows); g$pass <- as.logical(g$pass); g$git_head <- G_HEAD; cc_wcsv(g, file.path(GA, "equivalence_gate.csv")); print(g[, c("part", "check", "max_abs_diff", "pass")], row.names = FALSE)
  if (!all(g$pass)) { say("GATE A FAILED"); quit(save = "no", status = 1L) }
  say("GATE A PASSED (", nrow(g), " checks)"); quit(save = "no")
}

## ======================================================================================================================================
if (opt$mode == "refit_timing") {
  dir.create(LF, recursive = TRUE, showWarnings = FALSE); log <- file.path(LF, "run.log"); say <- function(...) cc_say(log, ...)
  common_guards(); cc_check_fixed(opt)
  nd <- cc_load_net(opt$quick); base_rss <- cc_vmrss_gb(); say(sprintf("refit timing (group half rep 1, label-free): data loaded, RSS %.2f GB", base_rss))
  t0 <- Sys.time(); nb <- d5_subset_groups(nd, d5_group_half(nd$n_groups, 1L)); X <- nb$X
  stopifnot(!any(apply(X, 1, function(x) diff(range(x)) == 0)))
  R <- cc_fit_three(X, opt$workers); secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  h <- vapply(R, cc_pack_hash, ""); peak <- cc_vmhwm_gb()
  man <- utils::read.delim(file.path(D5$out, "net4/refit/REFIT_MANIFEST.tsv"), stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
  ok_clr <- !SYN && opt$quick == 0L && identical(h[["clr"]], man$sha256[man$rep == 1 & man$method == "clr"])
  ok_att <- !SYN && opt$quick == 0L && identical(h[["attention"]], man$sha256[man$rep == 1 & man$method == "attention"])
  d <- data.frame(unit = "group", rep = 1L, n_chips = ncol(X), workers = opt$workers, seconds_total = round(secs, 1), peak_rss_gb = round(peak, 2), rss_after_load_gb = round(base_rss, 2),
                  clr_hash_matches_dream5_manifest = ok_clr, attention_hash_matches_dream5_manifest = ok_att, clr_sha256 = h[["clr"]], attention_sha256 = h[["attention"]], a3r_sha256 = h[["a3r"]],
                  declared_mem_6030_if_one_lane_gb = ceiling(1.4 * peak), git_head = G_HEAD)
  cc_wcsv(d, file.path(LF, "refit_timing.csv")); print(d, row.names = FALSE)
  if (!SYN && opt$quick == 0L && !(ok_clr && ok_att)) stop("refit timing: CLR/attention hashes differ from the DREAM5 refit manifest")
  say(sprintf("refit timing done: %.0f s, peak RSS %.2f GB", secs, peak)); quit(save = "no")
}

## ======================================================================================================================================
if (opt$mode == "labelfree") {
  dir.create(file.path(LF, "cache"), recursive = TRUE, showWarnings = FALSE); log <- file.path(LF, "run.log"); say <- function(...) cc_say(log, ...)
  common_guards(); cc_check_fixed(opt); say("label-free stage; git ", G_HEAD, "; options ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  if (!SYN) cc_guard_prior_a3r(file.path(LF, "guard_no_prior_a3r.csv"), say)
  nd <- cc_load_net(opt$quick); G <- nrow(nd$X); say(sprintf("net4: %d genes x %d chips, %d replicate groups, %d series", G, ncol(nd$X), nd$n_groups, length(unique(nd$series))))
  t0 <- Sys.time(); B <- cc_fit_three(nd$X, opt$workers, base = TRUE); say(sprintf("fit, operator, P^9, three readouts: %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  C <- B$clr; ATT <- B$attention; S1 <- B$a3r; Pt <- B$Pt; P <- B$P; A <- B$A; mi_rm <- B$mi_rm
  eq <- list(); add <- function(check, v, tol, gating = TRUE) eq[[length(eq) + 1L]] <<- data.frame(check = check, value = v, tolerance = tol, gating = gating, pass = is.finite(v) && v <= tol, stringsAsFactors = FALSE)
  chk_cache <- function(S, name) { f <- file.path(D5$out, "net4", paste0(name, "__replicate.scores"))
    sidecar <- utils::read.delim(paste0(f, ".sha256"), stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
    stored <- d5_read_scores(f); v <- d5_pack(S); md <- max(abs(v - stored)); tmp <- file.path(LF, "cache", ".tmp_hash.scores"); d5_write_scores(v, tmp); same <- identical(d5_sha256(tmp), sidecar$sha256); unlink(tmp)
    list(maxdiff = md, same_hash = same) }
  if (opt$quick == 0L) {
    r <- chk_cache(C, "clr"); add("Ga4: new CLR scores vs results/dream5/net4/clr__replicate.scores, max |diff|", r$maxdiff, 1e-12); add("Ga4: new CLR scores, SHA-256 of the packed file identical (1 = yes; informative)", as.numeric(!r$same_hash), 1, FALSE)
    r <- chk_cache(ATT, "attention"); add("Ga4: new attention scores vs results/dream5/net4/attention__replicate.scores, max |diff|", r$maxdiff, 1e-12); add("Ga4: new attention scores, SHA-256 identical (informative)", as.numeric(!r$same_hash), 1, FALSE)
    A4b <- d5_clean_sym((Pt + t(Pt)) / 2)          # from the very P^9 object that was passed to sk_a3r
    r <- chk_cache(A4b, "attention"); add("Ga4b: (P9 + P9')/2 from the P9 fed to sk_a3r vs attention__replicate.scores, max |diff|", r$maxdiff, 1e-12); rm(A4b)
  } else say("quick run: Ga part 4 skipped (gene subset)")
  A_rz <- rz_operator(C); add("info: package operator vs rz_operator(CLR), max |diff| (stack path vs package path)", max(abs(A_rz - A)), 1e-12, FALSE); rm(A_rz)
  cc_wcsv(do.call(rbind, eq), file.path(LF, "equivalence_net4.csv")); print(do.call(rbind, eq), row.names = FALSE)
  if (!all(vapply(eq, function(r) r$pass || !r$gating, NA))) stop("GATE Ga (net4, parts 4 / 4b) FAILED: no truth is read, nothing is scored")
  ## operator checks and Gb
  UT <- which(upper.tri(matrix(FALSE, G, G))); q_att <- rz_stationary(P); reach <- sk_reach(Pt); top1_C <- sk_top1(C, UT)
  dg <- sk_diag("dream5_net4", 1L, C, A, P, reach, S1, top1_C, top1_C, mi_rm, q_att, UT)
  L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA; mu <- rowMeans(L, na.rm = TRUE); nn <- rowSums(!is.na(L)); sdv <- sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (nn - 1)); rm(L)
  q <- stats::quantile(sdv[is.finite(sdv)], c(0, .05, .5, .95, 1), names = FALSE)
  oc <- cbind(dg, data.frame(G = G, sigma_min = q[1], sigma_p05 = q[2], sigma_med = q[3], sigma_p95 = q[4], sigma_max = q[5], git_head = G_HEAD)); cc_wcsv(oc, file.path(LF, "operator_checks.csv")); print(t(oc))
  dec <- grepl("^decoy", nd$gene_names); real <- which(!dec); Pr <- Pt[real, , drop = FALSE]; rsum <- rowSums(Pr); dshare <- rowSums(Pr[, dec, drop = FALSE]) / rsum
  Rm <- Pr > 0; Rm[cbind(seq_along(real), real)] <- FALSE
  dd <- data.frame(n_decoys = sum(dec), n_real = length(real), decoy_fraction_of_genes = mean(dec), mean_decoy_share_of_P9_row_mass_real_rows = mean(dshare), median_decoy_share_of_P9_row_mass = stats::median(dshare),
                   decoy_share_of_reachable_entries_real_rows = sum(Rm[, dec, drop = FALSE]) / sum(Rm), git_head = G_HEAD)
  cc_wcsv(dd, file.path(LF, "decoy_diagnostic.csv")); print(t(dd)); rm(Pr, Rm)
  tf <- get_truth_files(nd); T <- cc_build_truths(nd, tf, synthetic = SYN, ym = FALSE); main <- T$main; cn <- tf$counts[tf$counts$variant == "main", ]
  hashes <- if (SYN) rep(NA_character_, 2) else cc_truth_hash_guard()
  gb <- list(); addg <- function(gate, v, thr, ok) gb[[length(gb) + 1L]] <<- data.frame(gate = gate, value = v, threshold = thr, pass = ok, stringsAsFactors = FALSE)
  addg("Gb1 rows with an empty candidate list in the operator (CLR row has no positive entry)", dg$rows_no_candidate, "== 0", dg$rows_no_candidate == 0)
  addg("Gb2 unreachable share of off-diagonal P^9 pairs", dg$unreach_share, "<= 0.05", dg$unreach_share <= 0.05)
  addg("Gb3 rows with fewer than 2 reachable entries", dg$rows_lt2_reachable, "== 0", dg$rows_lt2_reachable == 0)
  addg("Gb4 all-zero rows of the a3r matrix S", dg$S_rows_all_zero, "== 0", dg$S_rows_all_zero == 0)
  if (!SYN) {
    addg("Gb5 members file SHA-256 matches the frozen hash (667cd567...)", 1, "match", TRUE); addg("Gb5 excluded-pairs file SHA-256 matches the frozen hash (b234b94d...)", 1, "match", TRUE)
    cov <- tf$src$coverage[tf$src$dataset == "dream5_net4"]; addg("Gb5 coverage of complex ORFs (complexes with >= 3 ORFs) in the net4 universe", cov, ">= 0.90", cov >= 0.90)
    addg("Gb5 merged complexes with positives (main)", sum(lengths(main$sets) >= 2L), ">= 100", sum(lengths(main$sets) >= 2L) >= 100)
    addg("Gb5 n_eff (frozen counts file, main)", cn$n_eff, ">= 10", cn$n_eff >= 10); addg("Gb5 positives (main)", sum(main$label), ">= 5000", sum(main$label) >= 5000)
    addg("Gb5 paralog pairs inside the universe (main)", main$n_removed, "== 49", main$n_removed == 49L)
    exp_cnt <- function(v) tf$counts[tf$counts$variant == v, ]
    for (v in c("main", "main_xsub_excl", "main_ribo_joined", "noribo", "nomerge", "size2")) { p <- T[[v]]; e <- exp_cnt(v)
      addg(sprintf("Gb5 structure of %s equals the frozen counts: positives (%d) and universe pairs (%d)", v, e$n_pos, e$n_pairs), sum(p$label), paste0("== ", e$n_pos, " and ", e$n_pairs), sum(p$label) == e$n_pos && length(p$label) == e$n_pairs) }
    for (pr in list(c("small_3to6", "small_3to6"), c("large_ge20", "large_ge20"))) { p <- T[[pr[1]]]; e <- exp_cnt(pr[2]); addg(sprintf("Gb5 stratum %s positives equal the frozen counts (%d)", pr[1], e$n_pos), sum(p$label), paste0("== ", e$n_pos), sum(p$label) == e$n_pos) }
    addg("Gb5 main universe pairs equal the frozen counts (1,282,352) and positives (13,438)", length(main$label), "== 1282352", length(main$label) == 1282352L && sum(main$label) == 13438L)
  }
  gbd <- do.call(rbind, gb); gbd$git_head <- G_HEAD; cc_wcsv(gbd, file.path(LF, "feasibility_gate.csv")); print(gbd, row.names = FALSE)
  if (!all(gbd$pass)) stop("GATE Gb (label-free feasibility) FAILED: outcome 'not evaluable'; nothing is scored")
  say("Ga parts 4/4b and Gb passed")
  ## caches (the three operator-derived score matrices, packed upper triangle, and the MI row mean for the hub stratifier) and their manifest
  files <- cc_cache_files(ROOT); vals <- list(d5_pack(C), d5_pack(ATT), d5_pack(S1), mi_rm); rows <- list()
  for (k in seq_along(files)) { d5_write_scores(vals[[k]], files[k])
    rows[[k]] <- data.frame(file = files[k], bytes = file.size(files[k]), sha256 = d5_sha256(files[k]), n_values = length(vals[[k]]), dims = G, finished_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), git_head = G_HEAD) }
  man <- do.call(rbind, rows); cc_wtsv(man, file.path(LF, "score_manifest.tsv")); print(man[, c("file", "bytes", "sha256")], row.names = FALSE)
  say(sprintf("label-free stage done; peak RSS %.2f GB", cc_vmhwm_gb())); quit(save = "no")
}

## ======================================================================================================================================
## shared by truth / expjack / summarize: production guards
production_guards <- function(stage) {
  if (SYN) return(invisible(TRUE))
  cc_check_fixed(opt); cc_design_guard(); cc_code_guard(); cc_truth_hash_guard(); if (stage != "summarize") cc_cache_guard(ROOT)
  invisible(TRUE)
}
est_files <- c(file.path(ROOT, "estimate.csv"), file.path(ROOT, "estimate_final.csv"))
truth_outputs <- c("truth_counts.csv", "metrics.csv", "contrasts.csv", "jack_gene_replicates.csv", "jack_complex_replicates.csv", "per_complex_ap.csv", "per_complex_summary.csv", "strata.csv", "rp_module_coherence.csv")

if (opt$mode == "truth") {
  dir.create(TR, recursive = TRUE, showWarnings = FALSE); log <- file.path(TR, "run.log"); say <- function(...) cc_say(log, ...)
  production_guards("truth")
  if (any(file.exists(file.path(TR, truth_outputs))) || any(file.exists(est_files))) stop("guard: a truth-stage output, estimate.csv or estimate_final.csv already exists: one run only")
  say("TRUTH STAGE start; git ", G_HEAD, "; options ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))
  if (!SYN) cc_guard_prior_a3r(file.path(TR, "guard_no_prior_a3r_6020.csv"), say)
  nd <- cc_load_net(opt$quick); G <- nrow(nd$X); tf <- get_truth_files(nd); rm_x <- nd$X; nd$X <- matrix(0, G, 1)   # the truth stage never reads expression values
  rm(rm_x); invisible(gc(verbose = FALSE))
  ## secondary item 8: the replicate-averaged input (321 x 5,950), the same three readouts (label-free fit), scored on the main truth below
  nda <- cc_load_net(opt$quick, "avg"); Ra <- cc_fit_three(nda$X, opt$workers); Savg <- stats::setNames(Ra[c("clr", "attention", "a3r")], paste0("avg_", c("clr", "attention", "a3r"))); rm(nda, Ra)
  say(sprintf("replicate-averaged input fitted (321 groups); RSS %.2f GB", cc_vmrss_gb())); invisible(gc(verbose = FALSE))
  T <- cc_build_truths(nd, tf, synthetic = SYN, ym = TRUE); tn <- names(T)
  halves <- cc_draw_halves(T, opt$jack, opt$seed)      # RNG: set.seed(20261114), then all truths in table order
  say("truths: ", paste(sprintf("%s (%d pairs, %d pos)", tn, vapply(T, function(p) length(p$label), 0L), vapply(T, function(p) sum(p$label), 0L)), collapse = "; "))
  cnt <- do.call(rbind, lapply(tn, function(v) { p <- T[[v]]; e <- if (v %in% tf$counts$variant) tf$counts[tf$counts$variant == v, ] else NULL
    data.frame(truth = v, role = p$role, n_complexes = length(p$sets), n_genes_in = length(p$genes_in), n_pairs = length(p$label), n_pos = sum(p$label), base_rate = mean(p$label), n_pairs_removed = p$n_removed,
               frozen_n_pairs = if (is.null(e) || v %in% c("small_3to6", "large_ge20")) NA else e$n_pairs, frozen_n_pos = if (is.null(e)) NA else e$n_pos, structure_md5 = cc_truth_struct_hash(p), note = p$note, stringsAsFactors = FALSE) }))
  man <- if (SYN) NULL else cc_cache_guard(ROOT)
  C <- d5_unpack(d5_read_scores(cc_cache_files(ROOT)[1]), G); ATT <- d5_unpack(d5_read_scores(cc_cache_files(ROOT)[2]), G); A3 <- d5_unpack(d5_read_scores(cc_cache_files(ROOT)[3]), G)
  mi_rm <- d5_read_scores(cc_cache_files(ROOT)[4]); H <- outer(mi_rm, mi_rm, "+")
  Sall <- list(clr = C, attention = ATT, a3r = A3); rm(C, ATT, A3)
  for (m in c("clr", "attention", "a3r")) { Sall[[paste0("hub_", m)]] <- rz_hub_stratify(Sall[[m]], H); invisible(gc(verbose = FALSE)) }
  rm(H); say(sprintf("six score matrices ready (raw and hub-stratified); RSS %.2f GB", cc_vmrss_gb()))
  ## pooled AUPR and the paired delete-half gene jackknife, every truth, both scopes
  met <- list(); rep_l <- list(); con <- list()
  for (v in tn) { t1 <- Sys.time(); invisible(gc(verbose = FALSE)); p <- T[[v]]
    r <- cc_score_truth(Sall, p, halves[[v]], G, nc)
    for (m in names(r$aupr)) { sc <- if (startsWith(m, "hub_")) "hub_mi" else "raw"; rd <- sub("^hub_", "", m)
      met[[length(met) + 1L]] <- data.frame(truth = v, role = p$role, scope = sc, readout = rd, aupr = r$aupr[[m]], n_pairs = length(p$label), n_pos = sum(p$label), base_rate = mean(p$label), enrichment = r$aupr[[m]] / mean(p$label), stringsAsFactors = FALSE)
      rep_l[[length(rep_l) + 1L]] <- data.frame(truth = v, scope = sc, readout = rd, rep = seq_len(nrow(r$jack)), aupr = r$jack[, m], stringsAsFactors = FALSE) }
    for (sc in c("raw", "hub_mi")) { pf <- if (sc == "raw") "" else "hub_"
      for (cn_ in names(CONTRASTS)) { a <- CONTRASTS[[cn_]]; d <- r$aupr[[paste0(pf, a[1])]] - r$aupr[[paste0(pf, a[2])]]; reps <- r$jack[, paste0(pf, a[1])] - r$jack[, paste0(pf, a[2])]
        con[[length(con) + 1L]] <- cbind(data.frame(truth = v, role = p$role, scope = sc, row = cn_, contrast = cc_contrast_names[[cn_]], stringsAsFactors = FALSE), cc_interval(d, reps)) } }
    say(sprintf("truth %-18s scored and jackknifed (%.0f s)", v, as.numeric(difftime(Sys.time(), t1, units = "secs")))) }
  ra <- cc_score_truth(Savg, T$main, halves$main, G, nc); rm(Savg)       # same halves as the main truth (paired draws), raw scope only
  for (m in names(ra$aupr)) { rd <- sub("^avg_", "", m)
    met[[length(met) + 1L]] <- data.frame(truth = "main", role = "descriptive_avg_input", scope = "avg_input", readout = rd, aupr = ra$aupr[[m]], n_pairs = length(T$main$label), n_pos = sum(T$main$label), base_rate = mean(T$main$label), enrichment = ra$aupr[[m]] / mean(T$main$label), stringsAsFactors = FALSE)
    rep_l[[length(rep_l) + 1L]] <- data.frame(truth = "main", scope = "avg_input", readout = rd, rep = seq_len(nrow(ra$jack)), aupr = ra$jack[, m], stringsAsFactors = FALSE) }
  for (cn_ in names(CONTRASTS)) { a <- CONTRASTS[[cn_]]; con[[length(con) + 1L]] <- cbind(data.frame(truth = "main", role = "descriptive_avg_input", scope = "avg_input", row = cn_, contrast = cc_contrast_names[[cn_]], stringsAsFactors = FALSE),
    cc_interval(ra$aupr[[paste0("avg_", a[1])]] - ra$aupr[[paste0("avg_", a[2])]], ra$jack[, paste0("avg_", a[1])] - ra$jack[, paste0("avg_", a[2])])) }
  say("replicate-averaged input scored on the main truth (descriptive)")
  met <- do.call(rbind, met); repdf <- do.call(rbind, rep_l); con <- do.call(rbind, con)
  ## complex-level jackknife (labelled sensitivity; raw scope, main truth)
  k_all <- cc_excl_keys(tf$excl, "all", nd); cmain <- T$main$sets; nC <- length(cmain)
  set.seed(opt$cseed); chalves <- lapply(seq_len(opt$cjack), function(r) sort(sample.int(nC, floor(nC / 2)))); set.seed(NULL)
  cone <- function(keep) { p <- cc_pairs(cmain[keep], nd, k_all); ix <- cbind(p$i, p$j)
    c(vapply(c("clr", "attention", "a3r"), function(m) ap_fast(Sall[[m]][ix], p$label), 0), n_pairs = length(p$label), n_pos = sum(p$label)) }
  invisible(gc(verbose = FALSE)); t1 <- Sys.time()
  cj <- if (nc > 1L) parallel::mclapply(chalves, cone, mc.cores = nc) else lapply(chalves, cone)
  if (any(vapply(cj, function(x) inherits(x, "try-error") || length(x) != 5L, NA))) stop("a complex-jackknife worker failed")
  cj <- do.call(rbind, cj); say(sprintf("complex-level jackknife: %d replicates (%.0f s)", nrow(cj), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  cjl <- do.call(rbind, lapply(c("clr", "attention", "a3r"), function(m) data.frame(rep = seq_len(nrow(cj)), readout = m, aupr = cj[, m], n_pairs = cj[, "n_pairs"], n_pos = cj[, "n_pos"])))
  pm <- met[met$truth == "main" & met$scope == "raw", ]; pa <- stats::setNames(pm$aupr, pm$readout)
  for (cn_ in names(CONTRASTS)) { a <- CONTRASTS[[cn_]]; con <- rbind(con, cbind(data.frame(truth = "main", role = "complex_level_jackknife_sensitivity", scope = "raw", row = cn_, contrast = cc_contrast_names[[cn_]], stringsAsFactors = FALSE),
    cc_interval(pa[[a[1]]] - pa[[a[2]]], cj[, a[1]] - cj[, a[2]]))) }
  ## per-complex AP (macro; complexes with >= 10 positive pairs; the regulon_ap_table definitions of the 26 Oct test) and r-protein module coherence
  opm <- cc_paralog_map(tf$excl, nd); regm <- list(members = cmain, class = stats::setNames(rep("complex", nC), names(cmain))); co <- comember_matrix(regm, G)
  pcl <- list(); for (m in c("clr", "attention", "a3r")) { x <- regulon_ap_table(Sall[[m]], regm, co, opm, nd$gene_names, min_pos = 10L); pcl[[m]] <- cbind(readout = m, x) }
  pc <- do.call(rbind, pcl); rm(co)
  pcs <- do.call(rbind, lapply(list(c("a3r", "attention"), c("attention", "clr"), c("a3r", "clr")), function(ab) { a <- pcl[[ab[1]]]; b <- pcl[[ab[2]]]; ok <- !is.na(a$ap) & !is.na(b$ap); d <- a$ap[ok] - b$ap[ok]
    data.frame(comparison = paste(ab[1], "minus", ab[2]), n_complexes_ge10pos = sum(ok), mean_d_ap = mean(d), median_d_ap = stats::median(d), n_improved = sum(d > 0), n_worse = sum(d < 0), stringsAsFactors = FALSE) }))
  rp <- cc_sets(tf$mem, "rp_module", nd); co_rp <- comember_matrix(list(members = rp), G)
  rpc <- do.call(rbind, lapply(c("clr", "attention", "a3r"), function(m) data.frame(readout = m, n_genes = length(rp[[1]]), coherence_auroc = regulon_coherence(Sall[[m]], list(members = rp), co_rp, opm, nd$gene_names), stringsAsFactors = FALSE)))
  ## write everything (staging directory, then rename)
  st <- file.path(TR, paste0("_staging", Sys.getpid())); dir.create(st)
  outs <- list(truth_counts.csv = cnt, metrics.csv = met, contrasts.csv = con, jack_gene_replicates.csv = repdf, jack_complex_replicates.csv = cjl, per_complex_ap.csv = pc,
               per_complex_summary.csv = pcs, strata.csv = con[con$truth %in% c("small_3to6", "mid_7to19", "large_ge20") & con$scope == "raw", ], rp_module_coherence.csv = rpc)
  for (f in names(outs)) utils::write.csv(outs[[f]], file.path(st, f), row.names = FALSE)
  for (f in names(outs)) file.rename(file.path(st, f), file.path(TR, f)); unlink(st, recursive = TRUE)
  print(con[con$truth == "main" & con$scope == "raw", c("contrast", "role", "d", "se_gene", "gene_lo", "gene_hi", "skew")], row.names = FALSE)
  say(sprintf("TRUTH STAGE done; peak RSS %.2f GB", cc_vmhwm_gb())); quit(save = "no")
}

## ======================================================================================================================================
if (opt$mode == "expjack") {
  dir.create(file.path(EJ, "reps"), recursive = TRUE, showWarnings = FALSE); log <- file.path(EJ, if (opt$collect == 1L) "run.log" else sprintf("run_lane%d.log", opt$lane)); say <- function(...) cc_say(log, ...)
  tasks <- rbind(data.frame(unit = "series", b = seq_len(opt$exp_n), stringsAsFactors = FALSE), data.frame(unit = "group", b = seq_len(opt$exp_n), stringsAsFactors = FALSE))
  repf <- function(u, b) file.path(EJ, "reps", sprintf("%s%03d.csv", u, b)); hashf <- function(u, b) file.path(EJ, "reps", sprintf("%s%03d_hash.csv", u, b))
  if (opt$collect == 1L) {
    production_guards("summarize")
    miss <- tasks[!file.exists(mapply(repf, tasks$unit, tasks$b)), ]; say("collect: ", nrow(tasks) - nrow(miss), " of ", nrow(tasks), " replicates present")
    for (u in c("series", "group")) { x <- tasks[tasks$unit == u & file.exists(mapply(repf, tasks$unit, tasks$b)), ]
      if (nrow(x)) cc_wcsv(do.call(rbind, lapply(seq_len(nrow(x)), function(k) utils::read.csv(repf(x$unit[k], x$b[k]), stringsAsFactors = FALSE))), file.path(EJ, paste0("replicates_", u, ".csv")))
      h <- tasks[tasks$unit == u & file.exists(mapply(hashf, tasks$unit, tasks$b)), ]
      if (nrow(h)) { hd <- do.call(rbind, lapply(seq_len(nrow(h)), function(k) utils::read.csv(hashf(h$unit[k], h$b[k]), stringsAsFactors = FALSE))); hd$git_head <- G_HEAD; cc_wcsv(hd, file.path(EJ, paste0("refit_hash_check_", u, ".csv"))) } }
    allh <- do.call(rbind, lapply(c("series", "group"), function(u) { f <- file.path(EJ, paste0("refit_hash_check_", u, ".csv")); if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) }))
    if (!is.null(allh)) cc_wcsv(allh, file.path(EJ, "refit_hash_check.csv"))
    if (nrow(miss)) { say("MISSING replicates: ", paste(sprintf("%s%d", miss$unit, miss$b), collapse = " ")); quit(save = "no", status = 2L) }
    if (!SYN && !all(allh$clr_hash_ok %in% TRUE & allh$attention_hash_ok %in% TRUE)) { say("HASH MISMATCH in at least one replicate"); quit(save = "no", status = 3L) }
    say("collect done: all 200 replicates present, all hashes match"); quit(save = "no")
  }
  production_guards("expjack")
  nd <- cc_load_net(opt$quick); G <- nrow(nd$X); tf <- get_truth_files(nd); T <- cc_build_truths(nd, tf, synthetic = SYN, ym = TRUE)
  if (!SYN) { tc <- utils::read.csv(file.path(TR, "truth_counts.csv"), stringsAsFactors = FALSE)
    if (!all(vapply(names(T), function(v) identical(tc$structure_md5[tc$truth == v], cc_truth_struct_hash(T[[v]])), NA))) stop("the truth built here differs from the truth-stage structure (truth_counts.csv)")
    say("truth structure identical to the truth stage (md5 of i, j, label for every truth)") }
  gman <- file.path(D5$out, "net4/refit/REFIT_MANIFEST.tsv"); sman <- file.path(D5$out, "net4/refit_series/REFIT_SERIES_MANIFEST.tsv")
  mine <- which((seq_len(nrow(tasks)) - 1L) %% opt$nlanes == opt$lane)
  say(sprintf("lane %d of %d: %d replicates; git %s", opt$lane, opt$nlanes, length(mine), G_HEAD))
  for (k in mine) { u <- tasks$unit[k]; b <- tasks$b[k]; if (file.exists(repf(u, b)) && file.exists(hashf(u, b))) next
    t0 <- Sys.time(); nb <- if (u == "group") d5_subset_groups(nd, d5_group_half(nd$n_groups, b)) else d5_subset_series(nd, b); X <- nb$X
    if (any(apply(X, 1, function(x) diff(range(x)) == 0))) stop("replicate ", u, b, ": zero-variance genes in the kept half")
    R <- cc_fit_three(X, opt$workers); rm(nb, X)
    hs <- list(clr = cc_pack_hash(R$clr), attention = cc_pack_hash(R$attention), a3r = cc_pack_hash(R$a3r))
    ok <- c(clr = NA, attention = NA)
    if (!SYN && opt$quick == 0L) { mf <- if (u == "group") gman else sman; man <- utils::read.delim(mf, stringsAsFactors = FALSE, colClasses = c(sha256 = "character"))
      for (m in names(ok)) { w <- man$sha256[man$net == 4 & man$rep == b & man$method == m]; ok[[m]] <- length(w) == 1L && identical(hs[[m]], w) } }
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cc_wcsv(data.frame(unit = u, rep = b, clr_hash_ok = ok[["clr"]], attention_hash_ok = ok[["attention"]], clr_sha256 = hs$clr, attention_sha256 = hs$attention, a3r_sha256 = hs$a3r, seconds = round(secs, 1), peak_rss_gb = round(cc_vmhwm_gb(), 2)), hashf(u, b))
    if (!SYN && opt$quick == 0L && !all(ok)) stop("refit replicate ", u, b, " differs from the committed DREAM5 hash (", paste(names(ok)[!ok], collapse = ","), "); not scoring")
    rows <- list(); for (v in names(T)) { r <- cc_score_truth(R[c("clr", "attention", "a3r")], T[[v]], NULL, G, 1L)
      rows[[v]] <- data.frame(unit = u, rep = b, truth = v, readout = names(r$aupr), aupr = unname(r$aupr), stringsAsFactors = FALSE) }
    cc_wcsv(do.call(rbind, rows), repf(u, b)); rm(R); invisible(gc(verbose = FALSE))
    say(sprintf("%s %3d done (%.0f s, peak RSS %.2f GB)", u, b, secs, cc_vmhwm_gb())) }
  say("lane done"); quit(save = "no")
}

## ======================================================================================================================================
if (opt$mode == "summarize") {
  dir.create(ROOT, recursive = TRUE, showWarnings = FALSE); log <- file.path(ROOT, "summarize.log"); say <- function(...) cc_say(log, ...)
  production_guards("summarize")
  if (any(file.exists(est_files))) stop("guard: estimate.csv or estimate_final.csv already exists; estimate.csv is written once and never rewritten")
  stopifnot(all(file.exists(file.path(TR, c("metrics.csv", "jack_gene_replicates.csv", "contrasts.csv")))), file.exists(file.path(EJ, "replicates_group.csv")))
  met <- utils::read.csv(file.path(TR, "metrics.csv"), stringsAsFactors = FALSE); jg <- utils::read.csv(file.path(TR, "jack_gene_replicates.csv"), stringsAsFactors = FALSE)
  con <- utils::read.csv(file.path(TR, "contrasts.csv"), stringsAsFactors = FALSE)
  ejg <- utils::read.csv(file.path(EJ, "replicates_group.csv"), stringsAsFactors = FALSE)
  ejs <- if (file.exists(file.path(EJ, "replicates_series.csv"))) utils::read.csv(file.path(EJ, "replicates_series.csv"), stringsAsFactors = FALSE) else NULL
  series_ok <- !is.null(ejs) && length(unique(ejs$rep)) == opt$exp_n && all(table(ejs$rep) == table(ejs$rep)[1])
  gene_reps <- function(truth, scope, a, b) { x <- jg[jg$truth == truth & jg$scope == scope, ]; xa <- x[x$readout == a, ]; xb <- x[x$readout == b, ]; xa$aupr[order(xa$rep)] - xb$aupr[order(xb$rep)] }
  exp_diff <- function(df, truth, a, b) { if (is.null(df)) return(NULL); x <- df[df$truth == truth, ]; xa <- x[x$readout == a, ]; xb <- x[x$readout == b, ]; xa$aupr[order(xa$rep)] - xb$aupr[order(xb$rep)] }
  pt <- function(truth, scope, rd) met$aupr[met$truth == truth & met$scope == scope & met$readout == rd]
  rows <- list()
  for (cn_ in names(CONTRASTS)) { a <- CONTRASTS[[cn_]]; d <- pt("main", "raw", a[1]) - pt("main", "raw", a[2]); gi <- cc_interval(d, gene_reps("main", "raw", a[1], a[2]))
    sg <- exp_diff(ejg, "main", a[1], a[2]); ss <- exp_diff(ejs, "main", a[1], a[2]); se_g <- stats::sd(sg); se_s <- if (series_ok) stats::sd(ss) else NA_real_
    tot_s <- sqrt(gi$se_gene^2 + se_s^2); tot_g <- sqrt(gi$se_gene^2 + se_g^2)
    rows[[cn_]] <- data.frame(row = cn_, contrast = cc_contrast_names[[cn_]], truth = "main", scope = "raw", d = d, se_gene = gi$se_gene, se_series = se_s, se_group = se_g,
      total_series_lo = d - 1.96 * tot_s, total_series_hi = d + 1.96 * tot_s, total_group_lo = d - 1.96 * tot_g, total_group_hi = d + 1.96 * tot_g, gene_lo = gi$gene_lo, gene_hi = gi$gene_hi,
      gene_pct_lo = gi$pct_lo, gene_pct_hi = gi$pct_hi, gene_skewness = gi$skew, n_gene_reps = gi$n_rep, n_series_reps = if (is.null(ejs)) 0L else length(ss), n_group_reps = length(sg), stringsAsFactors = FALSE) }
  est <- do.call(rbind, rows); P <- est[est$row == "P", ]
  basis <- if (series_ok) "total (gene + series) 95% interval" else "gene + group interval (series component missing; group-level used; deviation logged)"
  L <- if (series_ok) P$total_series_lo else P$total_group_lo; U <- if (series_ok) P$total_series_hi else P$total_group_hi
  sent <- cc_sentences(L, U); sg_ <- cc_sentences(P$gene_lo, P$gene_hi); sgr <- cc_sentences(P$total_group_lo, P$total_group_hi)
  est$interval_basis_for_sentences <- ""; est$flag <- ""; est$interval_basis_for_sentences[est$row == "P"] <- basis
  if (!series_ok) est$flag[est$row == "P"] <- "series component missing; group-level used"
  for (s in names(sent)) est[[s]] <- ifelse(est$row == "P", sent[[s]], "")
  for (s in setdiff(names(sg_), "S_fixed")) est[[paste0("side_gene_only_", s)]] <- ifelse(est$row == "P", paste("[gene jackknife only, comparable to the E. coli test]", sg_[[s]]), "")
  for (s in setdiff(names(sgr), "S_fixed")) est[[paste0("side_group_level_", s)]] <- ifelse(est$row == "P", paste("[group-level experiment unit]", sgr[[s]]), "")
  est$note <- CC_NOTE_CLR; est$git_head <- G_HEAD
  ## estimate_final.csv: every secondary and sensitivity line (gene-only interval; experiment SDs where the refits scored the truth)
  fin <- list()
  for (k in seq_len(nrow(con))) { x <- con[k, ]; a <- CONTRASTS[[x$row]]
    sgd <- if (x$scope == "raw" && x$role != "complex_level_jackknife_sensitivity") stats::sd(exp_diff(ejg, x$truth, a[1], a[2])) else NA_real_
    ssd <- if (x$scope == "raw" && x$role != "complex_level_jackknife_sensitivity" && !is.null(ejs)) stats::sd(exp_diff(ejs, x$truth, a[1], a[2])) else NA_real_
    fin[[k]] <- data.frame(line_type = if (x$role == "complex_level_jackknife_sensitivity") "complex_level_jackknife (sensitivity; SE_C)" else "gene_jackknife_interval", truth = x$truth, role = x$role, scope = x$scope, row = x$row, contrast = x$contrast,
      d = x$d, se = x$se_gene, lo = x$gene_lo, hi = x$gene_hi, pct_lo = x$pct_lo, pct_hi = x$pct_hi, skewness = x$skew, n_rep = x$n_rep, descriptive_sd_group_refits = sgd, descriptive_sd_series_refits = ssd, stringsAsFactors = FALSE) }
  fin <- do.call(rbind, fin)
  pcs <- utils::read.csv(file.path(TR, "per_complex_summary.csv"), stringsAsFactors = FALSE); rpc <- utils::read.csv(file.path(TR, "rp_module_coherence.csv"), stringsAsFactors = FALSE)
  extra <- rbind(data.frame(line_type = "per_complex_macro_AP (>= 10 positives)", truth = "main", role = "descriptive", scope = "raw", row = "", contrast = pcs$comparison, d = pcs$mean_d_ap, se = NA, lo = NA, hi = NA, pct_lo = NA, pct_hi = NA, skewness = NA, n_rep = pcs$n_complexes_ge10pos,
                            descriptive_sd_group_refits = NA, descriptive_sd_series_refits = NA, stringsAsFactors = FALSE),
                 data.frame(line_type = "rp_module_coherence_AUROC (descriptive)", truth = "rp_module", role = "descriptive", scope = "raw", row = "", contrast = rpc$readout, d = rpc$coherence_auroc, se = NA, lo = NA, hi = NA, pct_lo = NA, pct_hi = NA, skewness = NA, n_rep = rpc$n_genes,
                            descriptive_sd_group_refits = NA, descriptive_sd_series_refits = NA, stringsAsFactors = FALSE))
  fin <- rbind(fin, extra); fin$git_head <- G_HEAD
  f1 <- est_files[1]; f2 <- est_files[2]
  for (f in c(f1, f2)) if (file.exists(f)) stop("guard: ", f, " exists")
  t1 <- paste0(f1, ".tmp", Sys.getpid()); t2 <- paste0(f2, ".tmp", Sys.getpid()); utils::write.csv(est, t1, row.names = FALSE); utils::write.csv(fin, t2, row.names = FALSE)
  if (!file.link(t1, f1)) stop("could not create estimate.csv without overwriting"); file.remove(t1); if (!file.link(t2, f2)) stop("could not create estimate_final.csv without overwriting"); file.remove(t2)
  ## tables.md
  fm <- function(x) sprintf("%+.4f", x); f4 <- function(x) sprintf("%.4f", x)
  hdr <- "| Row | Contrast | Point estimate | Total 95% interval (gene + series) | Total interval (gene + group) | Gene-only interval (E. coli-comparable) | Gene replicates: 2.5% / 97.5% percentiles (centred on d), skewness |\n|---|---|---|---|---|---|---|"
  ln <- vapply(seq_len(nrow(est)), function(k) { e <- est[k, ]; sprintf("| %s | %s | %s | [%s, %s] | [%s, %s] | [%s, %s] | [%s, %s]; %.3f |", e$row, e$contrast, fm(e$d), fm(e$total_series_lo), fm(e$total_series_hi), fm(e$total_group_lo), fm(e$total_group_hi),
                                                                  fm(e$gene_lo), fm(e$gene_hi), fm(e$gene_pct_lo), fm(e$gene_pct_hi), e$gene_skewness) }, "")
  md <- c("# Complex estimation on DREAM5 net4: main table (written by the summariser; no interpretation)", "", hdr, ln, "", paste0("> ", CC_NOTE_CLR), "",
          sprintf("SE (main truth, raw): gene %s; series %s; group %s (P row)", f4(P$se_gene), f4(P$se_series), f4(P$se_group)), "",
          paste0("## Fixed sentences (row P; basis: ", basis, ")"), "", paste0("- ", sent), "", "Side-by-side lines (never the reported sentences):", "",
          paste0("- ", vapply(setdiff(names(sg_), "S_fixed"), function(s) paste("[gene jackknife only, comparable to the E. coli test]", sg_[[s]]), "")),
          paste0("- ", vapply(setdiff(names(sgr), "S_fixed"), function(s) paste("[group-level experiment unit]", sgr[[s]]), "")), "",
          "See estimate_final.csv for every secondary and sensitivity line (all gene-jackknife intervals; no p-value anywhere).")
  writeLines(md, file.path(ROOT, "tables.md"))
  print(est[, c("row", "contrast", "d", "se_gene", "se_series", "se_group", "total_series_lo", "total_series_hi", "gene_lo", "gene_hi")], row.names = FALSE); print(sent)
  say("summarize done: estimate.csv and estimate_final.csv written"); quit(save = "no")
}
stop("unknown --mode: ", opt$mode)
