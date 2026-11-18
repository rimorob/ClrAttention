#!/usr/bin/env Rscript
# Scoring and reporting driver of docs/plans/tcaclr_confirm_plan.md (draft 3a), one arm per call. Run from the repo root.
#
#   Rscript analysis/tcaclr_confirm/score_confirm.R --stage <gates|score|expjack|collect|report> --arm A|B --scores DIR --truth_dir DIR \
#       --members FILE --excluded FILE --sha FILE --out DIR [--expr FILE] [--perturbed FILE] [--counts FILE] ...
#
# Stages (each writes CSVs under --out; the hash guard runs first in every stage and the job stops on a mismatch):
#   gates    SHA-256 guard, gene universe, score files, truth counts and n_eff, power class, AP-implementation identity
#   score    pooled AUPR (ap_fast), paired delete-half gene jackknife (500, seed 20261117), complex-level sensitivity (500, seed 20261119),
#            per-complex AP, contributions, slopes, hub-stratified read-outs; every method and read-out present in --scores
#   expjack  experiment jackknife (100 replicates, seed 20261118): CLR, aCLR, tcaCLR, |r| (PIDC is not refitted under draft 3b; --pidc_refit 1 is kept for experiments) refitted on
#            each half of the profiles in --expr (replicate b: seed 20261118 + b - 1, draft 3b); --lane/--nlanes split the replicates; resumable (one file per replicate)
#   collect  gathers the replicate files into expjack/replicates_<name>.csv
#   report   estimate.csv, h3.csv, pairwise.csv, h2.csv, secondaries.csv, readouts.csv, descriptors, fixed sentences (sentences.txt)
# Arm B: same estimand, intervals and descriptors; no p-value, no intersection-union / Holm statement, no claim sentence.
# --validation 1 permits non-registered seeds / replicate counts / truth order (consumed-data validation only; every claim sentence is
# then labelled as a validation). Registered values are TC in tc_lib.R.
suppressPackageStartupMessages(source("analysis/tcaclr_confirm/tc_lib.R"))
opt <- tc_options(commandArgs(TRUE), list(stage = "", arm = "A", scores = "", truth_dir = "", members = "", excluded = "", sha = "", counts = "", perturbed = "",
  out = "", expr = "", exp_units = "", exp_name = "plan", exp_scheme = "plan", exp_primary = "", exp_import = "", exp_seed = TC$exp_seed, exp_n = TC$exp_n,
  jack = TC$jack_n, jack_seed = TC$jack_seed, cjack = TC$cjack_n, cjack_seed = TC$cjack_seed, order = "plan", validation = 0L, allow_missing = 0L,
  workers = 8L, threads = 8L, lane = 0L, nlanes = 1L, pidc_refit = 0L, hash_refit = 0L, upstream_gates = "", smoke_ref = "", claim_label = "",
  context = "", kind = "", fixed_sentence = "", overwrite = 0L, names_new = "tcaCLR", names_ref = "aCLR", quick = 0L, neff_override = -1))
VAL <- opt$validation == 1L; ARM <- opt$arm; stopifnot(ARM %in% c("A", "B"))
if (!nzchar(opt$out)) stop("--out is required")
G_HEAD <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) NA_character_)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(opt$out, sprintf("run_%s.log", opt$stage)); say <- function(...) tc_say(LOG, ...)
nc <- max(1L, min(opt$workers, parallel::detectCores() - 1L))
SC <- file.path(opt$out, "score"); EJ <- file.path(opt$out, "expjack"); dir.create(SC, showWarnings = FALSE); dir.create(EJ, showWarnings = FALSE)

# ---- registered constants: refuse anything else unless --validation 1 -------------------------------------------------------------
if (!VAL) {
  bad <- c(opt$jack != TC$jack_n, opt$jack_seed != TC$jack_seed, opt$cjack != TC$cjack_n, opt$cjack_seed != TC$cjack_seed, opt$exp_n != TC$exp_n,
           opt$exp_seed != TC$exp_seed, opt$order != "plan", opt$exp_scheme != "plan", opt$quick != 0L, opt$neff_override != -1)
  if (any(bad)) stop("production mode refuses non-registered settings (--jack, --jack_seed, --cjack, --cjack_seed, --exp_n, --exp_seed, --order, --exp_scheme, --quick, --neff_override); use --validation 1 on consumed data")
}
# ---- the guard, in front of everything --------------------------------------------------------------------------------------------
guard <- tc_hash_guard(opt$sha, must_cover = c(opt$members, opt$excluded))
if (nzchar(opt$upstream_gates)) { ug <- utils::read.csv(opt$upstream_gates); if (!"pass" %in% names(ug) || !all(as.logical(ug$pass))) stop("upstream gate file reports a failure: ", opt$upstream_gates) }
say("stage ", opt$stage, "; arm ", ARM, "; mode ", if (VAL) "VALIDATION (consumed data)" else "production", "; git ", G_HEAD, "; hash guard passed for ", nrow(guard), " files")
st <- tc_scores_open(opt$scores); genes <- st$genes; G <- st$G
tr <- tc_read_truth(opt$members, opt$excluded)
perturbed <- NULL; if (nzchar(opt$perturbed)) { pg <- tc_perturbed_symbols(readLines(opt$perturbed)); pg <- pg[nzchar(pg)]; perturbed <- stats::na.omit(match(pg, genes)); say("perturbed genes (symbols) in the universe: ", length(perturbed), " of ", length(pg)) }
core_present <- intersect(TC$methods, st$methods)
if (!all(c("tcaCLR", "aCLR") %in% core_present)) stop("tcaCLR and aCLR score files are required")
if (!VAL && opt$allow_missing == 0L && length(core_present) < length(TC$methods)) stop("production mode needs all seven methods; missing: ", paste(setdiff(TC$methods, core_present), collapse = ", "))
rw <- grep("^rewire_[0-9]+$", st$methods, value = TRUE)
readouts_present <- c(intersect(TC$readouts, st$methods), sort(rw))
hub_ok <- "mi_rowmean" %in% st$methods
score_methods <- c(core_present, readouts_present)
build_truths <- function() tc_build_truths(tr, genes, opt$order, perturbed)

## ======================================================================================================================================
if (opt$stage == "gates") {
  g <- list(); add <- function(gate, what, ok, value = "") g[[length(g) + 1L]] <<- data.frame(gate = gate, check = what, value = as.character(value), pass = ok, stringsAsFactors = FALSE)
  for (k in seq_len(nrow(guard))) add("1 hash", guard$file[k], guard$ok[k], guard$got[k])
  add("1 universe", "genes in score directory", TRUE, G)
  add("1 universe", "no duplicated gene names", !anyDuplicated(genes), "")
  T <- build_truths(); add("2 truth", "all member genes inside the score universe (tc_sets stops otherwise)", TRUE, length(T))
  tc <- do.call(rbind, lapply(names(T), function(v) { p <- T[[v]]; data.frame(truth = v, role = p$role, n_complexes = length(p$sets), n_genes_in = length(p$genes_in), n_pairs = length(p$label), n_pos = sum(p$label),
    base_rate = mean(p$label), n_pairs_removed = p$n_removed, structure_md5 = tc_truth_struct_hash(p), note = p$note, stringsAsFactors = FALSE) }))
  if (nzchar(opt$counts)) { cn <- utils::read.csv(opt$counts, stringsAsFactors = FALSE)
    for (v in intersect(tc$truth, cn$variant)) { r <- cn[cn$variant == v, ]   # size strata: the freeze counts pairs inside the stratum's own genes, the scored stratum keeps the main universe; positives only
      ok <- (v %in% c("small_3to6", "mid_7to19", "large_ge20", "mid_7to10", "mid_11to19") || tc$n_pairs[tc$truth == v] == r$n_pairs[1]) && tc$n_pos[tc$truth == v] == r$n_pos[1]
      add("2 truth", paste0("counts match the freeze: ", v), ok, sprintf("%d/%d pairs, %d/%d pos", tc$n_pairs[tc$truth == v], r$n_pairs[1], tc$n_pos[tc$truth == v], r$n_pos[1])) } }
  pc <- tc_pair_complex(T$main$sets, T$main, G); neff_frac <- tc_neff(pc, sum(T$main$label)); neff <- tc_neff_overlap(pc); if (VAL && opt$neff_override > 0) neff <- opt$neff_override
  cls <- if (ARM == "B") "arm_B_estimate" else if (neff >= TC$neff_test) "test_0.03" else if (neff >= TC$neff_degraded) "test_0.04" else "estimate_only"
  add("2 truth", "n_eff of the main truth, overlap-counted (the gate's definition)", TRUE, sprintf("%.3f", neff)); add("2 truth", "n_eff of the main truth, fractional assignment (descriptive)", TRUE, sprintf("%.3f", neff_frac)); add("2 truth", "power class", TRUE, cls)
  # score files: present, right length, finite
  for (m in score_methods) { v <- tc_scores_packed(st, m); add("5 score", paste0(m, ": length and finiteness"), all(is.finite(v)) || m %in% c("GENIE3", "GRNBoost2", "PIDC"), sprintf("nonfinite %d", sum(!is.finite(v)))); rm(v) }
  # the sort-once replicate AP equals ap_fast on the subset (bit-identical)
  set.seed(1L); p <- T$main; v <- tc_scores_packed(st, "tcaCLR"); s <- v[tc_pair_index(p$i, p$j)]; pr <- tc_ap_prep(s, p$label); worst <- abs(tc_ap_sub(pr, seq_along(s)) - ap_fast(s, p$label))
  for (r in 1:3) { keep <- sample(p$genes_in, floor(length(p$genes_in) / 2)); inh <- logical(G); inh[keep] <- TRUE; ix0 <- which(inh[p$i] & inh[p$j]); ixs <- which(inh[p$i[pr$o]] & inh[p$j[pr$o]])
    worst <- max(worst, abs(tc_ap_sub(pr, ixs) - ap_fast(s[ix0], p$label[ix0]))) }
  set.seed(NULL); add("5 score", "sort-once replicate AP equals ap_fast (max abs difference, tolerance 1e-15)", worst <= 1e-15, format(worst))
  gt <- do.call(rbind, g); tc_wcsv(gt, file.path(opt$out, "gates.csv")); tc_wcsv(tc, file.path(opt$out, "truth_counts_gates.csv"))
  print(gt, row.names = FALSE); if (!all(gt$pass)) stop("a gate failed: nothing is scored")
  say("gates passed; power class ", cls); quit(save = "no")
}

## ======================================================================================================================================
if (opt$stage == "score") {
  if (!VAL && file.exists(file.path(SC, "score_complete.flag")) && opt$overwrite == 0L) stop("guard: the score stage has already run in ", SC, " (one run only; --overwrite 1 to repeat)")
  unlink(file.path(SC, "score_complete.flag"))
  T <- build_truths(); tn <- names(T); genes_main <- T$main$genes_in
  halves <- tc_draw_halves(T, opt$jack, opt$jack_seed)
  say("truths: ", paste(sprintf("%s (%d pairs, %d pos)", tn, vapply(T, function(p) length(p$label), 0L), vapply(T, function(p) sum(p$label), 0L)), collapse = "; "))
  idx <- lapply(T, function(p) tc_pair_index(p$i, p$j))
  pc <- tc_pair_complex(T$main$sets, T$main, G); neff_frac <- tc_neff(pc, sum(T$main$label)); neff <- tc_neff_overlap(pc); say(sprintf("main truth: n_eff %.3f overlap-counted (gate definition), %.3f fractional", neff, neff_frac))
  cnt <- do.call(rbind, lapply(tn, function(v) { p <- T[[v]]; data.frame(truth = v, role = p$role, n_complexes = length(p$sets), n_genes_in = length(p$genes_in), n_pairs = length(p$label), n_pos = sum(p$label),
    base_rate = mean(p$label), n_pairs_removed = p$n_removed, structure_md5 = tc_truth_struct_hash(p), note = p$note, stringsAsFactors = FALSE) }))
  opm <- tc_paralog_map(tr$excl); cmain <- T$main$sets; nC <- length(cmain); size_c <- lengths(cmain)
  regm <- list(members = cmain, class = stats::setNames(rep("complex", nC), names(cmain))); co <- comember_matrix(regm, G)
  met <- list(); rep_l <- list(); per_complex <- list(); cscore <- list(); main_score <- list(); pv_cj <- list()
  record <- function(m, scope, S_or_v, is_packed) {
    for (v in tn) { p <- T[[v]]; s <- if (is_packed) S_or_v[idx[[v]]] else S_or_v[cbind(p$i, p$j)]
      ap <- ap_fast(s, p$label); jr <- tc_jack_method(s, p, halves[[v]], G, nc)
      met[[length(met) + 1L]] <<- data.frame(truth = v, role = p$role, scope = scope, method = m, aupr = ap, n_pairs = length(p$label), n_pos = sum(p$label), base_rate = mean(p$label), stringsAsFactors = FALSE)
      rep_l[[length(rep_l) + 1L]] <<- data.frame(truth = v, scope = scope, method = m, rep = seq_along(jr), aupr = jr, stringsAsFactors = FALSE)
      if (scope == "raw" && v == "main") main_score[[m]] <<- s }
  }
  mi_rm <- if (hub_ok) tc_read_bin(st$files[["mi_rowmean"]]) else NULL
  for (m in score_methods) { t1 <- Sys.time(); invisible(gc(verbose = FALSE)); v <- tc_scores_packed(st, m)
    record(m, "raw", v, TRUE)
    if (m %in% c("tcaCLR", "aCLR", "CLR")) pv_cj[[m]] <- v
    if (m %in% core_present || (hub_ok && m %in% c("tcaCLR", "aCLR", "CLR"))) { S <- tc_unpack(v, G)
      if (m %in% core_present) {
        x <- regulon_ap_table(S, regm, co, opm, genes, min_pos = TC$min_pos_complex); per_complex[[m]] <- cbind(method = m, x)
        cscore[[m]] <- data.frame(method = m, complex = names(cmain), size = size_c, mean_within_score = tc_complex_mean_score(S, cmain, genes, opm), stringsAsFactors = FALSE) }
      if (hub_ok && m %in% c("tcaCLR", "aCLR", "CLR")) { Hh <- outer(mi_rm, mi_rm, "+"); Sh <- rz_hub_stratify(S, Hh); rm(Hh); record(m, "hub_mi", Sh, FALSE); rm(Sh) }
      rm(S) }
    rm(v); say(sprintf("method %-12s scored on %d truths with %d jackknife replicates (%.0f s); RSS peak %.2f GB", m, length(tn), opt$jack, as.numeric(difftime(Sys.time(), t1, units = "secs")), tc_vmhwm_gb())) }
  rm(co); met <- do.call(rbind, met); repdf <- do.call(rbind, rep_l); pcdf <- do.call(rbind, per_complex); csdf <- do.call(rbind, cscore)
  ## per-complex contributions to d = AP(tcaCLR) - AP(comparator) on the main truth; descriptors (tcaCLR minus each comparator)
  contrib <- list(); desc <- list()
  for (cm in setdiff(core_present, "tcaCLR")) { r <- tc_complex_contrib(main_score$tcaCLR, main_score[[cm]], T$main, cmain, G, pc); ds <- tc_descriptors(r$d_c, size_c, r$d)
    contrib[[cm]] <- data.frame(comparator = cm, complex = names(cmain), size = size_c, d_c = unname(r$d_c), stringsAsFactors = FALSE)
    desc[[cm]] <- data.frame(comparator = cm, d = r$d, share_top3 = ds$share_top3, top3 = ds$top3, share_large = ds$share_large, k_carry = ds$k_carry, stringsAsFactors = FALSE) }
  contrib <- do.call(rbind, contrib); desc <- do.call(rbind, desc)
  ## slope of the mean within-complex score on log size, per method (H2 continuous form), raw and scaled by the SD of the method's pair scores
  sl <- do.call(rbind, lapply(core_present, function(m) { x <- csdf[csdf$method == m, ]; s <- tc_slope(x$mean_within_score, x$size); sdm <- stats::sd(main_score[[m]])
    data.frame(method = m, slope_raw = s[["slope"]], se_raw = s[["se"]], pair_score_sd = sdm, slope_scaled = s[["slope"]] / sdm, n_complexes = sum(is.finite(x$mean_within_score)), stringsAsFactors = FALSE) }))
  ## complex-level delete-half sensitivity (tcaCLR, aCLR, CLR; main truth)
  k_all <- tc_excl_keys(tr$excl, "all", genes); cjm <- intersect(c("tcaCLR", "aCLR", "CLR"), names(pv_cj))
  set.seed(opt$cjack_seed); chalves <- lapply(seq_len(opt$cjack), function(r) sort(sample.int(nC, floor(nC / 2)))); set.seed(NULL)
  cone <- function(keep) { p <- tc_pairs(cmain[keep], genes, k_all); ix <- tc_pair_index(p$i, p$j)
    r <- c(vapply(cjm, function(m) ap_fast(pv_cj[[m]][ix], p$label), 0), n_pairs = length(p$label), n_pos = sum(p$label)); invisible(gc(FALSE, full = FALSE)); r }
  invisible(gc(verbose = FALSE)); t1 <- Sys.time()
  cj <- if (nc > 1L) parallel::mclapply(chalves, cone, mc.cores = nc) else lapply(chalves, cone)
  if (any(vapply(cj, function(x) inherits(x, "try-error") || length(x) != length(cjm) + 2L, NA))) stop("a complex-jackknife worker failed")
  cj <- do.call(rbind, cj); say(sprintf("complex-level jackknife: %d replicates (%.0f s)", nrow(cj), as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  cjl <- do.call(rbind, lapply(cjm, function(m) data.frame(rep = seq_len(nrow(cj)), method = m, aupr = cj[, m], n_pairs = cj[, "n_pairs"], n_pos = cj[, "n_pos"])))
  ## write (staging, then rename)
  outs <- list(truth_counts.csv = cnt, metrics.csv = met, jack_gene_replicates.csv = repdf, jack_complex_replicates.csv = cjl, per_complex_ap.csv = pcdf, complex_scores.csv = csdf,
               complex_contrib.csv = contrib, descriptors.csv = desc, h2_slopes.csv = sl, neff.csv = data.frame(truth = "main", n_eff = neff, n_eff_fractional = neff_frac, n_pos = sum(T$main$label), n_complexes = nC, n_truth_genes = length(genes_main),
                                   n_truth_genes_perturbed = if (is.null(perturbed)) NA_integer_ else length(intersect(genes_main, perturbed))))
  for (f in names(outs)) tc_wcsv(outs[[f]], file.path(SC, f))
  if (nzchar(opt$smoke_ref)) { ref <- utils::read.csv(opt$smoke_ref, stringsAsFactors = FALSE); mm <- merge(met[met$scope == "raw", ], ref, by = c("truth", "method"), suffixes = c("", "_ref"))
    dd <- max(abs(mm$aupr - mm$aupr_ref)); tc_wcsv(data.frame(gate = "5 score: AUPR against the registered reference (tolerance 1e-12)", n = nrow(mm), max_abs_diff = dd, pass = dd <= 1e-12), file.path(SC, "smoke_gate.csv"))
    if (!is.finite(dd) || dd > 1e-12) stop("score gate failed: AUPR differs from the registered reference by ", format(dd)) }
  writeLines(c(format(Sys.time()), G_HEAD), file.path(SC, "score_complete.flag"))
  print(met[met$truth == "main" & met$scope == "raw", c("method", "aupr")], row.names = FALSE)
  say(sprintf("SCORE STAGE done; peak RSS %.2f GB", tc_vmhwm_gb())); quit(save = "no")
}

## ======================================================================================================================================
if (opt$stage %in% c("expjack", "collect")) {
  rd <- file.path(EJ, opt$exp_name, "reps"); dir.create(rd, recursive = TRUE, showWarnings = FALSE)
  repf <- function(b) file.path(rd, sprintf("rep%03d.csv", b)); hashf <- function(b) file.path(rd, sprintf("rep%03d_hash.csv", b))
  if (opt$stage == "collect") {
    have <- file.exists(vapply(seq_len(opt$exp_n), repf, "")); say("collect: ", sum(have), " of ", opt$exp_n, " replicates present")
    if (!all(have)) { say("MISSING: ", paste(which(!have), collapse = " ")); quit(save = "no", status = 2L) }
    x <- do.call(rbind, lapply(seq_len(opt$exp_n), function(b) utils::read.csv(repf(b), stringsAsFactors = FALSE))); tc_wcsv(x, file.path(EJ, paste0("replicates_", opt$exp_name, ".csv")))
    hh <- file.exists(vapply(seq_len(opt$exp_n), hashf, "")); if (all(hh)) tc_wcsv(do.call(rbind, lapply(seq_len(opt$exp_n), function(b) utils::read.csv(hashf(b), stringsAsFactors = FALSE))), file.path(EJ, paste0("refit_hashes_", opt$exp_name, ".csv")))
    say("collect done"); quit(save = "no")
  }
  if (!nzchar(opt$expr)) stop("--expr is required for the expjack stage")
  X <- if (grepl("\\.rds$", opt$expr)) readRDS(opt$expr) else { d <- utils::read.delim(opt$expr, check.names = FALSE, stringsAsFactors = FALSE); m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- d[[1]]; m }
  storage.mode(X) <- "double"; xrow <- match(genes, rownames(X)); if (anyNA(xrow)) stop("score-universe genes missing from the expression matrix: ", sum(is.na(xrow)))
  if (anyDuplicated(rownames(X))) stop("duplicated gene names in the expression matrix")
  units <- if (nzchar(opt$exp_units)) { u <- utils::read.delim(opt$exp_units, stringsAsFactors = FALSE); if (nrow(u) != ncol(X)) stop("--exp_units must have one row per column of --expr"); u$unit } else seq_len(ncol(X))
  ids <- sort(unique(units)); halves <- tc_exp_halves(length(ids), opt$exp_n, opt$exp_seed, opt$exp_scheme)
  T <- build_truths(); tn <- names(T)
  if (file.exists(file.path(SC, "truth_counts.csv"))) { tc <- utils::read.csv(file.path(SC, "truth_counts.csv"), stringsAsFactors = FALSE)
    if (!all(vapply(tn, function(v) identical(tc$structure_md5[tc$truth == v], tc_truth_struct_hash(T[[v]])), NA))) stop("the truth built here differs from the score-stage truth (truth_counts.csv)")
    say("truth structure identical to the score stage (md5 of i, j, label for every truth)") }
  mine <- which((seq_len(opt$exp_n) - 1L) %% opt$nlanes == opt$lane); say(sprintf("lane %d of %d: %d replicates; %d units of %d profiles; scheme %s; git %s", opt$lane, opt$nlanes, length(mine), length(ids), ncol(X), opt$exp_scheme, G_HEAD))
  for (b in mine) { if (file.exists(repf(b))) next
    t0 <- Sys.time(); cols <- which(units %in% ids[halves[[b]]]); Xb <- X[, cols, drop = FALSE]
    if (any(apply(Xb, 1, function(x) diff(range(x)) == 0))) stop("replicate ", b, ": zero-variance genes in the kept half")
    R <- tc_fit_three(Xb, opt$threads); R$absr <- tc_absr(Xb)
    hs <- if (opt$hash_refit == 1L) vapply(c("CLR", "aCLR", "tcaCLR"), function(m) { con <- rawConnection(raw(0), "wb"); writeBin(tc_pack(R[[m]]), con, size = 8L); h <- digest::digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE); close(con); h }, "") else NULL
    pidc <- if (opt$pidc_refit == 1L) tc_pidc_refit(Xb, opt$threads) else NULL
    rows <- list()
    for (v in tn) { p <- T[[v]]; ii <- xrow[p$i]; jj <- xrow[p$j]; ix <- cbind(ii, jj)
      for (m in names(R)) rows[[length(rows) + 1L]] <- data.frame(rep = b, truth = v, method = m, aupr = ap_fast(R[[m]][ix], p$label), stringsAsFactors = FALSE)
      if (!is.null(pidc)) rows[[length(rows) + 1L]] <- data.frame(rep = b, truth = v, method = "PIDC", aupr = ap_fast(pidc(ii, jj), p$label), stringsAsFactors = FALSE) }
    if (!is.null(hs)) tc_wcsv(data.frame(rep = b, method = names(hs), sha256 = unname(hs), stringsAsFactors = FALSE), hashf(b))
    tc_wcsv(do.call(rbind, rows), repf(b)); secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); rm(R, Xb); invisible(gc(verbose = FALSE))
    say(sprintf("%s %3d done (%d profiles kept, %.0f s, peak RSS %.2f GB)", opt$exp_name, b, length(cols), secs, tc_vmhwm_gb())) }
  say("lane done"); quit(save = "no")
}

## ======================================================================================================================================
if (opt$stage == "report") {
  est_f <- file.path(opt$out, "estimate.csv")
  if (file.exists(est_f) && opt$overwrite == 0L && !VAL) stop("guard: estimate.csv exists (written once; --overwrite 1 to repeat)")
  rd <- function(f) utils::read.csv(file.path(SC, f), stringsAsFactors = FALSE)
  met <- rd("metrics.csv"); jg <- rd("jack_gene_replicates.csv"); cj <- rd("jack_complex_replicates.csv"); desc <- rd("descriptors.csv"); pcdf <- rd("per_complex_ap.csv")
  contrib <- rd("complex_contrib.csv"); sl <- rd("h2_slopes.csv"); neffd <- rd("neff.csv"); neff <- neffd$n_eff[1]; if (VAL && opt$neff_override > 0) neff <- opt$neff_override
  cls <- if (ARM == "B") "arm_B_estimate" else if (neff >= TC$neff_test) "test_0.03" else if (neff >= TC$neff_degraded) "test_0.04" else "estimate_only"
  p_allowed <- ARM == "A" && cls != "estimate_only"
  ## experiment sources: expjack/replicates_<name>.csv, or imports "name=path:map,..." (map net4 renames clr/attention/a3r)
  srcs <- list()
  for (f in list.files(EJ, pattern = "^replicates_.*\\.csv$")) { nm <- sub("^replicates_(.*)\\.csv$", "\\1", f); srcs[[nm]] <- utils::read.csv(file.path(EJ, f), stringsAsFactors = FALSE) }
  if (nzchar(opt$exp_import)) for (it in strsplit(opt$exp_import, ",")[[1]]) { nm <- sub("=.*$", "", it); path <- sub("^[^=]*=", "", it); x <- utils::read.csv(path, stringsAsFactors = FALSE)
    x$method <- c(clr = "CLR", attention = "aCLR", a3r = "tcaCLR")[x$readout]; x <- x[!is.na(x$method), c("rep", "truth", "method", "aupr")]; srcs[[nm]] <- x }
  prim <- if (nzchar(opt$exp_primary)) opt$exp_primary else if (length(srcs)) names(srcs)[1] else ""
  if (nzchar(prim) && !prim %in% names(srcs)) stop("primary experiment source not found: ", prim)
  p_allowed <- p_allowed && nzchar(prim)     # no experiment component, no p-value
  say("experiment sources: ", if (length(srcs)) paste(names(srcs), collapse = ", ") else "none", "; primary: ", if (nzchar(prim)) prim else "none")
  aupr <- function(tr_, sc_, m) { x <- met$aupr[met$truth == tr_ & met$scope == sc_ & met$method == m]; if (length(x)) x[1] else NA_real_ }
  greps <- function(tr_, sc_, m) { x <- jg[jg$truth == tr_ & jg$scope == sc_ & jg$method == m, ]; x$aupr[order(x$rep)] }
  ereps <- function(src, tr_, m) { x <- srcs[[src]]; x <- x[x$truth == tr_ & x$method == m, ]; x$aupr[order(x$rep)] }
  contrast <- function(tr_, sc_, a, b, src = prim) {
    d <- aupr(tr_, sc_, a) - aupr(tr_, sc_, b); ra <- greps(tr_, sc_, a); rb <- greps(tr_, sc_, b)
    gi <- tc_interval(d, ra - rb)
    se_e <- NA_real_; ne <- 0L
    if (sc_ == "raw" && nzchar(src) && !(a %in% TC$no_refit) && !(b %in% TC$no_refit)) { ea <- ereps(src, tr_, a); eb <- ereps(src, tr_, b); if (length(ea) > 1L && length(ea) == length(eb)) { se_e <- stats::sd(ea - eb); ne <- length(ea) } }
    tt <- tc_total(d, gi$se_gene, se_e)
    needs_exp <- !(a %in% TC$no_refit) && !(b %in% TC$no_refit)     # CLR, aCLR, tcaCLR, |r|, PIDC are refitted; GENIE3 and GRNBoost2 are gene-only by registration
    p_ok <- p_allowed && sc_ == "raw" && (tt$exp_included || !needs_exp)
    data.frame(truth = tr_, scope = sc_, method_a = a, method_b = b, contrast = paste(tc_name(a), "-", tc_name(b)), d = d, se_gene = gi$se_gene, se_exp = se_e, n_exp_reps = ne, se_total = tt$se_total,
               exp_included = tt$exp_included, total_lo = tt$total_lo, total_hi = tt$total_hi, gene_lo = gi$gene_lo, gene_hi = gi$gene_hi, gene_pct_lo = gi$pct_lo, gene_pct_hi = gi$pct_hi,
               gene_skew = gi$skew, n_gene_reps = gi$n_rep, z = tt$z, p_one_sided = if (p_ok) tt$p_one_sided else NA_real_, exp_missing = needs_exp && !tt$exp_included, stringsAsFactors = FALSE)
  }
  have_m <- unique(met$method[met$scope == "raw"]); core <- intersect(TC$methods, have_m); comps <- intersect(TC$comparators, core)
  all_truths <- unique(met$truth)
  ## primary and the registered rows
  h1 <- contrast("main", "raw", "tcaCLR", "aCLR")
  h3rows <- do.call(rbind, lapply(comps, function(cm) contrast("main", "raw", "tcaCLR", cm)))
  pair_m <- core; prs <- utils::combn(pair_m, 2); pairwise <- do.call(rbind, lapply(seq_len(ncol(prs)), function(k) contrast("main", "raw", prs[1, k], prs[2, k])))
  pairwise$note <- ifelse(pairwise$method_a %in% TC$no_refit | pairwise$method_b %in% TC$no_refit, TC$note_ge, "")
  ## every truth and scope: tcaCLR - aCLR (P), tcaCLR - CLR (C1), aCLR - CLR (C2), gene interval and, where refitted, the experiment component
  ctr <- do.call(rbind, lapply(unique(met$scope), function(sc_) do.call(rbind, lapply(all_truths, function(v) do.call(rbind, lapply(list(c("P", "tcaCLR", "aCLR"), c("C1", "tcaCLR", "CLR"), c("C2", "aCLR", "CLR")), function(ab)
    if (all(ab[2:3] %in% met$method[met$scope == sc_ & met$truth == v])) cbind(row = ab[1], contrast(v, sc_, ab[2], ab[3])) ))))))
  ## alternative experiment sources for the primary contrast (reported side by side)
  alt <- if (length(srcs)) do.call(rbind, lapply(names(srcs), function(s) do.call(rbind, lapply(list(c("P", "tcaCLR", "aCLR"), c("C1", "tcaCLR", "CLR"), c("C2", "aCLR", "CLR")), function(ab)
    if (all(ab[2:3] %in% unique(srcs[[s]]$method))) { x <- cbind(row = ab[1], contrast("main", "raw", ab[2], ab[3], s)); x$exp_source <- s; x$is_primary <- s == prim; x })))) else NULL
  ## complex-level sensitivity (labelled; never decides)
  pa <- stats::setNames(met$aupr[met$truth == "main" & met$scope == "raw"], met$method[met$truth == "main" & met$scope == "raw"])
  cjrow <- function(a, b) { ra <- cj$aupr[cj$method == a][order(cj$rep[cj$method == a])]; rb <- cj$aupr[cj$method == b][order(cj$rep[cj$method == b])]; i <- tc_interval(pa[[a]] - pa[[b]], ra - rb)
    data.frame(method_a = a, method_b = b, d = i$d, se_complex = i$se_gene, cl_lo = i$gene_lo, cl_hi = i$gene_hi, cl_pct_lo = i$pct_lo, cl_pct_hi = i$pct_hi, cl_skew = i$skew, n_rep = i$n_rep, stringsAsFactors = FALSE) }
  cjm <- unique(cj$method); cjt <- do.call(rbind, lapply(list(c("tcaCLR", "aCLR"), c("tcaCLR", "CLR"), c("aCLR", "CLR")), function(ab) if (all(ab %in% cjm)) cjrow(ab[1], ab[2])))
  ## descriptors for H1 (tcaCLR minus aCLR)
  dh <- desc[desc$comparator == "aCLR", ]
  pa_ap <- pcdf[pcdf$method == "tcaCLR", ]; pb_ap <- pcdf[pcdf$method == "aCLR", ]; ord <- match(pa_ap$regulon, pb_ap$regulon)
  br <- tc_breadth(pa_ap$ap, pb_ap$ap[ord]); br_sent <- tc_breadth_sentence(br, dh$k_carry, dh$d)
  ## secondaries and size classes
  strata_names <- intersect(c("small_3to6", "mid_7to10", "mid_11to19", "mid_7to19", "large_ge20"), all_truths)
  sec_names <- intersect(c("noribo", "large_ge20", "small_3to6", "main_xsub_excl", "main_ribo_joined", "pert_neither", "pert_rest", "nomerge", "size2", "main_paralog_incl"), all_truths)
  sec <- do.call(rbind, lapply(unique(c(strata_names, sec_names)), function(v) contrast(v, "raw", "tcaCLR", "aCLR")))
  sec$basis <- ifelse(sec$exp_included, "gene + experiment", "gene only"); sec$lo <- ifelse(sec$exp_included, sec$total_lo, sec$gene_lo); sec$hi <- ifelse(sec$exp_included, sec$total_hi, sec$gene_hi)
  sec$sentence <- vapply(seq_len(nrow(sec)), function(k) if (anyNA(c(sec$d[k], sec$lo[k], sec$hi[k]))) sprintf("%s: interval not evaluable (too few replicates with a positive pair).", sec$truth[k]) else tc_secondary_sentence(sec$truth[k], sec$d[k], sec$lo[k], sec$hi[k], sec$basis[k]), "")
  ## H2
  lab_of <- function(v) { x <- sec[sec$truth == v, ]; if (!nrow(x)) return(NULL); lo <- if (x$exp_included) x$total_lo else x$gene_lo; hi <- if (x$exp_included) x$total_hi else x$gene_hi; if (anyNA(c(x$d, lo, hi))) return(NULL); tc_h2_labels(x$d, lo, hi) }
  h2 <- NULL; h2_sent <- "H2: not evaluable (a size class is missing)."
  ls_ <- lab_of("small_3to6"); lm_ <- lab_of("mid_7to19"); ll_ <- lab_of("large_ge20")
  if (!is.null(ls_) && !is.null(lm_) && !is.null(ll_)) { hv <- tc_h2_verdict(ls_, lm_, ll_); h2_sent <- hv$sentence }
  h2tab <- do.call(rbind, lapply(strata_names, function(v) { x <- sec[sec$truth == v, ]; l <- lab_of(v)
    data.frame(size_class = v, registered_in_prediction = v %in% c("small_3to6", "mid_7to19", "large_ge20"), d = x$d, se_gene = x$se_gene, se_exp = x$se_exp, lo = x$lo, hi = x$hi, basis = x$basis,
               label_positive = l[["positive"]], label_near_zero = l[["near_zero"]], label_not_positive = l[["not_positive"]], labels = paste(names(l)[l], collapse = "; "), stringsAsFactors = FALSE) }))
  ## H3
  h3 <- NULL; h3_sent <- NULL; ind_sent <- NULL
  if (nrow(h3rows)) {
    h3rows$note <- ifelse(h3rows$method_b %in% TC$no_refit, TC$note_ge, "")
    if (p_allowed) { hh <- tc_h3(data.frame(comparator = h3rows$method_b, d = h3rows$d, p = h3rows$p_one_sided, stringsAsFactors = FALSE))
      h3rows$won_iut <- hh$tab$won_iut; h3rows$holm_p <- hh$tab$holm_p; h3rows$won_holm <- hh$tab$won_holm
      lost <- tc_name(hh$tab$comparator[!hh$tab$won_iut])
      miss <- h3rows$method_b[h3rows$exp_missing]
      h3_sent <- if (length(miss)) sprintf("tcaCLR cannot be declared the best of the %d methods run: the experiment component is missing for %s, so those comparisons have no p-value.", hh$k, paste(tc_name(miss), collapse = ", "))
                 else if (hh$best) sprintf("tcaCLR is the best of the %d methods run (every one of the %d comparisons won at one-sided 0.025; intersection-union, no correction).", hh$k, nrow(hh$tab))
                 else sprintf("tcaCLR is not shown to be the best of the %d methods run: it does not exceed %s at one-sided 0.025 (intersection-union).", hh$k, paste(lost, collapse = ", "))
      if (any(h3rows$method_b %in% TC$no_refit)) h3_sent <- c(h3_sent, paste0("Note: ", TC$note_ge))
      ind_sent <- vapply(seq_len(nrow(h3rows)), function(k) if (h3rows$exp_missing[k]) sprintf("tcaCLR versus %s: no statement (experiment component missing).", tc_name(h3rows$method_b[k])) else if (h3rows$won_holm[k]) sprintf("tcaCLR exceeds %s (Holm-adjusted one-sided p %s, %d comparisons).", tc_name(h3rows$method_b[k]), tc_fmt_p(h3rows$holm_p[k]), nrow(h3rows)) else
                         sprintf("tcaCLR is not shown to exceed %s (Holm-adjusted one-sided p %s, %d comparisons).", tc_name(h3rows$method_b[k]), tc_fmt_p(h3rows$holm_p[k]), nrow(h3rows)), "")
      ind_sent <- c(ind_sent, if (any(h3rows$method_b %in% TC$no_refit)) paste0("Note: ", TC$note_ge))
    } else { h3rows$won_iut <- NA; h3rows$holm_p <- NA_real_; h3rows$won_holm <- NA }
    h3 <- h3rows }
  ## read-outs (intervals only, gene resampling only; none selected)
  ro_names <- intersect(c(TC$readouts, sort(grep("^rewire_[0-9]+$", have_m, value = TRUE))), have_m)
  ro <- NULL
  mkro <- function(label, ap_r, reps_r, ap_t, reps_t, ap_a, reps_a, scope) {
    one <- function(refap, refrep, refname) { i <- tc_interval(ap_r - refap, reps_r - refrep); data.frame(readout = label, scope = scope, aupr = ap_r, versus = refname, d = i$d, se_gene = i$se_gene, lo = i$gene_lo, hi = i$gene_hi, stringsAsFactors = FALSE) }
    rbind(one(ap_t, reps_t, "tcaCLR"), one(ap_a, reps_a, "aCLR")) }
  at <- aupr("main", "raw", "tcaCLR"); rt <- greps("main", "raw", "tcaCLR"); aa <- aupr("main", "raw", "aCLR"); ra <- greps("main", "raw", "aCLR")
  ro <- do.call(rbind, c(lapply(setdiff(ro_names, grep("^rewire_", ro_names, value = TRUE)), function(m) mkro(m, aupr("main", "raw", m), greps("main", "raw", m), at, rt, aa, ra, "raw")),
    lapply(grep("^rewire_", ro_names, value = TRUE), function(m) mkro(m, aupr("main", "raw", m), greps("main", "raw", m), at, rt, aa, ra, "raw")),
    if (length(grep("^rewire_", ro_names))) { rw_ <- grep("^rewire_", ro_names, value = TRUE); mr <- rowMeans(do.call(cbind, lapply(rw_, function(m) greps("main", "raw", m)))); ma <- mean(vapply(rw_, function(m) aupr("main", "raw", m), 0))
      list(mkro(sprintf("rewire_mean_of_%d", length(rw_)), ma, mr, at, rt, aa, ra, "raw")) },
    lapply(intersect(c("tcaCLR", "aCLR", "CLR"), unique(met$method[met$scope == "hub_mi"])), function(m) mkro(paste0("hub_stratified_", m), aupr("main", "hub_mi", m), greps("main", "hub_mi", m), at, rt, aa, ra, "hub_mi"))))
  if (!is.null(ro)) ro$basis <- "gene resampling only (read-outs are not refitted); intervals only, none may be cited in place of H1"
  ## sentences
  ctx <- if (nzchar(opt$context)) opt$context else "complexes"
  kind <- if (nzchar(opt$kind)) opt$kind else if (ARM == "A") "preregistered confirmatory test" else "preregistered secondary estimate; not a confirmatory test"
  fixed <- if (nzchar(opt$fixed_sentence)) opt$fixed_sentence else if (ARM == "A") "Arm A (primary): the claim rule of plan section 2 applies; the interval sentences above are descriptive of the same interval."
           else "Arm B (secondary estimate): no p-value and no claim rule apply to this arm (plan section 2)."
  Lh <- if (h1$exp_included) h1$total_lo else h1$gene_lo; Uh <- if (h1$exp_included) h1$total_hi else h1$gene_hi
  isent <- tc_interval_sentences(Lh, Uh, opt$names_new, opt$names_ref, ctx, kind = kind, fixed = fixed)
  flag <- if (!h1$exp_included) "experiment component missing; gene-only interval; no p-value, no claim" else ""
  if (nzchar(flag)) p_allowed <- FALSE
  claim <- NA_character_
  if (ARM == "A") {
    claim_txt <- if (nzchar(opt$claim_label)) opt$claim_label else if (VAL) "VALIDATION RUN on consumed data, not a claim: tcaCLR exceeds aCLR in pooled co-membership AUPR on the consumed truth" else TC$claim_a
    claim <- if (nzchar(flag)) paste0("The claim is not made: ", flag, "; ", TC$ecoli_fail, ".") else tc_claim_sentence(h1$p_one_sided, h1$d, cls, dh$share_large, claim_txt)
    if (VAL && !startsWith(claim, "VALIDATION RUN")) claim <- paste0("VALIDATION RUN on consumed data, not a claim: ", claim)
  }
  pert_sent <- if (all(c("pert_neither", "pert_rest") %in% sec$truth)) { a <- sec[sec$truth == "pert_neither", ]; b <- sec[sec$truth == "pert_rest", ]; if (anyNA(c(a$lo, b$lo))) "Perturbed-gene split: interval not evaluable." else tc_pert_sentence(a$d, a$lo, a$hi, b$d, b$lo, b$hi) } else NA_character_
  ## estimate.csv: row H1 with the mandatory descriptors, then the H3 rows
  cjp <- cjt[cjt$method_a == "tcaCLR" & cjt$method_b == "aCLR", ]
  h1row <- cbind(data.frame(arm = ARM, mode = if (VAL) "validation" else "production", row = "H1", power_class = cls, n_eff = neff, exp_source = prim, flag = flag, stringsAsFactors = FALSE), h1,
    data.frame(complex_se = cjp$se_complex, complex_lo = cjp$cl_lo, complex_hi = cjp$cl_hi, share_top3 = dh$share_top3, top3_complexes = dh$top3, share_ge20 = dh$share_large, k_carry = dh$k_carry,
               n_complexes_ge10pos = br$n_complexes, n_improved = br$n_improved, frac_improved = br$frac_improved, sign_test_p = br$p_sign_one_sided, breadth_pass = br$pass,
               replicate_skew = h1$gene_skew, stringsAsFactors = FALSE),
    data.frame(S_zero = isent[["S_zero"]], S_min = isent[["S_min"]], S_ecoli = isent[["S_ecoli"]], S_fixed = isent[["S_fixed"]], claim = claim, breadth_sentence = br_sent, h2_sentence = h2_sent, git_head = G_HEAD, stringsAsFactors = FALSE))
  tc_wcsv(h1row, est_f)
  if (!is.null(h3)) tc_wcsv(h3, file.path(opt$out, "h3.csv"))
  tc_wcsv(pairwise, file.path(opt$out, "pairwise.csv")); tc_wcsv(ctr, file.path(opt$out, "contrasts.csv")); tc_wcsv(sec, file.path(opt$out, "secondaries.csv")); tc_wcsv(h2tab, file.path(opt$out, "h2.csv"))
  tc_wcsv(cjt, file.path(opt$out, "complex_level_sensitivity.csv")); if (!is.null(alt)) tc_wcsv(alt, file.path(opt$out, "estimate_by_experiment_source.csv"))
  if (!is.null(ro)) tc_wcsv(ro, file.path(opt$out, "readouts.csv")); tc_wcsv(br, file.path(opt$out, "breadth.csv"))
  ## the fixed sentences, in one file
  sl_txt <- c(paste0("# Fixed sentences (generated by analysis/tcaclr_confirm/score_confirm.R, arm ", ARM, ", ", if (VAL) "VALIDATION on consumed data" else "production", ", git ", G_HEAD, ")"), "",
    "## H1 (tcaCLR minus aCLR, main truth)", if (ARM == "A") c(claim, ""), paste0("- ", isent), "",
    "## Breadth (sign test over complexes with 10 or more positives)", paste0("- ", br_sent), "", "## H2 (size classes)", paste0("- ", h2_sent), "",
    if (!is.null(h3_sent)) c("## H3 (comparators)", paste0("- ", h3_sent), paste0("- ", ind_sent), paste0("- Trap: ", TC$trap), ""),
    if (is.null(h3_sent) && !is.null(h3)) c(paste0("## Comparators (estimates with intervals; ", if (ARM == "B") "Arm B: no p-value, no best-of statement" else "no p-value in this report", ")"), paste0("- Trap: ", TC$trap),
      if (any(h3$method_b %in% TC$no_refit)) paste0("- Note: ", TC$note_ge), ""),
    "## Secondaries (intervals and script-chosen sentences; no p-values)", paste0("- ", sec$sentence), if (!is.na(pert_sent)) paste0("- ", pert_sent), "",
    "## Read-outs of the same operator", "- Exploratory, one table (readouts.csv), intervals only; none may be cited in place of H1; nothing is selected among them.")
  sens <- switch(cls, test_0.03 = sprintf("Power gate: n_eff = %.1f (overlap-counted); the registered sensitivity is +0.03.", neff), test_0.04 = sprintf("Power gate: n_eff = %.1f (overlap-counted), between 15 and 27; the test runs with its sensitivity stated as +0.04.", neff),
                 estimate_only = sprintf("Power gate: n_eff = %.1f (overlap-counted), below 15; estimate only, the arm is not a test.", neff), sprintf("Arm B: n_eff = %.1f (overlap-counted); estimate only by registration.", neff))
  link <- if (is.null(neffd$n_truth_genes_perturbed) || is.na(neffd$n_truth_genes_perturbed[1])) NULL else sprintf("Built-in data-truth link: %d of the %d truth genes are themselves perturbed in the panel.", neffd$n_truth_genes_perturbed[1], neffd$n_truth_genes[1])
  sl_txt <- append(sl_txt, c("", "## Gates and disclosure", paste0("- ", sens), if (!is.null(link)) paste0("- ", link)))
  writeLines(sl_txt, file.path(opt$out, "sentences.txt"))
  print(h1row[, c("row", "d", "se_gene", "se_exp", "total_lo", "total_hi", "p_one_sided", "complex_lo", "complex_hi", "share_top3", "n_improved", "n_complexes_ge10pos")], row.names = FALSE)
  cat(sl_txt, sep = "\n"); say("REPORT STAGE done"); quit(save = "no")
}
stop("unknown --stage: ", opt$stage)
