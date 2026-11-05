#!/usr/bin/env Rscript
# DREAM5 step-5 completeness: which required score files and refit hashes exist, and the disk used.
#   Rscript analysis/dream5_status.R            (read-only)
source("analysis/dream5_score_lib.R")
req <- rbind(d5_required(), data.frame(net = rep(D5$nets, each = 2), method = c("genie3_all", "grnboost2_all"), input = "replicate", stringsAsFactors = FALSE))
req$have <- file.exists(mapply(d5_cache_path, req$net, req$method, req$input)) & file.exists(paste0(mapply(d5_cache_path, req$net, req$method, req$input), ".sha256"))
req$kind <- ifelse(req$method %in% D5$primary, "required", "descriptive")
print(reshape(req[, c("net", "method", "have")], idvar = "method", timevar = "net", direction = "wide"), row.names = FALSE)
nrep <- length(Sys.glob(file.path(d5_refit_dir(4L), "rep*.sha256.tsv"))); cat(sprintf("net4 refit replicates with hashes: %d / %d\n", nrep, D5$n_refit))
cat("required files missing:", sum(!req$have & req$kind == "required"), "\n")
system("du -sh results/dream5")
