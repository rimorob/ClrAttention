# shared setup for the DREAM5 tests; run from the repo root (tests/dream5/run_tests.R)
setwd(Sys.getenv("D5_ROOT", unset = "."))   # repo root (set by run_tests.R)
Sys.setenv(OMP_NUM_THREADS = "2", OPENBLAS_NUM_THREADS = "1")
source("analysis/dream5_lib.R"); source("analysis/dream5_score_lib.R")
TESTSW <- function(...) { s <- d5_edge_switches(ties_at_cap = "keep_tied", ties_in_list = "end_of_block", tftf = "keep_both", universe = "all"); mod <- list(...); for (k in names(mod)) s[[k]] <- mod[[k]]; s }
