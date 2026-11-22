. jobs_pc/common.sh
# mem: 4
# cpus: 1
# after: 7910_spectral_summary
# HELD: SEALED HELD-OUT DATASET. DO NOT SUBMIT until the lead has (1) read the development ratios, (2) written the FROZEN variant list into the ledger of docs/plans/spectral_embedding_pass_plan.md and committed it, and (3) set FROZEN below (comma-separated variant names). Scores Jurkat ONCE, frozen variants only (--heldout YES is required by sp_score.R, which also refuses to run without --variants). Needs the Jurkat scGPT export and results/scgpt_ext/Jurkat/ (comparison columns) to exist, which is the lead's own step. No development path opens Jurkat. Output results/spectral_pass/Jurkat/ and results/spectral_pass/summary_jurkat/.
set -e
FROZEN="cos_kE,cos_k10"   # frozen 22 Nov 2026, see docs/plans/spectral_embedding_pass_plan.md section 10
[ -n "$FROZEN" ] || { echo "FROZEN variant list not set: the Jurkat held-out job is not released" >&2; exit 1; }
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/spectral_pass/sp_score.R --datasets Jurkat --heldout YES --variants "$FROZEN" --out results/spectral_pass --truth real --threads 1
[ -f results/scgpt_ext/Jurkat/metrics_ext.csv ] && nice -n 19 Rscript analysis/spectral_pass/sp_summary.R --datasets Jurkat --out results/spectral_pass/summary_jurkat
