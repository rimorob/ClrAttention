. jobs_pc/common.sh
# mem: 2
# cpus: 1
# after: 7900_spectral_hesc_score,7901_spectral_hhep_score,7902_spectral_k562_score,7903_spectral_hepg2_score,7904_spectral_rpe1_score
# Spectral pass: descriptive ratio tables (variant / scGPT token cosine, / tcaCLR2, / scGPT attention, per truth and per dataset; no test, no threshold). Reads the metric tables of this pass and results/scgpt_ext/<dataset>/metrics_ext.csv for the five development datasets only. Output results/spectral_pass/summary/.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 19 Rscript analysis/spectral_pass/sp_summary.R --datasets hESC,hHep,K562,HepG2,RPE1 --out results/spectral_pass/summary
