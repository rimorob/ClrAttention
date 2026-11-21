. jobs_pc/common.sh
# mem: 4
# cpus: 2
# scGPT extension, Jurkat: gene set (CollecTRI read ONLY for regulon sizes), 5 seeded subsamples of 1,024 non-targeting controls (fewer than 1,024 controls: all cells, one subsample; fewer than 512: aborts), extractor inputs, gene maps, per-cell UMI totals and gene-set counts (results/scfm_jurkat/cell_totals.tsv, expr/), so that the raw h5ad can be deleted afterwards. Reads control cells only; no scoring. Output results/scfm_jurkat/, ~/scfm_models/step2/inputs/hJurkat_s*.npz, results/scfm_step2/gene_map/hJurkat_s*.tsv. NOT TO BE SUBMITTED before the download (README) and the control-label check; the raw h5ad is deleted only after this job, the scoring inputs exist (results/scfm_jurkat/expr, cell_totals.tsv) and the deletion is logged in docs/disk_budget.md.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
$HOME/venvs/scfm/bin/python -I analysis/scfm_k562_prep.py --h5 data/jurkat_perturbseq/GSE264667_jurkat_raw_singlecell_01.h5ad --tag hJurkat --out results/scfm_jurkat --control-label non-targeting
