. jobs_pc/common.sh
# mem: 4
# cpus: 2
# scGPT extension, HepG2: gene set (CollecTRI read ONLY for regulon sizes), 5 seeded subsamples of 1,024 non-targeting controls (fewer than 1,024 controls: all cells, one subsample; fewer than 512: aborts), extractor inputs, gene maps, per-cell UMI totals and gene-set counts (results/scfm_hepg2/cell_totals.tsv, expr/), so that the raw h5ad can be deleted afterwards. Reads control cells only; no scoring. Output results/scfm_hepg2/, ~/scfm_models/step2/inputs/hHepG2_s*.npz, results/scfm_step2/gene_map/hHepG2_s*.tsv. Measured (full HepG2 input, scratch outputs): 11 s, peak RSS 0.76 GB; 4,976 controls, 1,253 genes (298 TFs + 1,000 variable, 45 overlap).
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
$HOME/venvs/scfm/bin/python -I analysis/scfm_k562_prep.py --h5 data/hepg2_perturbseq/GSE264667_hepg2_raw_singlecell_01.h5ad --tag hHepG2 --out results/scfm_hepg2 --control-label non-targeting
