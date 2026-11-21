. jobs_pc/common.sh
# mem: 4
# cpus: 2
# after: 7601_scgpt_ext_hepg2_extract
# scGPT extension, HepG2: freeze (SHA-256 manifest of the extraction in results/scfm_hepg2/sha256_manifest.csv; COMMIT it and record its hash in the change log BEFORE any scoring job runs) and export of the scGPT matrices to flat float32 files for the scorer (verifies the manifest first; exit 2 on a mismatch). No truth file. Output ~/scfm_models/step4_hepg2.
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
$HOME/venvs/scfm/bin/python -I analysis/scfm_k562_manifest.py --dir $HOME/scfm_models/step3/hepg2 --out results/scfm_hepg2/sha256_manifest.csv | tee results/scfm_hepg2/manifest_hash.txt
SHA=$(sha256sum results/scfm_hepg2/sha256_manifest.csv | cut -d' ' -f1)
$HOME/venvs/scfm/bin/python -I analysis/scfm_score_export.py --k562 --tag hepg2 --manifest-sha $SHA --out $HOME/scfm_models/step4_hepg2
