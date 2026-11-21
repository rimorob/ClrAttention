. jobs_pc/common.sh
# mem: 8
# cpus: 2
# after: 7700_scgpt_ext_rpe1_prep
# scGPT extension, RPE1: scGPT whole-human extraction on GPU DEVICE 1 (GTX 1060; the extractor refuses device 0), 5 subsamples, settings unchanged (last layer, 51 bins, binning seed 20261021), then the truth-free QC gate. No truth file. Same extraction as HepG2 (about 14 min).
set -e
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=1
REPO=/home/boris/Projects/ClrAttention; O=$HOME/scfm_models/step3/rpe1; mkdir -p $O
nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | awk -F', ' '$1==1 && $2+0 > 1500 {exit 1}'   # device 1 must be free
du -sh $HOME/scfm_models | tee /dev/stderr | awk "{ if (\$1+0 > 30 && \$1 ~ /G/) exit 1 }"
cd $HOME/scfm_models/lamin_instance
$HOME/venvs/scfm/bin/python -I $REPO/analysis/scfm_extract_scgpt.py --cases hRPE1_s1,hRPE1_s2,hRPE1_s3,hRPE1_s4,hRPE1_s5 --out $O/scgpt_A --device 1
$HOME/venvs/scfm_alt/bin/python -I $REPO/analysis/scfm_qc_dir.py $O/scgpt_A 5 scgpt
