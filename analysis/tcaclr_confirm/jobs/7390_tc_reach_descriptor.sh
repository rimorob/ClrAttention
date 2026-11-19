. jobs_pc/common.sh
# mem: 16
# cpus: 8
# Registered reachability descriptor for both arms (structure only: expression and complex membership; no score, no truth comparison). The scorer does not emit it.
set -e
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript analysis/tcaclr_confirm/reach_descriptor.R --expr results/tcaclr_confirm/armA/expression.tsv --members docs/tcaclr_confirm_selections/armA_complex_portal_members.tsv --out results/tcaclr_confirm/result_A/reach_descriptor --threads 8
Rscript analysis/tcaclr_confirm/reach_descriptor.R --expr results/tcaclr_confirm/armB/expression.tsv --members docs/tcaclr_confirm_selections/armB_cyc2008_members.tsv --out results/tcaclr_confirm/result_B/reach_descriptor --threads 8
