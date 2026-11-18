. jobs_pc/common.sh
# mem: 7
# cpus: 32
# after: 7330_tc_genie3_armA
# Arm A: GRNBoost2 by tools/grnboost2.py unchanged (legacy early stopping, seed 20261016, every gene a regulator), 32 workers (estimate 2.3 to 5.9 h).
set -e
~/venvs/clr-baselines/bin/python tools/grnboost2.py results/tcaclr_confirm/armA/expression.tsv results/tcaclr_confirm/armA/grnboost2_edges.tsv.part --workers 32 --seed 20261016
mv results/tcaclr_confirm/armA/grnboost2_edges.tsv.part results/tcaclr_confirm/armA/grnboost2_edges.tsv
