. jobs_pc/common.sh
# mem: 5
# cpus: 32
# after: 7331_tc_genie3_armB
# Arm B: GRNBoost2 by tools/grnboost2.py unchanged (legacy early stopping, seed 20261016, every gene a regulator), 32 workers (estimate 2.3 to 5.9 h).
set -e
~/venvs/clr-baselines/bin/python tools/grnboost2.py results/tcaclr_confirm/armB/expression.tsv results/tcaclr_confirm/armB/grnboost2_edges.tsv.part --workers 32 --seed 20261016
mv results/tcaclr_confirm/armB/grnboost2_edges.tsv.part results/tcaclr_confirm/armB/grnboost2_edges.tsv
