. jobs_pc/common.sh
# mem: 2
# cpus: 1
# after: 7550_sota_ext_m3d,7550_sota_ext_ecomac,7550_sota_ext_k562
# Extended SOTA comparison, assemble: bold-tie marks over all eight methods (aCLR, tcaCLR, tcaCLR2, CLR, GENIE3, PIDC, GRNBoost2, |r|), over the original six, and in the legacy convention (with Hub only in raw rows); paired differences with 95% gene-jackknife intervals (tcaCLR-CLR, tcaCLR-aCLR, tcaCLR2-CLR, tcaCLR2-aCLR, tcaCLR2-tcaCLR) for the 18 rows raw and hub-controlled; reproduction gate against docs/whitepapers/sota_table.md (printed values and marks of the old methods; stops on any difference). Output results/sota_extended/table_extended.md/.csv, marks_by_row.csv, paired_differences.csv, gate_assemble.csv. tcaCLR2 is exploratory; intervals resample genes only. Seconds.
set -e
nice -n 19 Rscript analysis/sota_extended.R --matrix assemble
