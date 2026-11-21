. jobs_pc/common.sh
# mem: 2
# cpus: 1
# after: 7610_scgpt_ext_hepg2_score,7612_scgpt_ext_beeline_score,7613_scgpt_ext_k562_score,7710_scgpt_ext_rpe1_score,7810_scgpt_ext_jurkat_score
# scGPT extension: the pre-specified statistics over six datasets (hESC, hHep, K562, HepG2, RPE1, Jurkat): P1, S1, S2 with Holm, descriptive tables, leave-one-dataset-out, fresh-only blocks. Reads the metric tables only. Output results/scgpt_ext/stats_n6/.
set -e
nice -n 19 Rscript analysis/scgpt_ext/stats_ext.R --in results/scgpt_ext --datasets hESC,hHep,K562,HepG2,RPE1,Jurkat --out results/scgpt_ext/stats_n6
