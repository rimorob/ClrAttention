. jobs_pc/common.sh
# mem: 2
# cpus: 1
# after: 7610_scgpt_ext_hepg2_score,7612_scgpt_ext_beeline_score,7613_scgpt_ext_k562_score
# scGPT extension: statistics over the four datasets already on disk (hESC, hHep, K562, HepG2), n = 4: d_i, P1/S1/S2, descriptive tables; at n = 4 NO interval and NO decision is possible (smallest p 0.0625). Reads the metric tables only. Output results/scgpt_ext/stats_n4/.
set -e
nice -n 19 Rscript analysis/scgpt_ext/stats_ext.R --in results/scgpt_ext --datasets hESC,hHep,K562,HepG2 --out results/scgpt_ext/stats_n4
