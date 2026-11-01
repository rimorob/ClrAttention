# EcoMAC / M3D overlap audit (label-free)

Script `analysis/ecomac_m3d_overlap.R` (28 s, capped scope, nice 19); outputs in `results/ecomac_m3d_overlap/` (`run.log` has every number below). No truth file was opened. Date 2026-11-01.

## Conclusion

EcoMAC is not an independent replication of M3D. About 38% of EcoMAC arrays (828 of 2,198) have a one-to-one counterpart in M3D, and those same pairs cover 91% of M3D's chips (828 of 907) and 92% of its experiments (431 of 466). EcoMAC is close to a superset of M3D plus 1,370 arrays that M3D does not contain, including every Genome 2.0 array. M3D carries no GEO accessions, so shared CEL files cannot be proved from identifiers; the claim rests on the similarity pattern below, and it cannot distinguish identical CEL files from independent replicates of identical conditions by the same labs.

## Data and gene mapping

- M3D: `data/E_coli_v4_Build_6`, 907 chips (`..._chips907probes4297.tab`) and 466 replicate-averaged experiments (`avg_..._exps466probes4297.tab`), 4,297 genes. Probe-set names `symbol_bNNNN_k` were mapped to b-numbers (second-to-last token), as in `analysis/baselines.R`.
- EcoMAC: `data/ecomac/EcoMAC.tsv`, 4,189 genes x 2,198 arrays; b-numbers from column 2 of `gene_names.tsv`; metadata from `samples.tsv` and `array_platform.csv`.
- No duplicated b-numbers in either. Shared genes: 4,056 (96.8% of EcoMAC, 94.4% of M3D). All similarity numbers use these 4,056.

## (a) Identifiers

M3D chip names, experiment names and descriptions contain no GSE/GSM/ArrayExpress accessions (0 matches), only names such as `dinI_U_N0025_r1`. EcoMAC has GSM ids for 1,178 arrays and GSE for most of the rest (`Not found` 214; literature citations or in-house for others). A direct accession join is therefore impossible. Weak metadata evidence: 29 author surnames are shared between EcoMAC `Author` and M3D `experimenter` (748 EcoMAC arrays, 630 M3D chips; `shared_author_surnames.csv`, for example Dwyer, Covert, Laubacher, Allen, Kohanski, Haddadin), which says only that the same labs contributed to both.

## (c) Normalization

- EcoMAC: described in `tools/get_ecomac.sh` (from the CovertLab redistribution) as 2,198 RMA-normalized arrays, two Affymetrix platforms (GPL199 E. coli Antisense 844 arrays, GPL3154 Genome 2.0 735 arrays; 619 unknown/not in GEO). Every array has nearly the same mean (9.105 to 9.107) and sd (1.503 to 1.511): a joint quantile-type normalization across both platforms. Sorted-value ranges across arrays are not exactly zero (max 1.73 at the tails).
- M3D: the local README says only "normalized"; the dump itself does not state RMA. Array means vary 7.13 to 8.62 and sd 1.00 to 1.75 (max sorted-value range 3.39), so it is not jointly quantile-normalized as released. I did not verify the published M3D method (Faith et al. 2008) here; whether it is RMA is not determined from the local files. Both are log2-like scales (EcoMAC 3.3 to 14.9, M3D 2.9 to 15.2).
- Consequence: the same physical array has different values in the two files, so raw cross-compendium Pearson r is capped below 1 (maximum 0.9934 across all 2,198 x 907 pairs), whereas duplicate arrays inside EcoMAC reach r = 1 (48 GPL199 arrays have a twin with r > 0.995, max 1, in another series).

## (b) Similarity

Rules fixed before computing (script header): Pearson on log2 after per-array standardisation, near-duplicate if r > 0.99, sensitivity 0.98 and 0.995; secondary rule on compendium-centred profiles (each gene centred on its own compendium mean) flagged if above the maximum seen in the time-matched null, which is the 518 EcoMAC arrays from series GSE20000 or later (released 2010+, after the September 2009 build; they cannot be in Build 6). Platform-matched: the null is also split GPL199 (137 arrays, 4 series) and GPL3154 (381).

Max Pearson from each EcoMAC array to any M3D chip (907):

| Group | n | median | max | > 0.98 | > 0.99 | > 0.995 |
|---|---|---|---|---|---|---|
| All EcoMAC | 2198 | 0.901 | 0.9934 | 815 (37.1%) | 406 (18.5%) | 0 |
| GPL199, any date | 844 | 0.989 | 0.9934 | 611 | 298 | 0 |
| GPL3154 | 735 | 0.749 | 0.888 | 0 | 0 | 0 |
| Null: GSE >= 20000, all | 518 | 0.828 | 0.963 | 0 | 0 | 0 |
| Null: GSE >= 20000, GPL199 | 137 | 0.916 | 0.963 | 0 | 0 | 0 |

Spearman > 0.99: 808 (36.8%), > 0.995: 337 (zero in the null). Against the 466 averaged experiments: > 0.98 703, > 0.99 157. Within-compendium context: M3D chip to nearest chip of another experiment median 0.976 (263 of 907 above 0.99, 70 above 0.995); EcoMAC to nearest array of another series median 0.926 (128 above 0.99).

Two things are therefore established. (1) The pre-stated primary threshold (0.99) flags 406 arrays (18.5%) and 0.98 flags 815 (37.1%); 0.995 flags none because of the normalization cap, so it is uninformative here, not evidence of absence. (2) Raw r alone is a weak discriminator, since distinct but similar conditions also reach 0.98 within a compendium; the null and the structure below carry the argument.

Centred Pearson (removes the shared mean-expression profile): the null maximum is 0.7036 (99.9th percentile 0.7027). Pre-build GPL199 arrays are bimodal against it: 659 of 707 above 0.75, versus none of the 137 null GPL199 arrays above 0.66; the GPL3154 arrays reach at most 0.73 (8 arrays above 0.70). Rule above null max: 904 EcoMAC arrays (41.1%), 0 in the null.

Post hoc additions (made after inspecting the distributions above): the mutual-nearest-neighbour requirement (the EcoMAC array is the best match of its own best M3D chip). 828 arrays are mutual nearest neighbours, all of them within the centred rule, and 0 in the post-build null. This is the headline call.

Structure that supports "same arrays": calls are whole series, one-to-one. Of 127 EcoMAC series/sources, 35 have a called array and 29 are called in >= 90% of their arrays, each array to a different M3D chip (GSE6836 266 of 284, GSE10158 43/43, GSE1121 43/43, MSB 3:97 62/62, ASAP-AffyData1 41/41, GSE11183 26/26, Cell 130(5) 17/17, J Bacteriol 186(10) 14/14 and others; `call_by_series.csv`). The matched pairs' median centred r is 0.825 versus 0.925 between M3D replicates inside M3D, which fits the same arrays under two different normalizations (a lower ceiling) and not noise.

Called arrays by group: GPL199 612 of 844 (86.6% of the 707 with GSE < 20000), non-GEO/unknown 216 of 619, GPL3154 0 of 735, GSE >= 20000 0 of 518. The remaining 95 pre-build GPL199 arrays that did not meet the rule have median raw r 0.976 (probably similar conditions, or near-ties that failed the mutual test); they are neither confirmed overlap nor confirmed independent. The 79 M3D chips with no EcoMAC mutual counterpart are not explained.

## What it implies

- Do not call EcoMAC an independent replication in the paper. Wording supported by the evidence: EcoMAC is a larger, partly overlapping compendium (about 38% of its arrays, and about 91% of M3D's chips, are the same or matched hybridizations), adds Genome 2.0 and post-2009 data, and agrees on direction as a robustness check with different batch handling (ComBat by platform).
- The paper's two "may share arrays" sentences (Results footnote and Limitations) can now be replaced by the numbers above; the Methods line "not yet an independent replication" is correct.
- Cleanest independent test available from the data in hand: re-score on EcoMAC arrays with no M3D counterpart (1,370 arrays: all 735 GPL3154, 518 post-2009, and unmatched others; note they form a different population from M3D's, so any change mixes platform and overlap). That is a new analysis, not done here.
- Caveat on this audit: identical CEL files cannot be proved without shared accessions; a lab re-measuring the same condition would give a similar signature. Either reading is a dependence for the purposes of replication. The series-number cutoff (GSE >= 20000 as post-build) is an approximation of release date; the 2009-09 build date is from the file timestamps.

Files: `overlap_rules.csv`, `similarity_summary.csv`, `per_array_best_match.csv` (all 2,198 arrays with best M3D chip and call), `call_by_series.csv`, `per_series_overlap.csv`, `centred_r_histogram.csv`, `scale_summary.csv`, `call_by_platform.csv`, `call_by_era.csv`, `m3d_within_nn_summary.csv`, `shared_author_surnames.csv`, `summary.txt`, `run.log`.
