# Manifests of the scGPT extension (frozen before scoring; see docs/plans/scgpt_extension_plan.md)

Each dataset's manifest (the SHA-256 of every file the scoring reads: control-cell subsamples, gene set, per-cell totals, scGPT extraction outputs) is recorded here and committed before that dataset is scored.

| Dataset | Manifest SHA-256 | Frozen | Notes |
|---|---|---|---|
| K562 | 787cbf87...080b (as frozen for the earlier K562 arm; recomputed byte-identical on 21 Nov 2026) | earlier arm | raw single-cell file kept |
| HepG2 | e9ef6a4bb97762fd0ec1c0bc985b16ab4aefa6effe52641e248aff40eef5f3c5 | 21 Nov 2026, 14:20 | 5 files verified, 0 mismatches; raw single-cell file kept (also the source of the HepG2 pseudobulk) |
| RPE1 | b6638a394b885827c0cf3dcec911362c5908f336f2b65e240fae5210be081d73 | 21 Nov 2026, 15:27 | 5 files verified, 0 mismatches; 11,485 non-targeting controls, five subsamples; raw single-cell file deleted after extraction (below) |

Raw single-cell downloads:

| File | SHA-256 | Size |
|---|---|---|
| rpe1_raw_singlecell_01.h5ad (figshare 20029387) | 9b05ef1f81526216fa008d677e9e0d03dce9a2f7a95499a4fb81e505e9d88ef1 | 8.70 GB |
| GSE264667_jurkat_raw_singlecell_01.h5ad (GEO GSE264667; 12,013 non-targeting controls) | ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3 | 9.37 GB |

Deleted after control-cell extraction (approved by the PI, 21 Nov 2026): data/perturbseq/rpe1_raw_singlecell_01.h5ad, 8.70 GB, sha256 as above; the extracted control cells (results/scfm_rpe1/expr and cell_totals.tsv, 25 MB) are kept and reproduce the scoring bit for bit (tested with --expr files). Re-download: figshare article 20029387.
