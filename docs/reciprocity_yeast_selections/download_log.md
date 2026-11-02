# Yeast reciprocity preregistration: protocol step 2 (download) log

- Preregistration fixed in 02cbcfd (2026-11-01). Download run as scheduler job 0980 (`Rscript analysis/yeast_prep.R download`), 2026-11-01 21:10:27 to 2026-11-01 21:30:06.
- Downloaded: GSE42526/GSE42527 final gene-expression matrices, protocols and series matrices; GSE45115 series matrix; GPL11232 annotation (GEO text API); SGD_features.tab; 755 GSE45115 wild-type sample tables and the 20 sample tables of the first 10 GSE42527 mutants (C-locale order of deleted gene), via the GEO text API. **No truth file was downloaded.**
- Integrity: R's 60 s timeout interrupted three first attempts (GSE42527 protocols, GSE42526 final matrix, GSE42526 series matrix); `fetch()` retried them. Final sizes equal the server-reported lengths (GSE42526 final 67,507,834 B; GSE42526 series matrix 125,602,415 B), every .gz passes `gzip -t`, and all 775 sample tables contain `!sample_table_end`.
- Not obtained: the authors' growth-corrected matrix (`deleteome_all_mutants_svd_transformed.txt`). The URL returned the site's HTML page, so the file was deleted and its hash removed. This affects one descriptive sensitivity only; the URL is unverified.
- Disk: data/yeast 920M (budget 3 GB); sample tables 476M; free space after: 321G.
- SHA-256 of every file: `download_SHA256SUMS.txt` (copy of data/yeast/SHA256SUMS, 784 files).
