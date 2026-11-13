# Complex estimation on DREAM5 net4: main table (written by the summariser; no interpretation)

| Row | Contrast | Point estimate | Total 95% interval (gene + series) | Total interval (gene + group) | Gene-only interval (E. coli-comparable) | Gene replicates: 2.5% / 97.5% percentiles (centred on d), skewness |
|---|---|---|---|---|---|---|
| P | a3r - attention | +0.0962 | [+0.0675, +0.1250] | [+0.0714, +0.1211] | [+0.0748, +0.1177] | [+0.0775, +0.1190]; 0.415 |
| C1 | a3r - CLR | +0.0385 | [+0.0099, +0.0671] | [+0.0229, +0.0541] | [+0.0276, +0.0495] | [+0.0281, +0.0514]; 0.387 |
| C2 | attention - CLR | -0.0577 | [-0.0979, -0.0175] | [-0.0815, -0.0340] | [-0.0757, -0.0397] | [-0.0749, -0.0396]; -0.105 |

> a3r − CLR has a smaller SE and may exclude 0 where a3r − attention does not; this is never a basis for 'a3r beats CLR' (E. coli S3).

SE (main truth, raw): gene 0.0110; series 0.0098; group 0.0064 (P row)

## Fixed sentences (row P; basis: total (gene + series) 95% interval)

- The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test).
- The whole interval lies above the smallest effect of interest (+0.005).
- The whole interval lies above the E. coli estimate (+0.027).
- Exploratory estimation, not a confirmatory test, a replication or a validation. a3r is not adopted; this study cannot change that by itself (PI decisions 12 and 13 Nov).

Side-by-side lines (never the reported sentences):

- [gene jackknife only, comparable to the E. coli test] The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test).
- [gene jackknife only, comparable to the E. coli test] The whole interval lies above the smallest effect of interest (+0.005).
- [gene jackknife only, comparable to the E. coli test] The whole interval lies above the E. coli estimate (+0.027).
- [group-level experiment unit] The total interval excludes no difference: a3r exceeded attention on yeast complexes (preregistered exploratory estimation; not a confirmatory test).
- [group-level experiment unit] The whole interval lies above the smallest effect of interest (+0.005).
- [group-level experiment unit] The whole interval lies above the E. coli estimate (+0.027).

See estimate_final.csv for every secondary and sensitivity line (all gene-jackknife intervals; no p-value anywhere).
