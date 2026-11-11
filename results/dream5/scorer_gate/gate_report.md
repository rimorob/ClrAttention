# DREAM5 scorer gate: PASS

Run: 2026-11-10 16:43:14. Pass rule: |ours - published| <= 0.001 for the six per-network AUPR/AUROC values and the three aggregate scores of a team (R, Python port, and matrix path where it ran); implementations vs each other 1e-9.
Inputs: pred_dir `data/dream5/published_extracted/flat_pred`, results `data/dream5/published_extracted/published_table_tidy.csv`. Matrix path: run; 0 of 108 team-networks skipped (file repeats a pair).
Teams compared: 36 (Bayesian1, Bayesian2, Bayesian3, Bayesian4, Bayesian5, Bayesian6, Community, Correlation1, Correlation2, Correlation3, Meta1, Meta2, Meta3, Meta4, Meta5, MI1, MI2, MI3, MI4, MI5, Other1, Other2, Other3, Other4, Other5, Other6, Other7, Other8, Regression1, Regression2, Regression3, Regression4, Regression5, Regression6, Regression7, Regression8). Teams fully matching the published table: 36; required 3.
Max |R - port| = 8.82e-12; max |R - matrix path| = 3.89e-15 (AUPR, AUROC, P_AUPR, P_AUROC and their -log10).
Teams in the results table without a prediction directory: none.

Columns: max abs difference to the published table over the six per-network AUPR/AUROC values and over the three aggregate scores, for R, the Python port and the matrix path; the last column lists failing quantities (impl:quantity[net]).

| team | matches | R: AUPR/AUROC | R: scores | port: AUPR/AUROC | port: scores | matrix: AUPR/AUROC | matrix: scores | failing |
|---|---|---|---|---|---|---|---|---|
| Bayesian1 | yes | 0.000295 | 0.000419 | 0.000295 | 0.000419 | 0.000295 | 0.000419 |  |
| Bayesian2 | yes | 0.000331 | 0.00037 | 0.000331 | 0.00037 | 0.000331 | 0.00037 |  |
| Bayesian3 | yes | 0.000461 | 0.000368 | 0.000461 | 0.000368 | 0.000461 | 0.000368 |  |
| Bayesian4 | yes | 0.000354 | 0.000401 | 0.000354 | 0.000401 | 0.000354 | 0.000401 |  |
| Bayesian5 | yes | 0.000338 | 0.000407 | 0.000338 | 0.000407 | 0.000338 | 0.000407 |  |
| Bayesian6 | yes | 0.000432 | 0.000428 | 0.000432 | 0.000428 | 0.000432 | 0.000428 |  |
| Community | yes | 0.000457 | 2.69e-05 | 0.000457 | 2.69e-05 | 0.000457 | 2.69e-05 |  |
| Correlation1 | yes | 0.000408 | 0.000487 | 0.000408 | 0.000487 | 0.000408 | 0.000487 |  |
| Correlation2 | yes | 0.00032 | 0.000427 | 0.00032 | 0.000427 | 0.00032 | 0.000427 |  |
| Correlation3 | yes | 0.000579 | 0.000363 | 0.000579 | 0.000363 | 0.000579 | 0.000363 |  |
| Meta1 | yes | 0.000476 | 0.000301 | 0.000476 | 0.000301 | 0.000476 | 0.000301 |  |
| Meta2 | yes | 0.000483 | 0.000392 | 0.000483 | 0.000392 | 0.000483 | 0.000392 |  |
| Meta3 | yes | 0.000376 | 0.000487 | 0.000376 | 0.000487 | 0.000376 | 0.000487 |  |
| Meta4 | yes | 0.000466 | 0.000415 | 0.000466 | 0.000415 | 0.000466 | 0.000415 |  |
| Meta5 | yes | 0.000467 | 0.000479 | 0.000467 | 0.000479 | 0.000467 | 0.000479 |  |
| MI1 | yes | 0.000513 | 0.000367 | 0.000513 | 0.000367 | 0.000513 | 0.000367 |  |
| MI2 | yes | 0.000496 | 0.000493 | 0.000496 | 0.000493 | 0.000496 | 0.000493 |  |
| MI3 | yes | 0.00049 | 3.89e-05 | 0.00049 | 3.89e-05 | 0.00049 | 3.89e-05 |  |
| MI4 | yes | 0.000339 | 0.00039 | 0.000339 | 0.00039 | 0.000339 | 0.00039 |  |
| MI5 | yes | 0.000439 | 0.000404 | 0.000439 | 0.000404 | 0.000439 | 0.000404 |  |
| Other1 | yes | 0.000359 | 0.000353 | 0.000359 | 0.000353 | 0.000359 | 0.000353 |  |
| Other2 | yes | 0.000461 | 0.000465 | 0.000461 | 0.000465 | 0.000461 | 0.000465 |  |
| Other3 | yes | 0.000469 | 0.000422 | 0.000469 | 0.000422 | 0.000469 | 0.000422 |  |
| Other4 | yes | 0.000493 | 0.000288 | 0.000493 | 0.000288 | 0.000493 | 0.000288 |  |
| Other5 | yes | 0.000431 | 0.000478 | 0.000431 | 0.000478 | 0.000431 | 0.000478 |  |
| Other6 | yes | 0.000411 | 0.000202 | 0.000411 | 0.000202 | 0.000411 | 0.000202 |  |
| Other7 | yes | 0.00039 | 0.000165 | 0.00039 | 0.000165 | 0.00039 | 0.000165 |  |
| Other8 | yes | 0.000454 | 0.000165 | 0.000454 | 0.000165 | 0.000454 | 0.000165 |  |
| Regression1 | yes | 0.000411 | 0.000467 | 0.000411 | 0.000467 | 0.000411 | 0.000467 |  |
| Regression2 | yes | 0.00047 | 0.000487 | 0.00047 | 0.000487 | 0.00047 | 0.000487 |  |
| Regression3 | yes | 0.000409 | 0.00048 | 0.000409 | 0.00048 | 0.000409 | 0.00048 |  |
| Regression4 | yes | 0.000399 | 0.000326 | 0.000399 | 0.000326 | 0.000399 | 0.000326 |  |
| Regression5 | yes | 0.000394 | 0.000197 | 0.000394 | 0.000197 | 0.000394 | 0.000197 |  |
| Regression6 | yes | 0.000402 | 0.000491 | 0.000402 | 0.000491 | 0.000402 | 0.000491 |  |
| Regression7 | yes | 0.000411 | 0.000474 | 0.000411 | 0.000474 | 0.000411 | 0.000474 |  |
| Regression8 | yes | 0.000277 | 0.000437 | 0.000277 | 0.000437 | 0.000277 | 0.000437 |  |
