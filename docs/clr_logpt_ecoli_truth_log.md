# E. coli complex truth build log (CLR(log P^t) preregistration, protocol step 4)

Built 2026-10-26 by `analysis/build_ecoli_complex_truth.R`. Build only: no method was scored, no AUPR computed, no scores read.
Logged here rather than in the preregistration's own change log because editing that file changes its last-commit hash and would invalidate the K562 reproduction check (`gate_repro`); the PI may copy this entry into the change log (and rerun the repro check) if wanted.

## Source chosen: EBI Complex Portal taxon 83333 (the preregistered availability fallback)
Why: EcoCyc `proteins.dat` is not obtainable without a BioCyc login/licence (https://ecocyc.org/ftp/ redirects to https://biocyc.org/ftp/?orgid=ECOLI; the SRI flat-file directory answers 401). Chosen on availability alone, before any scoring.
- Complexes: https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/83333.tsv, downloaded 2026-10-26, server Last-Modified 2026-01-14, 324 complexes, SHA-256 `46dd5690d8016ebb98d77a695de4f1842a6113ef401d6956791e128bdfffc8b1`.
- UniProt to b-number table (accession, ordered-locus names; query `organism_id:83333`), https://rest.uniprot.org/uniprotkb/stream, UniProt release 2026_03 (2 Sep 2026), downloaded 2026-10-26, SHA-256 `31dd3f9a09b2cd475a58c1d0d2b1411b113329eed91291fe7bb90052e8d019b6`.
- Raw copies: `data/ecoli_complexes/` (with `SHA256SUMS`).

## Rules applied (as preregistered)
Subcomplexes (CPX- participants) flattened recursively; non-protein participants (CHEBI, RNA URS, others) dropped; isoform/chain suffixes mapped to the parent accession; accession to b-number via UniProt ordered-locus names (1 of 993 flattened proteins unmapped); per dataset: in-universe genes (rows of M3D 466-experiment matrix and EcoMAC with non-constant rows, 4,297 and 4,189 genes), size 3-200, drop exact duplicates and proper subsets, complete-linkage hclust at Jaccard >= 0.5 (cutree h = 0.5), union per cluster. Pairs by `comembership_pairs()` with the RegulonDB operon map.
Ribosomal genes: `^rpl[A-Y]$|^rpm[A-J]$|^rps[A-V]$` symbols in `ecoli_k12_genes.tsv` (54 genes, all in both universes).

## Evaluability (truth structure only)
| dataset | merged complexes | genes | pairs (operon excl.) | positives | complexes with positives | operon-included positives | fallback needed |
|---|---|---|---|---|---|---|---|
| M3D | 129 | 561 | 156,267 | 1,150 | 54 | 1,828 | no (>= 200 and >= 10) |
| EcoMAC | 124 | 544 | 146,902 | 1,148 | 53 | 1,807 | no |

The operon-excluded main truth is therefore the scored endpoint in both compendia. Full counts for every variant (operon-included, size >= 2, no merge, ribosome-excluded >= 30% r-protein genes, r-protein module): `docs/clr_logpt_selections/ecoli_complexes_counts.csv`.

## Files
- `docs/clr_logpt_selections/ecoli_complexes_members.tsv` (long format: variant, dataset, complex, bnumber; variants main, size2, nomerge, noribo, rp_module)
- `docs/clr_logpt_selections/ecoli_complexes.sha256`: `72663d6986d0db8d3fd406be2ac77fa1656595086471400f5d72fc4426ce1c90`
