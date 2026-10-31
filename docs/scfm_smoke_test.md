# Learned-attention comparator, step 1: smoke test (2026-10-31)

Scope: technical sanity of scPRINT medium-v1.5 and scGPT whole-human attention extraction on the BEELINE hESC case. **Nothing was scored against any truth set, and no ground-truth file was opened** (the preregistration, `docs/preregistration_scfm.md`, is not yet committed with the foundation-model readouts frozen; scoring before that would invalidate the comparator). No leakage audit was done (step 2). Scripts: `analysis/scfm_smoke_scprint.py`, `analysis/scfm_smoke_scgpt.py`, `analysis/scfm_smoke_compare.py`. Outputs (about 2.8 GB, not in the repo): `~/scfm_models/smoke_out/`.

## Result at a glance

| Item | Result |
|---|---|
| `lamin init` (local sqlite, outside the repo) | pass, after one environment fix (httpx) |
| `import scprint` | pass, after one environment fix (setuptools<81) |
| Checkpoint SHA-256 | pass: `96ff7c5d9e1ad5076e818bde3c38f498c95cd354cc42029d90691b5ba2964ff8` (re-computed 2026-10-31, equals the value recorded 22 Oct; HF jkobject/scPRINT commit 8e1ad1a1bc40c1fc2f4f1700fe27923eea9a2bed) |
| Pascal (GTX 1060, sm_61) on the torch 2.3.1 sm_60 binary | **pass**: fp32 matmul/softmax/LayerNorm correct (matmul error vs CPU 6e-4 on 2048x2048 N(0,1)), both models ran end to end on device 1 |
| scPRINT loads and runs | pass, but needs a **shim** (checkpoint and installed package disagree on the architecture; below) |
| Attention extraction (all layers, all heads, per cell, fp32; special tokens removed, rows renormalised) | pass; no NaN/Inf, rows sum to 1 |
| Determinism (two identical runs) | pass: bitwise identical, both models |
| Expression-permuted control at smoke scale | runs; the permutation moves the mean attention no more than changing which 64 cells are used, i.e. no data dependence detectable at this scale (details below) |
| scGPT whole-human smoke test | pass (checkpoint already on disk), same checks |

## Environment changes (the only two; `docs/scfm_env_freeze.txt` regenerated)

1. `lamin init` crashed in the lamin hub client (`TypeError: Client.__init__() got an unexpected keyword argument 'proxy'`): `gotrue 2.9.1` needs httpx >= 0.26 but the venv had httpx 0.24.1. Fix: `uv pip install httpx==0.27.2` (httpcore 0.17.3 -> 1.0.9).
2. `import scprint` failed with `No module named 'pkg_resources'` (from `ctxcore`): setuptools 84.2.1 no longer ships it. Fix: `uv pip install setuptools==80.10.2`.

The freeze diff vs the 22 Oct file is exactly: httpx, httpcore, setuptools as above, plus `gdown==6.4.2` and `pysocks==1.7.1` (installed on 22 Oct when the scGPT weights were fetched, after the first freeze).

## `lamin init`: what it wrote and fetched

```
mkdir -p ~/scfm_models/lamin_instance && cd ~/scfm_models/lamin_instance
~/venvs/scfm/bin/lamin init --storage ~/scfm_models/lamin_instance/storage --schema bionty --name scfm_smoke
```

Output: `connected lamindb: anonymous/scfm_smoke`. Written:

- `~/scfm_models/lamin_instance/storage/c98543da8fa75c7f90d49be92ccb2dea.lndb` (the sqlite database, 2 MB; no ontologies populated).
- Under `~/.lamin/` (user-level, unavoidable): `current_user.env` (anonymous user), `current_instance.env`, `instance--anonymous--scfm_smoke.env`, `auto_connect`, and `bionty/versions/{sources_local.yaml,.current_sources.yaml,.lamindb_current_sources.yaml,bionty_base.lock}`. The `.env` files make every later Python process in any directory auto-connect to this instance (that is why a bare `import scprint` works now).
- Nothing in the repo. Network: the init contacts the lamin hub (anonymous lookup); no data download.

Separately, to map gene symbols to Ensembl IDs, `bionty.Gene.public(organism="human")` was queried (Ensembl release 110, 75,719 terms; read remotely from bionty's asset store, nothing cached by bionty); the table was saved by the script to `~/scfm_models/lamin_instance/ensembl110_human_genes.tsv` (5.9 MB) so step 2 can reuse it offline.

Consequence: the bionty `Organism` registry is empty, so scPRINT's high-level helpers that need it (`scprint.utils.load_genes`, the scdataloader `Collator`, `GNInfer`) fail until `scdataloader.utils.populate_my_ontology()` is run (a large ontology fetch). The smoke scripts bypass these helpers (gene tokens are the integer positions in the checkpoint's own `genes` list), so this is not needed for the preregistered R0 readout. R1 (the package's own gene-network output through `GNInfer`) would need it.

## scPRINT: shims (exact)

`scPrint.load_from_checkpoint(...)` fails twice (jsonargparse rejects `label_counts.cell_type_ontology_term_id`, then the class rejects extra hyperparameters). More importantly the checkpoint does not match the installed scprint 1.0.13 architecture, so the smoke script builds the model directly (`load_scprint()` in `analysis/scfm_smoke_scprint.py`):

- Drop the checkpoint hyperparameters the installed class does not accept: `compress_class_dim` (None), `cell_specific_blocks` (False), `depth_atinput` (False), `label_counts` (duplicate of `classes`). Set `precpt_gene_emb=None` (the gene embeddings are in the state dict) and `transformer="normal"` (no flash-attn on Pascal).
- The checkpoint has **no `depth_encoder` weights** and `expr_decoder.fc.0.weight` is 256x257 (sequencing depth is concatenated at the decoder input), whereas scprint 1.0.13 injects a depth token at the encoder input and has a 256x256 decoder layer. So the checkpoint was written by an earlier architecture variant of the code than the PyPI releases that are installable (1.0.13 and 1.6.1, the two inspected, both use the depth token). Handling: `expr_decoder.fc.0.weight` is dropped (decoder is not used), `strict=False` load (missing: `depth_encoder.*`, that one weight; unexpected: `grad_reverse_discriminator_loss.*`, training-only), and depth is never fed (`full_depth=None`), so the input is the 9 class tokens (base cell embedding + 8 classes, matching the 9-row `class_encoder`) followed by the gene tokens. All encoder weights (gene encoder, expression encoder, class encoder, 8 transformer blocks, final norm) load unchanged.
- Fidelity check of the shimmed encoder (functional, not a truth): the loaded organism classification head, read at class token index 1+i (the prefix without a depth token), predicts NCBITaxon:9606 for 64/64 hESC cells; the assay head gives one class (EFO:0008931, I believe Smart-seq2; not looked up) for 64/64 cells. The cell-type and tissue heads gave diffuse predictions (cell type CL:0001064 in 62/64 cells; not interpreted: the ontology tables are not populated, the input is the 1,161-gene case set rather than a typical 2,200-gene input, and the scale differs, below). **Residual risk, to be settled before the preregistration is frozen:** that the encoder reproduces the authors' original forward pass exactly cannot be verified from the heads alone. Options for step 2: (a) ask/inspect the scPRINT repo at the commit that produced `medium-v1.5` (git tags are non-semantic; v1.5 is not a PyPI release) and install that; (b) compare our attention with the package's own `GNInfer` output on one case.

## Data and input scale

- Case: `results/beeline_baselines/hESC_1000/expression.tsv` (the exported case matrix; only its gene list and values are used): 1,323 genes x 758 cells. Source: `data/beeline/BEELINE-data/inputs/scRNA-Seq/hESC/ExpressionData.csv`, 17,735 x 758.
- **Transform check (preregistration asks for it to be logged):** max 19.35, 51% exact zeros; `2^x - 1` gives per-cell column sums with median 1.13e6 (range 0.83e6 to 5.7e6), i.e. the values are log2(TPM + 1) of a TPM matrix (natural-log back-transform gives 2.4e7, not plausible). So the back-transform is `x_lin = 2^x - 1` (TPM).
- scPRINT (`normalization="log"`) applies `log2(1 + x)` internally, so passing TPM reproduces BEELINE's log2(TPM+1) values exactly at the model input. **The pretraining input was raw UMI-scale counts, not TPM** (values here reach 5e5); this scale mismatch is a step-2 issue (raw counts from GEO are the plan's sensitivity arm).
- scGPT: TPM normalised to 1e4 per cell, log1p, then per-cell quantile binning into 51 bins over the case genes (zeros stay bin 0), `<cls>` first with value 0. Binning is rank-based, so only the order of non-zero values matters. scGPT's `_digitize` breaks bin-edge ties at random: **the script seeds numpy before binning** (without it two identical runs differed by 0.5% relative L1; with it they are bitwise identical). Note that the pretraining binned over the cell's full gene set before gene selection; here binning is over the 1.2k case genes.
- Cells: 64 (seeded, seed 1; 256 for one timing run); control: seed 2 (a different 64 cells), and permutation of each gene independently across the 64 cells (seed 1001, same cells). Token set = case genes including zeros.

## Gene coverage (hESC_1000, 1,323 genes)

| Model | Used | Dropped | Why |
|---|---|---|---|
| scPRINT | 1,161 | 162 | 102 symbol not in Ensembl 110 (33 ERCC spike-ins and outdated symbols such as PHF17, HIST1H1B, H1F0, MLL3, AES); 60 symbol with several Ensembl IDs (see below). 0 lost to the vocabulary after mapping |
| scGPT | 1,221 | 102 | symbols absent from its 60,697-symbol vocabulary (ERCC spike-ins, outdated symbols) |
| Intersection | 1,158 | | |

Mapping used here is deliberately crude (exact HGNC symbol match to Ensembl 110, then to the checkpoint's 44,756-ID vocabulary). **The 60 "ambiguous" are mostly real, important genes** (POU5F1, GATA4, SOX7, TCF7L1, ID3, ...): Ensembl 110 lists patch/haplotype copies with the same symbol (POU5F1 has seven IDs, all in the vocabulary). Step 2 must resolve these to the primary-assembly ID (via the Ensembl chromosome table or HGNC IDs) and map outdated symbols via HGNC history; dropping them here costs transcription factors the comparison cares about.

## Attention extraction (as in the preregistration, R0)

Per cell, per layer, per head, in fp32: `softmax(QK^T / sqrt(d))` from the captured Q and K, taken over gene keys only for gene queries. That equals "strip the special tokens (rows and columns) and renormalise each row to sum 1" exactly, and avoids 0/0: a first version that renormalised after the softmax produced NaN rows in scPRINT's last layer (146 rows of layer 8 head 2) because those queries put all mass on the special tokens (float32 underflow). Averaged over cells and stored as (layers x heads x G x G). Orientation (S = A^T) is applied at scoring, not here. For scGPT the first layer's attention was checked equal to `torch.nn.MultiheadAttention`'s own weights (max abs difference 0).

| | scPRINT medium-v1.5 | scGPT whole-human |
|---|---|---|
| Layers x heads, d_model | 8 x 4, 256 (18.7M parameters in the checkpoint, 17.8M of them in the encoder side that is used) | 12 x 8, 512 |
| Output shape (64 cells) | 8 x 4 x 1161 x 1161 | 12 x 8 x 1221 x 1221 |
| NaN / Inf; row sums | 0 / 0; 1.000000 | 0 / 0; 1.000000 |
| Model load | 18 s | 4 s |
| Inference, 64 cells | 3.3 s (256 cells: 12.4 s, about 0.05 s per cell) | 9.6 s |
| GPU memory, device 1 (GTX 1060) | peak allocated 1.36 GiB, reserved 1.7 GiB | peak allocated 2.5 GiB, reserved 3.7 GiB |
| Host RSS peak; wall incl. imports | 2.2 GB; 36 s | 3.0 GB; 43 s |
| Determinism, two identical runs | bitwise identical | bitwise identical (after seeding the binning) |

Device 0 (desktop GPU) was never used; 6 GB was enough for both models, so the 8 GB card is not needed even at 256 cells with batch size 8.

Sink check (important for step 3): the share of a gene query's attention that goes to the 9 special tokens is near 0 in some heads and near 1 in many (scPRINT, per layer x head, 64 cells): layer 1 head 4 0.9999; layer 2 heads 1-4 0.90, 0.97, 0.29, 0.80; layers 5 to 7 have 9 of 12 heads above 0.9; layer 8 head 3 1.000. In those heads the renormalised gene-gene rows are determined by tiny residual logits, so they are numerically fragile and likely noise-dominated. scGPT has almost no `<cls>` mass (< 0.02 except layer 1: heads 3-7 0.04, 0.19, 0.03, 0.07, 0.08). Relative row entropy of the stripped rows is 0.4 to 0.9 (scPRINT) and 0.78 to 0.98 (scGPT), i.e. no head is close to a point mass but scPRINT heads are far from uniform; the full table is in `~/scfm_models/smoke_out/metrics_run1.json` and `scgpt_metrics_run1.json`.

## Permutation control and variability at smoke scale (technical only)

`python -I analysis/scfm_smoke_compare.py ~/scfm_models/smoke_out` (output below). The yardstick is how far the mean attention moves when a different 64 cells are used.

```
== scPRINT medium-v1.5: A (8, 4, 1161, 1161)
different 64-cell subset (seed1 vs seed2): relL1 0.2448, Spearman 0.9039, top-1000 overlap 0.805
observed vs permuted across cells        : relL1 0.2151, Spearman 0.9289, top-1000 overlap 0.884
observed vs uniform 1/G: relL1 1.495
== scGPT whole-human: A (12, 8, 1221, 1221)
different 64-cell subset (seed1 vs seed2): relL1 0.0286, Spearman 0.9913, top-1000 overlap 0.883
observed vs permuted across cells        : relL1 0.0163, Spearman 0.9972, top-1000 overlap 0.928
observed vs uniform 1/G: relL1 0.693
```

(relL1 = sum of absolute differences of the all-layer all-head mean over its total.) Reading: in both models the attention is far from uniform and changes very little when the expression is shuffled across cells (permuted-vs-observed distance is smaller than or about equal to the distance between two cell subsets), so at 64 cells the learned networks are dominated by the gene-embedding prior, with at most a weak expression-dependent component that this test cannot separate from sampling variation. scPRINT moves more than scGPT in every layer (largest in the last layer, relL1 0.6). 64 vs 256 cells (scPRINT, same seed) differ by relL1 0.20, so the mean has not converged at 64 cells. This is a necessary, not a final, version of the prior-versus-data control: step 3 should use all cells of each dataset and several permutation seeds, and the preregistered random-weight and gene-embedding-only controls.

## How to reproduce (exact commands)

Scheduler: `tools/scheduler.sh` has no GPU resource (headers are `# mem`, `# cpus`, `# after`). No job was put in `jobs_pc/`. Each run was a single foreground job under a capped transient scope at nice 19, one at a time (never more than one GPU process, none on device 0):

```
source ~/venvs/scfm/bin/activate
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
cd ~/scfm_models/lamin_instance          # so the lamin instance auto-connects; scripts use absolute repo paths
R="systemd-run --user --scope -q -p MemoryMax=20G -p MemorySwapMax=0 -p CPUWeight=10 nice -n 19 python -I"
OUT=~/scfm_models/smoke_out; REPO=/home/boris/Projects/ClrAttention
for t in run1 run2; do $R $REPO/analysis/scfm_smoke_scprint.py --device 1 --ncells 64 --out $OUT --tag $t; done
$R $REPO/analysis/scfm_smoke_scprint.py --device 1 --ncells 64 --permute --out $OUT --tag perm
$R $REPO/analysis/scfm_smoke_scprint.py --device 1 --ncells 64 --seed 2 --out $OUT --tag seed2
$R $REPO/analysis/scfm_smoke_scprint.py --device 1 --ncells 256 --out $OUT --tag n256
# same four runs with analysis/scfm_smoke_scgpt.py (tags run1 run2 perm seed2; no n256)
python -I $REPO/analysis/scfm_smoke_compare.py $OUT
```

## What step 2 needs

1. **Input scale.** Raw counts from GEO for the 7 BEELINE datasets (scPRINT and scGPT were trained on counts; here TPM-scale values reach 5e5). Decide the sensitivity arm (TPM vs counts) before freezing; this smoke test cannot say which matters.
2. **Gene-ID mapping.** Primary-assembly Ensembl ID for symbols with patch copies (the 60 dropped here include POU5F1, GATA4, SOX7), HGNC alias/previous-symbol mapping for outdated symbols (BEELINE uses 2016-era symbols; 102 unmapped here), mouse tables for the 5 mouse datasets (scPRINT NCBITaxon:10090), the coverage and dropped-edge table per case, and G-intersection with scGPT. The saved Ensembl 110 table is `~/scfm_models/lamin_instance/ensembl110_human_genes.tsv`; a mouse table is not yet fetched.
3. **scPRINT code/checkpoint match.** Pin the scprint code version that matches `medium-v1.5` (or document the shim as the preregistered extraction and reconcile R1) and run `scdataloader.utils.populate_my_ontology()` only if R1 (`GNInfer`) is wanted.
4. **Cells.** Full cell sets (758 for hESC), several permutation seeds, and the random-weight controls (3 seeds, two versions). Cost at the measured rates: scPRINT about 0.05 s per cell per case, scGPT about 0.15 s per cell, so all 14 BEELINE cases plus controls are well under one GPU-hour on the 6 GB card.
5. **Sink heads.** Decide, before scoring and without looking at any truth, whether R0 averages over all heads and layers as preregistered (then heads with 90-100% special-token mass contribute near-noise rows); the preregistration says all heads, so this is only a reporting point (report special-token mass per head).
6. **Leakage audit** (GEO/SRA accessions of BEELINE and Replogle against cellxgene census and scBaseCount), then the preregistration change log entries: checkpoint hashes, the shim, the package versions, and the scGPT binning-seed fix.

## scGPT specifics

- Checkpoint `~/scfm_models/scgpt_human/` (`best_model.pt` 205 MB, `vocab.json`, `args.json`) was already on disk from 22 Oct (downloaded from the Google Drive source named in the scGPT repo; the source is still not pinned in the preregistration). SHA-256, verified 2026-10-31 and re-checked by the script at every run: best_model.pt `6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922`, vocab.json `acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f`, args.json `c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60`. No download was done in this step.
- Loading: `scgpt.model.TransformerModel` with the `args.json` architecture (12 layers, 8 heads, d 512, n_bins 51, `use_fast_transformer=False`), `scgpt.utils.load_pretrained(strict=False)` (renames `Wqkv` to `in_proj_`). All encoder, value-encoder and transformer weights load; `flag_encoder.weight` is unused (no batch labels) and the classifier head is absent from the checkpoint (not needed). No shim beyond that. The first layer's attention equals torch's MHA weights (checked).
- Not smoke-tested (small, deferred to step 3): the token-embedding cosine and hidden-state cosine readouts of the secondary arm. Both are trivial reads of `encoder.embedding.weight` and the layer outputs of the same forward pass.
