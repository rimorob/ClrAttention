# Contributions and AI-use disclosure (draft for the PI to correct)

*Drafted 2026-10-21 by Claude (Anthropic) from the project history, for the PI to edit. Items marked [confirm] are my reading of the history and need the PI's check. Before submission, check the current wording of the journal's policy (Nature Portfolio editorial policies on AI, authorship and images) and adapt the final section.*

## How Nature Portfolio journals treat AI use (summary to verify before submission)

- Large language models do not meet authorship criteria; accountability rests with the human authors.
- Use of an LLM or other AI tool beyond copy editing should be documented in the Methods (or an equivalent section).
- Generative-AI images are generally not accepted in published content (narrow exceptions). Our figures are made with ggplot from computed results; any AI-generated schematic or graphical-abstract draft is a sketch for redrawing, not a submitted figure.
- Code and data availability statements apply as usual; all code is in the repository with its history.

## Contributions by the PI (Boris) [confirm and extend]

**Scientific idea and framing**
- The central thesis: CLR is a context-normalized attention operator, the first attention mechanism widely adopted in computational biology. The paper's aim is to bring advances in NLP attention back to CLR.
- The framing of CLR's robustness: it is the standard for de novo, unsupervised inference, distinct from prior-informed and regression or classification methods such as GENIE3.
- Historical scoping: CLR is not the first attention mechanism; the paper cites kernel regression and fast weights.
- Prior work on CLR, including the original implementations and the tiling approach for genome-scale matrices [confirm the authorship statement for the 2007 work].

**Experimental design decisions that shaped the results**
- Use a larger, multi-platform E. coli compendium (EcoMAC), and inspect platform batch effects ("Christmas tree" gene-pair plots). Use platform-level ComBat batch correction.
- Human data must explore many network states, so a single condition is the wrong test. This led to the Perturb-seq perturbation compendia.
- Normal-score (Blom) transformation for zero-inflated single-cell data, which produced CLR-NS, our best single-cell variant.
- Gene filters based on multimodality and conditional multimodality (stepAIC on Horn-retained principal components), not variance. Conditional relationships require batch removal first.
- Top principal components as condition-specific sigma-factor programs, and the sigmulon truth definition (drop σ70 and σS).
- Regulation as a union of networks ("leaves of the Venn diagram"). This led to the depth sweep as a multiscale, fine-to-coarse analysis, and the insistence that depth be analysed on a log scale.
- Per-regulon (macro) mean and median AUPR instead of filtering large regulons, with no merging of overlapping regulons, because multi-membership is meaningful.
- Questioning resampling: replicate averaging and bootstrap duplicates.
- Requiring empirical tests of all attention-width options; requiring a learned-attention (foundation model) comparator; differential CLR (diffCLR) for perturbation targets, possibly as a second paper.
- Figure standards: a Shapley decomposition of components and a consistent ggplot visual language.

**Process and quality control**
- The standing rule that every design, method choice and interpretation is checked by an independent high-effort AI reviewer before being presented. Requiring scheduling rules (memory safety, ~70-80% utilization) and status reporting.
- Choice of venue and authorship (with Jim Collins); all final decisions on claims and wording.

## Contributions by AI tools

**Claude (Anthropic), Opus model, in Claude Code, working in the PI's repository and workstation under the PI's direction**
- **Software:** ported and extended the CLR code base (R/C++), including the packed, band-parallel MI→CLR path adapted from the PI's tiling C code. Wrote the analysis scripts (attention operators, depth selection, bootstraps and jackknifes, baselines, macro and hub-controlled evaluation, BEELINE and Perturb-seq pipelines), the job scheduler and the status tools.
- **Computation:** ran all analyses through the job queue and monitored and fixed failures. Every run is logged in `jobs_pc/runner.log` and the repository history.
- **Literature:** retrieved and summarized literature (for example on foundation models and benchmarks). All citations are to be verified by the authors.
- **Writing:** drafted the paper text, supplement, plans, preregistrations, abstracts and this file, for the PI to rewrite.

**Independent AI review (Claude subagents, high effort, separate contexts)**
- Statistical and methodological reviews of designs and results before they were reported to the PI. They caught overclaims, a resampling bug that doubled regulon counts, hub-bias confounding (and an over-correction of it), and scheduler bugs; they also specified preregistered designs (keff variants, macro metrics, comparator choice). Reviews are recorded in `docs/paper_plan.md` (entries marked "QC").

**Image generation (if used)**
- PaperBanana (Gemini) may be used only for draft schematics or graphical-abstract sketches for an illustrator. None is to be submitted as a figure.

## Suggested Methods statement (draft)

> Software development, execution of analyses, literature retrieval and drafting of text were assisted by Claude (Anthropic; Opus model) used through Claude Code, under the direction of the authors. Statistical designs and results were additionally checked by independent instances of the same model acting as reviewers before being accepted. All analysis code, logs and preregistrations are available in the repository. The authors reviewed, verified and take full responsibility for all code, results and text. No AI-generated images are included in the figures.
