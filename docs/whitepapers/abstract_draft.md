# Abstract drafts: completed CLR attention (Paper 1)

*Drafts for the PI's rewrite, 2026-10-21. They use only claims that `docs/paper_plan.md` supports. Numbers, review statuses and exclusions are in `abstract_points.md`. [PENDING] marks the sentence to be written once the learned-attention comparator (`docs/plans/scfm_comparator.md`) is scored. Bracketed options inside that sentence are the possible outcomes.*

*How words are counted: each whitespace-separated token is one word, so hyphenated compounds and numbers count as one. The [PENDING] tag and Markdown italics are not counted.*

*Journal limits, checked on the journals' author pages on 2026-10-21:*
- *Nature Machine Intelligence Article: "Abstract – up to 150 words, unreferenced."*
- *Nature Communications: the abstract "should be no more than 200 words long and contain no references"; the maximum title length is 15 words.*

---

## (i) Nature Machine Intelligence abstract (150 words, unreferenced)

> Attention normalizes each query's compatibilities within its own context, an idea older than the transformer (kernel regression, fast weights). Its first widely adopted form in computational biology, context likelihood of relatedness (CLR), remains standard for de novo, unsupervised network inference but stops at the scores. We bring back from language modelling a temperature-controlled softmax, iterated attention and depth chosen on held-out data without labels. On strong-evidence *Escherichia coli* regulons, attention raises pooled recovery 16% over CLR, in all 1,000 resamplings, and, controlling for hub genes, leads GENIE3 and PIDC; the gain depends on CLR's normalization. Per regulon, it is on par with CLR and GENIE3. Iterated attention trades pair-level resolution for neighbourhood aggregation, helping large, loose regulons and hurting tight co-expression cliques such as ribosomal complexes. Gains are small in human K562 cells, and on single-cell data attention adds nothing to rank-based CLR. [PENDING] Zero-shot attention from single-cell foundation models [falls short of] both.

**Word count: 150** with the longest outcome option ("falls short of"); 148 with "matches" or "exceeds". 8 sentences, mean 18.8 words. No references.

**What each sentence rests on (plan section, review status):**

| Sentence | Evidence | Status |
|---|---|---|
| 1. Attention is older than the transformer | Nadaraya 1964; Watson 1964; Schmidhuber 1991; Schlag et al. 2021; Vaswani et al. 2017 (cited in the Introduction) | framing, plan §10 |
| 2. CLR, first widely adopted attention-like form; still standard for de novo, unsupervised inference | Faith, Hayete et al. 2007; DREAM5 (Marbach et al. 2012); strongest method tested on BEELINE | worked, §10 and §6 |
| 3. Softmax, iterated attention, label-free held-out depth | median t* = 9 over 50 lab splits; flat for t = 7–11 | worked, QC, §1 |
| 4a. Pooled recovery +16% over CLR, in all 1,000 resamplings | nested bootstrap +0.022 [+0.019, +0.026]; ahead in 1,000/1,000 | worked, QC, §1 |
| 4b. Leads GENIE3 and PIDC under hub control | within-hub-stratum AUPR leads on SC truths: M3D over GENIE3 and PIDC, EcoMAC over GENIE3 (reviewed); the EcoMAC PIDC cell is so far only in `sota_table.md` (1.57 against 1.40, not yet reviewed) | mixed, QC, §5 |
| 4c. Gain depends on CLR's normalization | uncalibrated propagation recovers about a third of the AUPR gain and none of the coherence gain | worked, QC, §2 |
| 5. Per-regulon parity | Holm p ≥ 0.21 against CLR, GENIE3 and PIDC | mixed, QC, §5 |
| 6. Resolution traded for aggregation | gain sits in regulons of ≥ 100 genes; 86% of the CORUM loss is ribosome and pre-rRNA complexes | QC, §5 |
| 7. Human gains small; single cell nothing | DoRothEA +0.006 [0.002, 0.010] under hub control; attention adds nothing over CLR-NS | mixed, QC, §5; did not work, QC, §6 |
| 8. [PENDING] learned attention | scPRINT (primary) and scGPT (secondary), zero-shot | not started, §5 |

**Notes for the rewrite**
- **"First widely adopted" is unhedged here** to save words. The Nature Communications version keeps "to our knowledge". Either way, the Introduction should separate CLR from topological overlap (see `abstract_points.md`, d).
- **GENIE3 and PIDC are named without a gloss** (tree-ensemble regression; partial information decomposition). An editor may ask for one.
- **The 16% is 0.022 / 0.136.** The current `paper.md` uses 15%, from one split. Use one basis everywhere.
- **"Held-out data"** means held-out laboratories for bulk data and random halves of cells for single-cell data. **"Strong-evidence"** means RegulonDB strong/confirmed.
- **"Both" in the last sentence** means rank-based CLR and estimated attention, which tie on single-cell data. If the comparator shows that CLR's normalization improves learned attention, add "and CLR's normalization improves it" and cut elsewhere: that is the strongest version of this paper for this journal.

---

## (ii) Nature Communications alternative (about 250 words)

The full text below is **260 words**, as asked for. Nature Communications caps abstracts at 200 words, so the spans in {braces} are marked for deletion. Deleting them gives the 200-word version that follows.

> {Inferring gene regulatory networks from expression data alone remains hard.} Context likelihood of relatedness (CLR), introduced in 2007, remains a standard for de novo, unsupervised network inference{, as distinct from prior-informed and regression-based methods}. CLR scores each gene's dependence on every other gene against that gene's own background. This per-query normalization is attention's central move, older than the transformer (kernel regression, fast weights); CLR was{, to our knowledge,} its first widely adopted form in computational biology. CLR, however, stopped at the scores. We bring back from language modelling what it lacked: a per-gene softmax temperature, iterated attention and depth chosen on held-out laboratories without labels. On strong-evidence *Escherichia coli* regulons, pooled recovery improves on CLR by 0.022 in precision–recall area ({95% interval 0.019–0.026; }about 16%), in all 1,000 resamplings of experiments. The lead over CLR, GENIE3 and PIDC survives control for hub genes, and without CLR's normalization the same propagation recovers only a third of the gain. Per regulon, the four methods are on par{, and CLR and GENIE3 remain better at detecting direct regulator–target edges}. Iterated attention trades pair-level resolution for neighbourhood aggregation, helping large, loose regulons and hurting tight co-expression cliques: in human K562 perturbation data it gains slightly on high-confidence transcription-factor targets and loses on protein complexes{, mostly ribosomal}. On single-cell data, a rank-based CLR improves on the original and attention adds nothing. [PENDING] Zero-shot attention from single-cell foundation models{, scored on the same cells and truths,} [falls short of both]. {In practice, single-step CLR should nominate direct edges and complexes, and attention should rank large regulatory programs.}

**Word count: 260** in full (13 sentences, mean 20 words).

**200-word version for submission** (with the {braced} spans deleted):

> Context likelihood of relatedness (CLR), introduced in 2007, remains a standard for de novo, unsupervised network inference. CLR scores each gene's dependence on every other gene against that gene's own background. This per-query normalization is attention's central move, older than the transformer (kernel regression, fast weights); CLR was its first widely adopted form in computational biology. CLR, however, stopped at the scores. We bring back from language modelling what it lacked: a per-gene softmax temperature, iterated attention and depth chosen on held-out laboratories without labels. On strong-evidence *Escherichia coli* regulons, pooled recovery improves on CLR by 0.022 in precision–recall area (about 16%), in all 1,000 resamplings of experiments. The lead over CLR, GENIE3 and PIDC survives control for hub genes, and without CLR's normalization the same propagation recovers only a third of the gain. Per regulon, the four methods are on par. Iterated attention trades pair-level resolution for neighbourhood aggregation, helping large, loose regulons and hurting tight co-expression cliques: in human K562 perturbation data it gains slightly on high-confidence transcription-factor targets and loses on protein complexes. On single-cell data, a rank-based CLR improves on the original and attention adds nothing. [PENDING] Zero-shot attention from single-cell foundation models [falls short of both].

**Word count: 200** with "falls short of both"; 198 with "matches both" or "exceeds both".

What each new sentence in the longer version rests on:
- **The direct-edge sentence:** TF-node AUPR ratio on M3D, CLR 6.75, GENIE3 6.54, attention 5.17 (mixed, QC).
- **"High-confidence transcription-factor targets":** DoRothEA A–B. CollecTRI shows no gain under hub control, so the abstract must not say "transcription-factor targets" unqualified.
- **The closing practical sentence:** the edge-detector versus module-ranker result. It is supported now and does not depend on the pending comparator.

---

## (iii) Three candidate titles

All three are at most 15 words (the Nature Communications limit) and use no abbreviations.

1. **Context matters: completing the attention in context likelihood of relatedness** (10 words). This echoes the construction of "Size matters:". It names the result an ML reader should remember: the per-query normalization CLR got right is what makes the rest of attention work. Risk: the wordplay, and "completing" asks the reader to accept the framing up front.
2. **Reach versus resolution: what attention adds to gene network inference, and what it costs** (14 words). This leads with the reframed mechanism and states the losses in the title, which protects against any reading of overclaim. It suits Nature Communications best.
3. **Attention before transformers: bringing language-model machinery back to gene regulatory inference** (11 words). This is the Nature Machine Intelligence hook and states the PI's goal. It claims only that CLR (2007) predates the transformer (2017), not that it is the first attention mechanism.

---

## (iv) Register note: the PI's paper with Jim Collins

**Citation.** Hayete, B., Gardner, T. S. & Collins, J. J. Size matters: network inference tackles the genome scale. *Mol. Syst. Biol.* **3**, 77 (2007). [doi:10.1038/msb4100118](https://doi.org/10.1038/msb4100118); PMID 17299414; PMCID PMC1828748.
- Retrieved from PubMed and the PMC full text.
- The title is "tackles the genome scale", not "comes of age" as recalled.
- It is a News & Views on Bonneau and colleagues' Inferelator study of the archaeon *Halobacterium salinarum* NRC-1 (cMonkey biclustering followed by LASSO regression).

**The register, measured.** Counts are on the PMC text. Italicized words lost in extraction (organism names, *in vivo*, *n* and *p*) were restored for counting, so the figures are approximate.
- **Length and rhythm.** About 1,250 words in 49 sentences; mean 26 words per sentence (median 26, range 8–49). Fifteen sentences run to 30 words or more and only three to 12 or fewer. Long, clause-rich sentences carry the argument; short ones serve as pivots ("This choice of an organism has two practical advantages.").
- **Framing.** It opens at the level of the field and its pressure: "The growing importance of microarray data challenges biologists, and especially the systems biology community, to come up with genome-scale analysis methods that can convert the large quantity of available high-throughput data into high-quality systems-level insights." Only after naming the problem class does it narrow to the method.
- **Historical positioning.** Methods are tools in "the bioinformatics toolbox", each judged by its "biologically plausible simplifying assumptions". The field is described as moving "from design to practical use, and from simple organisms to higher eukaryotes".
- **Numbers.** Concrete counts make the scale tangible and carry the argument:
  - "about 2400 question marks, each with thousands of possible answers";
  - 120 regulators against 268 chips;
  - "a nearly six-fold reduction in dimensionality".
- **Explanation by function.** Each method gets one plain sentence ("LASSO works by selecting a small set of the most likely regulators of a given gene, and simultaneously determines a quantitative influence function...").
- **Confidence.** Verdicts are stated flatly ("an impressive result", "a dramatic leap"), and hedges ("may") appear only in the outlook. The piece is evaluation-minded: "the obligatory next question: how much do we trust this network?" and "Biological yardsticks in general will gain in importance".
- **First-person plural.** Sparse: five uses, all in the last two paragraphs, all for judgement or outlook ("we believe", "we see the immediate future of network inference as promising and bright").
- **Occasional wit.** "Small genome or not", "more oomph from their microarrays", "the once and future model organisms".

**What the drafts match.**
- **Framing first.** Each draft names the problem class before the method: the definition of attention in (i), and the field-level opening sentence in (ii).
- **Methods placed in history and explained by function:**
  - kernel regression, fast weights, the transformer and CLR in 2007;
  - "scores each gene's dependence on every other gene against that gene's own background".
- **Numbers do the arguing.** 16%, 1,000 of 1,000, a third of the gain.
- **One short pivot sentence.** "CLR, however, stopped at the scores." echoes his short turns.
- **"We" only for what we did.** Verdicts are stated without hedging; the only hedges are "to our knowledge" and the pending sentence.
- **His evaluation habit appears as content:** hub control, per-regulon parity, and losses stated next to gains. Title 1 borrows his "Size matters:" construction.

**Where the drafts depart.**
- **Sentence length.** Sentences average 19–20 words, not 26: a 150- or 200-word abstract cannot spend 26-word sentences and still state five results.
- **No idiom or humour.** An abstract has no room for it. His one idiomatic move worth keeping in the main text is the "obligatory next question" turn, which would suit the Introduction's hand-off to the hub-control and per-regulon analyses.
