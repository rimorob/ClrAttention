# Paper framing: is CLR "the earliest implemented form of attention"?

Literature check written 2026-09-23. It is meant to help choose a claim we
can defend, not to argue for one we have already chosen.

## 1. The strong claim does not survive as stated

"CLR is the earliest implemented form of attention" would be contested
immediately, for three independent reasons.

- **Kernel smoothing (Nadaraya 1964; Watson 1964).** Kernel regression is
  now routinely presented as the prototype of attention. A query is compared
  with keys through a kernel, and the normalized similarities weight an
  average of values. The standard deep-learning textbook (Dive into Deep
  Learning, "Attention Pooling by Similarity") introduces attention this
  way. Shalizi's notebook on attention argues that neural "attention"
  re-invented this classical estimator around 2015.
- **Fast-weight programmers (Schmidhuber 1991–92).** Schlag, Irie and
  Schmidhuber (ICML 2021) showed that linear transformers are formally
  equivalent to the 1991 fast-weight programmers. Schmidhuber publicly
  claims 1991 as the origin of "linearized self-attention".
- **Neural attention by name (Bahdanau, Cho & Bengio 2014).** Neural
  attention under that name dates from 2014, and the transformer from 2017.

CLR (Faith, Hayete et al. 2007) postdates the first two and predates the third. A
priority claim over all attention would be read as overreach and would
distract from the results.

## 2. What is defensible, and is actually the more interesting claim

CLR has the defining structure of attention: *scores*. For every gene (the
query), CLR computes a data-dependent compatibility with every other gene
(the keys), and normalizes it **per query** against that query's own
background. That row-wise z-score is the same role softmax plays in
attention: it turns raw compatibilities into weights that are comparable
within one query's context. Three differences from both kernel smoothing
and transformer attention make CLR a distinct instance, not a relabelling:

- **The compatibility is an estimated statistical dependence.** It is
  mutual information on continuous measurements, so it is sign-blind and
  handles non-monotone relationships. It is not a fixed kernel on a
  distance, and it is not a learned dot product of embeddings. No
  "meaning" vectors are needed, which matters when meanings are
  tissue-dependent.
- **The normalization is contextual and calibrated.** Row-wise background
  subtraction, and here a permutation null with FDR control, give the
  scores statistical semantics that learned attention weights lack.
- **CLR stopped at the scores.** It never applied the attention *readout*:
  there was no aggregation of values, and no depth. The operator was
  used only to rank edges.

A defensible framing is therefore:

> CLR (2007) computed attention scores over continuous molecular
> measurements, per-query context-normalized dependence weights, which
> predates the neural-attention era by seven years. What it lacked was the
> rest of the attention machinery. Adding the components that attention
> research made standard substantially improves regulon recovery on the
> benchmark CLR was originally validated on. Those components are softmax
> normalization with a temperature, the aggregation readout, and depth
> (stacked attention). The depth is chosen on held-out data.

This keeps the historical point, "CLR was an early form of attention",
honest: an early, independently arrived-at instance of the attention-score
pattern in computational biology. It does not claim priority over kernel
smoothing or fast weights.

## 3. What reviewers will raise, and how to answer

- **"Multi-hop attention mass is just network propagation or a
  random-walk diffusion kernel."** This is mathematically true. It is the
  lazy random walk of Zhou & Schölkopf (2004), and it is related to the
  diffusion kernels of Kondor & Lafferty (2002). It is also the propagation
  that Cowen et al. (2017) describe as a universal amplifier in
  bioinformatics. The answer should not deny this; the novelty has to live
  elsewhere.
  1. The contribution is the *kernel*. Using a calibrated, context-normalized
     MI attention matrix as the transition kernel beats the identical
     construction on |Pearson| at matched density (`pearson_top10` control).
  2. The contribution is also the attention design choices. The softmax
     temperature, set by an effective-neighbour count, is what made the
     operator work; hard-thresholded FDR operators leave most genes
     isolated.
  3. Depth is selected on held-out labs, not on the gold standard.
- **"Learned attention already does GRN inference."** Transformer-based
  GRN methods exist: GRNFormer (Bioinformatics 2025/26), prior-informed
  transformers (Advanced Science 2025), and attention maps from
  single-cell foundation models. A 2026 systematic evaluation (Kendiukhov,
  BMC Genomics) found that attention-derived edge scores from scGPT and
  Geneformer add no incremental value over trivial gene-level features for
  perturbation-target prediction. That makes a useful contrast. Our
  attention is *estimated, not learned*. It needs no training corpus and no
  gene embeddings, and it carries a calibrated null.
- **"The gains may be within noise."** This is what `analysis/followup.R`
  settles: a paired delete-half jackknife over genes, and a paired
  bootstrap over regulons. The claim should be made only for comparisons
  whose intervals exclude zero.
- **"Which parts of the transformer toolkit did not help?"** Report them,
  because honest negatives strengthen the paper:
  - Value transforms (a signed value, and a conditional-expectation W_V)
    did not help on this benchmark, because 99% of the selected edges are
    positively correlated.
  - Diffusing profiles and then correlating them gave little gain.
  - Stouffer and Euclidean combination are equivalent.

## 4. Suggested title directions

- "Context likelihood of relatedness as attention: calibrated
  mutual-information attention with depth improves regulon recovery"
- "CLR was an attention mechanism: completing it with softmax, readout
  and depth"

## Sources

- Faith, Hayete et al. 2007, *PLoS Biology* 5:e8 (CLR; Faith and Hayete contributed equally).
- Nadaraya 1964; Watson 1964 (kernel regression). Dive into Deep Learning,
  "Attention Pooling by Similarity":
  https://d2l.ai/chapter_attention-mechanisms-and-transformers/attention-pooling.html
- Shalizi, notebook on attention and transformers:
  https://bactra.org/notebooks/nn-attention-and-transformers.html
- Schlag, Irie & Schmidhuber 2021, "Linear Transformers Are Secretly Fast
  Weight Programmers", ICML: https://arxiv.org/abs/2102.11174
- Bahdanau, Cho & Bengio 2014 (neural attention); Vaswani et al. 2017
  (transformer).
- Veličković et al. 2018, Graph Attention Networks, ICLR:
  https://arxiv.org/abs/1710.10903
- Kondor & Lafferty 2002, Diffusion kernels on graphs, ICML.
- Zhou & Schölkopf 2004 (lazy random-walk regularization; already cited as
  `zhouscholkopf2004`).
- Cowen et al. 2017, "Network propagation: a universal amplifier of genetic
  associations", *Nat Rev Genet*: https://www.nature.com/articles/nrg.2017.38
- Kendiukhov 2026, "Systematic evaluation of single-cell foundation model
  interpretability …", *BMC Genomics*:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC13390336/
- GRNFormer, *Bioinformatics*:
  https://academic.oup.com/bioinformatics/article/42/4/btag144/8540455
- Weng et al. 2025, prior-knowledge transformer for GRN inference,
  *Advanced Science*:
  https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202409990
