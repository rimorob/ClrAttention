# Graphical abstract: prompt and brief

*Drafted 2026-10-21. Nature Portfolio journals do not accept generative-AI images in published content (narrow exceptions), so the published graphical abstract must be drawn by a person or built in standard tools. Use this prompt for a PaperBanana/NanoBanana draft, and the same text as the brief for an illustrator. Update panel 4 once the foundation-model comparison is in.*

## Message in one sentence

CLR, a 2007 network-inference method, is a context-normalized attention mechanism. Completing it with a softmax and multi-step depth, without any learning, recovers large regulatory programs better than CLR and established inference methods. Context normalization is the key ingredient, and depth trades pair-level resolution for module-level recovery.

## Prompt (for the image model)

> A clean, publication-quality graphical abstract for a Nature-family computational biology paper, wide 16:9 layout, white background, flat vector style, no 3D, no photorealism, minimal text in a sans-serif font, four panels left to right connected by thin arrows.
>
> Panel 1, "Expression compendium": a compact heat-map grid (genes as rows, many experiments as columns) in muted blue-to-orange, with small icons beneath for "bacteria" and "human cells" to show both E. coli and human data.
>
> Panel 2, "CLR = context-normalized attention": one gene node in the centre with lines to many other genes; each line's thickness shows mutual information. A small bell curve behind the centre gene labelled "this gene's background" shows that each gene's scores are standardized against its own context. A side label in a small grey box: "dependence ≈ query·key; context normalization ≈ scaled logits". Caption: "CLR (2007): attention before transformers".
>
> Panel 3, "Completing the attention": the same gene with a softmax that concentrates weight on a few partners (thick arrows to about 10 genes, faint to the rest), then a multi-step walk drawn as concentric halos expanding over a small network, labelled "depth t" with a short logarithmic scale bar (t = 1, 3, 9, 27). A tiny note: "no learned parameters".
>
> Panel 4, "What it recovers": two small network clusters side by side. Left, a large, loose regulon (a hub transcription factor with many targets) highlighted in aqua with a check mark, labelled "large regulatory programs ↑". Right, a tight protein complex (small dense clique) in grey with a neutral mark, labelled "tight complexes: CLR better". Below, a slim bar strip comparing "CLR attention", "CLR", "GENIE3", "PIDC" with the first bar slightly longer, and a placeholder bar labelled "learned attention (scFM)".
>
> Colour palette (the paper's ggplot theme, analysis/plot_theme.R): blue #2a78d6 for dependence (mutual information), orange #eb6834 for context normalization (CLR), violet #6250d6 for the softmax, magenta #e87ba4 for depth; aqua #1baf7a for the CLR-attention bar and highlighted regulon; grey #8a8984 for competitors and neutral elements; ink #0b0b0b for text and #e7e6e2 for gridlines. Consistent thin dark-grey outlines. Generous white space. No logos, no real author names, no journal branding.

## Notes for the illustrator

- Keep the analogy to transformer attention visible but secondary (panel 2 side box). The audience is Nature Machine Intelligence readers who know transformers.
- Do not imply that attention wins on complexes or on single cells: panel 4 shows the trade-off honestly.
- Panel 4's bar lengths are schematic; take the final values from the method table (`docs/whitepapers/sota_table.md`) once reviewed.
- Colours match the ggplot theme in `analysis/plot_theme.R` (orange = context normalization throughout the paper).
