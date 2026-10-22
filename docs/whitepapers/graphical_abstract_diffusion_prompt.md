# Graphical abstract: prompt for diffusion image models (draft only)

*Drafted 2026-10-22. Companion to `graphical_abstract_prompt.md` (the brief) and to the code-drawn figure `figures/graphical_abstract.{png,pdf}` (`analysis/fig_graphical_abstract.R`).*

**Status: this produces a draft for an illustrator. It must not be submitted.** Nature Portfolio journals do not accept generative-AI images in published content (narrow exceptions do not apply here). Use the output only for layout and mood exploration. The publishable version is the code-drawn figure, or an illustrator's redraw of it. Diffusion models also cannot render reliable text or true bar lengths: every label and number must come from the code-drawn figure (`docs/whitepapers/sota_table.md`).

## Positive prompt

```
scientific graphical abstract, wide 16:9 infographic, white background, flat vector illustration, thin dark grey outlines, generous white space, minimal sans-serif labels.
Left two thirds: four small rounded panels in a row connected by grey chevrons, grouped under two labelled bands with square brackets: a light grey band over panels 1-2 and a pale violet band over panels 3-4; a thin light grey strip spans the full width under all four panels.
Panel 1: a compact heatmap grid of small squares in shades of a single blue; below it a rod-shaped bacterium icon with one flagellum and a round human cell icon with a nucleus.
Panel 2: a star graph, one central solid black node connected by blue lines of varying thickness to ten small white nodes; below it a smooth orange bell curve with short blue tick marks along its base, one tick far in the right tail.
Panel 3: the same central black node with about ten thick violet arrows pointing outward and many faint thin grey lines to other nodes.
Panel 4: a small node-link network with three concentric soft pink-magenta dashed halos expanding from a central black node; below it a short horizontal magenta scale bar with four evenly spaced ticks and one dot on the third tick.
Right third: one tall panel with two equal horizontal bar charts stacked vertically, each with four bars and a dashed vertical reference line; top chart has a small green hub-and-spoke icon, bars green (longest), grey, grey, orange; bottom chart has a small grey pentagon clique icon, bars green (shortest), orange, grey, grey, with a red arrow pointing left from the orange bar's length to the end of the green bar.
Palette: blue #2a78d6, orange #eb6834, violet #6250d6, magenta #e87ba4, green #1baf7a, grey #8a8984, red #e34948 (only the loss arrow), near-black #0b0b0b, light grey #e7e6e2 borders.
```

## Negative prompt

```
photorealistic, 3D render, isometric, glossy, gradient mesh, drop shadows, bokeh, glow, neon, dark background, textured paper, watercolor, sketch, hand-drawn, clip art, cartoon characters, faces, people, microscope photos, DNA double helix, brain, robot, circuit board, logos, watermark, journal branding, author names, dense paragraphs, gibberish text, misspelled labels, rainbow colormap, red-green contrast, cluttered, overlapping elements, extra panels, 3D bar charts
```

## Style keywords

flat vector, Nature-journal figure style, editorial infographic, Swiss minimal, clean line art, schematic network diagram, scientific illustration, consistent stroke weight, colour-blind-safe palette, high information density with white space.

## Settings

- **Aspect ratio:** 16:9 (e.g. 1920x1080 or 2400x1350; Midjourney `--ar 16:9`; SDXL 1344x768 then upscale).
- **Seed:** fix a seed once a layout works (e.g. try 8 to 16 seeds, keep the one whose four panels are evenly spaced), then vary only the prompt so changes are attributable. Record seed, model and version alongside each draft.
- **Guidance / stylization:** moderate CFG (5 to 7 on SDXL) or low stylization (`--stylize 50` to 100 on Midjourney); high values add decorative effects the negative prompt is fighting.
- **Text:** ask for no text, or treat any text as placeholder; overlay real labels afterwards. Labels and numbers to use are in `graphical_abstract_prompt.md` and the code-drawn figure.
- **Iterate per panel** if the model merges panels: generate each panel at 4:5 and compose them in a vector editor.

## Notes for the illustrator

- The output is a mood/layout draft only and cannot be submitted to Nature Portfolio journals; redraw everything as vector art, using the code-drawn figure as the reference for content and numbers.
- Do not credit CLR with the softmax or the walk. Published CLR (Faith et al. 2007) is only the per-gene z-scoring of MI (panels 1-2, grey band); the softmax and the multi-step walk are this work (violet band). The transformer analogy sits in its own strip under both bands, not in a panel.
- Colour key used in every paper figure (`analysis/plot_theme.R`): blue = dependence (MI); orange = CLR normalization (also the CLR bar); violet = softmax; magenta = depth/walk; green (aqua) = the full method of this work, bars only; grey = competitors; red = the loss marker only; query gene i in neutral black.
- Results must show the trade-off with equal weight: *E. coli* RegulonDB regulons (this work 1.68, GENIE3 1.51, PIDC 1.45, CLR 1.44) and K562 CORUM complexes (this work 3.51, CLR 6.39, GENIE3 6.32, |r| 7.00; 1.8x lower than CLR). Units: precision-recall area as a multiple of random. Bar lengths come from `docs/whitepapers/sota_table.md`, never from the generated image.
