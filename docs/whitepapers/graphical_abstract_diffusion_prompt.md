# Graphical abstract: prompt for diffusion image models (draft only)

*Drafted 2026-10-22. Companion to `graphical_abstract_prompt.md` (the brief) and to the code-drawn figure `figures/graphical_abstract.{png,pdf}` (`analysis/fig_graphical_abstract.R`).*

**Status: this produces a draft for an illustrator. It must not be submitted.** Nature Portfolio journals do not accept generative-AI images in published content (narrow exceptions do not apply here). Use the output only for layout and mood exploration. The publishable version is the code-drawn figure, or an illustrator's redraw of it. Diffusion models also cannot render reliable text or true bar lengths: every label and number must come from the code-drawn figure (`docs/whitepapers/sota_table.md`).

## Positive prompt

```
scientific graphical abstract, wide 16:9 infographic, white background, four rounded panels in a row left to right connected by small grey chevron arrows, flat vector illustration, thin dark grey outlines, generous white space, minimal sans-serif labels.
Panel 1: a compact heatmap grid of small squares in shades of a single blue, below it a simple rod-shaped bacterium icon with one flagellum and a round human cell icon with a nucleus.
Panel 2: a star graph, one central node with an orange ring connected by blue lines of varying thickness to twelve small white nodes; below it a smooth orange bell curve with short blue tick marks along its base, one tick far in the right tail; a light grey rounded box at the bottom.
Panel 3: the same central orange-ringed node with about ten thick violet arrows pointing outward and many faint thin grey lines to other nodes; below it a small node-link network with three concentric soft pink-magenta dashed halos expanding from one node; a short horizontal pink scale bar with four evenly spaced ticks; a rounded outlined pill at the bottom.
Panel 4: left, a hub node in green with many spoke edges to small green nodes and a green circular check mark badge; right, a small dense grey pentagon clique with all edges drawn and a grey circular minus badge; below, a slim horizontal bar chart with four bars (green longest, two grey, one orange) and one empty dashed placeholder bar.
Palette: blue #2a78d6, orange #eb6834, violet #6250d6, magenta #e87ba4, green #1baf7a, grey #8a8984, near-black #0b0b0b text, light grey #e7e6e2 borders.
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

- The output is a mood/layout draft only and cannot be submitted to Nature Portfolio journals; redraw everything as vector art.
- Keep colour semantics from `analysis/plot_theme.R`: orange = context normalization (CLR) throughout; blue = mutual information; violet = softmax; magenta = depth; green (aqua) = CLR attention; grey = competitors and neutral elements.
- Panel 4 must stay honest: CLR attention leads on large *E. coli* regulons (M3D strong-evidence AUPR / base rate: attention 1.68, GENIE3 1.51, PIDC 1.45, CLR 1.44), while CLR is better on tight complexes (K562 CORUM: attention 3.51, CLR 6.39). The "learned attention (scFM)" bar stays an empty dashed placeholder marked "pending" until that comparison is run.
- Bar lengths come from the table, not from the generated image.
