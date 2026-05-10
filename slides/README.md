# Slides

Judges-facing deck delivered at the **TBXT Hit Identification
Hackathon** (Pillar VC, Boston, 2026-05-09).

| File | What it is |
|---|---|
| `slides.pdf` | The rendered 19-slide deck (open this) |
| `slides.md` | Marp-style Markdown source |
| `architecture.png` | Pipeline overview graphic (referenced in slide 6) |
| `renders/` | 2D structures + 3D pose images for the 4 picks |

## Structure

The deck is organized around the three judging axes:

1. **Scientific rationale** — target biology, six-signal pipeline,
   pipeline architecture, seven-criterion strict filter
2. **Compound quality** — picks at a glance + a per-pick deep-dive
   slide for each of the four picks (with 2D + 3D pose renders)
3. **Hit-ID judgment** — cross-validation, tradeoffs we made

Closes with reproducibility, a conclusion slide, and a Q&A slide
carrying the InChIKey table for the four picks.

## Render to PDF

The repo includes a small renderer that converts the Marp-style
Markdown to PDF via headless Chromium (no Node / Marp CLI needed):

```bash
pip install markdown playwright
playwright install chromium
python tools/render_slides.py slides/slides.md
# → slides/slides.pdf
```

Alternatively, with the official Marp CLI:

```bash
npx @marp-team/marp-cli slides/slides.md -o slides/slides.pdf
```

## Underlying methodology

For the technical detail behind the deck, see:

- [`docs/methodology.md`](../docs/methodology.md) — six-signal
  scoring pipeline
- [`docs/filter_chain.md`](../docs/filter_chain.md) — seven-criterion
  strict gate
- [`docs/tier_definitions.md`](../docs/tier_definitions.md) —
  four-tier ranking
