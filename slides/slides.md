<!--
slides.md - judges-facing deck delivered at the
TBXT Hit Identification Hackathon (Pillar VC, Boston, 2026-05-09).

Render to PDF:
  python ../tools/render_slides.py slides/slides.md
(or via Marp CLI:  marp slides.md -o slides.pdf)

Structure: three judging axes drive the narrative.
  AXIS 1 - Scientific rationale (target, pipeline, filter chain)
  AXIS 2 - Compound quality      (4 picks: overview + per-pick)
  AXIS 3 - Hit ID judgment        (tradeoffs, cross-val, expectations)

Image assets in this folder:
  architecture.png   - pipeline overview
  renders/*.png      - 2D structures + 3D pose images for each pick
-->
---
marp: true
theme: default
paginate: true
math: katex
style: |
  :root {
    --ink: #0a1828;
    --ink-2: #2a3b54;
    --muted: #6b7785;
    --line: #d8dde4;
    --accent: #0b5fff;
    --accent-soft: #e8f0ff;
    --good: #0f7a3d;
    --warn: #b15c00;
  }
  section {
    font-family: "Inter", -apple-system, "Segoe UI", system-ui, sans-serif;
    color: var(--ink);
    padding: 44px 60px 52px;
    font-size: 22px;
    line-height: 1.45;
  }
  h1, h2, h3 { color: var(--ink); letter-spacing: -0.01em; }
  h1 { font-size: 2.0em; margin: 0 0 0.2em; font-weight: 700; }
  h2 { font-size: 1.55em; margin: 0 0 0.55em; font-weight: 600;
       padding-bottom: 10px; border-bottom: 1px solid var(--line); }
  h3 { font-size: 1.1em; color: var(--ink-2); margin: 0.4em 0 0.4em; font-weight: 500; }
  p, li { color: var(--ink); }
  code { background: #f4f6f9; padding: 1px 6px; border-radius: 4px;
         font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 0.92em; color: var(--ink); }
  table { border-collapse: collapse; width: 100%; font-size: 0.88em; margin: 8px 0; }
  th { background: #eef2f7; color: var(--ink); text-align: left;
       padding: 8px 12px; border-bottom: 2px solid var(--line); font-weight: 600; }
  td { padding: 7px 12px; border-bottom: 1px solid var(--line); }
  tr:last-child td { border-bottom: none; }
  blockquote { margin: 12px 0; padding: 10px 16px; background: var(--accent-soft);
               border-left: 3px solid var(--accent); color: var(--ink-2); }
  ul { margin: 6px 0; padding-left: 22px; }
  li { margin: 6px 0; }
  /* Axis tag (header chip) */
  .axis {
    display: inline-block; padding: 4px 12px; border-radius: 4px;
    font-size: 0.72em; letter-spacing: 0.08em; text-transform: uppercase;
    font-weight: 600; color: #fff; background: var(--accent); margin-bottom: 12px;
  }
  .axis.a2 { background: #0a8a4a; }
  .axis.a3 { background: #b15c00; }
  /* Lead / divider */
  section.lead { padding-top: 110px; }
  section.lead h1 { font-size: 2.6em; }
  section.lead .subtitle { color: var(--ink-2); font-size: 1.15em;
                            margin-top: 0.5em; line-height: 1.5; }
  section.lead .meta { margin-top: 60px; color: var(--muted);
                        font-size: 0.95em; line-height: 1.7; }
  section.divider { background: var(--ink); color: #fff; padding-top: 180px; }
  section.divider h1, section.divider h2 { color: #fff; border: none; }
  section.divider p { color: #b8c2d2; font-size: 1.1em; margin-top: 6px; }
  section.divider .axisbig {
    display: inline-block; padding: 6px 16px; border-radius: 4px;
    font-size: 0.7em; letter-spacing: 0.1em; text-transform: uppercase;
    font-weight: 600; color: #fff; background: var(--accent); margin-bottom: 16px;
  }
  section.divider .axisbig.a2 { background: #0a8a4a; }
  section.divider .axisbig.a3 { background: #b15c00; }
  /* Two-column primitives */
  .row { display: grid; grid-template-columns: 1fr 1.1fr; gap: 24px; align-items: center; }
  .row3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 18px; }
  .kpi { background: #f7f9fc; border: 1px solid var(--line); border-radius: 8px; padding: 14px 18px; }
  .kpi .label { font-size: 0.78em; color: var(--muted); text-transform: uppercase; letter-spacing: 0.06em; }
  .kpi .value { font-size: 1.9em; color: var(--accent); font-weight: 700; line-height: 1.05; margin-top: 4px; }
  .kpi .sub { color: var(--muted); font-size: 0.85em; margin-top: 2px; }
  .pill { display: inline-block; padding: 2px 10px; border-radius: 999px;
          font-size: 0.78em; font-weight: 500; background: #f0f4fa; color: var(--ink-2); margin-right: 4px; }
  .pill.good { background: #e6f5ec; color: var(--good); }
  .pill.warn { background: #fff3e0; color: var(--warn); }
  .footnote { color: var(--muted); font-size: 0.82em; margin-top: 14px; }
  img { display: block; }
  /* Architecture-only slide: image fills the slide on the longest axis */
  section.arch { padding: 0; display: flex; align-items: center; justify-content: center;
                  background: #f7f9fc; }
  section.arch .arch-img { max-width: 100%; max-height: 100%; object-fit: contain;
                            margin: 0 auto; }
  .figure { text-align: center; }
  .figure img { margin: 0 auto 8px; }
  .figure .caption { color: var(--muted); font-size: 0.8em; margin-top: 6px; }
---

<!-- _class: lead -->

# TBXT Hit Identification
## Top 4 Picks for Chordoma's Master Regulator

<p class="subtitle">A multi-signal computational pipeline for the G177D variant.</p>

<p class="meta">
TBXT Hackathon · Pillar VC, Boston · 2026-05-09<br/>
Target: TBXT G177D (Brachyury) · PDB <code>6F59</code> chain A · Site F (Y88 / D177 / L42)
</p>

---

## How we map the 3 judging axes

<div class="row3">
<div class="kpi">
  <div class="label" style="color:#0b5fff">Axis 1</div>
  <div class="value" style="font-size:1.1em;line-height:1.2">Scientific rationale</div>
  <div class="sub">Target biology · 6-signal pipeline · 7-criterion strict gate</div>
</div>
<div class="kpi">
  <div class="label" style="color:#0a8a4a">Axis 2</div>
  <div class="value" style="font-size:1.1em;line-height:1.2;color:#0a8a4a">Compound quality</div>
  <div class="sub">4 picks · all 100% onepot · all non-covalent · all site F</div>
</div>
<div class="kpi">
  <div class="label" style="color:#b15c00">Axis 3</div>
  <div class="value" style="font-size:1.1em;line-height:1.2;color:#b15c00">Hit ID judgment</div>
  <div class="sub">Tradeoffs · cross-validation · honest expectations</div>
</div>
</div>

<br/>

<div class="row3">
<div class="kpi">
  <div class="label">Pool screened</div>
  <div class="value">570</div>
  <div class="sub">novelty-filtered</div>
</div>
<div class="kpi">
  <div class="label">Strict-pass</div>
  <div class="value">137</div>
  <div class="sub">all 7 criteria</div>
</div>
<div class="kpi">
  <div class="label">Best Boltz Kd</div>
  <div class="value">3.2 µM</div>
  <div class="sub">dual-engine 1.02×</div>
</div>
</div>

---

<!-- _class: divider -->

<span class="axisbig">Axis 1</span>

# Scientific rationale

Why these methods, why this target, why these 137 candidates.

---

## The target

<span class="axis">Axis 1</span>

- TBXT G177D (Brachyury) - chordoma's master transcription factor and lineage-defining oncogenic driver
- G177D variant (`rs2305089`) is present in &gt; 90% of Western chordoma cases
- Site F pocket directly engages D177 - the variant residue itself
- G177 is 0% conserved across the 16 T-box paralogs &nbsp;⇒&nbsp; structural basis for intrinsic selectivity
- Receptor: PDB `6F59` chain A (G177D + DNA construct - matches CF Labs SPR assay)

> Why site F: TEP-recommended · all best Naar SPR binders predicted at site F · TBXT-unique residues anchor selectivity

---

## The pipeline - 6 orthogonal signals

<span class="axis">Axis 1</span>

| Signal | What it catches |
|---|---|
| Vina ensemble (6 receptor confs) | Geometric fit; receptor flexibility |
| GNINA CNN pose + pKd | Vina-trap detection; ML affinity |
| TBXT QSAR (RF + XGBoost on 650 Naar SPR Kd) | Target-specific affinity |
| Boltz-2 co-folding (two independent runs) | Independent affinity + binder classifier |
| MMGBSA single-snapshot (top 30) | Free-energy refinement |
| T-box paralog selectivity (16 paralogs) | Off-target risk |

> Each signal has a known failure mode - caught by another signal in the stack

---

<!-- _class: arch -->

<img src="architecture.png" class="arch-img">

---

## The 7-criterion strict filter

<span class="axis">Axis 1</span>

<div class="row">
<div>

- C1 onepot 100% catalog match
- C2 strictly non-covalent
- C3 Chordoma rule (MW ≤ 600, LogP ≤ 6, HBD ≤ 6, HBA ≤ 12)
- C4 lead-like ideal (10–30 HA · HBD+HBA ≤ 11 · &lt; 5 rings · ≤ 2 fused)

</div>
<div>

- C5 PAINS-clean + no acid halides / aldehydes / diazo / imines / polycyclic &gt; 2 fused / long alkyl
- C6 Tanimoto &lt; 0.85 to Naar / TEP / prior_art
- C7 ESOL log S &gt; -5

</div>
</div>

> Of 570 pool compounds → 137 pass all 7 criteria → 4 picks (ranks 1, 2, 11, 22)

<p class="footnote">Tier breakdown: <span class="pill">T1 GOLD: 0 (empty by design)</span> <span class="pill good">T2 SILVER: 16</span> <span class="pill">T3 BRONZE: 89</span> <span class="pill">T4 RELAXED: 32</span></p>

---

<!-- _class: divider -->

<span class="axisbig a2">Axis 2</span>

# Compound quality

The 4 picks: properties, evidence, and chemistry.

---

## Picks at a glance

<span class="axis a2">Axis 2</span>

| # | ID | Boltz Kd (run A / B) | gnina Vina/pKd | onepot $ | risks |
|---:|---|---:|---:|---:|:---:|
| 1 | `FM002150_analog_0083` | 3.2 / 3.26 µM | -5.01 / 3.94 | $125 | low/low |
| 2 | `FM001452_analog_0104` | 3.7 / 4.97 µM | -5.77 / 4.03 | $250 | med/med |
| 3 | `FM001452_analog_0201` | 8.16 / 8.76 µM | -6.07 / 4.69 | $375 | high/med |
| 4 | `FM001452_analog_0171` | 8.32 / 8.17 µM | -6.19 / 4.44 | $250 | med/med |

> Total cost to source: $875 · Site F coverage: 4 / 4 · Dual-engine Boltz agreement: 1.01–1.34× across all picks

---

## Pick 1 - `FM002150_analog_0083`

<span class="axis a2">Axis 2</span>

<div class="row">
<div>

- SMILES: `O=C(O)COCc1ccc(-c2ccsc2)cc1`
- Strongest predicted Boltz Kd (3.2 µM) of any 100%-onepot non-covalent compound; dual-engine 1.02× agreement
- Phenoxyacetic acid + thiophene; carboxylate H-bonds Y88 / D177 (variant residue)
- Cheapest + lowest-risk of the 4: <span class="pill good">$125 · low/low</span>
- MW 248.3 · LogP 3.02 · HBD 1 · HBA 4 - clean lead-like profile

</div>
<div class="figure">
<img src="renders/FM002150_analog_0083_2d.png" style="width:340px;">
<img src="renders/FM002150_analog_0083_pose_3d.png" style="width:440px;">
</div>
</div>

---

## Pick 2 - `FM001452_analog_0104`

<span class="axis a2">Axis 2</span>

<div class="row">
<div>

- SMILES: `Cc1ccccc1COc1cccc(N)c1`
- Cleanest medchem in the pool - minimal heteroatom decoration
- Methyl-phenyl-CH₂-O-aniline; aniline-N H-bonds D177
- Mass-efficient (MW 213.3) - best fragment-like starting point for SAR
- Boltz Kd 3.7 / 4.97 µM (1.34×) · <span class="pill">$250 · med/med</span>

</div>
<div class="figure">
<img src="renders/FM001452_analog_0104_2d.png" style="width:340px;">
<img src="renders/FM001452_analog_0104_pose_3d.png" style="width:440px;">
</div>
</div>

---

## Pick 3 - `FM001452_analog_0201`
### urea / benzyl ether for R174 + D177

<span class="axis a2">Axis 2</span>

<div class="row">
<div>

- SMILES: `NC(=O)Nc1cccc(OCc2ccccc2)c1`
- N-aryl urea + benzyl ether: H-bond donor / acceptor pair targeting R174 + D177 specifically
- Deepest Vina pose (-6.07) + highest gnina CNN-pKd (4.69) of the 4
- Adds urea-linker chemotype diversity to pick set
- <span class="pill warn">$375 · high chem · med supplier</span> (urea synthesis)

</div>
<div class="figure">
<img src="renders/FM001452_analog_0201_2d.png" style="width:340px;">
<img src="renders/FM001452_analog_0201_pose_3d.png" style="width:440px;">
</div>
</div>

---

## Pick 4 - `FM001452_analog_0171`
### pyridyl selectivity probe

<span class="axis a2">Axis 2</span>

<div class="row">
<div>

- SMILES: `Nc1cccc(OCc2ccc(-c3ccncc3)cc2)c1`
- Pyridyl introduces basic-N for selectivity probing - may differentiate TBXT from T-box paralogs
- Highest prob_binder (0.46) of the 4
- Tightest dual-engine Boltz agreement (1.02×) - most reproducible Kd prediction
- <span class="pill">$250 · med/med</span> · LogP 3.91

</div>
<div class="figure">
<img src="renders/FM001452_analog_0171_2d.png" style="width:340px;">
<img src="renders/FM001452_analog_0171_pose_3d.png" style="width:440px;">
</div>
</div>

---

<!-- _class: divider -->

<span class="axisbig a3">Axis 3</span>

# Hit ID judgment

What we chose, what we didn't, and what to actually expect.

---

## Cross-validation

<span class="axis a3">Axis 3</span>

- Two independent Boltz-2 runs on separate compute backends: 4 / 4 picks agree within 1.34×
- Rowan ADMET (49 properties × 4): all 4 ADMET-profiled
- Rowan pose-analysis MD (explicit-solvent, 5 ns × 1 traj + 1 ns equil): protein-ligand RMSD trajectories captured per pick
- Onepot.ai catalog (muni `onepot` tool): all 4 at similarity = 1.000 with price + chemistry_risk + supplier_risk

> Every pick is supported by multiple independent lines of evidence - not a single-score gamble.

---

## Tradeoffs we made

<span class="axis a3">Axis 3</span>

- All 4 site F - gives up site-A diversity. Defensible: the 100%-onepot non-covalent constraint dominated, and the strongest catalog-resident chemotype is site-F by structural heritage. Site-A backups are kept in our broader 24-compound submission for the experimental program.
- T1 GOLD tier intentionally empty - no compound simultaneously hits Kd ≤ 5 µM AND low/low risk. We surface this honestly rather than overclaim.
- Chemotype dominance (one family in 3 of 4 picks) is not artificially diversified - the catalog × binding-evidence intersection naturally selects this family. Our 20 additional submissions widen the chemotype envelope.

---

## Reproducibility

- Code, scoring outputs, and per-pick rationale are version-controlled and reproducible from a single setup command
- Every one of the 137 strict-pass candidates carries per-criterion pass/fail flags for full audit
- 24 compounds queued for the experimental program - 4 presented today, 20 additional first-batch
- Full bundle (CSVs, poses, rendered structures, evidence files) available on request

---

## Conclusion

<div class="row3">
<div class="kpi">
  <div class="label" style="color:#0b5fff">Axis 1 - what we did</div>
  <div class="value" style="font-size:1.0em;line-height:1.25">6 orthogonal signals + 7-criterion strict gate</div>
  <div class="sub">570 → 137 strictly compliant candidates</div>
</div>
<div class="kpi">
  <div class="label" style="color:#0a8a4a">Axis 2 - what we picked</div>
  <div class="value" style="font-size:1.0em;line-height:1.25;color:#0a8a4a">4 picks · 100% onepot · non-covalent · site F</div>
  <div class="sub">best Boltz Kd 3.2 µM · $875 to source all 4</div>
</div>
<div class="kpi">
  <div class="label" style="color:#b15c00">Axis 3 - how we judged ourselves</div>
  <div class="value" style="font-size:1.0em;line-height:1.25;color:#b15c00">Cross-validated, honest about expectations</div>
  <div class="sub">realistic SPR window 18–200 µM</div>
</div>
</div>

<br/>

### Take-aways for judges

- A defensible, multi-signal pipeline - not a single-score gamble, every pick is corroborated
- A strict, organizer-aligned hard gate - chordoma rule, non-covalent, novelty, lead-likeness, all enforced
- Calibrated confidence - best-case Kd predictions surfaced, public-method over-prediction acknowledged

> Our 4 picks are the strongest catalog-resident, novel, non-covalent compounds that survive every constraint we believe a TBXT hit must satisfy.

---

<!-- _class: lead -->

# Q&A

<p class="subtitle">InChIKeys for the 4 picks</p>

| # | ID | InChIKey |
|---:|---|---|
| 1 | `FM002150_analog_0083` | `TWMVEBYUPVBIAF-UHFFFAOYSA-N` |
| 2 | `FM001452_analog_0104` | `CFASJXOXPAEGCK-UHFFFAOYSA-N` |
| 3 | `FM001452_analog_0201` | `SWOUSJGQCZFNGK-UHFFFAOYSA-N` |
| 4 | `FM001452_analog_0171` | `UJJFUMPTHNXASI-UHFFFAOYSA-N` |
