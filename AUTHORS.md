# Authors and Contributions

This file recognizes everyone who contributed to the TBXT Hit
Identification Hackathon project (Pillar VC, Boston, 2026-05-09).
Roles below are grouped by the principal nature of each person's
contribution.

## Project lead

**Anand Sahu** &nbsp;`anandsahuofficial@gmail.com`

Conceptualized the project strategy and the multi-signal pipeline,
designed the pipeline architecture, made the methodology and
prioritization decisions, integrated the multi-signal consensus,
selected the final picks, coordinated the team, and delivered the
live demo at the hackathon.

Specific responsibilities:

- Idea conception and strategy design: target framing (G177D variant
  + site F as the chordoma-selective pocket), the six-signals-with-
  orthogonal-failure-modes design principle, the strict gate +
  four-tier ranking methodology
- 6-signal multi-method consensus design (Vina ensemble, GNINA CNN
  pose + pKd, TBXT-specific QSAR on 650 measured Naar SPR Kd,
  Boltz-2 co-folding, MMGBSA implicit-solvent refinement, and
  T-box paralog selectivity)
- Ran the full simulation matrix across the pipeline (docking,
  co-folding, MMGBSA, paralog selectivity) on top of the team's
  independent runs
- Overnight HPC variant pipelines (onepot generation, full-pool Boltz,
  receptor ensemble dock, MMGBSA + alchemical FEP, site-G dock)
- Onepot.ai catalog validation strategy and integration with the
  muni.bio CLI `onepot` tool plus Rowan platform engagement
- Convergence audit and the 100%-onepot non-covalent swap of all 4
  final picks under the strict catalog-membership criterion
- All judges-facing submission deliverables and the live demo

## Simulation & data generation

These team members ran simulations and generated the data that built
our case across the 6-signal pipeline:

- **Rabia** - independent onepot-style enumeration as a parallel run
  for cross-checking the lead's reaction enumeration
- **Mark M** - multi-seed GNINA dock at site F (10 seeds × 570
  compounds) validating pose stability of the locked picks
- **Jack** - independent local Boltz-2 full-pool run on the 570
  compounds; cross-validated the alternate-engine Boltz Kd
  predictions to within 1.01–1.13× on all 4 final picks
- **Jake Weiss** - additional simulation runs supporting the
  multi-engine cross-validation evidence base

## Analysis & event-time build

These team members analyzed results and built the submission during
the event:

- **Pridhi**
- **Kemal Özkırşehirli**
- **Zankhana Mehta**
- **Zeynep Gülen Erkoç**
- **Lihua Yu**
- **Arup**

## AI implementation assistance

Substantial portions of code, scripts, and prose were drafted with AI
assistance (Claude) under **Anand Sahu's** guidance. All scientific
decisions, parameter choices, judgment calls, and final outputs were
reviewed and accepted by Anand, who is responsible for all outcomes
of this work. Earlier commits' `Co-Authored-By: Claude` trailers were
removed at the project lead's request to consolidate the audit trail
under the human author.

This is the modern norm for compute-heavy projects with AI tooling
and is documented here for transparency.

## Acknowledgments

Thanks to the platform, tooling, and venue partners that made this
work possible:

- **[muni.bio](https://muni.bio)** - `onepot` catalog membership tool,
  CLI, and tooling credits ([CLI docs](https://muni.bio/cli) ·
  [docs](https://muni.bio/docs))
- **[Rowan](https://labs.rowansci.com)** - ADMET, docking, and
  pose-analysis MD platform
- **[onepot.ai](https://www.onepot.ai)** - virtual catalog and
  one-pot synthesis library
- **[Pillar VC](https://www.pillar.vc/)** - Boston event venue
- **[TBXT Hackathon](https://tbxtchallenge.org)** - organizers,
  mentors, and TEPs

Most importantly, thank you to the organizers and the chordoma
research community for treating this disease as a target worth this
much care.
