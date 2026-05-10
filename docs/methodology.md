# Methodology - Six-Signal Consensus Pipeline

This document explains the scoring half of the pipeline: the six
orthogonal signals that turn a virtual catalog of 570 compounds into
a ranked list of TBXT G177D binders. The downstream filter chain is
documented separately in [`filter_chain.md`](filter_chain.md), and
tier rules in [`tier_definitions.md`](tier_definitions.md).

![Pipeline architecture](architecture.png)

## Design principle

Every public scoring method has a known failure mode. Vina rewards
contact-maximizing decoys; CNN re-rankers can hallucinate plausibility
where there is none; co-folding models prioritize binders that look
like their training distribution; QSAR over-fits the chemotypes in
its training set; MMGBSA is sensitive to single-snapshot artifacts;
selectivity scores penalize legitimate on-target binders that happen
to share family-conserved residues.

The pipeline therefore scores every compound on **six signals chosen
to fail in different directions**, and only promotes compounds where
the signals agree. No pick depends on a single score, and no signal
is allowed to veto a candidate that the rest of the stack endorses.

## Target

- **Protein:** human TBXT (Brachyury), 435 aa, T-box DBD residues 42–219.
- **Variant:** G177D (`rs2305089`), allele frequency ~0.42, present in
  &gt; 90% of Western chordoma cases.
- **Receptor structure:** PDB `6F59` chain A - TBXT G177D + DNA. Matches
  the CF Labs SPR construct used by the hackathon experimental program.
- **Pocket of interest:** site F (Y88 / D177 / L42). The variant
  residue D177 is part of the pocket itself, which makes the pocket
  intrinsically variant-selective.
- **Selectivity reference:** the other 16 human T-box paralogs (TBX1,
  TBX2, …, MGA). G177 is 0% conserved across the family, so a
  D177-engaging ligand is structurally precluded from cross-reacting.

## Compound pool

The 570-compound pool is the union of three sources:

1. **Onepot enumeration** - 7 chordoma-relevant one-pot reactions
   (amide, sulfonamide, urea, reductive amination, Suzuki, ether,
   Buchwald–Hartwig) applied to a building-block library, then
   deduplicated and salted.
2. **Naar prior-art neighborhood** - Tanimoto-similar analogs of
   weak Naar SPR hits, capped at < 0.85 similarity to ensure novelty.
3. **TEP-suggested fragments** - the hackathon TEPs nominated a small
   list of fragment scaffolds; we expanded each with single-step SAR.

After deduplication the pool is fixed at exactly 570 unique SMILES,
which becomes the substrate for every downstream signal.

## The six signals

### 1. Vina ensemble (geometric fit + receptor flexibility)

- **Tool:** AutoDock Vina 1.2.x, exhaustiveness 32, num_modes 9.
- **Receptor handling:** 6 receptor conformations sampled from a
  short MD relaxation of `6F59:A` to capture pocket flexibility
  (Y88 / L42 sidechain rotamers, D177 carboxylate orientation).
- **Per-compound score:** `min(vina_affinity)` across the 6
  conformations. The minimum (most favorable) is used so that one
  compatible receptor conformation is enough to qualify.
- **Failure mode:** Vina rewards contact maximization - a fragment
  with a long flexible tail that fills volume can score better than
  a more rigid binder with a real H-bond. Caught downstream by GNINA
  CNN re-ranking and Boltz-2 affinity classification.

### 2. GNINA CNN pose + pKd (Vina-trap detection)

- **Tool:** GNINA 1.x, default CNN scoring (`--cnn_scoring rescore`)
  on Vina-generated poses. Run via a CUDA-12.8 / cuDNN-9 Singularity
  container to handle the SCC HPC's older glibc.
- **Per-compound output:** `cnn_score` (probability the pose is real)
  and `cnn_affinity` (predicted pKd from the CNN). Best-pose values
  per compound.
- **Pose-stability augmentation:** a 10-seed multi-start dock at
  site F across all 570 compounds gives a pose-stability σ for the
  top-pose centroid. σ &lt; 0.05 Å indicates a stable binding mode.
- **Failure mode:** the CNN is trained on protein-ligand crystal
  structures and can over-confidently endorse poses that look
  PDB-like even when the real binding free energy is poor. Caught
  by the QSAR signal (target-specific) and Boltz-2 (independent
  ML pipeline).

### 3. TBXT-specific QSAR (the only on-target signal)

- **Training data:** 650 measured Naar SPR Kd values against TBXT,
  released for the hackathon. The only target-specific affinity data
  available.
- **Models:** Random Forest + XGBoost regression on Morgan-2 (radius 2,
  nBits 2048) ECFP4 fingerprints. Both report on a held-out 20%
  validation split with R² ≈ 0.4–0.55 - modest, as expected for
  small-molecule SAR with limited data.
- **Per-compound output:** ensemble-mean predicted log Kd plus
  per-model standard deviation as an uncertainty estimate.
- **Failure mode:** the training set's chemotype coverage is narrow
  relative to the 570-compound pool, so QSAR predictions for novel
  scaffolds carry higher uncertainty. Caught by docking + co-folding,
  which are scaffold-agnostic.

### 4. Boltz-2 generative co-folding (dual-engine cross-validation)

- **Tool:** Boltz-2 (Wohlwend et al.), receptor + ligand co-folded
  from sequence + SMILES with no docking prior.
- **Engine 1:** local RTX-4090 single-GPU run on the 570-compound
  pool against `6F59:A`.
- **Engine 2:** independent re-run on the BU SCC cluster A100 nodes,
  same compounds, different infrastructure.
- **Per-compound output:** predicted Kd, plus a `prob_binder`
  classification head (binder / non-binder).
- **Cross-validation:** the dual-engine agreement is the primary
  rigor check on the affinity prediction. Final 4 picks agree within
  1.34× across the two engines.
- **Failure mode:** Boltz prefers binders that resemble its training
  distribution and can systematically under- or over-predict for
  out-of-distribution chemotypes. Cross-engine agreement controls
  for run-to-run variance but not for systematic bias; selectivity
  + MMGBSA + QSAR provide orthogonal counterweights.

### 5. MMGBSA implicit-solvent refinement (top 30)

- **Tool:** OpenMM 8.x with GBn2 implicit solvent.
- **Protocol:** top 30 candidates from the consensus-rerun stage are
  re-docked, then refined with a 100 ps MM-GBSA minimization +
  short equilibration. ΔG = ⟨E_complex⟩ − ⟨E_protein⟩ − ⟨E_ligand⟩.
- **Per-compound output:** ΔG and per-component decomposition
  (electrostatic, van der Waals, GB solvation).
- **Failure mode:** single-snapshot MMGBSA is sensitive to the
  starting pose. We mitigate by feeding the GNINA top pose (already
  CNN-validated) and by short-MD averaging rather than single-frame
  scoring.

### 6. T-box paralog selectivity (off-target risk)

- **Tool:** local sequence alignment of TBXT chain A residues 42–219
  against the 16 other human T-box paralogs (TBX1, TBX2, TBX3, TBX4,
  TBX5, TBX6, TBX10, TBX15, TBX18, TBX19, TBX20, TBX21, TBX22, EOMES,
  MGA, T-Box-3 isoforms).
- **Per-compound output:** for each paralog, a "selectivity score"
  is computed as the difference between the candidate's predicted
  binding score against TBXT vs. against the paralog at the
  homology-mapped pocket equivalent.
- **G177 conservation:** 0% across the 16 paralogs. Compounds whose
  binding signal depends on D177 contact are intrinsically selective
  by structure, not just by score.
- **Failure mode:** homology-based pocket mapping is approximate for
  paralogs with low overall sequence identity. Used as a relative
  ranking signal rather than a hard pass/fail gate.

## Consensus aggregation

Each compound carries a six-element score vector. We z-score each
signal across the 570-compound pool, invert the QSAR/Boltz/MMGBSA
signals so that "more positive = better" everywhere, and compute an
unweighted mean. Compounds in the top quartile of this consensus
score are then re-run through Boltz-2 and GNINA at higher exhaustiveness
to confirm - this is the "consensus rerun" step in the architecture
diagram.

We chose unweighted mean over a learned weighting because the
training data for a learned weighting (target-specific calibration
sets) does not exist for TBXT. An unweighted mean is interpretable
and robust to single-signal outliers.

## What this pipeline does not do

- **No free-energy perturbation (FEP) at scale.** Alchemical FEP was
  attempted as a variant on the top 8 picks (`v4 MMGBSA + alchemical
  FEP` overnight HPC variant), but the results are advisory only -
  single-pose FEP is unreliable enough that we did not let it
  re-rank.
- **No experimental MD pose-validation in the final.** We attempted
  Rowan's pose-analysis MD on the top 4 (5 ns explicit-solvent +
  1 ns equilibration); none of the runs completed in the available
  credit/compute window. This is documented honestly rather than
  silently dropped.
- **No direct chordoma-cell phenotypic data.** The pipeline is
  affinity-and-selectivity-driven; cellular activity is the
  experimental program's job.

## Outputs consumed by the filter chain

After scoring + consensus, every compound carries:

- 6 normalized signal scores + a consensus mean
- `prob_binder` (Boltz)
- `pose_sigma` (GNINA 10-seed)
- `predicted_logKd` ± uncertainty (QSAR)
- `vina_min`, `cnn_score`, `cnn_affinity`
- `mmgbsa_dG` (top 30 only)
- `paralog_selectivity_min` (most-similar paralog's score gap)

These are the inputs to the 7-criterion strict filter chain, which
is documented separately in [`filter_chain.md`](filter_chain.md).
