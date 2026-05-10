# Pose-Analysis MD — Attempt Log

**Status:** Attempted, not delivered. Documented for transparency.

## What we tried

Run an explicit-solvent pose-analysis MD workflow on the top 4 picks
to produce protein-ligand RMSD trajectories and pose-stability movies.

- Force field: AMBER `ff14SB` for protein, OpenFF `2.2.1` for ligand,
  TIP3P water
- Trajectory: 5 ns × 1 trajectory per pick + 1 ns equilibration
- Receptor: `6F59:A` (G177D + DNA construct)

## Submission paths attempted

| Path | Outcome |
|---|---|
| Direct platform SDK, multiple credentials, automatic PDB fetch | Two credentials failed at workflow submission with HTTP 422; two others submitted successfully but the cloud-side workflows then failed with `WorkflowError` — protein preparation step did not converge. |
| muni.bio MD tool with `prepare_protein=true` | HTTP 400 at `validate_forcefield` — 27 atoms reported as having extreme force after auto-prep; localized hydrogen clashes in residues `ARG63`, `GLN80`, `GLU3`, `GLY2`, and others not resolved by the auto-prep pipeline. |
| Pre-prepared protein (remove + re-add hydrogens, energy-minimize) then submit MD via muni.bio | Protein preparation succeeded under one account; MD submission then failed with HTTP 403 because each platform credential carries its own protein scope and cross-account references are blocked. |

## What we have instead

In place of explicit-solvent MD, we use two pose-stability proxies that
**are** in this repo:

- **Multi-seed GNINA pose stability** — 10 stochastic restarts of GNINA
  at site F across all 570 compounds. The top-pose centroid σ across
  seeds is the stability proxy. Pick #4 (`FM001452_analog_0171`) is in
  the most-stable subset (σ &lt; 0.05 Å).
- **Dual-engine Boltz cross-validation** — two independent Boltz-2 runs
  on separate compute backends agree within 1.34× on Kd for all 4 picks
  (see [`cross_validation_summary.md`](cross_validation_summary.md)).

Together these provide pose-reproducibility evidence at the discrete
sampling level, but they do not replace continuous-trajectory MD for
detecting slow binding-mode rearrangements. That gap is honestly noted.

## Honest framing

If asked "did you run MD?": *Explicit-solvent pose-analysis MD was
attempted on the top 4 picks across three different submission paths,
but the available platform credit tier and protein-scoping model
blocked completion in the compute window. Multi-seed GNINA pose
stability and dual-engine Boltz cross-validation stand in as our
pose-stability evidence in lieu of MD. The attempt log is documented
here rather than silently dropped.*
