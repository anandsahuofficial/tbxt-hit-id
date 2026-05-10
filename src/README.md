# `src/` - Pipeline Source

This is the working pipeline. Modules are grouped by stage; the
[`docs/methodology.md`](../docs/methodology.md) describes the
scientific design and per-signal failure modes, this README maps
each module to the stage it implements.

## Layout

```
src/
├── enumeration/    - building one-pot virtual catalogs from reactions + reagents
├── pipeline/       - the six orthogonal scoring signals
├── filters/        - seven-criterion strict gate (PAINS, motifs, lead-likeness, novelty, solubility)
├── ranking/        - four-tier classifier
└── viz/            - pose / structure rendering helpers
```

## Stage 0 - enumeration

Build the candidate pool. Either drop in your own pool CSV, or use
these to enumerate one from a building-block library + reaction set.

| Module | Role |
|---|---|
| [`enumeration/enumerate_analogs.py`](enumeration/enumerate_analogs.py) | Generate Tanimoto-similar analogs of seed compounds |
| [`enumeration/generate_proposals.py`](enumeration/generate_proposals.py) | Score-driven proposal generator from a small starting set |
| [`enumeration/enumerate_onepot.py`](enumeration/enumerate_onepot.py) | Apply 7 chordoma-relevant one-pot reactions to a building-block library |
| [`enumeration/onepot_query.py`](enumeration/onepot_query.py) | Query a virtual catalog by SMILES (muni.bio `onepot` integration) |
| [`enumeration/onepot_reachability.py`](enumeration/onepot_reachability.py) | Per-compound reachability score via the 7 one-pot reactions |
| [`enumeration/retrosynth_audit.py`](enumeration/retrosynth_audit.py) | Audit retrosynthetic accessibility for a candidate set |

## Stage 1 - receptor preparation

| Module | Role |
|---|---|
| [`pipeline/prep_receptor.py`](pipeline/prep_receptor.py) | Strip non-protein, add hydrogens, assign Gasteiger charges, write PDBQT |
| [`pipeline/prep_ensemble.py`](pipeline/prep_ensemble.py) | Sample 6 receptor conformations from short MD relaxation |
| [`pipeline/define_pockets.py`](pipeline/define_pockets.py) | Define site F + site A grid centroids and box dimensions |

## Stage 2 - six orthogonal scoring signals

These implement the scoring half of the pipeline (see
[`docs/methodology.md`](../docs/methodology.md) for the design rationale).

| Signal | Module | Output |
|---|---|---|
| Vina ensemble | [`pipeline/dock_vina.py`](pipeline/dock_vina.py) + [`pipeline/dock_ensemble.py`](pipeline/dock_ensemble.py) | min(vina_kcal) across 6 receptor confs |
| GNINA CNN | [`pipeline/dock_gnina.py`](pipeline/dock_gnina.py) | CNN pose score + CNN affinity (pKd) |
| GNINA pose stability | [`pipeline/dock_gnina_multiseed.py`](pipeline/dock_gnina_multiseed.py) | 10-seed top-pose centroid σ |
| TBXT QSAR | [`pipeline/parse_naar_spr.py`](pipeline/parse_naar_spr.py) → [`pipeline/train_qsar.py`](pipeline/train_qsar.py) | Predicted log Kd ± uncertainty |
| Boltz-2 co-folding | [`pipeline/run_boltz.py`](pipeline/run_boltz.py) (+ [`pipeline/_boltz_safeload.py`](pipeline/_boltz_safeload.py)) | Predicted Kd + prob_binder |
| MMGBSA refinement | [`pipeline/run_mmgbsa.py`](pipeline/run_mmgbsa.py) (+ [`pipeline/run_mmgbsa_md.py`](pipeline/run_mmgbsa_md.py)) | ΔG (kcal/mol) for top-30 picks |
| Paralog selectivity | [`pipeline/paralog_selectivity.py`](pipeline/paralog_selectivity.py) | Per-paralog binding-score gap (16 T-box paralogs) |
| Alchemical FEP (advisory) | [`pipeline/run_fep.py`](pipeline/run_fep.py) | Single-pose ΔΔG; not used to re-rank |

## Stage 3 - consensus

| Module | Role |
|---|---|
| [`pipeline/merge_signals.py`](pipeline/merge_signals.py) | Z-score and align signal directions |
| [`pipeline/consensus.py`](pipeline/consensus.py) | Unweighted-mean consensus + top-quartile flag |
| [`pipeline/analyze_poses.py`](pipeline/analyze_poses.py) | Per-pick pose summary (interaction map, distances) |

## Stage 4 - strict gate (7 criteria)

| Module | Role |
|---|---|
| [`filters/strict_gate.py`](filters/strict_gate.py) | Orchestrates all 7 criteria - outputs a CSV with C1–C7 pass/fail flag columns |
| [`filters/pains_and_motifs.py`](filters/pains_and_motifs.py) | C5 - PAINS A/B/C + reactive-motif SMARTS |
| [`filters/onepot_membership.py`](filters/onepot_membership.py) | C1 - onepot.ai 100% catalog match via muni.bio |

The C2–C4 + C7 criteria are implemented inline inside
`strict_gate.py` because they are short SMARTS / descriptor checks
that don't merit their own modules.

## Stage 5 - tier classification + ranking

| Module | Role |
|---|---|
| [`ranking/tier_classify.py`](ranking/tier_classify.py) | Assigns T1_GOLD / T2_SILVER / T3_BRONZE / T4_RELAXED per [`docs/tier_definitions.md`](../docs/tier_definitions.md); sorts strict-first then by Boltz Kd |

## Stage 6 - visualization

| Module | Role |
|---|---|
| [`viz/render_poses.py`](viz/render_poses.py) | Render 2D structure + 3D pose images (the PNGs in `slides/renders/`) |

## Running individually

Most modules support `python -m src.<group>.<module> --help`. The
end-to-end orchestrator that ties them together is
[`../examples/reproduce_top4.sh`](../examples/reproduce_top4.sh).

## Importing as a package

Although the modules are designed to be invoked from the command
line, the package is importable:

```python
from src.ranking.tier_classify import classify_dataframe, summarize
from src.filters.strict_gate import _evaluate_row, COVALENT_SMARTS
```

## What's not in `src/`

- `src/` contains the scoring + filtering + ranking pipeline. The
  *evaluation* of that pipeline (`results/`), the *narrative* of
  what was learned (`slides/`), and the *operational scaffolding*
  (`setup/`, `tools/`) live in their own top-level directories.
- Lead-only / working-repo-internal scripts (member-data uploaders,
  internal git tooling, on-day playbooks) are not migrated.
