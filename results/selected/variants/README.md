# Methodological variants

Five overnight HPC variants of the pipeline were run as
cross-validation, each stressing a different methodological choice.
The outputs in this folder are evidence that the pipeline runs
end-to-end across these variants — not the post-strict-gate ranked
list (that is in [`../../top4.csv`](../../top4.csv) and
[`../../all_candidates_tiered.csv`](../../all_candidates_tiered.csv)).

| File | Variant | What it tested |
|---|---|---|
| [`variant1_onepot_friendly_top50.csv`](variant1_onepot_friendly_top50.csv) | One-pot reaction enumeration | Top-50 from a parallel reaction-graph enumeration of building blocks across the 7 chordoma-relevant one-pot reactions |
| [`variant2_full_pool_boltz.json`](variant2_full_pool_boltz.json) | Full-pool Boltz on HPC | Boltz-2 inference on the full 570-pool from a separate HPC backend (cross-validates the local Boltz run) |
| [`variant3_ensemble_consensus_local.json`](variant3_ensemble_consensus_local.json) | Receptor ensemble (local) | Vina + GNINA across 4 receptor conformations sampled from short MD relaxation, run on the local GPU |
| [`variant3_ensemble_consensus_scc.json`](variant3_ensemble_consensus_scc.json) | Receptor ensemble (HPC) | Same as above, independently re-run on HPC for cross-validation |
| [`variant4_mmgbsa_md.json`](variant4_mmgbsa_md.json) | MMGBSA + short MD | Implicit-solvent MMGBSA refinement with short MD averaging on the top-30 consensus picks |
| [`variant4_fep_alchemical.json`](variant4_fep_alchemical.json) | Alchemical FEP | Single-pose alchemical free-energy perturbation on a small test set; advisory only — single-pose FEP is not reliable enough to re-rank |
| [`variant5_site_g_results.json`](variant5_site_g_results.json) | Alternate pocket (site G) | Re-dock at a TBXT-suggested alternate pocket (site G, centroid E48 / E50 / G81 / Y210) to test pocket-choice robustness |

## How variants relate to the final picks

The variants score the **same 570-compound pool** under different
methodological assumptions. The final picks in
[`../../top4.csv`](../../top4.csv) emerge from the post-scoring
seven-criterion strict gate (see
[`../../../docs/filter_chain.md`](../../../docs/filter_chain.md))
applied to the consensus-aggregated scores, with the additional
post-hoc constraint that all picks must be 100% catalog-matched at
onepot.ai (the post-swap subset). Some compound IDs in the variant
JSONs reflect intermediate-stage candidate identifiers and are not
identical to the final FM_* picks; this is expected — the variants
are evidence of pipeline rigor, not the final ranked output.
