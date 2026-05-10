# Filter Chain - Seven-Criterion Strict Gate

This document defines the seven hard criteria every candidate must
satisfy to be eligible for the final ranked list. The scoring half
of the pipeline is documented in
[`methodology.md`](methodology.md); tier rules in
[`tier_definitions.md`](tier_definitions.md).

The filter is applied **after** scoring, not before, so we can
report a fair "total candidates considered" → "passed all gates"
ratio (570 → 137).

## Why a strict gate (not a weighted score)

Hackathon organizers and the experimental program define hard
constraints (catalog membership, non-covalent, drug-like). Treating
those as soft penalty terms in a ranking score would let a high-
confidence compound that fails one constraint outrank a fully
compliant compound - which is how submissions get disqualified.

A strict pass/fail gate makes "compliant or not" a tier-0 decision,
and ranking happens only over the surviving population.

## The seven criteria

### C1 - onepot 100% catalog match

- **Rule:** Tanimoto similarity = 1.000 against the onepot.ai
  virtual catalog, queried via the muni.bio CLI `onepot` tool.
- **Source of truth:** muni.bio `onepot` returns canonical SMILES,
  catalog ID, list price, chemistry risk, and supplier risk per
  hit. Anything below 1.000 is rejected at this stage.
- **Why:** the experimental prize program will only synthesize
  compounds present in the onepot.ai catalog; near-matches require
  re-routing through medchem and are out of scope for the hackathon
  submission window.
- **Failure mode:** none - the muni.bio API is the authoritative
  source. The only ambiguity is canonicalization, which RDKit
  handles deterministically.

### C2 - strictly non-covalent

- **Rule:** the compound contains no SMARTS patterns matching
  reactive electrophiles known to form covalent bonds with cysteine,
  serine, lysine, or histidine in protein active sites.
- **Excluded patterns:** Michael acceptors (`[CX3]=[CX3]C(=O)`),
  α,β-unsaturated carbonyls (`C=CC(=O)`), boronic acids (`[B]`),
  acrylamides, vinyl sulfones (`C=CS(=O)(=O)`), aldehydes
  (`[CX3H1](=O)[#6]`), alkyl halides (`[CX4][F,Cl,Br,I]`), epoxides
  (`C1CO1`), aziridines, and isocyanates.
- **Why:** the hackathon explicitly disqualifies covalent binders
  to keep the experimental program tractable (covalent binders
  require kinact/Ki kinetics, not equilibrium SPR Kd).
- **Notable enforcement:** the `[B]` (any boron) pattern was used
  rather than `[B;!H0]` because boronic acids are reversible
  covalent binders - even reversible covalency is excluded.

### C3 - Chordoma chemistry rule

- **Rule:** simultaneously
  - MW ≤ 600 Da
  - LogP ≤ 6 (RDKit Crippen)
  - HBD ≤ 6
  - HBA ≤ 12
- **Why:** the chordoma research community has converged on these
  bounds as the reasonable envelope for compounds intended to reach
  intracranial / intraspinal tumor sites without CNS-penetration
  red flags. A lighter version of Lipinski for a target where
  CNS-leaning permeability is desirable.

### C4 - Lead-like ideal

- **Rule:** simultaneously
  - 10 ≤ heavy atoms ≤ 30
  - HBD + HBA ≤ 11
  - LogP &lt; 5
  - rings &lt; 5
  - fused rings ≤ 2
- **Why:** lead-likeness anchors the picks to chemistry that has
  room to grow during hit-to-lead - adding a methyl or a small
  H-bond donor without immediately violating Lipinski. Compounds
  that already saturate the Chordoma rule (C3) leave no headroom
  for SAR.

### C5 - PAINS + forbidden motifs

- **Rule:** no SMARTS hit against any of:
  - PAINS A, B, C lists (Baell & Holloway 2010)
  - acid halides (`C(=O)[F,Cl,Br,I]`)
  - aldehydes
  - diazo (`[N+]#N`)
  - imines (`C=N` not in aromatic ring)
  - polycyclic systems &gt; 2 fused rings
  - long alkyl chains &gt; C8
- **Why:** PAINS catches assay-interference chemotypes that produce
  false-positive signals across orthogonal assays; the additional
  forbidden motifs catch synthetically tractable but pharmacologically
  fragile chemistry.

### C6 - Novelty (Tanimoto < 0.85 to organizer DBs)

- **Rule:** maximum Tanimoto similarity against three reference
  sets is &lt; 0.85:
  - Naar SPR-measured TBXT binders (653 compounds)
  - TEP-suggested fragment list (curated by hackathon TEPs)
  - `prior_art_canonical` - public TBXT/Brachyury inhibitor
    literature compounds, canonicalized via RDKit
- **Why:** the experimental program rewards novel chemotypes, not
  Naar-similar analogs. 0.85 is the hackathon-defined threshold
  above which a compound is considered a Naar lookalike.

### C7 - ESOL log S > -5 (predicted solubility)

- **Rule:** ESOL-predicted log(S) &gt; -5, where S is in mol/L.
  This corresponds to ~10 µg/mL aqueous solubility - sufficient for
  DMSO @ 10 mM stock dilution into aqueous SPR buffer at 50 µM
  working concentration without precipitation.
- **Why:** SPR assays fail silently when the analyte precipitates,
  producing false-negatives that destroy the experimental program's
  signal. ESOL is a fast proxy that catches obvious solubility
  failures before wet-lab time is spent.
- **Caveat:** ESOL is a regression on a small training set and is
  unreliable for unusual scaffolds. The T3 BRONZE tier intentionally
  relaxes this criterion (see [`tier_definitions.md`](tier_definitions.md)).

## Composition of the strict-pass set (137 compounds)

After applying all seven criteria to the 570-compound pool:

| Filter point | Compounds remaining |
|---|---:|
| Initial pool | 570 |
| After C1 (onepot 100%) | ~290 |
| After C2 (non-covalent) | ~280 |
| After C3 (Chordoma rule) | ~250 |
| After C4 (lead-like ideal) | ~210 |
| After C5 (PAINS + forbidden) | ~180 |
| After C6 (Tanimoto < 0.85 novelty) | ~155 |
| After C7 (ESOL solubility) | **137** |

The "approximate" counts at intermediate steps are because criteria
are applied as a single set logically, not sequentially - the table
is illustrative of the contribution of each gate, not a literal
serial filter trace.

## Two-pass tightening

The filter chain is applied twice:

1. **Pre-scoring (C1–C2 only):** removes catalog-absent and
   covalent compounds before expensive scoring runs, saving roughly
   half the GPU time.
2. **Post-scoring (all 7):** applied to the full 570 after scoring
   so that the final 137 are guaranteed to be both compliant **and**
   to have full per-signal scores attached.

The 137 strict-pass set is the substrate for tier classification
([`tier_definitions.md`](tier_definitions.md)).

## What the filter chain does not do

- **Does not enforce stereochemistry.** Single SMILES per compound
  is treated as the canonical form; chiral submissions are deferred
  to medchem.
- **Does not enforce salt form.** Salt stripping happens once at
  pool ingest; salts are not re-attached for SPR.
- **Does not enforce supplier delivery time.** muni.bio reports
  supplier risk, but it is surfaced as metadata, not as a hard gate.
- **Does not penalize molecular weight against catalog price.**
  Cost-and-risk are reported per pick but are not part of the
  pass/fail gate.

These deferred concerns are addressed at tier-classification time
(BRONZE picks call out solubility risk; risk pills on slides surface
chem and supplier risk).
