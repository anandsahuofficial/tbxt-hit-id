# Tier Definitions — Four-Tier Ranking

After every candidate has cleared the [seven-criterion strict gate](filter_chain.md)
and carries a [six-signal score vector](methodology.md), it is
classified into one of four tiers. Tier is the primary ranking key;
within a tier, compounds are ordered by predicted Boltz-2 Kd
ascending.

## Tier rules

| Tier | Hard gate | Lead-like | Solubility | Boltz Kd | Risk profile |
|---|:---:|:---:|:---:|---:|:---:|
| **T1 GOLD** | ✓ all 7 | ✓ | ✓ ESOL > -5 | ≤ 5 µM | low chem AND low supplier |
| **T2 SILVER** | ✓ all 7 | ✓ | ✓ ESOL > -5 | ≤ 10 µM | any |
| **T3 BRONZE** | ✓ all 7 | ✓ | borderline (DMSO @ 10 mM ok) | ≤ 50 µM | any |
| **T4 RELAXED** | ✓ all 7 | ✓ | borderline | ≤ 100 µM | any |

In the 570-compound pool:

| Tier | Count | Notes |
|---:|---:|---|
| T1 GOLD | **0** | empty by design — see "T1 GOLD is empty" below |
| T2 SILVER | **16** | the cleanest tier; most picks pulled from here |
| T3 BRONZE | **89** | the bulk of the pool |
| T4 RELAXED | **32** | tail tier, used only when broader chemotype diversity is needed |
| **Total** | **137** | strict-pass set |

## T1 GOLD is empty (by design, not by accident)

T1 GOLD requires three things simultaneously:

1. Boltz-2 predicted Kd ≤ 5 µM
2. Low chemistry risk (single-step or well-precedented synthesis)
3. Low supplier risk (in-stock, fast turnaround)

In our 570-compound pool, no compound satisfies all three. The
strongest predicted binders (Boltz Kd 3.2 µM and 3.7 µM) are both
in T2 SILVER because at least one of the two risk axes is medium.

This is **surfaced honestly** rather than relaxed. We could have
loosened any of the three criteria to populate T1, but doing so
would mean the tier no longer carries a real signal.

The empty T1 is a methodological feature: an honest acknowledgment
that the 570-compound novelty-filtered pool has no compound that is
simultaneously ultra-potent, ultra-clean to synthesize, and immediately
sourceable. Future pools (different building blocks, different
catalog) might populate T1; this one does not.

## T2 SILVER (the working tier — 16 compounds)

The SILVER tier is the day-to-day working tier. SILVER compounds:

- Pass every hard gate (C1–C7)
- Meet the lead-like ideal (C4)
- Are predicted soluble (ESOL > -5)
- Have a Boltz-2 Kd ≤ 10 µM
- Carry any risk profile (low or medium)

All SILVER picks are catalog-resident, non-covalent, novel, and
within the chordoma-rule chemistry envelope. The two strongest
hackathon picks (`FM002150_analog_0083` Boltz Kd 3.2 µM,
`FM001452_analog_0104` Boltz Kd 3.7 µM) are both T2.

## T3 BRONZE (the broader tier — 89 compounds)

BRONZE is structurally identical to SILVER except that the
solubility criterion is relaxed to "DMSO-soluble at 10 mM" rather
than "aqueous-soluble at 50 µM". For SPR this still works — the DMSO
stock dilutes 1:200 into running buffer — but precipitation risk is
non-trivial for the lower end of the BRONZE solubility distribution.

BRONZE also relaxes the Kd ceiling to 50 µM, capturing weaker but
still useful starting points. Two of the four hackathon picks
(`FM001452_analog_0201` Boltz Kd 8.16 µM and `FM001452_analog_0171`
Boltz Kd 8.32 µM) are formally BRONZE — included because they
contribute chemotype diversity (urea linker, pyridyl selectivity
probe) that the SILVER picks alone do not provide.

## T4 RELAXED (the tail — 32 compounds)

RELAXED compounds pass every hard criterion but carry a Boltz Kd
between 50 and 100 µM. They exist in the ranking for one reason
only: chemotype coverage. If the SILVER and BRONZE tiers happen to
cluster around a single scaffold family, T4 picks can be promoted to
broaden the experimental program submission.

In the actual hackathon submission no T4 picks were used; T4 is
available for the experimental program first-batch submission if the
team wants to widen scaffold diversity.

## Within-tier ranking

Within each tier, compounds are sorted by Boltz-2 predicted Kd
ascending. Ties (same Kd to two decimals) are broken by:

1. Lower chemistry risk first
2. Lower supplier risk second
3. Lower onepot.ai cost third

This is a deterministic ranking — re-running the pipeline produces
the same ordering.

## Tier transitions

A compound's tier is fixed by the rules above and does not move.
There is no manual promotion or demotion. The only way to change
a tier population is to change the underlying scoring or filter
rules — and any such change must be applied to the full 570 pool
to remain auditable.

## Why four tiers (and not three or five)

- **One tier** would collapse the honest distinction between
  "ready for SPR today" (SILVER) and "ready for SPR if you accept
  solubility risk" (BRONZE) — losing decision-relevant information.
- **Two tiers** (pass / fail) hides the Kd potency dimension.
- **Three tiers** (GOLD / SILVER / BRONZE) drops the RELAXED tier
  and forces weaker-but-novel compounds out of the ranking entirely.
- **Five or more tiers** introduces distinctions the pipeline cannot
  reliably support — the signal-to-noise of public methods at the
  Kd extremes (sub-µM and > 100 µM) is too low to draw further
  meaningful sub-categories.

Four is the smallest number that preserves all decision-relevant
information without forcing the data to support distinctions it
cannot.

## How tiers feed the final picks

The hackathon submission pulled:

- 2 picks from T2 SILVER (ranks 1, 2)
- 2 picks from T3 BRONZE (ranks 11, 22) — chosen for chemotype
  diversity, not for Kd

The 20 additional picks (ranks 5–24) for the experimental program
first batch were drawn:

- 13 from T2 SILVER (ranks 5–17)
- 7 from T3 BRONZE (ranks 18–24)

with no T4 RELAXED picks included in the 24-compound submission
(but available in the curated CSVs for follow-up batches).
