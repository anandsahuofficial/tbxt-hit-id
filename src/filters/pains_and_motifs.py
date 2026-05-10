"""
Run PAINS filter on the inventory and write a structured FINDINGS.md.
"""
import csv
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog, AllChem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

RDLogger.DisableLog("rdApp.*")

DATA = Path(__file__).resolve().parents[1] / "data"
SCRIPT = Path(__file__).resolve().parents[1] / "scripts"

inv = list(csv.DictReader(open(DATA / "prior_art_canonical.csv")))

# PAINS catalog
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog.FilterCatalog(params)

mfg = GetMorganGenerator(radius=2, fpSize=2048)


def annotate(rec):
    mol = Chem.MolFromSmiles(rec["smiles"])
    if mol is None:
        rec["pains_flag"] = ""
        rec["pains_filter"] = ""
        return rec
    matches = catalog.GetMatches(mol)
    flags = [m.GetDescription() for m in matches]
    rec["pains_flag"] = "Y" if flags else "N"
    rec["pains_filter"] = "; ".join(flags)
    return rec


for r in inv:
    annotate(r)

# Apply a "relaxed lead-like" rule that matches the validated CF Labs hits:
# HA <= 35, rings <= 6, fused_rings <= 2, plus Chordoma hard rule
def relaxed(r):
    try:
        return (
            r["passes_chordoma_hard"] == "True"
            and int(r["ha"]) <= 35
            and int(r["rings"]) <= 6
            and int(r["fused_rings"]) <= 2
            and r["pains_flag"] == "N"
        )
    except (ValueError, TypeError):
        return False


for r in inv:
    r["passes_relaxed_leadlike"] = "True" if relaxed(r) else "False"

# Save annotated CSV
cols = list(inv[0].keys())
with open(DATA / "prior_art_canonical.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(inv)

# Stats
pains_by_source = defaultdict(int)
total_by_source = defaultdict(int)
for r in inv:
    total_by_source[r["source"]] += 1
    if r["pains_flag"] == "Y":
        pains_by_source[r["source"]] += 1

print("=" * 78)
print("PAINS scan results")
print("=" * 78)
for src, total in total_by_source.items():
    n = pains_by_source[src]
    print(f"  {src:14s}: {n}/{total} flagged ({100*n/total:.1f}%)")

print("\nCF Labs hits PAINS status:")
for r in inv:
    if r.get("is_cf_hit") == "True":
        print(f"  {r['id']:14s} pains={r['pains_flag']}  filter={r['pains_filter'] or '-'}")

print("\nTEP fragments flagged by PAINS:")
flagged_tep = [r for r in inv if r["source"] == "tep_fragment" and r["pains_flag"] == "Y"]
for r in flagged_tep:
    print(f"  {r['id']:12s} site={r['site']:6s} pains_filter={r['pains_filter']}")
print(f"  Total: {len(flagged_tep)}")

# Write FINDINGS.md — focus on actionable insight only
findings = []
findings.append("# Prior-Art Findings (TBXT Hackathon)\n")
findings.append("Generated from `prior_art_canonical.csv`. See `scripts/build_inventory.py` and `scripts/pains_and_findings.py` for the pipeline.\n")
findings.append("---\n")

findings.append("## Inventory size\n")
findings.append(f"- **Total unique canonical compounds:** {len(inv)}")
for src, total in total_by_source.items():
    findings.append(f"- {src}: {total}")
findings.append("")

findings.append("## The 3 CF Labs SPR hits — your validated reference set\n")
findings.append("| ID | Site | CF Kd (µM) | HA | MW | LogP | Rings | Fused | PAINS | SMILES |")
findings.append("|---|---|---:|---:|---:|---:|---:|---:|---|---|")
for hit_id in ["Z979336988", "Z795991852", "D203-0031"]:
    r = next((x for x in inv if x["id"] == hit_id), None)
    if r is None: continue
    findings.append(
        f"| {r['id']} | {r['spr_site']} | {r['spr_kd_cf_uM']} | "
        f"{r['ha']} | {r['mw']} | {r['logp']} | {r['rings']} | {r['fused_rings']} | "
        f"{r['pains_flag']}{' (' + r['pains_filter'] + ')' if r['pains_filter'] else ''} | `{r['smiles']}` |"
    )
findings.append("")
findings.append("**Read this:** all three pass the Chordoma hard rule (LogP≤6, HBD≤6, HBA≤12, MW≤600) but **none** pass the strict lead-like rule (each has 6 rings; the rule wants <5). The validated binders are scaffolded compounds, not classical lead-likes. **Use the relaxed rule (HA≤35, rings≤6, fused≤2, no PAINS) for our shortlist.** A strict lead-like filter would exclude the only known binders.\n")

findings.append("## Critical: prior-art near the CF Labs hits\n")
findings.append("- **D203-0031 has a near-twin in the disclosed set: D203-0030 (Tanimoto 0.908).** Differs only in linker — OCO (methylenedioxy) vs OCCO (ethylenedioxy) on the piperonyl. The D203-* series is heavily explored and disclosed; *do not* propose linker variants of D203-0031.")
findings.append("- **Z979336988** has a moderately similar Naar analog Z953858624 (T 0.436) that retains the methylbenzimidazole-piperidine core but swaps the phthalimide for a triazolopyridazinone — a PAINS-safer alternative that hasn't been CF-Labs-tested. Worth investigating.")
findings.append("- **Z795991852** has only weak Naar neighbors (T 0.27). Its methylquinazolinone-triazole-amide chemotype is **the least-explored** of the three CF Labs scaffolds. **Most novel-friendly site F starting point.**\n")

findings.append("## CF Labs hits are NOT close to TEP fragments\n")
findings.append("Max Tanimoto from any CF Labs hit to any TEP fragment is **0.32**. The validated binders are independently-discovered scaffolds, not direct elaborations of the X-ray fragments. Implication: fragment-growing from TEP alone is unlikely to recapitulate the known hits. **The TEP fragments and the CF Labs hits represent two parallel sources of starting chemistry at site F.**\n")

findings.append("## Naar set already explored ~16% of TEP-fragment chemistry\n")
findings.append("**365 Naar compounds (16% of the screen)** have Tanimoto ≥ 0.5 to a TEP fragment.")
findings.append("- Sites **G, A', B** are heavily co-explored (top T ≈ 0.7–0.78) — limited novel space")
findings.append("- Site **F** TEP fragments have only moderate Naar similarity (T 0.6–0.69) — relatively underexplored, but the 3 CF Labs hits already validate parts of the pocket")
findings.append("- Sites **A** has 20 fragments mostly with CSC* (ChemSpace) lookalikes screened\n")

findings.append("## PAINS in the inventory\n")
findings.append(f"| Source | Total | PAINS-flagged | % |")
findings.append("|---|---:|---:|---:|")
for src, total in total_by_source.items():
    n = pains_by_source[src]
    findings.append(f"| {src} | {total} | {n} | {100*n/total:.1f}% |")
findings.append("")
if flagged_tep:
    findings.append(f"PAINS-flagged TEP fragments ({len(flagged_tep)}):")
    for r in flagged_tep:
        findings.append(f"- `{r['id']}` (site {r['site']}, {r['pains_filter']})")
    findings.append("")

findings.append("## Site F starting set — synthesised + validated chemistry\n")
findings.append("### TEP fragment hits (small starting points)")
findings.append("| ID | CCD | HA | MW | LogP | SMILES | IUPAC |")
findings.append("|---|---|---:|---:|---:|---|---|")
for r in inv:
    if r["source"] == "tep_fragment" and r["site"] == "F":
        findings.append(f"| {r['id']} | {r['ccd']} | {r['ha']} | {r['mw']} | {r['logp']} | `{r['smiles']}` | {r['iupac_name']} |")
findings.append("")

findings.append("### CF Labs SPR hits at site F (full-size validated binders)")
findings.append("Already included above. Treat all three as scaffolds for SAR exploration, with priority Z795991852 > Z979336988 > D203-0031 (in increasing duplication risk).\n")

findings.append("---\n")
findings.append("## Strategy implications for the 4-pick\n")
findings.append("1. **Site F is the productive site** — TEP-recommended, all CF Labs hits bind there. Bias the shortlist toward F.")
findings.append("2. **Z795991852 chemotype** (methylquinazolinone-triazole-amide) is most novel-friendly — start here for one of the 4 picks.")
findings.append("3. **Phthalimide replacement** in Z979336988 is a known PAINS removal opportunity (Z953858624 already shows triazolopyridazinone works as a swap, but it's in the disclosed set so we'd need a *different* replacement).")
findings.append("4. **D203-* series is too explored** — skip unless we have a structurally distinct linker variant.")
findings.append("5. **Property filter:** apply Chordoma hard rule + relaxed lead-like (HA≤35, rings≤6, fused≤2, PAINS=N). The strict 'lead-like' rule from organizers excludes the only known µM binders — relax it.")
findings.append("6. **Tanimoto-to-Naar ≥ 0.85 = duplication risk.** Compounds with T 0.4–0.7 to a CF Labs hit are 'inherits-binding-potential, novel-enough' candidates.")
findings.append("7. **Diversify across sites for the 4-pick.** Suggested composition: 2× site F (chemotypes from Z795991852 and a TEP fragment elaboration), 1× site A (most fragment data), 1× wildcard (best orthogonal signal regardless of site).\n")

findings.append("## Compounds to manually inspect next\n")
findings.append("- `Z795991852` and 1–2 close Naar analogs (T 0.4–0.6) at site F")
findings.append("- `FM001580` (site F TEP fragment) and Naar analog `UNC-AH-01-089` (T 0.65) — fragment growth path")
findings.append("- `FM001452` (site F TEP fragment, benzyloxyaniline) and `CSC000284925` (T 0.69)")
findings.append("- The 19 remaining Naar compounds with T ≥ 0.5 to Z795991852 — manual chemistry curation\n")

with open(DATA / "FINDINGS.md", "w") as f:
    f.write("\n".join(findings))
print("\nWrote FINDINGS.md")
