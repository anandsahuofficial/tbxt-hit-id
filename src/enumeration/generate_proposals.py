"""
Strategy 7 — Generative chemistry via BRICS recombination across the prior-art reservoir.

Approach:
  1. Decompose all 42 TEP fragments + 135 Naar disclosed Sheet compounds + 4 priority
     scaffolds into BRICS fragments (synthetically-accessible bond breakages).
  2. Use BRICS.BRICSBuild to recombine into novel molecules, with shuffling enabled.
  3. Filter through:
     - Valid SMILES + sanitization
     - Chordoma hard rule + relaxed lead-like (HA ≤ 35, rings ≤ 6, fused ≤ 2, no PAINS)
     - Tanimoto < 0.85 to all 2274 known compounds (= novel relative to disclosed set)
     - QSAR predicted pKd ≥ 4.0 (predicted µM-range binder via TBXT-specific QSAR)
  4. Output top 200 by QSAR pKd as generative proposals.

Why BRICS-recombination counts as "generative":
  - Outputs are novel SMILES not in any existing database (we filter T < 0.85)
  - Fragments come from both validated TBXT binders (Naar) and crystallographically-
    bound TEP fragments — giving the model a pharmacophore-rich starting set
  - BRICSBuild's bond-rule constraints ensure synthesizable connections (mostly
    amide/ester/aryl-aryl/CH-CH bonds — overlap with Onepot's 7 reactions)

This is not pocket-conditioned (no 3D), but the QSAR scoring is target-specific
TBXT, which is the closest substitute available without installing pocket-conditioned
generative models like Pocket2Mol or DiffSBDD (both require external GitHub installs
+ checkpoint downloads).
"""
import csv
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem, BRICS, FilterCatalog, Lipinski, Descriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")

DATA = Path(__file__).resolve().parents[1] / "data"
QSAR_DIR = DATA / "qsar"
OUT = DATA / "generative"
OUT.mkdir(exist_ok=True)

MIN_QSAR_PKD = 4.0      # ≥ 100 µM predicted Kd
MAX_TANIMOTO_TO_KNOWN = 0.85
MAX_HA = 35
MAX_RINGS = 6
MAX_FUSED = 2
TOP_N = 200
MAX_GENERATE = 30000     # cap BRICSBuild attempts
TIMEOUT_S = 300           # seconds

# PAINS catalog
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

mfg = GetMorganGenerator(radius=2, fpSize=2048)


def featurize(smiles_list):
    """Same featurization as train_qsar.py (Morgan FP + 8 descriptors)."""
    fps = np.zeros((len(smiles_list), 2048), dtype=np.uint8)
    descs = np.zeros((len(smiles_list), 8), dtype=np.float32)
    valid = np.zeros(len(smiles_list), dtype=bool)
    for i, smi in enumerate(smiles_list):
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        try: Chem.SanitizeMol(m)
        except Exception: continue
        fp = mfg.GetFingerprint(m)
        fps[i] = np.frombuffer(fp.ToBitString().encode(), dtype=np.uint8) - ord('0')
        descs[i] = [
            Descriptors.MolWt(m), Descriptors.MolLogP(m),
            Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m),
            Descriptors.TPSA(m), Lipinski.NumRotatableBonds(m),
            m.GetRingInfo().NumRings(), m.GetNumHeavyAtoms(),
        ]
        valid[i] = True
    return np.hstack([fps, descs]), valid


def passes_filters(mol):
    try:
        mw = Descriptors.MolWt(mol); logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol); hba = Lipinski.NumHAcceptors(mol)
        ha = mol.GetNumHeavyAtoms()
    except Exception:
        return False, "desc_fail"
    if logp > 6 or hbd > 6 or hba > 12 or mw > 600: return False, "chordoma_hard"
    if ha > MAX_HA: return False, "too_big"
    rings = mol.GetRingInfo()
    if rings.NumRings() > MAX_RINGS: return False, "too_many_rings"
    # fused rings
    aring = rings.AtomRings()
    fused = 0
    for i in range(len(aring)):
        for j in range(i + 1, len(aring)):
            if len(set(aring[i]) & set(aring[j])) >= 2:
                fused += 1
                break
    if fused > MAX_FUSED: return False, "too_many_fused"
    if _pains_catalog.HasMatch(mol): return False, "pains"
    return True, "ok"


def main():
    print("Loading prior-art reservoir...")
    seed_smiles = []
    seed_sources = []
    # TEP fragments
    with open(DATA / "tep" / "tep_fragments_full.csv") as f:
        for r in csv.DictReader(f):
            if r["smiles"]:
                seed_smiles.append(r["smiles"])
                seed_sources.append(f"tep:{r['fragment_id']}")
    # Naar Sheet
    with open(DATA / "naar" / "naar_sheet_export.csv") as f:
        next(f)
        for row in csv.reader(f):
            if len(row) >= 2 and row[1]:
                seed_smiles.append(row[1])
                seed_sources.append(f"naar_sheet:{row[0]}")
    # Priority scaffolds (4 strong scaffolds)
    priority = {
        "Z795991852": "Cn1c(=O)c2ccccc2n2c(COC(=O)c3cccc(NC(=O)C4Cc5ccccc5O4)c3)nnc12",
        "Z979336988": "Cc1ccc2[nH]c(C3CCCN(C(=O)c4cccc(CN5C(=O)c6ccccc6C5=O)c4)C3)nc2c1",
        "FM001580":   "O=C(O)c1ccccc1OC(F)(F)F",
        "FM001452":   "Nc1cccc(OCc2ccccc2)c1",
        "FM002150":   "O=C(O)COCc1ccccc1",
    }
    for k, v in priority.items():
        seed_smiles.append(v)
        seed_sources.append(f"priority:{k}")
    print(f"  Seed reservoir: {len(seed_smiles)} compounds")

    # Decompose into BRICS fragments
    print("\nBRICS-decomposing seed reservoir...")
    all_fragments = set()
    for smi in seed_smiles:
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        try:
            frags = BRICS.BRICSDecompose(m, returnMols=False, minFragmentSize=3)
            all_fragments.update(frags)
        except Exception:
            continue
    fragment_mols = []
    for fsmi in all_fragments:
        m = Chem.MolFromSmiles(fsmi)
        if m is not None:
            fragment_mols.append(m)
    print(f"  Unique BRICS fragments: {len(fragment_mols)}")

    # BRICSBuild — generate novel molecules
    print("\nBRICS-building novel molecules...")
    t0 = time.time()
    generated = []
    seen = set()
    try:
        for i, m in enumerate(BRICS.BRICSBuild(fragment_mols, scrambleReagents=True, maxDepth=5)):
            if i >= MAX_GENERATE: break
            if time.time() - t0 > TIMEOUT_S:
                print(f"  Timeout reached after {i} attempts")
                break
            if m is None: continue
            try: Chem.SanitizeMol(m)
            except Exception: continue
            smi = Chem.MolToSmiles(m)
            if smi in seen: continue
            seen.add(smi)
            generated.append(smi)
            if len(generated) % 1000 == 0:
                print(f"  ...{len(generated)} unique generated")
    except Exception as e:
        print(f"  BRICSBuild ended: {e}")
    print(f"  Total generated unique SMILES: {len(generated)} in {time.time()-t0:.1f}s")

    # Apply property filters
    print("\nApplying property + PAINS filters...")
    survivors = []
    n_invalid = n_filter = 0
    fail_reasons = {}
    for smi in generated:
        m = Chem.MolFromSmiles(smi)
        if m is None: n_invalid += 1; continue
        ok, why = passes_filters(m)
        if not ok:
            n_filter += 1
            fail_reasons[why] = fail_reasons.get(why, 0) + 1
            continue
        survivors.append(smi)
    print(f"  After property filter: {len(survivors)}/{len(generated)}")
    print(f"  Filter reasons: {fail_reasons}")

    # Tanimoto vs known compounds
    print("\nLoading prior-art canonical inventory for novelty filter...")
    known_fps = []
    with open(DATA / "prior_art_canonical.csv") as f:
        for r in csv.DictReader(f):
            if r["smiles"]:
                m = Chem.MolFromSmiles(r["smiles"])
                if m is not None: known_fps.append(mfg.GetFingerprint(m))
    print(f"  Known fingerprints: {len(known_fps)}")

    print("Filtering by Tanimoto-to-known < {}...".format(MAX_TANIMOTO_TO_KNOWN))
    novel = []
    for smi in survivors:
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        fp = mfg.GetFingerprint(m)
        sims = DataStructs.BulkTanimotoSimilarity(fp, known_fps)
        max_sim = max(sims)
        if max_sim < MAX_TANIMOTO_TO_KNOWN:
            novel.append({"smiles": smi, "max_tanimoto_to_known": round(max_sim, 3)})
    print(f"  After novelty filter: {len(novel)}")

    # QSAR scoring
    print("\nLoading QSAR ensemble (RF + XGB)...")
    with open(QSAR_DIR / "model_rf.pkl", "rb") as f: rf = pickle.load(f)
    xgb_model = xgb.XGBRegressor()
    xgb_model.load_model(str(QSAR_DIR / "model_xgb.json"))

    print("Scoring with QSAR...")
    novel_smis = [r["smiles"] for r in novel]
    X, valid = featurize(novel_smis)
    rf_pred = rf.predict(X)
    xgb_pred = xgb_model.predict(X)
    ens_pred = 0.5 * (rf_pred + xgb_pred)

    for i, r in enumerate(novel):
        r["qsar_pkd_rf"] = round(float(rf_pred[i]), 3)
        r["qsar_pkd_xgb"] = round(float(xgb_pred[i]), 3)
        r["qsar_pkd_ens"] = round(float(ens_pred[i]), 3)
        r["qsar_kd_uM"] = round(10 ** (-float(ens_pred[i])) * 1e6, 2)

    # Apply QSAR threshold
    novel_qsar = [r for r in novel if r["qsar_pkd_ens"] >= MIN_QSAR_PKD]
    novel_qsar.sort(key=lambda r: -r["qsar_pkd_ens"])
    print(f"  After QSAR pKd ≥ {MIN_QSAR_PKD}: {len(novel_qsar)}")

    # Cap at TOP_N
    final = novel_qsar[:TOP_N]
    print(f"\nFinal proposals: {len(final)}")

    # Assign IDs + descriptors + write
    for i, r in enumerate(final):
        r["id"] = f"gen_{i+1:04d}"
        m = Chem.MolFromSmiles(r["smiles"])
        r["mw"] = round(Descriptors.MolWt(m), 1)
        r["ha"] = m.GetNumHeavyAtoms()
        r["logp"] = round(Descriptors.MolLogP(m), 2)
        r["hbd"] = Lipinski.NumHDonors(m)
        r["hba"] = Lipinski.NumHAcceptors(m)
        r["rings"] = m.GetRingInfo().NumRings()

    cols = ["id", "smiles", "qsar_pkd_ens", "qsar_kd_uM", "qsar_pkd_rf", "qsar_pkd_xgb",
            "max_tanimoto_to_known", "mw", "ha", "logp", "hbd", "hba", "rings"]
    out_path = OUT / "generative_proposals.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in final:
            w.writerow({k: r.get(k, "") for k in cols})
    print(f"Wrote {out_path}")

    # Summary
    print("\nTop 10 proposals by QSAR ensemble pKd:")
    for r in final[:10]:
        print(f"  {r['id']:10s}  pKd={r['qsar_pkd_ens']:.2f}  ({r['qsar_kd_uM']:>6.1f} µM)  "
              f"T_known={r['max_tanimoto_to_known']:.2f}  HA={r['ha']}  rings={r['rings']}  "
              f"smi={r['smiles'][:60]}")


if __name__ == "__main__":
    main()
