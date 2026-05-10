"""
TBXT-specific QSAR model — Strategy 1.

Inputs: data/qsar/naar_kd_dataset.csv  (653 compounds with SMILES + measured pKd)
Features: Morgan fingerprints (r=2, 2048 bits) + 8 RDKit descriptors
Target: pKd = -log10(Kd in M)

Models trained:
  - Random Forest (scikit-learn)
  - XGBoost
  - Ensemble = mean(RF, XGB)

Validation:
  - 5-fold cross-validation (random split)
  - Spearman ρ, Pearson r, MAE, RMSE per fold
  - Out-of-fold predictions for the full set

Held-out testing on:
  - 3 CF Labs hits (Z795991852 10 µM, Z979336988 30 µM, D203-0031 17 µM)
    Note: these ARE in our training set (curated Sheet IDs), so this is leakage check
  - 3 site F TEP fragments (FM001580, FM001452, FM002150) — never SPR-tested at CF Labs
  - 503-compound analog pool (compare rank vs GNINA)

Outputs:
  data/qsar/cv_results.json
  data/qsar/model_rf.pkl, data/qsar/model_xgb.pkl
  data/qsar/predictions_analogs.csv
  data/qsar/predictions_tep_fragments.csv
  data/qsar/QSAR_VALIDATION.md
"""
import csv
import json
import math
import pickle
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy.stats import spearmanr, pearsonr
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")

DATA = Path(__file__).resolve().parents[1] / "data"
QSAR = DATA / "qsar"

# Outlier handling: drop compounds with absurd pKd (>= 7 = sub-100nM, unrealistic for this target)
MAX_PKD_FOR_TRAINING = 7.0
MIN_PKD_FOR_TRAINING = 2.0

mfg = GetMorganGenerator(radius=2, fpSize=2048)


def featurize(smiles_list):
    """SMILES → (fingerprints array, descriptor array, valid mask)."""
    fps = np.zeros((len(smiles_list), 2048), dtype=np.uint8)
    descs = np.zeros((len(smiles_list), 8), dtype=np.float32)
    valid = np.zeros(len(smiles_list), dtype=bool)
    for i, smi in enumerate(smiles_list):
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        fp = mfg.GetFingerprint(m)
        fps[i] = np.frombuffer(fp.ToBitString().encode(), dtype=np.uint8) - ord('0')
        descs[i] = [
            Descriptors.MolWt(m),
            Descriptors.MolLogP(m),
            Lipinski.NumHDonors(m),
            Lipinski.NumHAcceptors(m),
            Descriptors.TPSA(m),
            Lipinski.NumRotatableBonds(m),
            m.GetRingInfo().NumRings(),
            m.GetNumHeavyAtoms(),
        ]
        valid[i] = True
    X = np.hstack([fps, descs])
    return X, valid


def cv_eval(X, y, model_factory, n_folds=5, seed=42):
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    fold_metrics = []
    oof_pred = np.zeros_like(y)
    for fi, (tr_idx, te_idx) in enumerate(kf.split(X)):
        m = model_factory()
        m.fit(X[tr_idx], y[tr_idx])
        pred = m.predict(X[te_idx])
        oof_pred[te_idx] = pred
        rho = spearmanr(y[te_idx], pred).correlation
        r = pearsonr(y[te_idx], pred)[0]
        mae = float(np.mean(np.abs(y[te_idx] - pred)))
        rmse = float(np.sqrt(np.mean((y[te_idx] - pred) ** 2)))
        fold_metrics.append({"fold": fi, "n_test": int(te_idx.size),
                             "spearman": float(rho), "pearson": float(r),
                             "mae_pkd": mae, "rmse_pkd": rmse})
    return oof_pred, fold_metrics


def main():
    # Load data
    print("Loading naar_kd_dataset.csv ...")
    rows = list(csv.DictReader(open(QSAR / "naar_kd_dataset.csv")))
    print(f"  Total compounds: {len(rows)}")

    # Filter outliers
    rows = [r for r in rows if MIN_PKD_FOR_TRAINING <= float(r["pkd"]) <= MAX_PKD_FOR_TRAINING]
    print(f"  After pKd ∈ [{MIN_PKD_FOR_TRAINING}, {MAX_PKD_FOR_TRAINING}] filter: {len(rows)}")

    smiles = [r["smiles"] for r in rows]
    y = np.array([float(r["pkd"]) for r in rows], dtype=np.float32)
    ids = [r["compound_id"] for r in rows]

    print("\nFeaturizing ...")
    X, valid = featurize(smiles)
    X = X[valid]; y = y[valid]; ids = [c for c, v in zip(ids, valid) if v]
    print(f"  Valid molecules: {valid.sum()} / {len(valid)}")
    print(f"  Feature shape: {X.shape}")

    # Cross-validate
    print("\n--- 5-fold CV: Random Forest ---")
    rf_oof, rf_metrics = cv_eval(X, y, lambda: RandomForestRegressor(
        n_estimators=400, max_depth=None, min_samples_leaf=2, n_jobs=-1, random_state=42))
    rf_overall = {
        "spearman": float(spearmanr(y, rf_oof).correlation),
        "pearson": float(pearsonr(y, rf_oof)[0]),
        "mae_pkd": float(np.mean(np.abs(y - rf_oof))),
        "rmse_pkd": float(np.sqrt(np.mean((y - rf_oof) ** 2))),
    }
    print(f"  OOF Spearman ρ: {rf_overall['spearman']:.3f}")
    print(f"  OOF Pearson r:  {rf_overall['pearson']:.3f}")
    print(f"  OOF MAE pKd:    {rf_overall['mae_pkd']:.3f}")
    print(f"  OOF RMSE pKd:   {rf_overall['rmse_pkd']:.3f}")

    print("\n--- 5-fold CV: XGBoost ---")
    xgb_oof, xgb_metrics = cv_eval(X, y, lambda: xgb.XGBRegressor(
        n_estimators=500, learning_rate=0.05, max_depth=6,
        subsample=0.8, colsample_bytree=0.8, reg_alpha=0.1, reg_lambda=1.0,
        tree_method="hist", n_jobs=-1, random_state=42))
    xgb_overall = {
        "spearman": float(spearmanr(y, xgb_oof).correlation),
        "pearson": float(pearsonr(y, xgb_oof)[0]),
        "mae_pkd": float(np.mean(np.abs(y - xgb_oof))),
        "rmse_pkd": float(np.sqrt(np.mean((y - xgb_oof) ** 2))),
    }
    print(f"  OOF Spearman ρ: {xgb_overall['spearman']:.3f}")
    print(f"  OOF Pearson r:  {xgb_overall['pearson']:.3f}")
    print(f"  OOF MAE pKd:    {xgb_overall['mae_pkd']:.3f}")
    print(f"  OOF RMSE pKd:   {xgb_overall['rmse_pkd']:.3f}")

    # Ensemble
    ens_oof = 0.5 * (rf_oof + xgb_oof)
    ens_overall = {
        "spearman": float(spearmanr(y, ens_oof).correlation),
        "pearson": float(pearsonr(y, ens_oof)[0]),
        "mae_pkd": float(np.mean(np.abs(y - ens_oof))),
        "rmse_pkd": float(np.sqrt(np.mean((y - ens_oof) ** 2))),
    }
    print(f"\n--- Ensemble (RF + XGB mean) ---")
    print(f"  OOF Spearman ρ: {ens_overall['spearman']:.3f}")
    print(f"  OOF Pearson r:  {ens_overall['pearson']:.3f}")
    print(f"  OOF MAE pKd:    {ens_overall['mae_pkd']:.3f}")
    print(f"  OOF RMSE pKd:   {ens_overall['rmse_pkd']:.3f}")

    # Train final models on full data
    print("\nTraining final models on all data ...")
    rf_final = RandomForestRegressor(
        n_estimators=400, max_depth=None, min_samples_leaf=2,
        n_jobs=-1, random_state=42)
    rf_final.fit(X, y)
    xgb_final = xgb.XGBRegressor(
        n_estimators=500, learning_rate=0.05, max_depth=6,
        subsample=0.8, colsample_bytree=0.8, reg_alpha=0.1, reg_lambda=1.0,
        tree_method="hist", n_jobs=-1, random_state=42)
    xgb_final.fit(X, y)

    # Save
    with open(QSAR / "model_rf.pkl", "wb") as f: pickle.dump(rf_final, f)
    xgb_final.save_model(str(QSAR / "model_xgb.json"))
    with open(QSAR / "feature_config.json", "w") as f:
        json.dump({"radius": 2, "n_bits": 2048, "n_descs": 8}, f)
    print(f"  Saved: model_rf.pkl, model_xgb.json")

    # Save CV results
    with open(QSAR / "cv_results.json", "w") as f:
        json.dump({
            "n_total": len(rows), "n_valid": int(valid.sum()),
            "training_pkd_range": [MIN_PKD_FOR_TRAINING, MAX_PKD_FOR_TRAINING],
            "rf_overall": rf_overall, "xgb_overall": xgb_overall, "ens_overall": ens_overall,
            "rf_folds": rf_metrics, "xgb_folds": xgb_metrics,
        }, f, indent=2)

    # Save OOF predictions
    with open(QSAR / "oof_predictions.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "smiles", "pkd_true", "pkd_pred_rf", "pkd_pred_xgb", "pkd_pred_ens"])
        for i, (cid, smi) in enumerate(zip(ids, smiles)):
            if not valid[i] if i < len(valid) else True: continue
            w.writerow([cid, smi, f"{y[i]:.3f}", f"{rf_oof[i]:.3f}", f"{xgb_oof[i]:.3f}", f"{ens_oof[i]:.3f}"])

    # Predict on validation set: 3 CF Labs hits + 3 TEP fragments
    print("\n--- Predicting on reference set ---")
    ref_smis = {
        "Z795991852 (CF: 10 µM, pKd 5.00)":  "Cn1c(=O)c2ccccc2n2c(COC(=O)c3cccc(NC(=O)C4Cc5ccccc5O4)c3)nnc12",
        "Z979336988 (CF: 30 µM, pKd 4.52)":  "Cc1ccc2[nH]c(C3CCCN(C(=O)c4cccc(CN5C(=O)c6ccccc6C5=O)c4)C3)nc2c1",
        "D203-0031  (CF: 17 µM, pKd 4.77)":  "O=C(c1ccc2c(c1)OCO2)N1CCC(c2nc(O)c3nnn(Cc4ccc(F)cc4)c3n2)CC1",
        "FM001580 (TEP frag site F)":        "O=C(O)c1ccccc1OC(F)(F)F",
        "FM001452 (TEP frag site F)":        "Nc1cccc(OCc2ccccc2)c1",
        "FM002150 (TEP frag site F)":        "O=C(O)COCc1ccccc1",
    }
    Xr, _ = featurize(list(ref_smis.values()))
    rf_pred = rf_final.predict(Xr)
    xgb_pred = xgb_final.predict(Xr)
    ens_pred = 0.5 * (rf_pred + xgb_pred)
    print(f"  {'Compound':45s}  {'RF pKd':>7}  {'XGB pKd':>7}  {'ENS pKd':>7}  {'Pred Kd µM':>10}")
    for label, rfp, xgbp, ensp in zip(ref_smis.keys(), rf_pred, xgb_pred, ens_pred):
        kd_uM = 10 ** (-ensp) * 1e6
        print(f"  {label:45s}  {rfp:>7.2f}  {xgbp:>7.2f}  {ensp:>7.2f}  {kd_uM:>10.2f}")

    # Predict on the full analog pool
    print("\n--- Predicting on 503-compound analog pool ---")
    pool = list(csv.DictReader(open(DATA / "analogs" / "all_candidates.csv")))
    pool_smis = [r["smiles"] for r in pool]
    Xp, vp = featurize(pool_smis)
    pred_rf = rf_final.predict(Xp)
    pred_xgb = xgb_final.predict(Xp)
    pred_ens = 0.5 * (pred_rf + pred_xgb)

    out_path = QSAR / "predictions_analogs.csv"
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "parent_id", "smiles", "qsar_pkd_rf", "qsar_pkd_xgb",
                    "qsar_pkd_ens", "qsar_kd_uM_ens", "tanimoto_to_parent", "max_tanimoto_to_naar"])
        for i, r in enumerate(pool):
            kd_uM = 10 ** (-pred_ens[i]) * 1e6
            w.writerow([r["id"], r["parent_id"], r["smiles"],
                        f"{pred_rf[i]:.3f}", f"{pred_xgb[i]:.3f}",
                        f"{pred_ens[i]:.3f}", f"{kd_uM:.2f}",
                        r["tanimoto_to_parent"], r["max_tanimoto_to_naar"]])
    print(f"  Wrote {out_path}")

    # Top-10 by ensemble pKd
    print(f"\n  Top 10 analogs by QSAR ensemble pKd:")
    sorted_idx = np.argsort(-pred_ens)
    for i in sorted_idx[:10]:
        kd_uM = 10 ** (-pred_ens[i]) * 1e6
        print(f"    {pool[i]['id']:30s}  pKd={pred_ens[i]:.2f}  ({kd_uM:>7.2f} µM)  parent={pool[i]['parent_id']}")


if __name__ == "__main__":
    main()
