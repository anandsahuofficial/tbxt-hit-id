"""
MMGBSA-style binding-energy estimation — FIXED version.

Builds three SEPARATE OpenMM systems (complex / apo / ligand) and computes the
potential energy of each in the minimized-complex geometry. This avoids the
ghost-atom bonded-energy bug in scripts/run_mmgbsa.py.

ΔE_bind = E_complex - E_protein - E_ligand    (kcal/mol)

Single-snapshot — no MD averaging. Suitable for ranking; ~5 min/compound.

Usage:
  python run_mmgbsa_fixed.py --smiles-csv <input.csv> --pose-dir <dir>
"""
import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, VerletIntegrator, unit
from openmm import app as omm_app
from openmm.app import ForceField, Modeller, PDBFile
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("mmgbsa_fixed")

ROOT = Path(__file__).resolve().parents[2]
DOCK = ROOT / "data/dock"
APO_PDB = DOCK / "receptor" / "6F59_apo.pdb"
OUT = ROOT / "data/mmgbsa"
OUT.mkdir(exist_ok=True)


def read_pdbqt_top_pose(pdbqt_path):
    """Read first MODEL of a Vina/GNINA PDBQT into a PDB-format string."""
    lines = pdbqt_path.read_text().splitlines()
    pdb_lines = []
    in_model = False
    for ln in lines:
        if ln.startswith("MODEL "):
            if in_model: break
            in_model = True; continue
        if ln.startswith("ENDMDL"): break
        if not in_model: continue
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            pdb_lines.append(ln[:54] + "  1.00  0.00          " + ln[77:79].strip().rjust(2))
    return "\n".join(pdb_lines) + "\nEND\n"


def build_off_mol(smiles):
    """Build OpenFF Molecule from SMILES with a 3D conformer."""
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        return None
    rdmol = Chem.AddHs(rdmol)
    if AllChem.EmbedMolecule(rdmol, randomSeed=42, useRandomCoords=True) != 0:
        return None
    AllChem.UFFOptimizeMolecule(rdmol, maxIters=200)
    return Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)


def make_force_field(off_mol):
    """Build an OpenMM ForceField with amber14 + GBn2 + SMIRNOFF for the ligand."""
    gen = SMIRNOFFTemplateGenerator(molecules=[off_mol], forcefield="openff-2.2.0")
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    ff.registerTemplateGenerator(gen.generator)
    return ff


def build_complex(off_mol, receptor_pdb_path):
    """Returns (system, topology, positions, n_protein_atoms)."""
    ff = make_force_field(off_mol)
    pdb = PDBFile(str(receptor_pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    n_protein_atoms = sum(1 for _ in pdb.topology.atoms())
    lig_topology = off_mol.to_topology().to_openmm()
    lig_positions = off_mol.conformers[0].to_openmm()
    modeller.add(lig_topology, lig_positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=omm_app.NoCutoff,
                             constraints=omm_app.HBonds)
    return system, modeller.topology, modeller.positions, n_protein_atoms


def build_protein_only(receptor_pdb_path):
    """Apo system using the same protein FF (no ligand template generator needed)."""
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    pdb = PDBFile(str(receptor_pdb_path))
    system = ff.createSystem(pdb.topology, nonbondedMethod=omm_app.NoCutoff,
                             constraints=omm_app.HBonds)
    return system, pdb.topology, pdb.positions


def build_ligand_only(off_mol):
    """Vacuum ligand system using the SMIRNOFF generator."""
    gen = SMIRNOFFTemplateGenerator(molecules=[off_mol], forcefield="openff-2.2.0")
    ff = ForceField("implicit/gbn2.xml")
    ff.registerTemplateGenerator(gen.generator)
    lig_topology = off_mol.to_topology().to_openmm()
    lig_positions = off_mol.conformers[0].to_openmm()
    system = ff.createSystem(lig_topology, nonbondedMethod=omm_app.NoCutoff,
                             constraints=omm_app.HBonds)
    return system, lig_topology, lig_positions


def energy(system, topology, positions):
    """Potential energy of a system at given positions (kcal/mol)."""
    integrator = VerletIntegrator(1.0 * unit.femtosecond)
    sim = omm_app.Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def minimize(system, topology, positions, n_steps=200):
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                          1.0 * unit.femtosecond)
    sim = omm_app.Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    sim.minimizeEnergy(maxIterations=n_steps)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    return state.getPositions(asNumpy=True), state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def mmgbsa_one(cid, smiles, pose_pdbqt_path):
    log.info(f"  Building system for {cid}")
    OUT.joinpath(f"{cid}_pose.pdb").write_text(read_pdbqt_top_pose(Path(pose_pdbqt_path)))

    off_mol = build_off_mol(smiles)
    if off_mol is None:
        return None

    try:
        system, topology, positions, n_p = build_complex(off_mol, APO_PDB)
    except Exception as e:
        log.warning(f"  complex build failed: {e}")
        return None

    try:
        min_pos, e_complex = minimize(system, topology, positions, n_steps=200)
    except Exception as e:
        log.warning(f"  minimization failed: {e}")
        return None

    n_total = system.getNumParticles()
    n_l = n_total - n_p

    # Split minimized geometry: protein atoms = first n_p; ligand atoms = remainder.
    protein_positions = [min_pos[i] for i in range(n_p)]
    ligand_positions  = [min_pos[i] for i in range(n_p, n_total)]

    # E_protein in apo system at the same protein geometry
    try:
        sys_p, top_p, _ = build_protein_only(APO_PDB)
        e_protein = energy(sys_p, top_p, protein_positions)
    except Exception as e:
        log.warning(f"  apo energy failed: {e}")
        return None

    # E_ligand in vacuum system at the same ligand geometry
    try:
        sys_l, top_l, _ = build_ligand_only(off_mol)
        e_ligand = energy(sys_l, top_l, ligand_positions)
    except Exception as e:
        log.warning(f"  ligand-only energy failed: {e}")
        return None

    delta_e = e_complex - e_protein - e_ligand
    return {
        "cid": cid,
        "e_complex_kcal": round(e_complex, 2),
        "e_protein_kcal": round(e_protein, 2),
        "e_ligand_kcal":  round(e_ligand, 2),
        "delta_e_kcal":   round(delta_e, 2),
        "n_protein_atoms": n_p,
        "n_ligand_atoms":  n_l,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles-csv", required=True)
    ap.add_argument("--pose-dir", default=str(ROOT / "data/full_pool_gnina_F" / "poses"),
                    help="dir containing <id>_F.pdbqt poses")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit: rows = rows[:args.limit]
    log.info(f"MMGBSA (FIXED) on {len(rows)} compounds")

    pose_dir = Path(args.pose_dir)
    results = []
    for i, row in enumerate(rows, 1):
        cid, smi = row["id"], row["smiles"]
        pose_path = pose_dir / f"{cid}_F.pdbqt"
        if not pose_path.exists():
            log.warning(f"[{i}/{len(rows)}] {cid}: pose not found at {pose_path}")
            continue
        log.info(f"[{i}/{len(rows)}] {cid}")
        try:
            r = mmgbsa_one(cid, smi, pose_path)
            if r is None:
                results.append({"cid": cid, "status": "failed"}); continue
            r["status"] = "ok"; r["smiles"] = smi
            results.append(r)
            log.info(f"  ΔE = {r['delta_e_kcal']:.2f} kcal/mol")
        except Exception as e:
            log.error(f"  error: {e}")
            import traceback; traceback.print_exc()
            results.append({"cid": cid, "status": f"error: {e}"})

    cols = ["cid", "smiles", "status", "delta_e_kcal", "e_complex_kcal",
            "e_protein_kcal", "e_ligand_kcal", "n_protein_atoms", "n_ligand_atoms"]
    out_csv = OUT / "mmgbsa_summary.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in cols})
    log.info(f"\nWrote {out_csv}")

    ok = [r for r in results if r.get("status") == "ok"]
    ok.sort(key=lambda r: r["delta_e_kcal"])
    log.info(f"\nResults sorted by ΔE (most negative = strongest predicted binding):")
    log.info(f"  {'compound':30s}  {'ΔE_kcal':>10}  {'E_complex':>10}  {'E_lig':>10}")
    for r in ok:
        log.info(f"  {r['cid']:30s}  {r['delta_e_kcal']:>10.2f}  "
                 f"{r['e_complex_kcal']:>10.2f}  {r['e_ligand_kcal']:>10.2f}")


if __name__ == "__main__":
    main()
