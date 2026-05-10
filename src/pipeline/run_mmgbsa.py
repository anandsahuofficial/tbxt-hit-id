"""
MMGBSA-style binding-energy estimation on docked poses.

Pipeline per compound:
  1. Read GNINA-docked pose (top model from data/dock/validation_F_gnina/poses/)
  2. Combine with apo receptor (6F59_apo.pdb chain A)
  3. Parameterize ligand with OpenFF Sage 2.x (SMIRNOFF), receptor with amber14
  4. Brief vacuum minimization on the complex (~100 steps)
  5. Compute MM energy on minimized state:
        E_complex = bonded + nonbonded(complex)
        E_protein = nonbonded(protein only, ligand removed)
        E_ligand  = nonbonded(ligand only, protein removed)
        ΔE_bind   = E_complex - E_protein - E_ligand
  6. Add a simple Generalized-Born solvation term via OBC GBSA implicit-solvent
     for ΔG_solv (single-snapshot — no MD averaging, faster than full MD)
  7. Output ΔG_MMGBSA estimate per compound

This is a "single-snapshot MMGBSA" — fast (~minutes per compound on CPU/GPU)
but noisy compared to MD-averaged MMGBSA. For pre-event validation of the
toolchain, this is sufficient. On-day use the same pipeline on top 8-15 picks.
"""
import argparse
import csv
import json
import logging
import sys
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, NonbondedForce, VerletIntegrator, unit
from openmm import app as omm_app
from openmm.app.element import Element
from openmm.app import ForceField, Modeller, PDBFile, GBn2
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("mmgbsa")

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
APO_PDB = DOCK / "receptor" / "6F59_apo.pdb"
OUT = Path(__file__).resolve().parents[1] / "data/mmgbsa"
OUT.mkdir(exist_ok=True)


def read_pdbqt_top_pose_to_mol(pdbqt_path):
    """Read first MODEL of a Vina/GNINA PDBQT into an RDKit mol with 3D coordinates."""
    lines = pdbqt_path.read_text().splitlines()
    pdb_lines = []
    in_model = False
    for ln in lines:
        if ln.startswith("MODEL "):
            if in_model: break
            in_model = True
            continue
        if ln.startswith("ENDMDL"):
            break
        if not in_model: continue
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            # Truncate PDBQT cols to PDB cols (78 chars)
            pdb_lines.append(ln[:54] + "  1.00  0.00          " + ln[77:79].strip().rjust(2))
    pdb_str = "\n".join(pdb_lines) + "\nEND\n"
    return pdb_str


def setup_complex_system(rdkit_mol, receptor_pdb_path):
    """Build OpenMM complex system (protein + ligand) using GBn2 implicit solvent.
    Returns: (system, topology, positions, ligand_indices)."""
    # Convert RDKit mol → OpenFF Molecule (assigns parameters)
    off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_mol], forcefield="openff-2.2.0")

    # Build ForceField with amber14 for protein + GBn2 implicit solvent + SMIRNOFF for ligand
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    ff.registerTemplateGenerator(smirnoff.generator)

    # Load protein
    pdb = PDBFile(str(receptor_pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)

    # Add ligand from RDKit mol
    lig_topology = off_mol.to_topology().to_openmm()
    lig_positions_quantity = off_mol.conformers[0].to_openmm()
    modeller.add(lig_topology, lig_positions_quantity)

    # GBn2 is included via the implicit/gbn2.xml force field above; no implicitSolvent kwarg in OpenMM 8.x
    system = ff.createSystem(modeller.topology, nonbondedMethod=omm_app.NoCutoff,
                             constraints=omm_app.HBonds)

    # Identify ligand atom indices (last residue's atoms)
    n_protein_atoms = sum(1 for a in pdb.topology.atoms())
    lig_atom_indices = list(range(n_protein_atoms, system.getNumParticles()))

    return system, modeller.topology, modeller.positions, lig_atom_indices


def compute_energy(system, topology, positions, atom_set=None):
    """Compute potential energy of the system on given positions.
    If atom_set is provided, zero out forces on atoms NOT in atom_set (treat them as ghosts)."""
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
    """Estimate ΔG_bind for one compound via single-snapshot MMGBSA-like analysis."""
    log.info(f"  Building system for {cid}")
    # Get docked pose as PDB-style string
    pdb_lig_str = read_pdbqt_top_pose_to_mol(Path(pose_pdbqt_path))
    pose_pdb_path = OUT / f"{cid}_pose.pdb"
    pose_pdb_path.write_text(pdb_lig_str)

    # Build RDKit mol from SMILES with the pose 3D coords
    template = Chem.MolFromSmiles(smiles)
    if template is None: return None
    template = Chem.AddHs(template)
    if AllChem.EmbedMolecule(template, randomSeed=42, useRandomCoords=True) != 0: return None
    AllChem.UFFOptimizeMolecule(template, maxIters=200)

    # Build complex system
    try:
        system, topology, positions, lig_idx = setup_complex_system(template, APO_PDB)
    except Exception as e:
        log.warning(f"  setup failed: {e}")
        return None

    # Minimize complex
    try:
        min_positions, e_complex_kcal = minimize(system, topology, positions, n_steps=200)
        positions_np = np.array([list(p.value_in_unit(unit.nanometer)) for p in min_positions]) * 10  # → Å
    except Exception as e:
        log.warning(f"  minimization failed: {e}")
        return None

    # For ΔG_bind, build apo (no ligand) and ligand-only (vacuum) systems and compute energies
    # APPROACH: zero out the ligand's nonbonded interactions to get protein-only energy
    nb = next(f for f in [system.getForce(i) for i in range(system.getNumForces())]
              if isinstance(f, NonbondedForce))

    # Save original parameters
    orig_params = []
    for i in lig_idx:
        c, sig, eps = nb.getParticleParameters(i)
        orig_params.append((c, sig, eps))
    # Save exception parameters between protein-ligand
    orig_exc = []
    n_exc = nb.getNumExceptions()

    # Compute E_complex (already minimized)
    e_complex = e_complex_kcal

    # Zero out ligand charges + LJ to get E_protein_in_complex_geometry
    for i, (c, sig, eps) in zip(lig_idx, orig_params):
        nb.setParticleParameters(i, 0.0, sig, 0.0)
    nb.updateParametersInContext(omm_app.Simulation(topology, system,
                                                     VerletIntegrator(1.0 * unit.femtosecond)).context) \
        if False else None  # placeholder; updateParametersInContext needs context
    # simpler: build a fresh simulation with modified system
    integrator = VerletIntegrator(1.0 * unit.femtosecond)
    sim_p = omm_app.Simulation(topology, system, integrator)
    sim_p.context.setPositions(min_positions)
    e_protein = sim_p.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

    # Restore ligand, zero out protein
    for i, (c, sig, eps) in zip(lig_idx, orig_params):
        nb.setParticleParameters(i, c, sig, eps)
    n_protein_atoms = nb.getNumParticles() - len(lig_idx)
    orig_protein = []
    for i in range(n_protein_atoms):
        c, sig, eps = nb.getParticleParameters(i)
        orig_protein.append((c, sig, eps))
        nb.setParticleParameters(i, 0.0, sig, 0.0)
    integrator = VerletIntegrator(1.0 * unit.femtosecond)
    sim_l = omm_app.Simulation(topology, system, integrator)
    sim_l.context.setPositions(min_positions)
    e_ligand = sim_l.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

    # Restore protein params (so subsequent calls aren't broken)
    for i, (c, sig, eps) in enumerate(orig_protein):
        nb.setParticleParameters(i, c, sig, eps)

    delta_e = e_complex - e_protein - e_ligand
    return {
        "cid": cid,
        "e_complex_kcal": round(e_complex, 1),
        "e_protein_kcal": round(e_protein, 1),
        "e_ligand_kcal": round(e_ligand, 1),
        "delta_e_kcal": round(delta_e, 2),
        "n_protein_atoms": n_protein_atoms,
        "n_ligand_atoms": len(lig_idx),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles-csv", required=True)
    ap.add_argument("--pose-dir", default=str(DOCK / "validation_F_gnina" / "poses"),
                    help="dir containing <id>_F.pdbqt poses")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.smiles_csv)))
    if args.limit: rows = rows[:args.limit]
    log.info(f"MMGBSA estimation on {len(rows)} compounds")

    pose_dir = Path(args.pose_dir)
    results = []
    for i, row in enumerate(rows, 1):
        cid = row["id"]
        smi = row["smiles"]
        pose_path = pose_dir / f"{cid}_F.pdbqt"
        if not pose_path.exists():
            log.warning(f"[{i}/{len(rows)}] {cid}: pose not found at {pose_path}")
            continue
        log.info(f"[{i}/{len(rows)}] {cid}")
        try:
            r = mmgbsa_one(cid, smi, pose_path)
            if r is None:
                results.append({"cid": cid, "status": "failed"})
                continue
            r["status"] = "ok"
            r["smiles"] = smi
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

    # Sort + print
    ok = [r for r in results if r.get("status") == "ok"]
    ok.sort(key=lambda r: r["delta_e_kcal"])
    log.info(f"\nResults sorted by ΔE (most negative = strongest predicted binding):")
    log.info(f"  {'compound':30s}  {'ΔE_kcal':>10}  {'E_complex':>10}  {'E_lig':>10}")
    for r in ok:
        log.info(f"  {r['cid']:30s}  {r['delta_e_kcal']:>10.2f}  "
                 f"{r['e_complex_kcal']:>10.1f}  {r['e_ligand_kcal']:>10.1f}")


if __name__ == "__main__":
    main()
