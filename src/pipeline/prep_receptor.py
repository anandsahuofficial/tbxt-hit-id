"""
Prepare TBXT G177D receptor for Vina docking.

Steps:
  1. Read 6F59 (G177D + DNA, dimer); take chain A only (one protein protomer)
  2. Strip DNA, ligands, ions, waters
  3. Add hydrogens at pH 7.5 via PDBFixer
  4. Save clean PDB + Vina-ready PDBQT (via meeko's prepare_receptor / openbabel)
  5. Print residue range, atom count, and residue 88/177 verification

Output: data/dock/receptor/6F59_apo.pdb, 6F59_apo.pdbqt, prep_log.txt
"""
import sys
from pathlib import Path

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, Select, PDBIO

DOCK = Path(__file__).resolve().parents[1] / "data/dock"
RCPT = DOCK / "receptor"
RCPT.mkdir(parents=True, exist_ok=True)

INPUT_PDB = RCPT / "6F59.pdb"


class ProteinChainASelector(Select):
    """Keep only chain A protein residues (no DNA, no waters, no ligands)."""

    def accept_residue(self, residue):
        # Standard amino acids
        STD_AA = {
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"
        }
        return residue.parent.id == "A" and residue.get_resname() in STD_AA


def main():
    log = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("6F59", str(INPUT_PDB))

    # Inspect chains
    log.append("=" * 70)
    log.append("Input 6F59 structure summary")
    log.append("=" * 70)
    for model in structure:
        for chain in model:
            residues = list(chain)
            res_types = {}
            for r in residues:
                res_types[r.get_resname()] = res_types.get(r.get_resname(), 0) + 1
            n_protein = sum(v for k, v in res_types.items() if k in {
                "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"})
            log.append(f"  chain {chain.id}: {len(residues)} residues, {n_protein} protein residues")
            unique = ", ".join(f"{k}({v})" for k, v in sorted(res_types.items())[:10])
            log.append(f"    residue types (first 10): {unique}{' ...' if len(res_types) > 10 else ''}")

    # Save chain A protein only
    io = PDBIO()
    io.set_structure(structure)
    chain_a_path = RCPT / "6F59_chainA_protein.pdb"
    io.save(str(chain_a_path), ProteinChainASelector())

    # Re-read to verify
    chain_a = parser.get_structure("6F59A", str(chain_a_path))
    chain_a_residues = list(chain_a[0]["A"].get_residues())
    log.append("")
    log.append("=" * 70)
    log.append("Chain A protein-only summary")
    log.append("=" * 70)
    log.append(f"  Residues: {len(chain_a_residues)}")
    log.append(f"  First residue: {chain_a_residues[0].get_resname()} {chain_a_residues[0].id[1]}")
    log.append(f"  Last  residue: {chain_a_residues[-1].get_resname()} {chain_a_residues[-1].id[1]}")

    # Verify residue 177 (should be ASP for G177D), residue 88 (TYR), residue 42 (LEU)
    for resnum, expected in [(42, "LEU"), (88, "TYR"), (177, "ASP")]:
        for r in chain_a_residues:
            if r.id[1] == resnum:
                log.append(f"  Residue {resnum}: {r.get_resname()} (expected {expected}, "
                           f"{'MATCH' if r.get_resname() == expected else 'MISMATCH'})")
                break
        else:
            log.append(f"  Residue {resnum}: NOT FOUND in chain A")

    # PDBFixer: protonation at pH 7.5
    log.append("")
    log.append("=" * 70)
    log.append("PDBFixer: clean + protonate at pH 7.5")
    log.append("=" * 70)
    fixer = PDBFixer(filename=str(chain_a_path))
    fixer.findMissingResidues()
    log.append(f"  Missing residues: {len(fixer.missingResidues)}")
    fixer.findNonstandardResidues()
    log.append(f"  Nonstandard residues: {len(fixer.nonstandardResidues)}")
    if fixer.nonstandardResidues:
        fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    log.append(f"  Missing atoms in residues: {len(fixer.missingAtoms)}")
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.5)

    apo_pdb_path = RCPT / "6F59_apo.pdb"
    with open(apo_pdb_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    log.append(f"  Wrote: {apo_pdb_path}")

    # Convert to PDBQT for Vina via openbabel
    pdbqt_path = RCPT / "6F59_apo.pdbqt"
    import subprocess
    log.append("")
    log.append("=" * 70)
    log.append("OpenBabel: convert to PDBQT")
    log.append("=" * 70)
    rc = subprocess.run(
        ["obabel", str(apo_pdb_path), "-O", str(pdbqt_path), "-xr", "-p", "7.5"],
        capture_output=True, text=True
    )
    log.append(f"  Return code: {rc.returncode}")
    if rc.stdout: log.append(f"  stdout: {rc.stdout.strip()}")
    if rc.stderr: log.append(f"  stderr: {rc.stderr.strip()[:300]}")
    if pdbqt_path.exists():
        log.append(f"  Wrote: {pdbqt_path} ({pdbqt_path.stat().st_size} bytes)")

    log_text = "\n".join(log)
    print(log_text)
    with open(RCPT / "prep_log.txt", "w") as f:
        f.write(log_text)


if __name__ == "__main__":
    main()
