"""
Generate publication-quality PyMOL renders of docked poses + RDKit 2D structures.

Usage:
  python scripts/team/render_poses.py \
      --top-4 data/tier_a/tier_a_candidates.csv \
      --pose-dir data/full_pool_gnina_F/poses \
      --receptor data/dock/receptor/6F59_apo.pdb \
      --out-dir data/slides/renders \
      --highlight-residues Y88,D177,L42,R174,M181,T183 \
      --max-picks 8

Produces per pick:
  data/slides/renders/<id>_pose_3d.png   (PyMOL ray-traced complex view)
  data/slides/renders/<id>_2d.png         (RDKit 2D structure)

PyMOL must be installed (conda install -c conda-forge pymol-open-source).
"""
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def render_2d(smiles, out_png):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return False
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(1200, 900), kekulize=True)
    img.save(out_png, dpi=(300, 300))
    return True


def render_3d_pymol(receptor_pdb, ligand_pdbqt, out_png, highlight_residues):
    """Render a 3D PyMOL view of the protein-ligand complex with anchor residues highlighted."""
    res_list = "+".join(r[1:] for r in highlight_residues)  # "88+177+42+174+181+183"

    pml = f"""
    load {receptor_pdb}, receptor
    load {ligand_pdbqt}, ligand
    bg_color white
    hide everything

    show cartoon, receptor and chain A
    color grey80, receptor and chain A

    # Highlight anchor residues
    show sticks, receptor and (resi {res_list})
    color cyan, receptor and (resi {res_list}) and name C*
    color blue, receptor and resi 174 and name N*
    color red, receptor and resi 177 and name O*
    label receptor and (resi {res_list}) and name CA, "%s%s" % (resn, resi)

    # Ligand
    show sticks, ligand
    color magenta, ligand and name C*
    hide everything, hydro

    # H-bonds (within 3.5 Å between ligand polar atoms and any receptor N/O)
    distance hbonds, ligand and (name N+O*), receptor and (name N+O*) and chain A and (resi {res_list}), 3.5
    hide labels, hbonds
    color yellow, hbonds

    # Camera
    zoom ligand, 6
    orient ligand
    set ray_shadows, 0
    set ambient, 0.4
    set specular, 0.3

    ray 2400, 1800
    png {out_png}, dpi=300
    """
    with tempfile.NamedTemporaryFile(suffix=".pml", delete=False, mode="w") as f:
        f.write(pml)
        pml_path = f.name
    try:
        rc = subprocess.run(
            ["pymol", "-cq", pml_path],
            capture_output=True, text=True, timeout=180
        )
        return rc.returncode == 0 and Path(out_png).exists()
    except Exception as e:
        print(f"  PyMOL failed: {e}")
        return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--top-4", required=True, help="CSV with top-N picks")
    p.add_argument("--pose-dir", required=True, help="Dir containing <id>_F.pdbqt files")
    p.add_argument("--receptor", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--highlight-residues", default="Y88,D177,L42,R174,M181,T183",
                   help="Comma-separated like Y88,D177,...")
    p.add_argument("--max-picks", type=int, default=8)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    highlight = args.highlight_residues.split(",")

    rows = list(csv.DictReader(open(args.top_4)))[:args.max_picks]

    print(f"Rendering {len(rows)} picks → {out_dir}")
    for r in rows:
        cid = r["id"]
        smi = r["smiles"]
        pose = Path(args.pose_dir) / f"{cid}_F.pdbqt"
        if not pose.exists():
            print(f"  ⚠ {cid}: pose file not found at {pose}")
            continue

        # 2D
        out_2d = out_dir / f"{cid}_2d.png"
        ok2d = render_2d(smi, str(out_2d))
        # 3D
        out_3d = out_dir / f"{cid}_pose_3d.png"
        ok3d = render_3d_pymol(args.receptor, str(pose), str(out_3d), highlight)
        print(f"  {cid:30s}  2D={'✓' if ok2d else '✗'}  3D={'✓' if ok3d else '✗'}")

    print(f"\nDone. Images at: {out_dir}")


if __name__ == "__main__":
    main()
