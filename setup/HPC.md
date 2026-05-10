# HPC + GPU Notes

This pipeline is designed to run on a single workstation with one
modern NVIDIA GPU (RTX 4090 or A100), but several stages benefit
from HPC. This file collects the operational notes that took us a
weekend to figure out.

## Single-workstation reference

| Stage | Time on RTX 4090 (24 GB) |
|---|---|
| Vina ensemble (570 cmpds × 6 receptor confs) | ~45 min |
| GNINA CNN re-rank on Vina poses | ~30 min |
| Boltz-2 co-folding (570 cmpds, 1 trajectory each) | ~3 h |
| MMGBSA on top 30 (100 ps refinement each) | ~1.5 h |
| QSAR train + predict | ~2 min |
| Paralog selectivity scan (16 paralogs) | ~5 min |
| **End-to-end** | **~6 h** |

For the full HPC variant matrix (5 overnight variants - receptor
ensemble × site, dual-engine Boltz, MMGBSA + alchemical FEP, site-G
re-dock), allocate ~12 h on a multi-GPU node.

## GPU memory

- **Boltz-2** is the memory-hungriest step. ≥ 16 GB GPU RAM required;
  24 GB is comfortable for 570 compounds in single-batch mode.
- **GNINA** runs comfortably in 8 GB but benefits from 16+ GB for
  the multi-seed (10-restart) pose-stability variant.
- **Vina** is CPU-only and parallelizes by receptor conformation -
  one CPU per conformation × 6 conformations = 6 cores recommended.

## GLIBC: GNINA on older HPC clusters

GNINA 1.x ships dynamically linked against `glibc ≥ 2.29`. Several
university HPC clusters still run RHEL/CentOS 7 with `glibc 2.17`
or 2.28 (Boston University's SCC was 2.28 at the time of the
hackathon). Symptom:

```
gnina: /lib64/libc.so.6: version `GLIBC_2.29' not found
```

**Fix:** run GNINA inside a Singularity / Apptainer container based
on Ubuntu 22.04 (which carries `glibc 2.35`). A working recipe:

```dockerfile
# Singularity definition (cuda12-cudnn9-ubuntu22.def)
Bootstrap: docker
From: nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04

%post
    apt-get update && apt-get install -y \
        wget git python3 python3-pip libboost-all-dev \
        && rm -rf /var/lib/apt/lists/*
    cd /opt && wget -q https://github.com/gnina/gnina/releases/download/v1.0/gnina \
        && chmod +x gnina

%environment
    export PATH="/opt:$PATH"
```

Build + bind the project / cluster scratch dirs:

```bash
sudo singularity build cuda12-cudnn9-ubuntu22.sif cuda12-cudnn9-ubuntu22.def
singularity exec --nv \
    --bind /projectnb,/scratch \
    cuda12-cudnn9-ubuntu22.sif \
    gnina --receptor data/receptor/6F59_chainA.pdbqt \
          --ligand   data/ligands/cmpd001.pdbqt \
          --out      out/cmpd001_gnina.pdbqt \
          --autobox_ligand data/receptor/site_F.pdb
```

The pipeline reads `GNINA_BIN` from the environment, so wrapper
scripts can point `GNINA_BIN=/opt/wrappers/gnina_singularity.sh`
to make the pipeline platform-agnostic.

## Boltz-2 model cache on disk-quota'd HPC homes

University HPC clusters often impose a 5–20 GB quota on `$HOME`,
which is much smaller than Boltz's ~7.6 GB model cache. First-run
Boltz will fail with `OSError: No space left on device` halfway
through download.

**Fix:** rsync the model cache to project space and symlink:

```bash
mkdir -p /projectnb/<your-project>/$USER/.boltz
rsync -av $HOME/.boltz/ /projectnb/<your-project>/$USER/.boltz/
rm -rf $HOME/.boltz
ln -s /projectnb/<your-project>/$USER/.boltz $HOME/.boltz
```

Or set the cache root via env: `export BOLTZ_CACHE_DIR=/path/to/space`.

## Resumable per-compound output

If your scoring run gets killed (job timeout, node reboot), you
don't want to start over. The pipeline writes per-compound CSV rows
incrementally and skips compounds whose row is already in the
output file. Concretely, pass `--resume` to any of:

```bash
python -m src.pipeline.dock_vina  --pool data/pool/candidate_pool_570.csv --resume
python -m src.pipeline.dock_gnina --pool data/pool/candidate_pool_570.csv --resume
python -m src.pipeline.run_boltz  --pool data/pool/candidate_pool_570.csv --resume
```

(The `src/` layer is pending; this README will be updated when those
modules land in the public repo.)

## SLURM submission templates

Per-stage SLURM scripts (single-GPU, multi-GPU array, CPU-only) are
not bundled with this repo to keep it cluster-neutral. The shape
of a typical Boltz run is:

```bash
#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --mem=32G

source ~/miniconda3/bin/activate tbxt-hit-id
export BOLTZ_CACHE_DIR=/projectnb/<proj>/$USER/.boltz
python -m src.pipeline.run_boltz --pool data/pool/candidate_pool_570.csv --resume
```

## Known gotchas

- **`MMGBSA` + AMBER force fields** require both `openmm` and
  `openmmforcefields`. The conda env spec installs both; if you
  only get OpenMM, MMGBSA will silently fall back to ff14SB-only and
  fail on non-standard residues.
- **`OPENBABEL` PDBQT conversion** writes per-atom partial charges
  using Gasteiger by default. For receptors with explicit metals or
  unusual residues, this can produce charge artifacts at the metal
  site. Inspect the PDBQT for `--` (negative-negative) pairs near
  metals.
- **First-time conda env solves** can take 30+ minutes if `bioconda`
  channels are slow. Use `mamba` (`pip install mamba` then
  `mamba env create -f environment.yml`) for ~5× faster solves.
