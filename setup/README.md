# Setup

Cold-start setup for `tbxt-hit-id`. Tested on Linux (Ubuntu 22.04,
RHEL 8/9) with NVIDIA CUDA 12.x. macOS works for everything except
GNINA (Linux-only binary).

## Quick start

```bash
# 1. Create the conda environment (one command — see ../environment.yml)
conda env create -f ../environment.yml
conda activate tbxt-hit-id

# 2. Fetch the receptor structure (PDB 6F59 chain A) from RCSB
bash setup/fetch_receptor.sh

# 3. Fetch the bulk data assets (compound pool, Naar SPR Kd, receptor
#    ensemble) from the companion HuggingFace dataset
bash setup/fetch_data.sh

# 4. (Optional) Skip expensive scoring by also pulling pre-computed
#    Boltz / GNINA / MMGBSA outputs (~600 MB)
bash setup/fetch_data.sh --include-poses

# 5. Reproduce the top 4 picks end-to-end
bash examples/reproduce_top4.sh
```

## What gets installed where

```
~/miniconda3/envs/tbxt-hit-id/   ← conda env (~3 GB)
<repo>/data/receptor/            ← 6F59:A PDB + PDBQT  (~200 KB)
<repo>/data/pool/                ← 570-compound SMILES CSV  (~80 KB)
<repo>/data/naar/                ← 650 measured SPR Kd values  (~40 KB)
<repo>/data/receptor/ensemble/   ← 6 relaxed receptor confs  (~3 MB)
<repo>/data/scored/  (optional)  ← pre-computed Boltz / GNINA / MMGBSA  (~600 MB)
```

`data/` is gitignored — bulk assets are intentionally not committed
to keep the repo small.

## Environment spec

The single source of truth for dependencies is
[`../environment.yml`](../environment.yml). The conda layer covers
everything with a clean conda recipe (RDKit, AutoDock Vina,
OpenMM, scikit-learn, etc.). The pip layer adds Boltz-2 and a few
PyPI-only packages.

For faster solves: install `mamba` first
(`conda install -n base -c conda-forge mamba`) and substitute
`mamba env create -f environment.yml` for the conda command above.

## GPU vs CPU

The conda env installs `pytorch-cpu` by default to keep the install
portable. If you have an NVIDIA GPU, [`fetch_data.sh`](fetch_data.sh)
detects `nvidia-smi` on first run and upgrades torch to the matching
CUDA 12.8 wheel automatically. To force this earlier:

```bash
pip install --force-reinstall --no-deps \
    "torch==2.8.0" "torchvision==0.23.0" \
    --index-url https://download.pytorch.org/whl/cu128
```

## Configuring the data source

[`fetch_data.sh`](fetch_data.sh) defaults to the public dataset
`anandsahuofficial/tbxt-hit-id-data` on Hugging Face. Override:

```bash
HF_USER=<user> HF_REPO=<repo> bash setup/fetch_data.sh
HF_TOKEN=hf_xxx bash setup/fetch_data.sh    # for private datasets
```

The bundle SHA-256 manifest is at `CHECKSUMS.sha256` in the dataset
repo; the script verifies every download against it.

## HPC

For HPC-specific notes — Singularity container for GNINA on
clusters with `glibc < 2.29`, Boltz model cache redirection on
disk-quota'd home directories, SLURM templates, resumable per-
compound output — see [`HPC.md`](HPC.md).

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `gnina: GLIBC_2.29 not found` | Old HPC `glibc` | See [`HPC.md`](HPC.md) Singularity recipe |
| `OSError: No space left on device` during Boltz first run | $HOME quota too small for ~7.6 GB Boltz cache | See [`HPC.md`](HPC.md) cache redirection |
| `conda env create` hangs > 30 min | Slow channel solve | `pip install mamba`, then `mamba env create -f ../environment.yml` |
| `torch.cuda.is_available() == False` on a GPU box | CPU torch installed (default) | See "GPU vs CPU" above |
| HF download returns HTML instead of binary | Wrong repo / private repo without `HF_TOKEN` / rate-limited | Check the URL printed by the script; set `HF_TOKEN` if private |
| `obabel: command not found` in `fetch_receptor.sh` | conda env not activated | `conda activate tbxt-hit-id` |
