# Setup

Cold-start setup for `tbxt-hit-id`. Two paths:

1. **Container (recommended)** — pull a single CUDA-enabled image
   from GHCR, get every binary + every Python dep + GNINA in one go.
2. **Native conda** — `conda env create -f environment.yml` and
   install GNINA / boltz weights yourself.

The container is the supported path; the conda path is for users
who prefer to live-edit deps or run on systems without
Apptainer/Docker.

---

## Path A — Container (recommended)

```bash
# 1. Pull the container (~6-12 GB; first pull takes 5-15 min)
bash setup/pull_container.sh
# → ./tbxt-hit-id.sif

# 2. Fetch the receptor (PDB 6F59:A)
bash setup/fetch_receptor.sh

# 3. Fetch the bulk data assets (compound pool, Naar SPR Kd, receptor ensemble)
bash setup/fetch_data.sh
# OPTIONAL: also pull pre-computed Boltz/GNINA/MMGBSA outputs (~600 MB)
bash setup/fetch_data.sh --include-poses

# 4. Run the pipeline
apptainer exec --nv --bind $PWD tbxt-hit-id.sif \
    bash examples/reproduce_top4.sh

# Or, demo mode (no GPU; needs --include-poses):
apptainer exec --bind $PWD tbxt-hit-id.sif \
    bash examples/reproduce_top4.sh --demo
```

The container image is built from [`Containerfile`](Containerfile)
and pushed to `ghcr.io/anandsahuofficial/tbxt-hit-id:latest` by the
[`.github/workflows/container.yml`](../.github/workflows/container.yml)
workflow on every push to `main`. To pull a specific tag:

```bash
bash setup/pull_container.sh --tag v1.0.0
```

To override the registry (mirror, fork, private registry):

```bash
GHCR_OWNER=<other> GHCR_REPO=<repo> bash setup/pull_container.sh
# or fully:
bash setup/pull_container.sh --image-ref docker://my.registry/foo:v1
```

### What's in the container

| Layer | Contents |
|---|---|
| Base | `nvidia/cuda:12.8.0-cudnn-runtime-ubuntu22.04` (Ubuntu 22.04 = glibc 2.35 — solves the GNINA glibc problem) |
| Python | Miniconda + `mamba` |
| Conda env (`tbxt-hit-id`) | Everything in [`environment.yml`](../environment.yml): RDKit, AutoDock Vina, OpenMM, OpenFF, scikit-learn, XGBoost, Biopython, pytorch (CUDA 12.8 wheel), boltz, etc. |
| Binaries | GNINA 1.0 (under `/opt/gnina`) |
| Source | The repo's `src/`, `tools/`, `examples/`, `setup/`, `docs/` baked under `/opt/tbxt-hit-id` |

The image is **single-arch (linux/amd64)** — that's the only
architecture CUDA + GNINA support.

### What's still on HuggingFace

The container ships **everything except the data**. Data lives in
the companion HF dataset (default: `anandsahuofficial/tbxt-hit-id-data`):

- `pool/candidate_pool_570.csv` — 570 SMILES
- `naar/naar_spr_kd_650.csv` — 650 measured affinities
- `receptor/6F59_chainA_ensemble.tar.gz` — 6 relaxed receptor conformations
- (optional) `scored/{boltz,gnina,mmgbsa}_outputs.tar.gz` — pre-computed signal CSVs

This split is intentional: containers are versioned with the code;
data is versioned independently and updates more often.

---

## Path B — Native conda

```bash
# 1. Create the env
conda env create -f environment.yml          # or: mamba env create -f ...
conda activate tbxt-hit-id

# 2. Install GNINA (Linux only — see HPC.md for the Singularity workaround
#    on clusters with glibc < 2.29)
wget -qO ~/bin/gnina https://github.com/gnina/gnina/releases/download/v1.0/gnina
chmod +x ~/bin/gnina

# 3. Steps 2-4 from Path A (fetch_receptor.sh, fetch_data.sh, reproduce_top4.sh)
```

The conda env defaults to `pytorch-cpu`; if a GPU is present,
[`fetch_data.sh`](fetch_data.sh) auto-upgrades to the CUDA 12.8 wheel
on first run. To force this earlier:

```bash
pip install --force-reinstall --no-deps \
    "torch==2.8.0" "torchvision==0.23.0" \
    --index-url https://download.pytorch.org/whl/cu128
```

---

## What gets installed where

```
./tbxt-hit-id.sif                 ← container (Path A only, ~6-12 GB)
~/miniconda3/envs/tbxt-hit-id/    ← conda env (Path B only, ~3 GB)
<repo>/data/receptor/             ← 6F59:A PDB + PDBQT  (~200 KB)
<repo>/data/pool/                 ← 570-compound SMILES CSV  (~80 KB)
<repo>/data/naar/                 ← 650 measured SPR Kd values  (~40 KB)
<repo>/data/receptor/ensemble/    ← 6 relaxed receptor confs  (~3 MB)
<repo>/data/scored/  (optional)   ← pre-computed Boltz/GNINA/MMGBSA  (~600 MB)
```

`data/` is gitignored — bulk assets stay on HF, not in the repo.

---

## Configuring the data source

[`fetch_data.sh`](fetch_data.sh) defaults to
`anandsahuofficial/tbxt-hit-id-data` on HuggingFace. Override:

```bash
HF_USER=<user> HF_REPO=<repo> bash setup/fetch_data.sh
HF_TOKEN=hf_xxx bash setup/fetch_data.sh    # for private datasets
```

The bundle SHA-256 manifest is at `CHECKSUMS.sha256` in the dataset
repo; the script verifies every download against it.

---

## HPC

For HPC-specific notes — Singularity container for GNINA on
clusters with `glibc < 2.29`, Boltz model cache redirection on
disk-quota'd home directories, SLURM templates, resumable per-
compound output — see [`HPC.md`](HPC.md). Note that **Path A
(container) already addresses the GNINA glibc issue**, so most HPC
users will want Path A.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `apptainer pull` returns HTTP 401 | Default expects public image; if you forked or moved the repo, the GHCR image visibility may not be set to public | Set the GHCR package visibility to public, or pass `--image-ref` to a public mirror |
| `apptainer exec --nv` fails: "no GPU detected" | Apptainer's `--nv` flag needs the host nvidia driver to match the container's CUDA runtime | Verify `nvidia-smi` on the host; check that the host driver supports CUDA 12.x |
| `gnina: GLIBC_2.29 not found` | Path B (native) on an old HPC | Switch to Path A (container) — that's literally what it solves |
| `OSError: No space left on device` during Boltz first run | $HOME quota too small for ~7.6 GB Boltz cache | See [`HPC.md`](HPC.md) cache redirection |
| `conda env create` hangs > 30 min | Slow channel solve in Path B | `mamba env create -f environment.yml` |
| `torch.cuda.is_available() == False` on a GPU box (Path B) | CPU torch installed by default | See "Path B" above |
| HF download returns HTML instead of binary | Wrong repo / private repo without `HF_TOKEN` / rate-limited | Check the URL printed by the script; set `HF_TOKEN` if private |
