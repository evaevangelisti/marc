#!/bin/bash
#SBATCH --job-name marc
#SBATCH -N1 --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=ELIX6_castrign
#SBATCH --partition=g100_usr_prod

set -e

module load anaconda3/2023.09-0
eval "$(conda shell.bash hook)"
conda activate "$WORK/envs/pymol_env"
trap 'conda deactivate' EXIT

module load profile/lifesc
module load autoload gromacs/2021.2

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

PDB="./data/raw/1gk8.pdb"
PROCESSED_DIR="./data/processed"
RESULTS_DIR="./results"

mkdir -p "$PROCESSED_DIR" "$RESULTS_DIR"

./scripts/clear_pdb.sh -i "$PDB" -o "$PROCESSED_DIR" || {
  echo "$PDB: cleaning failed" >&2
  exit 1
}

./scripts/simulate.sh -i "$PROCESSED_DIR/$(basename "$PDB" .pdb)_clean.pdb" -o "$RESULTS_DIR/$(basename "$PDB" .pdb)" || {
  echo "$PDB: simulation failed" >&2
  exit 1
}
