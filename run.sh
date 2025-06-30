#!/bin/bash
#SBATCH --job-name=marc
#SBATCH --partition=g100_usr_prod
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=48
#SBATCH --mem=350G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --array=0-5

set -e

module load anaconda3/2023.09-0
eval "$(conda shell.bash hook)"
conda activate "$WORK/envs/pymol_env"
trap 'conda deactivate' EXIT

module load profile/lifesc
module load gromacs/2021.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0-cuda-11.5.0

INPUT_DIR="./data/raw"
PROCESSED_DIR="./data/processed"
RESULTS_DIR="./results"

PDB_LIST="$INPUT_DIR/*.pdb"

main() {
  pdb="${PDB_LIST[$SLURM_ARRAY_TASK_ID]}"

  ./scripts/clear_pdb.sh -i "$pdb" -o "$PROCESSED_DIR" || {
    echo "$pdb: cleaning failed" >&2
    exit 1
  }

  ./scripts/simulate.sh -i "$PROCESSED_DIR/$(basename "$pdb" .pdb)_clean.pdb" -o "$RESULTS_DIR/$(basename "$pdb" .pdb)" -p || {
    echo "$pdb: simulation failed" >&2
    exit 1
  }
}

main "$@"
