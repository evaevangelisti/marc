#!/bin/bash
#SBATCH --job-name=marc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -e

module load anaconda3/2023.09-0
eval "$(conda shell.bash hook)"
conda activate "$WORK/envs/pymol_env"

module load profile/lifesc
module load gromacs/2021.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0-cuda-11.5.0

INPUT_DIR="./data/raw"
PROCESSED_DIR="./data/processed"
RESULTS_DIR="./results"

main() {
  rm -rf "$PROCESSED_DIR" "$RESULTS_DIR" logs
  mkdir -p "$PROCESSED_DIR" "$RESULTS_DIR" logs

  for pdb in "$INPUT_DIR"/*.pdb; do
    ./scripts/clear_pdb.sh -p "$pdb" -o "$PROCESSED_DIR" || {
      echo "$pdb: cleaning failed" >&2
      continue
    }

    ./scripts/simulate.sh -p "$PROCESSED_DIR/$(basename "$pdb" .pdb)_clean.pdb" -o "$RESULTS_DIR/$(basename "$pdb" .pdb)" || {
      echo "$pdb: simulation failed" >&2
      continue
    }
  done
}

main "$@"
conda deactivate
