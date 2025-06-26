#!/bin/bash
#SBATCH --job-name=marc
#SBATCH --partition=g100_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=350G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -e

module load anaconda3/2023.09-0
eval "$(conda shell.bash hook)"
conda activate "$WORK/envs/pymol_env"
trap 'conda deactivate' EXIT

module load profile/lifesc
module load autoload gromacs/2021.2

INPUT_DIR="./data/raw"
PROCESSED_DIR="./data/processed"
RESULTS_DIR="./results"

main() {
  rm -rf "$PROCESSED_DIR" "$RESULTS_DIR"
  mkdir -p "$PROCESSED_DIR" "$RESULTS_DIR"

  for pdb in "$INPUT_DIR"/*.pdb; do
    ./scripts/clear_pdb.sh -i "$pdb" -o "$PROCESSED_DIR" || {
      echo "$pdb: cleaning failed" >&2
      continue
    }

    ./scripts/simulate.sh -i "$PROCESSED_DIR/$(basename "$pdb" .pdb)_clean.pdb" -o "$RESULTS_DIR/$(basename "$pdb" .pdb)" || {
      echo "$pdb: simulation failed" >&2
      continue
    }
  done
}

main "$@"
