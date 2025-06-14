#!/bin/bash
#SBATCH --job-name=marc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load anaconda3/2023.09-0
conda activate "$WORK/envs/pymol_env"

module load profile/lifesc
module load gromacs/2021.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0-cuda-11.5.0

INPUT_DIR="./data/raw/microalgae"
CO2_PDB="./data/raw/co2.pdb"

PROCESSED_DIR="./data/processed"
RESULTS_DIR="./results"

main() {
  mkdir -p "$PROCESSED_DIR" "$RESULTS_DIR" logs

  for pdb in "$INPUT_DIR"/*.pdb; do
    echo "$pdb: starting cleaning"

    cleaned_pdb=$(./scripts/clear_pdb.sh -p "$pdb" -o "$PROCESSED_DIR") || {
      echo "$pdb: cleaning failed" >&2
      continue
    }

    echo "$pdb: starting simulation"

    ./scripts/simulate.sh -p "$cleaned_pdb" -c "$CO2_PDB" -o "$RESULTS_DIR" || {
      echo "$pdb: simulation failed" >&2
      continue
    }

    echo "$pdb: simulation completed successfully"
  done

  conda deactivate
}

main "$@"
