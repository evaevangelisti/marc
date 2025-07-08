#!/bin/bash

main() {
  source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../config" && pwd)/config.sh"

  pushd "$(pwd)" > /dev/null
  cd "$PROCESSED_DIR/$PDB"
  trap 'popd > /dev/null' EXIT

  preparation_id=$(sbatch --parsable --job-name="$PDB: preparation" "$ROOT/scripts/prepare_simulation.sh")

  first_simulation_id=$(sbatch --parsable --job-name="$PDB: simulation (first half)" --dependency=afterany:"$preparation_id" "$ROOT/scripts/simulate.sh")
  sbatch --job-name="$PDB: simulation (second half)" --dependency=afterany:"$first_simulation_id" "$ROOT/scripts/simulate.sh"
}

main "$@"
