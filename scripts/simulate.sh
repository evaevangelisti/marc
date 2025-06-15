#!/bin/bash

set -e

pdb=""
co2_molecules=1
output_dir=""

usage() {
  echo "usage: $0 [-h] [-p <pdb>] [-n <co2-molecules>] [-o <output-dir>]"
}

main() {
  while getopts ":hp:n:o:" opt; do
    case $opt in
      h) usage; exit 0 ;;
      p) pdb="$OPTARG" ;;
      n) co2_molecules="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '$pdb' must be a valid .pdb file" >&2; usage; exit 1
  fi

  if [[ -z "$output_dir" ]]; then
    output_dir="./$(basename "$pdb" .pdb)"
  fi

  mkdir -p "$output_dir"

  gmx pdb2gmx -f "$pdb" -o "$output_dir/processed.gro" -p "$output_dir/topol.top" -ff amber99sb-ildn -water tip3p
  gmx editconf -f "$output_dir/processed.gro" -o "$output_dir/box.gro" -c -d 1.0
  gmx solvate -cp "$output_dir/box.gro" -cs spc216.gro -o "$output_dir/solvated.gro" -p "$output_dir/topol.top"
  gmx genion -s "$output_dir/solvated.tpr" -o "$output_dir/ions.gro" -p "$output_dir/topol.top" -pname NA -nname CL -neutral

  TMP_DIR=$(mktemp -d)

  co2_pdb="$TMP_DIR/co2.pdb"
  cat <<EOF > "$co2_pdb"
HEADER    CARBON DIOXIDE
HETATM    1  C   CO2     0       0.000   0.000   0.000  0.00  0.00           C
HETATM    2  O   CO2     0       0.000   0.000  -1.208  0.00  0.00           O
HETATM    3  O   CO2     0       0.000   0.000   1.208  0.00  0.00           O
CONECT    1    2    3
CONECT    2    1
CONECT    3    1
TER
END
EOF

  gmx insert-molecules -f "$output_dir/ions.gro" -ci "$co2_pdb" -o "$output_dir/co2.gro" -p "$output_dir/topol.top" -nmol "$co2_molecules"

  rm -rf "$TMP_DIR"

  gmx grompp -f minimization.mdp -c "$output_dir/co2.gro" -p "$output_dir/topol.top" -o "$output_dir/minimization.tpr"
  gmx mdrun -v -deffnm "$output_dir/minimization"

  gmx grompp -f nvt.mdp -c "$output_dir/minimization.gro" -p "$output_dir/topol.top" -o "$output_dir/nvt.tpr"
  gmx mdrun -deffnm "$output_dir/nvt.tpr"

  gmx grompp -f npt.mdp -c "$output_dir/nvt.gro" -p "$output_dir/topol.top" -o "$output_dir/npt.tpr"
  gmx mdrun -deffnm "$output_dir/npt.tpr"

  gmx grompp -f production.mdp -c "$output_dir/npt.gro" -p "$output_dir/topol.top" -o "$output_dir/production.tpr"
  gmx mdrun -deffnm "$output_dir/production"

  gmx rms -s "$output_dir/production.tpr" -f "$output_dir/production.xtc" -o "$output_dir/rmsd.xvg"
  gmx rmsf -s "$output_dir/production.tpr" -f "$output_dir/production.xtc" -o "$output_dir/rmsf.xvg"
  gmx distance -s "$output_dir/production.tpr" -f "$output_dir/production.xtc" -select 'com of group "CO2" plus com of group "Rubisco"' -o "$output_dir/distances.xvg"
}

main "$@"
