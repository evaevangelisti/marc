#!/bin/bash

pdb=""
co2_pdb=""
output_dir="."

usage() {
  echo "usage: $0 [-h] [-p <pdb>] [-c <co2-pdb>] [-o <output-dir>]"
}

main() {
  while getopts ":hp:c:o:" opt; do
    case $opt in
      h) usage; exit 0 ;;
      p) pdb="$OPTARG" ;;
      c) co2_pdb="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '$pdb' must be a valid .pdb file" >&2
    exit 1
  fi

  if [[ -z "$co2_pdb" || ! -f "$co2_pdb" || "$co2_pdb" != *.pdb ]]; then
    echo "error: '$co2_pdb' must be a valid .pdb file" >&2
    exit 1
  fi

  pdb_output_dir="$output_dir"/"${pdb%.*}"
  mkdir -p "$pdb_output_dir"

  gmx pdb2gmx -f "$pdb" -o "$pdb_output_dir/processed.gro" -p "$pdb_output_dir/topol.top"
  gmx editconf -f "$pdb_output_dir/processed.gro" -o "$pdb_output_dir/box.gro" -c -d 1.0
  gmx solvate -cp "$pdb_output_dir/box.gro" -cs spc216.gro -o "$pdb_output_dir/solvated.gro" -p "$pdb_output_dir/topol.top"
  gmx genion -s "$pdb_output_dir/solvated.tpr" -o "$pdb_output_dir/ions.gro" -p "$pdb_output_dir/topol.top" -pname NA -nname CL -neutral
  gmx insert-molecules -f "$pdb_output_dir/ions.gro" -ci "$co2_pdb" -o "$pdb_output_dir/co2.gro" -p "$pdb_output_dir/topol.top" -nmol 9

  gmx grompp -f minimization.mdp -c "$pdb_output_dir/co2.gro" -p "$pdb_output_dir/topol.top" -o "$pdb_output_dir/minimization.tpr"
  gmx mdrun -v -deffnm "$pdb_output_dir/minimization"

  gmx grompp -f nvt.mdp -c "$pdb_output_dir/minimization.gro" -p "$pdb_output_dir/topol.top" -o "$pdb_output_dir/nvt.tpr"
  gmx mdrun -deffnm "$pdb_output_dir/nvt.tpr"

  gmx grompp -f npt.mdp -c "$pdb_output_dir/nvt.gro" -p "$pdb_output_dir/topol.top" -o "$pdb_output_dir/npt.tpr"
  gmx mdrun -deffnm "$pdb_output_dir/npt.tpr"

  gmx grompp -f production.mdp -c "$pdb_output_dir/npt.gro" -p "$pdb_output_dir/topol.top" -o "$pdb_output_dir/production.tpr"
  gmx mdrun -deffnm "$pdb_output_dir/production"

  gmx rms -s "$pdb_output_dir/production.tpr" -f "$pdb_output_dir/production.xtc" -o "$pdb_output_dir/rmsd.xvg"
  gmx rmsf -s "$pdb_output_dir/production.tpr" -f "$pdb_output_dir/production.xtc" -o "$pdb_output_dir/rmsf.xvg"
  gmx distance -s "$pdb_output_dir/production.tpr" -f "$pdb_output_dir/production.xtc" -select 'com of group "CO2" plus com of group "Rubisco"' -o "$pdb_output_dir/distances.xvg"
}

main "$@"
