#!/bin/bash

set -e

gmx_serial="autoload gromacs/2021.2 "
gmx_parallel="gromacs/2021.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0-cuda-11.5.0"

pdb=""
co2_molecules=1
output_dir=""
parallel=0

usage() {
  echo "usage: $0 [-h] [-i <pdb>] [-n <co2-molecules>] [-o <output-dir>] [-p]"
}

switch_module() {
  module purge

  module load profile/lifesc
  module load "$1"
}

main() {
  while getopts ":hi:n:o:p" opt; do
    case $opt in
      h) usage; exit 0 ;;
      i) pdb="$OPTARG" ;;
      n) co2_molecules="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      p) parallel=1 ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  config_dir="$(realpath "$(cd "$(dirname "$0")" && pwd)/../config")"

  if [[ ! -d "$config_dir" ]]; then
    echo "error: '$config_dir' does not exist" >&2
    exit 1
  fi

  required_files=(
    co2.itp
    co2.pdb
    ions.mdp
    minimization.mdp
    npt.mdp
    nvt.mdp
    production.mdp
  )

  for file in "${required_files[@]}"; do
    if [[ ! -f "$config_dir/$file" ]]; then
      echo "error: missing file $config_dir/$file" >&2
      exit 1
    fi
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '-i' must be a valid .pdb file" >&2; usage; exit 1
  fi

  pdb="$(realpath "$pdb")"

  if [[ ! "$co2_molecules" =~ ^[0-9]+$ || "$co2_molecules" -lt 1 ]]; then
    echo "error: '-n' must be a positive integer" >&2; usage; exit 1
  fi

  if [[ -z "$output_dir" ]]; then
    output_dir="./$(basename "$pdb" .pdb)"
  fi

  mkdir -p "$output_dir/intermediate"

  pushd "$(pwd)" > /dev/null
  cd "$output_dir/intermediate"
  trap 'popd > /dev/null' EXIT

  gmx pdb2gmx -f "$pdb" -o ./processed.gro -p ./topol.top -ff amber99sb-ildn -water tip3p

  gmx editconf -f "$config_dir/co2.pdb" -o "./co2.gro"

  cp "$config_dir/co2.itp" "./co2.itp"
  sed "/#include \"amber99sb-ildn.ff\/forcefield.itp\"/a\\
\\
; Include topology for co2\
\\
#include \"co2.itp\"\
\\
" ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top
  printf "CO2%-$((18 - ${#co2_molecules}))s%d\n" "" "$co2_molecules" >> ./topol.top

  gmx insert-molecules -f ./processed.gro -ci ./co2.gro -o ./complex.gro -nmol "$co2_molecules"

  gmx editconf -f ./complex.gro -o ./box.gro -bt dodecahedron -d 1.0

  gmx solvate -cp ./box.gro -cs spc216.gro -p ./topol.top -o ./solvated.gro

  gmx grompp -f "$config_dir/ions.mdp" -c ./solvated.gro -p ./topol.top -o ./ions.tpr
  echo SOL | gmx genion -s ./ions.tpr -o ./ions.gro -p ./topol.top -pname NA -nname CL -neutral

  gmx grompp -f "$config_dir/minimization.mdp" -c ./ions.gro -p ./topol.top -o ./minimization.tpr
  if [[ "$parallel" -eq 1 ]] && command -v srun &> /dev/null; then
    switch_module "$gmx_parallel"
    srun gmx-mpi mdrun -v -deffnm ./minimization
    switch_module "$gmx_serial"
  else
    gmx mdrun -v -deffnm ./minimization
  fi

  gmx make_ndx -f ./co2.gro -o ./co2.ndx << EOF
0 & ! a H*
q
EOF
  echo "CO2" | gmx genrestr -f ./co2.gro -n ./co2.ndx -o ./posre_co2.itp -fc 1000 1000 1000
  sed "/#include \"co2.itp\"/a\\
\\
; Ligand position restraints\
\\
#ifdef POSRES\
\\
#include \"posre_co2.itp\"\
\\
#endif\
\\
" ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top

  gmx make_ndx -f ./minimization.gro -o ./index.ndx << EOF
1 | 13
q
EOF

  gmx grompp -f "$config_dir/nvt.mdp" -c ./minimization.gro -r ./minimization.gro -p ./topol.top -n ./index.ndx -o ./nvt.tpr
  if [[ "$parallel" -eq 1 ]] && command -v srun &> /dev/null; then
    switch_module "$gmx_parallel"
    srun gmx-mpi mdrun -v -deffnm ./nvt
    switch_module "$gmx_serial"
  else
    gmx mdrun -v -deffnm ./nvt
  fi

  gmx grompp -f "$config_dir/npt.mdp" -c ./nvt.gro -t ./nvt.cpt -r ./nvt.gro -p ./topol.top -n ./index.ndx -o ./npt.tpr
  if [[ "$parallel" -eq 1 ]] && command -v srun &> /dev/null; then
    switch_module "$gmx_parallel"
    srun gmx-mpi mdrun -v -deffnm ./npt
    switch_module "$gmx_serial"
  else
    gmx mdrun -v -deffnm ./npt
  fi

  gmx grompp -f "$config_dir/production.mdp" -c ./npt.gro -t ./npt.cpt -p ./topol.top -n ./index.ndx -o ./production.tpr
  if [[ "$parallel" -eq 1 ]] && command -v srun &> /dev/null; then
    switch_module "$gmx_parallel"
    srun gmx-mpi mdrun -v -deffnm ./production
    switch_module "$gmx_serial"
  else
    gmx mdrun -v -deffnm ./production
  fi

  echo "Protein" | gmx trjconv -s ./production.tpr -f ./production.xtc -o ./production_centered.xtc -center -pbc mol -ur compact

  gmx rms -s ./production.tpr -f ./production_centered.xtc -o ../rmsd.xvg
  gmx rmsf -s ./production.tpr -f ./production_centered.xtc -o ../rmsf.xvg
  gmx distance -s ./production.tpr -f ./production_centered.xtc -select 'com of group "CO2" plus com of group "Rubisco"' -o ../distances.xvg
}

main "$@"
