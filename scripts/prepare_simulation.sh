#!/bin/bash

set -e

pdb=""
gas="co2"
molecules=1
output_dir=""

usage() {
  echo "usage: $0 [-h] [-i <pdb>] [-g <gas>] [-n <molecules>] [-o <output-dir>]"
}

main() {
  while getopts ":hi:g:n:o:" opt; do
    case $opt in
      h) usage; exit 0 ;;
      i) pdb="$OPTARG" ;;
      g) gas="${OPTARG,,}" ;;
      n) molecules="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  scripts_dir="$(cd "$(dirname "$0")" && pwd)"

  config_dir="$(realpath "$scripts_dir/../config")"

  if [[ ! -d "$config_dir" ]]; then
    echo "error: '$config_dir' does not exist" >&2
    exit 1
  fi

  required_files=(
    ions.mdp
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

  if [ ! -f "data/raw/gases/$gas/$gas.pdb" ] || [ ! -f "data/raw/gases/$gas/$gas.itp" ]; then
    echo "error: gas '$gas' not found in data/raw/gases/" >&2
    exit 1
  fi

  gas_pdb="$(realpath "data/raw/gases/$gas.pdb")"
  gas_itp="$(realpath "data/raw/gases/$gas.itp")"

  if [[ ! "$molecules" =~ ^[0-9]+$ || "$molecules" -lt 1 ]]; then
    echo "error: '-n' must be a positive integer" >&2; usage; exit 1
  fi

  if [[ -z "$output_dir" ]]; then
    output_dir="./$(basename "$pdb" .pdb)"
  fi

  mkdir -p "$output_dir/intermediate"

  pushd "$(pwd)" > /dev/null
  cd "$output_dir/intermediate"
  trap 'popd > /dev/null' EXIT

  gmx pdb2gmx -f "$pdb" -o ./processed.gro -p ./topol.top -ff amber99sb-ildn -water tip3p -ignh -merge all

  gmx editconf -f "$gas_pdb" -o "./$gas.gro"

  cp "$gas_itp" "./$gas.itp"
  sed "/#include \"amber99sb-ildn.ff\/forcefield.itp\"/a\\
\\
; Include topology for $gas\
\\
#include \"$gas.itp\"\
\\
" ./topol.top > ./topol.top.tmp && mv ./topol.top.tmp ./topol.top
  printf "%s%-$((21 - ${#gas} - ${#molecules}))s%d\n" "$gas" "" "$molecules" >> ./topol.top

  gmx insert-molecules -f ./processed.gro -ci "./$gas.gro" -o ./complex.gro -nmol "$molecules"

  gmx editconf -f ./complex.gro -o ./box.gro -bt dodecahedron -d 1.0

  gmx solvate -cp ./box.gro -cs spc216.gro -p ./topol.top -o ./solvated.gro

  gmx grompp -f "$config_dir/ions.mdp" -c ./solvated.gro -p ./topol.top -o ./ions.tpr
  echo SOL | gmx genion -s ./ions.tpr -o ./ions.gro -p ./topol.top -pname NA -nname CL -neutral
}

main "$@"
