#!/bin/bash

PML="./script.pml"

pdb=""
output_file=""

usage() {
  echo "usage: $0 [-h] [-p <pdb>] [-o <output>]"
}

main() {
  while getopts ":hp:o:" opt; do
    case "$opt" in
      h) usage; exit 0 ;;
      p) pdb="$OPTARG" ;;
      o) output_file="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '$pdb' must be a valid .pdb file" >&2; usage; exit 1
  fi

  if [[ -z "$output_file" ]]; then
    output_file="./$(basename "$pdb" .pdb)_clean.pdb"
  elif [[ -d "$output_file" ]]; then
    output_file="${output_file}/$(basename "$pdb" .pdb)_clean.pdb"
  elif [[ ! "$output_file" =~ \.pdb$ ]]; then
    output_file="${output_file}.pdb"
  fi

  mkdir -p "$(dirname "$output_file")"

  cat <<EOF > "$PML"
load $pdb
remove solvent
remove hetatm
save $output_file
quit
EOF

  pymol -cq "$PML"
  rm "$PML"
}

main "$@"
