#!/bin/bash

pdb=""
output_dir="."

usage() {
  echo "usage: $0 [-h] [-p <pdb>] [-o <output-dir>]"
}

main() {
  while getopts ":hp:o:" opt; do
    case "$opt" in
      h) usage; exit 0 ;;
      p) pdb="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '$pdb' must be a valid .pdb file" >&2
    exit 1
  fi

  mkdir -p "$output_dir"

  pymol -cq <<EOF
load "$pdb"
remove solvent
remove hetatm
save "$output_dir/${pdb%.*}_clean.pdb"
EOF

  echo "$output_dir/${pdb%.*}_clean.pdb"
}

main "$@"
