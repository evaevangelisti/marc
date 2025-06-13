#!/bin/bash

input_dir="."
output_dir="."

usage() {
  echo "usage: $0 [-h] [-i <input-dir>] [-o <output-dir>]"
}

main() {
  while getopts ":hi:o:" opt; do
    case $opt in
      h) usage; exit 0 ;;
      i) input_dir="$OPTARG" ;;
      o) output_dir="$OPTARG" ;;
      :) echo "missing argument for '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "unknown option: '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [ ! -d "$input_dir" ]; then
    echo "input directory '$input_dir' does not exist" >&2
    exit 1
  fi

  mkdir -p "$output_dir"

  for pdb in "$input_dir"/*.pdb; do
    if [ -f "$pdb" ]; then
      pymol -cq <<EOF
load $pdb
remove solvent
remove hetatm
save "$output_dir"/${pdb%.*}_clean.pdb
EOF
    fi
  done
}

main "$@"
