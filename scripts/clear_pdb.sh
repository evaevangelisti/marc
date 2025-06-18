#!/bin/bash

set -e

pdb=""
output_path=""

usage() {
  echo "usage: $0 [-h] [-p <pdb>] [-o <output>]"
}

main() {
  while getopts ":hp:o:" opt; do
    case "$opt" in
      h) usage; exit 0 ;;
      p) pdb="$OPTARG" ;;
      o) output_path="$OPTARG" ;;
      :) echo "error: missing argument for option '-$OPTARG'" >&2; usage; exit 1 ;;
      ?) echo "error: unknown option '-$OPTARG'" >&2; usage; exit 1 ;;
    esac
  done

  if [[ -z "$pdb" || ! -f "$pdb" || "$pdb" != *.pdb ]]; then
    echo "error: '$pdb' must be a valid .pdb file" >&2; usage; exit 1
  fi

  if [[ -z "$output_path" ]]; then
    output_path="./$(basename "$pdb" .pdb)_clean.pdb"
  elif [[ -d "$output_path" ]]; then
    output_path="${output_path}/$(basename "$pdb" .pdb)_clean.pdb"
  elif [[ ! "$output_path" =~ \.pdb$ ]]; then
    output_path="${output_path}.pdb"
  fi

  mkdir -p "$(dirname "$output_path")"

  tmp_dir=$(mktemp -d)
  trap 'rm -rf -- "$tmp_dir"' EXIT

  cat <<EOF > "$tmp_dir/script.pml"
load $pdb
remove solvent
remove hetatm
save $output_path
quit
EOF

  pymol -cq "$tmp_dir/script.pml"
}

main "$@"
