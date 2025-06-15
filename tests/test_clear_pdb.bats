#!/usr/bin/env bats

SCRIPT="./scripts/clear_pdb.sh"

TMP_DIR="./tmp"
INPUT="$TMP_DIR/input.pdb"
OUTPUT="$TMP_DIR/output.pdb"

setup() {
  if [[ ! -x "$SCRIPT" ]]; then
    echo "error: script '$SCRIPT' not found or not executable" >&2
    exit 1
  fi

  mkdir -p "$TMP_DIR"
  cat <<EOF > "$INPUT"
HEADER    TEST
ATOM      1  N   MET A   1      11.104  13.207  10.456  1.00  0.00           N
EOF
}

teardown() {
  rm -rf "$TMP_DIR"
}

@test "display usage information" {
  run "$SCRIPT" -h
  [ "$status" -eq 0 ]
  [[ "$output" =~ "usage:" ]]
}

@test "missing input file" {
  run "$SCRIPT" -o "$OUTPUT"
  [ "$status" -ne 0 ]
  [[ "$output" =~ must\ be\ a\ valid\ .pdb\ file ]]
}

@test "non-existent input file" {
  run "$SCRIPT" -p "$TMP_DIR/non_existent.pdb" -o "$OUTPUT"
  [ "$status" -ne 0 ]
  [[ "$output" =~ must\ be\ a\ valid\ .pdb\ file ]]
}

@test "invalid input file format" {
  touch "$TMP_DIR/input.txt"
  run "$SCRIPT" -p "$TMP_DIR/input.txt" -o "$OUTPUT"
  [ "$status" -ne 0 ]
  [[ "$output" =~ must\ be\ a\ valid\ .pdb\ file ]]
}

@test "default output behavior" {
  run "$SCRIPT" -p "$INPUT"
  [ "$status" -eq 0 ]
  [ -f "./$(basename "$INPUT" .pdb)_clean.pdb" ]
  rm -f "./$(basename "$INPUT" .pdb)_clean.pdb"
}

@test "output directory only" {
  mkdir -p "$TMP_DIR/output_dir"
  run "$SCRIPT" -p "$INPUT" -o "$TMP_DIR/output_dir/"
  [ "$status" -eq 0 ]
  [ -f "$TMP_DIR/output_dir/$(basename "$INPUT" .pdb)_clean.pdb" ]
}

@test "output file without extension" {
  run "$SCRIPT" -p "$INPUT" -o "$TMP_DIR/output"
  [ "$status" -eq 0 ]
  [ -f "$TMP_DIR/output.pdb" ]
}

@test "nested output directory" {
  run "$SCRIPT" -p "$INPUT" -o "$TMP_DIR/nested_dir/output.pdb"
  [ "$status" -eq 0 ]
  [ -f "$TMP_DIR/nested_dir/output.pdb" ]
}

@test "valid input and output files" {
  run "$SCRIPT" -p "$INPUT" -o "$OUTPUT"
  [ "$status" -eq 0 ]
  [ -f "$OUTPUT" ]
}
