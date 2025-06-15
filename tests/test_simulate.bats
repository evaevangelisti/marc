#!/usr/bin/env bats

SCRIPT="./scripts/simulate.sh"

TMP_DIR="./tmp"
INPUT="$TMP_DIR/input.pdb"
OUTPUT="$TMP_DIR/output"

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
  [ -d "./$(basename "$INPUT" .pdb)" ]
  rm -rf "./$(basename "$INPUT" .pdb)"
}

@test "valid input and output files" {
  run "$SCRIPT" -p "$INPUT" -o "$OUTPUT"
  [ "$status" -eq 0 ]
  [ -d "$OUTPUT" ]
}
