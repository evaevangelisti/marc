#!/bin/bash

PDB="1iwa"

# Paths
# ----------------------------------------------------------------

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

## Data paths
## --------------------------------

DATA_DIR="$ROOT/data"

### Raw paths

RAW_DIR="$DATA_DIR/raw/"

PROTEINS_DIR="$RAW_DIR/proteins"
GASES_DIR="$RAW_DIR/gases"

### Interim paths

RAW_DIR="$DATA_DIR/interim/"

### Processed paths

PROCESSED_DIR="$DATA_DIR/processed"

## Config paths
## --------------------------------

CONFIG_DIR="$ROOT/config"

### MDP paths

MDP_DIR="$CONFIG_DIR/mdp"

MDP_IONS="$MDP_DIR/ions.mdp"
MDP_MINIMIZATION="$MDP_DIR/minimization.mdp"
MDP_NPT="$MDP_DIR/npt.mdp"
MDP_NVT="$MDP_DIR/nvt.mdp"
MDP_PRODUCTION="$MDP_DIR/production.mdp"
