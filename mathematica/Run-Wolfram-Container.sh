#!/usr/bin/env bash

set -e

MATHEMATICA_FILE_NAME="${1:-ValidateMultifoldSumsOfPowers.txt}"

CURRENT_DIR="$(pwd)"

cd "$(dirname "$0")"

LICENSE_FILE="$PWD/.mathpass"
if [[ ! -f "$LICENSE_FILE" ]]; then
  echo "Wolfram license file not found: $LICENSE_FILE" >&2
  echo "Activate Wolfram Engine and save its mathpass file there." >&2
  exit 1
fi

SCRIPT_CODE="$(cat "$MATHEMATICA_FILE_NAME")"

echo "================================================= BEGIN: $MATHEMATICA_FILE_NAME ================================================="

docker run --rm \
  -v "$LICENSE_FILE:/home/wolframengine/.WolframEngine/Licensing/mathpass:ro" \
  wolframresearch/wolframengine:15.0 \
  wolframscript -code "$SCRIPT_CODE"

echo "================================================= END: $MATHEMATICA_FILE_NAME ================================================="

cd "$CURRENT_DIR"
