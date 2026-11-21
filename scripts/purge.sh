#!/usr/bin/env bash
# Run the custody purge with this machine's locations and attestor already filled in.
#
# Why this exists: the command it wraps is ~110 characters ending in a 40-character hash,
# and a terminal with bracketed paste enabled appends a stray "~" to a pasted line. That
# corrupted the hash twice; the pre-flight caught it both times and refused, but the fix is
# to make the command short enough to type rather than paste.
#
#   bash scripts/purge.sh              # dry run, deletes nothing
#   bash scripts/purge.sh --confirm    # deletes, then writes the attestation
#
# The revision hash is read from the Hub cache rather than typed, so it cannot be mistyped.
# Everything else the pre-flight needs is set below. See docs/data_custody.md.
set -euo pipefail

export MVA_DATA_ROOT="${MVA_DATA_ROOT:-$HOME/mva-data}"
export HF_HOME="${HF_HOME:-$MVA_DATA_ROOT/hf-cache}"
export CLAUDE_PROJECTS="${CLAUDE_PROJECTS:-/mnt/c/Users/kflor/.claude/projects}"
export MVA_ATTESTOR_NAME="${MVA_ATTESTOR_NAME:-Karen Flores}"
export MVA_ATTESTOR_ROLE="${MVA_ATTESTOR_ROLE:-Maintainer}"

PYTHON="${MVA_PYTHON:-$HOME/miniconda3/envs/mva-track2/bin/python}"
if [ ! -x "$PYTHON" ]; then
  echo "no interpreter at $PYTHON -- set MVA_PYTHON to the environment's python" >&2
  exit 2
fi

# The dataset's Hub revision, read from the cache the download left behind. Exactly one
# tree file must be present: zero means the data is already gone, more than one means the
# snapshot is ambiguous and the attestation would name the wrong one.
trees=("$MVA_DATA_ROOT"/raw/.cache/huggingface/trees/*.json)
if [ ! -e "${trees[0]}" ]; then
  echo "no Hub revision found under $MVA_DATA_ROOT/raw/.cache/huggingface/trees/." >&2
  echo "If the data is already purged, pass --revision-sha yourself against the value in" >&2
  echo "the previous attestation; do not guess it." >&2
  exit 2
fi
if [ "${#trees[@]}" -ne 1 ]; then
  echo "${#trees[@]} Hub revisions cached; cannot tell which snapshot was used." >&2
  exit 2
fi
revision="$(basename "${trees[0]}" .json)"

echo "dataset revision: $revision"
echo "attestor        : $MVA_ATTESTOR_NAME ($MVA_ATTESTOR_ROLE)"
echo

cd "$(dirname "$0")/.."
exec "$PYTHON" -m src.purge --revision-sha "$revision" "$@"
