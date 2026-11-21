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

# The dataset's Hub revision, in order of preference:
#   1. MVA_DATASET_REVISION, if the operator sets it deliberately
#   2. the cache the download left behind -- authoritative while it exists
#   3. DATA.md, which records it and is committed, so it outlives the data
# Three is needed because the cache lives inside C1: once that is purged, the only record
# left is the repository's own, and a resumed run must still be able to name the snapshot.
repo_root="$(cd "$(dirname "$0")/.." && pwd)"
source_of_revision=""
revision="${MVA_DATASET_REVISION:-}"
[ -n "$revision" ] && source_of_revision="MVA_DATASET_REVISION"

if [ -z "$revision" ]; then
  trees=("$MVA_DATA_ROOT"/raw/.cache/huggingface/trees/*.json)
  if [ -e "${trees[0]}" ]; then
    if [ "${#trees[@]}" -ne 1 ]; then
      echo "${#trees[@]} Hub revisions cached; cannot tell which snapshot was used." >&2
      exit 2
    fi
    revision="$(basename "${trees[0]}" .json)"
    source_of_revision="the Hub cache"
  fi
fi

if [ -z "$revision" ]; then
  revision="$(grep -oE '\b[0-9a-f]{40}\b' "$repo_root/DATA.md" | head -1 || true)"
  source_of_revision="DATA.md (the cache is gone with C1)"
fi

if [ -z "$revision" ]; then
  echo "no Hub revision available: not in MVA_DATASET_REVISION, no cache under" >&2
  echo "$MVA_DATA_ROOT/raw/.cache/huggingface/trees/, and none recorded in DATA.md." >&2
  echo "Do not guess it -- the attestation identifies the destroyed data by it alone." >&2
  exit 2
fi

echo "dataset revision: $revision"
echo "  from          : $source_of_revision"
echo "attestor        : $MVA_ATTESTOR_NAME ($MVA_ATTESTOR_ROLE)"
echo

cd "$repo_root"
exec "$PYTHON" -m src.purge --revision-sha "$revision" "$@"
