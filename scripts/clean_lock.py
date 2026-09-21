"""Strip the machine-specific and ToS-bearing parts out of the exported conda lock.

``conda env export`` writes three things this repository must not commit: the ``defaults``
channel, which conda merges in from base config and which carries Anaconda's Terms of
Service; a ``prefix:`` line holding the exporting machine's home directory; and no
indication of *when* it was exported, which is how the lock in this repository went stale
for ten days without anyone noticing.

Usage::

    conda env export > environment.lock.yml && python scripts/clean_lock.py

Run it from anywhere: the lock is located relative to this file, not to the caller.
"""

from __future__ import annotations

import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
LOCK = REPO_ROOT / "environment.lock.yml"

HEADER = """#
# Exact pins, builds included -- the reproducibility artifact CLAUDE.md requires, and the
# environment that produced the committed results.
#
#   conda env create -f environment.lock.yml
#
# Platform: linux-64. The pipeline is Linux/WSL-only regardless -- pysam and bcftools
# publish no Windows builds.
#
# environment.yml is pinned too (since 2026-09-21), so the two agree rather than splitting
# into a loose spec and an exact one. It stays the portable, readable form: use it if this
# lock will not solve on another machine, and say in the report that the environment
# differed.
#
# Three edits are applied to conda's own export, and all three matter:
#   * the `defaults` channel is removed. conda merges it in from base config, nothing here
#     comes from it, and it carries Anaconda's Terms of Service -- a judge should not hit
#     a licensing prompt to reproduce a CC-BY-4.0 result (see environment.yml).
#   * the `prefix:` line is removed. It embeds a local home directory and has no bearing
#     on the build.
#   * the generation date is stamped above, because a lock that predates a dependency is
#     worse than no lock: this one was stale for ten days after h5py and matplotlib
#     landed, and nothing said so. tests/test_reproducibility.py now asserts the lock
#     covers every declared dependency.
#
# Regenerate with:
#   conda env export > environment.lock.yml && python scripts/clean_lock.py
"""


def main() -> None:
    if not LOCK.is_file():
        raise SystemExit(
            f"{LOCK} not found. Export it first:\n"
            "  conda env export > environment.lock.yml && python scripts/clean_lock.py")

    stamp = f"# Generated: conda env export -n mva-track2  ({datetime.date.today()})\n"
    kept, removed = [], []
    for line in LOCK.read_text(encoding="utf-8").splitlines():
        if line.strip() == "- defaults":
            removed.append("defaults channel")
        elif line.startswith("prefix:"):
            removed.append("prefix line")
        elif line.startswith("#"):
            continue  # a previous run's header; this run writes a fresh one
        else:
            kept.append(line)

    LOCK.write_text(stamp + HEADER + "\n".join(kept) + "\n", encoding="utf-8")
    print(f"wrote {LOCK.relative_to(REPO_ROOT)}: {len(kept)} line(s), "
          f"removed {', '.join(removed) or 'nothing'}")


if __name__ == "__main__":
    main()
