"""Strip the machine-specific and ToS-bearing parts out of the exported conda lock."""
from pathlib import Path

LOCK = Path("/mnt/d/Labs/JupyterNotebooks/Hackathons/HKTN_RareDisease_MVA/environment.lock.yml")

HEADER = """# Exact pins, builds included -- the reproducibility artifact CLAUDE.md requires, and the
# environment that produced the committed results.
#
#   conda env create -f environment.lock.yml
#
# Platform: linux-64. The pipeline is Linux/WSL-only regardless -- pysam and bcftools
# publish no Windows builds.
#
# environment.yml is now pinned too (2026-09-21), so the two agree rather than splitting
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
#   * the generation date is recorded below rather than left implicit, because a lock that
#     predates a dependency is worse than no lock: this one was stale for ten days after
#     h5py and matplotlib landed, and nothing said so.
#
# Regenerate with: conda env export | python scripts/clean_lock.py
"""

import datetime

lines = LOCK.read_text(encoding="utf-8").splitlines()
stamp = f"# Generated: conda env export -n mva-track2  ({datetime.date.today()})\n#\n"
out, skipped = [], []
for line in lines:
    if line.strip() == "- defaults":
        skipped.append("defaults channel")
        continue
    if line.startswith("prefix:"):
        skipped.append("prefix line")
        continue
    out.append(line)

LOCK.write_text(stamp + HEADER + "\n".join(out) + "\n", encoding="utf-8")
print("removed:", ", ".join(skipped) or "nothing")
print("lines:", len(out))
