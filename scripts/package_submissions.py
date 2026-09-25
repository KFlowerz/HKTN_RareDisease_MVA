"""Assemble each track's submission package and check that nothing is missing.

The pipeline's output is organised by **layer**, not by track, and deliberately stays that
way: L0 feeds both tracks. Its causal-gene call is Track 1's answer and Track 2's premise,
so splitting ``results/l0_genomics`` between two track directories would claim a separation
that does not exist in the science.

What *is* per-track is the **submission**: which files get uploaded where, and whether they
are all present. This builds that, under ``results/submissions/<track>/``:

    results/submissions/
      track1/
        <name>.csv        the predictions file -- PATIENT DATA, never committed
        report.md         a copy of docs/report_track1.md, named for upload
        SUBMIT.md         the checklist, with what is ready and what is missing
      track2/
        report.md         a copy of docs/report_track2.md
        dossier/          a copy of the rendered dossier (index, exclusions, pages, figures)
        SUBMIT.md

Everything here is a **copy**, and the checklist says so. The sources stay canonical:
``docs/`` for the reports, ``results/l5/`` for the dossier. Editing a packaged copy is a
mistake, and the checklist names the source to edit instead.

Usage::

    python scripts/package_submissions.py            # both tracks
    python scripts/package_submissions.py --track 1  # one
    python scripts/package_submissions.py --name kflowerz_mane-clinvar-gnomad

Guardrail
    ``results/`` is gitignored and is custody location ``C3``. The Track 1 predictions file
    carries the subject's coordinates (decision D21) and **must not** be copied into the
    repository; this script only ever writes under ``results/`` and refuses a destination
    anywhere else.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "results"
SUBMISSIONS = RESULTS / "submissions"
DOCS = REPO_ROOT / "docs"


def _assert_inside_results(path: Path) -> None:
    """Refuse to write anywhere but ``results/``.

    The Track 1 package contains patient data. A destination outside the gitignored,
    custody-tracked output directory is a bug worth failing on rather than a path to fix
    after the fact.
    """
    resolved = path.resolve()
    if RESULTS.resolve() not in resolved.parents and resolved != RESULTS.resolve():
        raise SystemExit(
            f"refusing to write {resolved}: submission packages carry patient data and "
            "belong under results/ (gitignored, custody C3, decision D21)")


def _copy(source: Path, destination: Path) -> bool:
    if not source.exists():
        return False
    _assert_inside_results(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, destination, dirs_exist_ok=True)
    else:
        shutil.copy2(source, destination)
    return True


#: Files the **report itself** links to by relative path. Copied beside it because a
#: packaged report whose reference list is a dead link undercuts the one claim this
#: project defends hardest -- that every citation was resolved and retraction-checked.
#: None is patient data; all three are committed.
#:
#: The boundary is deliberate and stops here. These three link onward into the repository
#: -- to the config, the decision log, CLAUDE.md -- and following that chain would end in
#: packaging most of the repo. The submission carries the GitHub URL for exactly that, and
#: each SUBMIT.md says so rather than leaving the reader to find a dead link.
SUPPORTING = ("references.md", "references.bib", "glossary.md")


def _copy_supporting(out: Path) -> bool:
    """Copy what the reports link to, so those links resolve inside the package."""
    return all(_copy(DOCS / name, out / name) for name in SUPPORTING)


def _supporting_names() -> str:
    """The supporting filenames, for the checklist -- generated, so it cannot drift."""
    return ", ".join(f"`{name}`" for name in SUPPORTING)


def _tick(ready: bool) -> str:
    return "[x]" if ready else "[ ]"


def package_track1(name: str) -> Path:
    """Assemble the Track 1 package. Returns its directory."""
    out = SUBMISSIONS / "track1"
    _assert_inside_results(out)
    out.mkdir(parents=True, exist_ok=True)

    # The builder writes its CSV here directly, so this only checks rather than copies --
    # a build directory plus a copy of it is two things that can disagree.
    csv_destination = out / f"{name}.csv"
    csv_ready = csv_destination.exists()
    report_ready = _copy(DOCS / "report_track1.md", out / "report.md")
    supporting_ready = _copy_supporting(out)

    (out / "SUBMIT.md").write_text(f"""# Track 1 — Variant Prediction: what to upload

Packaged {date.today()}. Everything in this directory is a **copy**; edit the sources, not
these files, then re-run `python scripts/package_submissions.py`.

## Upload at the Space's Track 1 tab

| | File | Source of truth | Ready |
|---|---|---|---|
| Predictions CSV | `{csv_destination.name}` | `python scripts/build_track1_submission.py --name {name}` | {_tick(csv_ready)} |
| Report (`.md` or PDF) | `report.md` | `docs/report_track1.md` | {_tick(report_ready)} |
| Supporting files the report links to | {_supporting_names()} | `docs/` | {_tick(supporting_ready)} |
| GitHub URL | — | must start with `https://github.com/` | [ ] |
| Team / display name | — | your choice; keep it identical across teammates | [ ] |

## Before you upload

- [ ] **Confirm `proband_id`.** The builder defaults to `PROBAND01`; check it against the
      dataset README, because a wrong id is rejected outright and costs a slot.
- [ ] **Rename the CSV** to include your username and a short approach name, as the Space
      asks — for example `jane-doe_mane-clinvar-gnomad.csv`. Pass `--name` to this script.
- [ ] **Six submissions are available**, and only the highest-scoring one is featured.

## A known property of the copies

`report.md`'s own links all resolve inside this directory. The **supporting** files link
onward to things that are not packaged — `config/pipeline.yaml`, `mngmt/decisions.md`,
`CLAUDE.md` — because packaging those would mean packaging most of the repository. Follow
them at the GitHub URL above; that is what it is for.

## Do not

- **Do not copy `{csv_destination.name}` into the repository.** It carries the subject's
  variant coordinates. It lives here because `results/` is gitignored and is custody
  location `C3`; the repository must never hold it (decision D21).
- **Do not read the organizers' `groundtruth.py` or their ground-truth dataset.** It would
  contaminate every submission this project makes and cannot be undone.
""", encoding="utf-8")
    return out


def package_track2() -> Path:
    """Assemble the Track 2 package. Returns its directory."""
    out = SUBMISSIONS / "track2"
    _assert_inside_results(out)
    out.mkdir(parents=True, exist_ok=True)

    report_ready = _copy(DOCS / "report_track2.md", out / "report.md")
    supporting_ready = _copy_supporting(out)
    dossier_ready = _copy(RESULTS / "l5", out / "dossier")
    index = out / "dossier" / "index.html"

    (out / "SUBMIT.md").write_text(f"""# Track 2 — Drug Repurposing: what to upload

Packaged {date.today()}. Everything in this directory is a **copy**; edit the sources, not
these files, then re-run `python scripts/package_submissions.py`.

## Upload at the Space's Track 2 tab

| | File | Source of truth | Ready |
|---|---|---|---|
| Written report (`.md` or PDF) | `report.md` | `docs/report_track2.md` | {_tick(report_ready)} |
| Candidate dossier | `dossier/index.html` | `results/l5/` (regenerate, do not edit) | {_tick(dossier_ready)} |
| Supporting files the report links to | {_supporting_names()} | `docs/` | {_tick(supporting_ready)} |
| GitHub URL | — | must start with `https://github.com/` | [ ] |
| 3-minute video | — | script, shot list and privacy rules in `docs/video_script.md` — **not yet recorded** | [ ] |
| Methods description form | — | the Space's `.xlsx` template | [ ] |

## Before you upload

- [ ] **Re-run the claim check**: `python scripts/verify_report_claims.py`. It fails if the
      report and the artifacts have drifted apart after a re-run.
- [ ] **Open `dossier/index.html`** and read it as a judge would — exclusions first, then
      the two tiers.
- [ ] **Three submissions are available**, and only the latest is reviewed by the panel.
- [ ] **Record the AI disclosure** in the methods form: local open-weights model for the
      pipeline's reasoning step, and the development assistant with its plan and
      training setting (decision D17).

## Notes

- `report.md`'s own links all resolve inside this directory. The **supporting** files link
  onward to things that are not packaged — `config/pipeline.yaml`, `mngmt/decisions.md`,
  `CLAUDE.md` — because packaging those would mean packaging most of the repository. Follow
  them at the GitHub URL above.
- The dossier is a copy of `results/l5/`. Regenerate the source with
  `PYTHONHASHSEED=42 python -m src.pipeline --only l5_report`, then re-run this script.
- The dossier contains no patient identifiers — `src/l5_report/publish_guard.py` fails
  closed on every string it writes — but it is packaged here rather than committed because
  `results/` is where pipeline output lives.
""", encoding="utf-8")
    return out, index


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--track", choices=["1", "2"], help="package one track only")
    parser.add_argument("--name", default="track1_submission",
                        help="basename for the Track 1 CSV; include your username")
    args = parser.parse_args()

    if not RESULTS.is_dir():
        raise SystemExit("results/ not found. Run the pipeline first.")

    if args.track in (None, "1"):
        out = package_track1(args.name)
        files = sorted(p.name for p in out.iterdir())
        print(f"track 1 -> {out.relative_to(REPO_ROOT)}: {', '.join(files)}")

    if args.track in (None, "2"):
        out, index = package_track2()
        pages = len(list((out / "dossier" / "candidates").glob("*.html"))) \
            if (out / "dossier" / "candidates").is_dir() else 0
        print(f"track 2 -> {out.relative_to(REPO_ROOT)}: report.md, dossier/ "
              f"({pages} candidate pages)")
        if not index.exists():
            print("  WARNING: dossier/index.html missing -- run L5 first", file=sys.stderr)

    print()
    print("Read SUBMIT.md in each directory for the upload checklist.")
    print("results/ is gitignored (custody C3). Nothing here is committed.")


if __name__ == "__main__":
    main()
