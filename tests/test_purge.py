"""Custody register parsing and measurement.

These exercise the two read-only halves of ``src/purge.py``. Nothing here touches a real
custody location, and no fixture contains patient-derived content -- the register under
test is synthetic except where a case deliberately reads the committed one, which is read
and never written.
"""
from __future__ import annotations

import os
import textwrap
from pathlib import Path

import pytest

from src import purge

REGISTER = textwrap.dedent("""
    # Heading

    ## Custody register

    Status values: `planned` · `active` · `purged`

    | ID | Location | Holds | Purge method | Status |
    |---|---|---|---|---|
    | `C1` | `$TEST_DATA_ROOT/raw` (ext4, outside the repo) | The dataset as downloaded | `rm -rf`, then `fstrim` | **active** |
    | `C3` | `<repo>/results/` | Pipeline intermediates | `rm -rf` contents | **active** |
    | `C6` | Shell history, terminal scrollback | Only if a record is printed | Prevented rather than purged | n/a |
    | `C8` | The organizers' store — `SageBio/some-leaderboard`, a private HF dataset | The submitted CSV | **Not purgeable by this project.** | live once submitted |

    ## Next section

    | ID | Other | Table | Entirely | Different |
    |---|---|---|---|---|
    | `X1` | a | b | c | d |
    """)


@pytest.fixture()
def register(tmp_path, monkeypatch):
    monkeypatch.setenv("TEST_DATA_ROOT", str(tmp_path / "custody"))
    path = tmp_path / "data_custody.md"
    path.write_text(REGISTER, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- the register


def test_reads_every_row_and_stops_at_the_next_section(register):
    rows = purge.read_register(register)
    assert [row.id for row in rows] == ["C1", "C3", "C6", "C8"]


def test_expands_an_environment_placeholder(register, tmp_path):
    c1 = {row.id: row for row in purge.read_register(register)}["C1"]
    assert c1.path == (tmp_path / "custody" / "raw").resolve()


def test_expands_the_repo_placeholder(register):
    c3 = {row.id: row for row in purge.read_register(register)}["C3"]
    assert c3.path == (purge.REPO_ROOT / "results").resolve()


def test_a_row_naming_no_path_resolves_to_none(register):
    """C6 is a surface, not a directory. C8's backticks hold a Hub repo id, not a path."""
    rows = {row.id: row for row in purge.read_register(register)}
    assert rows["C6"].path is None
    assert rows["C8"].path is None


def test_unpurgeable_rows_are_not_marked_purgeable(register):
    rows = {row.id: row for row in purge.read_register(register)}
    assert rows["C1"].is_purgeable and rows["C3"].is_purgeable
    assert not rows["C6"].is_purgeable and not rows["C8"].is_purgeable


def test_emphasis_and_ticks_are_stripped_from_the_status(register):
    rows = {row.id: row for row in purge.read_register(register)}
    assert rows["C1"].status == "active"
    assert rows["C8"].status == "live once submitted"


def test_the_committed_register_needs_claude_projects_set(monkeypatch, tmp_path):
    """C7 resolves differently from Windows and WSL, and the wrong one reports 'absent'
    over a live store. The variable is mandatory so that ambiguity fails loudly."""
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "mva"))
    monkeypatch.setenv("HF_HOME", str(tmp_path / "mva" / "hf-cache"))
    monkeypatch.delenv("CLAUDE_PROJECTS", raising=False)
    with pytest.raises(ValueError, match="CLAUDE_PROJECTS"):
        purge.read_register()


def test_an_unset_variable_is_an_error_not_a_silent_skip(tmp_path, monkeypatch):
    """A row that cannot be located would be left behind while the attestation claims
    otherwise, so it must fail loudly rather than resolve to something plausible."""
    monkeypatch.delenv("TEST_DATA_ROOT", raising=False)
    path = tmp_path / "r.md"
    path.write_text(REGISTER, encoding="utf-8")
    with pytest.raises(ValueError, match="TEST_DATA_ROOT"):
        purge.read_register(path)


def test_a_duplicate_id_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setenv("TEST_DATA_ROOT", str(tmp_path))
    path = tmp_path / "r.md"
    path.write_text(REGISTER.replace("| `C3` |", "| `C1` |"), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate register id"):
        purge.read_register(path)


def test_an_unknown_status_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setenv("TEST_DATA_ROOT", str(tmp_path))
    path = tmp_path / "r.md"
    path.write_text(REGISTER.replace("| **active** |", "| probably gone |", 1), encoding="utf-8")
    with pytest.raises(ValueError, match="probably gone"):
        purge.read_register(path)


def test_a_register_without_the_heading_is_rejected(tmp_path):
    path = tmp_path / "r.md"
    path.write_text("# nothing here\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Custody register"):
        purge.read_register(path)


def test_the_committed_register_parses(monkeypatch, tmp_path):
    """The real register is the contract. If it stops parsing, the purge tool is blind."""
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "mva"))
    monkeypatch.setenv("HF_HOME", str(tmp_path / "mva" / "hf-cache"))
    monkeypatch.setenv("CLAUDE_PROJECTS", str(tmp_path / "claude" / "projects"))
    rows = purge.read_register()
    assert [row.id for row in rows] == ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]
    assert all(row.status in purge.CUSTODY_STATES or row.status in purge.UNPURGEABLE_STATES
               for row in rows)


@pytest.mark.parametrize("cell", [
    "see `.../CvYmNTTmt92D6ahmMBuLe1` for the deck",      # abbreviated link, not a path
    "at `https://claude.ai/artifact/XzuNPz` online",      # a URL, despite the // in it
    "the `SageBio/mva-hackathon-2026-leaderboard` set",   # a Hub repo id
    "Claude Code's local `.jsonl` store",                 # a suffix
])
def test_a_cell_that_only_looks_like_a_path_resolves_to_none(cell):
    """C9's links were first written as `.../<id>`, which a looser pattern read as a
    relative path and resolved against the repository."""
    assert purge._resolve_path(cell, "CX") is None


@pytest.mark.parametrize("cell", [
    "`$MVA_DATA_ROOT/raw` on ext4",
    "`<repo>/results/`",
    "`~/.claude/projects/slug`",
    "`./relative/path`",
    "`/absolute/path`",
])
def test_a_real_path_cell_still_resolves(cell, monkeypatch, tmp_path):
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path))
    assert purge._resolve_path(cell, "CX") is not None


def test_the_artifact_row_is_recorded_but_never_purged(monkeypatch, tmp_path):
    """C9 exists so the register is not read as exhaustive. It must never be acted on: it
    names no filesystem path and its status is outside the purgeable states."""
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "mva"))
    monkeypatch.setenv("HF_HOME", str(tmp_path / "mva" / "hf-cache"))
    monkeypatch.setenv("CLAUDE_PROJECTS", str(tmp_path / "claude" / "projects"))
    c9 = {row.id: row for row in purge.read_register()}["C9"]
    assert c9.path is None                  # artifact URLs are not paths
    assert not c9.is_purgeable
    assert c9.status == "retained, not patient-derived"


# --------------------------------------------------------------------------- measurement


def test_an_absent_location_measures_as_absent(tmp_path):
    state = purge.measure_path(tmp_path / "never-existed")
    assert state == purge.LocationState(exists=False, file_count=0, total_bytes=0)


def test_an_empty_directory_exists_but_holds_nothing(tmp_path):
    (tmp_path / "empty").mkdir()
    state = purge.measure_path(tmp_path / "empty")
    assert state.exists and state.file_count == 0 and state.total_bytes == 0


def test_counts_files_and_bytes_across_subdirectories(tmp_path):
    (tmp_path / "a" / "b").mkdir(parents=True)
    (tmp_path / "a" / "one").write_bytes(b"x" * 100)
    (tmp_path / "a" / "b" / "two").write_bytes(b"y" * 250)
    state = purge.measure_path(tmp_path / "a")
    assert state.file_count == 2
    assert state.total_bytes == 350


@pytest.mark.skipif(os.name == "nt", reason="symlink creation needs privilege on Windows")
def test_a_symlinked_file_is_counted_once_not_as_its_target(tmp_path):
    """The Hub cache trap: snapshots/ links to the bytes in blobs/. Following the links
    would report twice the payload; ignoring them would report an empty cache."""
    blobs = tmp_path / "blobs"
    snapshots = tmp_path / "snapshots"
    blobs.mkdir()
    snapshots.mkdir()
    (blobs / "deadbeef").write_bytes(b"z" * 10_000)
    (snapshots / "model.bin").symlink_to(blobs / "deadbeef")

    state = purge.measure_path(tmp_path)
    assert state.file_count == 2                   # the blob and the link
    assert state.total_bytes < 11_000              # the payload counted once, not twice
    assert state.total_bytes >= 10_000


@pytest.mark.skipif(os.name == "nt", reason="symlink creation needs privilege on Windows")
def test_a_symlinked_directory_is_not_descended_into(tmp_path):
    real = tmp_path / "real"
    (real / "deep").mkdir(parents=True)
    (real / "deep" / "file").write_bytes(b"q" * 500)
    (tmp_path / "link").symlink_to(real, target_is_directory=True)

    state = purge.measure_path(tmp_path)
    assert state.file_count == 2                   # the real file and the directory link
    assert state.total_bytes < 1_000


def test_measure_of_an_unpathed_row_is_absent(register):
    c6 = {row.id: row for row in purge.read_register(register)}["C6"]
    assert purge.measure(c6) == purge.LocationState(False, 0, 0)


# ------------------------------------------------------------------ refusing to delete


def _location(path, identifier="CX", method="rm -rf"):
    return purge.CustodyLocation(id=identifier, path=path, holds="x", method=method,
                                 status="active")


def test_refuses_a_path_at_a_filesystem_root(tmp_path):
    root = Path(tmp_path.anchor)
    with pytest.raises(ValueError, match="filesystem root"):
        purge._assert_deletable(_location(root))


def test_refuses_the_repository_itself(tmp_path):
    with pytest.raises(ValueError, match="repository working tree"):
        purge._assert_deletable(_location(purge.REPO_ROOT))


def test_refuses_a_repository_path_that_is_not_output(tmp_path):
    """A typo turning <repo>/results/ into <repo>/src/ must not delete the source."""
    with pytest.raises(ValueError, match="inside the repository"):
        purge._assert_deletable(_location(purge.REPO_ROOT / "src"))


def test_allows_the_two_repository_output_paths(tmp_path):
    purge._assert_deletable(_location(purge.REPO_ROOT / "results"))
    purge._assert_deletable(_location(purge.REPO_ROOT / "notebooks"))


def test_allows_a_path_outside_the_repository(tmp_path):
    purge._assert_deletable(_location(tmp_path / "custody" / "raw"))


# ------------------------------------------------------------------ deleting, on fixtures


def test_delete_removes_a_tree(tmp_path):
    target = tmp_path / "custody" / "raw"
    (target / "sub").mkdir(parents=True)
    (target / "sub" / "f").write_bytes(b"x" * 10)
    purge._delete(_location(target))
    assert not target.exists()


def test_results_is_emptied_but_keeps_its_gitkeep(tmp_path, monkeypatch):
    results = tmp_path / "repo" / "results"
    (results / "l0").mkdir(parents=True)
    (results / "l0" / "calls.json").write_text("{}", encoding="utf-8")
    (results / ".gitkeep").write_text("", encoding="utf-8")
    monkeypatch.setattr(purge, "REPO_ROOT", tmp_path / "repo")

    purge._delete(_location(results))
    assert results.exists()
    assert [p.name for p in results.iterdir()] == [".gitkeep"]


def test_a_dry_run_deletes_nothing(tmp_path, monkeypatch, register):
    target = tmp_path / "custody" / "raw"
    target.mkdir(parents=True)
    (target / "f").write_bytes(b"x" * 10)
    monkeypatch.setattr(purge, "read_register", lambda *a, **k: [_location(target, "C1")])

    verdicts = purge.purge({}, confirm=False)
    assert target.exists()
    assert verdicts[0].before.file_count == 1
    assert verdicts[0].after.file_count == 1          # unchanged: nothing was removed


# ------------------------------------------------------------------ the attestation


def _clean_verdict(identifier="C1", path=Path("/tmp/x")):
    return purge.LocationVerdict(
        location=_location(path, identifier),
        before=purge.LocationState(True, 3, 300),
        after=purge.LocationState(False, 0, 0),
        verified_utc="2026-11-01T00:00:00Z",
    )


SHA = "a" * 40


def test_refuses_a_revision_that_is_not_a_full_hash():
    with pytest.raises(ValueError, match="40-character"):
        purge.write_attestation([_clean_verdict()], revision_sha="abc123")


def test_refuses_to_attest_when_a_location_still_holds_bytes():
    dirty = purge.LocationVerdict(
        location=_location(Path("/tmp/x"), "C1"),
        before=purge.LocationState(True, 3, 300),
        after=purge.LocationState(True, 1, 100),
        verified_utc="2026-11-01T00:00:00Z",
    )
    with pytest.raises(ValueError, match="still hold bytes"):
        purge.write_attestation([dirty], revision_sha=SHA)


@pytest.mark.parametrize("payload, what", [
    ('{"sha256": "x"}', "sha256"),
    ('{"a": "' + "b" * 64 + '"}', "digest"),
    ('{"a": "sample.vcf.gz"}', "genomic filename"),
    ('{"a": "chr15:40161020"}', "variant coordinate"),
    ('{"a": "HP:0001250"}', "HPO term"),
])
def test_the_outbound_scan_catches_patient_data(payload, what):
    with pytest.raises(ValueError, match="refusing to write"):
        purge._assert_carries_no_patient_data(payload)


def test_a_clean_payload_passes_the_outbound_scan():
    purge._assert_carries_no_patient_data(
        '{"id": "C1", "file_count": 13, "total_bytes": 317500000}')


# ------------------------------------------------------------------ the pre-flight

def test_preflight_rejects_a_short_revision(monkeypatch):
    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    assert "40-character" in (purge._preflight("abc123") or "")


def test_preflight_rejects_an_interpreter_that_cannot_run_pytest(monkeypatch):
    """The failure this exists for: a bare 'python' that is not on PATH, discovered only
    after deletion, leaving no data and no attestation."""
    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    monkeypatch.setenv("MVA_PYTHON", "definitely-not-an-interpreter")
    problem = purge._preflight(SHA)
    assert problem and "definitely-not-an-interpreter" in problem


def test_preflight_requires_a_named_attestor(monkeypatch):
    monkeypatch.delenv("MVA_ATTESTOR_NAME", raising=False)
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    assert "MVA_ATTESTOR_NAME" in (purge._preflight(SHA) or "")


def test_preflight_refuses_inside_a_live_assistant_session(monkeypatch):
    """C7 is the session's own transcript store. Deleted from inside a live session it is
    immediately rewritten, so it can never verify empty -- and by then C1 is gone."""
    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    monkeypatch.delenv("MVA_PYTHON", raising=False)
    monkeypatch.setenv("CLAUDECODE", "1")
    problem = purge._preflight(SHA)
    assert problem and "transcript store" in problem


def test_preflight_refuses_while_a_session_elsewhere_is_still_writing(monkeypatch, tmp_path):
    """The variable check only sees a purge launched from inside a session. One open in
    another window sets nothing in this shell, so recent writes are the signal."""
    store = tmp_path / "projects" / "d--Labs-JupyterNotebooks-Hackathons-HKTN-RareDisease-MVA"
    store.mkdir(parents=True)
    (store / "live.jsonl").write_text("{}", encoding="utf-8")   # written just now

    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "mva"))
    monkeypatch.setenv("HF_HOME", str(tmp_path / "mva" / "hf-cache"))
    monkeypatch.setenv("CLAUDE_PROJECTS", str(tmp_path / "projects"))
    for marker in ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT"):
        monkeypatch.delenv(marker, raising=False)

    problem = purge._preflight(SHA)
    assert problem and "still has it open" in problem


def test_preflight_allows_a_store_that_has_gone_quiet(monkeypatch, tmp_path):
    store = tmp_path / "projects" / "d--Labs-JupyterNotebooks-Hackathons-HKTN-RareDisease-MVA"
    store.mkdir(parents=True)
    old = store / "finished.jsonl"
    old.write_text("{}", encoding="utf-8")
    stale = old.stat().st_mtime - (purge.LIVE_SESSION_SECONDS + 60)
    os.utime(old, (stale, stale))

    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    monkeypatch.setenv("MVA_DATA_ROOT", str(tmp_path / "mva"))
    monkeypatch.setenv("HF_HOME", str(tmp_path / "mva" / "hf-cache"))
    monkeypatch.setenv("CLAUDE_PROJECTS", str(tmp_path / "projects"))
    for marker in ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT"):
        monkeypatch.delenv(marker, raising=False)

    assert purge._preflight(SHA) is None


def test_preflight_passes_in_this_environment(monkeypatch):
    monkeypatch.setenv("MVA_ATTESTOR_NAME", "A")
    monkeypatch.setenv("MVA_ATTESTOR_ROLE", "B")
    monkeypatch.delenv("MVA_PYTHON", raising=False)
    for marker in ("CLAUDECODE", "CLAUDE_CODE_ENTRYPOINT"):
        monkeypatch.delenv(marker, raising=False)
    assert purge._preflight(SHA) is None


# ------------------------------------------------------------------ the Hub CLI and journal


def test_hub_cli_is_found_beside_the_interpreter(monkeypatch, tmp_path):
    """The env's bin is not on PATH when its python is invoked by absolute path, which is
    how the purge runs. A bare 'hf' failed there and cost a partial purge."""
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    cli = fake_bin / "hf"
    cli.write_text("#!/bin/sh\n", encoding="utf-8")
    cli.chmod(0o755)
    monkeypatch.setattr(purge.sys, "executable", str(fake_bin / "python"))
    monkeypatch.setenv("PATH", "")
    assert purge.hub_cli() == str(cli)


def test_hub_cli_is_none_when_absent(monkeypatch, tmp_path):
    monkeypatch.setattr(purge.sys, "executable", str(tmp_path / "nowhere" / "python"))
    monkeypatch.setattr(purge.shutil, "which", lambda _name: None)
    assert purge.hub_cli() is None


def test_the_journal_preserves_a_before_state_across_a_crash(tmp_path, monkeypatch):
    """C1 was deleted, the run died on C2, and C1's 13 files/317 MB survived only in the
    console. A re-run would otherwise record 'before: absent' -- a false record."""
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "_purge_journal.json")
    journal = {}
    purge._write_journal(journal, "C1", purge.LocationState(True, 13, 317_500_000))

    recovered = purge._read_journal()
    assert recovered["C1"] == {"exists": True, "file_count": 13,
                               "total_bytes": 317_500_000}
    assert purge._state_from_dict(recovered["C1"]).file_count == 13


def test_a_journal_that_cannot_be_written_does_not_stop_the_purge(tmp_path, monkeypatch):
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "no" / "such" / "dir" / "j.json")
    monkeypatch.setattr(purge.Path, "mkdir",
                        lambda *a, **k: (_ for _ in ()).throw(OSError("read-only")))
    purge._write_journal({}, "C1", purge.LocationState(True, 1, 1))   # must not raise


def test_the_markdown_companion_names_no_content():
    payload = {
        "generated_utc": "2026-11-01T00:00:00Z",
        "hackathon_close_date": "2026-10-24",
        "deletion_deadline_utc": "2026-11-23T23:59:00Z",
        "dataset": {"repo_id": "r", "revision_sha": SHA, "identification_note": "n"},
        "locations": [{"id": "C1", "method": "rm -rf",
                       "before": {"file_count": 3, "total_bytes": 300},
                       "after": {"exists": False, "file_count": 0, "total_bytes": 0}}],
        "checks": [{"name": "custody_register_complete", "result": "pass", "detail": "d"}],
        "attested_by": {"name": "n", "role": "r", "statement": "s"},
    }
    rendered = purge._render_markdown(payload)
    assert "C1" in rendered and "absent" in rendered
    purge._assert_carries_no_patient_data(rendered)


# ------------------------------------------------------------------ kept placeholders
#
# The run of 2026-11-21 emptied results/ correctly, kept its .gitkeep as designed, and then
# refused to attest because the measurement counted that .gitkeep as data. It also removed
# notebooks/ whole, taking the tracked placeholder with it.


def _in_repo(tmp_path, monkeypatch, name):
    repo = tmp_path / "repo"
    target = repo / name
    target.mkdir(parents=True)
    (target / ".gitkeep").write_text("", encoding="utf-8")
    monkeypatch.setattr(purge, "REPO_ROOT", repo)
    return target


@pytest.mark.parametrize("name", ["results", "notebooks"])
def test_an_emptied_repo_location_with_its_placeholder_is_clean(tmp_path, monkeypatch, name):
    target = _in_repo(tmp_path, monkeypatch, name)
    (target / "out").mkdir()
    (target / "out" / "f.json").write_text("{}", encoding="utf-8")
    location = _location(target)

    purge._delete(location)
    assert [p.name for p in target.iterdir()] == [".gitkeep"]
    after = purge.measure(location)
    assert (after.exists, after.file_count, after.total_bytes) == (True, 0, 0)


def test_a_placeholder_that_holds_bytes_is_still_counted(tmp_path, monkeypatch):
    target = _in_repo(tmp_path, monkeypatch, "results")
    (target / ".gitkeep").write_bytes(b"x" * 5)
    state = purge.measure(_location(target))
    assert (state.file_count, state.total_bytes) == (1, 5)


def test_a_placeholder_outside_the_repo_is_counted(tmp_path):
    target = tmp_path / "custody" / "raw"
    target.mkdir(parents=True)
    (target / ".gitkeep").write_text("", encoding="utf-8")
    assert purge.measure(_location(target)).file_count == 1


def test_a_nested_placeholder_is_counted(tmp_path, monkeypatch):
    target = _in_repo(tmp_path, monkeypatch, "results")
    (target / "l0").mkdir()
    (target / "l0" / ".gitkeep").write_text("", encoding="utf-8")
    assert purge.measure(_location(target)).file_count == 1


def test_an_emptied_location_keeps_its_journalled_before_state(tmp_path, monkeypatch):
    """A re-run finds results/ present but empty; the journal's 242 files must survive."""
    target = _in_repo(tmp_path, monkeypatch, "results")
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "j.json")
    purge._write_journal({}, "C3", purge.LocationState(True, 242, 13_105_072))
    monkeypatch.setattr(purge, "read_register", lambda *a, **k: [_location(target, "C3")])
    monkeypatch.setattr(purge, "_discard_freed_blocks", lambda: None)

    verdict = purge.purge({}, confirm=True)[0]
    assert verdict.before.file_count == 242
    assert verdict.is_clean
    assert purge._read_journal()["C3"]["file_count"] == 242


def test_rejournalling_the_same_state_keeps_its_note(tmp_path, monkeypatch):
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "j.json")
    journal = {"C1": {"exists": True, "file_count": 13, "total_bytes": 1, "note": "recovered"}}
    purge._write_journal(journal, "C1", purge.LocationState(True, 13, 1))
    assert purge._read_journal()["C1"]["note"] == "recovered"


def test_a_changed_state_drops_the_stale_note(tmp_path, monkeypatch):
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "j.json")
    journal = {"C1": {"exists": True, "file_count": 13, "total_bytes": 1, "note": "recovered"}}
    purge._write_journal(journal, "C1", purge.LocationState(True, 14, 2))
    assert "note" not in purge._read_journal()["C1"]


def test_a_refilled_location_keeps_what_the_earlier_run_found(tmp_path, monkeypatch):
    """C7 held 13 files when the first run emptied it; this session then wrote 1 more."""
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "j.json")
    journal = {}
    purge._write_journal(journal, "C7", purge.LocationState(True, 13, 63_625_301))
    purge._write_journal(journal, "C7", purge.LocationState(True, 1, 600_000))

    entry = purge._read_journal()["C7"]
    assert entry["file_count"] == 1
    assert entry["earlier"]["file_count"] == 13

    verdict = purge.LocationVerdict(
        location=_location(tmp_path / "c7", "C7"),
        before=purge.LocationState(True, 1, 600_000),
        after=purge.LocationState(False, 0, 0),
        verified_utc="2026-11-21T00:00:00Z")
    row = purge._location_entry(verdict, entry)
    assert row["earlier_before"]["total_bytes"] == 63_625_301
    assert "note" in row


def test_an_unchanged_state_keeps_its_earlier_record(tmp_path, monkeypatch):
    monkeypatch.setattr(purge, "JOURNAL", tmp_path / "j.json")
    journal = {"C7": {"exists": True, "file_count": 1, "total_bytes": 2,
                      "earlier": {"exists": True, "file_count": 13, "total_bytes": 9}}}
    purge._write_journal(journal, "C7", purge.LocationState(True, 1, 2))
    assert purge._read_journal()["C7"]["earlier"]["file_count"] == 13
