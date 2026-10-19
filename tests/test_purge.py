"""Custody register parsing and measurement.

These exercise the two read-only halves of ``src/purge.py``. Nothing here touches a real
custody location, and no fixture contains patient-derived content -- the register under
test is synthetic except where a case deliberately reads the committed one, which is read
and never written.
"""
from __future__ import annotations

import os
import textwrap

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
    rows = purge.read_register()
    assert [row.id for row in rows] == ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
    assert all(row.status in purge.CUSTODY_STATES or row.status in purge.UNPURGEABLE_STATES
               for row in rows)


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
