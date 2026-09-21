"""Smoke tests: every layer imports, every ``run()`` is a callable stub.

These do not test behavior -- there is none yet. They guard the scaffold's two
invariants: the package tree is importable, and the config's two gate decisions
(``causal_gene`` at G1, ``therapeutic_endpoint`` at G2) hold values a person recorded,
never values someone quietly filled in.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

LAYER_MODULES = [
    "src.l0_genomics.run",
    "src.l1_target.run",
    "src.l2_channels.run",
    "src.l3_integrate.run",
    "src.l4_validate.run",
    "src.l5_report.run",
]

CHANNEL_MODULES = [
    "src.l2_channels.channel_a_kg",
    "src.l2_channels.channel_b_proximity",
    "src.l2_channels.channel_c_signature",
    "src.l2_channels.channel_d_phenotype",
    "src.l2_channels.channel_e_prior",
]

HELPER_MODULES = [
    "src.l3_integrate.harmonize",
    "src.l3_integrate.aggregate",
    "src.l3_integrate.claude_reasoning",
    "src.l4_validate.safety_triage",
    "src.l4_validate.benchmark",
]


#: Layers still unimplemented end to end. L0, L1 and the L2 runner are absent: each is
#: implemented and needs real input, so an empty config fails on the missing input rather
#: than on an unimplemented body. Layers move off this list as they are built, one at a
#: time and deliberately. (The L2 *runner* is built; most channels behind it are not --
#: see IMPLEMENTED_CHANNELS.)
#: **Every layer is now implemented.** L5 was the last scaffold; it landed with the
#: dossier, so ``STUB_LAYER_MODULES`` is empty and the stub test below parametrises over
#: nothing. That is the intended end state, not a gap in coverage -- the rule the stub
#: test enforced (an unbuilt layer raises rather than returning empty) is now enforced for
#: every layer by ``test_implemented_layer_requires_real_config``.
IMPLEMENTED = {"src.l0_genomics.run", "src.l1_target.run", "src.l2_channels.run",
               "src.l3_integrate.run", "src.l4_validate.run", "src.l5_report.run"}
STUB_LAYER_MODULES = [m for m in LAYER_MODULES if m not in IMPLEMENTED]

#: Channels already built. A channel moves off the stub list only when it produces a real
#: ranking from real inputs -- the same one-at-a-time rule the layers follow.
#: Channel C is built and listed here even though ``channels.signature`` is false in the
#: shipped config: it produces a real ranking from real inputs and then declines to
#: nominate it, because the ranking does not beat an unrelated gene's knockdown (D19).
#: "Implemented" is about the code, not about whether the channel is enabled.
IMPLEMENTED_CHANNELS = {"src.l2_channels.channel_b_proximity",
                        "src.l2_channels.channel_c_signature",
                        "src.l2_channels.channel_d_phenotype",
                        "src.l2_channels.channel_e_prior"}
STUB_CHANNEL_MODULES = [m for m in CHANNEL_MODULES if m not in IMPLEMENTED_CHANNELS]


@pytest.mark.parametrize("module_name", STUB_LAYER_MODULES)
def test_layer_run_is_a_stub(module_name: str) -> None:
    """Each unimplemented layer exposes a callable ``run`` that raises.

    Parametrises over nothing now that L5 is built. Kept rather than deleted: the next
    layer added to ``LAYER_MODULES`` is covered by it the moment it appears, which is
    what made the rule hold as the stack was built one layer at a time.
    """
    module = importlib.import_module(module_name)
    assert callable(module.run)
    with pytest.raises(NotImplementedError):
        module.run({})


@pytest.mark.parametrize("module_name", sorted(IMPLEMENTED))
def test_implemented_layer_requires_real_config(module_name: str) -> None:
    """A built layer fails on its missing input, never on an unimplemented body.

    The counterpart to the stub assertion above, and what replaces it now that no
    scaffold is left: every layer must refuse an empty config rather than quietly
    producing an empty result, because an empty result is a scientific claim.
    """
    module = importlib.import_module(module_name)
    assert callable(module.run)
    with pytest.raises((KeyError, FileNotFoundError, ValueError)):
        module.run({})


def test_l0_is_partially_implemented() -> None:
    """L0 requires real config: it is no longer a stub that ignores its input.

    The burden half reads a VCF from ``data_dir``, so an empty config fails on the
    missing key rather than on an unimplemented body. Asserting this keeps the smoke
    suite honest about which layers are built.
    """
    # import_module, not `from src.l0_genomics import run` -- the package's __init__
    # re-exports the function under that name, so the latter binds the function.
    l0 = importlib.import_module("src.l0_genomics.run")

    assert callable(l0.run)
    with pytest.raises(KeyError):
        l0.run({})


def test_l2_runner_requires_a_channel_map() -> None:
    """The L2 runner refuses a config with no channel toggles instead of running nothing."""
    l2 = importlib.import_module("src.l2_channels.run")

    assert callable(l2.run)
    with pytest.raises(ValueError, match="channels"):
        l2.run({})


def test_l1_requires_a_causal_gene() -> None:
    """L1 refuses to run without a gene rather than inventing one.

    ``causal_gene`` stays null until gate G1, and L0 writes only a *candidate*. A layer
    that quietly picked a gene -- from L0's artifact or anywhere else -- would fabricate
    the scientific premise of the whole run.
    """
    l1 = importlib.import_module("src.l1_target.run")

    assert callable(l1.run)
    with pytest.raises(ValueError, match="gate G1"):
        l1.run({})


@pytest.mark.parametrize("module_name", STUB_CHANNEL_MODULES)
def test_channel_generate_is_a_stub(module_name: str) -> None:
    """Each unimplemented L2 channel exposes a callable ``generate`` that raises."""
    module = importlib.import_module(module_name)
    assert callable(module.generate)
    with pytest.raises(NotImplementedError):
        module.generate({})


@pytest.mark.parametrize("module_name", sorted(IMPLEMENTED_CHANNELS))
def test_implemented_channel_requires_real_config(module_name: str) -> None:
    """An implemented channel fails on its missing input, not on an unimplemented body.

    The stub assertion above is what proves an unbuilt channel cannot quietly return an
    empty list. This is its counterpart: a built channel must still refuse an empty
    config rather than inventing a module to score against.
    """
    module = importlib.import_module(module_name)
    assert callable(module.generate)
    with pytest.raises((KeyError, FileNotFoundError, ValueError)):
        module.generate({})


@pytest.mark.parametrize("module_name", HELPER_MODULES)
def test_helper_modules_import(module_name: str) -> None:
    """L3/L4 helper modules import cleanly."""
    assert importlib.import_module(module_name) is not None


def test_pipeline_imports_and_orchestrates() -> None:
    """The orchestrator imports and knows the layer order."""
    from src import pipeline

    assert pipeline.LAYER_ORDER == (
        "l0_genomics",
        "l1_target",
        "l2_channels",
        "l3_integrate",
        "l4_validate",
        "l5_report",
    )
    for name in ("load_config", "seed_everything", "run", "main"):
        assert callable(getattr(pipeline, name))


def test_channel_registry_covers_config_channels() -> None:
    """Every channel toggle in the config maps to a registered channel module."""
    from src.l2_channels import CHANNEL_REGISTRY

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    assert set(config["channels"]) == set(CHANNEL_REGISTRY)


def test_causal_gene_is_a_recorded_finding() -> None:
    """``causal_gene`` holds a gene a person decided at G1, traceable to L0's evidence.

    This replaces an assertion that it was ``None``. Gate G1 passed on 2026-09-17 (D9),
    so "still null" is no longer the invariant -- but the reason for the original test
    has not gone away: a gene nobody recorded is indistinguishable from a gene someone
    guessed. So the value must be a panel gene, must be named in a G1 decision, and --
    where L0's artifact exists on this machine -- must be one of L0's candidates. The
    artifact is gitignored, so on a fresh clone that last check is skipped, not faked.
    """
    import json

    from src.l0_genomics.run import SAC_PANEL

    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    gene = config["causal_gene"]
    assert config["seed"] == 42

    assert gene in SAC_PANEL, f"causal_gene {gene!r} is not a panel gene"
    decisions = (REPO_ROOT / "mngmt" / "decisions.md").read_text(encoding="utf-8")
    assert "Gate G1" in decisions and f"`causal_gene: {gene}`" in decisions, (
        f"causal_gene is {gene!r} but no G1 decision in mngmt/decisions.md records it"
    )

    call = REPO_ROOT / config.get("results_dir", "results") / "l0_genomics" / "causal_gene_call.json"
    if not call.exists():
        pytest.skip("L0 has not been run on this machine; the candidate check needs its artifact")
    candidates = json.loads(call.read_text(encoding="utf-8"))["candidate_genes"]
    assert gene in candidates, (
        f"causal_gene {gene!r} is not among L0's current candidates {candidates} -- "
        "G1 is reopened, not silently updated (D9)"
    )


def test_therapeutic_endpoint_is_decided_and_documented() -> None:
    """``therapeutic_endpoint`` holds a valid value, recorded in the decision log.

    This replaces an assertion that it was ``None``. Gate G2 passed on 2026-09-08, so
    "still null" is no longer the invariant -- but the reason the original test existed
    has not gone away, so it is tightened rather than removed: the value must be one of
    the three admissible endpoints, and it must be traceable to a written decision.
    An endpoint nobody recorded is indistinguishable from an endpoint someone guessed.
    """
    config = yaml.safe_load((REPO_ROOT / "config" / "pipeline.yaml").read_text())
    endpoint = config["therapeutic_endpoint"]

    assert endpoint in {"chemoprevention", "symptomatic", "mitotic_fidelity"}, (
        f"{endpoint!r} is not an admissible therapeutic endpoint"
    )

    decisions = (REPO_ROOT / "mngmt" / "decisions.md").read_text(encoding="utf-8")
    assert endpoint in decisions, (
        f"therapeutic_endpoint is {endpoint!r} but mngmt/decisions.md does not record it"
    )


def test_no_patient_data_files_present() -> None:
    """No genomic file may exist anywhere in the repo -- including outside data/."""
    patterns = ("*.vcf", "*.vcf.gz", "*.bam", "*.cram", "*.fastq", "*.fastq.gz", "*.fq")
    offenders = [
        path
        for pattern in patterns
        for path in REPO_ROOT.rglob(pattern)
        if ".git" not in path.parts
    ]
    assert offenders == [], f"genomic files found in repo: {offenders}"


#: Constraint 7 forbids writing the subject's clinical features into any committed file.
#: These match the *form* of that leak -- an assertion about this individual -- rather than
#: a clinical vocabulary, because the vocabulary is unbounded and the repo legitimately
#: discusses MVA's textbook features. Attribution is what makes a feature identifying.
#:
#: The first version keyed on the word "proband" alone and missed the worst instance in the
#: repository, which named the subject through the phenotype document instead. Both forms
#: are matched now.
PHENOTYPE_ATTRIBUTION = (
    (r"\bproband'?s?\s+(?:has|had|presents|shows|exhibits|lacks)\b",
     "a clinical assertion about the proband"),
    (r"\bproband'?s\s+(?:actual\s+|documented\s+)?(?:HPO|phenotype)\b",
     "a reference to the proband's phenotype"),
    (r"\b(?:documented|supplied|subject'?s?)\s+phenotype\s+"
     r"(?:has|had|shows|includes|contains|lacks|records)\b",
     "a characterisation of the documented phenotype"),
    (r"\bphenotype\s+document\s+records\b",
     "a clinical fact attributed to the phenotype document"),
)

#: HPO ids are checked only outside ``tests/``. Tests are required to use invented terms
#: (CLAUDE.md constraint 7) and two real ids appear there deliberately, as the fixtures
#: proving ``literature.refuse_private`` blocks an HPO id from reaching Europe PMC. That
#: exemption is a real gap: this guard cannot tell an invented id from a subject's own.
#: What it does cover is the surface the leak actually used -- prose in docs and source.
HPO_ID = r"\bHP:\d{7}\b"

#: Ontology structural terms, not clinical features: "All", "Mode of inheritance",
#: "Phenotypic abnormality". Channel D needs these to walk the ontology.
HPO_STRUCTURAL = frozenset({"HP:0000001", "HP:0000005", "HP:0000118"})

#: The one sanctioned statement of the subject's clinical basis, kept by the maintainer's
#: decision because D4 is unintelligible without its premise and the fact is close to the
#: base rate for this genotype (roughly 75% of BUB1B MVA patients develop cancer). It is
#: matched against the whole line so a second occurrence elsewhere in the same file still
#: fails; the delivery plan used to restate it and now points at D4 instead.
SANCTIONED = ((
    "mngmt/decisions.md",
    "The proband has already had one: the phenotype document records a malignancy",
),)

#: This file necessarily contains the patterns themselves, so it cannot be scanned with
#: them. Nothing else is exempt.
GUARD_FILE = "tests/test_smoke.py"


def test_no_subject_phenotype_in_tracked_files() -> None:
    """Constraint 7: no clinical feature of the subject in any committed file.

    A specific combination of features identifies in a population of roughly fifty, and
    this repository becomes public before judging. Vigilance already failed here twice --
    once when a term reached a committed file despite the rule being written in three
    places, and again when the sweep for it keyed on the wrong word -- so the rule is
    enforced rather than remembered.
    """
    import re
    import subprocess

    tracked = subprocess.run(
        ["git", "ls-files"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()

    offenders = []
    for rel in tracked:
        if rel == GUARD_FILE:
            continue
        path = REPO_ROOT / rel
        try:
            content = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue

        checks = list(PHENOTYPE_ATTRIBUTION)
        if not rel.startswith("tests/"):
            checks.append((HPO_ID, "an HPO term id"))

        for pattern, description in checks:
            for found in re.finditer(pattern, content):
                if found.group(0) in HPO_STRUCTURAL:
                    continue
                start = content.rfind("\n", 0, found.start()) + 1
                end = content.find("\n", found.end())
                whole_line = content[start: end if end != -1 else len(content)]
                if any(rel == f and s in whole_line for f, s in SANCTIONED):
                    continue
                line = content.count("\n", 0, found.start()) + 1
                offenders.append(f"{rel}:{line} contains {description}")

    assert offenders == [], (
        "subject clinical data in committed files (CLAUDE.md constraint 7):\n  "
        + "\n  ".join(offenders)
    )
