"""L2 Channel D -- Phenotype / HPO-driven.

Purpose
    Nominate candidates from the *phenotype* rather than the mechanism: find diseases whose
    annotated features resemble the subject's HPO-coded phenotype, then the approved drugs
    indicated for those diseases or for the subject's features themselves.

    This is the channel that reaches **symptomatic** candidates -- the ones a purely
    mechanism-first pipeline structurally cannot see.

Inputs
    Phenotype terms from the dataset, already HPO-coded (verified 2026-09-07), parsed
    from ``config["data_dir"]`` at runtime (:mod:`src.l2_channels.phenotype`); the Human
    Phenotype Ontology; Monarch's disease-to-phenotype associations; Open Targets'
    clinical indications and drug identities (restricted zone, :mod:`.enrichment`).

    **Read the terms; do not assume them.** This proband's presentation is not the
    textbook MVA description, and an earlier version of this docstring listed the
    literature phenotype (microcephaly, developmental delay, seizures) as though it were
    the subject's. Building against remembered features rather than the supplied ones
    would generate candidates for manifestations this child does not have, while missing
    the ones documented. The dataset is the only admissible source for which terms apply.

Outputs
    Under ``results_dir/l2/channel_d_phenotype/``:
      - ``candidates.tsv`` -- ranked drugs; identity fields through the enrichment
        whitelist, a score, the route that produced it, support counts, and
        ``endpoint=symptomatic`` on every row. No HPO id, no disease.
      - ``evidence.json`` -- **patient-derived, never published**: per drug, the diseases and
        HPO terms that produced it; and the disease ranking. Kept under the gitignored
        ``results_dir`` because which diseases resemble this child, and which of the
        child's features a drug addresses, is a statement about the child.
      - ``channel.json`` -- method, parameters, counts, sources with the HPO version,
        licensing, caveats, seed. Counts only; no term.

Method
    1. The subject's terms are resolved to live HPO ids (alt ids and obsolete terms with a
       replacement are followed) and restricted to the phenotypic-abnormality subtree.
    2. Every Mondo disease in Monarch's association table is scored against them by
       one-sided semantic similarity with an empirical p-value (:mod:`.phenosim`)
       [kohler2009] doi:10.1016/j.ajhg.2009.09.003 [resnik1999] doi:10.1613/jair.514.
    3. Disease route: diseases with ``p <= p_max``, at most ``top_diseases`` of them, are
       joined to indications approved for exactly that Mondo disease. A drug's disease
       score is the best normalised similarity among its supporting diseases.
    4. Phenotype route: indications approved for an HPO term on the lineage of one of the
       subject's terms score :func:`.phenosim.term_ratio`, kept at ``>= min_term_ratio``.
    5. A drug's score is the higher of the two routes; ties break on support count, then
       ChEMBL id.

Guardrail
    Symptomatic relief is **not** disease modification, and the distinction must survive
    into the final report -- a drug nominated here addresses a manifestation, not the
    underlying chromosome missegregation. Every candidate carries
    ``endpoint="symptomatic"`` so L3 and L5 cannot present it as mechanistic.

    Phenotype terms are patient-derived: use HPO codes and aggregate counts only. Never
    write free-text clinical narrative, ages, dates, birth weights, or anything else that
    could contribute to re-identification of a child in a ~50-patient population.

    **The term set itself never enters the repository** -- not config, not source, not a
    test fixture, not a committed example. A specific combination of features is
    identifying even without a name attached, and the family's own publications set the
    boundary for what is public about them (``COMPLIANCE.md``). Parse from ``data_dir``
    on every run; a test that needs terms uses invented ones. Logs carry counts only.

    **Nothing about the subject is sent anywhere** (decision D10): the resources are
    downloaded whole and matched on this machine.

    **This channel annotates; it does not exclude.** Safety and paediatric suitability are
    L4's; Monarch's own terms say its data "should not be used for direct diagnostic use or
    medical decision-making", and nothing here is either.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from . import enrichment, phenosim, phenotype

LOGGER = logging.getLogger(__name__)

ENDPOINT = "symptomatic"
INDICATION_DATASETS = ("clinical_indication", "drug_molecule")

DEFAULTS = {"permutations": 1000, "p_max": 0.001, "top_diseases": 50,
            "min_term_ratio": 0.5, "indication_stage": enrichment.APPROVED,
            "approved_only": True, "open_targets_release": "26.06"}

CAVEATS = (
    "Symptomatic, not disease-modifying: a drug here is indicated for a disease that "
    "resembles the subject's presentation, or for one of its features. It says nothing "
    "about the underlying mitotic defect.",
    "Similarity measures overlap of curated annotations, which are uneven: well-studied "
    "diseases carry more of them. The p-value corrects for annotation volume against "
    "random queries; it does not correct for curation bias in which features were "
    "recorded at all.",
    "p-values are empirical and uncorrected for the number of diseases scored. The "
    "p_max cut with a top_diseases cap is a nomination threshold, not a significance "
    "claim; the smallest attainable p is 1/(permutations + 1).",
    "Indications are matched to the exact Mondo disease or to an HPO term on the lineage "
    "of a subject's term. An indication recorded against a broader or differently "
    "identified disease is missed, so absence here is absence of a matched indication, "
    "not absence of use.",
    "Only indications whose highest stage for that indication is approval are used. "
    "Approval is jurisdiction- and population-specific; paediatric approval is not "
    "implied, and is L4's to assess.",
    "The two routes' scores are not calibrated against each other. A disease-route score "
    "is a whole-phenotype similarity that passed a random-query null; a phenotype-route "
    "score is how closely one indication term matches one of the subject's terms, with no "
    "null. The best_route column says which produced each rank.",
    "Coverage is thin by construction: most diseases that resemble a rare presentation "
    "are rare themselves, and few carry an approved indication recorded against that exact "
    "disease. channel.json reports how many selected diseases contributed any drug.",
    "Monarch's data is provided as is, without warranty, and is not to be used for "
    "diagnosis or medical decision-making (Monarch terms of use).",
)


def _settings(config: dict) -> dict:
    return {**DEFAULTS, **((config.get("l2") or {}).get("channel_d") or {})}


def _write_tsv(path: Path, header, rows) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def query_terms(raw_terms: list, obo: dict, onto: phenosim.Ontology) -> tuple:
    """Resolve the subject's terms to ontology indices. Returns ``(indices, counts)``.

    Counts only are returned alongside -- never which term failed, since that is itself a
    fact about the subject.
    """
    root = onto.index.get(phenotype.PHENOTYPIC_ABNORMALITY)
    if root is None:
        raise ValueError(f"the ontology lacks {phenotype.PHENOTYPIC_ABNORMALITY}")
    indices, unresolved, outside = [], 0, 0
    for term in raw_terms:
        live = phenotype.resolve(term, obo)
        if live is None or live not in onto.index:
            unresolved += 1
            continue
        idx = onto.index[live]
        if not onto.is_under(idx, root):
            outside += 1
            continue
        indices.append(idx)
    indices = sorted(set(indices))
    return indices, {"terms_in_document": len(raw_terms), "terms_used": len(indices),
                     "terms_unresolved": unresolved, "terms_outside_phenotypic_abnormality":
                     outside}


def build_corpus(annotations: dict, obo: dict, onto: phenosim.Ontology) -> tuple:
    """Map Monarch's annotations onto ontology indices, within the abnormality subtree."""
    root = onto.index[phenotype.PHENOTYPIC_ABNORMALITY]
    under = set(onto.descendants[root].tolist())
    mapped, dropped = {}, 0
    for disease, terms in annotations.items():
        keep = set()
        for term in terms:
            live = phenotype.resolve(term, obo)
            idx = onto.index.get(live) if live else None
            if idx is None or idx not in under:
                dropped += 1
                continue
            keep.add(idx)
        if keep:
            mapped[disease] = keep
    return phenosim.Corpus.build(onto, mapped), {"annotations_dropped": dropped,
                                                 "diseases_with_annotations": len(mapped)}


def rank_drugs(disease_hits: list, term_hits: dict, indications: dict, molecules: dict, *,
               approved_only: bool) -> list:
    """Combine both routes into ranked drug entries.

    Args:
        disease_hits: selected ``DiseaseScore`` objects (Mondo ids, ``MONDO:`` form).
        term_hits: ``{HP id: ratio}`` for indication terms on the subject's lineage.
        indications: ``{Open Targets disease id: {chembl_id}}``.
        molecules: ``{chembl_id: identity}``.
    """
    support: dict = {}
    for hit in disease_hits:
        for drug in indications.get(hit.disease.replace(":", "_"), ()):
            support.setdefault(drug, {"diseases": [], "terms": []})["diseases"].append(hit)
    for term, ratio in term_hits.items():
        for drug in indications.get(term.replace(":", "_"), ()):
            support.setdefault(drug, {"diseases": [], "terms": []})["terms"].append(
                (term, ratio))

    entries = []
    for drug, found in support.items():
        identity = molecules.get(drug, {})
        record = enrichment.DrugRecord(chembl_id=drug, **{k: identity.get(k, "") for k in
                                                          ("name", "drug_type",
                                                           "clinical_stage")})
        if approved_only and not record.approved:
            continue
        disease_score = max((h.score for h in found["diseases"]), default=0.0)
        term_score = max((r for _, r in found["terms"]), default=0.0)
        entries.append({
            "record": record,
            "score": max(disease_score, term_score),
            "route": "disease" if disease_score >= term_score else "phenotype",
            "diseases": sorted(found["diseases"], key=lambda h: (h.p_value, -h.score,
                                                                 h.disease)),
            "terms": sorted(found["terms"], key=lambda t: (-t[1], t[0])),
        })
    entries.sort(key=lambda e: (-e["score"],
                                -(len(e["diseases"]) + len(e["terms"])),
                                e["record"].chembl_id))
    return entries


def generate(config: dict) -> None:
    """Generate symptomatic candidates from the subject's phenotype terms.

    Args:
        config: Parsed pipeline configuration; uses ``data_dir`` (phenotype file),
            ``reference_dir``, ``enrichment_dir``, ``seed``, ``results_dir`` and
            ``l2.channel_d``.

    Raises:
        ValueError / FileNotFoundError: On an unreadable phenotype, an empty corpus, or no
            candidate -- each a broken input or an unmet threshold, reported as such rather
            than as an empty table (see ``l2_channels/run.py``).
    """
    settings = _settings(config)
    min_p = 1 / (int(settings["permutations"]) + 1)
    if float(settings["p_max"]) < min_p:
        # Otherwise no disease can ever pass, and the run would fail later as an "unmet
        # threshold" -- a configuration error dressed as a result.
        raise ValueError(f"l2.channel_d.p_max={settings['p_max']} is below the smallest "
                         f"attainable p-value 1/(permutations + 1) = {min_p:.6f}; raise "
                         "permutations or p_max")
    out_dir = Path(config["results_dir"]) / "l2" / "channel_d_phenotype"
    out_dir.mkdir(parents=True, exist_ok=True)

    # The fingerprint of a reviewed document comes from the environment, never from
    # config: it is a handle on a patient file and must not be committed.
    raw_terms, document_counts = phenotype.read_terms(
        phenotype.find_document(config),
        reviewed_fingerprint=os.environ.get("MVA_PHENOTYPE_REVIEWED", ""))

    paths, reference_provenance = phenotype.ensure_sources(config)
    with phenotype.open_text(paths["hpo"]) as handle:
        obo = phenotype.parse_obo(handle)
    onto = phenosim.Ontology.build(obo["parents"])
    with phenotype.open_text(paths["monarch"]) as handle:
        annotations, labels, association_stats = phenotype.load_associations(handle)

    query, query_counts = query_terms(raw_terms, obo, onto)
    if not query:
        raise ValueError(f"none of the {query_counts['terms_in_document']} phenotype term(s) "
                         "resolved to a live phenotypic-abnormality term in HPO "
                         f"{obo['version']}; check the ontology release")
    corpus, corpus_counts = build_corpus(annotations, obo, onto)
    LOGGER.info("scoring %d disease(s) against %d phenotype term(s), %d permutation(s)",
                len(corpus.diseases), len(query), settings["permutations"])

    scores, null_info = phenosim.score_diseases(onto, corpus, query, seed=config["seed"],
                                                permutations=int(settings["permutations"]))
    selected = [s for s in scores if s.p_value <= float(settings["p_max"])][
        :int(settings["top_diseases"])]

    ot_paths, drug_provenance = enrichment.ensure_datasets(
        config, INDICATION_DATASETS, channel="channel_d")
    indications, indication_stats = enrichment.load_indications(
        ot_paths["clinical_indication"], stage=str(settings["indication_stage"]))
    molecules = enrichment.load_molecules(ot_paths["drug_molecule"])

    term_hits = {}
    for disease_id in indications:
        if not disease_id.startswith("HP_"):
            continue
        term = disease_id.replace("_", ":")
        idx = onto.index.get(phenotype.resolve(term, obo) or "")
        if idx is None:
            continue
        ratio = phenosim.term_ratio(onto, corpus, query, idx)
        if ratio >= float(settings["min_term_ratio"]):
            term_hits[term] = ratio

    entries = rank_drugs(selected, term_hits, indications, molecules,
                         approved_only=bool(settings["approved_only"]))
    if not entries:
        raise ValueError(
            f"no approved drug reached channel D: {len(selected)} disease(s) passed "
            f"p <= {settings['p_max']} and {len(term_hits)} indication term(s) matched the "
            "phenotype's lineage. That is an unmet threshold on these inputs, not a finding "
            "that no drug addresses the phenotype; see channel.json.")

    rows = []
    for position, entry in enumerate(entries, start=1):
        public = enrichment.publishable(entry["record"])
        rows.append([position, public["chembl_id"], public["name"], public["drug_type"],
                     public["clinical_stage"], f"{entry['score']:.4f}", entry["route"],
                     len(entry["diseases"]), len(entry["terms"]), ENDPOINT])
    _write_tsv(out_dir / "candidates.tsv",
               ["rank", "chembl_id", "drug_name", "drug_type", "clinical_stage", "score",
                "best_route", "n_supporting_diseases", "n_supporting_phenotype_terms",
                "endpoint"], rows)

    (out_dir / "evidence.json").write_text(json.dumps({
        "patient_derived": True,
        "redistributable": False,
        "notice": ("Derived from the subject's phenotype and from licence-restricted "
                   "indication data. Never publish, quote or commit this file; report "
                   "counts and candidate identities only."),
        "drugs": [{
            "chembl_id": e["record"].chembl_id,
            "diseases": [{"mondo_id": h.disease, "label": labels.get(h.disease, ""),
                          "score": round(h.score, 4), "p_value": h.p_value}
                         for h in e["diseases"]],
            "phenotype_terms": [{"hpo_id": t, "ratio": round(r, 4)} for t, r in e["terms"]],
        } for e in entries],
        "query_self_similarity": null_info["self_similarity"],
        "diseases_selected": [{"mondo_id": s.disease, "label": labels.get(s.disease, ""),
                               "score": round(s.score, 4), "p_value": s.p_value}
                              for s in selected],
        "seed": config["seed"],
    }, indent=2, sort_keys=True), encoding="utf-8")

    (out_dir / "channel.json").write_text(json.dumps({
        "channel": "channel_d_phenotype",
        "endpoint": ENDPOINT,
        "method": {
            "similarity": "one-sided semantic similarity, query -> disease: mean over query "
                          "terms of the best Resnik similarity to the disease's annotations, "
                          "normalised by the query's self-similarity [resnik1999] "
                          "doi:10.1613/jair.514 [kohler2009] doi:10.1016/j.ajhg.2009.09.003",
            "information_content": "-ln(fraction of diseases annotated with the term or a "
                                   "descendant); unannotated terms floored at one disease",
            "null": "random queries of equal size drawn without replacement from terms used "
                    "in the corpus; p = (1 + exceedances) / (permutations + 1)",
            "disease_route": "approved indications for exactly the selected Mondo disease",
            "phenotype_route": "approved indications for an HPO term on a query term's "
                               "lineage, scored IC(MICA) / IC(query term)",
            "ranking": "descending score, then support count, then ChEMBL id",
        },
        "parameters": {k: settings[k] for k in DEFAULTS},
        "counts": {**query_counts, **document_counts, **corpus_counts, **association_stats,
                   "diseases_scored": len(corpus.diseases),
                   "diseases_selected": len(selected),
                   "selected_diseases_contributing_a_drug": len(
                       {h.disease for e in entries for h in e["diseases"]}),
                   "indication_terms_matched": len(term_hits),
                   "drugs_ranked": len(entries), **indication_stats},
        # self_similarity is the mean information content of the subject's own terms:
        # a number someone holding the public files could recompute for a guessed feature
        # set and match. It belongs with the other patient-derived values, not here.
        "null": {k: v for k, v in null_info.items() if k != "self_similarity"},
        "sources": {"reference": reference_provenance, "indications": drug_provenance,
                    "hpo_version": obo["version"]},
        "attribution": (f"This channel uses the Human Phenotype Ontology (version "
                        f"{obo['version']}) and the Monarch Initiative knowledge graph "
                        f"(release {reference_provenance.get('monarch', {}).get('release', 'unknown')}"
                        "); disease labels are Mondo (CC BY 4.0)."),
        "licensing": {
            "published_fields": list(enrichment.PUBLISHABLE_FIELDS),
            "note": ("Indication pairs are ChEMBL-derived and stay in the enrichment zone "
                     "(decision D7); candidates.tsv carries only whitelisted identity fields "
                     "and counts. evidence.json is patient-derived and never published. "
                     "See src/l4_validate/sources.md, Table 6."),
        },
        "caveats": list(CAVEATS),
        "seed": config["seed"],
    }, indent=2, sort_keys=True), encoding="utf-8")

    LOGGER.info("channel D written: %d drug(s) from %d disease(s) and %d indication term(s)",
                len(entries), len(selected), len(term_hits))
