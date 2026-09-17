"""L2 channel D -- phenotype semantic similarity: pure method, no I/O.

Purpose
    Rank diseases by how well their annotated features explain a query set of HPO terms,
    and attach an empirical p-value to each score, so that a disease is not ranked highly
    merely because it carries many annotations.

Inputs
    An ontology (term -> parents), disease -> annotated terms, a query term set, a seed.

Outputs
    :class:`Ontology`, :class:`Corpus` and :func:`score_diseases`' per-disease scores and
    p-values. Nothing is read from or written to disk here.

Method
    Similarity between two terms is the information content (IC) of their most
    informative common ancestor [resnik1999] doi:10.1613/jair.514, with IC estimated from
    how many diseases carry a term or one of its descendants. A query is compared with a
    disease one-sidedly: for each query term, the best match among the disease's terms;
    averaged over the query. Semantic similarity of this kind between a query and diseases
    annotated with the HPO, with p-values assigned to the scores, is the approach of
    Köhler et al. (2009) [kohler2009] doi:10.1016/j.ajhg.2009.09.003.

    The null is this module's own implementation, not a reproduction of theirs: random
    queries of the same size are drawn, without replacement, from the terms used in the
    annotation corpus; a disease's p-value is the fraction of random queries scoring at
    least as high against it, with the +1 correction so it is never zero.

Guardrail
    A similarity score is a statement about annotation overlap, not about aetiology or
    treatment response. This module ranks diseases; it never names a drug.

    Query terms are patient data when the caller is channel D. Nothing here logs, prints or
    stores them; functions return indices and numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

DEFAULT_PERMUTATIONS = 1000


@dataclass
class Ontology:
    """An is-a DAG over term ids, with ancestor and descendant closures (each incl. self)."""

    ids: list
    index: dict
    parents: list                      # per term: list of parent indices
    ancestors: list = field(default_factory=list)     # per term: np.ndarray (incl. self)
    descendants: list = field(default_factory=list)   # per term: np.ndarray (incl. self)

    @classmethod
    def build(cls, parent_map: dict) -> "Ontology":
        """``parent_map``: ``{term: [parent, ...]}``. Parents absent from the map are ignored."""
        ids = sorted(parent_map)
        index = {t: i for i, t in enumerate(ids)}
        parents = [[index[p] for p in parent_map[t] if p in index] for t in ids]
        onto = cls(ids=ids, index=index, parents=parents)
        onto._close()
        return onto

    def _close(self) -> None:
        n = len(self.ids)
        anc: list = [None] * n

        def visit(i: int) -> set:
            # Iterative, so a deep ontology cannot exhaust the recursion limit.
            stack = [(i, False)]
            while stack:
                node, expanded = stack.pop()
                if anc[node] is not None:
                    continue
                pending = [p for p in self.parents[node] if anc[p] is None]
                if expanded or not pending:
                    result = {node}
                    for p in self.parents[node]:
                        result |= anc[p]
                    anc[node] = result
                else:
                    stack.append((node, True))
                    stack.extend((p, False) for p in pending)
            return anc[i]

        for i in range(n):
            visit(i)
        desc: list = [[] for _ in range(n)]
        for i, a in enumerate(anc):
            for j in a:
                desc[j].append(i)
        self.ancestors = [np.fromiter(sorted(a), dtype=np.int64) for a in anc]
        self.descendants = [np.asarray(sorted(d), dtype=np.int64) for d in desc]

    def is_under(self, term: int, root: int) -> bool:
        return root in set(self.ancestors[term].tolist())


@dataclass
class Corpus:
    """Diseases and their directly annotated terms, with IC estimated over the corpus."""

    diseases: list                 # disease ids, in matrix order
    ptr: np.ndarray                # CSR row pointers into ``terms``
    terms: np.ndarray              # concatenated term indices
    ic: np.ndarray                 # per ontology term

    @classmethod
    def build(cls, onto: Ontology, annotations: dict) -> "Corpus":
        """``annotations``: ``{disease: iterable of term indices}``; empty diseases are dropped.

        IC(t) = -ln(n_t / N), where n_t counts diseases annotated with t or a descendant.
        A term no disease carries is given n_t = 1 -- the rarest observable -- rather than
        an infinite IC that would let one unannotated query term dominate every score.
        """
        diseases = sorted(d for d, ts in annotations.items() if ts)
        if not diseases:
            raise ValueError("the annotation corpus is empty")
        counts = np.zeros(len(onto.ids), dtype=np.int64)
        ptr, flat = [0], []
        for d in diseases:
            direct = sorted(set(annotations[d]))
            flat.extend(direct)
            ptr.append(len(flat))
            covered = set()
            for t in direct:
                covered.update(onto.ancestors[t].tolist())
            counts[list(covered)] += 1
        n = len(diseases)
        ic = -np.log(np.maximum(counts, 1) / n)
        return cls(diseases=diseases, ptr=np.asarray(ptr, dtype=np.int64),
                   terms=np.asarray(flat, dtype=np.int64), ic=ic)

    def annotated_terms(self) -> np.ndarray:
        return np.unique(self.terms)


def term_similarity(onto: Ontology, corpus: Corpus, query_term: int) -> np.ndarray:
    """IC of the most informative common ancestor of ``query_term`` and every term."""
    sim = np.zeros(len(onto.ids))
    for a in onto.ancestors[query_term]:
        idx = onto.descendants[a]
        sim[idx] = np.maximum(sim[idx], corpus.ic[a])
    return sim


def raw_scores(onto: Ontology, corpus: Corpus, query: list) -> np.ndarray:
    """Per disease: mean over query terms of the best-matching annotated term's similarity."""
    total = np.zeros(len(corpus.diseases))
    starts = corpus.ptr[:-1]
    for q in query:
        sim = term_similarity(onto, corpus, q)
        total += np.maximum.reduceat(sim[corpus.terms], starts)
    return total / len(query)


@dataclass(frozen=True)
class DiseaseScore:
    disease: str
    score: float        # raw score / the query's self-similarity, in [0, 1]
    raw: float
    p_value: float


def score_diseases(onto: Ontology, corpus: Corpus, query: list, *, seed: int,
                   permutations: int = DEFAULT_PERMUTATIONS, pool=None) -> tuple:
    """Score every disease against ``query`` and attach an empirical p-value.

    Args:
        query: Term indices (deduplicated by the caller). Order does not matter: it is
            sorted here, so the result cannot depend on how the terms were listed.
        pool: Term indices random queries are drawn from; defaults to every term used in
            the corpus.

    Returns:
        ``(scores, info)``: ``DiseaseScore`` per disease sorted by ascending p, then
        descending score, then disease id; and a dict of what the null used.
    """
    query = sorted(set(int(q) for q in query))
    if not query:
        raise ValueError("the query holds no usable term")
    self_similarity = float(np.mean(corpus.ic[query]))
    if self_similarity <= 0:
        raise ValueError("every query term has zero information content in this corpus")

    observed = raw_scores(onto, corpus, query)
    pool = np.asarray(sorted(set(int(t) for t in (corpus.annotated_terms() if pool is None
                                                  else pool))), dtype=np.int64)
    if len(pool) < len(query):
        raise ValueError(f"the random-query pool ({len(pool)} terms) is smaller than the "
                         f"query ({len(query)} terms)")
    rng = np.random.Generator(np.random.PCG64(seed))
    exceed = np.zeros(len(corpus.diseases), dtype=np.int64)
    for _ in range(permutations):
        sample = rng.choice(pool, size=len(query), replace=False)
        exceed += raw_scores(onto, corpus, sample.tolist()) >= observed - 1e-12

    p = (exceed + 1) / (permutations + 1)
    scores = [DiseaseScore(disease=d, score=float(observed[i] / self_similarity),
                           raw=float(observed[i]), p_value=float(p[i]))
              for i, d in enumerate(corpus.diseases)]
    scores.sort(key=lambda s: (s.p_value, -s.score, s.disease))
    info = {"permutations": permutations, "pool_size": int(len(pool)),
            "query_size": len(query), "self_similarity": self_similarity,
            "min_attainable_p": 1 / (permutations + 1)}
    return scores, info


def term_ratio(onto: Ontology, corpus: Corpus, query: list, term: int) -> float:
    """How well one term matches the query: max over query terms of sim(q, term) / IC(q).

    Only terms on a query term's own lineage count -- the query term itself, an ancestor,
    or a descendant. A sibling that merely shares a distant ancestor scores 0 here, since
    an indication for a different feature is not an indication for this one.
    """
    best = 0.0
    for q in query:
        if corpus.ic[q] <= 0:
            continue
        lineage = (term in set(onto.ancestors[q].tolist())
                   or q in set(onto.ancestors[term].tolist()))
        if not lineage:
            continue
        mica = max(corpus.ic[a] for a in set(onto.ancestors[q].tolist())
                   & set(onto.ancestors[term].tolist()))
        best = max(best, float(mica / corpus.ic[q]))
    return best
