# Track 2 candidate dossier

The generated output of layer L5, committed here so the Track 2 report's references to it
resolve for a reader who has only this repository. Open **[index.html](index.html)** — or,
if you are reading on GitHub, download the folder, since GitHub serves `.html` as source
rather than rendering it.

Read it in this order; it is arranged to be read as an argument, not browsed as a table.

| | |
|---|---|
| **[index.html](index.html)** | Start here. The two tiers, and why they are never merged into one ranked list |
| **[exclusions.html](exclusions.html)** | The headline result: 1,880 of 1,963 drugs removed, each with the label section and the sentence that caused it |
| **[candidates/](candidates/)** | One page per survivor, 83 of them, with provenance and the adversarial pass against each |
| **[figures/](figures/)** | Four figures, including the channel-contribution plot showing that zero candidates are supported by two discriminating methods |
| **[tables/](tables/)** | The same data as TSV |

## What this is not

**It is not a treatment recommendation.** Exactly one candidate has any published evidence
behind it, that evidence is *in vitro* and about aneuploid cancer cell lines rather than
this disease, and the pipeline's own adversarial pass concludes it should not be
shortlisted. That objection is printed in full on its page rather than suppressed. The
other 82 rest on a network proximity score this pipeline computed, with nothing published
linking any of them to this disease, and they carry no generated rationale for that reason.

Nothing here has been tested in a laboratory or a clinic for this condition.

## Why it is safe to publish

Every string was written through `src/l5_report/publish_guard.py`, which fails closed on
anything that could identify the subject. The per-chromosome aneuploidy result is
deliberately absent — `report.json` lists the ten withheld fields by name, because in a
population of about fifty patients the pattern of involved chromosomes approaches a
fingerprint (decision D20). The four figures are drug counts and channel contributions;
none is genomic.

The contents were scanned again at the point of copying, for coordinates, HGVS
expressions, dbSNP and ClinVar identifiers, HPO terms, and the sample identifier. No
match.

## Regenerating it

```
PYTHONHASHSEED=42 python -m src.pipeline --only l5_report
```

Needs the gated dataset. After the custody purge (see
[../data_custody.md](../data_custody.md)) that is no longer possible from this machine,
which is why the output is committed rather than left to be rebuilt on demand.
