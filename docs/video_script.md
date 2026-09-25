# The 3-minute video — production guide and script

**Rare Disease, Real Kid: MVA Hackathon 2026 — Track 2 deliverable**

Everything needed to record the submission video: the rules that constrain it, what to put
on screen, the words to say, and the checks to run before uploading.

The script is **410 words**. Measured against speaking rate: 2:38 at 155 words per minute,
2:47 at 147, 2:55 at 140 — and **3:09 at 130**, which is over. So read it at a normal
conversational pace or slightly brisker; a slow, deliberate delivery will not fit, and the
fix for that is cutting the 2:28 beat, not speeding up.

Read it aloud once with a timer before recording.

---

## 1. Before you record — the hard rules

These are not style preferences. Two of them can do real harm.

> **A screen recording is the easiest way to leak this patient's data.** Everything else in
> this project is protected by a gitignore, a publish guard or a path assertion. A video
> is protected by nothing except what you point the camera at.

**Never on screen, in any frame:**

| Do not show | Why |
|---|---|
| `data/` — any file in it | The patient's genome |
| `results/l0_genomics/*.json` | Contains variant coordinates |
| `results/submissions/track1/*.csv` | The Track 1 predictions — coordinates by design |
| `results/l0_genomics/aneuploidy_burden.json` | Per-chromosome result, withheld under decision D20 |
| Any variant coordinate, HGVS expression, rsID or ClinVar accession | Identifying in a population of about fifty |
| Any HPO code or the phenotype document | A specific combination of features identifies a person |
| An editor sidebar, terminal scrollback or tab bar showing the above | The leak is usually in the background, not the foreground |

**Safe to show, and what the video should use:**

- `results/l5/index.html` and the rest of the dossier — every string in it passes the
  publish guard
- `results/l5/figures/*.png`
- `docs/report_track2.md`, `docs/how_it_works.md`, `docs/glossary.md`
- The GitHub repository page
- `config/pipeline.yaml`

**Before you hit record:** close every editor tab and terminal that has touched `data/`,
`results/l0_genomics/` or `results/submissions/`. Start a clean window.

**Three things the narration must never do**, because they are wrong and because the
project has spent the whole build being careful about them:

1. **Never call `BUB1B` confirmed, or call this a diagnosis.** It is a research premise
   fixed by a person at a recorded gate.
2. **Always say *secondary* prevention**, never "prevention" alone. The unqualified word
   promises something much larger.
3. **Never claim a drug works, or might work, for this patient.** The output is a
   shortlist to investigate.

---

## 2. What to show, beat by beat

One shot per beat. Cuts on the beat boundaries.

| Time | On screen | Notes |
|---|---|---|
| 0:00–0:18 | Title card: project name, track, your name | Hold it. Do not animate |
| 0:18–0:33 | `docs/how_it_works.md` section 3 — the data-flow diagram | Scroll slowly through it once |
| 0:33–1:12 | `results/l5/exclusions.html`, scrolling | The refusal is the headline; give it the most screen time |
| 1:12–1:22 | `results/l5/figures/exclusions_by_reason.png` | Let the bar chart sit still |
| 1:22–1:45 | `results/l5/index.html` — the two tiers | Show tier 1 and tier 2 headings together |
| 1:45–2:00 | `results/l5/figures/channel_contribution.png` | The failed bet, shown not buried |
| 2:00–2:28 | `docs/report_track2.md` section 8 | The innovation and the withholding |
| 2:28–2:45 | `config/pipeline.yaml`, highlighting the gene line | The one-line change for the second disease |
| 2:45–3:00 | The GitHub repository page | End on something a judge can click |

---

## 3. The script

Read as written. Timings are cumulative. Numbers spoken as words are marked — say them
that way; the exact figure is on screen.

---

**[0:00 — The problem]**

> Mosaic variegated aneuploidy. About fifty patients worldwide. No therapy that touches the
> cause. The patient is a child whose disease already raises cancer risk — so a careless
> drug suggestion is not merely useless, it is harmful. I built a pipeline that runs from
> one whole-genome VCF to a shortlist you can audit, aimed at **secondary prevention**:
> reducing a further cancer in someone already at risk.

**[0:18 — What it is]**

> Six layers. Five independent search methods, deliberately unlike each other, on the bet
> that weak signals agreeing beat any single score. Config-driven, seeded, and every number
> names the file that produced it.

**[0:33 — The headline is a refusal]**

> The headline is a refusal. Nearly two thousand approved drugs went in. **Eighteen hundred
> and eighty** were removed by a paediatric safety triage. **Eighty-three** survived.
>
> Safety is a hard gate, never a weight — a ranked list reads as a recommendation whatever
> the caption says.
>
> And the largest bucket is the one I want checked. **Twelve hundred and thirty-one** of
> those mean "no safety information found", not "dangerous". The pipeline fails closed.
> Every exclusion records the label, the section, and the sentence that caused it.

**[1:22 — What survived, and what it rests on]**

> Of the eighty-three, exactly **one** has a published paper behind it. The other
> **eighty-two** rest on a number this pipeline computed — no published evidence links them
> to this disease. Two separate tiers, never one ranked table, because merging them invites
> you to read a distance calculation as a finding.

**[1:45 — The bet that failed]**

> The architecture rests on independent methods converging. **Zero** candidates are
> supported by two discriminating methods. I show that as a figure, not a footnote. One
> search, calibrated against **forty** unrelated genes, carried no gene-specific signal — it
> ships disabled and nominates nothing.

**[2:00 — Innovation, and what I withhold]**

> The dataset ships no alignments, so I quantify chromosome disruption from allele ratios in
> the called variants themselves — down to about **ten percent** of cells, measured from
> this sample's own depth.
>
> I publish the method, not the per-chromosome result. In a population of fifty, which
> chromosomes are involved is close to a fingerprint. That cost the project its most
> striking figure, and the withholding is named, not quietly dropped.

**[2:28 — Does it generalise]**

> One config line changed to cystic fibrosis. The module layer returned a different split —
> **fifty-four/one-forty-six** against **eighty-four/one-sixteen** — so it responds to the
> gene rather than reciting an answer. Two layers of six demonstrated; I claim no more.

**[2:45 — Close]**

> Hypothesis generation, for one patient. Nothing tested in a laboratory or a clinic. The
> causal gene is a research premise, not a diagnosis. Everything I have said names the
> artifact it came from.

---

## 4. Numbers spoken, and where they come from

Every figure in the script, with the artifact that produced it. `scripts/verify_report_claims.py`
checks this table against `results/`, so a re-run that moves a number fails the check rather
than leaving a wrong figure in a recording.

| Spoken | Exact | Artifact |
|---|---:|---|
| "nearly two thousand" | 1,963 | `results/l4/validation.json` · `counts.candidates` |
| "eighteen hundred and eighty" | 1,880 | `results/l4/validation.json` · `counts.excluded` |
| "eighty-three" | 83 | `results/l4/validation.json` · `counts.survivors` |
| "twelve hundred and thirty-one" | 1,231 | `counts.excluded_by_reason.insufficient_evidence` |
| "one" / "eighty-two" | 1 / 82 | `results/l5/report.json` · `counts` |
| "zero" | 0 | `results/l3/integration.json` · `counts.by_convergence.discriminating` |
| "forty" | 40 | `config/pipeline.yaml` · `l2.channel_c.n_control_genes` |
| "ten percent" | 0.098 | `results/l0_genomics/aneuploidy_burden.json` · `sensitivity` |
| "fifty-four/one-forty-six" | 54 / 146 | `results/scalability_cftr/l1_target/module.json` |
| "eighty-four/one-sixteen" | 84 / 116 | `results/l1_target/module.json` |

All at seed 42.

---

## 5. Recording notes

- **1080p, 16:9.** Screen recordings at native resolution — a scaled-down terminal is
  unreadable at a judge's playback size.
- **Record narration separately** from the screen capture if you can, then lay it under.
  Live narration while scrolling almost always runs long.
- **No music under the voice.** If you want a bed, keep it well under the narration and cut
  it entirely for the refusal section.
- **Burn in captions or ship an `.srt`.** Several numbers carry the argument, and a
  mis-heard "eighteen hundred" undoes the point.
- **Do not speed up the audio to fit.** If it runs over, cut the section 2:28 generalisation
  beat down to one sentence — it is the least load-bearing.
- **Under 3:00, not over.** Aim to land at 2:50.

---

## 6. Before you upload

- [ ] Watch it once at full screen, looking **only at the background** — sidebars, tab
      titles, notification pop-ups, terminal scrollback.
- [ ] Confirm no frame shows a coordinate, an HGVS expression, an rsID, a ClinVar
      accession, an HPO code, or a per-chromosome result.
- [ ] Confirm the narration never says "prevention" without "secondary", never calls the
      gene confirmed, and never says a drug would help.
- [ ] Run `python scripts/verify_report_claims.py` — if a number moved, the recording is
      wrong and needs the affected beat re-cut.
- [ ] Check the runtime is under 3:00.
- [ ] Record the AI disclosure in the methods form as well: the pipeline's local
      open-weights reasoning model, and the development assistant with training disabled
      (decision D17).

---

*Hypothesis generation only. Not medical advice, not a clinical recommendation, and not a
claim that any drug named here is safe or effective for any person.*
