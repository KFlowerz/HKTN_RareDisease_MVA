# mngmt — project management artifacts

Planning and tracking documents. Nothing here is part of the pipeline; nothing here may contain
patient-derived content.

| File | What it is |
|---|---|
| [delivery-plan.html](delivery-plan.html) | Scaffold-to-submission delivery plan: schedule matrix, phase timeline with gates G0–G4, risk register, cut line, and rubric map. Updated as gates pass and findings land. |

## Keeping it current

The plan is a live document, not a snapshot — it records gate outcomes and *realized* risks as they
happen, so the assumptions behind a decision stay legible weeks later. Update it when a gate passes,
a risk fires, or the schedule moves.

It is also published as a private artifact for viewing:
https://claude.ai/code/artifact/fb60f334-420f-4cf4-85a2-c287c5819922 — republishing from this path
updates that same URL rather than creating a second copy.

## A note on the HTML

The file is written as a document *fragment* — `<title>`, `<style>` and content, with no
`<!doctype>`, `<html>`, `<head>` or `<body>` wrapper, because the publishing step supplies those.
Browsers apply the title and styles correctly when the file is opened directly, so it is readable
straight from the repo; do not add the wrapper tags, as that would break publishing.
