# Compliance checklist
- No patient data / identifiable intermediates in git (see .gitignore). No Git LFS for data.
- No recontact of the subject, family, or MVA Society contacts.
- Delete all data within 30 days of Hackathon close; email confirmation to
  RarediseaserealkidMVAhackathon2026@synapse.org.
- Manuscript embargo until organizers post their summary/preprint; code/outputs shareable anytime.
- **Required acknowledgement — reproduce verbatim.** The rules mandate this exact text in any
  publication, preprint, conference abstract, or public communication arising from participation:

  > "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership
  > with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment
  > Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to
  > the child and their family who generously contributed their data and their story to advance
  > research into this rare disease. We acknowledge their trust in making this Hackathon possible."

  The dataset must also be cited using the reference on the Hackathon Synapse page at the time of
  publication.
- **Third-party LLM services are tools, not recipients** (decision D17). A service may be used only
  where it does not train on inputs or outputs, takes no rights in either, and retains nothing beyond
  short-lived operational logs. The pipeline's reasoning step runs a **local** model, so no candidate
  data leaves the machine. Development assistance is disclosed — provider, plan and the training
  setting — in the methods description form. Assistant session transcripts are patient-derived
  custody (`C7`) and are purged on the same deadline as the dataset.
- Re-identification-avoidance: publish nothing that could re-identify the child/family.
- **The dataset's own phenotype document sets a specific boundary** (verified 2026-09-07): no
  publication or communication arising from this challenge may include information that could
  re-identify the patient or family **beyond what the family already shares publicly through their
  own blog posts**. That is narrower than a general privacy rule and is the operative test — the
  family controls what is public about them, and this project must not widen it.
- Clinical phenotype (including HPO term sets) is patient data. It is parsed from `data_dir` at
  runtime and **never** committed — not to config, not to code, not to a fixture. A specific
  combination of features is identifying in a population of roughly 50 patients worldwide.
- Outputs CC-BY-4.0; avoid NC / ShareAlike sources in redistributed outputs.
- Data custody is registered at creation, not reconstructed at deletion. Every location holding
  patient-derived bytes is listed in [docs/data_custody.md](docs/data_custody.md) the moment it is
  created — including derived intermediates under `results/`, and the Hugging Face cache, which is
  why `HF_HOME` is set inside the custody root before the first download.
- Deletion is evidenced by an attestation from [src/purge.py](src/purge.py), recording per-location
  before/after file counts and byte totals plus verification checks. It identifies the destroyed
  dataset by its Hub revision SHA only — never filenames, never per-file digests, since both can
  aid re-identification.
- Attach that attestation to the deletion-confirmation email to
  RarediseaserealkidMVAhackathon2026@synapse.org.
