# Compliance checklist
- No patient data / identifiable intermediates in git (see .gitignore). No Git LFS for data.
- No recontact of the subject, family, or MVA Society contacts.
- Delete all data within 30 days of Hackathon close; email confirmation to
  RarediseaserealkidMVAhackathon2026@synapse.org.
- Manuscript embargo until organizers post their summary/preprint; code/outputs shareable anytime.
- Required attribution block (Sage Bionetworks, MVA Society, Hugging Face, BEACON, AWS, Anthropic, family).
- Re-identification-avoidance: publish nothing that could re-identify the child/family.
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
