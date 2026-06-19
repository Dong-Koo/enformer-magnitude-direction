# Submission Checklist

Pre-submission checklist for:
**"Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction"**
DongKoo Lee · ceo@nrootm.com

Work through these steps in order. Each step's output feeds into the next.

---

## Phase 1: Identity and repository setup

- [ ] **ORCID** — Go to [orcid.org/register](https://orcid.org/register) and create a free account. Obtain your 16-digit ORCID iD (format: `0000-0000-0000-0000`). Add your affiliation ("Independent Researcher, Gwangju-si, Republic of Korea") and email (ceo@nrootm.com) to your ORCID profile.

- [ ] **Fill ORCID placeholder** — Replace `[ORCID-to-be-assigned]` and `[ORCID]` in all files:
  - `CITATION.cff` (field: `orcid`)
  - `.zenodo.json` (field: `creators[0].orcid`)
  - `docs/cover_letter.md` (signature block)
  - `docs/researcher_email.md` (all three email signatures)
  - `manuscript/manuscript.html` (author line)

---

## Phase 2: GitHub

- [ ] **Create GitHub account** if you do not have one: [github.com/join](https://github.com/join)

- [ ] **Create repository** — Suggested name: `eqtl-enformer-magnitude-direction`. Set visibility to **Public**. Do not initialize with a README (you already have one).

- [ ] **Upload all files** — Push or upload the full directory tree:
  - `README.md` (top-level)
  - `CITATION.cff`
  - `.zenodo.json`
  - `manuscript/manuscript.html`
  - `manuscript/figures/fig1_magnitude_direction.{png,pdf}`
  - `manuscript/figures/fig2_finetune_direction.{png,pdf}`
  - `manuscript/figures/fig3_model_comparison.{png,pdf}`
  - `code/01_download_data.py` through `05_make_figures.py`
  - `code/requirements.txt`, `code/environment.yml`
  - `data/processed/gene_level_results.csv`
  - `data/processed/classA_D10_eGenes.csv`, `classB_D1_eGenes.csv`
  - `docs/README.md`, `docs/cover_letter.md`, `docs/DATA_AVAILABILITY.md`
  - `docs/METHODS_REPRODUCIBILITY.md`, `docs/SUBMISSION_CHECKLIST.md`
  - `docs/researcher_email.md`, `docs/LICENSE`

- [ ] **Verify** — Confirm the repository is publicly viewable and the top-level `README.md` renders correctly on GitHub.

- [ ] **Copy GitHub URL** — Format: `https://github.com/[your-username]/eqtl-enformer-magnitude-direction`

- [ ] **Fill GitHub URL placeholder** — Replace `[GitHub URL]` and `[username]` in:
  - `README.md` (Quick start clone URL, citation block)
  - `CITATION.cff` (field: `repository-code`)
  - `docs/DATA_AVAILABILITY.md` (data statement paragraph)
  - `docs/researcher_email.md` (all three emails)
  - `manuscript/manuscript.html`

---

## Phase 3: Zenodo

- [ ] **Create Zenodo account** at [zenodo.org](https://zenodo.org) (free; sign in with GitHub recommended).

- [ ] **Connect GitHub** — In Zenodo: Settings → GitHub → flip the toggle ON for your `eqtl-enformer-magnitude-direction` repository.

- [ ] **Create a GitHub release** — On GitHub: Releases → Draft a new release → Tag: `v1.0.0` → Title: `v1.0.0 — initial release` → Publish release. Zenodo automatically imports `.zenodo.json` and mints a DOI.

- [ ] **Confirm Zenodo DOI** — Visit your Zenodo record and copy the DOI (format: `10.5281/zenodo.XXXXXXX`).

- [ ] **Fill Zenodo DOI placeholder** — Replace `[Zenodo DOI]` and badge URL in:
  - `README.md` (badge + data availability table)
  - `docs/DATA_AVAILABILITY.md` (data statement paragraph)
  - `docs/cover_letter.md` (reproducibility paragraph)
  - `docs/researcher_email.md` (all three emails)
  - `manuscript/manuscript.html`

---

## Phase 4: bioRxiv

- [ ] **Export manuscript to PDF** — Open `manuscript/manuscript.html` in a browser and print to PDF (File → Print → Save as PDF). Save as `manuscript_lee2026_enformer_direction.pdf`.

- [ ] **Submit to bioRxiv** — Go to [biorxiv.org/submit](https://www.biorxiv.org/submit):
  - Subject area: **Genomics**
  - Title: full paper title (copy from `README.md`)
  - Abstract: 2–3 sentence summary (copy from `README.md` Abstract section)
  - Authors: DongKoo Lee; affiliation: Independent Researcher, Gwangju-si, Republic of Korea; email: ceo@nrootm.com; ORCID: [your ORCID]
  - Upload PDF
  - Check the box: "This is not being submitted to a journal" (until you have a journal submission; update later if needed)

- [ ] **Confirm bioRxiv DOI** — bioRxiv usually assigns a DOI within 24–48 hours. DOI format: `10.1101/XXXXXXXX`. URL: `https://www.biorxiv.org/content/10.1101/XXXXXXXX`

- [ ] **Fill bioRxiv DOI placeholder** — Replace `[bioRxiv DOI]` in:
  - `README.md` (badge + citation block)
  - `CITATION.cff` (field: `doi`)
  - `.zenodo.json` (field: `related_identifiers` if you wish to add it)
  - `docs/cover_letter.md` (preprint line)
  - `docs/researcher_email.md` (all three emails)
  - `manuscript/manuscript.html`

- [ ] **Update Zenodo record** — If the bioRxiv DOI is now available, create a new Zenodo version (`v1.0.1`) with the updated files so the archive is current.

---

## Phase 5: Researcher outreach

Send only after the bioRxiv DOI is live and confirmed accessible.

- [ ] **Email Priority 1: Mostafavi Lab** (mostafavi.lab@ubc.ca, CC Alexander Sasse)
  — Template in `docs/researcher_email.md` → "Priority 1"
  — Subject: "Reanalysis of EnformerAssessment reveals magnitude–direction dissociation in individual-level predictions"

- [ ] **Email Priority 2: Enformer/Borzoi team** (Kundaje/Avsec labs)
  — Template in `docs/researcher_email.md` → "Priority 2"

- [ ] **Email Priority 3: PrediXcan/GTEx researchers** (Gamazon/GTEx leads)
  — Template in `docs/researcher_email.md` → "Priority 3"

- [ ] **Track responses** — Log any responses and note whether revision is needed. Do not commit to authorship changes based on email feedback alone.

---

## Phase 6: Journal submission

- [ ] **Final manuscript check** — Confirm all placeholders (`[GitHub URL]`, `[Zenodo DOI]`, `[bioRxiv DOI]`, `[ORCID]`) are replaced in the PDF. Confirm figure files are embedded.

- [ ] **Bioinformatics Advances** (1st choice, Oxford Academic)
  - Submission portal: [academic.oup.com/bioinformaticsadvances/pages/submission-guidelines](https://academic.oup.com/bioinformaticsadvances/pages/submission-guidelines)
  - Article type: Research Article
  - Upload: manuscript PDF, cover letter (`docs/cover_letter.md`), figures as separate files
  - Declare: no competing interests, no external funding, preprint posted on bioRxiv at [DOI]

- [ ] **If rejected → PLOS Computational Biology** (2nd choice)
  - Submission: [journals.plos.org/ploscompbiol/s/submit-now](https://journals.plos.org/ploscompbiol/s/submit-now)

- [ ] **If rejected → PLOS ONE** (3rd choice)

- [ ] **If rejected → GigaScience** (4th choice)

- [ ] **If rejected → BMC Genomics** (5th choice)

See `docs/cover_letter.md` for full rationale, reviewer suggestions, and reviewers to avoid.

---

## Quick reference: placeholders to replace

| Placeholder | Replace with | Files |
|-------------|-------------|-------|
| `[ORCID-to-be-assigned]` / `[ORCID]` | Your 16-digit ORCID iD | CITATION.cff, .zenodo.json, cover_letter.md, researcher_email.md, manuscript.html |
| `[username]` / `[GitHub URL]` | Your GitHub repo URL | README.md, CITATION.cff, DATA_AVAILABILITY.md, researcher_email.md, manuscript.html |
| `[Zenodo DOI]` | `10.5281/zenodo.XXXXXXX` | README.md, DATA_AVAILABILITY.md, cover_letter.md, researcher_email.md, manuscript.html |
| `[bioRxiv DOI]` | `10.1101/XXXXXXXX` | README.md, CITATION.cff, cover_letter.md, researcher_email.md, manuscript.html |
| `10.1101/XXXXXXXX` (in README badge) | Actual DOI digits | README.md badge URLs |
| `zenodo.XXXXXXX` (in README badge) | Actual Zenodo record number | README.md badge URLs |

---

> Last updated: June 2026
