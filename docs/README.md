# docs/ — Documentation Index

This directory contains supporting documentation for the manuscript and repository.

---

## File descriptions

| File | Description |
|------|-------------|
| `README.md` | This file. Index of docs/ contents and publication workflow. |
| `cover_letter.md` | Journal submission cover letter for Bioinformatics Advances (and fallback journals). Includes target journal rationale and suggested reviewers. |
| `DATA_AVAILABILITY.md` | Full data provenance statement for manuscript submission. Lists all input sources, processed output files, and confirms no controlled-access data are used. |
| `METHODS_REPRODUCIBILITY.md` | Step-by-step methods guide corresponding to each analysis script. Includes software versions, statistical decision rationale, and a copy-pasteable verification snippet. |
| `SUBMISSION_CHECKLIST.md` | Practical pre-submission checklist: ORCID → GitHub → Zenodo → bioRxiv → researcher emails → journal. |
| `researcher_email.md` | Ready-to-send outreach email templates for Mostafavi Lab, Enformer/Borzoi team, and PrediXcan/GTEx researchers. Send after bioRxiv DOI is assigned. |
| `LICENSE` | MIT license text covering code in this repository. Text, figures, and tables are CC BY 4.0. |

---

## Repository overview

**Paper:** "Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction: A Magnitude–Direction Dissociation in Individual-Level Expression Prediction"

**Author:** DongKoo Lee · ceo@nrootm.com · ORCID: [to be assigned]  
**Affiliation:** Independent Researcher, Gwangju-si, Republic of Korea

**Preprint:** [bioRxiv DOI — to be assigned]  
**Code archive:** [Zenodo DOI — to be assigned]

---

## Data sources

| Dataset | Source | Access |
|---------|--------|--------|
| Enformer per-gene Pearson R | [mostafavilabuw/EnformerAssessment](https://github.com/mostafavilabuw/EnformerAssessment) | Public |
| GTEx v8 Brain Cortex eGenes | [GTEx Portal API v2](https://gtexportal.org/api/v2) | Public |

No controlled-access or individual-level data are used. All inputs are publicly available summary statistics or aggregate model outputs. See `DATA_AVAILABILITY.md` for the full statement.

---

## Key results summary

| Finding | Value |
|---------|-------|
| eGene rate, lowest \|R\| decile | 31.6% |
| eGene rate, highest \|R\| decile | 82.5% |
| \|R\| model AIC advantage over R+R² | ΔAIC = 104 |
| Direction accuracy for eGenes | 54.0% (near-random) |
| PrediXcanR for R<0 eGenes | 0.295 (correct direction) |

---

## Publication workflow

Steps in order. Each step depends on the previous.

1. **ORCID** — Create a free researcher ID at [orcid.org](https://orcid.org). Obtain your 16-digit ORCID iD (format: 0000-0000-0000-0000). Fill all `[ORCID-to-be-assigned]` placeholders in `CITATION.cff`, `.zenodo.json`, `docs/cover_letter.md`, `docs/researcher_email.md`, and `manuscript/manuscript.html`.

2. **GitHub** — Create a public repository (suggested name: `eqtl-enformer-magnitude-direction`). Upload all files from this package. Ensure the top-level `README.md`, `CITATION.cff`, and `.zenodo.json` are present at the repo root.

3. **Zenodo** — Connect Zenodo to your GitHub account (zenodo.org → Settings → GitHub). Enable the new repository. Create a GitHub release (e.g., `v1.0.0`). Zenodo auto-imports `.zenodo.json` and mints a DOI. Fill all `[Zenodo DOI]` placeholders.

4. **bioRxiv** — Upload `manuscript/manuscript.html` (printed to PDF) at [biorxiv.org/submit](https://www.biorxiv.org/submit). Fill metadata (title, abstract, author, subject: Genomics). bioRxiv assigns a DOI within ~24 hours. Fill all `[bioRxiv DOI]` placeholders.

5. **Update manuscript** — Replace all placeholder strings (`[GitHub URL]`, `[Zenodo DOI]`, `[bioRxiv DOI]`, `[ORCID]`) in the manuscript and docs files. See `SUBMISSION_CHECKLIST.md` for the complete list of files to update.

6. **Researcher emails** — Send outreach emails using templates in `researcher_email.md`. Send only after the bioRxiv DOI is live and confirmed. Priority order: Mostafavi Lab → Enformer/Borzoi team → PrediXcan/GTEx researchers.

7. **Journal submission** — Submit `manuscript.pdf` + `docs/cover_letter.md` to *Bioinformatics Advances* (Oxford Academic). See `cover_letter.md` for fallback journal list and reviewer suggestions.

---

## Status note

> **Not peer reviewed.** This is a preprint in preparation for journal submission.
> Please cite the bioRxiv preprint version until peer-reviewed publication.
