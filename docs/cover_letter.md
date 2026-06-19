# Cover Letter — Journal Submission

---

Dear Editor,

I am submitting for consideration the manuscript titled **"Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction: A Magnitude–Direction Dissociation in Individual-Level Expression Prediction"** as a Research Article.

**What this study does.** This paper reports a computational reanalysis of two large, publicly available datasets: per-gene individual-level Pearson R correlations from the Enformer sequence-to-expression benchmark (Sasse et al. 2023, *Nature Genetics*, n = 6,808 genes, 839 ROSMAP individuals) and GTEx v8 Brain Cortex eGenes (n = 9,082). The analysis characterizes a fundamental dissociation between prediction magnitude and prediction direction in Enformer's individual-level performance, and tests whether this dissociation is attributable to training-tissue mismatch.

**What this study finds.** Enformer's prediction magnitude (|R|) is a strong, monotone predictor of cis-eQTL gene identity (eGene rate rises from 31.6% in the lowest |R| decile to 82.5% in the highest; logistic |R| model AIC = 8,274, outperforming the R+R² model by 104 AIC units). However, the *direction* of Enformer's individual-level predictions for eGenes is near-random (54.0% correct; binomial p = 0.017 against 50%; maximum D10 = 62.5%). A 2×2 factorial analysis confirms that within any |R| stratum, genes with R > 0 and R < 0 have statistically identical eGene rates (50.6% vs. 51.6%; χ² p = 0.571). The direction error is irreversible: R < 0 eGenes retain negative fine-tuned correlations after ROSMAP fine-tuning (mean FinetuneR = −0.035 vs. +0.063; Mann–Whitney p ≈ 10⁻¹⁰⁰). In contrast, PrediXcan correctly assigns positive R to R < 0 eGenes (PrediXcanR = 0.295), confirming that the directional information is accessible from genotype data but inaccessible to Enformer.

**Why this is novel.** This is, to our knowledge, the first characterization of a magnitude–direction dissociation in any sequence-to-expression model evaluated at individual resolution. The finding has immediate practical implications: Enformer can identify which genes are cis-regulated (via |R|), but cannot be used to infer the sign of individual genetic effects on brain expression. The PrediXcan contrast provides a mechanistic anchor: the failure is specific to sequence-based prediction and is most parsimoniously explained by the hypothesis that Enformer's ENCODE cell-line training corpus lacks the brain-specific regulatory programs needed to correctly orient effects in DLPFC neurons (a causal mechanism that would require purpose-trained model variants to confirm).

**Reproducibility.** All analysis code (Python 3, scipy, statsmodels, matplotlib) and processed gene-level tables are available on GitHub and archived on Zenodo. The results can be reproduced entirely from publicly available data without controlled-access approval.

A preprint has been posted on bioRxiv: [insert DOI upon posting]

The author declares no competing interests and no external funding.

Sincerely,

**DongKoo Lee**  
Independent Researcher  
Gwangju-si, Gyeonggi-do, Republic of Korea  
ORCID: 0009-0006-4538-1101  
ceo@nrootm.com

---

## Target journals and rationale

| Priority | Journal | Reason |
|----------|---------|--------|
| 1st | **Bioinformatics Advances** | Computational genomics reanalysis; emphasizes methodological rigor and reproducibility; accepts null/dissociation results |
| 2nd | **PLOS Computational Biology** | Broad computational biology audience; open access; accepts computational reanalyses with novel conceptual contributions |
| 3rd | **PLOS ONE** | Methodological soundness over novelty; CC BY aligned; no novelty threshold |
| 4th | **GigaScience** | Emphasizes data/code packages; strong fit for reproducibility-first submissions |
| 5th | **BMC Genomics** | Open access; accepts computational reanalyses without experimental validation |

**Avoid initially:** Nature Genetics, Nature Methods, Genome Biology — require either experimental validation, multi-cohort replication, or mechanistic follow-up.

---

## Reviewer suggestions

Suggested reviewers (authors of relevant papers):

1. **Alexander Sasse** (Mostafavi Lab, UBC) — corresponding author of the original benchmark
2. **Sarah Mostafavi** (UBC/University of Washington) — PI of EnformerAssessment
3. **Julien Bryois** (Roche, formerly EPFL) — cortical eQTL and single-cell DLPFC specialist
4. **Eric R. Gamazon** (Vanderbilt) — PrediXcan developer; directly relevant to the genotype model comparison
5. **Francois Aguet** (Broad Institute, GTEx) — GTEx analysis lead

Reviewers to avoid: none required to disclose.
