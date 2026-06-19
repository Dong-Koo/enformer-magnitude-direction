# Researcher Outreach Emails

Send AFTER bioRxiv DOI is obtained. Copy-paste and fill in the name/DOI.

---

## Priority 1: Mostafavi Lab (original data authors)

**To:** mostafavi.lab@ubc.ca  *(or find direct contact from lab website)*  
**CC:** Alexander Sasse (check EnformerAssessment repo for current affiliation)
**Subject:** Reanalysis of EnformerAssessment reveals magnitude–direction dissociation in individual-level predictions

Dear Dr. Mostafavi (and Dr. Sasse),

I am an independent computational genomics researcher. I recently reanalyzed the publicly
available EnformerAssessment per-gene correlation data together with GTEx v8 Brain Cortex
eGenes, and identified a finding I believe you will find interesting.

The main result is a magnitude–direction dissociation: Enformer's prediction *magnitude*
(|R|) monotonically predicts cis-eQTL gene identity (eGene rate: 31.6% at lowest |R| → 82.5%
at highest |R|; logistic model AIC advantage = 104 over the R+R² U-shape), but the *direction*
of Enformer's individual predictions for eGenes is near-random (54.0% correct; maximum D10 =
62.5%). A 2×2 factorial analysis confirms that within any |R| stratum, R>0 and R<0 genes have
identical eGene rates (p = 0.571). Critically, PrediXcan correctly predicts positive expression
for R<0 eGenes (PrediXcanR = 0.295), suggesting the directional information is accessible from
genotype but inaccessible to sequence-based prediction—consistent with training-tissue mismatch
(ENCODE cell lines vs. DLPFC neurons).

Preprint: [bioRxiv DOI]  
Code + processed data: [GitHub URL] · [Zenodo DOI]

I would be very grateful for any feedback, particularly on whether the tissue mismatch
interpretation (K562/GM12878 vs. DLPFC) is consistent with your observations during the
benchmark. I would also welcome any correction if I have misunderstood any aspect of the
EnformerAssessment data.

Best regards,  
DongKoo Lee  
Independent Researcher, Republic of Korea  
ORCID: [ORCID] · ceo@nrootm.com

---

## Priority 2: Enformer/Borzoi team (Google DeepMind / Kundaje Lab)

**Subject:** Computational reanalysis: magnitude–direction dissociation in Enformer individual-level prediction

Dear Dr. [Kundaje / Avsec],

I am reaching out regarding a computational reanalysis of Enformer's individual-level prediction
performance (Sasse et al. 2023 benchmark) crossed with GTEx v8 Brain Cortex eGenes.

The key finding: Enformer's prediction *magnitude* (|R|) is a strong, monotone detector of
cis-eQTL gene identity, but its *direction* is near-random for eGenes (54% correct). Direction
errors are stable to fine-tuning on 839 ROSMAP individuals (p ≈ 10⁻¹⁰⁰ for the sign group
difference) and the directional signal IS recoverable by PrediXcan (a genotype-based linear
model, R = 0.295 for wrongly-directed eGenes). The most parsimonious explanation is
training-tissue mismatch: brain-specific cis-regulatory effects are absent or inverted in
Enformer's ENCODE cell-line training data.

This suggests a concrete target for Borzoi or next-generation models: incorporating
brain-tissue-specific chromatin data should recover directional accuracy while retaining
magnitude-based sensitivity.

Preprint: [bioRxiv DOI]  
Code + data: [GitHub / Zenodo]

I would welcome any feedback, especially on whether the observation is consistent with the
model's known limitations, and whether any internal analyses have addressed direction accuracy
for tissue-specific eQTL genes.

Best regards,  
DongKoo Lee · ceo@nrootm.com · ORCID: [ORCID]

---

## Priority 3: PrediXcan / GTEX researchers

**Subject:** Magnitude–direction dissociation in Enformer vs. PrediXcan — preprint

Dear Dr. [Gamazon / GTEx lead],

I am an independent researcher and have just posted a preprint examining the relationship
between Enformer individual-level predictions and GTEx v8 Brain Cortex eGenes.

A key finding directly relevant to your work: while Enformer's prediction magnitude (|R|)
identifies cis-eQTL genes with high sensitivity (eGene rate 82.5% at highest |R| decile),
its *direction* is near-random (54% correct for eGenes). In contrast, PrediXcan correctly
predicts the direction of expression for genes where Enformer fails (PrediXcanR = 0.295 for
R<0 eGenes). This suggests that the information needed for correct directional prediction of
cis-effects is captured by genotype-based linear models but is inaccessible to sequence-based
models trained on non-brain tissue.

Preprint: [bioRxiv DOI]

I would appreciate your perspective on whether this comparison is fair given the different
training frameworks, and whether the ~54% direction accuracy for eGenes is consistent with
any theoretical expectations for sequence models trained on non-target tissues.

Best regards,  
DongKoo Lee · ceo@nrootm.com

---

## After receiving a response

- Acknowledge within 48 hours
- If feedback leads to revision: update preprint version and re-notify
- If reviewer suggests addition/correction: add to Limitations or Discussion
- Do NOT commit to authorship changes based on email feedback alone
