# 출판 완전 가이드 (복붙 입력 포함)

**현재 상태:** 로컬 git 커밋 완료 (`main` 브랜치, 29개 파일)  
**다음 단계:** 아래 순서대로 진행. 각 단계에서 입력할 내용 전부 포함.

---

## STEP 0. GitHub CLI 로그인 (터미널에서 실행)

터미널을 열고 아래 명령어 실행:

```bash
gh auth login
```

프롬프트 응답:
```
? What account do you want to log into?  → GitHub.com  (엔터)
? What is your preferred protocol?        → HTTPS       (엔터)
? Authenticate GitHub CLI?                → Login with a web browser  (엔터)
→ 브라우저가 열림 → GitHub 계정으로 로그인 → Authorize 클릭
```

로그인 완료 후 아래 명령어로 repo 생성 + 업로드:

```bash
cd /Users/dongkoo/work/KnowledgeWork/eqtl_enformer_ushape

gh repo create enformer-magnitude-direction \
  --public \
  --description "Enformer recognizes cis-regulated genes by prediction magnitude but fails to predict effect direction (Lee 2026)" \
  --source=. \
  --remote=origin \
  --push
```

완료 시 URL이 출력됩니다: `https://github.com/[your-username]/enformer-magnitude-direction`

---

## STEP 1. ORCID 계정 생성

**URL:** https://orcid.org/register

### 입력 내용

| 필드 | 입력값 |
|------|--------|
| First name | DongKoo |
| Last name | Lee |
| Email | ceo@nrootm.com |
| Confirm email | ceo@nrootm.com |
| Password | (본인 설정) |
| Email visibility | Everyone (공개 권장) |

→ 이메일 인증 링크 클릭  
→ 로그인 후 ORCID ID 확인: `https://orcid.org/0000-0000-0000-XXXX` 형태

### ORCID에 논문 추가

1. 로그인 → "Works" 섹션 → "Add works" → "Add manually"
2. 입력:

| 필드 | 입력값 |
|------|--------|
| Work type | Preprint |
| Title | Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction: A Magnitude–Direction Dissociation in Individual-Level Expression Prediction |
| Journal/Repository | bioRxiv |
| Publication year | 2026 |
| DOI | (bioRxiv DOI — STEP 3 완료 후 추가) |

→ "Save" 클릭

---

## STEP 2. Zenodo — GitHub 연동 + DOI 발급

**URL:** https://zenodo.org

### 계정 생성

→ "Sign up" → **"Sign up with GitHub"** 클릭 (GitHub 계정 연동이 가장 간편)  
→ GitHub 계정으로 로그인 → "Authorize zenodo" 클릭

### GitHub 연동

1. 로그인 후 우측 상단 이름 → "GitHub" 클릭
2. `enformer-magnitude-direction` 저장소 찾기 → 토글 ON

### DOI 발급

1. GitHub에서 저장소 → "Releases" → "Create a new release"
2. 입력:

| 필드 | 입력값 |
|------|--------|
| Tag version | v1.0.0 |
| Release title | Initial release — Lee 2026 preprint |
| Description | First public release accompanying the bioRxiv preprint. Contains manuscript, analysis code, processed data, and figures. |

3. "Publish release" 클릭  
→ Zenodo에서 자동으로 DOI 발급: `10.5281/zenodo.20754856`

### Zenodo 메타데이터 확인

Zenodo 대시보드 → 해당 업로드 → "Edit" → 자동 입력 확인:
- Title, creators, keywords는 `.zenodo.json`에서 자동 로드됨
- "Submit for review" 클릭

---

## STEP 3. bioRxiv 논문 업로드

**URL:** https://www.biorxiv.org/submit

### 계정 생성

→ "Register" → 입력:

| 필드 | 입력값 |
|------|--------|
| First name | DongKoo |
| Last name | Lee |
| Email | ceo@nrootm.com |
| Institution | Independent Researcher |
| Country | Republic of Korea |

### 논문 업로드

1. 로그인 → "Submit a manuscript"
2. **Category:** Genomics
3. **Title:**
```
Enformer Recognizes cis-Regulated Genes by Prediction Magnitude but Fails to Predict Effect Direction: A Magnitude–Direction Dissociation in Individual-Level Expression Prediction
```

4. **Authors:**
```
DongKoo Lee
Independent Researcher, Gwangju-si, Gyeonggi-do, Republic of Korea
ceo@nrootm.com
ORCID: [STEP 1에서 발급된 ORCID]
```

5. **Abstract** (복붙):
```
Enformer, a sequence-to-expression deep learning model, has been benchmarked 
for individual-level prediction in the ROSMAP cohort (839 individuals, DLPFC 
bulk RNA-seq). We cross-referenced per-gene Pearson R correlations from the 
Enformer benchmark (Sasse et al. 2023) with GTEx v8 Brain Cortex eGenes 
(n=9,082) across 6,808 autosomal genes and find a magnitude-direction 
dissociation: prediction magnitude (|R|) monotonically predicts cis-eQTL gene 
identity (eGene rate 31.6% to 82.5%; |R| model AIC=8274 vs R+R² AIC=8378, 
ΔAIC=104), yet direction is near-random for eGenes (54.0% correct; maximum 
D10=62.5%). A 2x2 factorial analysis confirms that within any |R| stratum, 
genes with R>0 and R<0 have identical eGene rates (50.6% vs 51.6%; p=0.571). 
Direction errors are irreversible by ROSMAP fine-tuning (p≈10⁻¹⁰⁰) but 
recoverable by PrediXcan (R=0.295 for R<0 eGenes). We attribute the direction 
failure to training-tissue mismatch: Enformer's ENCODE cell-line training corpus 
lacks brain-specific regulatory programs needed to correctly orient effects in 
DLPFC neurons.
```

6. **Manuscript file:** `manuscript/manuscript.html` → PDF로 변환 필요

#### HTML → PDF 변환 (터미널):
```bash
# Chrome headless로 PDF 변환
/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome \
  --headless --print-to-pdf=/Users/dongkoo/work/KnowledgeWork/eqtl_enformer_ushape/manuscript/manuscript.pdf \
  /Users/dongkoo/work/KnowledgeWork/eqtl_enformer_ushape/manuscript/manuscript.html
```

7. **Figures:** fig1, fig2, fig3 PNG 파일 각각 업로드

8. **Subject area:** Genomics  
   **Related links:**
   - GitHub: `https://github.com/Dong-Koo/enformer-magnitude-direction`
   - Zenodo: `https://zenodo.org/record/XXXXXXX`

9. "Submit" → 약 24-48시간 후 DOI 발급: `10.1101/XXXXXX.XXXXX`

---

## STEP 4. 파일 내 placeholder 업데이트

bioRxiv DOI, Zenodo DOI, ORCID, GitHub URL 발급 후 아래 명령어로 일괄 교체:

```bash
cd /Users/dongkoo/work/KnowledgeWork/eqtl_enformer_ushape

GITHUB_URL="https://github.com/Dong-Koo/enformer-magnitude-direction"
ZENODO_DOI="10.5281/zenodo.20754856"
BIORXIV_DOI="10.1101/XXXXXX.XXXXX"
ORCID="0000-0000-0000-XXXX"

# README
sed -i '' "s|\[GitHub URL\]|$GITHUB_URL|g" README.md
sed -i '' "s|\[Zenodo DOI\]|$ZENODO_DOI|g" README.md
sed -i '' "s|10\.1101/XXXXXXXX|$BIORXIV_DOI|g" README.md

# CITATION.cff
sed -i '' "s|10\.1101/XXXXXXXX|$BIORXIV_DOI|g" CITATION.cff
sed -i '' "s|\[to be assigned\]|$ORCID|g" CITATION.cff

# manuscript
sed -i '' "s|\[GitHub URL\]|$GITHUB_URL|g" manuscript/manuscript.html
sed -i '' "s|\[Zenodo DOI to be assigned\]|$ZENODO_DOI|g" manuscript/manuscript.html
sed -i '' "s|\[bioRxiv DOI to be assigned\]|$BIORXIV_DOI|g" manuscript/manuscript.html
sed -i '' "s|\[to be assigned\]|$ORCID|g" manuscript/manuscript.html

# commit
git add -A && git commit -m "Add DOIs and ORCID after preprint registration"
git push origin main
```

---

## STEP 5. 연구자 이메일 발송

bioRxiv DOI 발급 확인 후 `docs/researcher_email.md` 내용 발송:

1. **Priority 1:** Mostafavi Lab — mostafavi.lab@ubc.ca
2. **Priority 2:** Enformer/Kundaje Lab
3. **Priority 3:** Gamazon Lab (PrediXcan 개발자)

---

## STEP 6. 저널 투고

**1차 타겟: Bioinformatics Advances**  
URL: https://academic.oup.com/bioinformaticsadvances/pages/submission_guidelines

Cover letter 내용: `docs/cover_letter.md` 참조

---

## 전체 타임라인 요약

| 일정 | 작업 |
|------|------|
| Day 1 | STEP 0 (GitHub) + STEP 1 (ORCID) |
| Day 1 | STEP 2 (Zenodo release — GitHub 연동 즉시 가능) |
| Day 1-2 | STEP 3 (bioRxiv 제출) |
| Day 2-4 | bioRxiv DOI 발급 대기 |
| Day 4 | STEP 4 (placeholder 업데이트) + STEP 5 (이메일) |
| Day 5+ | STEP 6 (저널 투고) |
