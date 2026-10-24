# Project Decision Log & Audit Trail

**Project Title:** Genomic Dissection of Autozygosity Architecture, Inbreeding Dynamics, and Population Differentiation in Indigenous Tropical Cattle
**Target Programme:** European Master in Animal Biodiversity and Genomics (EMABG)
**Lead Researcher / Applicant:** Sohith Reddy Korlagunta
**Log Initialized:** 2026-10-24T14:14:00+05:30

---

## Log Entry 001: Audit Discovery & Baseline Intake
- **Date & Time:** 2026-10-24T14:14:00+05:30
- **Action:** Located and thoroughly ingested existing audit file: `/Users/sohith/Desktop/emabg-archive/EMABG_Applicant_Portfolio_Audit_Sohith.md`.
- **Source Files Verified:**
  - `EMABG_Applicant_Portfolio_Audit_Sohith.md` (387 lines, 37.8 KB)
  - `emabg_admissions_dossier.html` (67.9 KB)
  - `SOHITH_MORPHOPHEN_CV.pdf` (101.6 KB)
  - `thesis.txt` (92.6 KB, 35+ EMABG master's theses across NMBU, WUR, BOKU, SLU, APT, UGOE)
  - `project-skill.md` (398 lines, 37.4 KB)
- **Key Audit Findings:**
  - Eligibility confirmed (B.E. Bioinformatics, 9.3/10 CGPA, 3-time University Topper, English proficiency passed).
  - Decisive weakness: Thematic bias toward human neurodegenerative pharmacology/CADD, with zero documented farm animal genomics, population genetics, or livestock breeding credits/projects.
  - Recommended intervention: Execute a high-rigor, reproducible dry-lab animal population genomics project demonstrating PLINK/R/Python pipelines, runs of homozygosity (ROH), population structure, and diversity metrics directly matching EMABG coursework (`BIO321`, `HFA350`, `MAT20306`).

## Log Entry 002: Computational Environment & Infrastructure Audit
- **Date & Time:** 2026-10-24T14:14:30+05:30
- **Hardware Platform:** Apple MacBook Air (Apple M1, 8-core CPU, ARM64 architecture, Darwin 23.6.0).
- **Storage Status:** 36 GB available on `/System/Volumes/Data`.
- **Software Stack:**
  - Python 3.13.5 (`/opt/miniconda3/bin/python3`) with NumPy (2.3.5), SciPy, Pandas, Scikit-learn, Matplotlib.
  - R version 4.4.1 (2024-06-14, `/usr/local/bin/R`) with Bioconductor foundation.
  - Git 2.39+ (`/opt/homebrew/bin/git`).
  - Homebrew (`/opt/homebrew/bin/brew`) available for formula installations (e.g. `bcftools`, `plink2`).
- **Resource Boundary Decision:**
  - Raw uncompressed whole-genome sequencing (WGS) FASTQs for hundreds of samples require >1–5 TB of disk space and hundreds of CPU hours, which exceeds local laptop storage (36 GB).
  - High-density (HD) SNP array / joint VCF genotype datasets (~100–300 MB compressed, covering hundreds of thousands of markers across 1,000+ individuals) are perfectly matched to local memory and storage, enabling complete, un-truncated, end-to-end execution, full statistical power, and strict reproducibility.

## Log Entry 003: Pipeline Implementation & Algorithm Standardization
- **Date & Time:** 2026-10-24T14:25:00+05:30
- **Action:** Engineered modular, deterministic Python and bash pipeline components:
  - `scripts/01_qc_filter.py`: Quality control, autosome filtering, sample matching against metadata, matrix extraction.
  - `scripts/02_population_structure.py`: Heterozygosity, inbreeding deficiency ($F_{IS}$), LD thinning ($r^2 < 0.2$), SVD PCA, and pairwise Hudson's $F_{ST}$.
  - `scripts/03_detect_roh.py`: Consecutive sliding-window ROH scan with 1.0 Mb gap cutoff (Nandolo et al. 2018), 5 temporal length bins ($F_{ROH, k}$).
  - `scripts/04_autozygosity_islands.py`: $100\text{ kb}$ window-based autozygosity incidence scanning ($P_{ROH}$), empirical top 1% thresholding, functional gene annotation on `ARS-UCD1.2`.
  - `scripts/05_generate_figures_and_tables.py`: Vector PDF/PNG rendering (300 DPI) and CSV publication tables.

## Log Entry 004: Unit Testing & Integrity Smoke Test Verification
- **Date & Time:** 2026-10-24T14:32:00+05:30
- **Action:** Executed `scripts/test_pipeline_unit.py`.
- **Result:** All 4 unit test suites passed with zero errors (`Ran 4 tests in 0.016s - OK`):
  - Validated metadata integrity (1,082 records, 1,063 included, 19 excluded, 40 breeds).
  - Verified lineage classification logic across all geographical population categories.
  - Verified deterministic length binning across all 5 discrete temporal classes.
  - Confirmed synthetic ROH identification and negative control rejection.

## Log Entry 005: Full Execution, Numerical Synchronization, and Submission Compilation
- **Date & Time:** 2026-10-24T15:02:00+05:30
- **Action:**
  - Completed end-to-end execution of the full analytical workflow across 985 quality-filtered cattle individuals and 659,107 autosomal SNPs.
  - Identified 34,955 distinct runs of homozygosity segments and 72 fine-mapped autozygosity islands.
  - Confirmed primary hypotheses:
    1. African Taurine exhibit elevated autozygosity ($F_{ROH} = 0.0735 \pm 0.0717$) driven by ancient demographic bottlenecks ($1\text{--}4\text{ Mb}$ segments) rather than recent inbreeding.
    2. Sanga and Zenga composites demonstrate admixture buffering with the lowest autozygosity ($F_{ROH} = 0.0266 \pm 0.0275$) and highest observed heterozygosity ($H_o = 0.3257 \pm 0.0259$).
    3. Shared autozygosity megabase hotspot on BTA 7 (49.53–52.63 Mb, peak frequency 64.04% in Zebu) overlaps the master Th2 cytokine *IL4*.
  - Synchronized all text in `manuscript/manuscript.md`, `manuscript/cover_letter.md`, `manuscript/highlights.md`, and `README.md` to exact empirical table values.
  - Created complete `manuscript/supplementary_materials.md`.
  - Compiled professional Word drafts `manuscript/manuscript.docx` and `manuscript/supplementary_materials.docx` via Pandoc.
  - Updated `scripts/03_detect_roh.py` with modular functions and verified unit tests (`scripts/test_pipeline_unit.py`) pass in 0.015s.
- **Milestone Status:** Phase 10 Gate completed; project is **SUBMISSION-READY FOR AUTHOR REVIEW**.

