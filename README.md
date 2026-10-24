# African Cattle Genomics: High-Density Autozygosity Landscapes & Inbreeding Architecture

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![DOI: 10.5281/zenodo.8317391](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.8317391-blue.svg)](https://doi.org/10.5281/zenodo.8317391)
[![Assembly: ARS-UCD1.2](https://img.shields.io/badge/Assembly-ARS--UCD1.2-brightgreen.svg)](https://www.ncbi.nlm.nih.gov/assembly/GCF_002263795.1/)
[![Python: >=3.10](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

This repository contains the complete, reproducible computational population genomics pipeline, analysis scripts, pre-computed results, and publication-ready figures for investigating genomic autozygosity, recent vs. historical inbreeding dynamics ($F_{ROH}$), and population differentiation across divergent indigenous African cattle lineages.

---

## 1. Project Overview

Using high-density Illumina BovineHD BeadChip array data (~777k markers mapped to `ARS-UCD1.2`) from the *Genomic Reference Resource for African Cattle* (GRRFAC; Zenodo Record 8317391), this pipeline analyzes **1,063 indigenous African cattle** spanning **40 distinct breeds** across **10 African nations**, retaining **985 quality-verified genomes** and **659,107 autosomal SNPs** following rigorous quality control.

### Primary Biological Findings
1. **Lineage Inbreeding Dynamics:** African Taurine cattle (*Bos taurus africanus*, e.g., N'Dama, Muturu, Baoulé; $N=170$) display significantly higher autozygosity ($F_{ROH} = 0.0735 \pm 0.0717$) than African Zebu (*Bos taurus indicus*, $N=470$, $F_{ROH} = 0.0299 \pm 0.0411$; $P < 10^{-15}$), dominated by short-to-intermediate segments ($1 \le L_{ROH} < 4\text{ Mb}$) reflecting ancient ancestral bottlenecks.
2. **Heterozygosity Buffering in Composites:** Ancient Sanga and Zenga hybrid cattle (e.g., Ankole, Sheko, Kuri; $N=345$) exhibit the lowest autozygosity ($F_{ROH} = 0.0266 \pm 0.0275$) and highest individual heterozygosity ($H_o = 0.3257 \pm 0.0259$), demonstrating the buffering effect of historical taurine-indicine admixture.
3. **Adaptive Selection Sweeps:** Fine-scale autozygosity scanning identified a prominent 3.1 Mb autozygosity peak on **BTA 7: 49.53–52.63 Mb** (reaching 64.04% frequency in Zebu and 45.29% in Taurine), directly harboring ***IL4*** (Interleukin 4), a critical regulator of protective Th2 immune polarization against vector-borne parasites.

---

## 2. Repository Directory Structure

```text
african-cattle-genomics/
├── README.md                      <- Master reproduction guide and project overview
├── LICENSE                        <- MIT code license & CC BY 4.0 data attribution
├── .gitignore                     <- Git ignore rules (excludes >100MB raw files)
├── workflow/
│   └── run_pipeline.sh            <- Master one-command automated pipeline runner
├── scripts/
│   ├── download_data.sh           <- Automated Zenodo data downloader & MD5 check
│   ├── test_pipeline_unit.py      <- Unit test suite (metadata & algorithm checks)
│   ├── 01_qc_filter.py            <- Step 1: QC, variant filtering, genotype extraction
│   ├── 02_population_structure.py <- Step 2: Heterozygosity, LD pruning, PCA, F_ST
│   ├── 03_detect_roh.py           <- Step 3: Consecutive ROH scanning & length binning
│   ├── 04_autozygosity_islands.py <- Step 4: Autozygosity hotspot mapping & annotation
│   └── 05_generate_figures_and_tables.py <- Step 5: Publication figures & CSV tables
├── environment/
│   ├── environment.yml            <- Conda environment definition
│   └── requirements.txt           <- Pinned pip dependencies
├── data/
│   ├── metadata/                  <- Sample cohort metadata & data dictionary
│   │   ├── dataset_inventory.tsv  <- Dataset provenance & checksums
│   │   ├── sample_metadata.tsv    <- Metadata for all 1,082 individual cattle
│   │   └── data_dictionary.md     <- Variable definitions and lineages
│   ├── interim/
│   │   └── variant_metadata.tsv   <- Chromosome positions and missingness (659k SNPs)
│   └── raw/                       <- Raw VCF directory (downloaded via script)
├── results/                       <- Machine-readable tabular pipeline outputs
│   ├── qc_exclusions.tsv          <- Samples excluded during QC
│   ├── breed_diversity_metrics.tsv<- Per-breed Ho, He, and FIS
│   ├── pca_coordinates.tsv        <- Individual PC1–PC4 coordinates
│   ├── pca_variance_explained.tsv <- Variance explained per PC
│   ├── lineage_fst_matrix.tsv     <- Pairwise Hudson F_ST between lineages
│   ├── roh_segments.tsv           <- All detected individual ROH segments
│   ├── individual_inbreeding_froh.tsv <- Per-animal FROH and length-class metrics
│   ├── lineage_froh_summary.tsv   <- Mean and SD of FROH by lineage
│   ├── genome_wide_roh_frequency.tsv  <- Marker-by-marker autozygosity frequencies
│   └── autozygosity_islands.tsv   <- Fine-mapped candidate autozygosity islands
├── figures/                       <- Publication figures (300-DPI PNG and vector PDF)
│   ├── figure_1_study_design_provenance.png / .pdf
│   ├── figure_2_population_structure_pca.png / .pdf
│   ├── figure_3_froh_inbreeding_distributions.png / .pdf
│   ├── figure_4_temporal_roh_length_stratification.png / .pdf
│   └── figure_5_autozygosity_islands_manhattan.png / .pdf
├── tables/                        <- Publication-formatted summary tables
│   ├── table_1_breed_diversity_summary.csv
│   ├── table_2_lineage_froh_stratification.csv
│   └── table_3_top_autozygosity_islands.csv
├── docs/                          <- Preregistered analysis plan & statistical audit
│   ├── analysis_plan.md           <- Methodology protocol & QC thresholds
│   ├── statistical_review.md      <- Robustness analysis & statistical justifications
│   ├── decision_log.md            <- Computational & analytical decision logs
│   └── known_limitations_and_unresolved_items.md
└── references/
    ├── references.bib             <- Verified BibTeX reference library
    └── novelty_matrix.csv         <- Prior art comparison matrix
```

---

## 3. Quick Start: Reproduction in 4 Steps

### Step 1: Clone the Repository
```bash
git clone https://github.com/sohithpydev/african-cattle-genomics.git
cd african-cattle-genomics
```

### Step 2: Set Up the Environment

Using **Pip**:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r environment/requirements.txt
```

Or using **Conda**:
```bash
conda env create -f environment/environment.yml
conda activate cattle-genomics
```

### Step 3: Download the Raw Genomic Data
The raw high-density genotype dataset (~185 MB compressed VCF) is publicly hosted on Zenodo (Record 8317391). Run the automated downloader to fetch and verify the data:
```bash
bash scripts/download_data.sh
```
*Note: This downloads `grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz` and its index into `data/raw/` and verifies the MD5 checksum (`f7b023667b9f603421b2e589685b2cc6`).*

### Step 4: Run the Full Pipeline
To execute the complete end-to-end analysis with a single command:
```bash
bash workflow/run_pipeline.sh
```
*The pipeline will automatically run unit verification, perform QC filtering, compute population structure, scan runs of homozygosity, map autozygosity islands, and generate all publication-ready figures and tables.*

---

## 4. Modular Step-by-Step Execution

Each pipeline step can also be run independently:

| Step | Script | Description | Primary Outputs |
| :--- | :--- | :--- | :--- |
| **0.5** | `scripts/test_pipeline_unit.py` | Validates metadata integrity, lineage rules, and synthetic ROH detection | Unit test summary |
| **1** | `scripts/01_qc_filter.py` | Filters 1,063 samples to autosomes (BTA 1–29), extracts genotypes | `data/interim/clean_genotypes.npz`, `results/qc_exclusions.tsv` |
| **2** | `scripts/02_population_structure.py` | LD pruning ($r^2 < 0.2$), PCA, observed/expected heterozygosity, Hudson $F_{ST}$ | `results/pca_coordinates.tsv`, `results/breed_diversity_metrics.tsv` |
| **3** | `scripts/03_detect_roh.py` | Detects consecutive ROH ($\ge 1\text{ Mb}, \ge 50\text{ SNPs}$) & bins across 5 length classes | `results/roh_segments.tsv`, `results/individual_inbreeding_froh.tsv` |
| **4** | `scripts/04_autozygosity_islands.py` | Marker-by-marker autozygosity frequency, top 1% threshold, gene annotation | `results/genome_wide_roh_frequency.tsv`, `results/autozygosity_islands.tsv` |
| **5** | `scripts/05_generate_figures_and_tables.py` | Produces publication figures (Figures 1–5) and tables (Tables 1–3) | `figures/*.png`, `figures/*.pdf`, `tables/*.csv` |

---

## 5. Pre-Computed Results & Figures

For immediate inspection without re-running compute-intensive genotype processing:
- **`figures/`**: Contains Figures 1 through 5 in both raster (PNG, 300 DPI) and editable vector (PDF) formats.
- **`tables/`**: Contains Tables 1 through 3 in tidy CSV format.
- **`results/`**: Contains raw TSV output files for PCA coordinates, $F_{ST}$ matrices, ROH segment calls, and autozygosity island boundaries.

---

## 6. Data Provenance & References

- **Dataset Source:** *Genomic Reference Resource for African Cattle (GRRFAC)*, Centre for Tropical Livestock Genetics and Health (CTLGH) / ILRI.
- **Data Archive:** Zenodo Record 8317391, DOI: [10.5281/zenodo.8317391](https://doi.org/10.5281/zenodo.8317391).
- **Primary Data Descriptor:** *Nature Scientific Data* 11, 792 (2024), DOI: [10.1038/s41597-024-03589-2](https://doi.org/10.1038/s41597-024-03589-2).
- **Reference Assembly:** *Bos taurus* `ARS-UCD1.2` (GCA_002263795.2).

---

## 7. License

- The code and analysis workflows in this repository are licensed under the **[MIT License](LICENSE)**.
- Data derivatives, tables, and figures are distributed under **[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)**.
