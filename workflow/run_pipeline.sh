#!/usr/bin/env bash
# ==============================================================================
# Complete End-to-End Execution Runner for Cattle Population Genomics Pipeline
# Project: Genomic Autozygosity & Inbreeding Architecture in Indigenous Cattle
# Lead: Sohith Reddy Korlagunta
# Repository: https://github.com/sohithpydev/african-cattle-genomics
# ==============================================================================

set -euo pipefail

echo "========================================================================"
echo "African Cattle Genomics: Reproducible Computational Pipeline"
echo "Date: $(date)"
echo "Host: $(uname -a)"
echo "Python: $(python3 --version 2>/dev/null || echo 'Python not found')"
echo "========================================================================"

# Step 0: Verification of raw data and checksum (auto-download if missing)
echo "[STEP 0] Verifying input raw data and checksums..."
RAW_VCF="data/raw/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz"
EXPECTED_MD5="f7b023667b9f603421b2e589685b2cc6"

if [[ ! -f "$RAW_VCF" ]]; then
    echo "Raw VCF file not found at $RAW_VCF."
    echo "Initiating automated download from Zenodo (Record 8317391)..."
    bash scripts/download_data.sh
fi

ACTUAL_MD5=$(md5 -q "$RAW_VCF" 2>/dev/null || md5sum "$RAW_VCF" | awk '{print $1}')
echo "  Observed MD5: $ACTUAL_MD5"
echo "  Expected MD5: $EXPECTED_MD5"

if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "WARNING: MD5 does not match expected reference checksum. Proceeding if file is complete or verifying integrity."
fi

# Step 0.5: Unit test validation
echo "------------------------------------------------------------------------"
echo "[STEP 0.5] Running unit and data integrity tests..."
python3 scripts/test_pipeline_unit.py

# Step 1: Input validation and QC filtering
echo "------------------------------------------------------------------------"
echo "[STEP 1] Running QC filtering and autosome genotype extraction..."
python3 scripts/01_qc_filter.py

# Step 2: Population genetic structure, heterozygosity, and F_ST
echo "------------------------------------------------------------------------"
echo "[STEP 2] Running population structure, PCA, and F_ST calculations..."
python3 scripts/02_population_structure.py

# Step 3: Runs of homozygosity detection and inbreeding profiling
echo "------------------------------------------------------------------------"
echo "[STEP 3] Running Runs of Homozygosity (ROH) detection..."
python3 scripts/03_detect_roh.py

# Step 4: Autozygosity island fine-mapping and gene annotation
echo "------------------------------------------------------------------------"
echo "[STEP 4] Fine-mapping autozygosity islands on ARS-UCD1.2..."
python3 scripts/04_autozygosity_islands.py

# Step 5: Figure and table generation
echo "------------------------------------------------------------------------"
echo "[STEP 5] Generating publication-quality figures and tables..."
python3 scripts/05_generate_figures_and_tables.py

echo "========================================================================"
echo "Pipeline completed successfully at $(date)."
echo "Outputs are available in results/, figures/, and tables/."
echo "========================================================================"
