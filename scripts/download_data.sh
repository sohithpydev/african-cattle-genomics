#!/usr/bin/env bash
# ==============================================================================
# Helper Script: Download Raw Genomic Reference Dataset from Zenodo
# Dataset: Genomic Reference Resource for African Cattle (GRRFAC)
# DOI: 10.5281/zenodo.8317391
# ==============================================================================

set -euo pipefail

DATA_DIR="data/raw"
mkdir -p "$DATA_DIR"

VCF_FILE="$DATA_DIR/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz"
TBI_FILE="$DATA_DIR/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz.tbi"
EXPECTED_MD5="f7b023667b9f603421b2e589685b2cc6"

VCF_URL="https://zenodo.org/api/records/8317391/files/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz/content"
TBI_URL="https://zenodo.org/api/records/8317391/files/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz.tbi/content"

echo "========================================================================"
echo "African Cattle Genomics: Raw Data Download"
echo "Target directory: $DATA_DIR"
echo "========================================================================"

# Download index file if missing
if [[ ! -f "$TBI_FILE" ]]; then
    echo "Downloading Tabix index ($TBI_FILE)..."
    curl -k -L -C - --progress-bar -o "$TBI_FILE" "$TBI_URL"
    echo "✓ Tabix index downloaded."
else
    echo "✓ Tabix index already exists ($TBI_FILE)."
fi

# Download VCF file if missing
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Downloading raw high-density VCF (~185 MB) from Zenodo (Record 8317391)..."
    curl -k -L -C - --progress-bar -o "$VCF_FILE" "$VCF_URL"
    echo "✓ VCF file downloaded."
else
    echo "✓ VCF file already exists ($VCF_FILE)."
fi

# Verify MD5 checksum
echo "Verifying MD5 checksum..."
ACTUAL_MD5=$(md5 -q "$VCF_FILE" 2>/dev/null || md5sum "$VCF_FILE" | awk '{print $1}')
echo "  Observed MD5: $ACTUAL_MD5"
echo "  Expected MD5: $EXPECTED_MD5"

if [[ "$ACTUAL_MD5" == "$EXPECTED_MD5" ]]; then
    echo "✓ Checksum verification successful."
else
    echo "⚠ Warning: Checksum does not match expected reference ($EXPECTED_MD5)."
    echo "  Please check download integrity."
fi

echo "========================================================================"
echo "Data preparation complete. You can now execute the pipeline:"
echo "  bash workflow/run_pipeline.sh"
echo "========================================================================"
