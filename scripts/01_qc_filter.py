#!/usr/bin/env python3
"""
Step 1: Input Validation, Quality Control, and Genotype Extraction (Optimized)
=============================================================================
Project: Genomic Autozygosity & Inbreeding Architecture in Indigenous Tropical Cattle
Author: Sohith Reddy Korlagunta / Autonomous Pipeline
Target: EMABG Portfolio Project

Operations:
1. Parse sample metadata from data/metadata/sample_metadata.tsv.
2. Inspect input VCF (data/raw/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz).
3. Validate sample header IDs against metadata; confirm 1,063 passing samples.
4. Filter variants to bovine autosomes BTA 1 to BTA 29.
5. Compute sample-level and variant-level missingness.
6. Record exclusions in results/qc_exclusions.tsv.
7. Save filtered autosome genotype matrix to data/interim/clean_genotypes.npz.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import pysam

def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def main():
    log("Starting Step 1: Input Validation and QC Filtering...")

    # Define paths
    vcf_path = "data/raw/grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz"
    meta_path = "data/metadata/sample_metadata.tsv"
    out_qc_exclusions = "results/qc_exclusions.tsv"
    out_interim = "data/interim/clean_genotypes.npz"
    out_variants_meta = "data/interim/variant_metadata.tsv"

    os.makedirs("results", exist_ok=True)
    os.makedirs("data/interim", exist_ok=True)

    # 1. Load metadata
    log(f"Loading metadata from {meta_path}...")
    meta_df = pd.read_csv(meta_path, sep="\t")
    log(f"Metadata loaded: {len(meta_df)} total records.")

    included_meta = meta_df[meta_df["Study_Cohort"] == "INCLUDED"].copy()
    valid_sample_set = set(included_meta["VCF_ID"].astype(str).str.strip())
    log(f"Included study samples in metadata: {len(valid_sample_set)}.")

    # 2. Inspect VCF header
    if not os.path.exists(vcf_path):
        log(f"ERROR: VCF file not found at {vcf_path}!")
        sys.exit(1)

    vcf_size_mb = os.path.getsize(vcf_path) / (1024 * 1024)
    log(f"Inspecting VCF file ({vcf_size_mb:.2f} MB)...")

    # Read VCF sample names using pysam
    vcf_in = pysam.VariantFile(vcf_path)
    vcf_samples = list(vcf_in.header.samples)
    n_vcf_samples = len(vcf_samples)
    log(f"VCF sample header count: {n_vcf_samples}")

    # Check sample overlap
    matched_samples = [s for s in vcf_samples if s in valid_sample_set]
    log(f"Matching samples between VCF and metadata: {len(matched_samples)}")

    # Write initial QC exclusions record
    qc_records = []
    excluded_meta = meta_df[meta_df["Study_Cohort"] == "EXCLUDED"]
    for _, r in excluded_meta.iterrows():
        qc_records.append({
            "Sample_or_Variant_ID": str(r["VCF_ID"]),
            "Type": "Sample",
            "Reason": str(r["Exclusion_Reason"]),
            "Stage": "Provider_Pre_QC",
            "Action": "Excluded"
        })

    # 3. Parse autosomal variants and genotypes
    valid_chroms = {str(i): int(i) for i in range(1, 30)}
    valid_chroms.update({f"chr{i}": int(i) for i in range(1, 30)})

    log("Scanning autosomes (BTA 1 to 29) using high-speed integer dosage extraction...")
    var_chroms = []
    var_positions = []
    var_ids = []

    genotypes_blocks = []
    current_block = []
    block_size = 50000

    total_vars_scanned = 0
    autosomal_vars_count = 0
    non_autosomal_count = 0

    t0 = time.time()
    for record in vcf_in:
        total_vars_scanned += 1
        chrom = record.chrom
        chrom_clean = chrom.replace("chr", "").replace("BTA", "")

        if chrom_clean not in valid_chroms and chrom not in valid_chroms:
            non_autosomal_count += 1
            continue

        std_chrom = valid_chroms.get(chrom_clean, valid_chroms.get(chrom, None))
        if std_chrom is None:
            non_autosomal_count += 1
            continue

        autosomal_vars_count += 1
        var_chroms.append(std_chrom)
        var_positions.append(record.pos)
        var_ids.append(record.id if record.id else f"{std_chrom}_{record.pos}")

        # High-speed numpy fromiter extraction across sample proxies
        row = np.fromiter(
            ( (a[0] + a[1]) if (a is not None and len(a) == 2 and a[0] is not None and a[1] is not None and a[0] <= 1 and a[1] <= 1) else -1
              for a in (s.allele_indices for s in record.samples.values()) ),
            dtype=np.int8,
            count=n_vcf_samples
        )
        current_block.append(row)

        if len(current_block) == block_size:
            genotypes_blocks.append(np.array(current_block, dtype=np.int8))
            current_block = []
            elapsed = time.time() - t0
            rate = autosomal_vars_count / elapsed
            log(f"  Processed {autosomal_vars_count} autosomal variants ({elapsed:.1f}s, {rate:.0f} vars/s)...")

    if current_block:
        genotypes_blocks.append(np.array(current_block, dtype=np.int8))
        current_block = []

    vcf_in.close()

    total_elapsed = time.time() - t0
    log(f"Variant scanning complete in {total_elapsed:.1f}s: {total_vars_scanned} total scanned, {autosomal_vars_count} autosomal retained, {non_autosomal_count} non-autosomal excluded.")

    if non_autosomal_count > 0:
        qc_records.append({
            "Sample_or_Variant_ID": "Non_Autosomal_Markers",
            "Type": "Variant_Batch",
            "Reason": f"Excluded {non_autosomal_count} non-autosomal/unplaced markers",
            "Stage": "Autosome_Filtering",
            "Action": "Excluded"
        })

    # Concatenate blocks
    log("Concatenating genotype blocks into contiguous array...")
    gt_matrix = np.vstack(genotypes_blocks)
    del genotypes_blocks

    n_variants, n_samples = gt_matrix.shape
    log(f"Genotype matrix assembled: {n_variants} variants by {n_samples} samples.")

    # 4. Variant missingness filter (threshold: <= 10% missing)
    log("Applying variant call rate filter (missingness <= 10%)...")
    missing_mask = (gt_matrix == -1)
    var_missing_rates = np.mean(missing_mask, axis=1)
    pass_var_mask = (var_missing_rates <= 0.10)
    n_vars_pass = int(np.sum(pass_var_mask))
    n_vars_fail = n_variants - n_vars_pass
    log(f"Variant call rate filter: {n_vars_pass} pass, {n_vars_fail} fail.")

    if n_vars_fail > 0:
        qc_records.append({
            "Sample_or_Variant_ID": "High_Missingness_Variants",
            "Type": "Variant_Batch",
            "Reason": f"Missing rate > 0.10 (N={n_vars_fail})",
            "Stage": "Variant_QC",
            "Action": "Excluded"
        })

    # 5. Sample missingness filter (threshold: <= 5% missing among passing variants)
    log("Applying sample call rate filter (missingness <= 5%)...")
    gt_filtered = gt_matrix[pass_var_mask, :]
    sample_missing_rates = np.mean(gt_filtered == -1, axis=0)
    pass_sample_mask = (sample_missing_rates <= 0.05)
    n_samples_pass = int(np.sum(pass_sample_mask))
    n_samples_fail = n_samples - n_samples_pass
    log(f"Sample call rate filter: {n_samples_pass} pass, {n_samples_fail} fail.")

    final_gt = gt_filtered[:, pass_sample_mask]
    final_samples = np.array(vcf_samples)[pass_sample_mask]
    final_chroms = np.array(var_chroms, dtype=np.int8)[pass_var_mask]
    final_positions = np.array(var_positions, dtype=np.int32)[pass_var_mask]
    final_ids = np.array(var_ids)[pass_var_mask]

    log(f"Final clean dataset: {final_gt.shape[0]} variants x {final_gt.shape[1]} samples.")

    # Save QC exclusions table
    qc_df = pd.DataFrame(qc_records)
    qc_df.to_csv(out_qc_exclusions, sep="\t", index=False)
    log(f"Saved QC exclusions to {out_qc_exclusions}.")

    # Save variant metadata
    var_meta_df = pd.DataFrame({
        "Chrom": final_chroms,
        "Position": final_positions,
        "Variant_ID": final_ids,
        "Missing_Rate": np.round(var_missing_rates[pass_var_mask], 4)
    })
    var_meta_df.to_csv(out_variants_meta, sep="\t", index=False)
    log(f"Saved variant metadata to {out_variants_meta}.")

    # Save compressed clean genotypes
    log(f"Saving compressed interim genotypes to {out_interim}...")
    np.savez_compressed(
        out_interim,
        genotypes=final_gt,
        samples=final_samples,
        chromosomes=final_chroms,
        positions=final_positions,
        variant_ids=final_ids
    )
    log("Step 1 completed successfully.")

if __name__ == "__main__":
    main()
