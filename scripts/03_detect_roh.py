#!/usr/bin/env python3
"""
Step 3: Runs of Homozygosity (ROH) Detection & Inbreeding Architecture (Vectorized)
===================================================================================
Project: Genomic Autozygosity & Inbreeding Architecture in Indigenous Tropical Cattle
Author: Sohith Reddy Korlagunta / Autonomous Pipeline
Target: EMABG Portfolio Project

Operations:
1. Load clean genotype matrix and variant positions from data/interim/clean_genotypes.npz.
2. Detect runs of homozygosity (ROH) per individual across chromosomes BTA 1 to 29:
   - Minimum physical length: 1,000,000 bp (1.0 Mb).
   - Minimum marker count: 50 consecutive homozygous SNPs.
   - Maximum inter-marker gap: 1,000,000 bp (Nandolo et al. 2018).
   - Density requirement: at least 1 SNP per 50 kb.
3. Classify ROH segments into 5 discrete temporal length classes:
   - 1-2 Mb (ancient background autozygosity, ~25-50 generations)
   - 2-4 Mb (historical inbreeding, ~12-25 generations)
   - 4-8 Mb (intermediate inbreeding, ~6-12 generations)
   - 8-16 Mb (recent inbreeding, ~3-6 generations)
   - >=16 Mb (very recent consanguinity, <=3 generations)
4. Compute total F_ROH and stratified F_ROH_k per individual and breed.
5. Save results to results/roh_segments.tsv and results/individual_inbreeding_froh.tsv.
"""

import os
import sys
import time
import numpy as np
import pandas as pd

# ARS-UCD1.2 bovine autosomal genome length (chromosomes 1 to 29) in megabases
L_AUTOSOME_MB = 2489.38

def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def classify_length(length_mb):
    if length_mb < 2.0:
        return "1-2 Mb"
    elif length_mb < 4.0:
        return "2-4 Mb"
    elif length_mb < 8.0:
        return "4-8 Mb"
    elif length_mb < 16.0:
        return "8-16 Mb"
    else:
        return ">=16 Mb"

def detect_roh_chromosome(c_calls, p, chrom=1, s_id="TEST", min_length_bp=1000000, min_snps=50, max_gap_bp=1000000):
    """
    Detect runs of homozygosity on a single chromosome for one individual.
    """
    ind_roh = []
    is_het = (c_calls == 1)
    padded = np.concatenate(([True], is_het, [True]))
    diffs = np.diff(padded.astype(np.int8))
    starts = np.where(diffs == -1)[0]
    ends = np.where(diffs == 1)[0]

    lengths_bp = p[ends - 1] - p[starts]
    snp_counts = ends - starts

    cand_mask = (lengths_bp >= min_length_bp) & (snp_counts >= min_snps)
    cand_starts = starts[cand_mask]
    cand_ends = ends[cand_mask]

    for s_start, s_end in zip(cand_starts, cand_ends):
        seg_pos = p[s_start:s_end]
        if len(seg_pos) > 1:
            gaps = np.diff(seg_pos)
            large_gaps = np.where(gaps > max_gap_bp)[0]
            if len(large_gaps) > 0:
                split_indices = [0] + list(large_gaps + 1) + [len(seg_pos)]
                for k in range(len(split_indices) - 1):
                    sub_start = s_start + split_indices[k]
                    sub_end = s_start + split_indices[k+1]
                    sub_len_bp = p[sub_end - 1] - p[sub_start]
                    sub_snps = sub_end - sub_start
                    if sub_len_bp >= min_length_bp and sub_snps >= min_snps and sub_len_bp <= sub_snps * 50000:
                        len_mb = sub_len_bp / 1000000.0
                        ind_roh.append({
                            "Sample_ID": s_id,
                            "Chrom": int(chrom),
                            "Start_Pos": int(p[sub_start]),
                            "End_Pos": int(p[sub_end - 1]),
                            "Length_bp": int(sub_len_bp),
                            "Length_Mb": round(len_mb, 4),
                            "SNP_Count": int(sub_snps),
                            "Length_Class": classify_length(len_mb)
                        })
                continue

        len_bp = p[s_end - 1] - p[s_start]
        snp_cnt = s_end - s_start
        if len_bp <= snp_cnt * 50000:
            len_mb = len_bp / 1000000.0
            ind_roh.append({
                "Sample_ID": s_id,
                "Chrom": int(chrom),
                "Start_Pos": int(p[s_start]),
                "End_Pos": int(p[s_end - 1]),
                "Length_bp": int(len_bp),
                "Length_Mb": round(len_mb, 4),
                "SNP_Count": int(snp_cnt),
                "Length_Class": classify_length(len_mb)
            })
    return ind_roh

def detect_roh_individual(calls, chrom_arr, pos_arr, s_id="TEST", min_length_bp=1000000, min_snps=50, max_gap_bp=1000000):
    """
    Detect runs of homozygosity genome-wide across all chromosomes for one individual.
    """
    unique_chroms = np.unique(chrom_arr)
    ind_roh = []
    for c in unique_chroms:
        mask = (chrom_arr == c)
        ind_roh.extend(detect_roh_chromosome(calls[mask], pos_arr[mask], chrom=c, s_id=s_id,
                                             min_length_bp=min_length_bp, min_snps=min_snps, max_gap_bp=max_gap_bp))
    return ind_roh

def main():
    log("Starting Step 3: Runs of Homozygosity Detection & Inbreeding Profiling (Vectorized)...")

    in_genotypes = "data/interim/clean_genotypes.npz"
    pca_meta_path = "results/pca_coordinates.tsv"
    out_roh_segments = "results/roh_segments.tsv"
    out_individual_froh = "results/individual_inbreeding_froh.tsv"
    out_lineage_summary = "results/lineage_froh_summary.tsv"

    if not os.path.exists(in_genotypes):
        log(f"ERROR: {in_genotypes} not found! Run scripts/01_qc_filter.py first.")
        sys.exit(1)

    # 1. Load genotypes
    log(f"Loading genotypes from {in_genotypes}...")
    npz_data = np.load(in_genotypes)
    gt_matrix = npz_data["genotypes"] # (n_variants, n_samples)
    sample_ids = npz_data["samples"]
    chroms = npz_data["chromosomes"]
    positions = npz_data["positions"]

    n_variants, n_samples = gt_matrix.shape
    log(f"Scanning ROH across {n_samples} individuals over {n_variants} markers...")

    # Load metadata
    if os.path.exists(pca_meta_path):
        meta_df = pd.read_csv(pca_meta_path, sep="\t")
    else:
        meta_df = pd.read_csv("data/metadata/sample_metadata.tsv", sep="\t")
        meta_df["Lineage"] = "Unknown"

    sample_meta_map = {row["VCF_ID"]: row for _, row in meta_df.iterrows()}

    # 2. Pre-slice chromosomes and positions for vectorized speed
    unique_chroms = np.unique(chroms)
    c_masks = {c: np.where(chroms == c)[0] for c in unique_chroms}
    c_pos = {c: positions[c_masks[c]] for c in unique_chroms}

    all_roh_records = []
    ind_summaries = []

    t0 = time.time()
    min_length_bp = 1000000
    min_snps = 50
    max_gap_bp = 1000000

    for s_idx in range(n_samples):
        s_id = str(sample_ids[s_idx])
        ind_roh = []

        for c in unique_chroms:
            c_calls = gt_matrix[c_masks[c], s_idx]
            p = c_pos[c]
            ind_roh.extend(detect_roh_chromosome(c_calls, p, chrom=c, s_id=s_id,
                                                min_length_bp=min_length_bp, min_snps=min_snps, max_gap_bp=max_gap_bp))

        all_roh_records.extend(ind_roh)

        # Compute individual metrics
        s_meta = sample_meta_map.get(s_id, {})
        breed = s_meta.get("Breed", "Unknown")
        population = s_meta.get("Population", "Unknown")
        lineage = s_meta.get("Lineage", "Unknown")
        country = s_meta.get("Country_of_origin", "Unknown")
        ho = s_meta.get("Observed_Heterozygosity_Ho", np.nan)

        if ind_roh:
            ind_df = pd.DataFrame(ind_roh)
            total_n_roh = len(ind_df)
            total_s_roh_mb = ind_df["Length_Mb"].sum()
            total_froh = total_s_roh_mb / L_AUTOSOME_MB

            froh_1_2 = ind_df[ind_df["Length_Class"] == "1-2 Mb"]["Length_Mb"].sum() / L_AUTOSOME_MB
            froh_2_4 = ind_df[ind_df["Length_Class"] == "2-4 Mb"]["Length_Mb"].sum() / L_AUTOSOME_MB
            froh_4_8 = ind_df[ind_df["Length_Class"] == "4-8 Mb"]["Length_Mb"].sum() / L_AUTOSOME_MB
            froh_8_16 = ind_df[ind_df["Length_Class"] == "8-16 Mb"]["Length_Mb"].sum() / L_AUTOSOME_MB
            froh_ge16 = ind_df[ind_df["Length_Class"] == ">=16 Mb"]["Length_Mb"].sum() / L_AUTOSOME_MB

            n_roh_1_2 = len(ind_df[ind_df["Length_Class"] == "1-2 Mb"])
            n_roh_2_4 = len(ind_df[ind_df["Length_Class"] == "2-4 Mb"])
            n_roh_4_8 = len(ind_df[ind_df["Length_Class"] == "4-8 Mb"])
            n_roh_8_16 = len(ind_df[ind_df["Length_Class"] == "8-16 Mb"])
            n_roh_ge16 = len(ind_df[ind_df["Length_Class"] == ">=16 Mb"])
        else:
            total_n_roh = 0
            total_s_roh_mb = 0.0
            total_froh = 0.0
            froh_1_2 = froh_2_4 = froh_4_8 = froh_8_16 = froh_ge16 = 0.0
            n_roh_1_2 = n_roh_2_4 = n_roh_4_8 = n_roh_8_16 = n_roh_ge16 = 0

        ind_summaries.append({
            "Sample_ID": s_id,
            "Breed": breed,
            "Population": population,
            "Lineage": lineage,
            "Country": country,
            "Observed_Heterozygosity_Ho": round(float(ho), 4) if pd.notna(ho) else np.nan,
            "Total_N_ROH": total_n_roh,
            "Total_S_ROH_Mb": round(float(total_s_roh_mb), 3),
            "F_ROH_Total": round(float(total_froh), 5),
            "F_ROH_1_2Mb": round(float(froh_1_2), 5),
            "F_ROH_2_4Mb": round(float(froh_2_4), 5),
            "F_ROH_4_8Mb": round(float(froh_4_8), 5),
            "F_ROH_8_16Mb": round(float(froh_8_16), 5),
            "F_ROH_ge16Mb": round(float(froh_ge16), 5),
            "N_ROH_1_2Mb": n_roh_1_2,
            "N_ROH_2_4Mb": n_roh_2_4,
            "N_ROH_4_8Mb": n_roh_4_8,
            "N_ROH_8_16Mb": n_roh_8_16,
            "N_ROH_ge16Mb": n_roh_ge16
        })

    elapsed = time.time() - t0
    log(f"ROH scan completed in {elapsed:.2f}s! Total segments: {len(all_roh_records)}.")

    # 3. Save ROH segments table
    roh_df = pd.DataFrame(all_roh_records)
    roh_df.to_csv(out_roh_segments, sep="\t", index=False)
    log(f"Saved {len(roh_df)} ROH segments to {out_roh_segments}.")

    # 4. Save individual inbreeding table
    ind_df = pd.DataFrame(ind_summaries)
    ind_df.to_csv(out_individual_froh, sep="\t", index=False)
    log(f"Saved individual inbreeding metrics to {out_individual_froh}.")

    # 5. Summarize by Lineage
    log("Computing summary statistics across evolutionary lineages...")
    lineage_summary = ind_df.groupby("Lineage").agg({
        "Sample_ID": "count",
        "Total_N_ROH": ["mean", "median", "std"],
        "Total_S_ROH_Mb": ["mean", "median", "std"],
        "F_ROH_Total": ["mean", "median", "std"],
        "F_ROH_1_2Mb": "mean",
        "F_ROH_2_4Mb": "mean",
        "F_ROH_4_8Mb": "mean",
        "F_ROH_8_16Mb": "mean",
        "F_ROH_ge16Mb": "mean"
    }).round(4)
    lineage_summary.to_csv(out_lineage_summary, sep="\t")
    log(f"Saved lineage summary to {out_lineage_summary}.")

    log("Lineage Inbreeding Overview (Mean F_ROH):")
    for lineage, grp in ind_df.groupby("Lineage"):
        log(f"  {lineage} (N={len(grp)}): Mean F_ROH = {grp['F_ROH_Total'].mean():.4f} +/- {grp['F_ROH_Total'].std():.4f}, Mean S_ROH = {grp['Total_S_ROH_Mb'].mean():.1f} Mb")

    log("Step 3 completed successfully.")

if __name__ == "__main__":
    main()
