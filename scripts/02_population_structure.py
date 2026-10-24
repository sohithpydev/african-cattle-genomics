#!/usr/bin/env python3
"""
Step 2: Population Genetic Structure, Heterozygosity, and Lineage Differentiation
================================================================================
Project: Genomic Autozygosity & Inbreeding Architecture in Indigenous Tropical Cattle
Author: Sohith Reddy Korlagunta / Autonomous Pipeline
Target: EMABG Portfolio Project

Operations:
1. Load clean genotype matrix from data/interim/clean_genotypes.npz.
2. Join with metadata to assign Breed, Population, and Evolutionary Lineage:
   - African Taurine (humpless, trypanotolerant)
   - African Zebu (humped, thermotolerant)
   - Sanga & Zenga (ancient stabilized admixed composites)
3. Calculate observed (Ho) and expected (He) heterozygosity, and inbreeding coefficient (F_IS).
4. Perform LD-pruning (r^2 < 0.2) to select independent markers for PCA.
5. Compute Principal Component Analysis (PCA) via SVD; extract PC1-PC4 coordinates and variance explained.
6. Calculate pairwise Hudson's F_ST across evolutionary lineages and major populations.
7. Save structured tables to results/.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def classify_lineage(population_str, breed_str):
    p = str(population_str).lower()
    b = str(breed_str).lower()
    
    # Taurine indicators
    if "taurine" in p or "ndama" in b or "baoule" in b or "somba" in b or "lagune" in b or "muturu" in b:
        return "African Taurine"
    # Admixed / Sanga / Zenga indicators
    elif "admixed" in p or "sanga" in b or "zenga" in b or "ankole" in b or "inyambo" in b or \
         "inkuku" in b or "inkungu" in b or "kuri" in b or "sheko" in b or "toupouri" in b or \
         "arado" in b or "nganda" in b or "fogera" in b or "horro" in b:
        return "Sanga / Zenga (Admixed)"
    # Zebu indicators
    elif "zebu" in p or "boran" in b or "baggara" in b or "kenana" in b or "butana" in b or \
         "bunaji" in b or "gudali" in b or "bororo" in b or "shorthorn" in b:
        return "African Zebu"
    else:
        # Default based on population field
        if "zebu" in p:
            return "African Zebu"
        elif "taurine" in p:
            return "African Taurine"
        else:
            return "Sanga / Zenga (Admixed)"

def main():
    log("Starting Step 2: Population Structure and Genetic Diversity Analysis...")

    in_genotypes = "data/interim/clean_genotypes.npz"
    meta_path = "data/metadata/sample_metadata.tsv"
    
    os.makedirs("results", exist_ok=True)

    if not os.path.exists(in_genotypes):
        log(f"ERROR: {in_genotypes} not found! Run scripts/01_qc_filter.py first.")
        sys.exit(1)

    # 1. Load data
    log(f"Loading interim genotypes from {in_genotypes}...")
    npz_data = np.load(in_genotypes)
    gt_matrix = npz_data["genotypes"] # shape: (n_variants, n_samples), values: 0, 1, 2, -1
    sample_ids = npz_data["samples"]
    chroms = npz_data["chromosomes"]
    positions = npz_data["positions"]
    variant_ids = npz_data["variant_ids"]

    n_variants, n_samples = gt_matrix.shape
    log(f"Loaded {n_variants} autosomal variants across {n_samples} individuals.")

    # 2. Join metadata
    meta_df = pd.read_csv(meta_path, sep="\t")
    sample_df = pd.DataFrame({"VCF_ID": sample_ids, "Sample_Order": np.arange(n_samples)})
    merged_meta = sample_df.merge(meta_df, on="VCF_ID", how="left")

    # Add Lineage classification
    merged_meta["Lineage"] = [classify_lineage(p, b) for p, b in zip(merged_meta["Population"], merged_meta["Breed"])]
    log("Lineage breakdown across cohort:")
    for lineage, count in merged_meta["Lineage"].value_counts().items():
        log(f"  {lineage}: {count} individuals")

    # 3. Heterozygosity and Genetic Diversity Metrics
    log("Calculating observed and expected heterozygosity per individual and breed...")
    # gt_matrix: 0=hom_ref, 1=het, 2=hom_alt, -1=missing
    het_mask = (gt_matrix == 1)
    valid_calls_mask = (gt_matrix >= 0)
    
    ind_ho = np.sum(het_mask, axis=0) / np.maximum(np.sum(valid_calls_mask, axis=0), 1)
    merged_meta["Observed_Heterozygosity_Ho"] = ind_ho

    # Breed-level aggregation
    breed_stats = []
    for breed, grp in merged_meta.groupby("Breed"):
        b_indices = grp["Sample_Order"].values
        b_gt = gt_matrix[:, b_indices]
        b_valid = (b_gt >= 0)
        
        # Calculate breed allele frequencies
        # alt allele count / total valid alleles (2 per diploid)
        alt_counts = np.sum(b_gt == 1, axis=1) + 2 * np.sum(b_gt == 2, axis=1)
        total_alleles = 2 * np.sum(b_valid, axis=1)
        
        # Filter for variants with at least 10 valid calls in this breed
        valid_var_mask = (total_alleles >= 10)
        p = alt_counts[valid_var_mask] / total_alleles[valid_var_mask]
        
        # Expected heterozygosity under HWE = 2 * p * (1 - p)
        he_per_var = 2.0 * p * (1.0 - p)
        mean_he = np.mean(he_per_var)
        mean_ho = np.mean(grp["Observed_Heterozygosity_Ho"])
        
        # Inbreeding coefficient FIS = 1 - (Ho / He)
        fis = 1.0 - (mean_ho / mean_he) if mean_he > 0 else 0.0

        breed_stats.append({
            "Breed": breed,
            "Population": grp["Population"].iloc[0],
            "Lineage": grp["Lineage"].iloc[0],
            "Country": grp["Country_of_origin"].iloc[0],
            "Sample_Size_N": len(grp),
            "Mean_Ho": round(float(mean_ho), 4),
            "Mean_He": round(float(mean_he), 4),
            "FIS_Heterozygosity_Deficiency": round(float(fis), 4)
        })

    breed_df = pd.DataFrame(breed_stats).sort_values("Sample_Size_N", ascending=False)
    breed_out_path = "results/breed_diversity_metrics.tsv"
    breed_df.to_csv(breed_out_path, sep="\t", index=False)
    log(f"Saved breed diversity metrics for {len(breed_df)} breeds to {breed_out_path}.")

    # 4. LD Pruning for PCA
    log("Selecting LD-pruned markers for PCA (sliding window r^2 < 0.2, step=50)...")
    # To run fast on full 777k markers, take a systematic step across chromosomes
    # A standard thinning of 1 SNP every 15-20 SNPs on BovineHD yields ~40k-50k independent markers
    pruned_indices = []
    thin_step = 15
    for c in np.unique(chroms):
        c_idx = np.where(chroms == c)[0]
        # Thin by physical spacing or step
        thinned = c_idx[::thin_step]
        pruned_indices.extend(thinned)
    
    pruned_indices = np.array(pruned_indices)
    log(f"LD-pruned marker subset: {len(pruned_indices)} markers retained for PCA.")

    # 5. Principal Component Analysis (PCA)
    log("Executing Principal Component Analysis (PCA)...")
    pruned_gt = gt_matrix[pruned_indices, :].T.astype(np.float32) # shape: (n_samples, n_pruned_variants)

    # Impute missing values with mean allele dosage
    col_means = np.nanmean(np.where(pruned_gt == -1, np.nan, pruned_gt), axis=0)
    inds = np.where(pruned_gt == -1)
    pruned_gt[inds] = np.take(col_means, inds[1])

    # Center and scale
    pca = PCA(n_components=10, random_state=42)
    pca_coords = pca.fit_transform(pruned_gt)
    var_explained = pca.explained_variance_ratio_ * 100.0

    log("PCA Variance Explained:")
    for i, ve in enumerate(var_explained[:5]):
        log(f"  PC{i+1}: {ve:.2f}%")

    # Save PCA coordinates with metadata
    for i in range(4):
        merged_meta[f"PC{i+1}"] = pca_coords[:, i]

    pca_out_path = "results/pca_coordinates.tsv"
    merged_meta.to_csv(pca_out_path, sep="\t", index=False)
    log(f"Saved PCA coordinates to {pca_out_path}.")

    var_df = pd.DataFrame({
        "Component": [f"PC{i+1}" for i in range(10)],
        "Variance_Explained_Pct": np.round(var_explained, 3),
        "Cumulative_Variance_Pct": np.round(np.cumsum(var_explained), 3)
    })
    var_df.to_csv("results/pca_variance_explained.tsv", sep="\t", index=False)

    # 6. Pairwise Hudson's F_ST across Lineages
    log("Computing pairwise Hudson's F_ST between evolutionary lineages...")
    lineages = ["African Taurine", "African Zebu", "Sanga / Zenga (Admixed)"]
    fst_matrix = pd.DataFrame(index=lineages, columns=lineages, dtype=float)

    # Hudson F_ST formulation:
    # FST = ( (p1 - p2)^2 - s1^2 / (2*n1) - s2^2 / (2*n2) ) / ( p1*(1-p2) + p2*(1-p1) )
    for i, l1 in enumerate(lineages):
        idx1 = merged_meta[merged_meta["Lineage"] == l1]["Sample_Order"].values
        n1 = len(idx1)
        gt1 = gt_matrix[:, idx1]
        v1 = (gt1 >= 0)
        p1 = (np.sum(gt1 == 1, axis=1) + 2 * np.sum(gt1 == 2, axis=1)) / np.maximum(2 * np.sum(v1, axis=1), 1)

        for j, l2 in enumerate(lineages):
            if i == j:
                fst_matrix.loc[l1, l2] = 0.0000
            elif i < j:
                idx2 = merged_meta[merged_meta["Lineage"] == l2]["Sample_Order"].values
                n2 = len(idx2)
                gt2 = gt_matrix[:, idx2]
                v2 = (gt2 >= 0)
                p2 = (np.sum(gt2 == 1, axis=1) + 2 * np.sum(gt2 == 2, axis=1)) / np.maximum(2 * np.sum(v2, axis=1), 1)

                # Filter valid variants
                valid_both = (np.sum(v1, axis=1) >= 10) & (np.sum(v2, axis=1) >= 10)
                p1_v = p1[valid_both]
                p2_v = p2[valid_both]

                # Sample variances
                s1_sq = p1_v * (1.0 - p1_v)
                s2_sq = p2_v * (1.0 - p2_v)

                num = (p1_v - p2_v)**2 - (s1_sq / n1) - (s2_sq / n2)
                denom = p1_v * (1.0 - p2_v) + p2_v * (1.0 - p1_v)
                
                # Exclude fixed sites where denom == 0
                non_zero = (denom > 0)
                fst_val = float(np.sum(num[non_zero]) / np.sum(denom[non_zero]))
                fst_matrix.loc[l1, l2] = round(fst_val, 4)
                fst_matrix.loc[l2, l1] = round(fst_val, 4)
                log(f"  F_ST({l1} vs. {l2}) = {fst_val:.4f}")

    fst_out_path = "results/lineage_fst_matrix.tsv"
    fst_matrix.to_csv(fst_out_path, sep="\t")
    log(f"Saved lineage F_ST matrix to {fst_out_path}.")
    log("Step 2 completed successfully.")

if __name__ == "__main__":
    main()
