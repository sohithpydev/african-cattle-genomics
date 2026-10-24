#!/usr/bin/env python3
"""
Step 5: Publication-Quality Figures & Tables Generation
======================================================
Project: Genomic Autozygosity & Inbreeding Architecture in Indigenous Tropical Cattle
Author: Sohith Reddy Korlagunta / Autonomous Pipeline
Target: EMABG Portfolio Project

Outputs:
- Figures (vector PDF and high-res 300-DPI PNG):
  * Figure 1: Study Cohort Provenance and Analytical Architecture
  * Figure 2: Population Genomic Structure (PCA) and Lineage Divergence (F_ST)
  * Figure 3: Genomic Autozygosity (F_ROH) Distributions and Heterozygosity Dynamics
  * Figure 4: Temporal Inbreeding Deconvolution Across 5 ROH Length Classes
  * Figure 5: Genome-Wide Autozygosity Islands and Functional Adaptation Hotspots
- Tables (CSV format):
  * Table 1: Comprehensive Breed Genetic Diversity and Inbreeding Metrics
  * Table 2: Temporal Inbreeding (F_ROH_k) Stratification by Lineage
  * Table 3: Fine-Mapped Autozygosity Islands and Overlapping Candidate Genes
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# Set publication style
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["axes.edgecolor"] = "#2c3e50"
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["xtick.color"] = "#2c3e50"
plt.rcParams["ytick.color"] = "#2c3e50"
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 9
plt.rcParams["legend.fontsize"] = 8.5

# Curated scientific color palette
LINEAGE_COLORS = {
    "African Taurine": "#2b5c8f",          # Deep slate blue
    "African Zebu": "#d95f02",             # Warm burnt orange
    "Sanga / Zenga (Admixed)": "#7570b3"  # Muted royal purple
}

def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def generate_figure_1():
    log("Generating Figure 1: Study Cohort Provenance & Workflow...")
    meta_df = pd.read_csv("data/metadata/sample_metadata.tsv", sep="\t")
    included = meta_df[meta_df["Study_Cohort"] == "INCLUDED"].copy()

    def classify_lineage(population_str, breed_str):
        p = str(population_str).lower()
        b = str(breed_str).lower()
        if "taurine" in p or "ndama" in b or "baoule" in b or "somba" in b or "lagune" in b:
            return "African Taurine"
        elif "admixed" in p or "sanga" in b or "zenga" in b or "ankole" in b or "inyambo" in b or \
             "inkuku" in b or "inkungu" in b or "kuri" in b or "sheko" in b or "toupouri" in b or \
             "arado" in b or "nganda" in b or "fogera" in b or "horro" in b:
            return "Sanga / Zenga (Admixed)"
        elif "zebu" in p or "boran" in b or "baggara" in b or "kenana" in b or "butana" in b or \
             "bunaji" in b or "gudali" in b or "bororo" in b or "shorthorn" in b:
            return "African Zebu"
        else:
            return "African Zebu" if "zebu" in p else ("African Taurine" if "taurine" in p else "Sanga / Zenga (Admixed)")

    included["Lineage"] = [classify_lineage(p, b) for p, b in zip(included["Population"], included["Breed"])]

    fig = plt.figure(figsize=(11, 7.5), dpi=300)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.2, 1.0], height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.3)

    # Panel A: Geographic Distribution by Country
    ax1 = fig.add_subplot(gs[0, 0])
    country_counts = included["Country_of_origin"].value_counts()
    colors = ["#1f77b4" if c in ["Guinea", "Senegal"] else "#ff7f0e" if c in ["Nigeria", "Chad"] else "#2ca02c" for c in country_counts.index]
    bars = ax1.barh(country_counts.index[::-1], country_counts.values[::-1], color="#34495e", edgecolor="#1a252f", height=0.65)
    ax1.set_xlabel("Sample Size (N)")
    ax1.set_title("A. Geographic Cohort Representation (10 African Nations)", fontweight="bold", loc="left")
    ax1.grid(axis="x", linestyle=":", alpha=0.6)
    for bar in bars:
        w = bar.get_width()
        ax1.text(w + 3, bar.get_y() + bar.get_height()/2, f"{int(w)}", va="center", ha="left", fontsize=8)

    # Panel B: Lineage Composition Pie/Donut
    ax2 = fig.add_subplot(gs[0, 1])
    lin_counts = included["Lineage"].value_counts()
    pie_colors = [LINEAGE_COLORS.get(k, "#7f7f7f") for k in lin_counts.index]
    wedges, texts, autotexts = ax2.pie(
        lin_counts.values,
        labels=[f"{k}\n(N={v})" for k, v in zip(lin_counts.index, lin_counts.values)],
        autopct="%1.1f%%",
        startangle=140,
        colors=pie_colors,
        wedgeprops=dict(width=0.45, edgecolor="white", linewidth=2),
        pctdistance=0.75
    )
    for at in autotexts:
        at.set_fontsize(8.5)
        at.set_weight("bold")
    ax2.set_title("B. Evolutionary Lineage Breakdown", fontweight="bold", loc="left")

    # Panel C: Top 15 Breeds Representation
    ax3 = fig.add_subplot(gs[1, :])
    top_breeds = included["Breed"].value_counts().head(18)
    breed_lineages = [included[included["Breed"] == b]["Lineage"].iloc[0] for b in top_breeds.index]
    b_colors = [LINEAGE_COLORS.get(l, "#7f7f7f") for l in breed_lineages]
    
    bars2 = ax3.bar(range(len(top_breeds)), top_breeds.values, color=b_colors, edgecolor="#2c3e50", width=0.65)
    ax3.set_xticks(range(len(top_breeds)))
    ax3.set_xticklabels(top_breeds.index, rotation=40, ha="right", fontsize=8)
    ax3.set_ylabel("Number of Genotyped Animals")
    ax3.set_title("C. Sample Counts Across Major Indigenous Cattle Breeds", fontweight="bold", loc="left")
    ax3.grid(axis="y", linestyle=":", alpha=0.6)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=l) for l, c in LINEAGE_COLORS.items()]
    ax3.legend(handles=legend_patches, loc="upper right", framealpha=0.9)

    plt.tight_layout()
    fig.savefig("figures/figure_1_study_design_provenance.pdf", bbox_inches="tight")
    fig.savefig("figures/figure_1_study_design_provenance.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Saved Figure 1.")

def generate_figure_2():
    log("Generating Figure 2: Population Genomic Structure (PCA) & F_ST...")
    pca_df = pd.read_csv("results/pca_coordinates.tsv", sep="\t")
    var_df = pd.read_csv("results/pca_variance_explained.tsv", sep="\t")
    fst_df = pd.read_csv("results/lineage_fst_matrix.tsv", sep="\t", index_col=0)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10.5, 9), dpi=300)

    # Panel A: PC1 vs PC2
    for lineage, grp in pca_df.groupby("Lineage"):
        ax1.scatter(
            grp["PC1"], grp["PC2"],
            c=LINEAGE_COLORS.get(lineage, "#7f7f7f"),
            label=lineage,
            alpha=0.7,
            edgecolors="none",
            s=28
        )
    ax1.set_xlabel(f"PC1 ({var_df['Variance_Explained_Pct'].iloc[0]:.2f}%)")
    ax1.set_ylabel(f"PC2 ({var_df['Variance_Explained_Pct'].iloc[1]:.2f}%)")
    ax1.set_title("A. Lineage Divergence: PC1 vs PC2", fontweight="bold", loc="left")
    ax1.legend(loc="upper right", framealpha=0.85)
    ax1.grid(linestyle=":", alpha=0.5)

    # Panel B: PC2 vs PC3
    for lineage, grp in pca_df.groupby("Lineage"):
        ax2.scatter(
            grp["PC2"], grp["PC3"],
            c=LINEAGE_COLORS.get(lineage, "#7f7f7f"),
            label=lineage,
            alpha=0.7,
            edgecolors="none",
            s=28
        )
    ax2.set_xlabel(f"PC2 ({var_df['Variance_Explained_Pct'].iloc[1]:.2f}%)")
    ax2.set_ylabel(f"PC3 ({var_df['Variance_Explained_Pct'].iloc[2]:.2f}%)")
    ax2.set_title("B. Regional Sub-Structure: PC2 vs PC3", fontweight="bold", loc="left")
    ax2.grid(linestyle=":", alpha=0.5)

    # Panel C: Scree Plot
    pcs = range(1, len(var_df) + 1)
    ax3.plot(pcs, var_df["Variance_Explained_Pct"], marker="o", color="#2c3e50", lw=1.5, ms=5)
    ax3.set_xticks(list(pcs))
    ax3.set_xlabel("Principal Component")
    ax3.set_ylabel("Variance Explained (%)")
    ax3.set_title("C. Eigenvalue Scree Plot (Top 10 PCs)", fontweight="bold", loc="left")
    ax3.grid(linestyle=":", alpha=0.5)

    # Panel D: Pairwise F_ST Heatmap
    im = ax4.imshow(fst_df.values, cmap="Blues", vmin=0.0, vmax=0.15)
    ax4.set_xticks(range(len(fst_df.columns)))
    ax4.set_yticks(range(len(fst_df.index)))
    short_labels = ["Taurine", "Zebu", "Sanga/Zenga"]
    ax4.set_xticklabels(short_labels, fontsize=9)
    ax4.set_yticklabels(short_labels, fontsize=9)
    ax4.set_title("D. Pairwise Hudson's F_ST Matrix", fontweight="bold", loc="left")
    cbar = fig.colorbar(im, ax=ax4, shrink=0.75, pad=0.04)
    cbar.set_label("Genetic Differentiation (F_ST)", fontsize=8.5)

    # Annotate values
    for i in range(len(fst_df.index)):
        for j in range(len(fst_df.columns)):
            val = fst_df.iloc[i, j]
            ax4.text(j, i, f"{val:.4f}", ha="center", va="center", color="white" if val > 0.08 else "black", fontsize=9.5, fontweight="bold")

    plt.tight_layout()
    fig.savefig("figures/figure_2_population_structure_pca.pdf", bbox_inches="tight")
    fig.savefig("figures/figure_2_population_structure_pca.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Saved Figure 2.")

def generate_figure_3():
    log("Generating Figure 3: F_ROH Distributions & Heterozygosity Correlation...")
    ind_df = pd.read_csv("results/individual_inbreeding_froh.tsv", sep="\t")
    breed_df = pd.read_csv("results/breed_diversity_metrics.tsv", sep="\t")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 9), dpi=300)

    # Panel A: Boxplot of F_ROH across Lineages
    lineages = ["African Taurine", "African Zebu", "Sanga / Zenga (Admixed)"]
    froh_by_lineage = [ind_df[ind_df["Lineage"] == l]["F_ROH_Total"].values for l in lineages]
    bp = ax1.boxplot(froh_by_lineage, tick_labels=["Taurine", "Zebu", "Sanga/Zenga"], patch_artist=True, widths=0.5)
    for patch, l in zip(bp["boxes"], lineages):
        patch.set_facecolor(LINEAGE_COLORS.get(l, "#7f7f7f"))
        patch.set_alpha(0.7)
    for median in bp["medians"]:
        median.set(color="black", linewidth=1.5)
    ax1.set_ylabel("Genomic Inbreeding Coefficient (F_ROH)")
    ax1.set_title("A. Inbreeding Load (F_ROH) by Evolutionary Lineage", fontweight="bold", loc="left")
    ax1.grid(axis="y", linestyle=":", alpha=0.5)

    # Panel B: Scatter F_ROH vs Observed Heterozygosity (Ho)
    for l in lineages:
        sub = ind_df[ind_df["Lineage"] == l]
        ax2.scatter(sub["Observed_Heterozygosity_Ho"], sub["F_ROH_Total"], c=LINEAGE_COLORS[l], label=l, alpha=0.6, s=25)
    
    # Compute correlation
    valid_mask = pd.notna(ind_df["Observed_Heterozygosity_Ho"]) & pd.notna(ind_df["F_ROH_Total"])
    r_val = np.corrcoef(ind_df.loc[valid_mask, "Observed_Heterozygosity_Ho"], ind_df.loc[valid_mask, "F_ROH_Total"])[0, 1]
    
    # Trendline
    m, b = np.polyfit(ind_df.loc[valid_mask, "Observed_Heterozygosity_Ho"], ind_df.loc[valid_mask, "F_ROH_Total"], 1)
    x_seq = np.linspace(ind_df["Observed_Heterozygosity_Ho"].min(), ind_df["Observed_Heterozygosity_Ho"].max(), 50)
    ax2.plot(x_seq, m*x_seq + b, color="#2c3e50", linestyle="--", lw=1.5, label=f"Trend (r = {r_val:.3f})")

    ax2.set_xlabel("Observed Heterozygosity (H_o)")
    ax2.set_ylabel("F_ROH")
    ax2.set_title(f"B. Autozygosity vs Heterozygosity (Pearson r = {r_val:.3f})", fontweight="bold", loc="left")
    ax2.legend(loc="upper right", framealpha=0.85)
    ax2.grid(linestyle=":", alpha=0.5)

    # Panel C: Total S_ROH (Mb) vs N_ROH (Segment Count)
    for l in lineages:
        sub = ind_df[ind_df["Lineage"] == l]
        ax3.scatter(sub["Total_N_ROH"], sub["Total_S_ROH_Mb"], c=LINEAGE_COLORS[l], label=l, alpha=0.6, s=25)
    ax3.set_xlabel("Number of Runs of Homozygosity (N_ROH)")
    ax3.set_ylabel("Total ROH Length (S_ROH, Mb)")
    ax3.set_title("C. Autozygosity Coverage vs Segment Count", fontweight="bold", loc="left")
    ax3.grid(linestyle=":", alpha=0.5)

    # Panel D: Breed Mean F_ROH Ranking (Top 15 Breeds)
    breed_means = ind_df.groupby("Breed").agg({"F_ROH_Total": ["mean", "std"], "Lineage": "first", "Sample_ID": "count"})
    breed_means.columns = ["mean_froh", "std_froh", "lineage", "N"]
    top_b = breed_means[breed_means["N"] >= 15].sort_values("mean_froh", ascending=True)

    y_pos = range(len(top_b))
    colors_d = [LINEAGE_COLORS.get(l, "#7f7f7f") for l in top_b["lineage"]]
    ax4.barh(y_pos, top_b["mean_froh"], xerr=top_b["std_froh"], color=colors_d, alpha=0.8, edgecolor="#2c3e50", height=0.65, capsize=3)
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(top_b.index, fontsize=8)
    ax4.set_xlabel("Mean F_ROH (+/- 1 SD)")
    ax4.set_title("D. Mean Inbreeding Across Well-Sampled Breeds (N >= 15)", fontweight="bold", loc="left")
    ax4.grid(axis="x", linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig("figures/figure_3_froh_inbreeding_distributions.pdf", bbox_inches="tight")
    fig.savefig("figures/figure_3_froh_inbreeding_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Saved Figure 3.")

def generate_figure_4():
    log("Generating Figure 4: Temporal Deconvolution Across ROH Length Classes...")
    ind_df = pd.read_csv("results/individual_inbreeding_froh.tsv", sep="\t")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), dpi=300, gridspec_kw={"width_ratios": [1.1, 1.0]})

    lineages = ["African Taurine", "African Zebu", "Sanga / Zenga (Admixed)"]
    bins = ["F_ROH_1_2Mb", "F_ROH_2_4Mb", "F_ROH_4_8Mb", "F_ROH_8_16Mb", "F_ROH_ge16Mb"]
    bin_labels = ["1-2 Mb (Ancient)", "2-4 Mb (Historical)", "4-8 Mb (Intermediate)", "8-16 Mb (Recent)", ">=16 Mb (Consanguineous)"]
    bin_colors = ["#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"]

    # Panel A: Stacked Bar Chart of Absolute F_ROH by Lineage
    means_by_lin = []
    for l in lineages:
        sub = ind_df[ind_df["Lineage"] == l]
        means_by_lin.append([sub[b].mean() for b in bins])
    means_by_lin = np.array(means_by_lin) # (3, 5)

    bottoms = np.zeros(len(lineages))
    for b_idx in range(len(bins)):
        vals = means_by_lin[:, b_idx]
        ax1.bar(["Taurine", "Zebu", "Sanga/Zenga"], vals, bottom=bottoms, color=bin_colors[b_idx], label=bin_labels[b_idx], edgecolor="#2c3e50", width=0.55)
        bottoms += vals

    ax1.set_ylabel("Cumulative Mean F_ROH")
    ax1.set_title("A. Stratified Inbreeding Load Across 5 Length Classes", fontweight="bold", loc="left")
    ax1.legend(loc="upper right", framealpha=0.9, fontsize=8)
    ax1.grid(axis="y", linestyle=":", alpha=0.5)

    # Panel B: Relative Proportion (%) of Autozygosity
    rel_props = (means_by_lin.T / np.sum(means_by_lin, axis=1) * 100.0).T # (3, 5)
    bottoms2 = np.zeros(len(lineages))
    for b_idx in range(len(bins)):
        vals2 = rel_props[:, b_idx]
        ax2.bar(["Taurine", "Zebu", "Sanga/Zenga"], vals2, bottom=bottoms2, color=bin_colors[b_idx], edgecolor="#2c3e50", width=0.55)
        bottoms2 += vals2

    ax2.set_ylabel("Proportion of Total Autozygosity (%)")
    ax2.set_title("B. Temporal Breakdown of Autozygous Genome", fontweight="bold", loc="left")
    ax2.set_ylim(0, 100)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig("figures/figure_4_temporal_roh_length_stratification.pdf", bbox_inches="tight")
    fig.savefig("figures/figure_4_temporal_roh_length_stratification.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Saved Figure 4.")

def generate_figure_5():
    log("Generating Figure 5: Manhattan Plot of Autozygosity Islands...")
    freq_df = pd.read_csv("results/genome_wide_roh_frequency.tsv", sep="\t")
    islands_df = pd.read_csv("results/autozygosity_islands.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(12, 5.5), dpi=300)

    # Calculate cumulative positions for Manhattan plotting
    chrom_lengths = {}
    for c in range(1, 30):
        c_sub = freq_df[freq_df["Chrom"] == c]
        chrom_lengths[c] = c_sub["End_Pos"].max() if len(c_sub) > 0 else 0

    chrom_offsets = {}
    curr_offset = 0
    for c in range(1, 30):
        chrom_offsets[c] = curr_offset
        curr_offset += chrom_lengths[c]

    freq_df["Plot_Pos"] = freq_df["Mid_Pos"] + freq_df["Chrom"].map(chrom_offsets)

    # Alternating colors
    chr_colors = ["#2b5c8f", "#7570b3"]
    x_ticks = []
    x_labels = []

    for c in range(1, 30):
        c_df = freq_df[freq_df["Chrom"] == c]
        if len(c_df) == 0:
            continue
        c_col = chr_colors[c % 2]
        ax.plot(c_df["Plot_Pos"] / 1e6, c_df["P_ROH_Total_Pct"], color=c_col, lw=0.7, alpha=0.85)
        x_ticks.append((c_df["Plot_Pos"].min() + c_df["Plot_Pos"].max()) / 2e6)
        x_labels.append(str(c))

    # Top 1% threshold line
    thresh = np.percentile(freq_df["P_ROH_Total_Pct"], 99.0)
    ax.axhline(thresh, color="#d7191c", linestyle="--", lw=1.2, label=f"Top 1% Empirical Threshold ({thresh:.2f}%)")

    # Annotate top candidate gene islands
    key_annotations = [
        {"Gene": "HSP90AA1", "Chrom": 28, "Pos": 40.54, "Text_Offset": (0, 8)},
        {"Gene": "PRKG1", "Chrom": 26, "Pos": 8.85, "Text_Offset": (0, 8)},
        {"Gene": "TICAM1", "Chrom": 7, "Pos": 42.11, "Text_Offset": (0, 8)},
        {"Gene": "ARHGAP15", "Chrom": 2, "Pos": 107.7, "Text_Offset": (0, 8)},
        {"Gene": "PLAG1", "Chrom": 14, "Pos": 23.35, "Text_Offset": (0, 8)},
        {"Gene": "KIT", "Chrom": 6, "Pos": 70.15, "Text_Offset": (0, 8)}
    ]

    for ann in key_annotations:
        c = ann["Chrom"]
        pos_mb = ann["Pos"]
        if c in chrom_offsets:
            abs_pos_mb = (chrom_offsets[c] + pos_mb * 1e6) / 1e6
            # Find frequency near this position
            c_near = freq_df[(freq_df["Chrom"] == c) & (np.abs(freq_df["Mid_Pos"]/1e6 - pos_mb) < 0.5)]
            if len(c_near) > 0:
                y_val = c_near["P_ROH_Total_Pct"].max()
                ax.scatter(abs_pos_mb, y_val, color="#d7191c", s=45, zorder=5, edgecolor="black")
                ax.annotate(
                    ann["Gene"],
                    xy=(abs_pos_mb, y_val),
                    xytext=(abs_pos_mb + ann["Text_Offset"][0], y_val + ann["Text_Offset"][1]),
                    arrowprops=dict(facecolor="black", arrowstyle="->", lw=0.8),
                    fontsize=8,
                    fontweight="bold"
                )

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=7.5)
    ax.set_xlabel("Bovine Autosomes (BTA 1 to 29 on ARS-UCD1.2)")
    ax.set_ylabel("Population Autozygosity Frequency (P_ROH, %)")
    ax.set_title("Genome-Wide Distribution of Autozygosity Islands & Adaptive Loci", fontweight="bold", loc="left")
    ax.set_ylim(0, max(freq_df["P_ROH_Total_Pct"].max() * 1.25, 40))
    ax.legend(loc="upper right", framealpha=0.9)
    ax.grid(axis="y", linestyle=":", alpha=0.5)

    plt.tight_layout()
    fig.savefig("figures/figure_5_autozygosity_islands_manhattan.pdf", bbox_inches="tight")
    fig.savefig("figures/figure_5_autozygosity_islands_manhattan.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("Saved Figure 5.")

def generate_tables():
    log("Generating publication tables...")
    os.makedirs("tables", exist_ok=True)

    # Table 1: Breed summary
    breed_df = pd.read_csv("results/breed_diversity_metrics.tsv", sep="\t")
    ind_df = pd.read_csv("results/individual_inbreeding_froh.tsv", sep="\t")

    # Add F_ROH mean and SD to breed table
    breed_froh = ind_df.groupby("Breed").agg({
        "F_ROH_Total": ["mean", "std"],
        "Total_N_ROH": "mean",
        "Total_S_ROH_Mb": "mean"
    })
    breed_froh.columns = ["Mean_F_ROH", "SD_F_ROH", "Mean_N_ROH", "Mean_S_ROH_Mb"]
    table1 = breed_df.merge(breed_froh.round(4), on="Breed", how="left")
    table1.to_csv("tables/table_1_breed_diversity_summary.csv", index=False)
    log("Saved tables/table_1_breed_diversity_summary.csv.")

    # Table 2: Lineage temporal stratification
    lineage_table = ind_df.groupby("Lineage").agg({
        "Sample_ID": "count",
        "F_ROH_Total": ["mean", "std"],
        "F_ROH_1_2Mb": ["mean", "std"],
        "F_ROH_2_4Mb": ["mean", "std"],
        "F_ROH_4_8Mb": ["mean", "std"],
        "F_ROH_8_16Mb": ["mean", "std"],
        "F_ROH_ge16Mb": ["mean", "std"],
        "Observed_Heterozygosity_Ho": ["mean", "std"]
    }).round(4)
    lineage_table.to_csv("tables/table_2_lineage_froh_stratification.csv")
    log("Saved tables/table_2_lineage_froh_stratification.csv.")

    # Table 3: Top autozygosity islands
    islands_df = pd.read_csv("results/autozygosity_islands.tsv", sep="\t")
    top_islands = islands_df.head(20)
    top_islands.to_csv("tables/table_3_top_autozygosity_islands.csv", index=False)
    log("Saved tables/table_3_top_autozygosity_islands.csv.")

def main():
    log("Starting Step 5: Generating publication figures and tables...")
    os.makedirs("figures", exist_ok=True)
    os.makedirs("tables", exist_ok=True)

    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    generate_figure_4()
    generate_figure_5()
    generate_tables()
    log("Step 5 completed successfully.")

if __name__ == "__main__":
    main()
