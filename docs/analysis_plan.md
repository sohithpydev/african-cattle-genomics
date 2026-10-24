# Local Preregistered Analysis Plan & Analytical Protocol

**Protocol Identifier:** AP-EMABG-2026-001  
**Date of Preregistration:** October 24, 2026  
**Status:** Frozen Prior to Outcome-Level Exploration  
**Lead Analyst / Applicant:** Sohith Reddy Korlagunta  
**Repository:** `/Users/sohith/Desktop/emabg-archive`  

---

## 1. Research Question and Primary Estimands

### Primary Question:
How does the architecture of genomic autozygosity (runs of homozygosity segment count, length distribution, and cumulative coverage) differ across indigenous African Taurine, Zebu, and ancient admixed (Sanga/Zenga) cattle lineages?

### Primary Estimand:
The individual genomic inbreeding coefficient derived from runs of homozygosity ($F_{ROH}$), defined as:
$$F_{ROH} = \frac{\sum L_{ROH}}{L_{auto}}$$
where $\sum L_{ROH}$ is the total physical length (in Mb) of all identified runs of homozygosity in an individual, and $L_{auto}$ is the total physical length of the autosomal genome covered by markers (for Hereford assembly `ARS-UCD1.2`, $L_{auto} \approx 2,489.38\text{ Mb}$).

### Stratified Estimands ($F_{ROH, k}$ across 5 length classes):
Recombination breaks ancestral chromosome blocks over successive meiotic generations. The expected length of an IBD segment follows an exponential distribution with mean $100 / (2g)\text{ cM}$, where $g$ is the number of generations since the common ancestor. Assuming an average bovine recombination rate of $1\text{ cM} \approx 1\text{ Mb}$:
1. **$F_{ROH, 1-2\text{Mb}}$:** Ancient autozygosity ($g \approx 25$–$50$ generations ago; ~150–300 years ago).
2. **$F_{ROH, 2-4\text{Mb}}$:** Historical background autozygosity ($g \approx 12$–$25$ generations ago; ~70–150 years ago).
3. **$F_{ROH, 4-8\text{Mb}}$:** Intermediate autozygosity ($g \approx 6$–$12$ generations ago; ~35–70 years ago).
4. **$F_{ROH, 8-16\text{Mb}}$:** Recent inbreeding ($g \approx 3$–$6$ generations ago; ~18–35 years ago).
5. **$F_{ROH, \ge 16\text{Mb}}$:** Very recent consanguinity ($g \le 3$ generations ago; within past 1–3 generations).

### Secondary Estimand:
- **Autozygosity Islands (ROH Hotspots):** Marker-by-marker population homozygosity frequency ($P_{ROH}$), with islands defined as contiguous genomic intervals where $P_{ROH}$ exceeds the empirical top 1% threshold across the lineage genome.

---

## 2. Inclusion, Exclusion, and Quality Control Gates

### Sample-Level QC:
1. Retain only the $N=1,063$ individuals passing initial provider QC in Table S4.
2. Individual genotype completion rate $\ge 95\%$ (missingness $F_{MISS} \le 0.05$).
3. First-degree relatives or duplicates ($PHI \ge 0.498$) identified in Table S6 are strictly excluded.
4. All sample removals recorded with reason in `results/qc_exclusions.tsv`.

### Variant-Level QC:
1. Autosomes only: Retain chromosomes BTA 1 to BTA 29. Exclude sex chromosomes, mitochondrial DNA, and unplaced contigs.
2. Genotype call rate $\ge 90\%$ per variant across the cohort.
3. Minor Allele Frequency (MAF):
   - **For ROH detection:** No MAF pruning (filtering on MAF introduces ascertainment distortion, eliminates rare homozygous alleles, and artificially truncates runs of homozygosity).
   - **For PCA and $F_{ST}$:** MAF $\ge 0.05$ and linkage disequilibrium (LD) pruning ($r^2 < 0.2$ in 50-SNP sliding windows with 10-SNP steps) to prevent LD-driven distortion of eigenvectors.

---

## 3. Computational Algorithms & Implementation Parameters

### Runs of Homozygosity (ROH) Detection:
ROH will be identified using a rigorous sliding-window algorithm and consecutive homozygous scanning implemented via standardized, deterministic tools (PLINK / `scikit-allel` / R):
- Minimum physical length: $1,000\text{ kb}$ ($1.0\text{ Mb}$).
- Minimum number of consecutive homozygous SNPs: $50$ SNPs.
- Maximum gap between consecutive homozygous SNPs: $1,000\text{ kb}$ ($1.0\text{ Mb}$ to prevent false joins across centromeric or unmapped gaps, following Nandolo et al. 2018).
- Density threshold: at least 1 SNP per $50\text{ kb}$.
- Allowed heterozygous calls per run:
  - For segments $<4\text{ Mb}$: 0 heterozygous SNPs permitted.
  - For segments $\ge 4\text{ Mb}$: maximum 1 heterozygous SNP permitted (to account for sequencing/genotyping error or somatic mutation).
- Allowed missing calls: maximum 1 missing call per run.

### Population Structure & Fixation Index ($F_{ST}$):
- Principal Component Analysis (PCA) performed on LD-pruned, centered, and scaled genotype matrices.
- Pairwise Hudson's $F_{ST}$ calculated across the three evolutionary lineages (Taurine, Zebu, Sanga/Zenga) with jackknife standard errors.
- Observed ($H_o$) and expected ($H_e$) heterozygosity calculated per breed.

---

## 4. Sensitivity Analyses & Robustness Checks

To ensure all findings are robust and not artifacts of parameter tuning, the following sensitivity checks are prespecified:
1. **Window Size Sensitivity:** Benchmark ROH detection across minimum SNP count thresholds ($N_{SNP} \in \{30, 50, 75\}$).
2. **Heterozygous Tolerance Check:** Compare $F_{ROH}$ computed with 0 vs. 1 allowed heterozygous call across all length bins.
3. **Sample Imbalance Adjustment:** Re-evaluate breed-level means using balanced subsampling ($N=20$ per representative breed) to confirm that large-cohort breeds (e.g., N'Dama $N=140$, EASZ $N=114$) do not bias lineage comparisons.
4. **Correlation with Heterozygosity:** Test Pearson and Spearman rank correlation between $F_{ROH}$ and $(1 - H_o)$ as an internal negative control (a strong negative correlation is mathematically expected; deviation flags technical batch artifacts).

---

## 5. Planned Figures and Tables

- **Figure 1:** Geographic provenance map and study workflow showing cohort composition across 40 breeds and 10 countries.
- **Figure 2:** Population genomic structure (PCA of PC1 vs PC2 and PC2 vs PC3) separating Taurine, Zebu, and Sanga/Zenga lineages.
- **Figure 3:** Distribution of total genomic inbreeding ($F_{ROH}$) across breeds and evolutionary lineages.
- **Figure 4:** Temporal dissection of inbreeding: Proportion of $F_{ROH}$ contributed by each of the 5 segment length classes.
- **Figure 5:** Manhattan plot of genome-wide autozygosity islands (BTA 1–29) and candidate gene annotations.
- **Table 1:** Cohort summary table (breed names, sample sizes, countries, lineage classification, mean $H_o$, $H_e$, and mean $F_{ROH}$).
- **Table 2:** Lineage comparison of $F_{ROH}$ stratified by segment length classes.
- **Table 3:** Top shared autozygosity islands, genomic coordinates on `ARS-UCD1.2`, peak frequencies, and overlapping candidate genes.
- **Table S1:** Complete individual-level metrics for all 1,063 animals.
- **Table S2:** Full QC exclusion ledger.

---

## 6. Criteria for Null or Contradictory Results

- If Sanga/Zenga composites show *higher* autozygosity than purebred Taurine or Zebu, Hypothesis 2 is refuted. The manuscript will report this transparently, discussing whether recent closed-herd breeding or geographic isolation in the Great Lakes region counteracted ancestral hybrid vigor.
- If no autozygosity island exceeds the empirical top 1% threshold, candidate regions will be reported as exploratory, and no claims of positive selection will be made.
