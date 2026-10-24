# Data Dictionary & Metadata Specification

**Study:** Genomic Dissection of Autozygosity Landscapes, Recent vs. Historical Inbreeding, and Population Differentiation in Indigenous Tropical Cattle  
**Reference Dataset:** Genomic Reference Resource for African Cattle (GRRFAC)  
**Primary Archive:** Zenodo (Record 8317391, DOI: 10.5281/zenodo.8317391)  
**Published Descriptor:** *Nature Scientific Data* 11:792 (2024), DOI: 10.1038/s41597-024-03589-2  

---

## 1. Primary Metadata Table (`data/metadata/sample_metadata.tsv`)

| Field Name | Data Type | Description | Allowed Values / Examples |
| :--- | :--- | :--- | :--- |
| `S/NO` | Integer | Original serial numbering in GRRFAC Table S4 | 1 to 1082 |
| `VCF_ID` | String | Unique sample identifier matching genotype columns in the VCF header | e.g. `SU019201`, `ADGX1005311`, `Boran_01` |
| `Population` | String | Broad continental population grouping based on geographic and lineage ancestry | `East African Zebu`, `West African Taurine`, `East African Admixed`, `West African Zebu`, `Central African Admixed`, `Central African Zebu`, `West African Admixed` |
| `Breed` | String | Standardized breed classification of the individual animal | 40 distinct breeds (e.g. `Ndama`, `East Africa Shorthorn Zebu`, `Kuri`, `Toupouri`, `Ethiopia Boran`, `Kenana`, `Butana`, `Ankole`, `Sheko`) |
| `INDV3` | String | Alternative institutional individual accession code | e.g. `ILRI_00124` |
| `SOURCE` | String | Data origin (newly sequenced under GRRFAC or collated from public domain) | `CTLGH/GRRFAC`, `Public`, `ILRI` |
| `F_MISS` | Float | Variant missingness fraction reported in original provider QC | Range 0.000 to 0.233 |
| `Provider_QC` | String | QC status assigned by original dataset creators in Table S4 | `PASSED` ($N=1080$), `FILTERED` ($N=2$) |
| `Country_of_origin` | String | Country where the sample was collected | `Nigeria`, `Kenya`, `Chad`, `Ethiopia`, `Rwanda`, `Sudan`, `Uganda`, `Guinea`, `Senegal`, `Burkina Faso` |
| `Country_code` | String | ISO 3166-1 alpha-3 international country code | `NGA`, `KEN`, `TCD`, `ETH`, `RWA`, `SDN`, `UGA`, `GIN`, `SEN`, `BFA` |
| `Exclusion_Reason`| String | Specific reason for exclusion if flagged in Table S6 | `Duplicate / Relatedness (PHI >= 0.498)`, `High Missingness (>20%)`, `None (Included)` |
| `Study_Cohort` | String | Final inclusion status for downstream population genomics pipeline | `INCLUDED` ($N=1063$), `EXCLUDED` ($N=19$) |

---

## 2. Genomic Data Specification (`grrfac_public_sample1063_HDSNP_ARS-UCD1.2.vcf.gz`)

- **Format:** Variant Call Format (VCF) version 4.2, Block-GZIP (BGZF) compressed with corresponding Tabix index (`.tbi`).
- **Reference Assembly:** Hereford *Bos taurus* genome assembly `ARS-UCD1.2` (GCA_002263795.2).
- **Chromosomes:** Autosomes BTA 1 through BTA 29.
- **Marker Density:** ~777,962 high-density SNP loci from the Illumina BovineHD Genotyping BeadChip.
- **Genotype Fields:** `GT` (phased or unphased diploid calls: `0/0` homozygous reference, `0/1` heterozygous, `1/1` homozygous alternate, `./.` missing).

---

## 3. Evolutionary Lineage Classification

For comparative population genetics, the 40 breeds are grouped into three primary evolutionary lineages:

1. **African Taurine (*Bos taurus africanus*):**
   - Breeds: N'Dama (Guinea, Senegal), Baoulé, Somba, Lagune.
   - Phenotype: Humpless, compact stature, pronounced innate trypanotolerance.
2. **African Zebu (*Bos taurus indicus*):**
   - Breeds: East African Shorthorn Zebu (EASZ), Ethiopia Boran, Orma Boran, Kenana, Butana, Baggara, Red Bororo, Bunaji, Adamawa Gudali, Sokoto Gudali.
   - Phenotype: Prominent thoracic/cervico-thoracic hump, pendulous dewlap, heat tolerance, tick/vector resistance.
3. **Sanga & Zenga Composites (Ancient Admixed Lineages):**
   - Breeds: Ankole (Uganda, Rwanda), Inyambo (Rwanda), Inkuku, Inkungu, Sheko (Ethiopia - taurine/zenga intermediate), Kuri (Chad/Nigeria - lake-adapted gigantean horns), Toupouri, Nganda, Arado.
   - Phenotype: Stabilized historical hybrids (taurine $\times$ indicine) adapted to transitional agro-ecological zones.
