# Dietary reconstruction in sympatric muntjacs — analysis code and data

Analysis code accompanying:

> **Robustness of dietary reconstruction in sympatric muntjacs using fecal shotgun
> metagenomics and a curated plant reference database**
> Xia Qiu, Chao Ji, Yuanyuan Luo
> College of Life Sciences, China Jiliang University, Hangzhou, China

---

## Overview

This repository contains the complete bioinformatic workflow used to reconstruct
the diets of two sympatric muntjac species from fecal shotgun metagenomic data:

- **Reeves' muntjac** (*Muntiacus reevesi*), n = 52 samples
- **black muntjac** (*M. crinifrons*), n = 13 samples

collected in Jiulongshan National Nature Reserve, Zhejiang Province, China
(December 2024). Dietary taxa were assigned at the **plant family** level against a
**curated multi-barcode reference database** built from local vegetation surveys and
published dietary records.

Two independent reconstruction strategies were compared:

| Strategy | Principle | Thresholds |
|---|---|---|
| **DIR** (direct read-based) | map reads directly to reference barcodes with Bowtie2 | Loose / Medium / Strict |
| **ASM** (assembly-based) | assemble contigs with MEGAHIT, then BLASTn | pident 90% / 95% / 99% |

**ASM95 was selected as the primary dataset** for downstream ecological analyses
(see `04_reproduction_scripts/asm_threshold_mrra.py` for the quantitative justification).

---

## Repository structure

```
.
├── 01_upstream_pipeline/           # raw reads -> family-level matrices
│   ├── 01_DIR_read_mapping.sbatch
│   └── 02_ASM_assembly_blast_tpm.sbatch
│
├── 02_downstream_analysis/         # R/Python: matrices -> figures
│   ├── 03_figures_main.R
│   ├── 04_figure1_map.R
│   ├── 05_verify_ASM_reproduce.R
│   └── 06_build_TPM_matrix.py
│
├── 03_host_assignment/             # mitochondrial diagnostic-site voting
│   ├── mito_hostID.sbatch
│   ├── check_env.sh
│   ├── refs/muntjac_mito_combined.fna
│   ├── sites_all.tsv
│   ├── sites_transversion_only.tsv
│   ├── meta.json
│   └── README_mito_method.md
│
├── 04_reproduction_scripts/        # scripts that generate every number in the paper
│   ├── build_diet_tables.py
│   ├── build_downstream_data.py
│   ├── build_downstream_new.py
│   ├── compute_downstream_stats.R
│   ├── repro_all_values.py
│   ├── topn_scientific.py
│   └── asm_threshold_mrra.py
│
├── environment.yml
├── requirements.txt
├── sessionInfo.txt
└── LICENSE
```

---

## Sample naming

Samples are labelled with a **two-part prefix**:

| Prefix | Meaning |
|---|---|
| **D\*** | samples from the core 30-ha forest dynamics plot area |
| **E\*** | samples from the satellite plots |

The number following the letter is the field sample ID. Sample IDs are identical to
those in **Table S2** and in the archived sequencing data.

**65 samples were retained** for analysis after quality control (52 *M. reevesi* +
13 *M. crinifrons*); 46 of these yielded plant families above the ASM95 detection
threshold (40 *M. reevesi* + 6 *M. crinifrons*) and were used for the ecological
analyses.

---

## Two workflows: host assignment vs. dietary reconstruction

These are **independent** and use different reference sequences.

### A. Host assignment (`03_host_assignment/`)

Assigns each fecal sample to a host species. Whole-genome mapping rates against the
two nuclear reference assemblies are **not** reliable here: the two assemblies differ
substantially in chromosome number (2n = 46 vs. 2n = 4) and scaffold N50
(113 Mb vs. 647 Mb), so mapping rate is governed by reference quality rather than by
the true source species.

Instead, host identity is assigned by the **bases carried by individual reads at
interspecific diagnostic sites** between the two mitochondrial genomes:

| | |
|---|---|
| *M. reevesi* reference | `NC_004069.1` (16,354 bp) |
| *M. crinifrons* reference | `NC_004577.1` (16,357 bp) |
| Diagnostic sites | **947 substitutions** (864 transitions, 83 transversions) |
| Assignment rule | base matching a reference = vote; exact two-sided binomial test, p < 0.05 |

A 150 bp read spans ~8.7 diagnostic sites, giving ~58-fold separation between correct
and incorrect assignment. Assignments were verified using the 83 transversion sites
alone as an independent, more conservative criterion.

### B. Dietary reconstruction (`01_upstream_pipeline/`)

Uses the **curated plant barcode reference database** (not the mitogenomes) to assign
host-depleted reads to plant families. See `01_upstream_pipeline/`.

---

## Running the pipeline

### Requirements

```bash
conda env create -f environment.yml
conda activate muntjac_diet
```

### Step 1 — host assignment

```bash
cd 03_host_assignment
bash check_env.sh                 # verify tools are available
sbatch mito_hostID.sbatch         # parallelisable: PARALLEL=8 sbatch mito_hostID.sbatch
```

Outputs `host_call_allsites.tsv` and `host_call_transversion.tsv`.

### Step 2 — dietary reconstruction

```bash
cd 01_upstream_pipeline
sbatch 01_DIR_read_mapping.sbatch
sbatch 02_ASM_assembly_blast_tpm.sbatch
```

Both scripts expect host-depleted reads in a single directory; set `CLEAN_DIR` at the
top of each script.

### Step 3 — downstream statistics

```bash
cd 04_reproduction_scripts
python build_diet_tables.py            # family-level diet tables (mRRA, FOO, wPOO)
python build_downstream_data.py        # RRA matrices + rarefied counts
Rscript compute_downstream_stats.R     # beta diversity, IndVal, niche metrics
```

---

## Reproducing specific results

| Result | Script |
|---|---|
| **§3.1** QC and host-removal statistics | `01_upstream_pipeline/01_DIR_read_mapping.sbatch` (flagstat) |
| **§3.2** cross-method Jaccard similarity | `04_reproduction_scripts/repro_all_values.py` |
| **§3.2** top-N ranking method | `04_reproduction_scripts/topn_scientific.py` |
| **§3.3** family-level diet composition (mRRA, FOO) | `04_reproduction_scripts/build_diet_tables.py` |
| **§3.4** alpha diversity | `04_reproduction_scripts/compute_downstream_stats.R` |
| **§3.4** beta diversity (PERMANOVA, PERMDISP, SIMPER) | `04_reproduction_scripts/compute_downstream_stats.R` |
| **§3.4** indicator taxa (IndVal) | `04_reproduction_scripts/compute_downstream_stats.R` |
| **§3.4** niche breadth, overlap, specialisation | `04_reproduction_scripts/compute_downstream_stats.R` |
| **ASM90/95/99 threshold choice** | `04_reproduction_scripts/asm_threshold_mrra.py` |
| **Figures** | `02_downstream_analysis/03_figures_main.R` |

---

## Data

| Data | Repository | Accession |
|---|---|---|
| Raw sequencing reads (65 samples, paired-end 150 bp) | GSA (NGDC) | *PRJCA to be inserted* |
| Curated plant reference database, family-level matrices, diet tables, statistical outputs | Zenodo | *DOI to be inserted* |

The curated plant reference database (`seq_Table.fasta`, 1,589 sequences) and its
sequence-to-species lookup table (`Reference Table.xlsx`) are archived on Zenodo
together with all derived matrices.

---

## Three important notes for reproducibility

### 1. The taxonomic dictionary requires three manual mappings

The species-to-family dictionary (`official_species_family_dictionary.csv`) does **not**
contain three genera that occur in the reference database. These must be mapped
explicitly, otherwise the family counts will be **28 / 31 instead of 29 / 32**:

```python
FAMILY_FIX = {
    "itea":        "Iteaceae",
    "indocalamus": "Poaceae",
    "ilex":        "Aquifoliaceae",
}
# plus: reference sequence "seq_695" -> Schisandraceae
```

### 2. Alpha diversity requires the rarefied count matrices

Alpha diversity (richness, Shannon, Pielou) is **not** computed from the raw abundance
matrices. It is computed from `Counts_rarefied_ASM{90,95}.csv`, in which each sample's
reads have been rarefied to the minimum total sequencing depth
(R = 3,424,157 reads) using `numpy.random.default_rng(1)`.

**If the rarefied matrices (or the fixed random seed) are not used, the alpha-diversity
values in §3.4 cannot be reproduced.**

### 3. `sites_hit` counts site × chromosome, not independent sites

The mitochondrial voting output reports `sites_hit = 1894`. This is **947 diagnostic
sites recorded in two coordinate frames** (one per reference mitogenome), not 1,894
independent sites. The two references differ in length (16,354 vs. 16,357 bp), so each
site is listed once per coordinate system.

---

## Reference sequence species names

Species names follow the current NCBI taxonomy. Note in particular that
`MH500770.1` is ***Muntiacus vaginalis*** (northern red muntjac), not *M. muntjak*.

---

## License

Code: **MIT** (see `LICENSE`)
Data: **CC-BY-4.0**

---

## Contact

Xia Qiu — College of Life Sciences, China Jiliang University
