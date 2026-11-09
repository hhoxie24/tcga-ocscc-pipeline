# TCGA Oral Cavity SCC: Transcriptomic Analysis Pipeline

An R / R Markdown pipeline for exploring oral cavity squamous cell carcinoma
(OCSCC) patient data from The Cancer Genome Atlas head and neck cohort (TCGA-HNSC). 
It covers data acquisition from the GDC, sample and clinical metadata, mRNA and
miRNA differential expression, and weighted gene co-expression network analysis
(WGCNA) with module-trait correlation and pathway enrichment.

Parameters are set in a JSON configuration file, every module writes a
timestamped log, and one parameter that requires manual inspection
(WGCNA soft-thresholding power) is recorded explicitly rather than
hard-coded.

## Modules

| Module | Script | Purpose |
|---|---|---|
| 0 | `Module_0_DataAcquisition.Rmd` | Download clinical, mRNA, miRNA and somatic mutation data from the GDC; restrict to OCSCC patients |
| 1 | `Module_1_Metadata.Rmd` | Patient- and sample-level metadata: stage, HPV status, mutation flags and burden |
| 2 | `Module_2_DifferentialExpression.Rmd` | limma-voom differential expression (mRNA, miRNA) and miRNA-mRNA validated target pairs |
| 3A | `Module_3A_WGCNA_Diagnostic.Rmd` | WGCNA input matrix and soft-threshold diagnostic; stops for the power to be recorded |
| 3B | `Module_3B_WGCNA_Network.Rmd` | Signed co-expression network, module-trait correlations, eigengene trajectories, enrichment, hub genes |

## Cohort

- **Source:** TCGA-HNSC, via the NCI Genomic Data Commons (GDC) and TCGAbiolinks.
- **Selection:** patients whose `tissue_or_organ_of_origin` is an oral cavity site
  (list in Module 0). The list includes `Lip, NOS` and `Mandible`; lip is a
  separate ICD-O-3 topography group (C00) from oral cavity (C02-C06), and
  mandible is a bone site. Edit the list in Module 0 to apply a stricter definition.
- **Size (GDC Data Release 45.0):** 252 patients; 264 RNA-seq samples
  (247 primary tumor, 17 normal adjacent tissue, NAT) after excluding one
  metastatic (sample type `06`) sample. Omics layers contain different subsets of
  the 252 patients.
- **HPV status** is taken from the TCGA HNSC publication supplement and is
  available for 108 of 252 patients.
- GDC data releases change over time; a fresh download may differ from these counts.
  Module 0 logs the release in use.

## Repository structure

```
tcga-ocscc-pipeline/
├── config/
│   ├── config.json            # parameters (tracked; placeholder path)
│   ├── config_local.json      # your project_dir (not tracked; create this)
│   └── runs/<run_id>.json     # recorded WGCNA soft-thresholding power
├── scripts/
│   ├── helpers/
│   │   ├── config.R           # config loading, paths, logging
│   │   ├── de_helpers.R       # limma-voom and volcano-plot functions
│   │   └── run_record.R       # run record and input fingerprint checks
│   └── Module_*.Rmd
├── CHANGELOG.md
├── LICENSE
└── README.md
```

Data, results and logs are written to a separate project directory.

```
<project_dir>/
├── data/raw/          # GDC downloads
├── data/external/     # HPV table (user-supplied)
├── data/processed/    # .rds objects passed between modules
├── results/           # tables and figures (DE/, WGCNA/)
└── logs/              # 1 timestamped log per module run
```

## Requirements

Tested with R 4.4.2 on macOS. Packages (versions used):

| Source | Packages |
|---|---|
| Bioconductor | TCGAbiolinks 2.34.1, SummarizedExperiment 1.36.0, maftools 2.22.0, edgeR 4.4.2, limma 3.62.2, org.Hs.eg.db 3.20.0, clusterProfiler 4.15.1.002, ReactomePA 1.50.0, multiMiR 1.28.0, miRBaseConverter 1.30.0 |
| CRAN | WGCNA 1.74, dplyr 1.2.0, tidyr 1.3.2, ggplot2 4.0.2, ggrepel 0.9.8, ggpubr 0.6.3, cowplot 1.2.0, here 1.0.2, jsonlite 2.0.0, digest 0.6.39, rmarkdown 2.30, knitr 1.51 |

The clusterProfiler version used is a development build.

## Setup

1. Clone the repository.
2. Create `config/config_local.json` with an absolute path for outputs:
   ```json
   {"paths": {"project_dir": "/absolute/path/to/project"}}
   ```
3. **External data:** download the HPV status table from the supplementary
   information of the TCGA HNSC publication (Cancer Genome Atlas Network, 2015;
   see References), save it as tab-delimited text with columns `Barcode` and
   `Final_HPV_Status`, and place it at
   `<project_dir>/data/external/HPV_summary_TCGA2015NatureSupp1_2.txt`.
   The table is not redistributed here.

## Running

Knit each module in order (RStudio Knit button or `rmarkdown::render()`):

1. Module 0, Module 1, Module 2
2. **Module 3A**, then:
   - inspect `results/WGCNA/SoftThreshold.pdf` and `soft_threshold_candidates.csv`;
   - in `config/runs/<run_id>.json` (created by 3A), set `selected_value`,
     `selected_by` and `selected_on` (YYYY-MM-DD).
3. **Module 3B**

Logs are written to `<project_dir>/logs/`, each beginning with the module name,
start time, project directory and git commit.

The soft-thresholding power is chosen by reading the scale-free topology
diagnostic, so it cannot be computed or defaulted. Module
3B stops if no power has been recorded, or if the current input no longer matches
the fingerprint, for example after changing the gene-type filter. A changed
input always requires a new reading of the diagnostic.

## Configuration (`config/config.json`)

| Key | Default | Meaning |
|---|---|---|
| `run_id` | `ocscc_264` | names the run record file |
| `compute.seed`, `compute.wgcna_threads` | 123, 4 | reproducibility and threading |
| `mutations.flag_genes` | TP53, CDKN2A, FAT1, NOTCH1, PIK3CA, CASP8, NSD1, HRAS, TTN, MUC16 | genes given per-patient mutation flags |
| `mutations.nonsynonymous_classes` | maftools default classes | variant classes counted as mutations |
| `de.adj_p`, `de.abs_logfc` | 0.05, 1 | significance thresholds |
| `de.gene_types` | `["protein_coding"]` | GENCODE gene types analysed; `null` for all |
| `mirna.mirbase_version` | `v21` | miRBase version of GDC miRNA IDs |
| `mirna.target_table` | `mirtarbase` | multiMiR table |
| `mirna.support_types` | `["Functional MTI"]` | miRTarBase support categories kept; other categories are `Functional MTI (Weak)`, `Non-Functional MTI`, `Non-Functional MTI (Weak)`; `null` keeps all |
| `wgcna.top_n_by_variance` | 10000 | genes used for the network |
| `wgcna.network_type` | `signed` | network and soft-threshold fit type |
| `wgcna.min_module_size`, `wgcna.merge_cut_height` | 30, 0.25 | module detection |
| `wgcna.max_block_size` | `null` | `null` builds one block; set a number to split |
| `wgcna.hub_kme_threshold` | 0.8 | hub gene module membership cutoff |
| `wgcna.stage_traits_tumor_only` | `true` | stage traits compare tumors only (NAT = `NA`) |
| `wgcna.expression_traits` | `["RARG"]` | genes whose expression is used as continuous traits |
| `wgcna.soft_power.*` | powers 1-20, R² 0.85 | candidate powers and selection criterion |

## Analysis notes and limitations

- Stage contrasts are fairly weak at this sample size (Stage I n = 20).
- HPV status is known for 108 of 252 patients and 8 are HPV-positive.
- GDC miRNA counts are per precursor. Target pairs assign each precursor to the
  validated targets of every mature miRNA it produces, so pairs are not resolved
  by arm (5p/3p).
- Mutation, HPV and stage traits are tested among tumors only; tumor vs normal
  is captured by `Tumor_Status`.
- KEGG, Reactome, multiMiR results depend on database versions at run time.

## Planned

- Mature-level miRNA quantification (GDC isoform expression files) so that
  differential expression and target matching are resolved per mature miRNA.

## Acknowledgements

The results shown here are in whole or part based upon data generated by the
TCGA Research Network: https://www.cancer.gov/tcga.

Code refactoring, documentation, and review were assisted by Claude Opus 5.5
(Anthropic). All analytical decisions, parameter choices, and interpretation of results were made and verified by the author.

## References

**Data**

- Cancer Genome Atlas Network (2015). Comprehensive genomic characterization of
  head and neck squamous cell carcinomas. *Nature* 517(7536):576-582.
  doi:10.1038/nature14129

**Methods and software**

- Colaprico A, Silva TC, Olsen C, et al. (2015). TCGAbiolinks: an R/Bioconductor
  package for integrative analysis of TCGA data. *Nucleic Acids Research*.
  doi:10.1093/nar/gkv1507
- Mayakonda A, Lin D, Assenov Y, Plass C, Koeffler HP (2018). Maftools: efficient
  and comprehensive analysis of somatic variants in cancer. *Genome Research*.
  doi:10.1101/gr.239244.118
- Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). edgeR v4: powerful
  differential analysis of sequencing data with expanded functionality and
  improved support for small counts and larger datasets. *Nucleic Acids Research*
  53(2):gkaf018. doi:10.1093/nar/gkaf018
- Ritchie ME, Phipson B, Wu D, et al. (2015). limma powers differential expression
  analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*
  43(7):e47. doi:10.1093/nar/gkv007
- Law CW, Chen Y, Shi W, Smyth GK (2014). voom: precision weights unlock linear
  model analysis tools for RNA-seq read counts. *Genome Biology* 15(2):R29.
  doi:10.1186/gb-2014-15-2-r29
- Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation
  network analysis. *BMC Bioinformatics* 9:559. doi:10.1186/1471-2105-9-559
- Langfelder P, Horvath S (2012). Fast R functions for robust correlations and
  hierarchical clustering. *Journal of Statistical Software* 46(11):1-17.
- Wu T, Hu E, Xu S, et al. (2021). clusterProfiler 4.0: a universal enrichment
  tool for interpreting omics data. *The Innovation* 2(3):100141.
  doi:10.1016/j.xinn.2021.100141
- Xu S, Hu E, Cai Y, et al. (2024). Using clusterProfiler to characterize
  multiomics data. *Nature Protocols* 19(11):3292-3320.
  doi:10.1038/s41596-024-01020-z
- Yu G, He Q (2016). ReactomePA: an R/Bioconductor package for reactome pathway
  analysis and visualization. *Molecular BioSystems* 12(12):477-479.
  doi:10.1039/C5MB00663E
- Ru Y, Kechris KJ, Tabakoff B, et al. (2014). The multiMiR R package and database:
  integration of microRNA-target interactions along with their disease and drug
  associations. *Nucleic Acids Research* 42(17):e133. doi:10.1093/nar/gku631
- Huang HY, Lin YCD, Cui S, et al. (2022). miRTarBase update 2022: an informative
  resource for experimentally validated miRNA-target interactions.
  *Nucleic Acids Research* 50(D1):D222-D230. doi:10.1093/nar/gkab1079
- Xu T, Su N, Liu L, et al. (2018). miRBaseConverter: a comprehensive and
  high-efficiency tool for converting and retrieving the information of miRNAs
  in different miRBase versions. *BMC Bioinformatics* 19(Suppl 19):514.
  doi:10.1186/s12859-018-2531-5
- Carlson M (2024). org.Hs.eg.db: Genome wide annotation for Human. R package
  version 3.20.0.
- Lei. How to Build Gene Co-expression Networks from RNA-seq Data Using WGCNA –
  Complete Step-by-Step Guide For Absolute Beginners. NGS101.
  https://ngs101.com/how-to-build-gene-co-expression-networks-from-rna-seq-data-using-wgcna-complete-step-by-step-guide-for-absolute-beginners/

## License

MIT; see `LICENSE`. The license covers the code in this repository, not TCGA
data or the external HPV table.
