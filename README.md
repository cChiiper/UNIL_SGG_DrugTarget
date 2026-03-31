# UNIL_SGG_DrugTarget
Code accompanying the manuscript "Integration of genetic evidence to identify approved drug targets". 

## Authors
The scripts were written by Samuel Moix, with partial assistance from large language models (GitHub Copilot and OpenAI GPT models) for code suggestions and drafting, and include some scripts originally written by [Marie Sadler](https://github.com/masadler) that were reused or modified.

## Folders

### data
This folder contains files required to run the scripts. Other required files are publicly available online.

| File                      | Description                                                                       |
|---------------------------|-----------------------------------------------------------------------------------|
| `gene_universe_rbs.txt`   | Genes used for the cross-trait analysis.                                          |
| `HLA_uk10k_rsids.txt`     | UK10K rsIDs of SNPs in the extended HLA region.                                   |
| `no_pQTL_proteins.txt`    | IDs of proteins without instruments for Mendelian randomization.                  |
| `target_genes_TTD.csv`    | List of drug target genes for the 30 traits (`1` = target, `0` = not a target).   |


### scripts
The scripts folder contains the necessary code in R, Python, and Bash, with some workflows managed through Snakemake.

#### 01_MR
Contains scripts to perform Mendelian randomization (MR) to estimate the effect of protein expression on the 30 traits.
- 01_Get_IVs: Script to select instrumental variables
- 02_protein_to_trait: Script to run the MR analysis

#### 02_combine
Contains the script to test the various integration methods. 
- 01_aggregation.R: Script to combine the different prioritization methods (also includes OR calculation)
- 02_cross_trait.R: Script to perform cross trait analysis (also includes benchmark)

#### 03_benchmark
Contains standalone scripts to summarize key benchmarking metrics.
- 01_delta_diff.R: Computes paired mean-rank differences against a baseline method
- 02_roc_metrics.R: Computes AUROC values, ROC curve coordinates, and pdiff summaries
- 03_jaccard_index.R: Computes pairwise Jaccard indices between trait drug target sets

#### 04_rg_ldsc
Contains the scripts to compute genetic correlation from [LDSC](https://github.com/bulik/ldsc)
