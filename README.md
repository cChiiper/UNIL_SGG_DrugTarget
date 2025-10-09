Hey there, early bird!

You’ve landed here a bit early — some of the scripts are still being polished and others will be added soon.
If you’d like more details before the publication is out, feel free to contact me.
This page will be updated over the coming weeks.

🗓️ Last updated: October 9, 2025

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
- 01_Get_IVs: Contains scripts to select instrumental variables
- 02_protein_to_trait: Contains scripts to run the MR analysis

#### 02_combine
Contains the script to test the various integration method. 

#### 03_benchmark
This section is currently being updated and will be expanded over the next few days.

#### 04_rg_ldsc
Contains the scripts to compute genetic correlation from [LDSC](https://github.com/bulik/ldsc)
