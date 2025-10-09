###################################################################################
### Script to aggregate p-values from different sources and compute odds ratios ###
### Author: S. Moix, DBC, UNIL, 2025 with partial assistance from LLMs ############
### Date: October 2025 ############################################################
###################################################################################

################################################
### Load libraries #############################

library(dplyr)
library(progress)
library(conflicted)
library(argparse)

################################################
### Parser ####################################

# Create a parser
parser <- ArgumentParser(description = "Script to run SLURM job with varying TRAIT and DTRAIT")

# Add arguments
parser$add_argument("--trait", type = "character", required = TRUE, help = "The trait to analyze.")
parser$add_argument("--dtrait", type = "character", required = TRUE, help = "The drug trait to analyze.")

# Parse the arguments
args <- parser$parse_args()

################################################
### Parameters #################################
TRAIT <- args$trait
DTRAIT <- args$dtrait

print(paste("Trait is: ", TRAIT))
print(paste("Drug trait is: ", DTRAIT))

par_ecdf <- TRUE
### Computed range (0 - 100)
par_max_rank_th <- 5.1
### Threshold for drug targets
nb_db_for_DT <- 3
par_TTD <- TRUE
### Subset to GWAS present genes
par_GWAS_uni <- TRUE
par_exome_uni <- FALSE
par_PPP_uni <- FALSE
### Don't keep NAs
par_full <- FALSE

### Method to merge ranks
par_method <- "minimum" # "minimum", "product", "mean", "PCA", weightsum

### Selected data
# Select from : p_rank_pQTL_PPP, p_rank_exome, p_rank_eQTL_GWAS_blood, p_rank_GWAS, p_rank_pQTL_GWAS, p_rank_LM, p_rank_OT
INCLUDED_DATA <- c("p_rank_GWAS", "p_rank_exome", "p_rank_eQTL_GWAS_blood")

################################################
### Set working directories ####################

### Path to files from Sadler et al 2023 DOI: 10.1016/j.xgen.2023.100341 
path_to_sadler <- "/.../"
### Path to MR results with UK Biobank Pharma Proteomics Project data (See Supplementary Table 16 to directly get p-values or run pipeline in the  01_MR folder)
path_to_MR_data <- "/.../"
### Path to drug targets from Sadler et al 2023 DOI: 10.1016/j.xgen.2023.100341 converted to hgnc (if empty "" it will be mapped here)
path_to_hgnc_DT <- "/.../"
### Path to target genes from The Therapeutic Target Database (TTD)
path_to_TTD <- "/.../target_genes_TTD.csv" # See data folder
### Path to output folder
path_to_output <- "/.../"

################################################
### Prepare drug targets dataframe #############

if(path_to_hgnc_DT == ""){
  library(biomaRt)
  drug_targets <- read.table(paste0("/.../",DTRAIT,"_drug_all_db.tsv"), header = T, sep = "\t") # If path_to_hgnc_DT is empty, path should be provided here
  drug_targets$Sum <- rowSums(drug_targets[, c("DrugBank_DGIdb", "Ruiz_DGIdb", "Ruiz_STITCH", "DrugBank_STITCH", "ChEMBL_ChEMBL")])
  
  ### convert DT EnsemblId to Gene name 
  
  # Set up the Ensembl dataset using biomaRt
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Query to convert ENSEMBL IDs to gene names
  gene_info <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    filters = 'ensembl_gene_id',
    values = drug_targets$EnsemblId,
    mart = ensembl
  )
  drug_targets <- merge(drug_targets, gene_info, by.x = "EnsemblId", by.y = "ensembl_gene_id", all.x = TRUE)
}else{
  drug_targets <-read.table(paste0(path_to_hgnc_DT, DTRAIT,"_drug_all_db_with_hgnc.tsv"), 
                   sep = "\t", header = T)
}

if(par_TTD){
  ### Add genes that are present in TTD
  TTD_genes <- read.csv(path_to_TTD,
                        header = TRUE, sep = ",")
  
  drug_targets <- drug_targets %>% 
    select(-Sum)
  trait_genes <- TTD_genes %>%
    dplyr::filter(!!sym(DTRAIT) == 1) %>%
    dplyr::pull(Gene)
  
  ### We won't add genes that were not present
  drug_targets <- drug_targets %>%
    dplyr::mutate(TTD = ifelse(hgnc_symbol %in% trait_genes, 1, 0))
  
  # Recalculate the Sum column: sum of database-related columns
  drug_targets <- drug_targets %>%
    dplyr::mutate(Sum = DrugBank_DGIdb + Ruiz_DGIdb + Ruiz_STITCH + DrugBank_STITCH + ChEMBL_ChEMBL + TTD)
}

### As there can be multiple hgnc (different EnsemblId) for the same gene, ###
### we will keep the one with the highest sum ################################
drug_targets <- drug_targets %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::slice_max(order_by = Sum, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

drug_targets_genes <- drug_targets %>%
  dplyr::rename(Gene = hgnc_symbol) %>%
  dplyr::filter(Sum >= nb_db_for_DT) %>% 
  dplyr::filter(!is.na(Gene)) %>%
  dplyr::filter(Gene != "")

################################################
### Load data ##################################

### pQTL --> GWAS 
MR_forward <- read.table(file.path(path_to_MR_data, "MS_merged_results.tsv"), 
                         sep = "\t", header = T)

### Sadler 2023 data
exome_data <- read.table(file.path(path_to_sadler, "/Zenodo/gene_prioritization/main/Exome/", 
                                   paste0(TRAIT ,"_method_Exome.tsv")), sep = "\t", header = T)

eQTL_GWAS_blood_data <- read.table(file.path(path_to_sadler, "/Zenodo/gene_prioritization/main/eQTL_GWAS_blood/", 
                                             paste0(TRAIT ,"_method_eQTL_GWAS_blood.tsv")), sep = "\t", header = T)   

GWAS_data <- read.table(file.path(path_to_sadler, "/Zenodo/gene_prioritization/main/GWAS/", 
                                  paste0(TRAIT ,"_method_GWAS.tsv")), sep = "\t", header = T)

pQTL_GWAS_data <- read.table(file.path(path_to_sadler, "/Zenodo/gene_prioritization/main/pQTL_GWAS/", 
                                       paste0(TRAIT ,"_method_pQTL_GWAS.tsv")), sep = "\t", header = T)


################################################
### Subset data ################################

MR_forward <- MR_forward %>%
  dplyr::filter(outcome == TRAIT) %>%
  dplyr::filter(method == "Inverse variance weighted" | method == "Wald ratio") %>%
  dplyr::mutate(Gene = sub(":.*", "", exposure)) %>%
  dplyr::arrange(((b)))

### Drop duplicated genes !!! WARNING !!! Duplicated genes (i.e., same hgnc symbol) are simply dropped here 
# (i.e., keep only the first occurrence)
MR_forward <- MR_forward %>% 
  dplyr::group_by(Gene) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup()

################################################
### Add ranks and keep only ranks ##############

if(par_method != "PCA"){
  ### Trait Gene 
  MR_forward <- MR_forward %>%
    # dplyr::filter(pval < 0.05/nrow(MR_forward)) %>%
    # dplyr::mutate(pval = ifelse(pval < 0.05 / nrow(MR_forward), pval, NA)) %>%
    dplyr::mutate(p_rank_pQTL_PPP = min_rank(pval)) 
  MR_forward <- MR_forward %>% 
    dplyr::mutate(p_rank_pQTL_PPP = p_rank_pQTL_PPP/nrow(MR_forward)*100) %>%
    dplyr::select(c("Gene","p_rank_pQTL_PPP"))

  exome_data <- exome_data %>%
    # dplyr::filter(p_value < 0.05/nrow(exome_data)) %>%
    # dplyr::mutate(p_value = ifelse(p_value < 0.05 / nrow(exome_data), p_value, NA)) %>%
    dplyr::mutate(p_rank_exome = min_rank(p_value)) 
  exome_data <- exome_data %>%
    dplyr::mutate(p_rank_exome = p_rank_exome/nrow(exome_data)*100) %>%
    dplyr::select(c("Gene", "p_rank_exome"))

  eQTL_GWAS_blood_data <- eQTL_GWAS_blood_data %>%
    # dplyr::filter(p_value < 0.05/nrow(eQTL_GWAS_blood_data)) %>%
    # dplyr::mutate(p_value = ifelse(p_value < 0.05 / nrow(eQTL_GWAS_blood_data), p_value, NA)) %>%
    dplyr::mutate(p_rank_eQTL_GWAS_blood = min_rank(p_value))
  eQTL_GWAS_blood_data <- eQTL_GWAS_blood_data %>%
    dplyr::mutate(p_rank_eQTL_GWAS_blood = p_rank_eQTL_GWAS_blood/nrow(eQTL_GWAS_blood_data)*100) %>%
    dplyr::select(c("Gene", "p_rank_eQTL_GWAS_blood"))

  GWAS_data <- GWAS_data %>%
    # dplyr::filter(p_value < 0.05/nrow(GWAS_data)) %>%
    # dplyr::mutate(p_value = ifelse(p_value < 0.05 / nrow(GWAS_data), p_value, NA)) %>%
    dplyr::mutate(p_rank_GWAS = min_rank(p_value))
  GWAS_data <- GWAS_data %>%
    dplyr::mutate(p_rank_GWAS = p_rank_GWAS/nrow(GWAS_data)*100) %>%
    dplyr::select(c("Gene", "p_rank_GWAS"))

  pQTL_GWAS_data <- pQTL_GWAS_data %>%
    # dplyr::filter(p_value < 0.05/nrow(pQTL_GWAS_data)) %>%
    # dplyr::mutate(p_value = ifelse(p_value < 0.05 / nrow(pQTL_GWAS_data), p_value, NA)) %>%
    dplyr::mutate(p_rank_pQTL_GWAS = min_rank(p_value))
  pQTL_GWAS_data <- pQTL_GWAS_data %>%
    dplyr::mutate(p_rank_pQTL_GWAS = p_rank_pQTL_GWAS/nrow(pQTL_GWAS_data)*100) %>%
    dplyr::select(c("Gene", "p_rank_pQTL_GWAS"))
}

if(par_method == "PCA"){
  ### Trait Gene 
  MR_forward <- MR_forward %>% 
    dplyr::rename(p_rank_pQTL_PPP = pval) %>%
    dplyr::select(c("Gene","p_rank_pQTL_PPP"))

  exome_data <- exome_data %>%
      dplyr::rename(p_rank_exome = p_value) %>%
      dplyr::select(c("Gene", "p_rank_exome"))

  eQTL_GWAS_blood_data <- eQTL_GWAS_blood_data %>%
      dplyr::rename(p_rank_eQTL_GWAS_blood = p_value) %>%
      dplyr::select(c("Gene", "p_rank_eQTL_GWAS_blood"))

  GWAS_data <- GWAS_data %>%
    dplyr::rename(p_rank_GWAS = p_value) %>%
    dplyr::select(c("Gene", "p_rank_GWAS"))

  pQTL_GWAS_data <- pQTL_GWAS_data %>%
    dplyr::rename(p_rank_pQTL_GWAS = p_value) %>%
    dplyr::select(c("Gene", "p_rank_pQTL_GWAS"))

  if("p_rank_LM" %in% INCLUDED_DATA){
    pLM_data <- lm_protein_trait %>%
      dplyr::rename(p_rank_LM = p_value) %>%
      dplyr::select(c("Gene", "p_rank_LM"))
  }

  if("p_rank_OT" %in% INCLUDED_DATA){
    OTscore_data <- OT_TPF %>%
      dplyr::mutate(globalScore = 1 - globalScore) %>%
      dplyr::rename(p_rank_OT = globalScore) %>%
      dplyr::select(c("Gene", "p_rank_OT"))
  }
}

################################################
### Merge data #################################

merged_data <- merge(MR_forward, exome_data, by = "Gene", all = TRUE)
merged_data <- merge(merged_data, eQTL_GWAS_blood_data, by = "Gene", all = TRUE)
merged_data <- merge(merged_data, GWAS_data, by = "Gene", all = TRUE)
merged_data <- merge(merged_data, pQTL_GWAS_data, by = "Gene", all = TRUE)

### Define gene space 
GWAS_gene_space <- unique(GWAS_data$Gene)
Exome_gene_space <- unique(exome_data$Gene)
eQTL_GWAS_blood_gene_space <- unique(eQTL_GWAS_blood_data$Gene)
# pQTL_GWAS_gene_space <- unique(pQTL_GWAS_data$Gene)
pQTL_GWAS_gene_space <- unique(MR_forward$Gene)

if(par_full){
  gene_space_full <- Reduce(intersect, list(GWAS_gene_space, Exome_gene_space, eQTL_GWAS_blood_gene_space, pQTL_GWAS_gene_space))
  print(paste("Full gene space: ", length(gene_space_full)))
}

merged_data <- merged_data %>% 
  dplyr::select(c("Gene", all_of(INCLUDED_DATA))) %>% 
  dplyr::filter(if_any(all_of(INCLUDED_DATA), ~ !is.na(.))) # Remove rows not included

### Drop duplicated genes !!! WARNING !!! Duplicated genes (i.e., same hgnc symbol) are simply dropped here
# (i.e., keep only the first occurrence)
merged_data <- merged_data %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

if(par_full){
  merged_data <- merged_data[complete.cases(merged_data[, INCLUDED_DATA]), ]
  merged_data <- merged_data %>%
    dplyr::filter(Gene %in% gene_space_full)
}

################################################
### Change minimum GWAS ########################

if("p_rank_GWAS" %in% INCLUDED_DATA & length(INCLUDED_DATA) > 1){
  # 1) Identify the global minimum of p_rank_GWAS
  min_val <- min(merged_data$p_rank_GWAS, na.rm = TRUE)
  
  # 2) Compute the minimum of other p_rank columns for each row
  calculate_minimum <- function(row) {
    # Select columns that are in INCLUDED_DATA but not p_rank_GWAS
    other_values <- as.numeric(row[setdiff(INCLUDED_DATA, "p_rank_GWAS")])
    # Remove NA values
    other_values <- other_values[!is.na(other_values)]
    if (length(other_values) == 0) {
      return(NA)
    } else {
      return(min(other_values))
    }
  }
  merged_data$min_other_p_rank <- apply(merged_data, 1, calculate_minimum)
  
  # 3) Find rows that have the minimum p_rank_GWAS (and a non-NA computed minimum)
  idx_min <- which(!is.na(merged_data$min_other_p_rank) & merged_data$p_rank_GWAS == min_val)
  
  # 4) Spread those values evenly from 0 up to min_val, in order of their computed minimum
  if (length(idx_min) > 0) {
    # Order the indices by min_other_p_rank (lowest first)
    idx_min_ordered <- idx_min[order(merged_data$min_other_p_rank[idx_min])]
    
    # Create a sequence of equally spaced values from 0 to min_val
    n <- length(idx_min_ordered)
    seq_vals <- seq(0, min_val, length.out = n)
    
    # Assign the new values to p_rank_GWAS based on the ordering
    merged_data$p_rank_GWAS[idx_min_ordered] <- seq_vals
  }
}

################################################
### Get min rank and source ####################

if(par_method != "PCA"){
  merged_data <- merged_data %>%
  # Process the data row-by-row to perform calculations on each row individually.
  rowwise() %>%
  mutate(
    # Calculate 'min_rank' based on the specified 'par_method' parameter.
    min_rank = {
      if (par_method == "minimum") {
        # For the "minimum" method:
        # - Use do.call with pmin to compute the minimum value across all columns listed in INCLUDED_DATA.
        # - na.rm = TRUE ensures that NA values are ignored in the comparison.
        do.call(pmin, c(across(all_of(INCLUDED_DATA)), na.rm = TRUE))
      } else if (par_method == "product") {
        # For the "product" method:
        # - c_across collects values from the columns specified in INCLUDED_DATA into a vector.
        # - prod calculates the product of these values, ignoring NA values.
        # - log takes the natural logarithm of the product.
        # - The result is then divided by the count of non-NA values, which is calculated using sum(!is.na(...)).
        log(prod(c_across(all_of(INCLUDED_DATA)), na.rm = TRUE)) /
          sum(!is.na(c_across(all_of(INCLUDED_DATA))))
      } else if (par_method == "mean") {
        # For the "mean" method:
        # - Simply compute the mean of the values in the INCLUDED_DATA columns.
        # - na.rm = TRUE ensures that NA values are excluded from the calculation.
        mean(c_across(all_of(INCLUDED_DATA)), na.rm = TRUE)
      } else if (par_method == "weightsum") {
        # For the "weightsum" method:
        # 1) Drop NAs, invert to high-is-good (100 – x), sort descending.
        # 2) Weight by 1 / rank^2 and compute weighted average.
        # 3) Back-invert to low-is-good percentile.
        vals <- c_across(all_of(INCLUDED_DATA))
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
          NA_real_
        } else {
          inv_sorted <- sort(100 - vals, decreasing = TRUE)
          weights    <- 1 / (seq_along(inv_sorted)^2)
          S          <- sum(inv_sorted * weights) / sum(weights)   # weighted average
          100 - S                                               # back to original scale
        }
      } else {
        # If 'par_method' does not match any expected method, assign NA to 'min_rank'.
        NA_real_
      }
    },
    # Identify the column name corresponding to the minimum value across the INCLUDED_DATA columns.
    source_min = {
      # Collect all values from the columns in INCLUDED_DATA for the current row.
      vals <- c_across(all_of(INCLUDED_DATA))
      # Check if all values are NA; if so, assign NA_character_ to source_min.
      if (all(is.na(vals))) {
        NA_character_
      } else {
        # Determine the minimum value in the vector, ignoring any NA values.
        the_min <- min(vals, na.rm = TRUE)
        # Identify the index of the first occurrence of this minimum value.
        first_match <- which(vals == the_min)[1]
        # Use the index to return the corresponding column name from INCLUDED_DATA.
        INCLUDED_DATA[first_match]
      }
    }
  ) %>%
  # Remove the rowwise grouping to return to the normal data frame structure.
  ungroup()
}

if(par_method == "PCA"){
  # Convert p-values to z-statistics using two-tailed conversion: z = qnorm(1 - p).
  merged_data <- merged_data %>%
    dplyr::mutate(across(all_of(INCLUDED_DATA), ~ qnorm(1 - .)))

  # Replace machine precisions resulting INF | Not perfect but works more or less (could also take max)
  merged_data <- merged_data %>%
    dplyr::mutate(across(all_of(INCLUDED_DATA),
         ~ ifelse(. == Inf, qnorm(1 - .Machine$double.eps),
           ifelse(. == -Inf, qnorm(.Machine$double.eps), .))))

  # Replace NA values with the median of the respective column
  merged_data <- merged_data %>%
    mutate(across(all_of(INCLUDED_DATA), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

  # Run PCA on the specified columns
  pca_result <- prcomp(merged_data %>% dplyr::select(all_of(INCLUDED_DATA)), 
                      center = TRUE, scale. = TRUE)

  # Add the first principal component scores as a new column "min_rank"
  merged_data$min_rank <- pca_result$x[, 1]

  # Compute the row-wise mean of the original INCLUDED_DATA columns as a reference
  row_avg <- rowMeans(merged_data %>% dplyr::select(all_of(INCLUDED_DATA)), na.rm = TRUE)

  # Check if the correlation between min_rank and row_avg is positive;
  # if so, flip the sign of min_rank
  if(cor(merged_data$min_rank, row_avg, use = "complete.obs") > 0) {
    print("Changed sign of min_rank")
    merged_data$min_rank <- -merged_data$min_rank
  }

  # Place holder for the source of the minimum rank
  merged_data$source_min <- "PCA"
}

################################################
### ECDF transformation: #######################
if(par_ecdf){
  merged_data$adjusted_min_rank <- merged_data$min_rank
  ecdf_function <- ecdf(merged_data$adjusted_min_rank)
  merged_data <- merged_data %>%
    mutate(min_rank = 100 * (ecdf_function(adjusted_min_rank)))
}

################################################
### Prediction with Trait-Gene score ###########
prediction_data <- merged_data[,c("Gene","min_rank", "source_min")]

if(par_GWAS_uni){ # Same Gene universe
  prediction_data <- prediction_data[which(prediction_data$Gene %in% GWAS_gene_space),]
}
if(par_exome_uni){ # Same Gene universe
  prediction_data <- prediction_data[which(prediction_data$Gene %in% Exome_gene_space),]
}

if(par_PPP_uni){ # Same Gene universe
  prediction_data <- prediction_data[which(prediction_data$Gene %in% pQTL_GWAS_gene_space),]
}

prediction_data$Target <- ifelse(prediction_data$Gene %in% unique(drug_targets_genes$Gene), 1, 0)

# Define the range of thresholds
# thresholds <- seq(0, 100, by = 0.01)
thresholds <- c(0, sort(prediction_data$min_rank)[which(sort(prediction_data$min_rank) < par_max_rank_th)])


# Allocate memory for the results data frame
results <- data.frame(threshold = numeric(length(thresholds)),
                      odds_ratio = numeric(length(thresholds)),
                      lower_ci = numeric(length(thresholds)),
                      upper_ci = numeric(length(thresholds)),
                      p_value = numeric(length(thresholds)))


# Create a progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed",
  total = length(thresholds), clear = FALSE, width = 60
)

# Source of data
vec_pred_DT <- c()

# Loop over thresholds
for (i in seq_along(thresholds)) {
  threshold <- thresholds[i]
  
  # Add predicted_exp column based on threshold
  prediction_data$predicted_exp <- NA # Make sure it's empty
  prediction_data <- prediction_data %>%
    dplyr::mutate(predicted_exp = ifelse(min_rank <= threshold, 1, 0))
  
  # Create a contingency table
  contingency_table <- table(prediction_data$Target, prediction_data$predicted_exp)
  
  # Ensure contingency table has all cells (2x2 matrix)
  full_table <- matrix(0, nrow = 2, ncol = 2, 
                       dimnames = list(c("0", "1"), c("0", "1")))
  existing_rownames <- rownames(contingency_table)
  existing_colnames <- colnames(contingency_table)
  
  full_table[existing_rownames, existing_colnames] <- contingency_table
  
  
  # Perform Fisher's Exact Test
  fishers_test <- fisher.test(full_table)
  # Store results
  results$threshold[i] <- threshold
  results$odds_ratio[i] <- fishers_test$estimate
  results$lower_ci[i] <- fishers_test$conf.int[1]
  results$upper_ci[i] <- fishers_test$conf.int[2]
  results$p_value[i] <- fishers_test$p.value
  
  results$TP[i] <- full_table[2, 2]
  results$FP[i] <- full_table[1, 2]
  results$FN[i] <- full_table[2, 1]
  results$TN[i] <- full_table[1, 1]
  
  
  # Source of the DT
  temp_pred_DT <- prediction_data$Gene[which(prediction_data$Target == 1 & prediction_data$predicted_exp == 1)]
  new_pred_DT <- temp_pred_DT[which(!temp_pred_DT %in% vec_pred_DT)]
  vec_pred_DT <- temp_pred_DT
  
  results$new_TP_gene[i] <- new_pred_DT[1] # Take the first one if there are multiple that are added at the same time
  if (length(new_pred_DT) > 1) {
    print(paste(
      "Warning there are", 
      length(new_pred_DT),
      "new TP at once (", 
      paste(new_pred_DT, collapse = ", "),  # Separates elements with commas
      "), only the first one is considered for DT source"
    ))
  }
  
  # Update the progress bar
  pb$tick()
}

### Zoom parameter
par_optimal_F1 <- FALSE

### Filling NA of new TP Genes and adding DT source:
results <- merge(results, drug_targets_genes[, c("Gene", "Sum")], 
                 by.x = "new_TP_gene", by.y = "Gene", 
                 all.x = TRUE)
results <- results[order(results$threshold), ]

### Adding metrics 
results$precision <- results$TP / (results$TP + results$FP)
results$recall <- results$TP / (results$TP + results$FN)
results$TPR <- results$TP / (results$TP + results$FN)
results$FPR <- results$FP / (results$FP + results$TN)
results$StN <- results$TP / results$FP
results$F1 <- 2 * (results$precision * results$recall) / (results$precision + results$recall)

# Remove infinity and 0 cases
plot_results_pos <- results %>%
  dplyr::filter(!odds_ratio %in% c(Inf, 0))

# In case no TP are found in the given threshold range
found_match_bool <- !is.na(match(TRUE, !is.na(results$new_TP_gene)))

# Add first occurence of new TP gene to plotting (removed cause INF)
if(found_match_bool){
  plot_results_pos$new_TP_gene[1] <- results$new_TP_gene[which(!is.na(results$new_TP_gene))[1]]
  plot_results_pos$Sum[1] <- results$Sum[which(!is.na(results$Sum))[1]]
}else{
  print("No TP found in the given threshold range. No plot_results_pos will be created.")
}

### Start from 0
# Extract the first row
first_row <- plot_results_pos[1, ]

# Modify the threshold to 0
first_row$threshold <- 0
first_row$new_TP_gene <- NA
first_row$Sum <- NA

# Append the new row to the original dataframe
plot_results_pos <- rbind(first_row, plot_results_pos)

# Ensure the data is sorted by threshold
plot_results_pos <- plot_results_pos[order(plot_results_pos$threshold), ]

################################################
### Make dataframe for DT source ###############

target_source <- merge(prediction_data, drug_targets[,c("hgnc_symbol","Sum")], by.x = "Gene", by.y = "hgnc_symbol", all.x = T)

# Check for duplicate Genes in target_source
dup_target_source <- target_source[duplicated(target_source$Gene) | 
                                     duplicated(target_source$Gene, fromLast = TRUE), ]
cat("Duplicated rows in target_source after clean up:\n")
print(dup_target_source)

# Drop redundant columns
target_source <- target_source[,c("Gene","min_rank","source_min","Sum")]

# First, add the always-present columns
target_source <- target_source %>%
  mutate(
    p_rank_GWAS            = ifelse(Gene %in% GWAS_gene_space, 1, 0),
    p_rank_exome           = ifelse(Gene %in% Exome_gene_space, 1, 0),
    p_rank_eQTL_GWAS_blood = ifelse(Gene %in% eQTL_GWAS_blood_gene_space, 1, 0),
    p_rank_pQTL_PPP        = ifelse(Gene %in% pQTL_GWAS_gene_space, 1, 0)
  )

# Then, conditionally add p_rank_LM if pLM_gene_space exists
if (exists("pLM_gene_space")) {
  target_source <- target_source %>%
    mutate(p_rank_LM = ifelse(Gene %in% pLM_gene_space, 1, 0))
}

# And conditionally add p_rank_OT if OT_gene_space exists
if (exists("OT_gene_space")) {
  target_source <- target_source %>%
    mutate(p_rank_OT = ifelse(Gene %in% OT_gene_space, 1, 0))
}

################################################
### Export results and plot data ###############

### Export parameters naming
export_suffixes <- paste0("_", as.character(nb_db_for_DT), "DTS_", DTRAIT, "T")
if(par_TTD){
  export_suffixes <- paste0(export_suffixes, "_TTD")
}

if(par_full){
  export_suffixes <- paste0(export_suffixes, "_full")
}

if(par_GWAS_uni){
  export_suffixes <- paste0(export_suffixes, "_GWuni")
}
if(par_exome_uni){
  export_suffixes <- paste0(export_suffixes, "_Exomeuni")
}

if(par_PPP_uni){
  export_suffixes <- paste0(export_suffixes, "_PPPuni")
}

if("p_rank_pQTL_PPP" %in% INCLUDED_DATA){
  export_suffixes <- paste0(export_suffixes, "_pQTL")
}
if("p_rank_exome" %in% INCLUDED_DATA){
  export_suffixes <- paste0(export_suffixes, "_exome")
}
if("p_rank_eQTL_GWAS_blood" %in% INCLUDED_DATA){
  export_suffixes <- paste0(export_suffixes, "_eQTL")
}
if("p_rank_GWAS" %in% INCLUDED_DATA){
  export_suffixes <- paste0(export_suffixes, "_GWAS")
}
if("p_rank_pQTL_GWAS" %in% INCLUDED_DATA){
  export_suffixes <- paste0(export_suffixes, "_pQTLdC")
}

if(length(INCLUDED_DATA) > 1){
  export_suffixes <- paste0(export_suffixes, "_", par_method)
}

### Export target_source, plot_results_pos, results as rds files
saveRDS(target_source, file = paste0(path_to_output, TRAIT, export_suffixes, "_target_source.rds"))
if(found_match_bool){
  saveRDS(plot_results_pos, file = paste0(path_to_output, TRAIT, export_suffixes, "_plot_results_pos.rds"))
}
saveRDS(results, file = paste0(path_to_output, TRAIT, export_suffixes, "_results.rds"))

print("Results saved in folder, with TRAIT and suffix:")
print(paste0(path_to_output, TRAIT, export_suffixes))

  