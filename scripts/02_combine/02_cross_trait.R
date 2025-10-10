###################################################################################
### Script to get bent correlation matrix and compare cross-trait predictions #####
### Author: S. Moix, DBC, UNIL, 2025 with partial assistance from LLMs ############
### Date: October 2025 ############################################################
###################################################################################

################################################
### Load libraries #############################
library(dplyr)
library(tidyr)
library(mbend)
library(MASS)
library(conflicted)

################################################
### Set path and parameters ####################
path_to_data <- "/.../" # Main folder that contains the data (i.e., results from 01_aggregation.R and other files)
path_to_acronyms <- "/.../Exome_data.csv" # File with trait acronyms (https://github.com/masadler/DrugTargetMethodComparison/tree/main/data)
path_to_rg <- "/.../" # Folder with genetic correlation results (see 04_rg_ldsc)


traits <- c("LOAD_consortium", "AD", "AF","ASTHMA",
            "BIPOLAR", "BMD", "CAD_consortium",
            "COPD_FIG_UKBB", "eGFR",
            "EPILEPSY", "ENDOMETRIOSIS_FIG_UKBB", "GLAUCOMA",
            "IBD_consortium", "IBS",
            "MDD", "MS", "OA",
            "PD", "PNEUMONIA_FIG_UKBB", "PSORIASIS_FIG_UKBB",
            "RA_consortium", "SCZ_consortium", "SBP",
            "STROKE", "T1D", "T2D",
            "TC", "VTE_UKBB") 


# Parameters (modify these as needed)
nb_DT_ttest <- 3 # Minimum number of drug targets to perform t-test (defines lenient or moderate set...)
par_threshold <- 100 # Percentile threshold for filtering genes (e.g., 100 for all genes, 50 for top 50%)
### Suffixes from below  depend on the parameters used in 01_aggregation.R
export_suffixes <- "TTD_GWuni_pQTL_exome_eQTL_GWAS_minimum" # Suffix for the ranking results files
par_VIP_removal <- TRUE # Whether to remove VIP genes from the drug target list
nb_sim_random <- 100 # Number of random simulations to run to obtain the random control results

################################################
### Load data ##################################

acronym_dict <- read.csv(path_to_acronyms, header = TRUE, sep = ",") 
acronym_dict <- acronym_dict %>%
  mutate(Acronym = gsub("SBP_meta", "SBP", Acronym),
         Acronym = gsub("DBP_meta", "DBP", Acronym))

acronym_map <- setNames(as.character(acronym_dict$Abbreviation), acronym_dict$Acronym)
ext_VIP_gene_list <- scan(file.path(path_to_data, "extended_VIP_gene_list.txt"), what = "character") # See manuscript supplemental data
gene_universe <- readLines(file.path(path_to_data, "gene_universe_rbs.txt")) # See data folder

################################################
### Get bent matrix from rg ####################

set.seed(12345)
# Shuffle genes
gene_universe <- as.vector(sample(gene_universe))

################################################
### Load genetic correlation data ##############

grg_summary <- read.table(file.path(path_to_rg, "pw_genetic_correlation_results.tsv"),
                          sep = "\t", header = TRUE)


# Clean trait names by removing ".sumstats.gz" and compute -log10(p)
grg_summary <- grg_summary %>%
  mutate(trait1 = gsub("\\.sumstats\\.gz$", "", p1),
         trait2 = gsub("\\.sumstats\\.gz$", "", p2),
         logp = -log10(p))

grg_summary <- grg_summary %>%
  mutate(
    trait1 = gsub("SBP_meta", "SBP", trait1),
    trait2 = gsub("SBP_meta", "SBP", trait2),
    trait1 = gsub("DBP_meta", "DBP", trait1),
    trait2 = gsub("DBP_meta", "DBP", trait2)
  )

grg_summary <- grg_summary %>%
  dplyr::filter(trait1 %in% traits & trait2 %in% traits)


### Populate heritability data frame based on the mean across traits
variance_df <- data.frame(trait = traits,
                         h2 = NA,
                         h2_se = NA,
                         stringsAsFactors = FALSE)

for(var_trait in traits){
  temp_df_1 <- grg_summary %>%
    dplyr::select(c("trait1", "h2_p1", "h2_se_p1")) %>%
    dplyr::filter(trait1 == var_trait) %>%
    dplyr::rename(trait = trait1, h2 = h2_p1, h2_se = h2_se_p1)
  
  temp_df_2 <- grg_summary %>%
    dplyr::select(c("trait2", "h2_p2", "h2_se_p2")) %>%
    dplyr::filter(trait2 == var_trait) %>%
    dplyr::rename(trait = trait2, h2 = h2_p2, h2_se = h2_se_p2)
  
  temp_df <- rbind(temp_df_1, temp_df_2)
  
  temp_df <- temp_df %>%
    dplyr::summarise(h2 = mean(h2),
                     h2_se = mean(h2_se))
  
  variance_df$h2[variance_df$trait == var_trait] <- temp_df$h2
  variance_df$h2_se[variance_df$trait == var_trait] <- temp_df$h2_se
}


grg_summary <- grg_summary %>%
  dplyr::select(c("trait1", "trait2", "gcov", "gcov_se"))

################################################
### Bend matrix and make mvrnorm ###############

# Create an empty square matrix with these traits as both row and column names
gcov_matrix <- matrix(NA, nrow = length(traits), ncol = length(traits),
                      dimnames = list(traits, traits))

# Fill the matrix with the correlation values and their symmetric counterparts
for(i in seq_len(nrow(grg_summary))) {
  t1 <- grg_summary$trait1[i]
  t2 <- grg_summary$trait2[i]
  r_val <- grg_summary$gcov[i]
  
  gcov_matrix[t1, t2] <- r_val
  gcov_matrix[t2, t1] <- r_val  # ensure symmetry
}

### Slower but safer 
for (i in seq_along(traits)) {
  trait_name <- traits[i]
  gcov_matrix[trait_name, trait_name] <- variance_df$h2[variance_df$trait == trait_name]
}

# Print the resulting correlation matrix
print(gcov_matrix)

### Same for gcov_se
gcov_se_matrix <- matrix(NA, nrow = length(traits), ncol = length(traits),
                         dimnames = list(traits, traits))

# Fill the matrix with the standard error values and their symmetric counterparts
for(i in seq_len(nrow(grg_summary))) {
  t1 <- grg_summary$trait1[i]
  t2 <- grg_summary$trait2[i]
  r_val <- grg_summary$gcov_se[i]
  
  gcov_se_matrix[t1, t2] <- r_val
  gcov_se_matrix[t2, t1] <- r_val  # ensure symmetry
}

for (i in seq_along(traits)) {
  trait_name <- traits[i]
  gcov_se_matrix[trait_name, trait_name] <- variance_df$h2_se[variance_df$trait == trait_name]
}

print(gcov_se_matrix)

# Bend correlation matrix to make it positive-definite
# bent_gcov_matrix_results <- mbend::bend(gcov_matrix)

# Calculate the weight matrix as the inverse squared standard errors.
weight_matrix <- 1 / (gcov_se_matrix^2)

# Bend the covariance matrix using the weight matrix.
bent_gcov_matrix_results <- mbend::bend(gcov_matrix, weight_matrix)
bent_gcov_matrix <- bent_gcov_matrix_results$bent


# Create the diagonal scaling matrix D.
# The diagonal elements of D are the inverse square roots of the diagonal of the bent matrix.
D <- diag(1 / sqrt(diag(bent_gcov_matrix))) # diag((diag(bent_gcov_matrix)**-0.5))

# Compute the standardized (or scaled) bent matrix.
# This rescales the bent covariance matrix so that its diagonal becomes 1.
final_matrix <- D %*% bent_gcov_matrix %*% D


rownames(final_matrix) <- rownames(bent_gcov_matrix)
colnames(final_matrix) <- colnames(bent_gcov_matrix)



################################################
### Compare RBS (Random Best Self) #############

### Summary of the approach: ###################
# trait_DT stands for the trait that who's drug targets we want to predict
# trait_score is the priorization scores based on the minimum approach for a given trait (for each trait_DT we have 28 trait_scores)
#       1) A t-test between drug targets and non-drug targets is performed for each trait_score (groups based on the drug target list for trait_DT)
#       2) The "Best" trait_score is selected as the one with the lowest p-value among those with a negative difference in means
#       3) The "Self" trait_score is also selected (i.e., where trait_score == trait_DT)
#       4) A random analysis is performed 
#           - For each replicate (nb_sim_random), new random z-stats are generated based on the bent correlation matrix
#           - For each random_trait_score, a t-test is performed between drug targets and non-drug targets
#           - The best random_trait_score is selected as the one with the lowest p-value among those with a negative difference in means
#           - The mean across each replicate is taken to obtain the random control results (one-sample t-test against mu = 0)
#       5) Results are for each trait_DT, for the "Best", "Self" and "Random" measures
################################################

n_sim <- length(gene_universe)
mu <- rep(0, ncol(bent_gcov_matrix))

drug_target_df <- read.table(file.path(path_to_data, "merged_drug_all_db_with_hgnc_TTD.tsv"), # Drug target information
                                  sep = "\t", header = TRUE)

panel_data <- list()
for(trait_DT in traits){
  
  # Build results based on the min_rank t-tests for each trait_score
  trait_results <- data.frame(trait_score = character(),
                              diff_mean = numeric(),
                              lower_ci = numeric(),
                              upper_ci = numeric(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)
  
   # Load drug target data for the current trait_DT
  drug_target_ttest <- drug_target_df %>%
    dplyr::filter(trait == trait_DT, Sum >= nb_DT_ttest) %>%
    dplyr::pull(hgnc_symbol)
  
  if(par_VIP_removal){
    drug_target_ttest <- setdiff(drug_target_ttest, ext_VIP_gene_list)
  }
    
  for(trait_score in traits){
    # Here just load the minimum ranking results for the given trait_score (thus number of DTS doesn't matter, here 2 is used)
    file_name <- paste0(trait_score, "_", 2, "DTS_", trait_score, "T_", export_suffixes, "_target_source.rds")
    file_path <- file.path(path_to_data, file_name)
    
    if (!file.exists(file_path)) {
      message(paste("File not found for trait:", trait_score, "and", trait_DT))
      next  # Skip if file doesn't exist
    }
    
    # Load the data and filter by par_threshold
    min_scores_ttest <- readRDS(file_path) %>%
      dplyr::select(Gene, min_rank) %>%
      dplyr::filter(min_rank <= par_threshold)
    
   
    # Create grouping variable based on whether Gene is in drug_target_ttest
    min_scores_ttest <- min_scores_ttest %>%
      dplyr::mutate(Group = dplyr::if_else(Gene %in% drug_target_ttest, "Drug Target", "Non-Drug Target"))
    
    if(length(unique(min_scores_ttest$Group)) < 2){
      warning(paste("Not enough groups for t-test for trait_score:", trait_score, "and trait_DT:", trait_DT))
      next
    }
    
    # Run the t-test for min_rank values
    ttest_result <- try(t.test(min_rank ~ Group, data = min_scores_ttest), silent = TRUE)
    if(inherits(ttest_result, "try-error")){
      message("t-test failed for trait_score: ", trait_score, " and trait_DT: ", trait_DT)
      next
    }
    
    estimates <- ttest_result$estimate
    diff_mean <- estimates[["mean in group Drug Target"]] - estimates[["mean in group Non-Drug Target"]]
    ci <- ttest_result$conf.int
    p_val <- ttest_result$p.value
    
    trait_results <- rbind(trait_results,
                           data.frame(trait_score = trait_score,
                                      diff_mean = diff_mean,
                                      lower_ci = ci[1],
                                      upper_ci = ci[2],
                                      p_value = p_val,
                                      stringsAsFactors = FALSE))
  }
  
  # From trait_results, select the "Best" measure: among those with a negative diff_mean, choose the one with the lowest p-value.
  valid_best <- trait_results %>% dplyr::filter(diff_mean < 0)
  if(nrow(valid_best) > 0){
    best_row <- valid_best %>% dplyr::arrange(p_value) %>% dplyr::slice(1)
    best_row$measure_type <- "Best"
  } else {
    best_row <- data.frame(trait_score = NA,
                           diff_mean = 0,
                           lower_ci = 0,
                           upper_ci = 0,
                           p_value = NA,
                           measure_type = "Best",
                           stringsAsFactors = FALSE)
  }
  
  # Get the "Self" measure (where trait_score equals trait_DT)
  self_row <- trait_results %>% dplyr::filter(trait_score == trait_DT)
  if(nrow(self_row) == 0){
    self_row <- data.frame(trait_score = trait_DT,
                           diff_mean = 0,
                           lower_ci = 0,
                           upper_ci = 0,
                           p_value = NA,
                           measure_type = "Self",
                           stringsAsFactors = FALSE)
  } else {
    self_row$measure_type <- "Self"
  }
  
  # Now perform the random analysis based on random_df_percentile
  
  # Initialize an empty data frame to store simulation results:
  all_random_results <- data.frame(trait_score = character(), 
                                   diff_mean = numeric(), 
                                   lower_ci = numeric(), 
                                   upper_ci = numeric(), 
                                   p_value = numeric(), 
                                   measure_type = character(),
                                   stringsAsFactors = FALSE)
  
  ### Run many random simulations and take average of best
  for(i in 1:nb_sim_random){
    
    # For each replicate, initialize a temporary container for this simulation's results:
    random_trait_results <- data.frame(trait_score = character(),
                                       diff_mean = numeric(),
                                       lower_ci = numeric(),
                                       upper_ci = numeric(),
                                       p_value = numeric(),
                                       stringsAsFactors = FALSE)

      # Generate new random z-stats each time
      z_stats <- MASS::mvrnorm(n = n_sim, mu = mu, Sigma = final_matrix)
      colnames(z_stats) <- colnames(bent_gcov_matrix)
      rownames(z_stats) <- gene_universe
      random_df <- as.data.frame(z_stats)
      
      # Convert z_stats into percentiles
      random_df_percentile <- random_df %>%
        dplyr::mutate(across(everything(), ~ {
          f <- ecdf(.x)
          100 * f(.x)
        }))                                  
    
    # Loop over all 28 trait_scores:
    for(random_trait_score in traits){
      # Check if the trait_score exists in your generated data:
      if(!random_trait_score %in% colnames(random_df_percentile)){
        next
      }
      
      # Subset and create temporary results
      random_values <- random_df_percentile[, random_trait_score, drop = FALSE]
      random_temp <- data.frame(Gene = rownames(random_values),
                                score = random_values[[1]])
      
      # Apply filtering and grouping
      random_temp <- random_temp %>%
        dplyr::filter(score <= par_threshold) %>%
        dplyr::mutate(Group = dplyr::if_else(Gene %in% drug_target_ttest, "Drug Target", "Non-Drug Target"))
      
      if(length(unique(random_temp$Group)) < 2){
        message(paste("Not enough groups for t-test for random trait_score:", random_trait_score, "and trait_DT:", trait_DT))
        next
      }
      
      # Perform t-test
      ttest_random <- try(t.test(score ~ Group, data = random_temp), silent = TRUE)
      if(inherits(ttest_random, "try-error")){
        message(paste("t-test failed for random trait_score:", random_trait_score, "and trait_DT:", trait_DT))
        next
      }
      
      estimates_random <- ttest_random$estimate
      diff_mean_random <- estimates_random[["mean in group Drug Target"]] - estimates_random[["mean in group Non-Drug Target"]]
      ci_random <- ttest_random$conf.int
      p_val_random <- ttest_random$p.value
      
      # Append result for this trait_score into temporary container
      random_trait_results <- rbind(random_trait_results,
                                    data.frame(trait_score = random_trait_score,
                                               diff_mean = diff_mean_random,
                                               lower_ci = ci_random[1],
                                               upper_ci = ci_random[2],
                                               p_value = p_val_random,
                                               stringsAsFactors = FALSE))
    } # End inner loop over random_trait_score
    
    # After processing all trait_scores in this replicate, select the best one
    if(nrow(random_trait_results) > 0){
      valid_random_best <- random_trait_results %>% dplyr::filter(diff_mean < 0)
      if(nrow(valid_random_best) > 0){
        best_random_row <- valid_random_best %>% dplyr::arrange(p_value) %>% dplyr::slice(1)
        best_random_row$measure_type <- "Random"
      } else {
        best_random_row <- data.frame(trait_score = NA,
                                      diff_mean = 0,
                                      lower_ci = 0,
                                      upper_ci = 0,
                                      p_value = NA,
                                      measure_type = "Random",
                                      stringsAsFactors = FALSE)
      }
    } else {
      best_random_row <- data.frame(trait_score = NA,
                                    diff_mean = 0,
                                    lower_ci = 0,
                                    upper_ci = 0,
                                    p_value = NA,
                                    measure_type = "Random",
                                    stringsAsFactors = FALSE)
    }
    
    # Overwrite trait_score to "Random" in best row
    best_random_row$trait_score <- "Random"
    
    # Append only this replicate's best outcome to the master results
    all_random_results <- rbind(all_random_results, best_random_row)
    
  }
  
  # Compute overall statistics from the simulations:
  mean_diff_mean <- mean(all_random_results$diff_mean)
  
  # Use a one-sample t-test to compare against 0:
  t_test_result <- t.test(all_random_results$diff_mean, mu = 0)
  agg_p_value <- t_test_result$p.value
  agg_ci_lower <- t_test_result$conf.int[1]
  agg_ci_upper <- t_test_result$conf.int[2]
  
  # Combine these into a summary result data frame:
  agg_results <- data.frame(trait_score = "Random",
                            diff_mean = mean_diff_mean,
                            lower_ci = agg_ci_lower,
                            upper_ci = agg_ci_upper,
                            p_value = agg_p_value,
                            measure_type = "Random",
                            stringsAsFactors = FALSE)
  
  # Combine Best, Self, and Random rows (three entries) for the current trait_DT
  panel_df <- rbind(best_row, self_row, agg_results)
  panel_df$trait_DT <- trait_DT
  
  # Store the panel data for the current trait_DT in the list
  panel_data[[trait_DT]] <- panel_df
  print(paste0("Done for: ", trait_DT))
}

# Combine all panel data into one data frame with 3 rows per trait_DT
all_panel_data <- do.call(rbind, panel_data)

################################################
### Get pdiff ##################################

# First, add an estimated standard error to each row
all_panel_data <- all_panel_data %>%
  mutate(se = (upper_ci - lower_ci) / (2 * 1.96))

# Compute pdiff only from Best and Self rows
diff_data <- all_panel_data %>%
  dplyr::filter(measure_type %in% c("Best", "Self")) %>%
  dplyr::select(trait_DT, measure_type, diff_mean, se) %>%
  tidyr::pivot_wider(names_from = measure_type, values_from = c(diff_mean, se)) %>%
  dplyr::mutate(pdiff = 2 * pnorm(-abs((diff_mean_Self - diff_mean_Best) /
                                         sqrt(se_Self^2 + se_Best^2))))


# Merge the computed pdiff (at the trait_DT level) into all_panel_data
all_panel_data$shape_symbol <- ifelse(is.na(all_panel_data$p_value), 92, 16)
all_panel_data <- merge(all_panel_data, diff_data[, c("trait_DT", "pdiff")], by = "trait_DT", all.x = TRUE)

################################################
### Change labels ##############################

# Update trait_score and trait_DT labels based on acronym_map (if available)
all_panel_data$trait_score <- sapply(as.character(all_panel_data$trait_score), function(x) {
  if(x %in% names(acronym_map)) acronym_map[x] else x
})

all_panel_data$trait_DT <- sapply(as.character(all_panel_data$trait_DT), function(x) {
  if(x %in% names(acronym_map)) acronym_map[x] else x
})


################################################
### Adjusted pdiff Best vs Self ################

adjusted_best <- all_panel_data %>%
  dplyr::filter(measure_type %in% c("Best", "Random")) %>% 
  dplyr::select(trait_DT, measure_type, diff_mean, se) %>%
  pivot_wider(
    names_from   = measure_type,
    values_from  = c(diff_mean, se)
  ) %>%
  transmute(
    trait_DT,
    measure_type = "Adjusted_best",
    diff_mean    = diff_mean_Best - diff_mean_Random,
    se           = se_Best
  )


adjusted_best <- bind_rows(all_panel_data[which(all_panel_data$measure_type == "Self"),
                                          c("trait_DT", "measure_type", "diff_mean","se")], adjusted_best)

pdiff_adj <- adjusted_best %>%
  dplyr::filter(measure_type %in% c("Self", "Adjusted_best")) %>%
  dplyr::select(trait_DT, measure_type, diff_mean, se) %>%
  pivot_wider(
    names_from   = measure_type,
    values_from  = c(diff_mean, se)
  ) %>%
  transmute(
    trait_DT,
    pdiff_adjusted = 2 * pnorm(
      -abs(
        (diff_mean_Self - diff_mean_Adjusted_best) /
          sqrt(se_Self^2 + se_Adjusted_best^2)
      )
    )
  )

# Merge the new pdiff_adjusted back in
all_panel_data <- all_panel_data %>%
  left_join(pdiff_adj, by = "trait_DT")

################################################
### Export data as rds #########################

VIP_export <- ifelse(par_VIP_removal, "No_VIP", "VIP")
saveRDS(all_panel_data, file.path(path_to_data, paste0("all_panel_data_", VIP_export, "_TH", par_threshold, "_nbDT", nb_DT_ttest, "_nbsimR", nb_sim_random,".rds")))
print(paste("File saved to:", file.path(path_to_data, paste0("all_panel_data_", VIP_export, "_TH", par_threshold, "_nbDT", nb_DT_ttest, "_nbsimR", nb_sim_random,".rds"))))
