################################################################################
### Minimal delta-difference benchmarking script ###############################
### These scripts have been adapted from original scripts with #################
### OpeanAI's GPT-5.4 model and rerun to check consistency. ####################
### Date: 2026-03-31 ###########################################################
### The output summary file contains the mean differences presented in #########
### the manuscript. ############################################################
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

################################################
### Editable parameters ########################

path_to_data <- "/.../"

traits <- c("LOAD_consortium", "AD", "AF", "ASTHMA",
            "BIPOLAR", "BMD", "CAD_consortium",
            "COPD_FIG_UKBB", "eGFR",
            "EPILEPSY", "ENDOMETRIOSIS_FIG_UKBB", "GLAUCOMA",
            "IBD_consortium", "IBS",
            "MDD", "MS", "OA",
            "PD", "PNEUMONIA_FIG_UKBB", "PSORIASIS_FIG_UKBB",
            "RA_consortium", "SCZ_consortium", "SBP",
            "STROKE", "T1D", "T2D",
            "TC", "VTE_UKBB", "LDL", "DBP")

file_nb_dts <- 2 # This is just used for file naming, it could be made cleaner
nb_DT_ttest <- c(2, 3) # Number of drug target threshold
par_threshold <- 100 # This should be kept to 100 to include all genes (this is kept as a legacy parameter for exploration)
par_VIP_removal <- FALSE # While we use VIP removal for cross trait analysis (i

baseline_method <- "GWAS"
method_files <- c(
  GWAS = "TTD_GWuni_GWAS",
  minimum = "TTD_GWuni_pQTL_exome_eQTL_GWAS_minimum",
  mean = "TTD_GWuni_pQTL_exome_eQTL_GWAS_mean",
  product = "TTD_GWuni_pQTL_exome_eQTL_GWAS_product",
  PCA = "TTD_GWuni_pQTL_exome_eQTL_GWAS_PCA",
  weightsum = "TTD_GWuni_pQTL_exome_eQTL_GWAS_weightsum"
)

################################################
### Files #####################################

drug_target_file <- file.path(path_to_data, "merged_drug_all_db_with_hgnc_TTD.tsv")
vip_file <- file.path(path_to_data, "extended_VIP_gene_list.txt")
summary_file <- file.path(path_to_data, "delta_diff_summary.tsv")
detail_file <- file.path(path_to_data, "delta_diff_detail.tsv")

################################################
### Helpers ###################################

build_target_source_path <- function(trait, method_name) {
  suffix <- unname(method_files[[method_name]])
  file.path(
    path_to_data,
    paste0(trait, "_", file_nb_dts, "DTS_", trait, "T_", suffix, "_target_source.rds")
  )
}

read_target_source <- function(trait, method_name, column_name) {
  readRDS(build_target_source_path(trait, method_name)) %>%
    dplyr::select(Gene, min_rank) %>%
    dplyr::rename(!!column_name := min_rank)
}

write_tsv <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

empty_summary_df <- function() {
  data.frame(
    trait = character(),
    method = character(),
    baseline_method = character(),
    file_nb_dts = integer(),
    nb_DT_ttest = integer(),
    par_threshold = numeric(),
    n_targets_before_vip = integer(),
    n_targets_after_vip = integer(),
    n_genes_tested = integer(),
    n_complete_pairs = integer(),
    diff_mean = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    p_value = numeric(),
    method_file = character(),
    baseline_file = character(),
    stringsAsFactors = FALSE
  )
}

empty_detail_df <- function() {
  data.frame(
    trait = character(),
    method = character(),
    baseline_method = character(),
    file_nb_dts = integer(),
    nb_DT_ttest = integer(),
    Gene = character(),
    min_rank_method = numeric(),
    min_rank_baseline = numeric(),
    delta = numeric(),
    paired_complete = logical(),
    stringsAsFactors = FALSE
  )
}

################################################
### Load inputs ################################

if (!baseline_method %in% names(method_files)) {
  stop("baseline_method must be one of names(method_files).")
}

if (!file.exists(drug_target_file)) {
  stop("Drug target file not found: ", drug_target_file)
}

drug_target_df <- read.table(
  drug_target_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

vip_genes <- character(0)
if (par_VIP_removal) {
  if (file.exists(vip_file)) {
    vip_genes <- scan(vip_file, what = "character", quiet = TRUE)
  } else {
    warning("VIP removal enabled but file not found: ", vip_file)
  }
}

################################################
### Compute delta differences ##################

methods_to_compare <- setdiff(names(method_files), baseline_method)
summary_results <- list()
detail_results <- list()

for (current_trait in traits) {
  baseline_path <- build_target_source_path(current_trait, baseline_method)

  for (method_name in methods_to_compare) {
    method_path <- build_target_source_path(current_trait, method_name)

    if (!file.exists(method_path) || !file.exists(baseline_path)) {
      message("Skipping ", current_trait, " / ", method_name, ": missing input file(s).")
      next
    }

    method_data <- read_target_source(current_trait, method_name, "min_rank_method")
    baseline_data <- read_target_source(current_trait, baseline_method, "min_rank_baseline")
    merged_data <- merge(method_data, baseline_data, by = "Gene", all.x = TRUE)

    for (current_nb in nb_DT_ttest) {
      drug_target_ttest <- drug_target_df %>%
        dplyr::filter(trait == current_trait, Sum >= current_nb, hgnc_symbol != "") %>%
        dplyr::pull(hgnc_symbol) %>%
        unique()

      n_targets_before_vip <- length(drug_target_ttest)

      if (par_VIP_removal && length(vip_genes) > 0) {
        drug_target_ttest <- setdiff(drug_target_ttest, vip_genes)
      }

      dt_data <- merged_data %>%
        dplyr::filter(
          (min_rank_method <= par_threshold | min_rank_baseline <= par_threshold) &
            Gene %in% drug_target_ttest
        )

      if (nrow(dt_data) < 2) {
        message(
          "Skipping ", current_trait, " / ", method_name,
          " / nb_DT=", current_nb, ": not enough drug target rows."
        )
        next
      }

      complete_pairs <- stats::complete.cases(
        dt_data$min_rank_method,
        dt_data$min_rank_baseline
      )

      if (sum(complete_pairs) < 2) {
        message(
          "Skipping ", current_trait, " / ", method_name,
          " / nb_DT=", current_nb, ": not enough complete pairs."
        )
        next
      }

      diff_mean <- mean(dt_data$min_rank_method, na.rm = TRUE) -
        mean(dt_data$min_rank_baseline, na.rm = TRUE)

      ttest_result <- try(
        stats::t.test(
          dt_data$min_rank_method,
          dt_data$min_rank_baseline,
          paired = TRUE,
          alternative = "less"
        ),
        silent = TRUE
      )

      if (inherits(ttest_result, "try-error")) {
        message(
          "t-test failed for ",
          current_trait, " / ", method_name, " / nb_DT=", current_nb
        )
        next
      }

      upper_ci <- as.numeric(ttest_result$conf.int[2])
      lower_ci <- diff_mean - abs(diff_mean - upper_ci)

      summary_results[[length(summary_results) + 1]] <- data.frame(
        trait = current_trait,
        method = method_name,
        baseline_method = baseline_method,
        file_nb_dts = file_nb_dts,
        nb_DT_ttest = current_nb,
        par_threshold = par_threshold,
        n_targets_before_vip = n_targets_before_vip,
        n_targets_after_vip = length(drug_target_ttest),
        n_genes_tested = nrow(dt_data),
        n_complete_pairs = sum(complete_pairs),
        diff_mean = diff_mean,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        p_value = ttest_result$p.value,
        method_file = method_path,
        baseline_file = baseline_path,
        stringsAsFactors = FALSE
      )

      detail_results[[length(detail_results) + 1]] <- dt_data %>%
        dplyr::mutate(
          trait = current_trait,
          method = method_name,
          baseline_method = baseline_method,
          file_nb_dts = file_nb_dts,
          nb_DT_ttest = current_nb,
          delta = min_rank_method - min_rank_baseline,
          paired_complete = complete_pairs
        ) %>%
        dplyr::select(
          trait,
          method,
          baseline_method,
          file_nb_dts,
          nb_DT_ttest,
          Gene,
          min_rank_method,
          min_rank_baseline,
          delta,
          paired_complete
        )
    }
  }
}

summary_df <- if (length(summary_results) > 0) {
  dplyr::bind_rows(summary_results)
} else {
  empty_summary_df()
}

detail_df <- if (length(detail_results) > 0) {
  dplyr::bind_rows(detail_results)
} else {
  empty_detail_df()
}

write_tsv(summary_df %>% dplyr::arrange(method, trait, nb_DT_ttest), summary_file)
write_tsv(detail_df %>% dplyr::arrange(method, trait, nb_DT_ttest, Gene), detail_file)

message("Wrote delta summary: ", summary_file)
message("Wrote delta detail: ", detail_file)
