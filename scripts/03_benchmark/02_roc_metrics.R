################################################################################
### Minimal ROC benchmarking script ############################################
### These scripts have been adapted from original scripts with #################
### OpeanAI's GPT-5.4 model and rerun to check consistency. ####################
### Date: 2026-03-31 ###########################################################
### The output summary file contains the AUC values presented in the ###########
### the manuscript. ############################################################
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
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
            "TC", "VTE_UKBB", "DBP", "LDL")

file_nb_dts <- 2 # This is just used for file naming, it could be made cleaner
par_nb_DT <- 3
boot_n <- 2000

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

summary_file <- file.path(path_to_data, "roc_auc_summary.tsv")
curve_file <- file.path(path_to_data, "roc_curve_points.tsv")
pdiff_file <- file.path(path_to_data, "roc_pdiff.tsv")

################################################
### Helpers ###################################

build_target_source_path <- function(trait, method_name) {
  suffix <- unname(method_files[[method_name]])
  file.path(
    path_to_data,
    paste0(trait, "_", file_nb_dts, "DTS_", trait, "T_", suffix, "_target_source.rds")
  )
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

empty_auc_df <- function() {
  data.frame(
    trait = character(),
    method = character(),
    file_nb_dts = integer(),
    par_nb_DT = integer(),
    n_genes = integer(),
    n_cases = integer(),
    n_controls = integer(),
    AUC = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    file_path = character(),
    stringsAsFactors = FALSE
  )
}

empty_curve_df <- function() {
  data.frame(
    trait = character(),
    method = character(),
    file_nb_dts = integer(),
    par_nb_DT = integer(),
    threshold = numeric(),
    specificity = numeric(),
    sensitivity = numeric(),
    FPR = numeric(),
    TPR = numeric(),
    stringsAsFactors = FALSE
  )
}

################################################
### Compute ROC metrics ########################

if (!baseline_method %in% names(method_files)) {
  stop("baseline_method must be one of names(method_files).")
}

summary_results <- list()
curve_results <- list()

for (method_name in names(method_files)) {
  for (current_trait in traits) {
    file_path <- build_target_source_path(current_trait, method_name)

    if (!file.exists(file_path)) {
      message("Skipping ", current_trait, " / ", method_name, ": missing input file.")
      next
    }

    current_df <- readRDS(file_path)

    if (!all(c("Gene", "min_rank", "Sum") %in% names(current_df))) {
      message("Skipping ", current_trait, " / ", method_name, ": required columns missing.")
      next
    }

    current_df <- current_df %>%
      dplyr::select(Gene, min_rank, Sum) %>%
      dplyr::mutate(target_label = ifelse(!is.na(Sum) & Sum >= par_nb_DT, 1L, 0L))

    if (length(unique(current_df$target_label)) < 2) {
      message(
        "Skipping ", current_trait, " / ", method_name,
        ": ROC needs both cases and controls."
      )
      next
    }

    roc_obj <- pROC::roc(
      response = current_df$target_label,
      predictor = current_df$min_rank,
      levels = c(0, 1),
      direction = ">",
      quiet = TRUE
    )

    auc_val <- as.numeric(pROC::auc(roc_obj))
    ci_vals <- as.numeric(
      pROC::ci.auc(
        roc_obj,
        method = "bootstrap",
        boot.n = boot_n,
        conf.level = 0.95
      )
    )

    summary_results[[length(summary_results) + 1]] <- data.frame(
      trait = current_trait,
      method = method_name,
      file_nb_dts = file_nb_dts,
      par_nb_DT = par_nb_DT,
      n_genes = nrow(current_df),
      n_cases = sum(current_df$target_label == 1L),
      n_controls = sum(current_df$target_label == 0L),
      AUC = auc_val,
      CI_lower = ci_vals[1],
      CI_upper = ci_vals[3],
      file_path = file_path,
      stringsAsFactors = FALSE
    )

    coords_df <- as.data.frame(
      pROC::coords(
        roc_obj,
        x = "all",
        ret = c("threshold", "specificity", "sensitivity"),
        transpose = FALSE
      )
    )

    curve_results[[length(curve_results) + 1]] <- coords_df %>%
      dplyr::mutate(
        trait = current_trait,
        method = method_name,
        file_nb_dts = file_nb_dts,
        par_nb_DT = par_nb_DT,
        FPR = 1 - specificity,
        TPR = sensitivity
      ) %>%
      dplyr::select(
        trait,
        method,
        file_nb_dts,
        par_nb_DT,
        threshold,
        specificity,
        sensitivity,
        FPR,
        TPR
      )
  }
}

summary_df <- if (length(summary_results) > 0) {
  dplyr::bind_rows(summary_results)
} else {
  empty_auc_df()
}

curve_df <- if (length(curve_results) > 0) {
  dplyr::bind_rows(curve_results)
} else {
  empty_curve_df()
}

pdiff_df <- summary_df %>%
  dplyr::mutate(se = (CI_upper - CI_lower) / (2 * 1.96)) %>%
  dplyr::left_join(
    summary_df %>%
      dplyr::filter(method == baseline_method) %>%
      dplyr::transmute(
        trait,
        AUC_baseline = AUC,
        se_baseline = (CI_upper - CI_lower) / (2 * 1.96),
        baseline_file = file_path
      ),
    by = "trait"
  ) %>%
  dplyr::mutate(
    baseline_method = baseline_method,
    pdiff = dplyr::case_when(
      method == baseline_method ~ NA_real_,
      is.na(AUC_baseline) ~ NA_real_,
      TRUE ~ 2 * stats::pnorm(-abs((AUC - AUC_baseline) / sqrt(se^2 + se_baseline^2)))
    )
  ) %>%
  dplyr::select(
    trait,
    method,
    baseline_method,
    AUC,
    CI_lower,
    CI_upper,
    AUC_baseline,
    se,
    se_baseline,
    pdiff,
    file_path,
    baseline_file
  )

write_tsv(summary_df %>% dplyr::arrange(method, trait), summary_file)
write_tsv(curve_df %>% dplyr::arrange(method, trait, FPR, TPR), curve_file)
write_tsv(pdiff_df %>% dplyr::arrange(method, trait), pdiff_file)

message("Wrote ROC summary: ", summary_file)
message("Wrote ROC curve points: ", curve_file)
message("Wrote ROC pdiff: ", pdiff_file)
