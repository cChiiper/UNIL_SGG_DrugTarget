################################################################################
### Minimal Jaccard-index computing script #####################################
### These scripts have been adapted from original scripts with #################
### OpeanAI's GPT-5.4 model and rerun to check consistency. ####################
### Date: 2026-03-31 ###########################################################
### The output summary file contains pairwise trait Jaccard indices ############
### used in the manuscript. ####################################################
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

nb_DT_values <- c(2, 3)

################################################
### Files #####################################

drug_target_file <- file.path(path_to_data, "merged_drug_all_db_with_hgnc_TTD.tsv")
summary_file <- file.path(path_to_data, "jaccard_pairwise_summary.tsv")

################################################
### Helpers ###################################

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

jaccard_stats <- function(genes_1, genes_2) {
  overlap_genes <- intersect(genes_1, genes_2)
  union_genes <- union(genes_1, genes_2)

  overlap_n <- length(overlap_genes)
  union_n <- length(union_genes)

  if (union_n == 0) {
    return(list(jaccard_index = NA_real_, overlap = overlap_n, union = union_n))
  }

  list(
    jaccard_index = overlap_n / union_n,
    overlap = overlap_n,
    union = union_n
  )
}

empty_summary_df <- function() {
  data.frame(
    trait1 = character(),
    trait2 = character(),
    nb_DT = integer(),
    jaccard_index = numeric(),
    overlap = integer(),
    union = integer(),
    n_trait1 = integer(),
    n_trait2 = integer(),
    stringsAsFactors = FALSE
  )
}

################################################
### Load inputs ################################

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

required_columns <- c("trait", "hgnc_symbol", "Sum")
missing_columns <- setdiff(required_columns, names(drug_target_df))

if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

################################################
### Compute pairwise Jaccard indices ###########

pairwise_results <- list()

for (current_nb_DT in nb_DT_values) {
  filtered_dt <- drug_target_df %>%
    dplyr::filter(
      trait %in% traits,
      Sum >= current_nb_DT,
      !is.na(hgnc_symbol),
      hgnc_symbol != ""
    ) %>%
    dplyr::group_by(trait, hgnc_symbol) %>%
    dplyr::slice_max(order_by = Sum, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  gene_sets <- setNames(vector("list", length(traits)), traits)
  for (current_trait in traits) {
    gene_sets[[current_trait]] <- filtered_dt$hgnc_symbol[filtered_dt$trait == current_trait] %>%
      unique()
  }

  trait_pairs <- utils::combn(traits, 2, simplify = FALSE)

  for (current_pair in trait_pairs) {
    trait1 <- current_pair[1]
    trait2 <- current_pair[2]

    genes_trait1 <- gene_sets[[trait1]]
    genes_trait2 <- gene_sets[[trait2]]
    current_stats <- jaccard_stats(genes_trait1, genes_trait2)

    pairwise_results[[length(pairwise_results) + 1]] <- data.frame(
      trait1 = trait1,
      trait2 = trait2,
      nb_DT = current_nb_DT,
      jaccard_index = current_stats$jaccard_index,
      overlap = as.integer(current_stats$overlap),
      union = as.integer(current_stats$union),
      n_trait1 = as.integer(length(genes_trait1)),
      n_trait2 = as.integer(length(genes_trait2)),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- if (length(pairwise_results) > 0) {
  dplyr::bind_rows(pairwise_results)
} else {
  empty_summary_df()
}

write_tsv(
  summary_df %>% dplyr::arrange(nb_DT, trait1, trait2),
  summary_file
)

message("Wrote Jaccard summary: ", summary_file)
