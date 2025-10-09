############################################################################################################
### Script to run TwoSampleMR analysis #####################################################################
### Author: S. Moix, DBC, UNIL, 2025 #######################################################################
### Date: October 2025 #####################################################################################
############################################################################################################


################################################
### Libraries ##################################
library(TwoSampleMR) #TwoSampleMR version 0.5.6
library(ggplot2)
library(dplyr)
library(gridExtra)
library(data.table)
library(tidyr)
library(conflicted)
library(argparse)

################################################
### Get arguments with argparse ################

# Create argument parser
parser <- ArgumentParser(description = 'strand-files')

# Add arguments
parser$add_argument('-mr_f', '--mr_data', help = 'MR data frame', required = TRUE)
parser$add_argument('-hla_f', '--HLA_snps', help = 'HLA SNPs file', required = TRUE)
parser$add_argument('-o_f', '--mr_result', help = 'MR result file', required = TRUE)
parser$add_argument('-plot_f', '--mr_plots', help = 'MR plots file', required = TRUE)
parser$add_argument('-exp_par', '--expo_trait', help = 'Exposure trait', required = TRUE) 
parser$add_argument('-out_par', '--outcome_trait', help = 'Outcome trait', required = TRUE)

# Parse arguments
args <- parser$parse_args()

arg_mr_data <- args$mr_data
arg_HLA_snps <- args$HLA_snps
arg_mr_result <- args$mr_result
arg_mr_plots <- args$mr_plots
arg_expo_trait <- args$expo_trait
arg_outcome_trait <- args$outcome_trait

################################################
### Rsids to exclude files #####################
HLA_rsids <- read.table(arg_HLA_snps, header = F)
names(HLA_rsids) = c('SNP')

################################################
### TwoSampleMR ################################
EXPOSURE <- arg_expo_trait
OUTCOME <- arg_outcome_trait
    
### Read exposure data
exposure_dat <- read_exposure_data(
    file.path(arg_mr_data),
                                   sep = " ",
                                   snp_col = "SNP", 
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2",
                                   eaf_col = "Freq_EXPO",
                                   beta_col = "b_EXPO",
                                   se_col = "se_EXPO",
                                   pval_col = "p_EXPO",
                                   samplesize_col = "N_EXPO")

exposure_dat$exposure <- EXPOSURE

### Read outcome data
outcome_dat <- read_outcome_data(file.path(arg_mr_data),
                                 sep = " ",
                                 snp_col = "SNP", 
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "Freq_OUT",
                                 beta_col = "b_OUT",
                                 se_col = "se_OUT",
                                 pval_col = "p_OUT",
                                 samplesize_col = "N_OUT")

outcome_dat$outcome <- OUTCOME

### Harmonize data
dat <- harmonise_data(exposure_dat, outcome_dat)
# In cases where eaf = 0 NA's are introduced
if(sum(rowSums(is.na(dat)) > 0) != 0){
  print("Removing the following SNPs (eaf = 0)")
  print(dat[rowSums(is.na(dat)) > 0, "SNP"]) 
}

# Here we remove such cases
dat <- dat[rowSums(is.na(dat)) == 0,]

### allele frequency filter
print(paste0("Starting with ", nrow(dat), " SNPs"))
dat <- dat[abs(dat$eaf.exposure - dat$eaf.outcome) < 0.05,]
print(paste0("Number of SNPs after frequency check: ", nrow(dat)))

### Steiger filter: remove SNPs with larger outcome than exposure effects
dat$zval_steiger <- (abs(dat$beta.exposure)-abs(dat$beta.outcome))/sqrt(dat$se.exposure**2 + dat$se.outcome**2)
dat <- dat[dat$zval_steiger > -1.96, ]
print(paste0("Number of SNPs after steiger filter: ", nrow(dat)))

### With TwoSampleMR function
# dat <- steiger_filtering(dat)
# dat <- dat[(dat$steiger_dir) | (!dat$steiger_dir & dat$steiger_pval > 0.05), ]
# print(paste0("Number of SNPs after Steiger filter: ", nrow(dat)))

### Remove SNP from the HBB locus (chr11 5246696-5248301)
# dat <- dplyr::filter(dat, !SNP %in% HBB_rsids$SNP)
# print(paste0("Number of SNPs without HBB locus: ", nrow(dat)))
    
### Remove SNPs from HLA locus (chr6 25MB-37MB)
dat <- dplyr::filter(dat, !SNP %in% HLA_rsids$SNP)
print(paste0("Number of SNPs without HLA locus: ", nrow(dat)))

### MR analysis
mr_results <- mr(dat)

### MR results to export
mr_add_res_table <- mr_results %>% select(c("exposure", "outcome", "method", "nsnp", 
                                            "b", "se", "pval"))
if(nrow(dat) >= 2){ # If only one SNP can't compute heterogeneity
  q_stat <- mr_heterogeneity(dat)[,c("method","Q","Q_pval")]
} else{
  q_stat <- data.frame("method" = c("MR Egger","Inverse variance weighted"), 
                       "Q" = c(NA, NA),
                       "Q_pval" = c(NA,NA))
}
mr_add_res_table <- merge(mr_add_res_table, q_stat, by="method", all.x = T)
mr_add_res_table <- mr_add_res_table[,c("exposure", "outcome", "method", "nsnp", 
                                        "b", "se", "pval","Q","Q_pval")]

################################################
### Plotting ###################################

### Avoid regenerating existing plots
# Scatter plot
p1 <- mr_scatter_plot(mr_results, dat)

# Forest plot
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)

# Leave-one-out plot
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)

# Funnel plot
p4 <- mr_funnel_plot(res_single)

### Plot all plots together
mr_plots <- gridExtra::grid.arrange(grobs = c(p1,p2,p3,p4), ncol=2)

################################################
### Save results ###############################
fwrite(mr_add_res_table, arg_mr_result, col.names = T, row.names = F, quote = F, sep = "\t")
ggsave(arg_mr_plots, mr_plots, width=3, height=3, units="in", scale=3)  
