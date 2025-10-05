import pandas as pd
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Harmonize SNPs based on input GWAS and SNP list.")
parser.add_argument('--outcome_gwas', help='Path to outcome GWAS .ma file', required=True)
parser.add_argument('--expo_snps', help='Path to exposure SNPs file', required=True)
parser.add_argument('--harmo_snps', help='Path for output harmonized SNPs file', required=True)

args = parser.parse_args()

# Load the data
df = pd.read_csv(args.outcome_gwas, sep = '\t')
snps_df = pd.read_csv(args.expo_snps, sep='\t', names=['SNP'])

# Harmonize the SNPs
snps_df = snps_df[snps_df.SNP.isin(df.SNP)]

# Save the harmonized SNPs
snps_df.to_csv(args.harmo_snps, index=False, header=False)
