import numpy as np
import pandas as pd

protein = snakemake.wildcards["prot"]
gene_start = snakemake.params["gene_start"]

snp_annotation = pd.read_csv(snakemake.input["snp_info"], compression='gzip', sep = '\t', usecols=["rsid", "POS38"])
pQTL_df = pd.read_csv(snakemake.params["path"], compression='gzip', sep = ' ', 
                      usecols=["GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P" ])

 

print("Gene start")
print(gene_start)

print("SNP annotation")
print(snp_annotation.head(5))

print("pQTL")
print(pQTL_df.head(5))

merged_df = pd.merge(snp_annotation, pQTL_df, left_on='POS38', right_on='GENPOS')

# Drop the 'GENPOS' column as it is redundant
merged_df.drop('GENPOS', axis=1, inplace=True)

# Sort the merged dataframe by POS38 in ascending order
merged_df = merged_df.sort_values(by='POS38')

# Format to .ma format
merged_df['z'] =  merged_df['BETA']/merged_df['SE']
merged_df['b'] = merged_df['z'] / np.sqrt(merged_df['N'] + merged_df['z']**2)
merged_df['se'] = 1/np.sqrt(merged_df['N'] + merged_df['z']**2)

merged_df['p'] = 10**(-merged_df['LOG10P'].astype(np.float64))
merged_df = merged_df.rename(columns = {'ALLELE1':'A1', 'ALLELE0':'A2', 'A1FREQ':'freq', "rsid": "SNP"}) 
final_df = merged_df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N', "POS38"]].drop_duplicates()
final_df = final_df.dropna()

# Save the resulting dataframe
final_df.to_csv(snakemake.output["format_sumstat"], index = False, sep = '\t')

