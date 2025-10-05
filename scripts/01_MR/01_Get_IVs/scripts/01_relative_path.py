############################################################################################################
### Script to add the relative path column to the protein_annotation_3k.tsv file of the UKB PPP ############¨
############################################################################################################

### Example
# PCSK9_Q8NBP7_OID20235_v1_Cardiometabolic/discovery_chrX_PCSK9:Q8NBP7:OID20235:v1:Cardiometabolic.gz
# UKBPPP_ProteinID would be PCSK9:Q8NBP7:OID20235:v1 and Panel would be Cardiometabolic

# Import the pandas library
import pandas as pd

# Define the path to the input file
path_to_protein_annotation = "/.../protein_annotation_3k.tsv"

# Read the input file
df = pd.read_csv(path_to_protein_annotation, sep='\t')

# Function to generate the rel_path based on the given columns
def generate_rel_path(row):
    protein_id_formatted = row['UKBPPP_ProteinID'].replace(':', '_')
    rel_path = f"{protein_id_formatted}_{row['Panel']}/discovery_chr{row['chr']}_{row['UKBPPP_ProteinID']}:{row['Panel']}.gz"
    return rel_path

# Apply the function to generate the rel_path
df['rel_path'] = df.apply(generate_rel_path, axis=1)

# Save the resulting dataframe to a new tsv file
df.to_csv('/.../protein_annotation_3k_RP.tsv', sep='\t', index=False)
