#step one for missing gene accession
import pandas as pd

# File paths
peptideatlas_file = 'query_guest_20241002-152803-981.tsv'
massive_kb_proteins_file = 'MassIVE-KB_HPP_proteins.tsv'
output_file = 'PeptideAtlas_missing_in_MassIVE-KB.tsv'

# Load PeptideAtlas dataset
peptideatlas_df = pd.read_csv(peptideatlas_file, sep='\t')

# Load MassIVE-KB protein Gene Accessions as a single column
with open(massive_kb_proteins_file, 'r') as file:
    massive_kb_protein_accessions = set(line.strip() for line in file)

# Filter PeptideAtlas data for proteins not in MassIVE-KB
missing_proteins_df = peptideatlas_df[~peptideatlas_df['accession'].isin(massive_kb_protein_accessions)]

# Write the missing proteins to a new TSV file
missing_proteins_df.to_csv(output_file, sep='\t', index=False)

print(f"Proteins missing in MassIVE-KB have been written to {output_file}")
