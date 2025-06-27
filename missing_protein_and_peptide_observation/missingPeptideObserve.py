import pandas as pd
import requests
from bs4 import BeautifulSoup

# Paths to files
peptideatlas_missing_file = 'PeptideAtlas_missing_in_MassIVE-KB.tsv'
massive_kb_peptides_file = 'MassIVE-KB_HPP_peptides.tsv'
output_peptides_file = 'PeptideAtlas_peptides.tsv'
output_observations_file = 'PeptideAtlas_peptide_observations.tsv'

# Load the PeptideAtlas missing protein data and MassIVE-KB peptides
missing_proteins_df = pd.read_csv(peptideatlas_missing_file, sep='\t')

# Load the MassIVE-KB peptide sequences into a set for faster lookup
with open(massive_kb_peptides_file, 'r') as f:
    massive_kb_peptide_sequences = set(line.strip() for line in f)

# Check if a peptide sequence is in MassIVE-KB
def is_peptide_in_massive_kb(sequence):
    return sequence in massive_kb_peptide_sequences

# Function to generate canonical URL for a given gene accession
def get_canonical_link(gene_accession):
    base_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein"
    atlas_build_id = "572"
    return f"{base_url}?atlas_build_id={atlas_build_id}&apply_action=QUERY&protein_name={gene_accession}"

# Function to fetch the distinct peptides for a given protein from the canonical page
def fetch_distinct_peptides(canonical_url, gene_accession):
    try:
        response = requests.get(canonical_url)
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find the peptide table by its ID 'individual_spectra'
        table = soup.find('table', {'id': 'individual_spectra'})
        if not table:
            print(f"No peptide table found for {gene_accession}")
            return
        
        # Process each row in the peptide table
        for row in table.find_all('tr')[1:]:  # Skip the header row
            cols = row.find_all('td')
            if len(cols) < 17:
                continue  # Ensure the row has enough columns
            
            sequence = cols[4].text.strip()  # Sequence column
            is_unique = cols[16].text.strip()  # is_unique column
            n_experiments = cols[13].text.strip()  # N Experiments column

            # Skip if peptide is in MassIVE-KB or is not unique
            if is_peptide_in_massive_kb(sequence) or is_unique == 'N':
                continue
            
            # Extract peptide information
            accession = cols[0].text.strip()  # Peptide Accession
            peptide_info = [gene_accession] + [col.text.strip() for col in cols]  # Save the peptide details

            # Fetch and process observations for the current peptide
            num_experiments_observed = fetch_peptide_observations(accession, sequence, gene_accession)

            # Add the calculated number of experiments (num_experiments_observed) to the peptide information
            peptide_info.append(str(num_experiments_observed))  # Add NumExperiments

            # Save peptide information to output file
            with open(output_peptides_file, 'a') as peptides_file:
                peptides_file.write('\t'.join(peptide_info) + '\n')
    except requests.exceptions.ConnectionError:
      print("Timeout occurred at fetching "+ gene_accession)
# Function to fetch and parse the observations for a given peptide
def fetch_peptide_observations(peptide_accession, sequence, gene_accession):
    observations_url = f"https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis={peptide_accession}&action=QUERY"
    try:
        response = requests.get(observations_url)
        soup = BeautifulSoup(response.content, 'html.parser')

        # Find the observations table
        # Find all tables with the class 'PA_sort_table'
        tables = soup.find_all('table', {'class': 'PA_sort_table'})

        # Check if there are at least two tables found
        if len(tables) >= 2:
            table = tables[1]  # Select the second occurrence (index 1)

        if not table:
            print(f"No observations table found for peptide {peptide_accession}")
            return 0

        # Locate the table
        table2 = soup.find('table', style="margin-left:15px;")
        
        if not table2:
            print(f"No overview table found for peptide {peptide_accession}")
            return 0

        # Locate the table under the div with id 'getpeptide_overview_div'
        table2 = soup.find('div', {'id': 'getpeptide_overview_div'}).find('table')

        # Select the 8th row (index 7) for the '# Experiments' info
        experiments_row = table2.find_all('tr')[7]

        # Extract the number of experiments from the <b> tag inside the corresponding row
        num_experiments = experiments_row.find_all('td')[1].b.text

        # Process each row in the observation table
        for row in table.find_all('tr')[1:]:  # Skip the header row
            cols = row.find_all('td')
            if len(cols) < 7:
                continue  # Ensure the row has enough columns

            experiment_info = [gene_accession, peptide_accession, sequence] + [col.text.strip() for col in cols]
            
            # Save observation details to output file
            with open(output_observations_file, 'a') as observations_file:
                observations_file.write('\t'.join(experiment_info) + '\n')
            
            # num_experiments += 1
        
        return num_experiments  # Return the count of experiments observed

    except requests.exceptions.ConnectionError:
        print("Timeout occurred at fetching "+ peptide_accession)

# Main processing loop for missing proteins
with open(output_peptides_file, 'w') as peptides_file, open(output_observations_file, 'w') as observations_file:
    for index, protein_row in missing_proteins_df.iterrows():
        gene_accession = protein_row['accession']
        
        # Get the canonical URL for the current gene accession
        canonical_url = get_canonical_link(gene_accession)
        
        # Fetch the canonical page and parse the distinct peptides
        fetch_distinct_peptides(canonical_url, gene_accession)

print("Processing complete.")
# hi professor, for peptides that satisfy all requirements and have over 500+ experiments for example: https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis=PAp00764553&action=QUERY
# do I have to put every experiment down to the file? cus the web truncate automatically to 500 columns and it would be much more complicated to download and operate, and other than that I have the program running and have the two output files as follows. If every experiment is required ill take more time to create and deleting temp tsv.