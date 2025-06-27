import requests
from bs4 import BeautifulSoup


# Function to create the correct PeptideAtlas URL for a given peptide
def create_peptide_url(peptide_id):
    base_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide"
    url = f"{base_url}?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis={peptide_id}&action=QUERY"
    print(f"Created PeptideAtlas URL: {url}")
    return url

# Function to capture the USI from the Spectrum page
def get_usi_from_spectrum(spectrum_url):
    try:
        print(f"Fetching Spectrum URL: {spectrum_url}")
        response = requests.get(spectrum_url)
        response.raise_for_status()  # Raises an error for bad responses
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Locate the USI form and extract the value
        usi_input = soup.find('input', {'name': 'USI'})
        if usi_input:
            usi = usi_input['value']
            print(f"USI found: {usi}")
            if usi.startswith('mzspec:PXD') or usi.startswith('mzspec:MSV'):
                return usi
        else:
            print("No USI found in this spectrum.")
    except requests.exceptions.RequestException as e:
        print(f"Error fetching spectrum URL: {spectrum_url}, Error: {e}")
    
    return None

# Function to extract USIs from the individual spectra table
def process_peptide_spectra(protein_id, peptide_id):
    print(f"Processing peptide: {peptide_id} for protein: {protein_id}")
    url = create_peptide_url(peptide_id)
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raises an error for bad responses
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Parse individual spectra table
        spectra_table = soup.find('table', {'id': 'individual_spectra'})
        if not spectra_table:
            print(f"No spectra table found for peptide: {peptide_id}")
            return
        
        for row in spectra_table.find_all('tr')[1:]:  # Skip the header row
            spectrum_link = row.find('a', href=True)
            if spectrum_link:
                spectrum_url = "https://db.systemsbiology.net" + spectrum_link['href']
                usi = get_usi_from_spectrum(spectrum_url)
                if usi:
                    peptide_identifier = usi.split(':')[-1].split('/')[0]
                    if '[' in peptide_identifier or ']' in peptide_identifier:
                        print(f"Modified peptide identified: {peptide_identifier}")
                        save_usi(protein_id, usi, modified=True)
                    else:
                        print(f"Unmodified peptide identified: {peptide_identifier}")
                        save_usi(protein_id, usi, modified=False)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching peptide URL: {url}, Error: {e}")

# Function to save the USI to the corresponding file
def save_usi(protein_id, usi, modified=False):
    file_name = f"{protein_id}_{'modified' if modified else 'unmodified'}_peptides_usi.txt"
    print(f"Saving USI to file: {file_name}")
    
    with open(file_name, 'a') as file:
        file.write(f"{usi}\n")

# Process the PA2024_HPP_peptides.txt file
def process_peptide_file(file_path):
    print(f"Processing peptide file: {file_path}")
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            protein_id, peptide_id, sequence = line.strip().split('\t')
            print(f"Processing line: {line.strip()}")
            process_peptide_spectra(protein_id, peptide_id)

# Example usage
file_path = 'PA2024_HPP_peptides.txt'
process_peptide_file(file_path)

print("Processing complete and USIs stored.")
