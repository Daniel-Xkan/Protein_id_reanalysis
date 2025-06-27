import os
import sys
from multiprocessing import Pool, cpu_count

import requests
from bs4 import BeautifulSoup


# Function to create the correct PeptideAtlas URL for a given peptide
def create_peptide_url(peptide_id):
    base_url = "https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide"
    url = f"{base_url}?_tab=3&atlas_build_id=572&searchWithinThis=Peptide+Name&searchForThis={peptide_id}&action=QUERY"
    return url

# Function to capture the USI from the Spectrum page
def get_usi_from_spectrum(spectrum_url):
    try:
        response = requests.get(spectrum_url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')

        usi_input = soup.find('input', {'name': 'USI'})
        if usi_input:
            usi = usi_input['value']
            if usi.startswith('mzspec:PXD') or usi.startswith('mzspec:MSV'):
                return usi
    except requests.exceptions.RequestException:
        pass
    return None

# Function to extract USIs from the individual spectra table
def process_peptide_spectra(data):
    protein_id, peptide_id, worker_id = data
    url = create_peptide_url(peptide_id)

    try:
        response = requests.get(url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')

        spectra_table = soup.find('table', {'id': 'individual_spectra'})
        if not spectra_table:
            return

        for row in spectra_table.find_all('tr')[1:]:  # Skip the header row
            spectrum_link = row.find('a', href=True)
            if spectrum_link:
                spectrum_url = "https://db.systemsbiology.net" + spectrum_link['href']
                usi = get_usi_from_spectrum(spectrum_url)
                if usi:
                    peptide_identifier = usi.split(':')[-1].split('/')[0]
                    modified = '[' in peptide_identifier or ']' in peptide_identifier
                    save_usi(protein_id, usi, modified)
    except requests.exceptions.RequestException:
        pass

    print(f"Worker {worker_id} finished processing peptide {peptide_id}")

# Function to save the USI to the corresponding file
def save_usi(protein_id, usi, modified=False):
    file_name = f"{protein_id}_{'modified' if modified else 'unmodified'}_peptides_usi.txt"
    with open(file_name, 'a') as file:
        file.write(f"{usi}\n")

# Function to process peptide file split into parts
def process_peptide_file(file_path, worker_id):
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            protein_id, peptide_id, _ = line.strip().split('\t')
            process_peptide_spectra((protein_id, peptide_id, worker_id))

# Split the large file into smaller parts
def split_file(file_path, num_splits=50):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    header = lines[0]
    data_lines = lines[1:]

    split_size = len(data_lines) // num_splits
    for i in range(num_splits):
        split_file_name = f"{file_path.split('.txt')[0]}_part_{i + 1}.txt"
        with open(split_file_name, 'w') as split_file:
            split_file.write(header)  # Write the header
            start_index = i * split_size
            end_index = start_index + split_size if i < num_splits - 1 else len(data_lines)
            split_file.writelines(data_lines[start_index:end_index])

# Main function for parallel processing
def main(file_path, num_splits=50):
    split_file(file_path, num_splits)
    files = [f"{file_path.split('.txt')[0]}_part_{i + 1}.txt" for i in range(1, num_splits + 1)]

    # Use a Pool with imap_unordered for multiprocessing
    with Pool(processes=cpu_count()) as pool:
        pool.imap_unordered(process_peptide_file_with_worker, enumerate(files, 1))
        pool.close()
        pool.join()

def process_peptide_file_with_worker(data):
    worker_id, file = data
    process_peptide_file(file, worker_id)

# Usage
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_peptide_file>")
        sys.exit(1)

    peptide_file = sys.argv[1]
    main(peptide_file)
