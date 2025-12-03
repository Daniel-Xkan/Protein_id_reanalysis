import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any
import argparse
from tqdm import tqdm

import urllib.parse
#cd usi_acquisition/usi_filter
#python3 filter_pride_api.py ../new_usi/all_usi.txt
def read_usi_file(filepath: str) -> List[str]:
    """Read USI identifiers from a file."""
    with open(filepath, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def fetch_spectrum(usi: str) -> Dict[str, Any]:
    """Fetch spectrum data from PRIDE API for a given USI."""
    base_url = "https://www.ebi.ac.uk/pride/proxi/archive/v0.1/spectra"
    
    params = {
        'resultType': 'full',
        'usi': usi  # Don't encode - requests will handle it correctly
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        return {'usi': usi, 'data': response.json(), 'error': None}
    except Exception as e:
        return {'usi': usi, 'data': None, 'error': str(e)}


def check_mass_difference(spectrum_data: List[Dict[str, Any]], threshold: float = 4.0) -> bool:
    """Check if the difference between selected ion m/z and isolation window target m/z exceeds threshold."""
    if not spectrum_data or len(spectrum_data) == 0:
        return False
    
    attributes = spectrum_data[0].get('attributes', [])
    
    selected_ion_mz = None
    isolation_window_mz = None
    
    for attr in attributes:
        if attr.get('accession') == 'MS:1000744':  # selected ion m/z
            selected_ion_mz = float(attr.get('value'))
        elif attr.get('accession') == 'MS:1000827':  # isolation window target m/z
            isolation_window_mz = float(attr.get('value'))
    
    if selected_ion_mz is not None and isolation_window_mz is not None:
        difference = abs(selected_ion_mz - isolation_window_mz)
        return difference > threshold
    
    return False


def process_usi_file(input_file: str, output_file: str = 'pride_mass_differences.txt', 
                     threshold: float = 4.0, max_workers: int = 10):
    """Process USI file and filter spectra based on mass difference threshold."""
    usi_list = read_usi_file(input_file)
    results_to_write = []
    
    print(f"Processing {len(usi_list)} USIs with threshold {threshold}...")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_usi = {executor.submit(fetch_spectrum, usi): usi for usi in usi_list}
        
        for future in tqdm(as_completed(future_to_usi), total=len(usi_list), desc="Fetching spectra"):
            result = future.result()
            usi = result['usi']
            
            if result['error']:
                print(f"Error fetching {usi}: {result['error']}")
                with open('error.txt', 'a') as error_file:
                    error_file.write(f"{usi}\t{result['error']}\n")
                continue
            
            if check_mass_difference(result['data'], threshold):
                results_to_write.append(usi)
    
    # Write results to output file
    with open(output_file, 'w') as f:
        for usi in results_to_write:
            f.write(f"{usi}\n")
    
    print(f"Done! Found {len(results_to_write)} USIs exceeding threshold.")
    print(f"Results written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter USIs based on mass difference from PRIDE API')
    parser.add_argument('input_file', help='Path to file containing USI identifiers')
    parser.add_argument('--output', default='pride_mass_differences.txt', help='Output file path')
    parser.add_argument('--threshold', type=float, default=4.0, help='Mass difference threshold')
    parser.add_argument('--workers', type=int, default=10, help='Number of parallel workers')
    
    args = parser.parse_args()
    
    process_usi_file(args.input_file, args.output, args.threshold, args.workers)
