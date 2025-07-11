import argparse
from pathlib import Path
import requests
import urllib.parse
import re
import html
import os
import time
import sys
##################tier seperation##################
def read_tier_assignments(file_path):
    """
    Reads USI tier assignments from a TSV file.
    
    Tier 1: MSV column is not empty and Importing != 'import'
    Tier 2: MSV is empty but Importing = 'import'
    Tier 3: All others
    
    Returns a dictionary with {usi: tier_number}
    """
    usi_tiers = {}
    print(f"Reading tier assignments from {file_path}")
    try:
        with open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            
            # Find the column indices
            usi_idx = header.index('USI') if 'USI' in header else None
            msv_idx = header.index('MSV') if 'MSV' in header else None
            importing_idx = header.index('Importing') if 'Importing' in header else None
            
            if any(idx is None for idx in [usi_idx, msv_idx, importing_idx]):
                print(f"Warning: Required columns not found in {file_path}")
                print(f"Expected: USI, MSV, Importing")
                print(f"Found: {header}")
                return usi_tiers
            
            for line in f:
                if not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) <= max(usi_idx, msv_idx, importing_idx):
                    continue
                
                usi = fields[usi_idx]
                msv = fields[msv_idx].strip()
                importing = fields[importing_idx].strip()
                
                # Assign tier based on the criteria
                if 'MSV' in msv and importing != 'Import':
                    tier = 1
                elif 'MSV' not in msv and importing == 'Import':
                    tier = 2
                else:
                    tier = 3
                
                usi_tiers[usi] = tier
                
        print(f"Loaded {len(usi_tiers)} USI tier assignments")
        
    except Exception as e:
        print(f"Error reading tier assignments: {str(e)}")
    
    return usi_tiers
###################################################
QUERY_USI_BASE_URL = 'https://proteomics3.ucsd.edu/ProteoSAFe/QuerySpectrum?id='
LORIKEET_BASE_URL = (
    'http://proteomics3.ucsd.edu/ProteoSAFe/DownloadResultFile?invoke='
    'annotatedSpectrumImageText&block=0&format=JSON&peptide=*..*&uploadfile=True'
)

# Fallback URL with task parameter if the regular URL fails
LORIKEET_FALLBACK_BASE_URL = LORIKEET_BASE_URL + '&task=4f2ac74ea114401787a7e96e143bb4a1'
HTML_TAG_RE = re.compile(r'<[^>]+>')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert a list of USIs to an MGF with all available fields.'
    )
    parser.add_argument(
        '--usi_list', '-i', type=Path, required=True,
        help='Input text file with one USI per line'
    )
    parser.add_argument(
        '--out', '-o', type=Path, required=True,
        help='Output MGF path'
    )
    parser.add_argument(
        '--attempts', type=int, default=3,
        help='Maximum number of attempts (default: 3)'
    )
    parser.add_argument(
        '--tier', '-t', type=int, choices=[1, 2, 3],
        help='Process only USIs belonging to the specified tier (1, 2, or 3)'
    )
    parser.add_argument(
        '--tier_file', type=Path, default="all_usi.tsv",
        help='Path to TSV file containing tier assignments'
    )
    return parser.parse_args()

def main():
    args = parse_args()
    failed_usis = []
    successful_count = 0
    
    # Load tier assignments if tier filtering is requested
    usi_tiers = {}
    if args.tier is not None:
        if not args.tier_file:
            print("Error: --tier_file must be provided when using --tier")
            sys.exit(1)
        usi_tiers = read_tier_assignments(args.tier_file)
        if not usi_tiers:
            print("Warning: No tier assignments loaded, processing all USIs")
    
    with open(args.out, 'w') as out_f:
        with open(args.usi_list) as f:
            lines = [line.strip() for line in f if line.strip()]
            total_usis = len(lines)
            processed_usis = 0
            skipped_usis = 0
            
            for i, usi in enumerate(lines, 1):
                # Skip USIs that don't match the specified tier
                if args.tier is not None:
                    if usi not in usi_tiers:
                        print(f"Skipping USI {i} (not found in tier assignments): {usi[:60]}{'...' if len(usi) > 60 else ''}")
                        skipped_usis += 1
                        continue
                    if usi_tiers[usi] != args.tier:
                        print(f"Skipping USI {i} (tier {usi_tiers[usi]}, requested tier {args.tier}): {usi[:60]}{'...' if len(usi) > 60 else ''}")
                        skipped_usis += 1
                        continue
                
                processed_usis += 1
                try:
                    print(f"Processing USI {i}/{total_usis}: {usi[:60]}{'...' if len(usi) > 60 else ''}")
                    process_usi(usi, out_f, processed_usis, args.attempts)
                    successful_count += 1
                    print(f"✓ Successfully processed USI {i}")
                except Exception as e:
                    print(f"✗ Failed to process USI {i}: {str(e)}")
                    error_msg = str(e)
                    failed_usis.append((i, usi, error_msg))
                    # Write error to file immediately
                    error_file = f"{os.path.splitext(args.usi_list)[0]}_tier{args.tier if args.tier is not None else '_alltier'}_error.txt"
                    with open(error_file, 'a') as err_f:
                        if os.path.getsize(error_file) == 0:  # Add header if file is empty
                            err_f.write(f"=== Failed USIs for {args.usi_list} ===\n\n")
                        err_f.write(f"USI {i}: {usi}\n")
                        err_f.write(f"Error: {error_msg}\n\n")
                    continue
    
    # Print summary
    print(f"\n=== Processing Summary ===")
    print(f"Total USIs: {total_usis}")
    if args.tier is not None:
        print(f"Skipped (not in tier {args.tier}): {skipped_usis}")
        print(f"Attempted processing: {processed_usis}")
    print(f"Successful: {successful_count}")
    print(f"Failed: {len(failed_usis)}")
    
    # Update the summary in the error file
    error_file = f"{os.path.splitext(args.usi_list)[0]}_tier{args.tier if args.tier is not None else '_alltier'}_error.txt"
    # if failed_usis:
        # Read existing content
    error_content = ""
    if os.path.exists(error_file):
        with open(error_file, 'r') as err_f:
            error_content = err_f.read()
    
    # Write updated file with summary at the beginning
    with open(error_file, 'w') as err_f:
        err_f.write(f"=== Processing Summary ===\n")
        err_f.write(f"Total USIs: {total_usis}\n")
        if args.tier is not None:
            err_f.write(f"Skipped (not in tier {args.tier}): {skipped_usis}\n")
            err_f.write(f"Attempted processing: {processed_usis}\n")
        err_f.write(f"Successful: {successful_count}\n")
        err_f.write(f"Failed: {len(failed_usis)}\n\n")
        err_f.write(f"=== Failed USIs for {args.usi_list} ===\n\n")
        err_f.write(error_content.replace(f"=== Failed USIs for {args.usi_list} ===\n\n", ""))
    
    if failed_usis:
        print(f"\n=== Failed USIs ===")
        for idx, usi, error in failed_usis:
            print(f"USI {idx}: {usi}")
            print(f"  Error: {error[:100]}{'...' if len(error) > 100 else ''}")


def process_usi(usi, out_f, i, attempts=3):
    r = requests.get(QUERY_USI_BASE_URL + urllib.parse.quote_plus(usi), timeout=30)
    if r.status_code != 200:
        raise ValueError('Error querying USI {}:\n\n{}'.format(
            usi, parse_servlet_error(r.text)
        ))
    
    json_data = r.json()
    if not json_data.get('row_data') or len(json_data['row_data']) == 0:
        raise ValueError(f'No spectrum data found for USI: {usi}')
        
    spectrum = json_data['row_data'][0]
    if spectrum['resolved_nativeids'] != '':
        nativeid = 'nativeid=' + spectrum['resolved_nativeids']
    else:
        nativeid = spectrum['nativeid']
    (
        spectrum['peaks'], spectrum['precursor_mz'],
        spectrum['precursor_charge']
    ) = call_lorikeet(
        spectrum['file_descriptor'].replace('f.',''), nativeid, attempts
    )
    write_spectrum_to_mgf(usi, spectrum, out_f, i)


def call_lorikeet(filename, nativeid, attempts=3, retry_delay=2):
    url = '{}&file=FILE-%3E{}&{}'.format(
        LORIKEET_BASE_URL, urllib.parse.quote_plus(filename), nativeid
    )
    
    for attempt in range(attempts):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                json_data = r.json()
                return(
                    json_data['peaks'], 
                    json_data['precursor']['mz'],
                    json_data['precursor']['charge']
                )
            elif r.status_code != 200 and attempt == 0:
                # Try fallback URL on first failure before second attempt
                print(f"  Using fallback URL for retry...")
                url = '{}&file=FILE-%3E{}&{}'.format(
                    LORIKEET_FALLBACK_BASE_URL, urllib.parse.quote_plus(filename), nativeid
                )
                time.sleep(retry_delay)
                continue
            elif r.status_code !=200 and attempt < attempts - 1:
                print(f"  Server error (attempt {attempt + 1}/{attempts}), retrying in {retry_delay}s...")
                time.sleep(retry_delay)
                retry_delay *= 2  # Exponential backoff
                continue
            else:
                raise ValueError('Lorikeet error for URL {}:\n\n{}'.format(
                    url, parse_servlet_error(r.text)
                ))
        except requests.exceptions.RequestException as e:
            if attempt < attempts - 1:
                print(f"  Network error (attempt {attempt + 1}/{attempts}), retrying in {retry_delay}s...")
                time.sleep(retry_delay)
                retry_delay *= 2
                continue
            else:
                raise ValueError(f'Network error for URL {url}: {str(e)}')
    
    raise ValueError(f'Failed to retrieve data after {attempts} attempts for URL: {url}')


def parse_servlet_error(error_message):
    return html.unescape(re.sub(
        HTML_TAG_RE, '',
        error_message.replace('<p>', '\n\n').replace('<h3>', '\n\n')
    ))


def write_spectrum_to_mgf(usi, spectrum, out_f, i):
    out_f.write(
        'BEGIN IONS\n'
        'PEPMASS={pepmass:.5f}\n'
        'CHARGE={charge}\n'
        'CHARGE_ID={charge_id}\n'
        'MSLEVEL=2\n'
        'COLLISION_ENERGY=0.0\n'
        'FILENAME=\n'
        'SEQ={seq}\n'
        'PROTEIN=\n'
        'SCANS={spec_index}\n'
        'SCAN={spec_index}\n'
        'PROVENANCE_FILENAME={filename}\n'
        'PROVENANCE_SCAN={scan}\n'
        'PROVENANCE_USI={prov_usi}\n'
        '{peaks}\n'
        'END IONS\n\n'.format(
            pepmass=spectrum['precursor_mz'],
            charge=spectrum['precursor_charge'],
            charge_id=spectrum['charge'],
            seq=spectrum['peptide'],
            spec_index=i,
            filename=spectrum['file_descriptor'][2:],
            scan=spectrum['nativeid'].replace('scan=', ''),
            prov_usi=usi,
            peaks='\n'.join(['{mz:.6f} {intensity:.6f}'.format(
                mz=z[0], intensity=z[1]
            ) for z in spectrum['peaks']])
        )
    )

if __name__ == '__main__':
    main()