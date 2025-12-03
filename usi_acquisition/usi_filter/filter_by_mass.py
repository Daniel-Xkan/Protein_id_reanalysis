import argparse
import requests
import re
import html
import time
import sys
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from multiprocessing import Pool, cpu_count

#!/usr/bin/env python3
"""
Simple helper to query a USI and print the precursor mz using the ProteoSAFe endpoints.
Save as a script and run: python get_mz.py --usi "mzspec:..." 
"""
import urllib.parse

QUERY_USI_BASE_URL = 'https://massive.ucsd.edu/ProteoSAFe/QuerySpectrum?id='
LORIKEET_BASE_URL = (
    'http://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?invoke='
    'annotatedSpectrumImageText&block=0&format=JSON&peptide=*..*&uploadfile=True'
)
LORIKEET_FALLBACK_BASE_URL = LORIKEET_BASE_URL + '&task=4f2ac74ea114401787a7e96e143bb4a1'
HTML_TAG_RE = re.compile(r'<[^>]+>')

def parse_args():
    p = argparse.ArgumentParser(description='Get precursor mz for a USI')
    p.add_argument('--usi', '-u', help='USI string (e.g. mzspec:...)')
    p.add_argument('--file', '-f', help='File containing USI strings (one per line)')
    p.add_argument('--attempts', type=int, default=3, help='Retries for lorikeet call')
    p.add_argument('--output', '-o', default='mass_differences.txt', help='Output file for mass differences')
    p.add_argument('--tolerance', '-t', type=float, default=4, help='Mass difference tolerance (Da)')
    args = p.parse_args()
    
    if not args.usi and not args.file:
        p.error('Either --usi or --file must be provided')
    if args.usi and args.file:
        p.error('Cannot specify both --usi and --file')
    
    return args

def parse_servlet_error(error_message):
    return html.unescape(re.sub(HTML_TAG_RE, '', error_message.replace('<p>', '\n\n').replace('<h3>', '\n\n')))

def call_lorikeet(filename, nativeid, attempts=3, retry_delay=2):
    url = '{}&file=FILE-%3E{}&{}'.format(
        LORIKEET_BASE_URL, urllib.parse.quote_plus(filename), nativeid
    )
    for attempt in range(attempts):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                j = r.json()
                return j['peaks'], j['precursor']['mz'], j['precursor']['charge']
            elif r.status_code != 200 and attempt == 0:
                # try fallback once
                url = '{}&file=FILE-%3E{}&{}'.format(
                    LORIKEET_FALLBACK_BASE_URL, urllib.parse.quote_plus(filename), nativeid
                )
                time.sleep(retry_delay)
                continue
            elif attempt < attempts - 1:
                time.sleep(retry_delay)
                retry_delay *= 2
                continue
            else:
                raise ValueError('Lorikeet error for URL {}:\n\n{}'.format(url, parse_servlet_error(r.text)))
        except requests.exceptions.RequestException as e:
            if attempt < attempts - 1:
                time.sleep(retry_delay)
                retry_delay *= 2
                continue
            raise ValueError(f'Network error for URL {url}: {e}')
    raise ValueError(f'Failed after {attempts} attempts for URL: {url}')

def get_precursor_mz_for_usi(usi, attempts=3):
    resp = requests.get(QUERY_USI_BASE_URL + urllib.parse.quote_plus(usi), timeout=30)
    if resp.status_code != 200:
        raise ValueError(f'Query error for USI {usi}:\n\n{parse_servlet_error(resp.text)}')
    j = resp.json()
    if not j.get('row_data'):
        raise ValueError(f'No spectrum data found for USI: {usi}')
    spectrum = j['row_data'][0]
    if spectrum.get('resolved_nativeids'):
        nativeid = 'nativeid=' + spectrum['resolved_nativeids']
    else:
        nativeid = spectrum.get('nativeid', '')
    filename = spectrum['file_descriptor'].replace('f.', '') if spectrum.get('file_descriptor') else ''
    _, mz, _ = call_lorikeet(filename, nativeid, attempts=attempts)
    # Derive charge (parse from USI if possible; otherwise query lorikeet once more)
    try:
        charge = int(usi.rsplit('/', 1)[1])
    except Exception:
        _, _, charge = call_lorikeet(filename, nativeid, attempts=attempts)

    mh_actual = mz * charge - (charge - 1) * 1.007276  # monoisotopic MH+
    return mh_actual

def compute_theoretical_mass_from_usi(usi):
    try:
        peptide_part = usi.split(':')[-1]
        peptide = peptide_part.split('/')[0]
        aa_mono = {
            'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185,
            'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064,
            'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764,
            'S':87.032028,'T':101.047679,'W':186.079313,'Y':163.063329,'V':99.068414
        }
        # expanded set of common modification names -> monoisotopic mass shifts
        mod_masses = {
            # alkylation / cysteine
            'Carbamidomethyl': 57.021464,
            'Carboxymethyl': 58.005479,
            'Dicarbamidomethyl': 114.042927,
            'Pyro-carbamidomethyl': 39.994915,

            # oxidation / common variable mods
            'Oxidation': 15.994915,
            'Ox': 15.994915,

            # phosphorylation
            'Phospho': 79.966331,
            'Phosphorylation': 79.966331,

            # acetylation / N-term
            'Acetyl': 42.010565,
            'Acetylation': 42.010565,

            # deamidation
            'Deamidation': 0.984016,
            'Deamidated': 0.984016,
            'Deamidated:18O(1)': 2.988261,

            # ubiquitin remnant
            'GlyGly': 114.042927,
            'GlyGlyRemnant': 114.042927,
            'GG': 114.042927,
            'LRGG': 383.228103,

            # methylation series
            'Methyl': 14.015650,
            'Dimethyl': 28.031300,
            'Dimethyl:2H(4)': 32.056407,
            'Dimethyl:2H(4)13C(2)': 34.063117,
            'Dimethyl:2H(6)13C(2)': 36.075670,
            'Trimethyl': 42.046950,

            # small additions / others
            'Formyl': 27.994915,
            'Sulfation': 79.956815,
            'Pyro-glu': -17.026549,
            'PyroGlu': -17.026549,
            'Gln->pyro-Glu': -17.026549,
            'Glu->pyro-Glu': -18.010565,
            'Propionyl': 56.026215,
            'Propionamide': 71.037114,
            'Carbamyl': 43.005814,
            'Dethiomethyl': -48.003371,
            'Ammonia-loss': -17.026549,
            'Methylthio': 45.987721,
            'Cysteinyl': 119.004099,
            'Nethylmaleimide': 125.047679,
            'Thiazolidine': 87.998285,

            # common labeling reagents (mass added to peptide backbone)
            'TMT6plex': 229.162932,
            'TMT6plex114': 229.162932,
            'TMT10plex': 229.162932,
            'TMT11plex': 229.162932,
            'TMT16plex': 229.162932,
            'TMTpro': 304.207146,

            'iTRAQ4plex': 144.102063,
            'iTRAQ4plex114': 144.102063,
            'iTRAQ8plex': 304.205360,
            'iTRAQ8plex:13C(6)15N(2)': 304.205360,

            'DiLeu4plex117': 145.125595,

            # isotope labels
            'Label:13C(6)': 6.020129,
            'Label:13C(6)15N(1)': 7.017164,
            'Label:13C(6)15N(2)': 8.014199,
            'Label:13C(6)15N(4)': 10.008269,
            'Label:2H(4)': 4.025107,

            # ADP-Ribosyl
            'ADP-Ribosyl': 541.061110,

            # numeric modifications (use exact values from your data)
            '+141.11544': 141.11544,
            '+141.1154': 141.1154,
            '+1431.83104': 1431.83104,
            '+1541.85014': 1541.85014,
            '+1555.95614': 1555.95614,
            '+186.1127': 186.1127,
            '+186.1165': 186.1165,
            '+229.162931': 229.162931,
            '+271.1735': 271.1735,
            '+271.1736': 271.1736,
            '+454.18121': 454.18121,
            '+46.03274': 46.03274,
            '+46.0328': 46.0328,
            '+471.20776': 471.20776,
            '+64.10696': 64.10696,
            '+75.04729': 75.04729,
            '+85.05549': 85.05549,
            '0.0233': 0.0233,
            'Xlink:BuUrBu[85]': 85.0,
        }

        total = 18.010564  # + H2O
        i = 0
        # handle leading N-terminal modification like [Dimethyl]-...
        n = len(peptide)
        if i < n and peptide[i] == '[':
            j = peptide.find(']', i)
            if j != -1:
                mod_name = peptide[i+1:j]
                total += mod_masses.get(mod_name, 0.0)
                i = j + 1
                if i < n and peptide[i] == '-':
                    i += 1

        while i < n:
            ch = peptide[i]
            if ch in aa_mono:
                total += aa_mono[ch]
                i += 1
                # residue-level modification like K[Dimethyl]
                if i < n and peptide[i] == '[':
                    j = peptide.find(']', i)
                    if j != -1:
                        mod_name = peptide[i+1:j]
                        total += mod_masses.get(mod_name, 0.0)
                        i = j + 1
                continue
            # skip common separators
            if ch in '-.()':
                i += 1
                continue
            # stray bracket (treat as N-term or apply to next residue)
            if ch == '[':
                j = peptide.find(']', i)
                if j != -1:
                    mod_name = peptide[i+1:j]
                    total += mod_masses.get(mod_name, 0.0)
                    i = j + 1
                    if i < n and peptide[i] == '-':
                        i += 1
                    continue
            # unknown character - skip
            i += 1
        return total + 1.007276
    except Exception:
        return None

def main():
    args = parse_args()
    if args.file:
        with open(args.file, 'r') as f:
            usi_list = [line.strip() for line in f if line.strip()]
        
        results = []
        
        
        results_lock = Lock()
        
        def process_usi(usi):
            try:
                mh_actual = get_precursor_mz_for_usi(usi, attempts=args.attempts)
                theoretical_mass = compute_theoretical_mass_from_usi(usi)
                
                if theoretical_mass is not None:
                    mass_diff = abs(mh_actual - theoretical_mass)
                    if mass_diff > args.tolerance:
                        return {
                            'usi': usi,
                            'mh_actual': mh_actual,
                            'theoretical_mass': theoretical_mass,
                            'mass_diff': mass_diff
                        }
            except Exception as e:
                tqdm.write(f'Error processing {usi}: {e}')
            return None
        
        # Use threading instead of multiprocessing
        num_threads = min(10, len(usi_list))
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = {executor.submit(process_usi, usi): usi for usi in usi_list}
            
            for future in tqdm(as_completed(futures), total=len(usi_list), desc="Processing USIs", unit="USI"):
                result = future.result()
                if result is not None:
                    with results_lock:
                        results.append(result)
                        # Write immediately to file
                        with open(args.output, 'a') as out:
                            # Write header if this is the first result
                            if len(results) == 1:
                                out.write('USI\tMH_Actual\tTheoretical_Mass\tMass_Difference\n')
                            out.write(f"{result['usi']}\t{result['mh_actual']:.6f}\t{result['theoretical_mass']:.6f}\t{result['mass_diff']:.6f}\n")
        
        print(f"\nSummary written to {args.output}")
        print(f"Total USIs processed: {len(usi_list)}")
        print(f"USIs with mass difference > {args.tolerance} Da: {len(results)}")
        sys.exit(0)
    else: 
        try:
            mh_actual = get_precursor_mz_for_usi(args.usi, attempts=args.attempts)
            print(mh_actual)
        except Exception as e:
            print(f'Error: {e}', file=sys.stderr)
            if 'mz' not in locals():
                sys.exit(1)

        theoretical_mass = compute_theoretical_mass_from_usi(args.usi)
        if theoretical_mass is not None:
            print('theoretical_mass', theoretical_mass)

if __name__ == '__main__':
    main()