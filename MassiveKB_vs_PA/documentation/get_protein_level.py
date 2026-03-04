import pandas as pd
import os
import re

def main(mismatch_tsv_path, matched_tsv_path, mgf_path):
    usage = (
        "Usage:\n"
        "  python get_protein_level.py <mismatched.tsv> <matched.tsv> <spectra.mgf>\n\n"
        "Or call main(mismatch_tsv_path, matched_tsv_path, mgf_path) from Python.\n"
        "Arguments:\n"
        "  mismatch_tsv_path   Path to the mismatched TSV file\n"
        "  matched_tsv_path    Path to the matched TSV file\n"
        "  mgf_path            Path to the MGF file containing spectra and provenance\n"
    )

    # Print usage/help if any argument requests it
    if any(str(x).lower() in ("-h", "--help", "help", "usage") for x in (mismatch_tsv_path, matched_tsv_path, mgf_path)):
        print(usage)
        return

    # Basic validation of inputs
    if not mismatch_tsv_path or not matched_tsv_path or not mgf_path:
        raise ValueError("Missing required arguments.\n\n" + usage)

    # Normalize and check file existence
    mismatch_tsv_path = os.path.abspath(mismatch_tsv_path)
    matched_tsv_path = os.path.abspath(matched_tsv_path)
    mgf_path = os.path.abspath(mgf_path)

    for path, desc in [(mismatch_tsv_path, "mismatched TSV"), (matched_tsv_path, "matched TSV"), (mgf_path, "MGF")]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{desc} file not found: {path}")
# This section reads the matched and mismatched TSV files, extracts relevant information, and creates summary TSV files for both matched and mismatched scans. 
# The summary files include scan information, scores, protein IDs, peptide sequences (with modifications removed), and reasons for mismatch if applicable.

# mismatched and matched files have the same columns, so we can read only the necessary columns to save memory
    columns_to_read = ['rowid', 'Scan', 'Annotation','Score','ProtsAll', 'AnnotationOther', 'ChargeOther','ScoreOther','ProtsAllOther',
                    'MinNTermAdd', 'minNTermSubtract', 'MinCTermAdd', 'minCTermSubtract']
    # mismatch_df = pd.read_csv('PRISMAL_REPROCESS_PA_UniproKB_Normal123_mismatched.tsv', sep='\t', usecols=columns_to_read)
    mismatch_df = pd.read_csv(mismatch_tsv_path, sep='\t', usecols=columns_to_read)
    # matched_df = pd.read_csv('PRISMAL_REPROCESS_PA_UniproKB_Normal123_matched.tsv', sep='\t', usecols=columns_to_read)
    matched_df = pd.read_csv(matched_tsv_path, sep='\t', usecols=columns_to_read)


    def remove_modifications(peptide_string):
        """Remove modification annotations from peptide sequence"""
        if pd.isna(peptide_string):
            return ''
        # Remove modification patterns like +42.011, -17.027, +28.011, etc.
        cleaned = re.sub(r'[+-][\d.]+', '', str(peptide_string))
        return cleaned

    # Create the new dataframe with required columns
    def extract_protein_ids(prots_string):
        """Extract protein IDs from the ProtsAll format"""
        if pd.isna(prots_string):
            return ''
        # Extract protein IDs from format like (sp|Q9H0U3-2|MAGT1_HUMAN,1,13,0,0)
        protein_ids = re.findall(r'\|([A-Z0-9-]+)\|', str(prots_string))
        return ';'.join(protein_ids) if protein_ids else ''

    def extract_fist_protein_id(prots_string):
        """Extract the first protein ID from the ProtsAll format"""
        if pd.isna(prots_string):
            return ''
        match = re.search(r'\|([A-Z0-9-]+)\|', str(prots_string))
        if match:
            return match.group(1).split('-')[0]
        return ''

    matched_summary_df = pd.DataFrame({
        'scan': matched_df['Scan'],
        'DB_search_score': matched_df['Score'],
        'precursor_score': matched_df['ScoreOther'],
        'DB_proteins': matched_df['ProtsAll'].apply(lambda x: extract_fist_protein_id(x) if 'Annotation' in matched_df.columns else ''),
        'precursor_proteins': matched_df['ProtsAllOther'].apply(lambda x: extract_fist_protein_id(x) if 'AnnotationOther' in matched_df.columns else ''),
        'precursor_better': (matched_df['Score'] < matched_df['ScoreOther']).astype(int),
        'charge': matched_df['ChargeOther'],
        'peptide': matched_df['AnnotationOther'],
        'peptide_demod': matched_df['AnnotationOther'].apply(remove_modifications),
        'peptide_length': matched_df['AnnotationOther'].apply(remove_modifications).apply(len),
        'reason_mismatch': '',
        'MinNTermAdd': matched_df['MinNTermAdd'],
        'minNTermSubtract': matched_df['minNTermSubtract'],
        'MinCTermAdd': matched_df['MinCTermAdd'],
        'minCTermSubtract': matched_df['minCTermSubtract']
    })

    # Save to TSV file
    matched_summary_df.to_csv('matched_summary.tsv', sep='\t', index=False)

    print(f"Created match_summary.tsv with {len(matched_summary_df)} rows")


    mismatched_summary_df = pd.DataFrame({
        'scan': mismatch_df['Scan'],
        'DB_search_score': mismatch_df['Score'],
        'precursor_score': mismatch_df['ScoreOther'],
        'DB_proteins': mismatch_df['ProtsAll'].apply(lambda x: extract_fist_protein_id(x) if 'Annotation' in mismatch_df.columns else ''),
        'precursor_proteins': mismatch_df['ProtsAllOther'].apply(lambda x: extract_fist_protein_id(x) if 'AnnotationOther' in mismatch_df.columns else ''),
        'precursor_better': (mismatch_df['Score'] < mismatch_df['ScoreOther']).astype(int),
        'charge': mismatch_df['ChargeOther'],
        'peptide': mismatch_df['AnnotationOther'],
        'peptide_demod': mismatch_df['AnnotationOther'].apply(remove_modifications),
        'peptide_length': mismatch_df['AnnotationOther'].apply(remove_modifications).apply(len),
        'reason_mismatch': '',
        'MinNTermAdd': mismatch_df['MinNTermAdd'],
        'minNTermSubtract': mismatch_df['minNTermSubtract'],
        'MinCTermAdd': mismatch_df['MinCTermAdd'],
        'minCTermSubtract': mismatch_df['minCTermSubtract']
    })



    # Extract peptide information from matched_df
    print("Processing matched_df to create identified_peptide_list...")

    identified_peptide_list = set()

    for _, row in matched_df.iterrows():
        # Get demodified peptide
        peptide_demod = remove_modifications(row.get('Annotation'))
        
        # Extract mass shift from the modified peptide
        mod_matches = re.findall(r'([+-][\d.]+)', str(row.get('Annotation')))
        mass_shift = mod_matches[0] if mod_matches else '+0'
        
        # Get charge (assuming there's a 'Charge' column in matched_df)
        charge = row.get('Charge', 0)
        
        # Add to set as tuple
        identified_peptide_list.add((peptide_demod, mass_shift, charge))

    def analyze_peptide_mismatches(df):
        """
        Function to analyze peptide mismatches and assign peptide types
        """
        def is_identified(peptide_demod, mass_shift, charge, identified_peptide_list):
            """Check if the peptide is in the identified peptide list"""
            return (peptide_demod, mass_shift, charge) in identified_peptide_list

        def count_missed_cleavages(peptide_seq):
            """Count missed cleavages for trypsin (K/R not at C-terminus)"""
            if pd.isna(peptide_seq) or len(peptide_seq) == 0:
                return 0
            missed = 0
            for aa in peptide_seq[:-1]:
                if aa in ['K', 'R']:
                    missed += 1
            return missed

        def is_cterm_tryptic(peptide_seq):
            """Check if C-terminus is tryptic (ends with K or R)"""
            if pd.isna(peptide_seq) or len(peptide_seq) == 0:
                return False
            return peptide_seq[-1] in ['K', 'R']

        def check_HLA(length):
            """Assign peptide type based on length"""
            if 8 <= length <= 12:
                return 'HLA1'
            elif 13 <= length <= 24:
                return 'HLA2'
            else:
                return 'Others'

        def check_C57(peptide_seq):
            """Check if at least one Cysteine carries a +57 modification"""
            if pd.isna(peptide_seq) or len(peptide_seq) == 0:
                return False
            if 'C' not in peptide_seq:
                return True
            return bool(re.search(r'C\+57(?:\.\d+)?', peptide_seq))

        def get_mismatch_reason(row):
            is_ident = is_identified(
                row['peptide_demod'],
                row.get('mass_shift', '+0'),
                row.get('Charge', 0),
                identified_peptide_list
            )
            if is_ident:
                return 'IDENTIFIED'

            standard_mods = [1, 16, 42, 43, -17]
            tolerance = 0.5
            mod_matches = re.findall(r'([+-][\d.]+)', str(row['peptide']))
            for mod_str in mod_matches:
                mod_value = float(mod_str)
                if not any(abs(mod_value - std_mod) <= tolerance for std_mod in standard_mods):
                    return 'non_standard_modification'
                
            if not check_C57(row['peptide']):
                return 'C_without_+57_modification' 
            
            mod_count = row['peptide'].count('+') + row['peptide'].count('-')
            if mod_count >= 2:
                return 'multiple_modifications'

            if row['peptide_length'] > 40:
                return 'peptide_length_>_40'

            if is_cterm_tryptic(row['peptide_demod']):
                if row['MinNTermAdd'] >= 3:
                    return 'lost_3_aa_Nterm'
                if count_missed_cleavages(row['peptide_demod']) >= 3:
                    return 'miss_3_cleavages'
                return 'regular_tryptic'

            if 8 <= row['peptide_length'] <= 12:
                return 'HLA1'
            elif 13 <= row['peptide_length'] <= 24:
                return 'HLA2'


            return 'Others'

        if df.empty:
            df['reason_mismatch'] = pd.Series(dtype='object')
            return df

        reasons = df.apply(get_mismatch_reason, axis=1)
        if isinstance(reasons, pd.DataFrame):
            if reasons.shape[1] == 1:
                reasons = reasons.iloc[:, 0]
            else:
                raise ValueError("Mismatch reason assignment produced multiple columns")
        df['reason_mismatch'] = reasons
        return df

    # Apply the analysis to the summary dataframe
    mismatched_summary_df = analyze_peptide_mismatches(mismatched_summary_df)


    print(f"Created mismatch_summary.tsv with peptide types and mismatch reasons")
    print("\nMismatch reason distribution:")
    print(mismatched_summary_df['reason_mismatch'].value_counts())
    mismatched_summary_df.to_csv('mismatch_summary.tsv', sep='\t', index=False)
    print(f"Created mismatch_summary.tsv with {len(mismatched_summary_df)} rows")

    # At this point, we have two TSV files: 'matched_summary.tsv' and 'mismatch_summary.tsv', each containing detailed information about the matched and mismatched scans, respectively.
    # This section, we will read these summary files and perform protein-level analysis by aggregating the peptide-level information to the protein level. We will create a new TSV file
    # that summarizes the evidence for each protein based on the peptides identified in both matched and mismatched scans.


    #mismatched df
    df = pd.read_csv('mismatch_summary.tsv', sep='\t')
    #matched df
    df2 = pd.read_csv('matched_summary.tsv', sep='\t')

    scan_to_dataset = {}

    with open(mgf_path, 'r') as file:
        dataset = None
        scan_number = None
        for line in file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                dataset = None
                scan_number = None
            elif line.startswith("PROVENANCE_FILENAME="):
                provenance = line.split('=')[1]
                if provenance.startswith('ProteomeCentral'):
                    parts = provenance.split('/')
                    dataset = parts[1] if len(parts) > 1 else None
                else:
                    dataset = provenance.split('/')[0]
            elif line.startswith("SCAN="):
                scan_number = int(line.split('=')[1])
            elif line.startswith("END IONS") and scan_number is not None and dataset is not None:
                scan_to_dataset[scan_number] = dataset

    # print(scan_to_dataset)

    prec_protein = set(df['precursor_proteins'].dropna().str.split(r'[;-]').str[0])
    # print(f"Number of unique precursor proteins: {len(prec_protein)}")
    print(f"Number of unique precursor proteins: {len(prec_protein)}")

    db_protein = set(df2['precursor_proteins'].dropna().str.split(r'[;-]').str[0])
    print(f"Number of unique DB proteins: {len(db_protein)}")



    # Create a dictionary to store the data for each precursor protein
    #iterate through filtered df each row, make a new df with columns: precursor_protein, precursor_id_number - this is the number of rows(peptides) that has that protein in precursor_proteins, [(evidence_peptides,scan_number) - this is a tuple with multiple peptides and scan ]

    protein_data = {}

    for idx, row in df.iterrows():
        # Create a dictionary to store the data for each precursor protein
        #iterate through filtered df each row, make a new df with columns: precursor_protein, precursor_id_number - this is the number of rows(peptides) that has that protein in precursor_proteins, [(evidence_peptides,scan_number) - this is a tuple with multiple peptides and scan ]

        proteins = row['precursor_proteins']
        if pd.isna(proteins):
            continue
        protein = re.split(r'[;-]', proteins)[0]
        peptide = row['peptide_demod']
        if len(peptide) <9:
            continue
        scan = row['scan']
        dataset = scan_to_dataset.get(scan, 'Unknown')
        reason_mismatch = row['reason_mismatch']
        evidence = False
        source = 'Missed'

        
        if protein not in protein_data:
            protein_data[protein] = {
                'precursor_protein': protein,
                'precursor_id_number': 0,
                'evidence_peptides_scans': []
            }
        
        protein_data[protein]['precursor_id_number'] += 1
        protein_data[protein]['evidence_peptides_scans'].append((peptide,dataset, scan,reason_mismatch,evidence,source))

    for protein, data in protein_data.items():
        # print(data)
        peptides_temp = []
        peptides_temp_atomic = []
        peptides_atomic_max_l_dict = {}
        for peptide, dataset, _, _, _,_ in protein_data[protein]['evidence_peptides_scans']:
            peptides_temp.append((peptide,dataset))

        # initialize peptide dict
        for p_d in peptides_temp:
            peptide = p_d[0]  # Extract the peptide sequence from the tuple
            dataset = p_d[1]  # Extract the dataset from the tuple
            peptides_atomic_max_l_dict[p_d] = len(peptide)  # Initialize the dictionary with the length of each peptide

        peptides_to_keep = set(peptides_temp)  # Create a set to track peptides to keep
        for pd1 in list(peptides_temp):  # Iterate over the peptides in peptides_temp
            peptide1 = pd1[0]  # Extract the peptide sequence from the tuple
            dataset1 = pd1[1]  # Extract the dataset from the tuple
            if len(peptide1) < 9:
                peptides_to_keep.discard(pd1)  # Mark p1 for removal if its length is less than 9
                continue
            for pd2 in list(peptides_temp):  # Compare p1 with other peptides in peptides_temp
                peptide2 = pd2[0]  # Extract the peptide sequence from the tuple
                dataset2 = pd2[1]  # Extract the dataset from the tuple

                if peptide1 in peptide2 and peptide1 != peptide2 and dataset1 == dataset2:
                    # Check if p1 is contained within p2 and is not the same as p2
                    if peptides_atomic_max_l_dict[pd1] < len(peptide2):
                        # print(f"Peptide {p1} is contained in {p2}")
                        peptides_atomic_max_l_dict[pd1] = len(peptide2)

                        if pd1 not in peptides_temp_atomic:
                            peptides_temp_atomic.append(pd1)
                        peptides_to_keep.discard(pd1)  # Mark p1 for removal if it is smaller and overlaps with p2
                    continue

        peptides_temp = list(peptides_to_keep)  # Update peptides_temp with the peptides to keep
        peptides_temp = list(set(peptides_temp))
        
        # print(peptides_temp)
    
        data['n_noncontained_le9'] = len(peptides_temp)
    #peptide,dataset, scan,reason_mismatch,evidence,source
        for peptide, dataset, _, _, evidence, _ in protein_data[protein]['evidence_peptides_scans']:
            if (peptide, dataset) in peptides_temp:
                for i, (pep, ds, scan, reason, _, _) in enumerate(protein_data[protein]['evidence_peptides_scans']):
                    if (pep, ds) == (peptide, dataset):
                        protein_data[protein]['evidence_peptides_scans'][i] = (pep, ds, scan, reason, True, 'Missed')
    # Create the new dataframe
    df_protein_summary = pd.DataFrame(protein_data.values())

    matched_protein_data = {}

    for idx, row in df2.iterrows():
        proteins = row['precursor_proteins']
        if pd.isna(proteins):
            continue

        protein = re.split(r'[;-]', proteins)[0]
        peptide = row['peptide_demod']
        if len(peptide) < 9:
            continue

        scan = row['scan']
        dataset = scan_to_dataset.get(scan, 'Unknown')
        reason_mismatch = row['reason_mismatch']
        evidence = False
        source = 'Matched'

        if protein not in matched_protein_data:
            matched_protein_data[protein] = {
                'precursor_protein': protein,
                'precursor_id_number': 0,
                'evidence_peptides_scans': []
            }

        matched_protein_data[protein]['precursor_id_number'] += 1
        matched_protein_data[protein]['evidence_peptides_scans'].append(
            (peptide, dataset, scan, reason_mismatch, evidence, source)
        )

    for protein, data in matched_protein_data.items():
        peptides_temp = []
        peptides_temp_atomic = []
        peptides_atomic_max_l_dict = {}

        for peptide, dataset, _, _, _, _ in data['evidence_peptides_scans']:
            peptides_temp.append((peptide, dataset))

        for p_d in peptides_temp:
            peptides_atomic_max_l_dict[p_d] = len(p_d[0])

        peptides_to_keep = set(peptides_temp)
        for pd1 in list(peptides_temp):
            peptide1, dataset1 = pd1
            if len(peptide1) < 9:
                peptides_to_keep.discard(pd1)
                continue

            for pd2 in list(peptides_temp):
                peptide2, dataset2 = pd2
                if peptide1 in peptide2 and peptide1 != peptide2 and dataset1 == dataset2:
                    if peptides_atomic_max_l_dict[pd1] < len(peptide2):
                        peptides_atomic_max_l_dict[pd1] = len(peptide2)
                        if pd1 not in peptides_temp_atomic:
                            peptides_temp_atomic.append(pd1)
                        peptides_to_keep.discard(pd1)

        peptides_temp = list(set(peptides_to_keep))
        data['n_noncontained_le9'] = len(peptides_temp)

        for peptide, dataset, _, _, _, _ in data['evidence_peptides_scans']:
            if (peptide, dataset) in peptides_temp:
                for i, (pep, ds, scan, reason, _, _) in enumerate(data['evidence_peptides_scans']):
                    if (pep, ds) == (peptide, dataset):
                        data['evidence_peptides_scans'][i] = (pep, ds, scan, reason, True, 'Matched')

    df_matched_protein_summary = pd.DataFrame(matched_protein_data.values())

    mismatch_summary = df_protein_summary.rename(columns={
        'precursor_id_number': 'mismatch_precursor_id_number',
        'evidence_peptides_scans': 'mismatch_evidence_peptides_scans',
        'n_noncontained_le9': 'mismatch_n_noncontained_le9'
    })

    matched_summary = df_matched_protein_summary.rename(columns={
        'precursor_id_number': 'matched_precursor_id_number',
        'evidence_peptides_scans': 'matched_evidence_peptides_scans',
        'n_noncontained_le9': 'matched_n_noncontained_le9'
    })

    combined_protein_summary = pd.merge(
        mismatch_summary,
        matched_summary,
        on='precursor_protein',
        how='outer'
    )

    def _ensure_list(x):
        return x if isinstance(x, list) else []

    combined_protein_summary['mismatch_evidence_peptides_scans'] = \
        combined_protein_summary['mismatch_evidence_peptides_scans'].apply(_ensure_list)
    combined_protein_summary['matched_evidence_peptides_scans'] = \
        combined_protein_summary['matched_evidence_peptides_scans'].apply(_ensure_list)

    combined_protein_summary['evidence_peptides_scans'] = (
        combined_protein_summary['mismatch_evidence_peptides_scans']
        + combined_protein_summary['matched_evidence_peptides_scans']
    )

    combined_protein_summary.to_csv('combined_protein_level.tsv', sep='\t', index=False)

    def _atomic_noncontained(peptide_records):
        peptide_dataset_pairs = []

        for rec in peptide_records:
            if not isinstance(rec, (tuple, list)) or len(rec) < 2:
                continue

            # Expected format: (peptide, dataset, scan, reason, evidence, source)
            if len(rec) >= 6:
                peptide, dataset = rec[0], rec[1]
            # Backward-compatibility for old format: (peptide, scan, reason, evidence, source)
            elif len(rec) == 5:
                peptide, scan = rec[0], rec[1]
                dataset = scan_to_dataset.get(scan, 'Unknown')
            else:
                continue

            if isinstance(peptide, str) and len(peptide) >= 9:
                peptide_dataset_pairs.append((peptide, dataset))

        peptides_to_keep = set(peptide_dataset_pairs)

        for p1, d1 in peptide_dataset_pairs:
            for p2, d2 in peptide_dataset_pairs:
                # Only compare within the same dataset
                if d1 != d2:
                    continue
                if p1 != p2 and p1 in p2 and len(p1) < len(p2):
                    peptides_to_keep.discard((p1, d1))
                    break

        return peptides_to_keep


    combined_n_noncontained = []
    updated_evidence_records = []

    for _, row in combined_protein_summary.iterrows():
        records = row['evidence_peptides_scans']
        if not isinstance(records, list):
            records = []

        atomic_pairs = _atomic_noncontained(records)
        combined_n_noncontained.append(len(atomic_pairs))

        updated_records = []
        for rec in records:
            if not isinstance(rec, (tuple, list)):
                continue

            if len(rec) >= 6:
                pep, dataset, scan, reason, _, source = rec[:6]
            elif len(rec) == 5:
                pep, scan, reason, _, source = rec
                dataset = scan_to_dataset.get(scan, 'Unknown')
            else:
                continue

            evidence = (pep, dataset) in atomic_pairs
            updated_records.append((pep, dataset, scan, reason, evidence, source))

        updated_evidence_records.append(updated_records)

    combined_protein_summary['combined_n_noncontained_le9'] = combined_n_noncontained
    combined_protein_summary['evidence_peptides_scans'] = pd.Series(
        updated_evidence_records,
        index=combined_protein_summary.index,
        dtype='object'
    )


    combined_protein_summary.to_csv('combined_protein_level.tsv', sep='\t', index=False)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 4:
        print("Usage: python get_protein_level.py <mismatched.tsv> <matched.tsv> <spectra.mgf>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])