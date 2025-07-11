import pandas as pd
import argparse

def create_dataset_tier_dictionary(file_path):
    """
    Read TSV file and create a dictionary mapping datasets to tiers.
    
    Tier 1: MSV column contains 'MSV'
    Tier 2: No 'MSV' in MSV column but 'import' in Action column
    Tier 3: Otherwise
    """
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    dataset_tier_dict = {}
    
    for _, row in df.iterrows():
        dataset = row['Dataset']
        msv_value = str(row['MSV']) if pd.notna(row['MSV']) else ''
        action_value = str(row['Action']) if pd.notna(row['Action']) else ''
        
        # Determine tier based on conditions
        if 'MSV' in msv_value:
            tier = 1
        elif 'import' in action_value.lower():
            tier = 2
        else:
            tier = 3
            
        dataset_tier_dict[dataset] = tier
    
    return dataset_tier_dict

# Read the file and create the dictionary
dataset_tier_dict = create_dataset_tier_dictionary('all_usi_by_datasets.tsv')


def calculate_tier_portions(usi_file_path, dataset_tier_dict):
    """
    Calculate and print the portion of tier 1 USIs vs (tier 2 + tier 3) USIs.
    """
    tier_counts = {1: 0, 2: 0, 3: 0}
    
    with open(usi_file_path, 'r') as file:
        for line in file:
            usi = line.strip()
            if usi:
                # Extract dataset from USI (format: mzspec:DATASET:...)
                parts = usi.split(':')
                if len(parts) >= 2:
                    dataset = parts[1]
                    tier = dataset_tier_dict.get(dataset, 3)  # Default to tier 3 if not found
                    tier_counts[tier] += 1
    
    total_tier1 = tier_counts[1]
    total_tier2_tier3 = tier_counts[2] + tier_counts[3]
    
    print(f"\nTier counts:")
    print(f"Tier 1: {total_tier1}")
    print(f"Tier 2: {tier_counts[2]}")
    print(f"Tier 3: {tier_counts[3]}")
    print(f"Total Tier 2 + Tier 3: {total_tier2_tier3}")
    
    total_all = total_tier1 + total_tier2_tier3
    if total_all > 0:
        portion = total_tier1 / total_all
        print(f"\nPortion of Tier 1 / Total: {portion:.4f}")
    else:
        print("\nNo USIs found.")

# Add argument parser
parser = argparse.ArgumentParser(description='Calculate tier portions from USI file')
parser.add_argument('-i', '--input', required=True, help='Input USI file path')
args = parser.parse_args()

# Call the function with the input file parameter
calculate_tier_portions(args.input, dataset_tier_dict)
