import os
import pandas as pd

# Directory containing the files
directory = os.path.dirname(__file__)

# Output file
output_file = os.path.join(directory, '1merged.txt')

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file ends with 'usi.txt'
        if filename.endswith('usi.txt'):
            file_path = os.path.join(directory, filename)
            # Open each file in read mode
            with open(file_path, 'r') as infile:
                # Write the contents of the file to the output file
                outfile.write(infile.read())
                outfile.write('\n')  # Add a newline character after each file's content``


# File paths
observations_file = os.path.join(directory, '1observations.csv')
peptides_file = os.path.join(directory, '1PA2024_HPP_peptides.txt')
difference_file = os.path.join(directory, 'difference.txt')

# Read the CSV and TXT files into dataframes
observations_df = pd.read_csv(observations_file)
peptides_df = pd.read_csv(peptides_file, sep='\t', header=None)

# Find the difference between the two dataframes
difference_df = pd.concat([observations_df, peptides_df]).drop_duplicates(keep=False)

# Write the difference to a file
difference_df.to_csv(difference_file, index=False, header=False)