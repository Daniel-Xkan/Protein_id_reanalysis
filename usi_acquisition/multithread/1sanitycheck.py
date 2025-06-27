import csv

# input_file = '1PeptideAtlas_peptide_observations.csv'
# output_file = '1observations.csv'

# # Read the input CSV file
# with open(input_file, mode='r', newline='') as infile:
#     reader = csv.DictReader(infile)
#     unique_sequences = {}
    
#     # Process each row in the input file
#     for row in reader:
#         sequence = row['Sequence']
#         if sequence not in unique_sequences:
#             unique_sequences[sequence] = {
#                 'Protein': row['Protein'],
#                 'Paref': row['Paref'],
#                 'Sequence': sequence
#             }

# # Write the output CSV file
# with open(output_file, mode='w', newline='') as outfile:
#     fieldnames = ['Protein', 'Paref', 'Sequence']
#     writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    
#     writer.writeheader()
#     for sequence in unique_sequences.values():
#         writer.writerow(sequence)
###############################################################################
# # Read the input TXT file
# txt_input_file = '1PA2024_HPP_peptides.txt'
# csv_output_file = '1PA2024.csv'

# with open(txt_input_file, mode='r') as txtfile:
#     lines = txtfile.readlines()
#     header = lines[0].strip().split()
#     data = [line.strip().split() for line in lines[1:]]

# # Write the output CSV file
# with open(csv_output_file, mode='w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(header)
#     writer.writerows(data)

###############################################################################
## Compare the differences between 1PAPA2024.csv and 1observations.csv
# observations_file = '1observations.csv'
# papa2024_file = '1PA2024.csv'
new_peptides_file = '1new_peptides.csv'

# # Read the observations CSV file
# with open(observations_file, mode='r', newline='') as obsfile:
#     obs_reader = csv.DictReader(obsfile)
#     observations_sequences = {row['Sequence'] for row in obs_reader}

# # Read the PAPA2024 CSV file
# with open(papa2024_file, mode='r', newline='') as papafile:
#     papa_reader = csv.DictReader(papafile)
#     papa_sequences = {row['Sequence'] for row in papa_reader}

# # Find sequences that are in observations but not in PAPA2024
# unique_sequences = observations_sequences - papa_sequences

# # Write the unique sequences to the new peptides CSV file
# with open(observations_file, mode='r', newline='') as obsfile:
#     obs_reader = csv.DictReader(obsfile)
#     fieldnames = obs_reader.fieldnames

#     with open(new_peptides_file, mode='w', newline='') as newfile:
#         writer = csv.DictWriter(newfile, fieldnames=fieldnames)
#         writer.writeheader()
        
#         for row in obs_reader:
#             if row['Sequence'] in unique_sequences:
#                 writer.writerow(row)
##########################################################################
# Convert the new peptides CSV file to a TXT file
new_peptides_txt_file = '1new_peptides.txt'

with open(new_peptides_file, mode='r', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    fieldnames = reader.fieldnames

    with open(new_peptides_txt_file, mode='w') as txtfile:
        txtfile.write('\t'.join(fieldnames) + '\n')
        for row in reader:
            txtfile.write('\t'.join(row[field] for field in fieldnames) + '\n')