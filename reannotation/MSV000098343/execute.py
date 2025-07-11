import os
import pandas as pd
import re

with open("mzml_files.txt", "r") as file:
    file_names = [line.strip() for line in file.readlines()]

# Function to parse filename and extract required fields
def parse_filename(filename): #TODO


    dataset = "MSV000098343"  # Replace with actual logic if needed
    run = filename

    # Further processing to determine Condition and Cohort
    cohort = ''  # Assuming the tissue is the cohort
    condition = '' 
    condition = ''   # Join conditions by comma
    time_point = ''   # Not provided in the filename
    experiment = ''   # Could be inferred if logic is provided
    channel = ''   # Not applicable to mass spectrometry data
    fraction = ''  # Not provided in the filename
    species = ""  # Assuming human cell lines or tissues
    disease = ''   # Not provided in the filename
    enzyme = ''   # Not applicable to mass spectrometry data
    notes = ""  # Placeholder for future notes

    return {
        "Dataset": dataset,
        "Run": run,
        "Cohort": cohort,
        "Subject": '' ,
        "Condition": condition,
        "TimePoint": time_point,
        "BioReplicate": '' ,
        "Experiment": experiment,
        "Batch": '' ,
        "Channel": channel,
        "Fraction": fraction,
        "Species": species,
        "Disease": disease,
        "Tissue": '' ,
        "Enzyme": enzyme,
        "Notes": notes,
    }


# Create a list to hold each parsed file data
data = []

# Process each filename and append parsed data to the list
for filename in file_names:
    parsed_data = parse_filename(filename)
    if parsed_data:
        data.append(parsed_data)

# Create a DataFrame and save to CSV
df = pd.DataFrame(data)
df.to_csv('MSV000098343_annotated.csv', index=False)

print("CSV file has been created: annotated_data.csv")
