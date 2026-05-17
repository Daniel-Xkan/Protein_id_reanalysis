# Reading PRISMAL-REPROCESS TSV Files

This document explains how to read and process the PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv file.

## Overview

The `read_prismal_tsv.py` script provides functionality to read TSV (Tab-Separated Values) files from the PRISMAL reprocessing pipeline, specifically designed for UniproKB protein data.

## Requirements

- Python 3.x
- pandas library

Install pandas if not already installed:
```bash
pip install pandas
```

## Usage

### Basic Usage (Default File)

The script will read `PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv` from the current directory:

```bash
python read_prismal_tsv.py
```

### Specifying a Custom File Path

You can provide a specific file path as an argument:

```bash
python read_prismal_tsv.py /path/to/your/file.tsv
```

## Output

The script displays:
- File name and location
- Number of rows and columns
- Column names (numbered list)
- First few rows of data
- Data types for each column
- Basic statistical summary of numeric columns

## Example Output

```
Reading TSV file: PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv

================================================================================
File: PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv
================================================================================

Number of rows: 5
Number of columns: 5

Column names:
  1. Protein_ID
  2. Gene_Name
  3. Peptide_Sequence
  4. Confidence_Score
  5. Coverage

First few rows:
  Protein_ID Gene_Name           Peptide_Sequence  Confidence_Score  Coverage
0     P12345     BRCA1          MGDVLPDNHYLSTQSAL              0.95      45.2
1     P67890      TP53         MEEPQSDPSVEPPLSQET              0.89      67.8
...
```

## Script Functions

### `read_prismal_tsv(file_path)`
- Reads the TSV file and returns a pandas DataFrame
- Handles file validation and error checking

### `display_file_info(df, file_path)`
- Displays comprehensive information about the loaded data
- Shows column names, data types, and statistics

### `main()`
- Entry point for the script
- Handles command-line arguments
- Manages error handling and user feedback

## Integration with Existing Workflows

This script follows the same patterns used in other parts of the Protein_id_reanalysis repository:
- Uses pandas for data handling (consistent with reanalysis.ipynb notebooks)
- Reads TSV files with `sep='\t'` parameter
- Uses `low_memory=False` for large files

## Error Handling

The script includes comprehensive error handling for:
- Missing files (FileNotFoundError)
- Malformed TSV files
- Invalid file paths
- General exceptions during file reading

## Further Processing

Once the file is loaded into a pandas DataFrame, you can perform additional operations:

```python
from read_prismal_tsv import read_prismal_tsv

# Load the file
df = read_prismal_tsv("PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv")

# Filter data
high_confidence = df[df['Confidence_Score'] > 0.9]

# Export to different format
df.to_csv("output.csv", index=False)

# Merge with other datasets
# merged = pd.merge(df, other_df, on='Protein_ID')
```

## Related Files

- `indexing.ipynb` - Protein indexing and k-mer generation
- `Evaluate_reanalysis/files_needed/reanalysis.ipynb` - Main reanalysis workflow
- `compare_reanalyze_peptide.py` - Peptide comparison scripts
