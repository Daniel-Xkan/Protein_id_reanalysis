#!/usr/bin/env python3
"""
Script to read PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv file

This script reads the TSV file and displays basic information about its contents.
"""

import pandas as pd
import os
import sys


def read_prismal_tsv(file_path):
    """
    Read the PRISMAL TSV file and return a pandas DataFrame.
    
    Args:
        file_path (str): Path to the TSV file
        
    Returns:
        pd.DataFrame: The loaded dataframe
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    try:
        # Read the TSV file
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
        return df
    except Exception as e:
        raise Exception(f"Error reading TSV file: {str(e)}")


def display_file_info(df, file_path):
    """
    Display basic information about the loaded dataframe.
    
    Args:
        df (pd.DataFrame): The dataframe to analyze
        file_path (str): Path to the original file
    """
    print(f"\n{'='*80}")
    print(f"File: {file_path}")
    print(f"{'='*80}")
    print(f"\nNumber of rows: {len(df)}")
    print(f"Number of columns: {len(df.columns)}")
    print(f"\nColumn names:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")
    
    print(f"\nFirst few rows:")
    print(df.head())
    
    print(f"\nData types:")
    print(df.dtypes)
    
    print(f"\nBasic statistics:")
    print(df.describe())
    print(f"\n{'='*80}\n")


def main():
    """Main function to read and display TSV file information."""
    # Default file name
    default_file = "PRISMAL-REPROCESS_PA_UniproKB_Normal123.tsv"
    
    # Check if file path is provided as argument
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    else:
        file_path = default_file
    
    try:
        print(f"Reading TSV file: {file_path}")
        df = read_prismal_tsv(file_path)
        display_file_info(df, file_path)
        
        return df
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print(f"\nUsage: python {sys.argv[0]} [path_to_tsv_file]")
        print(f"Default file: {default_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    df = main()
