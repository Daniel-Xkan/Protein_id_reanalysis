import re
import os
import sys

def correct_usi(usi):
    # Check if USI starts with specific dataset identifiers
    if "mzspec:PXD011967:" in usi:
        # Find the filename part (between PXD011967: and :scan:)
        match = re.search(r'mzspec:PXD011967:([^:]+):scan:', usi)
        if match:
            filename = match.group(1)
            # Remove all underscores from filename
            corrected_filename = filename.replace('_', '')
            # Replace the original filename with the corrected one
            return usi.replace(f"mzspec:PXD011967:{filename}:", f"mzspec:PXD011967:{corrected_filename}:")
    
    elif "mzspec:PXD002098:" in usi:
        # Find the filename part (between PXD002098: and :scan:)
        match = re.search(r'mzspec:PXD002098:([^:]+):scan:', usi)
        if match:
            filename = match.group(1)
            # Look for pattern where we need to add parentheses around the last digit before the final underscore
            pattern_match = re.search(r'(.+)(\d)(_\d+)$', filename)
            if pattern_match:
                prefix = pattern_match.group(1)
                digit = pattern_match.group(2)
                suffix = pattern_match.group(3)
                # Add parentheses around the digit
                corrected_filename = f"{prefix}({digit}){suffix}"
                # Replace the original filename with the corrected one
                return usi.replace(f"mzspec:PXD002098:{filename}:", f"mzspec:PXD002098:{corrected_filename}:")
    
    elif "mzspec:PXD004352:" in usi:
        # Find the filename part (between PXD004352: and :scan:)
        match = re.search(r'mzspec:PXD004352:([^:]+):scan:', usi)
        if match:
            filename = match.group(1)
            # Remove the last underscore in the filename
            if '_' in filename:
                # Find the last underscore and remove it
                last_underscore_pos = filename.rfind('_')
                corrected_filename = filename[:last_underscore_pos] + filename[last_underscore_pos+1:]
                # Replace the original filename with the corrected one
                return usi.replace(f"mzspec:PXD004352:{filename}:", f"mzspec:PXD004352:{corrected_filename}:")
            
    elif "mzspec:PXD014058:" in usi:
        # Find the filename part (between PXD014058: and :scan:)
        match = re.search(r'mzspec:PXD014058:([^:]+):scan:', usi)
        if match:
            filename = match.group(1)
            # Replace "-" with "&" in filename
            corrected_filename = filename.replace('-', '&')
            # Replace the original filename with the corrected one
            return usi.replace(f"mzspec:PXD014058:{filename}:", f"mzspec:PXD014058:{corrected_filename}:")
        
    elif "mzspec:PXD019643:" in usi:
        # Find the filename part (between PXD019643: and :scan:)
        match = re.search(r'mzspec:PXD019643:([^:]+):scan:', usi)
        if match:
            filename = match.group(1)
            # Replace double underscores with single underscore
            corrected_filename = filename.replace('__', '_')
            # Replace the original filename with the corrected one
            return usi.replace(f"mzspec:PXD019643:{filename}:", f"mzspec:PXD019643:{corrected_filename}:")
    
    # If no corrections needed or pattern not matched, return original USI
    return usi

def process_usi_file(input_file):
    # Get base filename without extension
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_usi_fixed.txt"
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line:  # Skip empty lines
                corrected_usi = correct_usi(line)
                outfile.write(corrected_usi + '\n')
    
    print(f"Processed USIs written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python USI_correction.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    process_usi_file(input_file)