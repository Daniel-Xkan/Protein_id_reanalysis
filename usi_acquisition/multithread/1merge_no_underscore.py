# Read through the file and replace double underscores with single underscores
input_file = "1merged.txt"
output_file = "1merged_fixed.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        fixed_line = line.replace("__", "_")
        outfile.write(fixed_line)

print(f"Processed file saved as {output_file}")