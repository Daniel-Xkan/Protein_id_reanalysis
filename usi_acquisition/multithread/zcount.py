def count_unique_proteins(file_path):
    unique_proteins = set()
    
    with open(file_path, 'r') as f:
        for line in f:
            # Extract the protein ID (the first column)
            protein_id = line.strip().split('\t')[0]
            unique_proteins.add(protein_id)
    
    print(f"Total unique proteins: {len(unique_proteins)}")

# Usage example
file_path = 'PA2024_HPP_peptides.txt'
count_unique_proteins(file_path)
