import os
import sys
from Bio import SeqIO

def parse_hmmsearch_output(hmmsearch_file):
    """Parse HMMER tblout file and extract hit details."""
    hits = []
    with open(hmmsearch_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            # Split the line
            fields = line.strip().split()

            # Ensure we have enough fields
            if len(fields) >= 1:
                try:
                    hits.append({
                        'target_name': fields[0],  # Protein sequence name
                    })
                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {line}")
                    print(f"Error details: {e}")

    print(f"Parsed {len(hits)} hits from {hmmsearch_file}")  # Debugging output
    return hits

def retrieve_sequences(fna_file, hits, output_file, mode='a'):
    """Retrieve sequences for hits from FNA file and write to output."""
    # Read all sequences from the FNA file
    all_sequences = list(SeqIO.parse(fna_file, 'fasta'))

    print(f"Total sequences in {fna_file}: {len(all_sequences)}")  # Debugging output

    # List to track retrieved sequences
    retrieved_sequences = []

    # Open output file in append or write mode
    with open(output_file, mode) as outfile:
        for hit in hits:
            # Debugging: check the hit ID
            print(f"Searching for sequence: {hit['target_name']}")  

            # Find matching sequence by ID
            matching_seq = next((seq for seq in all_sequences if seq.id == hit['target_name']), None)

            if matching_seq:
                # Write sequence to output file
                SeqIO.write(matching_seq, outfile, 'fasta')
                retrieved_sequences.append(matching_seq)
            else:
                print(f"Sequence {hit['target_name']} not found in {fna_file}")  # Debugging output

    print(f"Retrieved {len(retrieved_sequences)} sequences from {os.path.basename(fna_file)}")
    return retrieved_sequences

def process_directory(directory):
    """Process all files in a directory."""
    # Ensure output directory exists
    output_dir = os.path.join(directory, 'retrieved_sequences')
    os.makedirs(output_dir, exist_ok=True)

    # Prepare output files
    rpoB_output = os.path.join(output_dir, 'retrieved_rpoB.fasta')
    rpoC_output = os.path.join(output_dir, 'retrieved_rpoC.fasta')
    missing_rpoB = os.path.join(output_dir, 'missing_rpoB.txt')
    missing_rpoC = os.path.join(output_dir, 'missing_rpoC.txt')

    # Clear previous output files
    open(rpoB_output, 'w').close()
    open(rpoC_output, 'w').close()
    open(missing_rpoB, 'w').close()
    open(missing_rpoC, 'w').close()

    # Find all relevant files
    for filename in os.listdir(directory):
        # Match patterns for input files
        if filename.endswith('.fna'):
            base_name = filename.replace('.fna', '')

            # Check for corresponding hmmsearch files
            rpoB_hmm = os.path.join(directory, f'{base_name}_rpoB.hmmsearch')
            rpoC_hmm = os.path.join(directory, f'{base_name}_rpoC.hmmsearch')

            # Process rpoB
            if os.path.exists(rpoB_hmm):
                rpoB_hits = parse_hmmsearch_output(rpoB_hmm)
                if rpoB_hits:
                    retrieve_sequences(
                        os.path.join(directory, filename), 
                        rpoB_hits, 
                        rpoB_output
                    )
                else:
                    with open(missing_rpoB, 'a') as missing_file:
                        missing_file.write(f"{base_name}\n")

            # Process rpoC
            if os.path.exists(rpoC_hmm):
                rpoC_hits = parse_hmmsearch_output(rpoC_hmm)
                if rpoC_hits:
                    retrieve_sequences(
                        os.path.join(directory, filename), 
                        rpoC_hits, 
                        rpoC_output
                    )
                else:
                    with open(missing_rpoC, 'a') as missing_file:
                        missing_file.write(f"{base_name}\n")

    print("Sequence retrieval complete.")

# Allow running as a script
if __name__ == '__main__':
    # If directory provided as argument
    if len(sys.argv) > 1:
        process_directory(sys.argv[1])
    else:
        print("Please provide a directory path.")
        sys.exit(1)

