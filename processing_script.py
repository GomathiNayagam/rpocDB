import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
import tempfile

# Directories for input and output
genomes_dir = "genomes/"
output_dir = "rpoc_output/"
os.makedirs(output_dir, exist_ok=True)

# Path to HMM profile (only rpoC)
hmm_profile = "rpoC.hmm"

# Function to process each genome
def process_genome(genome_file):
    base_name = os.path.splitext(os.path.basename(genome_file))[0]

    # Temporary files for Prodigal outputs
    with tempfile.NamedTemporaryFile(delete=False, suffix=".faa") as temp_prot, \
         tempfile.NamedTemporaryFile(delete=False, suffix=".fna") as temp_nucl:
        prot_file = temp_prot.name
        nucl_file = temp_nucl.name

    # Run Prodigal to predict proteins and corresponding nucleotide CDS
    prodigal_cmd = [
        "prodigal",
        "-i", genome_file,
        "-a", prot_file,
        "-d", nucl_file,
        "-p", "meta"
    ]
    subprocess.run(prodigal_cmd, check=True)

    # Run HMMER for rpoC
    hmm_output_file = os.path.join(output_dir, f"{base_name}_rpoC.hmmsearch")
    hmm_cmd = [
        "hmmsearch",
        "--cpu", "6",
        "-E", "1e-20",  
        "--tblout", hmm_output_file,
        hmm_profile,
        prot_file
    ]
    subprocess.run(hmm_cmd, check=True)

    # Extract matched nucleotide CDS sequences (not proteins)
    extract_matching_nucleotide_seqs(hmm_output_file, nucl_file, f"{base_name}_rpoC")

    # Clean up temporary files
    os.remove(prot_file)
    os.remove(nucl_file)

# Function to extract nucleotide sequences based on HMM hits on proteins
def extract_matching_nucleotide_seqs(hmm_output_file, nucl_file, output_prefix):
    hits = parse_hmm_hits(hmm_output_file)

    output_nucl_file = os.path.join(output_dir, f"{output_prefix}_retrieved.fna")

    with open(output_nucl_file, 'w') as out_handle:
        for record in SeqIO.parse(nucl_file, "fasta"):
            if record.id in hits:
                SeqIO.write(record, out_handle, "fasta")

# Function to parse HMMER hits from --tblout format
def parse_hmm_hits(hmm_output_file):
    hits = set()
    with open(hmm_output_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                hits.add(line.split()[0])  # Protein ID that matched the HMM
    return hits

# Main
if __name__ == "__main__":
    genome_files = [os.path.join(genomes_dir, f) for f in os.listdir(genomes_dir) if f.endswith('.fna')]

    with ProcessPoolExecutor(max_workers=28) as executor:
        executor.map(process_genome, genome_files)
