import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
import tempfile

# Directories for input and output
genomes_dir = "no_rpoc_genomes/"
output_dir = "second_attempt_rpoc_output/"
os.makedirs(output_dir, exist_ok=True)

# Paths to HMM profiles
hmm_profiles = ["rpoB.hmm", "rpoC.hmm"]

# Function to process each genome
def process_genome(genome_file):
    base_name = os.path.splitext(os.path.basename(genome_file))[0]

    # Temporary files for Prodigal outputs
    with tempfile.NamedTemporaryFile(delete=False, suffix=".faa") as temp_prot, \
         tempfile.NamedTemporaryFile(delete=False, suffix=".fna") as temp_nucl:
        prot_file = temp_prot.name
        nucl_file = temp_nucl.name

    # Run Prodigal to get protein and nucleotide sequences
    prodigal_cmd = [
        "prodigal",
        "-i", genome_file,
        "-a", prot_file,
        "-d", nucl_file,
        "-p", "meta"
    ]
    subprocess.run(prodigal_cmd, check=True)

    # Run HMMER for each profile
    for hmm_profile in hmm_profiles:
        profile_name = os.path.basename(hmm_profile).split(".")[0]
        hmm_output_file = os.path.join(output_dir, f"{base_name}_{profile_name}.hmmsearch")

        hmm_cmd = [
            "hmmsearch",
            "--cpu", str(max(1, 28 // len(hmm_profiles))),  # Allocate CPUs dynamically
            "--tblout", hmm_output_file,
            hmm_profile,
            prot_file
        ]
        subprocess.run(hmm_cmd, check=True)

        # Extract sequences for this HMM profile
        extract_sequences(hmm_output_file, nucl_file, f"{base_name}_{profile_name}")

    # Cleanup temporary files
    os.remove(prot_file)
    os.remove(nucl_file)

# Function to extract sequences based on HMM hits
def extract_sequences(hmm_output_file, nucl_file, output_prefix):
    output_prot_file = os.path.join(output_dir, f"{output_prefix}_retrieved.faa")
    output_nucl_file = os.path.join(output_dir, f"{output_prefix}_retrieved.fna")

    # Parse HMM hits
    hits = parse_hmm_hits(hmm_output_file)

    # Write protein sequences
    with open(output_prot_file, 'w') as prot_handle:
        for record in SeqIO.parse(nucl_file, "fasta"):
            if record.id in hits:
                SeqIO.write(record, prot_handle, "fasta")

    # Write nucleotide sequences
    with open(output_nucl_file, 'w') as nucl_handle:
        for record in SeqIO.parse(nucl_file, "fasta"):
            if record.id in hits:
                SeqIO.write(record, nucl_handle, "fasta")

# Function to parse HMMER hits
def parse_hmm_hits(hmm_output_file):
    hits = set()
    with open(hmm_output_file, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                hits.add(line.split()[0])  # Adjust based on HMMER output format
    return hits

# Main function
if __name__ == "__main__":
    # Collect genome files
    genome_files = [os.path.join(genomes_dir, f) for f in os.listdir(genomes_dir) if f.endswith('.fna')]

    # Parallel processing with multiprocessing
    with ProcessPoolExecutor(max_workers=28) as executor:
        executor.map(process_genome, genome_files)

