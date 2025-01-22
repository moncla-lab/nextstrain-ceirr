from Bio import SeqIO
import glob

# List of fasta files
fasta_files = ["ha.fasta", "mp.fasta", "na.fasta", "np.fasta", "ns.fasta", "pa.fasta", "pb1.fasta", "pb2.fasta"]

# Dictionary to store headers and corresponding sequences for each file
seq_dicts = {}
all_headers = set()

# Read sequences from each fasta file
for fasta_file in fasta_files:
    seq_dict = {record.id: record for record in SeqIO.parse(fasta_file, "fasta")}
    seq_dicts[fasta_file] = seq_dict
    if not all_headers:
        all_headers = set(seq_dict.keys())
    else:
        all_headers.intersection_update(seq_dict.keys())

# Sort headers to maintain order
sorted_headers = sorted(all_headers)

# Write filtered sequences to new fasta files
for fasta_file in fasta_files:
    output_file = f"filtered_{fasta_file}"
    with open(output_file, "w") as out_f:
        for header in sorted_headers:
            SeqIO.write(seq_dicts[fasta_file][header], out_f, "fasta")
    print(f"Written: {output_file}")
