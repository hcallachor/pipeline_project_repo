#identify the longest contig from each assembly
from Bio import SeqIO
with open(snakemake.input[0]) as handle:
    longest_contig = []
    for record in SeqIO.parse(handle, "fasta"): #iterate through the sequences
        if len(record.seq) > len(longest_contig): #keep the longest_contig list updated every time a longer contig is identified
            longest_contig = record.seq
with open(snakemake.output[0], "w") as output_handle:
    output_handle.write(">longest_contig\n")
    output_handle.write(str(longest_contig) + "\n") #save the longest contig for each sample

#repeat for the other assemblies
with open(snakemake.input[1]) as handle:
    longest_contig = [] #reset the contig length to zero for the next sample
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > len(longest_contig):
            longest_contig = record.seq
with open(snakemake.output[1], "w") as output_handle:
    output_handle.write(">longest_contig\n")
    output_handle.write(str(longest_contig) + "\n")

with open(snakemake.input[2]) as handle:
    longest_contig = []
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > len(longest_contig):
            longest_contig = record.seq
with open(snakemake.output[2], "w") as output_handle:
    output_handle.write(">longest_contig\n")
    output_handle.write(str(longest_contig) + "\n")

with open(snakemake.input[3]) as handle:
    longest_contig = []
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > len(longest_contig):
            longest_contig = record.seq
with open(snakemake.output[3], "w") as output_handle:
    output_handle.write(">longest_contig\n")
    output_handle.write(str(longest_contig) + "\n")