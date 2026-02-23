#clean the reference data to only include the protein_id and the sequence, and remove the description
HCMV_data = []
num_records = 0
from Bio import SeqIO
with open(snakemake.input[0], "r") as in_handle:
    for record in SeqIO.parse(in_handle, 'fasta'):
        num_records += 1 #count the number of records in the input file
        if "protein_id=" in record.description:
               protein_id = record.description.split("protein_id=")[1].split("=")[0] #extract the protein_id for each record
               protein_id = protein_id.split("]")[0]
               record.description = "" #clear the description
               record.id = protein_id
        HCMV_data.append(record)
    with open(snakemake.output[0], "w") as out_handle:
            SeqIO.write(HCMV_data, out_handle, "fasta")

#write a report with the number of CDS in the HCMV genome
with open(snakemake.output[1], 'w') as f:
    f.write("The HCMV genome (GCF_000845245.1) has " + str(num_records) + " CDS." + "\n" + "\n")