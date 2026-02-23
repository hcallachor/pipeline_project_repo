#count the number of reads for each sample before vs after filtering with bowtie2
from Bio import SeqIO
with open(snakemake.input[0], "r") as handle: #count the number of reads from each paired end fastq file for the sample
    num_before = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[1], "r") as handle: #count the reads from the other paired end fastq
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[2], "r") as handle: #count the number of reads after filtering with bowtie2
    num_after = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_after += 1
with open(snakemake.output[0], "a") as report: #append info to report
    report.write("\n" + "Sample SRR5660030 had " + str((num_before)//2) + " read pairs before and " + str((num_after)//2) + " read pairs after Bowtie2 filtering." + "\n")

#repeat for other samples
with open(snakemake.input[3], "r") as handle:
    num_before = 0 #reset the count to zero for next sample
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[4], "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[5], "r") as handle:
    num_after = 0 #reset the count to zero for the next sample
    for record in SeqIO.parse(handle, "fastq"):
        num_after += 1
with open(snakemake.output[0], "a") as report:
    report.write("Sample SRR5660033 had " + str((num_before)//2) + " read pairs before and " + str((num_after)//2) + " read pairs after Bowtie2 filtering." + "\n")

with open(snakemake.input[6], "r") as handle:
    num_before = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[7], "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[8], "r") as handle:
    num_after = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_after += 1
with open(snakemake.output[0], "a") as report:
    report.write("Sample SRR5660044 had " + str((num_before)//2) + " read pairs before and " + str((num_after)//2) + " read pairs after Bowtie2 filtering." + "\n")

with open(snakemake.input[9], "r") as handle:
    num_before = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[10], "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        num_before += 1
with open(snakemake.input[11], "r") as handle:
    num_after = 0
    for record in SeqIO.parse(handle, "fastq"):
        num_after += 1
with open(snakemake.output[0], "a") as report:
    report.write("Sample SRR5660045 had " + str((num_before)//2) + " read pairs before and " + str((num_after)//2) + " read pairs after Bowtie2 filtering." + "\n" + "\n")