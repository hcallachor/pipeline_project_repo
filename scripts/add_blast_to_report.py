#add top 5 blast hits from each assembly to the pipeline report
with open(snakemake.output[0], "w") as report: #append data to report
    report.write("SRR5660030:\n") #first header
    report.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\n") #second header, tab delimited
    with open(snakemake.input[0], "r") as blast_file: #add tab delimited data from the blast results
        for i, line in enumerate(blast_file):
            if i <= 5: #only add top 5 hits
                report.write(line + "\n")
    report.write("SRR5660033:\n") #repeat for other samples
    report.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\n")
    with open(snakemake.input[1], "r") as blast_file:
        for i, line in enumerate(blast_file):
            if i <= 5:
                report.write(line + "\n")
    report.write("SRR5660044:\n")
    report.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\n")
    with open(snakemake.input[2], "r") as blast_file:
        for i, line in enumerate(blast_file):
            if i <= 5:
                report.write(line + "\n")
    report.write("SRR5660045:\n")
    report.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle\n")
    with open(snakemake.input[3], "r") as blast_file:
        for i, line in enumerate(blast_file):
            if i <= 5:
                report.write(line + "\n")