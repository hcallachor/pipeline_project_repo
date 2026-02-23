#create a tab delimited file with the following columns: sample, condition, path to kallisto abundance file
import csv

data = [
    ['sample', 'condition', 'path'], #header
    ['SRR5660030', '2dpi', 'results/SRR5660030'], #set the 2dpi samples as control
    ['SRR5660033', '6dpi', 'results/SRR5660033'], #set the 6dpi samples as experimental
    ['SRR5660044', '2dpi', 'results/SRR5660044'],
    ['SRR5660045', '6dpi', 'results/SRR5660045'],
]

tsv_file = snakemake.output[0] #create tab delimited file

with open(tsv_file, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t') #tab delimit the data
    writer.writerows(data)