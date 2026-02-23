#test for differential expression between the 2 timepoints (2dpi and 6dpi) using sleuth
library(sleuth)

input_file <- snakemake@input[[1]]
#read in the sleuth input file
stab=read.table(input_file, header=TRUE)
so=sleuth_prep(stab)

#fit the full and reduced models
so=sleuth_fit(so, ~condition, 'full')
so=sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test
so=sleuth_lrt(so, 'reduced', 'full')
library(dplyr)
sleuth_table=sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#filter the results to only include significant results
sleuth_significant=dplyr::filter(sleuth_table, qval < 0.05)

#write the significant results to a file
#header should include target_id, test_stat, pval, qval
output_file <- snakemake@output[[1]]
write.table(sleuth_significant[,c("target_id", "test_stat", "pval", "qval")], file=output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)