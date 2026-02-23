Hello thank you for using my pipeline!
This pipeline is intended to compare the transcriptomes of HCMV (Human Cytomegalovirus) from patient samples 2 and 6 days post infection.

In order to run this pipeline, you must install the Snakemake language, Python, R, Biopython, NCBI Datasets command line tools, NCBI BLAST+ command line tools, sleuth R package, SPAdes genome assembler, Kallisto, SAMtools, and Bowtie2.

The repository includes a completed Pipeline Report for the full samples of SRR5660030, SRR5660033, SRR5660044, and SRR5660045.
The data folder contains much shorter sample files from each SRR sample (each paired end fastq file contains only 10,000 reads), so you can test this pipeline in just a few minutes.
These sample files were created by using the command "head -n 40000 data.fastq > sampledata.fastq", which creates a new file with only the first 10,000 reads from the original data file.

If instead you would like to test the full SRR data (or different SRR samples), it can be retrieved using the wget command line tool and the NCBI website.
First, search for your sample by typing the SRR number (such as SRR5660030) into the NCBI search bar, then open the sequence read archive (SRA) result. In the "Runs" section of the page there will be a link to the run of your SRR. After opening this link, open the "Data access" tab, and copy the link next to "SRA Normalized" (not the webpage link).

Finally, change your directory to the location you want to store your data, and in the command line type wget and your link (such as wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030).
Once you have obtained this data, use "fasterq-dump ./" followed by the SRR number (fasterq-dump ./SRR5660030) to unpack the paired end fastq files.

If you would like to use SRR numbers that are different from those used for this pipeline, you must change the SRR numbers throughout the pipeline and in the scripts.

Finally, to run the code simply ensure you are in the "pipeline_project" directory, and enter into the command line "snakemake -c #", with # being the number of cores you would like to use. The data folder already contains the sample data necessary to run the pipeline.

If you would like to run the pipeline multiple times, you must delete all the outputs created by the pipeline.

Good luck have fun!