#define the Pipeline Report as the final output of the workflow
rule all:
    input:
        "PipelineReport.txt"
#step 2 retrieve reference transcriptome for HCMV with datasets command line tool
rule retrieve_transcriptome:
    output:
        "HCMV_data/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna",
        "HCMV_data/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
    shell:
        """
        datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report --filename HCMV_data/ncbi_dataset.zip
        unzip HCMV_data/ncbi_dataset.zip -d HCMV_data
        """
#step 2 clean the reference fasta file to have only the protein ID in the headers
#additionally, counts the number of transcripts in the reference and adds that information to the pipeline report
rule clean_ref:
    input:
        "HCMV_data/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna",
    output:
        "data/HCMV.fasta",
        "results/clean_ref_report.txt"
    script:
        "scripts/clean_ref.py"
#build a transcriptome index for HCMV with kallisto for quantification of TPM values in the next step
rule index:
    input:
        "data/HCMV.fasta"
    output:
        "HCMV_index.idx"
    shell:
        "kallisto index -i {output} {input}"
#quantify TPM with kallisto
rule quantify:
    input:
        index="HCMV_index.idx",
        reads1="data/SRR5660030_1.fastq",
        reads2="data/SRR5660030_2.fastq",
        reads3="data/SRR5660033_1.fastq",
        reads4="data/SRR5660033_2.fastq",
        reads5="data/SRR5660044_1.fastq",
        reads6="data/SRR5660044_2.fastq",
        reads7="data/SRR5660045_1.fastq",
        reads8="data/SRR5660045_2.fastq"
    output:
        "results/SRR5660030/abundance.tsv",
        "results/SRR5660033/abundance.tsv",
        "results/SRR5660044/abundance.tsv",
        "results/SRR5660045/abundance.tsv",
        touch("results/quant_done.txt")
    shell:
        """
        kallisto quant -i {input.index} -o results/SRR5660030 -b 30 -t 4 {input.reads1} {input.reads2}
        kallisto quant -i {input.index} -o results/SRR5660033 -b 30 -t 4 {input.reads3} {input.reads4}
        kallisto quant -i {input.index} -o results/SRR5660044 -b 30 -t 4 {input.reads5} {input.reads6}
        kallisto quant -i {input.index} -o results/SRR5660045 -b 30 -t 4 {input.reads7} {input.reads8}
        """
#prepare kallisto results for sleuth, python script creates a single tsv file with the data table needed for sleuth analysis
rule prepare_sleuth:
    input:
        "results/SRR5660030/abundance.tsv",
        "results/SRR5660033/abundance.tsv",
        "results/SRR5660044/abundance.tsv",
        "results/SRR5660045/abundance.tsv",
        "results/quant_done.txt" #use the touch file from the quantify step as input to ensure this step waits for all quantification to be done before running
    output:
        "results/sleuth_input.tsv"
    script:
        "scripts/prepare_sleuth.py"
#find differential expression between the 2 timepoints with sleuth
rule sleuth:
    input:
        "results/sleuth_input.tsv"
    output:
        "results/sleuth_significant.tsv"
    script:
        "scripts/sleuth_analysis.R"
#create a genome index for HCMV with bowtie2
rule bowtie_index:
    input:
        "HCMV_data/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna"
    output:
        "HCMV_index.1.bt2",
        "HCMV_index.2.bt2",
        "HCMV_index.3.bt2",
        "HCMV_index.4.bt2"
    shell:
        "bowtie2-build {input} HCMV_index"
#map reads from the SRR paired end reads to the HCMV genome index with bowtie2
rule bowtie_map:
    input:
        reads1="data/SRR5660030_1.fastq",
        reads2="data/SRR5660030_2.fastq",
        reads3="data/SRR5660033_1.fastq",
        reads4="data/SRR5660033_2.fastq",
        reads5="data/SRR5660044_1.fastq",
        reads6="data/SRR5660044_2.fastq",
        reads7="data/SRR5660045_1.fastq",
        reads8="data/SRR5660045_2.fastq",
        index="HCMV_index.1.bt2" #use one of the bowtie2 index files as input to ensure bowtie2 waits for the index to be built before mapping
    output:
        "results/SRR5660030/bowtie2_alignment.sam",
        "results/SRR5660033/bowtie2_alignment.sam",
        "results/SRR5660044/bowtie2_alignment.sam",
        "results/SRR5660045/bowtie2_alignment.sam"
    shell:
        """
        bowtie2 --no-unal -x HCMV_index -1 {input.reads1} -2 {input.reads2} -S results/SRR5660030/bowtie2_alignment.sam
        bowtie2 --no-unal -x HCMV_index -1 {input.reads3} -2 {input.reads4} -S results/SRR5660033/bowtie2_alignment.sam
        bowtie2 --no-unal -x HCMV_index -1 {input.reads5} -2 {input.reads6} -S results/SRR5660044/bowtie2_alignment.sam
        bowtie2 --no-unal -x HCMV_index -1 {input.reads7} -2 {input.reads8} -S results/SRR5660045/bowtie2_alignment.sam
        """

#convert bowtie2 sam files to fastq files for spades
rule sam_to_fastq:
    input:
        "results/SRR5660030/bowtie2_alignment.sam",
        "results/SRR5660033/bowtie2_alignment.sam",
        "results/SRR5660044/bowtie2_alignment.sam",
        "results/SRR5660045/bowtie2_alignment.sam"
    output:
        "results/SRR5660030/filtered_reads.fastq",
        "results/SRR5660033/filtered_reads.fastq",
        "results/SRR5660044/filtered_reads.fastq",
        "results/SRR5660045/filtered_reads.fastq"
    shell:
        """
        samtools fastq results/SRR5660030/bowtie2_alignment.sam > results/SRR5660030/filtered_reads.fastq
        samtools fastq results/SRR5660033/bowtie2_alignment.sam > results/SRR5660033/filtered_reads.fastq
        samtools fastq results/SRR5660044/bowtie2_alignment.sam > results/SRR5660044/filtered_reads.fastq
        samtools fastq results/SRR5660045/bowtie2_alignment.sam > results/SRR5660045/filtered_reads.fastq
        """
#count the number of unfiltered vs filtered reads, python script adds this information to the pipeline report
rule count_reads:
    input:
        "data/SRR5660030_1.fastq",
        "data/SRR5660030_2.fastq",
        "results/SRR5660030/filtered_reads.fastq",
        "data/SRR5660033_1.fastq",
        "data/SRR5660033_2.fastq",
        "results/SRR5660033/filtered_reads.fastq",
        "data/SRR5660044_1.fastq",
        "data/SRR5660044_2.fastq",
        "results/SRR5660044/filtered_reads.fastq",
        "data/SRR5660045_1.fastq",
        "data/SRR5660045_2.fastq",
        "results/SRR5660045/filtered_reads.fastq"
    output:
        "results/read_counts.txt"
    script:
        "scripts/count_reads.py"
#build assemblies for each SRR sample with SPAdes
#using kmer size 127
rule spades:
    input:
        "results/SRR5660030/filtered_reads.fastq",
        "results/SRR5660033/filtered_reads.fastq",
        "results/SRR5660044/filtered_reads.fastq",
        "results/SRR5660045/filtered_reads.fastq"
    output:
        "results/SRR5660030/spades_output/contigs.fasta",
        "results/SRR5660033/spades_output/contigs.fasta",
        "results/SRR5660044/spades_output/contigs.fasta",
        "results/SRR5660045/spades_output/contigs.fasta"
    shell:
        """
        spades.py -k 127 -t 2 --only-assembler -s results/SRR5660030/filtered_reads.fastq -o results/SRR5660030/spades_output
        spades.py -k 127 -t 2 --only-assembler -s results/SRR5660033/filtered_reads.fastq -o results/SRR5660033/spades_output
        spades.py -k 127 -t 2 --only-assembler -s results/SRR5660044/filtered_reads.fastq -o results/SRR5660044/spades_output
        spades.py -k 127 -t 2 --only-assembler -s results/SRR5660045/filtered_reads.fastq -o results/SRR5660045/spades_output
        """
#identify longest contig in each assembly to be used for blast analysis
rule longest_contig:
    input:
        "results/SRR5660030/spades_output/contigs.fasta",
        "results/SRR5660033/spades_output/contigs.fasta",
        "results/SRR5660044/spades_output/contigs.fasta",
        "results/SRR5660045/spades_output/contigs.fasta"
    output:
        "results/SRR5660030/spades_output/longest_contig.fasta",
        "results/SRR5660033/spades_output/longest_contig.fasta",
        "results/SRR5660044/spades_output/longest_contig.fasta",
        "results/SRR5660045/spades_output/longest_contig.fasta"
    script:
        "scripts/longest_contig.py"
#obtain data to make local database of Betaherpesvirinae subfamily for blast analysis
rule local_db_data:
    output:
        "BHV_data/ncbi_dataset/data/genomic.fna"
    shell:
        """
        datasets download virus genome taxon betaherpesvirinae --refseq --include genome --filename BHV_data/ncbi_dataset.zip
        unzip BHV_data/ncbi_dataset.zip -d BHV_data
        """
#create local blast database
rule make_blast_db:
    input:
        "BHV_data/ncbi_dataset/data/genomic.fna"
    output:
        "betaherpesvirinae.ndb",
    shell:
        "makeblastdb -in BHV_data/ncbi_dataset/data/genomic.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl"
#blast longest contigs against local BHV database
rule blast:
    input:
        "results/SRR5660030/spades_output/longest_contig.fasta",
        "results/SRR5660033/spades_output/longest_contig.fasta",
        "results/SRR5660044/spades_output/longest_contig.fasta",
        "results/SRR5660045/spades_output/longest_contig.fasta",
        "betaherpesvirinae.ndb" #use one of the blast db files as input to ensure blast waits for the database to be built before running
    output:
        "results/SRR5660030/blast_results.tsv",
        "results/SRR5660033/blast_results.tsv",
        "results/SRR5660044/blast_results.tsv",
        "results/SRR5660045/blast_results.tsv"
    shell:
        """
        blastn -task megablast -query results/SRR5660030/spades_output/longest_contig.fasta -db betaherpesvirinae -out results/SRR5660030/blast_results.tsv -max_hsps 1 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' 
        blastn -task megablast -query results/SRR5660033/spades_output/longest_contig.fasta -db betaherpesvirinae -out results/SRR5660033/blast_results.tsv -max_hsps 1 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'
        blastn -task megablast -query results/SRR5660044/spades_output/longest_contig.fasta -db betaherpesvirinae -out results/SRR5660044/blast_results.tsv -max_hsps 1 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'
        blastn -task megablast -query results/SRR5660045/spades_output/longest_contig.fasta -db betaherpesvirinae -out results/SRR5660045/blast_results.tsv -max_hsps 1 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'
        """
add blast results to pipeline report
rule add_blast_to_report:
    input:
        "results/SRR5660030/blast_results.tsv",
        "results/SRR5660033/blast_results.tsv",
        "results/SRR5660044/blast_results.tsv",
        "results/SRR5660045/blast_results.tsv"
    output:
        "results/blast_summary.txt",
    script:
        "scripts/add_blast_to_report.py"
#combine all report information into a single pipeline report text file
rule make_report:
    input:
        report1="results/clean_ref_report.txt",
        report2="results/sleuth_significant.tsv",
        report3="results/read_counts.txt",
        report4="results/blast_summary.txt"
    output:
        "PipelineReport.txt"
    shell:
        """
        cat {input.report1} > {output}
        cat {input.report2} >> {output}
        cat {input.report3} >> {output}
        cat {input.report4} >> {output}
        """