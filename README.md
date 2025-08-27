# Algorithms for synthetic data generation in fetal aneuploidy detection.

# Input File
1. Trim single-end .Fastq data to the first 35 nucleotides
2. Align trimmed reads to the hg19 reference index using Bowtie2
3. Convert SAM to BAM, then sort and index the alignment files
4. Remove duplicate reads using samtools rmdup -s (single-end mode)
5. Filter out only uniquely mapped reads and save them into the final *.unique* file
