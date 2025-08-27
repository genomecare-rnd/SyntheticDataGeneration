# Algorithms for Synthetic Data Generation in Fetal Aneuploidy Detection

## Formatting Input Files
1. **Trim** single-end `.Fastq` data to the first 35 nucleotides  
2. **Align** the trimmed reads to the *hg19* reference index using **Bowtie2**  
3. **Convert** SAM to BAM, then **sort** and **index** the alignment files  
4. **Remove duplicate reads** using `samtools rmdup -s` (single-end mode)  
5. **Filter** for uniquely mapped reads and save them into the final ***.unique*** file  
6. For multiple samples, generate a ***.unique*** file for each sample, and **record the fetal fraction (FF) and GC content values** for each through a separate process  
