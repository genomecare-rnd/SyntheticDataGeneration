# Algorithms for Synthetic Data Generation in Fetal Aneuploidy Detection

## Preparing Input Data
1. **Prepare** multiple single-end `.Fastq` files generated from NIPT sequencing, which will be used for synthetic data generation 
2. **Trim** single-end `.Fastq` data to the first 35 nucleotides  
3. **Align** the trimmed reads to the *hg19* reference index using **Bowtie2**  
4. **Convert** SAM to BAM, then **sort** and **index** the alignment files  
5. **Remove duplicate reads** using `samtools rmdup -s` (single-end mode)  
6. **Filter** for *uniquely mapped reads* and save them into the final `.unique` file  
7. For multiple samples, generate a `.unique` file for each sample, and **record the fetal fraction (FF) and GC content values** for each through a separate process  

## Generating Synthetic ACA (Trisomy) Dataset
- Refer to the implementation and annotations within the `make_synthetic_ACA.py` script

## Generating Synthetic SCA Dataset
- Refer to the implementation and annotations within the `make_synthetic_SCA.py` script 
