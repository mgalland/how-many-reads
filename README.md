# how-many-reads
Analyses related to RNA-seq experimental design: how many reads per sample, power, etc.

## Number of genes detected per million reads 
__Question__: how many genes can I detect for a certain number of reads?

Here, I want to calculate for different species, different tissues, how many genes I can detect at different thresholds for a given number of reads.

Parameters to explore to see how it affects the final curve:
* Take 10% to 100% of the initial reads.
* Test different detection threshold: 10 reads, 100 reads ,1000 reads.

__Workflow steps:__ 
* Trim the intial reads (adapter, base quality)
* Create the different read subsets (10%, 20%, etc.)
* Pseudo-align to reference transcriptome (Kallisto)
* Output a matrix of gene counts for each sample (initial and subsetted).
* Plot a rarefaction curve of the number of detected genes (Y) as a function of the number of reads (X) for each sample.


## Quantification accuracy
__Question__: how the quantification accuracy affected by the number of reads being sequenced?





## Number of differentially expressed genes per million reads


# datasets
1. Xu 2015: Moneymaker stem trichomes. (Illumina HiSeq 125 nt.)
2. 