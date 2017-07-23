# RNAi genes duplication

All supporting scripts that I use to analyze genomic data in my Master's Dissertation. 
The detailed explanation for each script is given below

## Variant Calling Script

I used [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) , [samtools v1.4](http://samtools.sourceforge.net/) and [GATK v3.5](https://software.broadinstitute.org/gatk/) for identifying variants of 15-20 genes from short reads genomic data of *Drosophila pseudoobscura* (12 strains), *Drosophila miranda* (12 strains) and *Drosophila athabsca* (28 strains). At the end of script I used fasta_formatter from [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) to organize the fasta files according to the genes.


