# RNAi genes duplication

![Drosophila](http://obbard.bio.ed.ac.uk/photo_gallery/flies/Drosophila_lacicola.JPG) 

All supporting scripts that I use to analyze genomic data in my Master's Dissertation. 
The detailed explanation for each script is given below

## Variant Calling 

I used [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) , [samtools v1.4](http://samtools.sourceforge.net/) and [GATK v3.5](https://software.broadinstitute.org/gatk/) to collect variants of 15-20 genes from short reads population genomic data of *Drosophila pseudoobscura* (12 strains), *Drosophila miranda* (12 strains) [see McGaugh *et.al*, 2012](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001422) and *Drosophila athabsca* (28 strains) [see Miller *et al.* ,2017](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx134/3738284/Patterns-of-Genome-Wide-Diversity-and-Population). At the end of script I used fasta_formatter from [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) to organize the fasta files according to the genes. 


