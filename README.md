# RNAi genes duplications in *Drosophila*

![Drosophila](http://obbard.bio.ed.ac.uk/photo_gallery/flies/Drosophila_lacicola.JPG) 

All supporting scripts that I use to analyze genomic data in my Master's Dissertation. 
The detailed explanations for each script is given below

## Variant Calling 

* I used [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) , [samtools v1.4](http://samtools.sourceforge.net/) and [GATK v3.5](https://software.broadinstitute.org/gatk/) to collect variants of 15-20 genes from short reads population genomic data of *Drosophila pseudoobscura* (12 strains), *Drosophila miranda* (12 strains) [see McGaugh *et.al*, 2012](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001422) and *Drosophila athabsca* (28 strains) [see Miller *et al.* ,2017](https://academic.oup.com/mbe/article/doi/10.1093/molbev/msx134/3738284/Patterns-of-Genome-Wide-Diversity-and-Population). 
Finally, I used fasta_formatter from [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) to organize the fasta files according to the genes. The script is available here: [varcall.sh](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Population_Genetic/varcall.sh)
* Any linkage information in heterozygous individual is recovered using [FastPhase](http://scheet.org/software.html). I created a R script to parse fasta files into format suitable for FastPhase, running FastPhase and reconvert the output back into Fasta format. The script is available here: [FastPhaseIntercovert.R](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Population_Genetic/FastaPhaseInterconvert.r). The original script credits to [Dr Darren Obbard](http://obbard.bio.ed.ac.uk/) and I modified the script so that it compatible with current version of FastPhase.
* To gain information on the effect of weakly deleterious mutations, I created a R script which removes variants with MAF less than 0.15. The scripts run by parsing fasta files into matrix, remove variants with MAF less than 0.15 and replace those variants with major variants at corresponding site and finally, convert back the matrix into Fasta format. The script is available here: [MAFRemover.R](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Population_Genetic/MAF_remover.R)

![VarCall](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Population_Genetic/varcall_workflow.png)

## Expression Analysis
* I collated transcriptome datasets from [ENA](http://www.ebi.ac.uk/ena) and [DDBJ](ddbj.nig.ac.jp) and used [custom parser](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Expression_Analysis/ENA_DDBJ_parser.sh) to process the text files into format that I set for expression script
* The [expression script](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Expression_Analysis/mapping_expression.sh) runs by mapping RNA-seq reads with reference transcriptome [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and uses [samtools v1.4](http://samtools.sourceforge.net/) to retain and count only the mapped reads for each genes (by using idxstats samtools)
* The rest of the scripts are used for [normalization](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Expression_Analysis/Normalization_Gene.R), [transform and plotting data](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Expression_Analysis/Normalization_Gene.R) using [ggverse packages](http://tidyverse.org/) and [mixed model analysis](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Expression_Analysis/mcmcGLMM.sh) using [mcmcGLMM](https://cran.r-project.org/web/packages/MCMCglmm/index.html)

## Phylogenetic Analysis

* I collected gene sequences using tBLASTN of local database. To automate the process, I write a script which will take sequences query, doing BLAST search and outputing fasta from BLAST hits output (and reverse complementing when the hits map in reverse direction). The script is available here: [Reverse-ComplementBLASThits.sh](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Phylogenetic_Analysis/Reverse-ComplementBLASThits.bash) and I also created script to manually select region of interest from sequence database (which I find quite handy when inspecting sequence and looking for the surrounding sequence). This script is available here: [BLAST search](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Phylogenetic_Analysis/BLAST_search.sh)
* For the species that only has genomic reads, I used a targetted assembly approach to assemble only reads that match with known RNAi protein. I identify the reads using [Diamond](https://github.com/bbuchfink/diamond) and then use single cell/bacterial assembly program, [Spades](Spades.bioinf.spbau.ru). The script is available here: [Diamond_Spades.sh](https://github.com/danangcrysnanto/RNAi-duplication/blob/master/Phylogenetic_Analysis/Diamond_Spades.sh).
