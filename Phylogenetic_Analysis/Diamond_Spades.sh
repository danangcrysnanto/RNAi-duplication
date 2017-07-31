
###I use this script for targetted assembly of Drosophila algonquin reads 
###I assembled only reads that match known RNAi genes protein sequence (identified using Diamond)
###And then used the scaffold from Spades assembly as reference in local BLAST

#make diamond database
 diamond makedb --in alldpseprot --db refdpse.dmnd
 
 #search diamond (blastx)
 diamond blastx --db refdpse.dmnd --daa algonhits.daa  --query   algonquin.fastq.gz --max-target-seqs 1 --evalue 1 --threads 10
 
 #look output
 diamond view  --daa  algonhits.daa > Blast-like_output.tsv
 (two strands)
 diamond view  --daa  dathf.daa > Blast-like_output.tsv
  diamond view  --daa  dathr.daa >> Blast-like_output.tsv
  
 #Get read names of the hits
 cut -f 1 Blast-like_output.tsv | sort | uniq > read_names.txt
 
 #Get the sequence from fastq file using zgrep
 zgrep -F -f read_names.txt -A 3 --no-group-separator  algonquin.fastq.gz| gzip > RNAiGene-like.fastq.gz
 
 #de-novo built with spades
 /media/DataDrive5/danang_working/SPAdes-3.9.0-Linux/bin/spades.py -s RNAiGene-like.fastq.gz -o RNAi-algon
 /media/DataDrive5/danang_working/SPAdes-3.9.0-Linux/bin/spades.py --s1 forward_RNAiGene-like.fastq.gz --s2 reverse_RNAiGene-like.fastq.gz -o RNAi-dath
 
 #convert to two lines fasta instead of multiple fasta lines
 cat scaffolds.fasta |fasta_formatter > algonref.fasta
 
 #Use the assembled sequences as  reference
 makeblastdb -in scaffolds.fasta -out algonref -dbtype 'nucl'
 
 #Doing tblastn
tblastn -query dathprot -db algonref -outfmt 6 -evalue 1e-40 > temp5
