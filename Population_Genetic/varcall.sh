while read i
do
echo "############################################################################"
echo "$i"
#Url-parser
strain=`echo "$i"| cut -f1`
acc=`echo "$i" | cut -f2`
url=`echo "$i" | cut -f3`
echo "strain: ${strain}"
echo "acc: $acc"
echo "url: $urlf"

#get the file
wget $url -O temp_${acc}.fastq.gz


#Bowtie2-mapping
bowtie2  --very-fast -x dpseref  -q -U temp_${acc}.fastq.gz -p 10 |#Mapping reads to references
samtools view -F 4 -L posgen.bed -hb | #taking only position of interest and take only the mapped reads
samtools sort -@ 6 > ${acc}_${strain}.bam #sorted
samtools index ${acc}_${strain}.bam


#Mark Duplicates (Picard Tools)
java -jar MarkDuplicates.jar INPUT=${acc}_${strain}.bam OUTPUT=${acc}_${strain}_dedup.bam  METRICS_FILE=${acc}_${strain}.matrix

#Indexing dedup bam
samtools index ${acc}_${strain}_dedup.bam

#Add Read Group Information (needed for GATK indel realignment)
java -jar AddOrReplaceReadGroups.jar I=${acc}_${strain}_dedup.bam O=${acc}_${strain}_final.bam RGID=${acc}_${strain}  RGLB=${acc}_${strain}  RGPL=illumina RGPU= ${acc}_${strain} RGSM= ${acc}_${strain}

#Indexing bam with read group 
samtools index ${acc}_${strain}_final.bam

#Target interval for indel realignment (GATK)
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R  dpse-all-chromosome-r3.04.fasta -I ${acc}_${strain}_final.bam -o ${acc}_${strain}.intervals

#Indel realigner (GATK)
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R  dpse-all-chromosome-r3.04.fasta -I ${acc}_${strain}_final.bam -targetIntervals ${acc}_${strain}.intervals -o ${acc}_${strain}_realigned.bam

#Indexing realigned indel bam
samtools index ${acc}_${strain}_realigned.bam

#Samtools mpileup
samtools mpileup -f  dpse-all-chromosome-r3.04.fasta --position posgen.bed -u  --skip-indels  ${acc}_${strain}_realigned.bam|bcftools call -m  > ${acc}_${strain}_realigned.vcf

#Select only SNP (GATK)
java -jar GenomeAnalysisTK.jar -T SelectVariants -R dpse-all-chromosome-r3.04.fasta -V ${acc}_${strain}_realigned.vcf -o  ${acc}_${strain}_snp.vcf --selectTypeToInclude SNP

#Low quality variant filtering (to be masked with NNNN in fasta alternate, GATK)
java -jar GenomeAnalysisTK.jar -T SelectVariants -R dpse-all-chromosome-r3.04.fasta -V ${acc}_${strain}_realigned.vcf -o  ${acc}_${strain}_masked.vcf -select "DP<5" 

#Constructing the alternate reference using FastaAlternate (GATK)
#Emit Heterezoygote site with IUPAC code (Y,R,K)
java -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R dpse-all-chromosome-r3.04.fasta -o ${acc}_${strain}_realigned.fas -L posgen.bed -V ${acc}_${strain}_snp.vcf --snpmask ${acc}_${strain}_masked.vcf --use_IUPAC_sample ${acc}_${strain}

rm temp*

done <urllist


###Converting multiple line fasta into single line sequence fasta 

for file in SRR*fas
do
name=`echo $file|cut -d'_' -f2,3`
cat $file|fasta_formatter > SRR_tidy_${name}
done


#####Getting fasta with strain identifier Dpse####### 
###I havel already know that there are 12 genes
for i in $(seq 12)
do
echo $i
for fasta in  SRR_tidy*
do
strain=`echo "$fasta"|cut -d'_' -f3| sed "s/.fas//g"`
echo ">Dpse_"$strain >> Dpse_gene_${i}.fas
grep -A 1 ">${i}\b" $fasta | tail -n 1 >> Dpse_gene_${i}.fas
done
done   

