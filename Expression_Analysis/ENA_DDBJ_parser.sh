####I use this script to parse ENA and DDBJ text files into tabular files
###So that the format suitable for Expression Analysis script
###This script is highly customable 



#Getting female dataset from pseudoobscura
grep "pseudoobscura" SourceOfRNAData.tsv |grep -v 'miranda'|grep -i "female"|grep -o 'ftp[^;]\+gz'|cut -f 1 |grep -v '.*_2.*'
#Getting male dataset pseduoobscura
grep "pseudoobscura" SourceOfRNAData.tsv |grep -v 'miranda'|grep -i "male"|grep -v -i "female"|grep -o 'ftp[^;]\+gz'|cut -f 1 |grep -v '.*_2.*'\
#Long identifier
grep "pseudoobscura" SourceOfRNAData.tsv |grep -v 'miranda'|grep -i "male"|grep -v -i "female"|cut -f 15
#Short identifier
grep "pseudoobscura" SourceOfRNAData.tsv |grep -v 'miranda'|grep -i "male"|grep -v -i "female"|cut -f 
#Getting identifier from ftp 

cat url2| cut -f 2| cut -d'/' -f 6 |sed 's/.fastq.gz//g'
#Getting only fasta identifier from fasta file
grep  '>' dpsegene|tr -d '>'
#Print sequence length from fasta
awk '$0 !~ />/{print length($0)}' dpsegene

#Additional miranda dataset from embryo
grep 'miranda' ENA_ObscuraGroup_19thMay2017.tsv|grep -v 'pseudoobscura' | grep -i 'transcriptomic' |cut -f 14|grep -i '_M_'
grep 'miranda' ENA_ObscuraGroup_19thMay2017.tsv|grep -v 'pseudoobscura' | grep -i 'transcriptomic' |grep '_F_'|cut -f 14 > temp2
grep 'miranda' ENA_ObscuraGroup_19thMay2017.tsv|grep -v 'pseudoobscura' | grep -i 'transcriptomic' |grep '_F_'|grep -o 'ftp[^;]\+gz'|cut -f 1 |grep -v '.*_2.*' > temp3
grep 'miranda' ENA_ObscuraGroup_19thMay2017.tsv|grep -v 'pseudoobscura' | grep -i 'transcriptomic' |grep '_F_'|cut -f 5 > temp4
paste temp2 temp3 temp4 > embryofemalemiranda

#Getting useful info from ENA, 1: ID, 2: species,3: Experiment,20: link

grep 'pseudoobscura' SourceOfRNAData.tsv|cut -f 5,6,15,20>alldpse

#Removing and only left the fist read
cut -f 4 alldpse|sed 's/;.*//g'|less

#Removing all Ago genes
grep -v -i  ".*Ago.*" dpsegene|grep '>' -A 1 --no-group-separator |less

DRR055236 Illumina HiSeq 2000 paired end sequencing; mRNA sequencing of male abdomens of Drosophila pseudoobscura (replicate 1)	ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA004/DRA004463/DRX049920/DRR055236_1.fastq.bz2

#Delete blank line and shorten Drosophila
cut -f 2 alldpse|tr -d ' '|sed 's/Drosophila/D_/g'

#Use tr for removing spaces
cut -f 3 alldpse|tr ' ' '_'

#Use awk to filter only D. pseudoobscura
awk '$2 ~ /D_pseudoobscura/ {print $0}' tidyurk

#Removing space (no url) in column4
awk '$4 !~ /^$/  {print $0}' temp3>temp4

#Only getting the forward (_1) strand
grep -o 'ftp[^;]\+gz' temp4 |cut -f 1 |grep -v '.*_2.*'

#Extract column containing non reverse url download with tab separated value
awk -v OFS='\t' '{print $1,$2,$3,$5}' temp6 


makeblastdb -in scaffolds.fasta -out bifasago -dbtype 'nucl'
#URL DDBJ generate automatically
paste -d"," <(paste -d'\0' <(cat urlddjb |cut -f3) <(cat urlddjb |cut -f1)) <(cut -f2 urlddjb)|sed 's/,/\/DRR/g'>urlnig
awk -v OFS='' '{print $0,"_1.fastq.bz2"}' urlnig|less

#Getting dobs rpl32
tblastn -query rpl32obs -db CombinedTranscripomes -evalue 1e-40 -outfmt 6 |grep 'DOBS'


#Correcting url with adding 0 in front of DRX and DDR

sed 's/DRX/DRX0/g' urlnig2|sed 's/DRR/DRR0/g'




