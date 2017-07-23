
for gene in search*
do

cat $gene | fasta_formatter | tr -d "-" > tempseq #Fasta wrapper and removing connector
blastn -query tempseq -db CombinedAfflow -outfmt 6 -evalue 1e-20 > tempsearch_${gene}

done

for i in tempsearch
do
 
#getting the fasta read (location) and start and stop at the respective location
loc=`cat $i |cut -f2`  
start=`cat $i|cut -f9` 
stop=`cat $i|cut -f10`

#Print fasta information to a file
grep "$loc" CombinedAfflow >> Afflowei_${gene}.fasta
#Append the sequence that has been cut at the respective location
grep -A 1 "$loc" CombinedAfflow | head -1 | cut -c $start-$stop >> Afflowei_${gene}.fasta
done

done

rm temp*



grep -A 1  ">Dlow_4_group3" CombinedAfflow.fas | tail -n 1| cut -c $s1-$s2 >> tempres
