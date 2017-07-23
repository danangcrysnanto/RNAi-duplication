#BLAST searching
#Iterating over all genes which  all has initial name as search

for gene in search*
do

cat $gene | fasta_formatter | tr -d "-" > tempseq #Fasta wrapper and removing connector
blastn -query tempseq -db CombinedAfflow -outfmt 6 -evalue 1e-20 > tempsearch_${gene} #Doing blast save in tempsearch

done



#Extracting sequence
for gene in tempsearch*
do

hits=`cat $gene|wc -l`
echo "total hits: ${hits}"

while read i
do
echo "############################################################################"
echo "$i"
loc=`echo "$i"| cut -f2`
start=`echo "$i" | cut -f9`
stop=`echo "$i" | cut -f10`
echo "location: ${loc}"
echo "position: ${start}-${stop}"

#Differentiate between match in + or - (need to reverse complement)
if [ "$start" -lt "$stop" ]; then
grep -A 1 "$loc" CombinedAfflow.fas | tail -n 1 | cut -c $start-$stop >> blast_${gene}.fasta
grep "$loc" CombinedAfflow.fas >> blast_${gene}.fasta
else
grep "$loc" CombinedAfflow.fas >> blastr_${gene}.fasta
grep -A 1 "$loc" CombinedAfflow.fas | tail -n 1 | cut -c $stop-$start | rev | tr  ATGC TACG >> blastr_${gene}.fasta
fi

done < ${gene}

echo "Finished!"
echo "############################################################################"

done

rm temp*










#Extracting the sequence

#Iterating over the result for location and start-stop index
#IFS=IFS=$'\n'


#for i in  tempsearch
#do
#sed -n 2p "$i"

#loc=`cat $i | awk '{print $1}`
#stop=`cat $i | awk '{print $9}'`
#start=`cat $i | awk '{print $10}'`
#echo "${start} - ${stop}"
#echo "$loc"
 

#getting the fasta read (location) and start and stop at the respective location
#loc= `cat $i | awk '{print $1}'`  
#stop= `cat $i | awk '{print $9}'` 
#start=`cat $i | awk '{print $10}'`
#echo "$loc"
#echo "$loc"
#echo "${start} - ${stop}"

#Print fasta information to a file
#grep "$loc" CombinedAfflow.fas >> Afflowei_${gene}.fasta
#Append the sequence that has been cut at the respective location
#grep -A 1 "$loc" CombinedAfflow.fas | tail -n 1 | cut -c $start-$stop >> Afflowei_${gene}.fasta
#done 

#done

#rm temp*
