while read line
do
   url=`echo "$line"| cut -f 4`
   species=`echo "$line"| cut -f2`
   id2=`echo "$line" | cut -f 3`
   name=`echo "$line" | cut -f 1`
   if [ "$species" == "Drosophila miranda" ] ; then
   refbow='refmir'
   elif [ "$species" == "Drosophila obscura" ] ; then
   refbow='refobs'
   elif [ "$species" == "Drosophila pseudoobscura" ] ; then
   refbow='refobs'
   fi
   echo $url
   echo $id2
   echo $refbow
   echo 'Result from:' $id2 >> dpsemale.tsv
   echo 'SSR id:' $name >> dpsemale.tsv
   wget $url -O tempmapped.fastq.gz
   bowtie2  --very-sensitive  -x dpseref  -q -U  tempmapped.fastq.bz2 -p 10 | samtools view -F 4 -u -bS | samtools sort -@ 6 > $name.bam
   samtools index $name.bam
   samtools idxstats $name.bam | grep -v '*'  >> dpsemale.tsv
   rm temp*

done<urlnig
