read -p 'Location: ' loc
read -p 'Start: ' start
read -p 'Stop: ' stop
echo ">${loc} ${start} - ${stop}"  

if [ "$start" -lt "$stop" ]; then
grep -A 1 "$loc" CombinedAfflow.fas | tail -n 1 | cut -c $start-$stop 
else
grep -A 1 "$loc" CombinedAfflow.fas | tail -n 1 | cut -c $stop-$start | rev | tr  ATGC TACG >> blastr_${gene}.fasta
fi
