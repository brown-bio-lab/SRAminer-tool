#!/bin/sh

input=keyword_IDlist.txt
i=1
touch sam_duplicates.fasta
while IFS= read -r line
do
qstart=$(awk 'FNR == '$i' {print $10}' keywords.txt)
qend=$(awk 'FNR == '$i' {print $11}' keywords.txt)
if [ $qstart -gt  $qend ]
then
        temp=$qstart
        qstart=$qend
        qend=$temp
fi
samtools faidx keywords_contigs.fa $line:$qstart-$qend >> sam_duplicates.fasta
i=$i+1
done < "$input"

