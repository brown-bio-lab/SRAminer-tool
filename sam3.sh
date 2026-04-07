#!/bin/sh

input=keyword_IDlist.txt
i=1
touch sam_duplicates.fa
while IFS= read -r line
do
qstart=$(awk 'FNR == '$i' {print $10}' keyword.txt)
qend=$(awk 'FNR == '$i' {print $11}' keyword.txt)
if [ $qstart -gt  $qend ]
then
	temp=$qstart
	qstart=$qend
	qend=$temp
fi
samtools faidx keywords_contigs.fa $line:$qstart-$qend >> sam_duplicates.fa
i=$i+1
done < "$input"



