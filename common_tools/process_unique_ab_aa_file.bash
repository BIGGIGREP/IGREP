#!/bin/bash

#assumes tab delimited file , also assumes that the last column in the file is counts in a specific experiment
#step1 => sort by the first column  (the antibody sequence)
#step2 => use an awk script to create counts of unique sesquences for every experiment and total counts 
	#a=>read the last column number. this corresponds to which experiment this file came from
	#b=>for every line read, only include fileds from the first column to the second to last column (i=2;i<NF). i=1 is accounted for right before, when ab_line=$1
	#c=>store the counts for that sequence in every experiment using the number from 'a'
	#d=> once a new antibody is observed ($1!=ab) then print the unique antibody seuqen ce to file and all of the counts obvserved in all experiemnts 
#step 3=> make a final sort by the 'total_counts' column (which will always be equal to the length of header_row(total_count_column)
#step 4=> append the header row to the file 	

inputfile=$1
outputfile=$2
headerfile=$3
numberexps=$4
totalcount=$5
dirname=$(dirname "${inputfile}")
sort -t$'\t' -k1 "$inputfile" | 
awk -v expnum="$numberexps" 'BEGIN{FS="\t";OFS="\t"}
	FNR==1{count=0;ab=$1;ab_line=$1;for(i=2;i<NF;i++)ab_line=ab_line""OFS""$i;for(i=0;i<expnum;i++)exp_list[i]=0}
	{exp_num=$NF-1}
	$1!=ab{line=ab_line"\t"count;for(i=0;i<expnum;i++)line=line""OFS""exp_list[i];print line;count=0;for(i=0;i<expnum;i++)exp_list[i]=0;ab=$1;ab_line=$1;for(i=2;i<NF;i++)ab_line=ab_line""OFS""$i}
	{count+=1;exp_list[exp_num]+=1}
	END{line=ab_line"\t"count;for(i=0;i<expnum;i++)line=line""OFS""exp_list[i];print line}' |
sort -T "$dirname" -t$'\t' -k"$totalcount","$totalcount"nr -k1,1 |cat "$headerfile" - > "$outputfile";
rm "$3";
rm "$1";