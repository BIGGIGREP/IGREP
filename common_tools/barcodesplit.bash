barcodefile=$1
sequencefile=$2
output_file=$3
output_extension=$4
lineloop=$5
prefix=$6

allsequencefiles=${sequencefile/;/ }

awk -v lines_per_seq="$lineloop" -v output_file="$output_file" -v prefix="$prefix" -v ext="$output_extension" 'BEGIN{seq_count=0;barcodes_count["nomatch"]=0;print "Searching files"}	 
	 NR==FNR{barcodes[toupper($2)]=$1;barcodes_count[$1]=0; next} 
	 {	 			 			
		if((FNR-1)%lines_per_seq==1){			
			seq_count+=1;
			seq[1]=toupper($1);
			barcode_file="nomatch"
			for (b in barcodes){							
				if( index(seq[1],b) > 0 ){					
					barcode_file = barcodes[b]
					break;
				}				
			}
		}
		else{
			seq[(FNR-1)%lines_per_seq]=$1;			
		}
		if(FNR%lines_per_seq==0){			
			barcodes_count[barcode_file]+=1;			
			for (i in seq){
				print seq[i]> output_file""prefix""barcode_file"."ext;
			}
			if(seq_count%100000==0){
				print seq_count;
			}
		}
		
	 }
	 END{
	 	unmatched=0;
	 	print "Total Sequences:","\t",seq_count>output_file""prefix"barcodesplitsummarytextfile.txt";
		for(barcode in barcodes_count){			
			if(barcode=="nomatch"){
				unmatched=barcodes_count[barcode]
			}
			else{
				print "Barcode file:","\t",output_file""prefix""barcode"."ext,"\t",barcodes_count[barcode]>output_file""prefix"barcodesplitsummarytextfile.txt";
			}
		}
		print "Unmatched:","\t",output_file""prefix"nomatch."ext,"\t",barcodes_count["nomatch"]>output_file""prefix"barcodesplitsummarytextfile.txt";
	 }
	 ' $barcodefile $allsequencefiles
	 
	 

