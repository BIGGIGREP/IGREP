import immunogrep_makeappfunctions as makeapp

import os
import immunogrep_useful_immunogrep_functions as useful
import subprocess

barcode_split_bash_program = os.path.join(os.path.relpath(__file__), "barcodesplit.bash")


# will a string of FASTA sequences into tab delimited file that can be used by fastx toolkit
def Write_FASTA_String_Barcodes_To_File(barcode_string,filename,include_rc = False):
	list_of_seqs = []
	barcode_array = barcode_string.split('>')[1:]
	for seq in barcode_array:		
		header_seq = seq.split('\n')
		list_of_seqs.append([header_seq[0],''.join(header_seq[1:]).replace(' ','')])
	
	with open(filename,'w') as out:
		for seq in list_of_seqs:		
			out.write(seq[0]+'\t'+seq[1]+'\n')
			if include_rc:
				out.write(seq[0]+'\t'+useful.Reverse_Complement(seq[1])+'\n')


def Run_Barcode_Awk_Splitting(barcode_file_path,input_file_path,outfile_path,file_ext,lines_per_seq,prefix=''):
	if not isinstance(input_file_path,list):
		input_file_path = [input_file_path]
	subprocess.call("bash '" + barcode_split_bash_program + "' '{0}' '{1}' '{2}' {3} {4} {5}".format(barcode_file_path,';'.join(input_file_path),outfile_path,file_ext,str(lines_per_seq),prefix), shell=True)	
	
	with open(outfile_path+prefix+"barcodesplitsummarytextfile.txt") as e:
		g = e.readlines()	
		
	result = {'barcodes':{}}
	print "Starting barcode splitting"
	for each_line in g:
		if each_line.strip().lower().startswith('total sequences:'):
			result['total'] = int(each_line.split('\t')[1].strip())
		elif each_line.strip().lower().startswith('unmatched:'):
			result['unmatched'] = int(each_line.split('\t')[2].strip())
		elif each_line.strip().lower().startswith('barcode file:'):
			each_line = each_line.split('\t')
			result['barcodes'][each_line[1].strip()] = int(each_line[2].strip())
	print "Barcode splitting complete"

	return result

