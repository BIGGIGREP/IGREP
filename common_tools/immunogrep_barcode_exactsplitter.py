#@APP title: "IGREP Barcode Splitter",tags:"immunogrep_ngsprocess"

import immunogrep_makeappfunctions as makeapp
import appsoma_api
import os
import immunogrep_useful_immunogrep_functions as useful
import subprocess

#will a string of FASTA sequences into tab delimited file that can be used by fastx toolkit
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
	os.system("bash barcodesplit.bash '{0}' '{1}' '{2}' {3} {4} {5}".format(barcode_file_path,';'.join(input_file_path),outfile_path,file_ext,str(lines_per_seq),prefix))			
	
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
	
	

appsoma_api.resource_pull("https://www.appsoma.com/programs/get/cchrysostomou/default/simple_barcode_splitter_awk.bash",'barcodesplit.bash')

if __name__ == '__main__':
		
	file_list = [makeapp.LocateIGFile(divlocid='file-select_0',allowedFileFormats=['FASTQ','FASTA'])]
	
	while True:
		response = appsoma_api.communicate_await()	
		#add a file button was pressed
		if response.get('add-file'):
			file_list.append(makeapp.LocateIGFile(divlocid='file-select_'+str(response.get('add-file')),allowedFileFormats=['FASTQ','FASTA']))
		#delete a file button was pressed
		elif response.get('delete-file'):		
			file_list[response.get('delete-file')] = None
		#form was submitted (submit button was pressed)
		elif 'barcodes' in response:
			var_settings = response
			all_files_defined = True
			file_locations_fastq = []
			file_locations_fasta = []
			for i in file_list:
				if i:
					f = i.Get_Data()
					if f:
						if f['File_type'].upper() == 'FASTQ':
							file_locations_fastq.append(f['Folder_name_fullpath']+'/'+f['File_name'])
						elif f['File_type'].upper() == 'FASTA':
							file_locations_fasta.append(f['Folder_name_fullpath']+'/'+f['File_name'])
					else:
						all_files_defined = False
			if all_files_defined:
				break
		#any other event probably with respect to file_list was pressed
		else:
			for i in file_list:
				if i:
					i.Process_Events(response)
					
	concatenate_files = var_settings.pop('cat',None)
	
	appsoma_api.html_append("Running barcode splitting")
	if file_locations_fastq:	
		working_directory = os.path.dirname(file_locations_fastq[0])
		barcode_path = working_directory+'/barcodelist.txt'
		Write_FASTA_String_Barcodes_To_File(var_settings.pop('barcodes',None),working_directory+'/'+'barcodelist.txt',include_rc=var_settings['rcsearch'])	
		if concatenate_files:
			outfile_path = working_directory+'/'+var_settings['prefix']+os.path.basename(file_locations_fastq[0])+'.'
			result = Run_Barcode_Awk_Splitting(barcode_path, file_locations_fastq , outfile_path,'fastq',4)
		else:
			for file in file_locations_fastq:
				outfile_path = os.path.dirname(file)+'/'+var_settings['prefix']+os.path.basename(file)+'.'
				result = Run_Barcode_Awk_Splitting(barcode_path, [file] , outfile_path,'fastq',4)
				#g = subprocess.check_output("bash barcodesplit.bash '{0}' '{1}' '{2}' fastq 4".format(barcode_path,file,outfile_path),shell=True)
				
		
	if file_locations_fasta:	
		working_directory = os.path.dirname(file_locations_fastq[0])
		barcode_path = working_directory+'/barcodelist.txt'
		Write_FASTA_String_Barcodes_To_File(var_settings.pop('barcodes',None),working_directory+'/'+'barcodelist.txt',include_rc=var_settings['rcsearch'])
		if concatenate_files:
			outfile_path = working_directory+'/'+var_settings['prefix']+os.path.basename(file_locations_fastq[0])+'.'
			result = Run_Barcode_Awk_Splitting(barcode_path,file_locations_fastq,outfile_path,'fasta',2)
			#g = subprocess.check_output("bash barcodesplit.bash '{0}' '{1}' '{2}' fasta 2".format(barcode_path,';'.join(file_locations_fastq),outfile_path),shell=True)			
		else:
			for file in file_locations_fastq:
				outfile_path = os.path.dirname(file)+'/'+var_settings['prefix']+os.path.basename(file)+'.'
				result = Run_Barcode_Awk_Splitting(barcode_path,[file],outfile_path,'fasta',2)
				#g = subprocess.check_output("bash barcodesplit.bash '{0}' '{1}' '{2}' fasta 2".format(barcode_path,file,outfile_path),shell=True)
	
	appsoma_api.html_append("Analysis complete")			
	print result			
	
				
