import appsoma_api
import os
import subprocess
import copy
from immunogrep_useful_immunogrep_functions import Reverse_Complement
from collections import defaultdict
import immunogrep_useful_immunogrep_functions as useful


#will a string of FASTA sequences into tab delimited file that can be used by fastx toolkit
def Write_FASTA_String_Barcodes_To_File(barcode_string,filename):
	list_of_seqs = []
	barcode_array = barcode_string.split('>')[1:]
	for seq in barcode_array:		
		header_seq = seq.split('\n')
		list_of_seqs.append([header_seq[0],''.join(header_seq[1:]).replace(' ','')])
	
	with open(filename,'w') as out:
		for seq in list_of_seqs:
			out.write(seq[0]+'\t'+seq[1]+'\n')

#CONVERTS lists of sequences into tab delimited file that can be used by fastx toolkit
#assumes 2D list, each row is unique/first index is unique sequence. within each row is 2d list containing sequence header, seqhence
#[[SEQUENCE_1,ACGT],[SEQUENCE2,AAAA]]
def Write_SEQ_list_Barcodes_To_File(barcode_list,filename):		
	with open(filename,'w') as out:
		for seq in barcode_list:
			out.write(seq[0]+'\t'+seq[1].strip().replace(' ','')+'\n')

def Run_FASTX_Barcode_Splitter(files,output_dir,settings={'orientation':'bol'},search_reverse_complement=True):	
	
	parameters = copy.deepcopy(settings)
	if 'orientation' not in parameters:
		raise Exception('"Orientation" is required in the parameters field')
	
	if not type(files) is list:
		files = [files]
	
	for i,each_file in enumerate(files):
		if each_file.endswith('.gz'):
			print "Unzipping file: {0}...".format(each_file)
			os.system("gunzip '{0}'".format(each_file))
			files[i] = each_file[:-3]			
			print "Unzipping complete"		
	
	suffix = parameters.pop('suffix') if 'suffix' in parameters else ''
	
	barcode_splitter_command = 'cat '+' '.join(files)+' | '
	
	if output_dir[-1] == '/':
		output_dir = output_dir[:-1]
	
		
	if 'prefix' in parameters and parameters['prefix'] != '':
		prefix = output_dir+'/'+parameters['prefix'] 
	else:
	 	prefix = output_dir+'/'	
	
	parameters.pop('prefix',None)
	
	additional_folders = os.path.dirname(prefix)
		
	if not os.path.isdir(additional_folders):
		os.system("mkdir '{0}'".format(additional_folders))
		
	orientation = parameters.pop('orientation',None)
	
	barcode_splitter_command += 'fastx_barcode_splitter.pl '
	for p in parameters:
		barcode_splitter_command+='--{0} {1} '.format(p,parameters[p])
	
	barcode_splitter_command +='--prefix '+prefix+ ' --suffix '+suffix+' --'+orientation
	
	output = useful.get_stdout(barcode_splitter_command).rstrip(' \n').split('\n')#output = subprocess.check_output(barcode_splitter_command,shell=True).rstrip(' \n').split('\n')

	if output[0].lower().startswith('error'):
		raise Exception("Error found in barcode split program: "+output[0])
	
	
	result = {'barcodes':defaultdict(int)}
	for line in output[1:-2]:
		line = line.split('\t')
		result['barcodes'][line[2]] = int(line[1])
	result['total'] = int(output[-1].split('\t')[1])
	result['unmatched'] = int(output[-2].split('\t')[1])
	
	
	if search_reverse_complement:		
		initial_file = []
		new_file = []
		map_barcode_to_file = {}
		with open(parameters['bcfile']) as file:
			lines=file.readlines()
			new_bcfile=open(settings['bcfile']+'rc','w')
			for l in lines:
				c = l.split('\t')
				initial_file.append(c[0].strip())
				new_file.append(c[0].strip()+'rev')
				map_barcode_to_file[c[0].strip()+'rev'] = prefix+c[0].strip()+suffix
				new_bcfile.write(c[0].strip()+'rev'+'\t'+Reverse_Complement(c[1].strip())+'\n')
			new_bcfile.close()
		
		os.system("cp '{0}' '{0}.temp'".format(prefix+'unmatched'+suffix))
		files = [prefix+'unmatched'+suffix+'.temp']
		parameters['bcfile'] +='rc'
		
		if orientation == 'eol':
			orientation='bol'
		elif orientation == 'bol':
			orientation= 'eol'
			
		barcode_splitter_command = 'cat "'+' '.join(files)+'" | '

		barcode_splitter_command += 'fastx_barcode_splitter.pl '
		for p in parameters:
			barcode_splitter_command+='--{0} {1} '.format(p,parameters[p])
		
		barcode_splitter_command +='--prefix '+prefix+' --suffix '+suffix+' --'+orientation
		
		#os.system(barcode_splitter_command)
		output = useful.get_stdout(barcode_splitter_command).rstrip(' \n').split('\n') #subprocess.check_output(barcode_splitter_command,shell=True).rstrip(' \n').split('\n')
		
		for i,line in enumerate(output[1:-2]):
			line = line.split('\t')
			result['barcodes'][map_barcode_to_file[line[0].strip()]] += int(line[1])
		
		result['unmatched'] = int(output[-2].split('\t')[1])
				
		cleanup_command = ''
		for i,each_bc_file in enumerate(initial_file):
			cleanup_command+= "mv '{0}{1}{3}' '{0}{1}{3}.temp';cat '{0}{1}{3}.temp' '{0}{2}{3}' > '{0}{1}{3}'; rm '{0}{1}{3}.temp';rm '{0}{2}{3}';".format(prefix,each_bc_file,new_file[i],suffix)				
		cleanup_command+="rm '{0}{1}'; ".format(prefix,'unmatched'+suffix+'.temp')
		os.system(cleanup_command)
		
		
	return result
	
	

def Run_Quality_Filter(files,output_dir,quality,percent,encoding='-Q33'):
		
	if not type(files) is list:
		files = [files]
		
	
	for i,each_file in enumerate(files):
		if each_file.endswith('.gz'):
			print "Unzipping file: {0}...".format(each_file)
			os.system("gunzip '{0}'".format(each_file))
			files[i] = each_file[:-3]			
			print "Unzipping complete"
				
	
	file_list = 'cat '+' '.join(['"'+f+'"' for f in files]) +' | '
	
	if output_dir[-1] == '/':
		outfile = output_dir+os.path.basename(files[0]).replace('.fastq','')+'.filtered.{0}.fastq'.format('q'+str(quality)+'p'+str(percent))
	else:
		outfile = output_dir+'/'+os.path.basename(files[0]).replace('.fastq','')+'.filtered.{0}.fastq'.format('q'+str(quality)+'p'+str(percent))

	print "Running filtering..."
	os.system('{3} fastq_quality_filter -v {4} -o "{0}" -q {1} -p {2}'.format(outfile,str(quality),str(percent),file_list,encoding))
	print "filtering complete"

	
	return outfile
