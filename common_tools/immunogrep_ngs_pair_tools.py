import os
import subprocess
import immunogrep_useful_functions as useful

from immunogrep_global_variables import program_folder
trimmomatic_location = os.path.join(program_folder, "trimmomatic-0.32.jar")
flash_location = os.path.join(program_folder, "flash")
pear_location = os.path.join(program_folder, "pear") 


def run_trimmomatic(files, output_directory=None, method='SE', phred=None, optional_parameters = {}):
	'''
		Wrapper function for running trimmomatic program within python
		Trimmomatic will remove low quality bases from the ends of NGS reads using an average quality score in a given window size
		
		Parameters
		----------
		files : string or list of strings
			List of input filenames (fastq or fastq.gz) for the MISEQ files. We either accept a single string or a list of two strings.
		working_directory : string, default none
			Pathname of desired output directory
		outfile : string, default empty string
			Desired filename name
		method : SE or PE, default 'SE'
			String representing whether to treat input files as single (SE) or paired-end files (PE)
		phred : integer, default None
			If None, then will rely on trimmomatic to guess the quality encoding. If a number, then will pass this value into the phred field.
		optional_parameters : dict, default empty parameters
			An optional dict of all parameters you would like to pass to trimmomatic 
			http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf		
	'''
	
	method = method.upper()
	if method not in ['SE', 'PE']:
		raise Exception('Incorrect value provided for parameter "method". Provided value: ' + method)
	
	if not isinstance(files, list):
		files = [files]
	if len(files) > 2:
		raise Exception(str(len(files)) + 'total files have been passed to function. We only except 1 or 2 filepaths representing the R1/R2 reads')
	
	for i,f in enumerate(files):
		if f.endswith('.gz'):
			print('Unzipping: ', f)
			files[i] = useful.gunzip_python(f)		
	
	output_directory = useful.get_parent_dir(files[0]) if not output_directory else os.path.abspath(output_directory)		
		
	
	return_file_names = []	
	command_loops = []
	
	if method == 'SE':
		# Trim each file at a time
		for f in files:
			input_file_names = []
			output_file_names = []
			input_file_names.append('"'+f+'"')			
			out = f[:-6] if f.endswith('.fastq') else f
			output_file_names.extend(['"'+out + '.trimmed.fastq"'])
			return_file_names.append(out+'.trimmed.fastq')
			command_loops.append([input_file_names, output_file_names])
	else:
		input_file_names = []
		output_file_names = []
		# trim all files simultaneously		
		for f in files:
			input_file_names.append('"'+f+'"')			
			out = f[:-6] if f.endswith('.fastq') else f
			output_file_names.extend(['"'+out + '.trimmed.fastq"', '"'+out + '.trimmed.unpaired.fastq"'])
			return_file_names.append(out+'.trimmed.fastq')
		command_loops.append([input_file_names, output_file_names])
	phred_var = '-phred'+str(phred) if phred else ''	
	
	# We should change the java folder to recognize /usr/local/bin...	
	for loops in command_loops:
		inputs = loops[0]
		outputs = loops[1]
		trim_command = 'java -jar {5} {0} {4} -threads 2 {1} {2} {3}'.format(method, ' '.join(inputs), ' '.join(outputs), ' '.join([key + ':' + str(value) for key, value in optional_parameters.iteritems()]), phred_var, trimmomatic_location)		
		worked = subprocess.call(trim_command, shell=True)
		if worked > 0:
			raise Exception('Trimmomatic failed')
	return return_file_names
	

def run_flash(r1file, r2file, working_directory, outfile='', parameters={}, suffix=''):
	r1_path = useful.get_parent_dir(r1file)  # '/'.join(r1file.split('/')[:-1])
	
	r2_path = useful.get_parent_dir(r2file)  # '/'.join(r2file.split('/')[:-1])
	
	if not parameters:
		print "PARAMETERS NOT PASSED INTO FLASH PROGRAM. USING DEFAULT IGSEQ PARAMETERS: R = 300, F = 400"
		parameters = {'r':300,'f':400}
	
	if r1file.endswith('.gz'):
		print "Unzipping R1 File.."
		
		#os.system("gunzip '{0}'".format(r1file))
		#r1file = r1file.replace('.gz','')
		r1file = useful.gunzip_python(r1file)
	
	if r2file.endswith('.gz'):		
		print "Unzipping R2 File.."
		#os.system("gunzip '{0}'".format(r2file))		
		#r2file = r2file.replace('.gz','')
		r2file = useful.gunzip_python(r2file)
		
	working_directory = os.path.abspath(working_directory)
	if r1_path != working_directory:
		os.rename(r1file,os.path.join(working_directory, os.path.basename(r1file)))
		#os.system("mv '{0}' '{1}/'".format(r1file,working_directory) )
	if r2_path != working_directory:	
		os.rename(r2file,os.path.join(working_directory, os.path.basename(r2file)))
		#os.system("mv '{0}' '{1}/'".format(r2file,working_directory) )
		
	if outfile == '':		
		outfile=os.path.basename(r1file)		
		if '_R1' in outfile:
			r_pos=outfile.index('_R1')
			outfile=outfile[:r_pos]
		elif '_R2' in outfile:
			r_pos=outfile.index('_R2')
			outfile=outfile[:r_pos]
		else:				
			outfile =  os.path.commonprefix([os.path.basename(r1file),os.path.basename(r2file)]) # useful.removeFileExtension(os.path.basename(r1file))
		#if len(outfile)<5:
		#try:
		#	r1_pos = outfile.index('_R1')
		#	outfile = outfile[:r1_pos]
		#except:
		#	outfile = outfile.replace('R1','')
	else:		
		outfile = os.path.basename(outfile)
	outfile = outfile.replace('.fastq','').replace('.fasta','')
	outfile+='.flashed'+suffix
		
			
	if os.path.isfile(os.path.join(working_directory, outfile)):# in resulting_files:		
		print('WARNING: FILE {0} ALREADY PRESENT IN FOLDER. FILE WILL BE OVERWRITTEN'.format(working_directory+'/'+outfile))
							
	r1file = os.path.join(working_directory,os.path.basename(r1file)) # working_directory+'/'+os.path.basename(r1file)
	r2file = os.path.join(working_directory,os.path.basename(r2file)) # working_directory+'/'+os.path.basename(r2file)

	flash_command = "{2} {0} {1}".format(r1file, r2file, flash_location)
	
	parameters['o'] = outfile
	parameters['d'] = working_directory

	for p,val in parameters.iteritems():
		flash_command+=' -{0} {1}'.format(p,str(val))
	
	flash_command+=' -q'#run on quiet command
	# os.system(flash_command)
	worked = subprocess.call(flash_command, shell=True)
	if worked > 0:
		raise Exception('Flash failed')
	os.rename(os.path.join(working_directory, outfile+'.extendedFrags.fastq'), os.path.join(working_directory, outfile))
	
	try:
		read_count_r1_file = useful.file_line_count(r1file)
	except Exception as e:
		read_count_r1_file = 1
		print("Could not get number of lines in read file: " + str(e))
	
	try:
		read_count_flashed_file = useful.file_line_count(os.path.join(working_directory, outfile))
	except Exception as e:
		read_count_flashed_file = 1
		print("Could not get number of lines in outfile read file: "+ str(e))
	resulting_counts=(
		os.path.join(working_directory, outfile),
		read_count_flashed_file/4,
		read_count_r1_file/4,
		float(100)*(read_count_flashed_file/float(read_count_r1_file))
	)
	
	return resulting_counts
	
	
def run_pear(r1file, r2file, working_directory, outfile='', parameters={}, suffix='', num_threads=1, memory='1G'):
	r1_path = useful.get_parent_dir(r1file)
	r2_path = useful.get_parent_dir(r2file)

	if r1file.endswith('.gz'):
		print("Unzipping R1 File..")
		r1file = useful.gunzip_python(r1file)
	
	if r2file.endswith('.gz'):		
		print("Unzipping R2 File..")
		r2file = useful.gunzip_python(r2file)
				
	working_directory = os.path.abspath(working_directory)
	if r1_path != working_directory:
		os.rename(r1file, os.path.join(working_directory, os.path.basename(r1file)))		
	if r2_path != working_directory:	
		os.rename(r2file, os.path.join(working_directory, os.path.basename(r2file)))		
		
	if outfile == '':		
		outfile = os.path.basename(r1file)		
		if '_R1' in outfile:
			r_pos = outfile.index('_R1')
			outfile = outfile[:r_pos]
		elif '_R2' in outfile:
			r_pos = outfile.index('_R2')
			outfile = outfile[:r_pos]
		else:				
			outfile = os.path.commonprefix([os.path.basename(r1file),os.path.basename(r2file)])
	else:		
		outfile = os.path.basename(outfile)
	
	outfile = outfile.replace('.fastq', '').replace('.fasta', '')
	
	outfile = os.path.join(working_directory, outfile)
	if os.path.isfile(os.path.join(working_directory, outfile)):# in resulting_files:
		print('WARNING: FILE {0} ALREADY PRESENT IN FOLDER. FILE WILL BE OVERWRITTEN'.format(working_directory+'/'+outfile))
							
	r1file = os.path.join(working_directory, os.path.basename(r1file))
	r2file = os.path.join(working_directory, os.path.basename(r2file))

	pear_command = "{2} -f {0} -r {1}".format(r1file, r2file, pear_location)
	

	parameters['o'] = outfile
	parameters['y'] = memory
	parameters['j'] = num_threads
	
	for p, val in parameters.iteritems():
		pear_command += ' -{0} {1}'.format(p, str(val))
			
	worked = subprocess.call(pear_command, shell=True)
	
	if worked > 0:
		raise Exception('Error in pear program')
	
	try:
		read_count_r1_file = useful.file_line_count(r1file)
	except Exception as e:
		read_count_r1_file = 1
		print("Could not get number of lines in read file: " + str(e))
	
	try:
		read_count_flashed_file =useful.file_line_count(outfile + '.assembled.fastq')
	except Exception as e:
		read_count_flashed_file = 1
		print("Could not get number of lines in outfile read file: " + str(e))

	resulting_counts = (
		outfile + '.assembled.fastq',
		read_count_flashed_file / 4,
		read_count_r1_file / 4,
		float(100) * (read_count_flashed_file / float(read_count_r1_file))
	)
	
	return resulting_counts
