
"""
	Set of functions for running igfft, parsing results, CDR3 calls, and ISOTYPE calls
"""

from Bio.Seq import Seq
import sys
import re

from Bio.Alphabet import generic_dna
import time
from time import gmtime, strftime
import json
import os

import subprocess
import traceback
from pprint import pprint
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict
import collections
import copy
from multiprocessing import Process
from multiprocessing import Queue

# Useful functions
import immunogrep_useful_functions as useful
from immunogrep_useful_functions import purge
# For global variables
from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import fasta_file_delimiter
from immunogrep_global_variables import descriptor_symbol  # this is the symbol we will use to designate 'comment' lines in the file 
from immunogrep_global_variables import translation_var  # this key will signify the translation/translator key 
from immunogrep_global_variables import filetype_var
from immunogrep_global_variables import program_folder
# For reading files
import immunogrep_read_file as readwrite
# use this for querying germlines sequences
import immunogrep_query_germline_functions as query_germlines
# For CDR3
import immunogrep_cdr3_tools as CDR3tools
# For isotyping
try:	
	import immunogrep_fft_align_tools as fftaligner
	from immunogrep_isotype_fft import defaultBarcodes
	isotypeworking = True
except:	
	isotypeworking = False

# Keep an interal folder of database files. Within this folder are subfolders by species and loci.
# A text file for each locus set is stored as a tab delimited file
# If this folder is not present, then a user MUST provide the location of their own files OR change the location of the database folder
global databaseFolder
databaseFolder = os.path.join(os.path.dirname(os.path.relpath(__file__)), "igfft_germline_files")

# Keep the location of the binary as a string
igfft_location = os.path.join(program_folder, 'igfft')


def modify_germline_DB_path(newpath):
	"""
		If the default database folder does not exist (subfolder within the module script location) then
		this function will change the default path for finding germline sequences defined by the databaseFolder
	"""
	if not os.path.isdir(newpath):
		raise Exception("The provided db path does not exist! " + newpath)
	global databaseFolder
	databaseFolder = newpath


def igfft_multiprocess(input_file, species, locus, file_type=None, output_dir=None, germline_source="default", num_processes=1, igfft_settings={}, parsing_settings={}, custom_germline_file_loc={}, germline_db_query={}, return_igrep_doc_line=False):
	"""
		Wrapper function for running igfft on multiple processes. This will take input file, split into multiple files, and then run 
		igfft and igfft parsing function on each file. igfft parsing function will always check for CDR3 and Isotyping using default settings.

		input_file can only be a FASTA or FASTQ file

		Parameters
		---------
		input_file : string, 
			Filepath of the FASTQ/FASTA dataset		
		species : list of strings
			Species databases used to run igfft and identify CDR3
		locus : list of strings
			Loci databases used to run igfft and identify CDR3
		file_type: string, default=None
			Allowed file types are "FASTQ" or "FASTA"
			If None, then program will attempt to guess filetype
		output_dir : string, default None
			Directory of the desired output file. If None, then output_dir will be equal to the directory of the input_file
		germline_source : string ("default", "db", or "custom")
			Settings for where we should look for germline files.
			..important::
			If germline_source is "custom" then the user must have defined the location of the v, d, and j files in the parameter custom_germline_file_loc
			For "default" and "db" we will use the species list and locus list parameters to extract germlines.
		num_processes : int, default 1
			Number of files to create and run in parallel
		igfft_settings : dict
			Additional alignment settings used for running igfft. See run_igfft function below. 
		parsing_settings : dict
			Additional alignment settings used for parsing igfft. I.e. settings for isotyping, idnentifying CDR3, handling insertions. See parse_alignment_file function below.
		custom_germline_file_loc : dict, default {}
			Only used when germline_source is "custom"
			Allowed keys in the dict are "v", "d", and "j"
			Values are list of file paths for all germlines to use as "v" set, "j" set and "d" set
		germline_db_query : dict, default {}
			Only used when germline_source is "db"
			Dict defines fields you want to use in a query that are in addition to the "SPECIES" and "LOCUS" settings
		return_igrep_doc_line : boolean, default False
			When True, then a comment line is added to the top of the annotated file after parsing. This line defines how to add select fields from the file into the database.
	
		Returns
		-------
		List of two strings:
			Element 1: location of the final parsed file containing genes, fr, cdr, cdr3, and isotype
			Element 2: location of the file created by igfft containing genes, fr, cdr

		Examples
		--------
		1) Standard usage: only define the  input file, species, and locus (will run as a single process, and will guess the input file type)
		>>> igfft_multiprocess(file.fastq, species="Homo sapiens", locus="IGH")
		2) Define file type, define multiple loci, run as 6 parallel processes
		>>> igfft_multiprocess(file.fasta, species="Homo sapiens", locus=["IGH","IGK","IGL"], num_processes=6, file_type="FASTA")
		3) Define multiple loci, run as 6 parallel processes, change a parameter in parsing settings (always remove insertins from alignment)
		>>> igfft_multiprocess(file.fasta, species="Homo sapiens", locus=["IGH","IGK","IGL"], num_processes=6, parsing_settings={'remove_insertions': 1})
		4) Define multiple loci, run as 6 parallel processes, change a paramter in annotation settings (when annotating only return first top hit), change a parameter in parsing settings (always remove insertins from alignment)
		>>> igfft_multiprocess(file.fasta, species="Homo sapiens", locus=["IGH","IGK","IGL"], num_processes=6, igfft_settings={'num_hits': 1}, parsing_settings={'remove_insertions': 1})
		5) Annotate V germlines of a FASTQ file using a custom v germline set
		>>> igfft_multiprocess(file.fastq, germline_source="custom", custom_germline_file_loc={'v': '~/common_tools/igfft_germline_files/homosapiens/v/homosapiens_igh_v.txt'})	
	"""
	if not file_type:
		# Attempt to guess filetype
		tempfile = readwrite.immunogrepFile(input_file)
		file_type = tempfile.getFiletype()
		if file_type is None:
			raise Exception("Could not predict the filetype for this file. Please manually define the file_type using file_type parameter in function")
	
	if file_type not in ['FASTQ', 'FASTA']:
		raise Exception("The provided file, {0}, does not appear to be a FASTQ or FASTA file. Predicted type: {1}. Please manually define if the file is a FASTQ or FASTA using the file_type file_type parameter".format(input_file, file_type))

	if not output_dir:
		output_dir = os.path.dirname(input_file)
	else:
		if not os.path.isdir(output_dir):
			os.path.makedirs(output_dir)
	basename_file = useful.removeFileExtension(os.path.basename(input_file))
	
	# Files created by program will be placed here
	outfile_queue = Queue()
	# this should be the final name provided to the alignment file 	
	of1_prefix = os.path.join(output_dir, basename_file)

	if 'isotype' not in parsing_settings:
		# Always run isotyping
		parsing_settings['isotype'] = {}
	if 'cdr3_search' not in parsing_settings:
		# Always run CDR3 search
		parsing_settings['cdr3_search'] = {'SPECIES': species, 'LOCUS': locus}
	if germline_source == "default":
		germline_location_settings = {'SPECIES': species, 'LOCUS': locus}
	elif germline_source == "db":
		germline_location_settings = {'SPECIES': species, 'LOCUS': locus}
		germline_location_settings.update(germline_db_query)
	elif germline_source == "custom":
		germline_location_settings = custom_germline_file_loc
		
	def calligfft(f1, output_file_prefix, index=0):
		"""
			This inside function actually runs the programs
			index => current file being analyzed/called by the function
			output_file_prefix => this defines how the result files wil be named WITH OUT EXTENSION (we always add the .alignment and .annoation)
			queue => this uses multithreading queues to keep track of what files were genreated in each process/thread
			index => this the number of the thread/process being run (makes it easy for sorting the results so that they mach the order they were input)
		"""
		# Filename created by igfft program	
		outfile_igfft = output_file_prefix + '.igfft.alignment'
		# Filename created by parsing program
		outfile_annotation = output_file_prefix + '.igfft.annotation'		
		# use deepcopy so we cant screw with original settings		
		fxn_params = copy.deepcopy(igfft_settings)
		fxn_parse = copy.deepcopy(parsing_settings)
		# Step 1: Run the igfft program and return the settings used for the program (command_val)
		command_val = run_igfft(f1, outfile=outfile_igfft, germline_source=germline_source, germlines=germline_location_settings, input_filetype=file_type.upper(), variable_parameters=fxn_params)		
		# Step 2: Run the parsing function which should also perform cdr3 analysis and isotype analysis		
		parse_alignment_file(outfile_igfft, outfile_annotation, fxn_parse, command_val=command_val, return_igrep_doc_line=return_igrep_doc_line)
		# add outputfiles to queue
		outfile_queue.put([outfile_annotation, outfile_igfft, index])		
		
	if num_processes == 1:
		# No need to muliprocess,so just run function 
		calligfft(input_file, of1_prefix)				
		final_files = outfile_queue.get()[0:2]
	else:			
		# First split files into mulitple small files 
		if file_type.upper() == 'FASTQ':
			split_files = useful.split_files_by_seq(input_file, num_processes, number_lines_per_seq=4, contains_header_row=False)
		else:
			split_files = useful.split_files_by_seq(input_file, num_processes, number_lines_per_seq=2, contains_header_row=False)
		proc = []
		for n, small_file in enumerate(split_files):
			# Annotate split files in parallel
			p = Process(target=calligfft, args=(small_file, small_file, n)) 	
			p.start()
			# Add a small delay between files 
			time.sleep(2)
			proc.append(p)
		for p in proc:
			p.join()

		# return all files generated by queue , they may be out of order depending upon upon when analysis completes; sort files so that they appear int he order they were submitted to run 
		results = sorted([outfile_queue.get() for i in range(len(split_files))], key=lambda x: x[2])	
		generated_annotated_files = [f[0] for f in results] 
		generated_alignment_files = [f[1] for f in results]

		# Merge the files back into a single file 
		print('annotation complete, merging files...')
		if return_igrep_doc_line:
			# This file has 2 header lines => database translator, and fields
			useful.merge_multiple_files(generated_annotated_files, 2, of1_prefix + '.igfft.annotation')
		else:
			# this file has 1 header lines => Fields
			useful.merge_multiple_files(generated_annotated_files, 1, of1_prefix + '.igfft.annotation') 
		useful.merge_multiple_files(generated_alignment_files, 1, of1_prefix + '.igfft.alignment')
		print('files merged.')
		# Delete all of the split files 		
		# delete split files and delete any files created from split files
		purge(split_files)		
		final_files = [of1_prefix + '.igfft.annotation', of1_prefix + '.igfft.alignment']

	return final_files


def run_igfft(dataset_path, germline_source="default", germlines=None, outfile=None, variable_parameters={}, input_filetype='FASTQ', header_field='header', sequence_field='sequence', quality_field='phred'):
	"""
		Python wrapper for running igfft binary from terminal

		Parameters
		----------
		dataset_path : string,
			Filepath of sequences to be annotated
		germline_source : string, default="default"
			Defines where the program will search for germlines
			If default: search for germlines in the default folder defined by igfft_germline_files
			If custom: the location of germline files will be explicity defined
			If db: Germlines will be downloaded from the IGREP database
		germlines : dict, default = None
			Defines settings for identifying germlines. Settings depends on the germline_source defined above
			..important::germlines=None
				if this field is not defined by the user, then we assume that the path is defined for v and j germlines in the dict variable_parameters
			If germline_source = default:
				germlines = {'SPECIES': list of species, 'LOCUS': list of loci}
			If germline_source = custom: 
				germlines = {'v': path to vgermline set, 'd': path to d germline set, 'j': path to j germline set}
			If germline_source = db:
				dict for querying database for germlines:				
				i.e. {_id: [list of germline ids]}
				i.e. {'SPECIES': 'Homo sapiens', 'LOCUS': 'IGH', 'PRODUCTIVITY': 'P'}
				..note::Allowed keys in query:
					The following keys are allowed in the query: ['_id', 'SPECIES', 'PRODUCTIVITY', 'LOCUS', 'MOLECULAR_COMPONENT', 'SOURCE', 'VERSION']
		outfile : string, default=None
			If defined, it will be the source of the output file
		variable_parameters : dict, default = {}
			Set of key, values for defining settings when running igfft. Keys for settings can be found by running ./igfft --defaults
		input_filetype : string, default="FASTQ"
			Format of the input file (dataset_path)
			Allowed values are None, FASTA, FASTQ, TAB, JSON
		header_field : string, default="header"
			If the format is not FASTQ/FASTA, we need to know what field in the file corresponds to the sequence header.
			This string defines the field name in the file that represents the DNA sequence header 
		sequence_field : string, default="sequence"
			If the format is not FASTQ/FASTA, we need to know what field in the file corresponds to the sequence.
			This string defines the field name in the file that represents the DNA sequence 
		phred_field : string, default="phred"
			If the format is not FASTQ, we need to know what field in the file corresponds to the sequence quality score.
			If there is no quality field, enter None

		Returns
		-------
		Settings used for annotating dataset

		Examples
		-------
		Annotate a FASTQ file, containing heavy and light chain sequences from humans, using the default germline database
		>>>run_igfft("folder/filetest.fastq", germlines={'SPECIES': 'Homo sapiens', 'LOCUS':['IGH', 'IGK', 'IGL']})
		Annotate a FASTA file, containing heavy and light chain sequences from humans, using the default germline database
		>>>run_igfft("folder/filetest.fasta", germlines={'SPECIES': 'Homo sapiens', 'LOCUS':['IGH', 'IGK', 'IGL']}, input_filetype='FASTA')
		Annotate a TAB file, containing heavy and light chain sequences from humans, using the default germline database. The DNA sequence is found under the column name 'dnaseq'.
		>>>run_igfft("folder/filetest.fasta", germlines={'SPECIES': 'Homo sapiens', 'LOCUS':['IGH', 'IGK', 'IGL']}, input_filetype=None (or 'TAB'), header_field='header', sequence_field='dnaseq')
		Annotate a FASTQ file, containing heavy and light chain sequences from humans, using the database
		>>>run_igfft("folder/filetest.fastq", germline_settings="db", germlines={'SPECIES': 'Homo sapiens', 'LOCUS':['IGH', 'IGK', 'IGL']})
		Annotate a FASTQ file but change the number of germlines to report. In this case, return the top two germlines.
		>>>run_igfft("folder/filetest.fastq", germlines={'SPECIES': 'Homo sapiens', 'LOCUS':['IGH', 'IGK', 'IGL']}, variable_parameters={"num_hits": 2})		
	"""
	
	# first go through the variable germlines to figure out the germline files to use 
	global databaseFolder
	default_germline_folder = databaseFolder
	
	if not germlines:
		# we assume that the user manually provided v and j parameters to variable parameters. therefore, this is a CUSTOM source 		
		custom_details = {}
		v = variable_parameters.pop('v', None)
		j = variable_parameters.pop('j', None)
		if v:			
			custom_details['v'] = v
		if j:			
			custom_details['j'] = j
		if custom_details == {}:
			# the user has not provided any information at all regarding a germline 
			raise Exception("You must provide at least one germline database. If you do not have a germline database then define the variable germlines to use default germline database in program or download germline from GG lab database. See documentation in function")
		else:
			germline_source = "custom"
			germlines = custom_details
		
	if not isinstance(germlines, dict):
		raise Exception('Germline definition must be a dict')
	
	if type(variable_parameters) is not dict:
		raise Exception('Incorrect format of variable_parameters. The correct format of the function variable_parameters is a dictionary. keys correspond to the parameters used in the program and values correspond to the settings for that parameter')

	# IGNORE d FOR NOW  (igfft requires lowercalse subtypes)
	subtypes = ['v', 'j']
	germline_files = {}
	if germline_source == 'custom':		
		germline_settings = {'GERMLINE-SOURCE': 'CUSTOM-FILES', 'PARAMS': germlines}							
		# the user is providing proper germline files for each 
		for key_subtypes, pathlist in germlines.iteritems():
			keys = key_subtypes.lower()
			if keys not in subtypes:
				print('Warning: the provided subtype will not be included: ' + keys)
				continue
			if not isinstance(pathlist, list):
				pathlist = [pathlist]
			updated_pathlist = []
			for paths in pathlist:
				if os.path.isfile(paths):
					updated_pathlist.append(paths)
				else:
					print("Warning: the provided germline path does not exist and will not be included in germlines: " + paths)										
			germline_files[keys] = updated_pathlist
			
	elif germline_source == 'default':		
		germline_files = {}
		germline_settings = {'GERMLINE-SOURCE': 'DEFAULT FOLDER', 'PARAMS': germlines}	
		# the user wants to use the default germlines from the program 		
		species = [germlines['SPECIES']] if not isinstance(germlines['SPECIES'], list) else germlines['SPECIES']
		locus = [germlines['LOCUS']] if not isinstance(germlines['LOCUS'], list) else germlines['LOCUS']
		for st in subtypes:
			germline_files[st] = []
			for s in species:
				s = s.replace(' ', '').lower()
				for l in locus:				
					l = l.replace(' ', '').lower()
					germline_file = os.path.join(default_germline_folder, s, st, '_'.join([s, l, st]) + '.txt')				
					if not os.path.isfile(germline_file):
						print("Warning: the provided germline path does not exist and will not be included in germlines: " + germline_file)
						continue
						# raise Exception()
					germline_files[st].append(germline_file)						
	elif germline_source == 'db':
		# the user wants to download germlines from database		
		germline_settings = {'GERMLINE-SOURCE': 'DATABASE-DOWNLOAD'}		
		[germline_files, germline_settings['PARAMS']] = download_germlines_from_db(germlines)				
	else:		
		raise Exception('Germline source is only allowed to be equal to the following strings: custom, default, db')			
	variable_parameters.update(germline_files)			

	# make sure igfft parameters are good in variable_parameters dictionary  
	default_parameters = GetDefaultParameters()

	default_parameters = {v: default_parameters[v] for v in default_parameters}
	if not outfile:
		# User did not define path, so define it manuall
		outfile = useful.removeFileExtension(dataset_path) + '.igfft.alignment'
			
	variable_parameters['o'] = outfile
	# this will remove subtype preffix for parameters (for example, if parameter is num_gaps, but we are only modifying vgene, then parameter will be vnum_gaps. this will remove the v part of the string)	
	ignore_seq_prefix = {var: var[1:] if var[0] in ['v', 'd', 'j'] and len(var) > 1 else var for var in variable_parameters}
	# remove any paramters that are already set as default in the program
	remove_vals_set_as_default = [var for var in variable_parameters if ignore_seq_prefix[var] in default_parameters and default_parameters[ignore_seq_prefix[var]] and variable_parameters[var] == default_parameters[ignore_seq_prefix[var]]] 
	
	for var in remove_vals_set_as_default:
		variable_parameters.pop(var)
				
	deletetemp = False
	# FASTA/FASTQ FILE FORMAT IS ALREADY SUPPORTED BY IGFFT
	if input_filetype is 'FASTA':
		tempFile = dataset_path
		filetype = 'FASTA'
	elif input_filetype is'FASTQ':
		tempFile = dataset_path
		filetype = 'FASTQ'
	else:
		# Read in the input file and convert it to a TAB file so that that can be read by igfft	
		deletetemp = True
		inputIFFile = readwrite.immunogrepFile(filelocation=dataset_path, filetype=input_filetype)
		date = str(datetime.now())
		remove_chars = [':', '-', ' ', '.']
		for char in remove_chars:
			date = date.replace(char, '')		
		# Temporary TAB file we will make for running program 
		tempFile = dataset_path + '_{0}'.format(date)
		output_temp_file = open(tempFile, 'w')		
		output_temp_file.write('HEADER\tSEQUENCE\tSEQUENCE_QUALITY\n')
		for line in inputIFFile.read():		
			if line:			
				db_info = '{0}{1}'.format(fasta_file_delimiter, json.dumps({idIdentifier: line[idIdentifier]})) if idIdentifier in line else ''
				quality = line[quality_field] if quality_field in line else ''
				if sequence_field in line:													
					output_temp_file.write('{0}{2}\t{1}\t{3}\n'.format(line[header_field].strip(), line[sequence_field].strip(), db_info, quality))
		output_temp_file.close()
		inputIFFile.IFclass.close()
		filetype = 'TAB'
	
	variable_parameters['i'] = filetype
	
	print('\n\nRunning IgFFT using the following settings:\n')
	pprint(variable_parameters)		
	running_program_text = "Input file location: {0}\n".format(tempFile)
	running_program_text += "Output file will be saved as: {0}\n".format(outfile)
	running_program_text += "IgFFT Run has started at {0}\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()))			
	# NOW RUN THE BINARY 
	print(running_program_text)
	command_string = '''"{1}" "{0}" '''.format(tempFile, igfft_location)	
	for var in variable_parameters:
		if var == 'o':
			# OUTPUT FILE ADD DOUBLE QUOTES
			command_string += '-{0} "{1}" '.format(var, variable_parameters[var])			
		elif var in ['v', 'd', 'j']:
			newsets = [os.path.abspath(g).replace(' ', '\ ') for g in variable_parameters[var]]
			command_string += '-{0} {1} '.format(var, ' '.join(newsets))
			variable_parameters[var] = [os.path.basename(n) for n in newsets]
		else:
			command_string += '-{0} {1} '.format(var, variable_parameters[var])					
	subprocess.call(command_string, shell=True)		

	# program complete 
	print("Analysis Completed at {0}\nAnalysis saved to: {1}\n\n\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()), outfile))	
	if os.path.isfile(tempFile + '.convert_tab.txt'):
		# c++ program produces this, just delete it
		os.remove(tempFile + '.convert_tab.txt')
	if deletetemp:
		# delete any temp files that were made 
		os.remove(tempFile)

	command = variable_parameters
	# Dont store output path 
	command.pop('o', None) 
	# Dont store input path 
	command.pop('i', None)
	command['germline-details'] = germline_settings	
	command['annotation'] = 'IGFFT'
	return command


def parse_alignment_file(annotatedFile, outfile=None, annotation_settings={}, command_val={}, return_results_dict=False, return_igrep_doc_line=False):
	"""
		Function  for parsing the results of the IGFFT annotation file. 
		It will go through the results and ensure results were found, identify CDR3 if requested, remove insertions if requested, and identify isotype sequences if requested

		Parameters
		----------
		annotatedFile : string
		Filepath of the annotated file output from igfft program
		outfile : string, default None
		Filepath of the file after being parsed. If path is none, then output file will be = annotatedFile+'.annotation'
		annotation_settings : dict, default {}
		These keys will define how the user wants to parse the annotation results. If empty, then default settings will be used.
		The user can modify the following fields 'min_cutoff'=min percent id to germline, 'min_len_cutoff'=min alignment length to germline, 'remove_insertions': 0 => never, 2=> only when stop codon, 1=> always
		'isotype'={'barcode-list':[], 'penalize_truncations':, 'maxmismatch':, 'min-len':, 'iso_search_direction' } => use isotpying if 'isotype' is a key
		'cdr3_search'={'SPECIES':[], 'LOCUS':[]}
		command_val : dict, default {}
		Dict defining the command used in the igfft program. This is just an optional way to keep track of annotation settings.		
		return_igrep_doc_line : boolean, default False
		If True, will return a comment line above TAB file defining how to insert results into the database
	"""

	# ##JUST SETTING UP SETTINGS AND PARAMETERS### #
	write_format = 'TAB'	
	numSeqs = useful.file_line_count(annotatedFile) - 1
	dna_alphabet = "ACTGUKMRYSWBVHDXN"
	bad_dna_char_pattern = re.compile('[^' + dna_alphabet + ']', re.IGNORECASE)
	
	print('Parsing and summarizing alignment file')	
	
	if not(outfile):
		outfile = useful.removeFileExtension(annotatedFile) + '.igfft.annotation'
	outfile_unknown_rec = outfile + ".unk_recombination"

	# This file will be used to note any errors that occur while parsing the file
	outfile_errors = outfile + ".igfft.error_log"
	
	if annotatedFile == outfile:
		os.rename(annotatedFile, annotatedFile + '.temp')		
		annotatedFile += '.temp'
	
	# Decides what min cutoff (percent identity to germline) to use to consider a successful sequence match
	if 'min_cutoff' in annotation_settings:
		min_per_algn_cutoff = annotation_settings['min_cutoff']
	else:
		min_per_algn_cutoff = 0.2
	
	# minimum number of sequences that must align to germlines to be considered antibody
	if 'min_len_cutoff' in annotation_settings:
		min_algn_len_cutoff = annotation_settings['min_len_cutoff']
	else:
		min_algn_len_cutoff = 50
	
	# remove base insertions detected in sequence alignment
	if 'remove_insertions' in annotation_settings:
		remove_insertions = annotation_settings['remove_insertions']
	else:
		# 0 -> dont remove insertions, 
		# 1 -> remove insertions always, 
		# 2-> remove insertions only if there is a stop codon
		remove_insertions = 0 						
	
	# Determine whether we will perform isotyping
	isotype_aligner = None

	if 'isotype' in annotation_settings and isotypeworking:
		iso_settings = annotation_settings['isotype']
		identifyIsotype = True
		isotype_barcodes = iso_settings['barcode-list'] if 'barcode-list' in iso_settings and iso_settings['barcode-list'] else copy.deepcopy(defaultBarcodes())	
		p_t = iso_settings['penalize_truncations'] if 'penalize_truncations' in iso_settings else True
		num_mismatch = iso_settings['maxmismatch'] if 'maxmismatch' in iso_settings else 2
		minimum_iso_alignment_length = iso_settings['min-len'] if 'min-len' in iso_settings else 15
		# only consider forward direction of barcodes
		search_rc_isotype = iso_settings['iso_search_direction'] if 'iso_search_direction' in iso_settings else 0 
		isotype_aligner = fftaligner.BarcodeAligner(isotype_barcodes, p_t, search_rc_isotype, num_mismatch, minimum_iso_alignment_length)
	else:
		if isotypeworking is False:
			print("Warning: PYFFTW is not installed. Therefore isotyping will not be possible")
			
		identifyIsotype = False
	
	if 'cdr3_search' in annotation_settings:
		# Set up settings for finding cdr3 from results
		identifyCDR3 = True
		cdr3_search_parameters = CDR3tools.findCDR3(annotation_settings['cdr3_search'])
		motifs = cdr3_search_parameters.get_motifs()
		unique_locus_sets = set([a[0].upper() for a in motifs])
		# Find the maximum length of the poossible motifs
		max_l_motif_len = 4 * max(sorted([int(k) for a in motifs for k in a[1]]))
		# Find the maximum lenght of the possible right motifs
		max_r_motif_len = 4 * max(sorted([int(k) for a in motifs for k in a[2]]))
		cdr3_search_name = annotation_settings['cdr3_search']
	else:
		identifyCDR3 = False
		cdr3_search_name = "None selected"
		print('No parameters were provided for identifying the CDR3. CDR3 analysis will not be included')
				
	# CHECK IF IGFFT FILE WAS ACTUALLY CREATED	
	if not(os.path.isfile(annotatedFile)):			   
		raise Exception('ERROR: IGFFT FILE WAS NOT CREATED.PLEASE MODIFY SETTINGS AND/OR RE-RUN')		
	
	foutfile = open(outfile, 'w')
	foutfile_unknown_recombination = open(outfile_unknown_rec, 'w')
	ferrorlog = open(outfile_errors, 'w')
	output_files = {'annotated': foutfile, 'UNK': foutfile_unknown_recombination, 'ERROR': ferrorlog}

	if return_igrep_doc_line:
		# Save the DB translator
		# in the first line of the file output a comment line definining how to convert this file into a database file 
		translator = DatabaseTranslator()
		translator[filetype_var] = write_format 
		translator_comment = descriptor_symbol
		translator_string = json.dumps(translator)
		foutfile.write(translator_comment + translator_string + '\n')

	# Save the settings used for parsing, again this is just for convenience to store the settings you used in file	
	command_val['Cdr3-Fr4Identification'] = json.dumps({'method': 'gglab', 'param': cdr3_search_name}) if identifyCDR3 else 'not selected'
	command_val['Min_Alignment_Idendity_Cutoff'] = str(min_per_algn_cutoff)
	command_val['Min_Alignment_Length_Cutoff'] = str(min_algn_len_cutoff)
	if identifyIsotype:
		command_val['Isotyping'] = {'Barcodes': isotype_barcodes, 'mismatch_cutoff': num_mismatch, 'penalize_truncations': p_t, 'minimum_length_cutoff': minimum_iso_alignment_length}					
	if remove_insertions == 0:
		command_val['Fix Mutations'] = 'Never' 
	elif remove_insertions == 1:
		command_val['Fix Mutations'] = 'Always'
	elif remove_insertions == 2:
		command_val['Fix Mutations'] = 'WhenStopCodon'
	commandString = json.dumps(command_val)
	
	isostring = '\t\tPerforming isotyping on sequences using provided barcodes: ' + '\n' if 'isotype' in annotation_settings else ''	
	parse_igfft_notification = '''
	The resulting annotation output will be parsed using the following settings:
	\tMinimum alignment length cutoff: {0},
	\tMinimum percent identity: {1},
	\tFix detected insertions: {2},
	\tIdentify CDR3: {3}\n{4}
	\n\n
	'''.format(str(min_algn_len_cutoff), str(min_per_algn_cutoff), command_val['Fix Mutations'], str(identifyCDR3), isostring)
	print(parse_igfft_notification)
		
	useDebug = False			
	# ## ALL PARAMETERS HAVE BEEN SET####

	vdj_hits = 0
	vj_hits = 0
	unknown_hits = 0	
	f1 = readwrite.immunogrepFile(filelocation=annotatedFile, filetype='TAB')	
	count = 0	
	startPer = 0	
	chainDic = {
		'IGH': ["heavy", "VDJ", "IGH"],
		'IGK': ["light", "VJ", "IGK"],
		'IGL': ["light", "VJ", "IGL"],
		'TRA': ["alpha", "VJ", "TRA"],
		'TRB': ["beta", "VDJ", "TRB"],
		'VH': ["heavy", "VDJ", "IGH"],
		'VK': ["light", "VJ", "IGK"],
		'VL': ["light", "VJ", "IGL"],
		'TA': ["alpha", "VJ", "TRA"],
		'TB': ["beta", "VDJ", "TRB"]
	}
	
	chainTypes = {
		'heavy': 'VDJ',
		'light': 'VJ',
		'alpha': 'VJ',
		'beta': 'VDJ'
	}
	
	debugMe = False
	useDebug = False
	igblast_read_line = True
	total_parsing_errors = 0
	loop_status_gen = useful.LoopStatusGen(numSeqs, 10) 
	tab_header_var = TABFileHeader()
	overlap_len = 10 
	# create header row
	if write_format is 'TAB':		
		output_files['annotated'].write('\t'.join([field for field in tab_header_var]) + '\n')
	
	# while count<numSeqs: #read through every sequence in the file
	for query_results in f1.read():		
		count += 1
		# allow the user to monitor what percent of the sequences have been processed
		loop_status_gen.next()
		if not(query_results):
			continue	
		[seqHeader, additionalInfo] = readwrite.GetAdditionalInfo(query_results['Header'])				
		query_results['Document_Header'] = query_results['Header']
		query_results['Header'] = seqHeader										
		additionalInfo.pop('Header', None)
		additionalInfo.pop('document_header', None)		
		query_results = defaultdict(str, dict(query_results.items() + additionalInfo.items()))		
		query_results['Quality_Score'] = query_results['Sequence quality']
		error_dic = defaultdict(str, {'Sequence': query_results['Sequence'], 'Header': query_results['Header'], 'Document_Header': query_results['Document_Header']})
		if idIdentifier in query_results:
			error_dic[idIdentifier] = query_results[idIdentifier]			
		if query_results['Sequence'] != "":
			seq = query_results['Sequence']			
		else:			 
			error_dic["Notes"] = "Sequence not found"
			error_dic["Errors"] = "Sequence not found"
			error_dic["Percent_Identity"] = None
			error_dic["Alignment_Length"] = None					
			Write_Seq_TAB(error_dic, tab_header_var, output_files['annotated'])
			continue			
		try:							
			# CHECK TO SEE IF SEQUENCE HAS WEIRD CHARACTERS IN IT
			if bad_dna_char_pattern.search(seq):
				notes = "Sequence contains unknown characters"
				print("Seq # " + str(count) + " contains unusual characters and therefore we are ignoring this sequence: " + seqHeader)				
				error_dic["Percent_Identity"] = None
				error_dic["Alignment_Length"] = None
				error_dic['Notes'] = notes
				error_dic['Errors'] = notes								
				Write_Seq_TAB(error_dic, tab_header_var, output_files['annotated'])				
				continue										
			notes = query_results['Notes']			
						
			algn_seq = query_results['Strand_Corrected_Sequence']
			if algn_seq == "":		
				error_dic = query_results
				error_dic['Notes'] = notes
				error_dic['Errors'] = notes							
				error_dic['Percent_Identity'] = None
				error_dic['Alignment_Length'] = None
				
				#Write_Seq_JSON(error_dic,tab_header_var,output_files['ERROR'])
				if write_format == "TAB":								
					Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
				else:																									   				
					Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])			
				continue		
			
			start_of_ab = int(query_results['Query_Start'])
			end_of_ab = int(query_results['Query_End'])
			query_translating_codon = int(query_results['Codon_Start'])
			query_translating_frame = int(query_results['Codon_Frame'])		
			
			if not(query_results['VGENE: Query_Start'] == '') and not(query_results['VGENE: Query_End'] == ''):
				v_start = int(query_results['VGENE: Query_Start'])
				v_end = int(query_results['VGENE: Query_End'])
				algnLenV = v_end - v_start + 1
				matchValsV = int(query_results['VGENE: Total_Matches'])
				vpresent = True						
			else:
				v_start = -1
				v_end = len(algn_seq) + 1
				algnLenV = 0
				matchValsV = 0
				vpresent = False
								
			if not(query_results['JGENE: Query_Start'] == '') and not(query_results['JGENE: Query_End'] == ''):
				j_start = int(query_results['JGENE: Query_Start'])
				j_end = int(query_results['JGENE: Query_End'])
				algnLenJ = j_end - j_start + 1
				matchValsJ = int(query_results['JGENE: Total_Matches'])
				jpresent = True  
				queries = [query_results['JGENE: Alignment_Sequence_Query'].split(';')[0]]; #fornow, just take the top germline
				germlines = [query_results['JGENE: Alignment_Sequence_Germline'].split(';')[0]];#fornow, just take the top germline
				jGermlineInfo = []
				for q,g in enumerate(germlines):
					jGermlineInfo.append(			
						{
							'query_seq':queries[q],
							'germline_seq':g,
							'start':j_start				
						}
					)
			else:
				j_start = -1
				fr4start = -1
				j_end = len(algn_seq) + 1
				algnLenJ = 0
				matchValsJ = 0
				jpresent = False
				jGermlineInfo = None
			
			algnLen = algnLenV + algnLenJ
			algnLen += 1 if algnLen <= 0 else 0
			matchVals = matchValsV + matchValsJ
			
			perIden = matchVals / float(algnLen)
			
			query_results['Percent_Identity'] = perIden
			query_results['Alignment_Length'] = algnLen
			
			if perIden < min_per_algn_cutoff or algnLen < min_algn_len_cutoff:							
				keep_fields = ['Sequence',idIdentifier,'Header','Document_Header','Strand_Corrected_Sequence','Locus','Quality_Score','5_Prime_Annotation','3_Prime_Annotation','Direction','Full_Length_Sequence.NT','Full_Length_Sequence.AA','Codon_Start','Query_Start','Query_End','Codon_Frame']
				error_dic = defaultdict(str,{a:query_results[a] for a in keep_fields})
				notes = 'Sequence results did not pass alignment filter settings after parsing output;' + notes								
				error_dic['Notes'] = notes
				error_dic['Errors'] = notes							
				error_dic['Percent_Identity'] = perIden
				error_dic['Alignment_Length'] = algnLen
				#Write_Seq_JSON(error_dic,chain_ind_fields,chain_dep_fields,'',output_files['ERROR'])
				#Write_Seq_JSON(error_dic,tab_header_var,output_files['ERROR'])
				if write_format == "TAB":								
					Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
				else:																									   				
					Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])		
				continue		
			
			
			newFrame = query_results['VGENE_Reading_Frames: FR1,CDR1,FR2,CDR2,FR3,CDR3'].split(',')
			if len(newFrame)<7:
				for ifr in range(len(newFrame),7):
					newFrame.append('')
			
			#newFrame.append('')
			if not(query_results['Locus'] == ''):
				locus = query_results['Locus'].split(',')
			else:
				locus = [k for k in chainDic.keys() for genes in ['Top_V-Gene_Hits','Top_J-Gene_Hits'] for hits in query_results[genes].split(',')[0].split(' ') if hits.upper().startswith(k) ]		
					
			locus = set(locus)
			locus_list = list(locus)
					
			if not(query_results['Chain'] == ''):
				chain = query_results['Chain'].split(',')			
				var_type = chainTypes[chain[0].lower()] if chain[0].lower() in chainTypes else 'UNK'									
			else:
				chain = chainDic[locus_list[0]][0] if len(locus_list) > 0 and locus_list[0] in chainDic else 'UNK'
				var_type = chainDic[locus_list[0]][1] if len(locus_list) > 0 and locus_list[0] in chainDic else 'UNK'
			
			query_results['Chain'] = chain
				
			
			query_results['Recombination_Type'] = var_type
			query_results['Locus'] = ','.join(locus_list)
			
			if  var_type == 'UNK':			
				print('This chain type was not recognized in the parsing script. Analysis information for this sequence will be placed in the file "{0}". Consider updating the variable chainTypes in the funcion "parse_alignment_file"'.format(outfile_unknown_rec))
				error_dic = query_results
				error_dic['Notes'] = "Chain type not recognized"
				error_dic['Errors'] = "Chain type not recognized. analysis info placed in the file '{0}'".format(outfile_unknown_rec)				
				unknown_hits+=1
				#Write_Seq_JSON(error_dic,chain_ind_fields,chain_dep_fields,'',output_files['ERROR'])
				Write_Seq_JSON(error_dic,tab_header_var,output_files['unk_recombination'])
				if write_format == "TAB":								
					Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
				else:																									   				
					Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])		
						 
			
			if not(query_results['VGENE: Query_FR3_Start::End'] == '') and not(query_results['VGENE: Query_FR3_Start::End'] == '::'):
				fr3present = True				
				fr3_pos = [int(p) for p in query_results['VGENE: Query_FR3_Start::End'].split('::')]				
				
				cdr3start = fr3_pos[1] + 1			
			else:
				fr3present = False
				cdr3start = -1
				fr3_pos = []
			
			
			if identifyCDR3:		
				locus = locus_list if locus<=unique_locus_sets else []				
				
				[cdr3_nt, cdr3_aa, fr4_nt, fr4_aa, cdr3Frame, fr4Frame, cdr3start, fr4start, cdr3notes] = CDR3_search(cdr3_search_parameters, algn_seq, cdr3start, end_of_ab, locus, max_l_motif_len, max_r_motif_len)																	

				query_results['CDR3_Sequence.NT'] = cdr3_nt
				query_results['CDR3_Sequence.AA'] = cdr3_aa
				query_results['FR4_Sequence.NT'] = fr4_nt
				query_results['FR4_Sequence.AA'] = fr4_aa						
				notes+=cdr3notes
				
				if cdr3_nt!='':														
					query_results['Query_CDR3_Start::End'] = str(cdr3start)+'::'+str(fr4start-1)				
					newFrame[5] = str(cdr3Frame)	
					
															
					if fr3present and cdr3start>0 and (cdr3start-1)!=fr3_pos[1]: #this means that after the cdr3 analysis search. our end position for the fr3 has changed. we need to update values
					
						#adjusting.....fr3...................					
						fr3_pos_new = cdr3start-1#ok fr3 actually ends at this position now...
						query_results['VGENE: Query_FR3_Start::End'] = str(fr3_pos[0])+"::"+str(fr3_pos_new) #update the start and end positions for the FR3 with respect to original sequence
						query_results['FR3_Sequence.NT'] = algn_seq[fr3_pos[0]:fr3_pos_new+1]#now update the nucleotide sequence
					
					
						#now we have the update the positions corresponding to the "aligned sequences" (algned query and germline)
						new_algn_pos = query_results['VGENE: Alignment_FR3_Start::End'].split('::')
						fr3_algn_pos = [int(p) for p in new_algn_pos]												
					
						query_algn_seq = query_results['VGENE: Alignment_Sequence_Query'][fr3_algn_pos[0]:fr3_algn_pos[1]] #this is the alignmetn string for the query
						
						index_count = fr3_algn_pos[1]-fr3_algn_pos[0]
						s_count = fr3_algn_pos[1] #this is our current character position with respect to alignment sequence
						if fr3_pos_new<fr3_pos[1]: #we have to move backwwords
							seqpos = fr3_pos[1] #the original position was fr3_pos[1] alon the query						
							while (index_count>0 and index_count<len(query_algn_seq) and seqpos>fr3_pos_new):
								if(query_algn_seq[index_count]!='-'): #if there is no deletion at that position, then decrease the sequence position by 1
									seqpos-=1
								index_count-=1
								s_count-=1								
						else: #then we traverse forwards
							seqpos = fr3_pos[1]																				
							while (index_count<len(query_algn_seq)-1 and index_count>=0 and (seqpos<fr3_pos_new or query_algn_seq[index_count]=='-')): #the second part of the if statmenet is provided to protect against a deletion detected right before the new position. so if its ac-g and the new start position is g, then this will ensure that the count starts at g and not -
								if(query_algn_seq[index_count]!='-'):
									seqpos+=1
								index_count+=1								
								s_count+=1						
						query_results['VGENE: Alignment_FR3_Start::End'] = str(fr3_algn_pos[0])+"::"+str(s_count)
						#fr3......adjusted........					
					
															
				else:
					query_results['Query_CDR3_Start::End'] = '::'
					
				if fr4start>-1:
					query_results['JGENE: Query_FR4_Start::End'] = str(fr4start)+'::'+str(end_of_ab)
					newFrame[6] = str(fr4Frame)
				else:
					query_results['JGENE: Query_FR4_Start::End'] = '::'			
			else:
				cdr3_nt = ''
				fr3_nt = ''
				cdr3Frame=0
				fr4Frame=0
			
	
			algn_info = []
			adjustedFrame =  query_translating_codon - start_of_ab
			
			
			if vpresent:#perfrom framework annotation translation and remove mutations 
				annotations = ['FR1','CDR1','FR2','CDR2','FR3']
				annotations_present = [(i,region) for i,region in enumerate(annotations) if (query_results['VGENE: Query_{0}_Start::End'.format(region)]!='' and  query_results['VGENE: Query_{0}_Start::End'.format(region)]!='::')]
				last_base = int(query_results['VGENE: Query_{0}_Start::End'.format(annotations_present[-1][1])].split('::')[1])+1
				
				query_results['Full_Length_Sequence.AA'] = TranslateSeq(query_results['Full_Length_Sequence.NT'],adjustedFrame)# Seq(query_results['Full_Length_Sequence.NT'][adjustedFrame:],generic_dna).translate().tostring()				   
				contains_stop_codons = '*' in query_results['Full_Length_Sequence.AA']					 
				query_results['VREGION.NUM_GAPS'] = int(query_results['VGENE: Total_Indel'])
				numDel = 0#if more than 0, then, we will be replacing field Full_Length_Sequence.NT with new_nt_seq
				if query_results['VREGION.NUM_GAPS'] > 0:
					indel_profile = query_results['VGENE_Indels: FR1,CDR1,FR2,CDR2,FR3,CDR3'].split(',')					 
					startingNt = query_translating_codon
					new_seq_nt = ''# query_results['Full_Length_Sequence.NT'][0:adjustedFrame]
					
					mutation_notes = ''
					
					for i,region in annotations_present:										
						newFrame[i] = str(startingNt%3+1)  
						update_annotated_seq = False				  
						
						translation_frame = adjustedFrame if region == query_results['5_Prime_Annotation'] else  0 #this is just the frame we will be translating the substring. not the frame of the region with respect to full length sequence					
						
						if int(indel_profile[i])>0:
							algn_query_substring_pos = query_results['VGENE: Alignment_{0}_Start::End'.format(region)].split('::')
							algn_query_substring = query_results['VGENE: Alignment_Sequence_Query'][int(algn_query_substring_pos[0]):int(algn_query_substring_pos[1])+1]
							algn_germline_substring = query_results['VGENE: Alignment_Sequence_Germline'][int(algn_query_substring_pos[0]):int(algn_query_substring_pos[1])+1]
																							 
							if remove_insertions == 1:
								seq_with_gaps = ''
								for j,char in enumerate(algn_germline_substring):
									if j<translation_frame:
										continue #only modify bases that are after translation frame
									if char != '-':									
										seq_with_gaps+=algn_query_substring[j]
									else:
										numDel+=1 #if more than 0, then, we will be replacing field Full_Length_Sequence.NT with new_seq_nt
										update_annotated_seq = True #f true, we will be replacing subsequence of annotated region
										
							elif remove_insertions == 2:
								if contains_stop_codons:
									seq_with_gaps = ''
									for j,char in enumerate(algn_germline_substring):
										if j<translation_frame:
											continue #only modify bases that are after translation frame
										if char != '-':										
											seq_with_gaps+=algn_query_substring[j]
										else:
											numDel+=1
											update_annotated_seq = True
											
								else:
									seq_with_gaps = algn_query_substring							
							else: #remove insertions = 0, or some weird value
								seq_with_gaps = algn_query_substring	
																		   
							if update_annotated_seq:
								query_results['{0}_Sequence.NT'.format(region)] = seq_with_gaps.replace('-','')
								mutation_notes += region+','
							query_results['{0}_Sequence.NT.Gapped'.format(region)] = algn_query_substring.replace('-','n') #this will give the dna sequence of the query containing both insertions and deletions (ignores any requested modifications in remove_insertions parameters)
							query_results['{0}_Sequence.AA.Gapped'.format(region)] = TranslateSeq(query_results['{0}_Sequence.NT.Gapped'.format(region)],translation_frame)#  Seq(query_results['{0}_Sequence.NT.Gapped'.format(region)][translation_frame:],generic_dna).translate().tostring() #amino acid of gapped dna sequence 
						else:						
							query_results['{0}_Sequence.NT.Gapped'.format(region)] = ''
							query_results['{0}_Sequence.AA.Gapped'.format(region)] = ''
						
						query_results['{0}_Sequence.AA'.format(region)] =  TranslateSeq(query_results['{0}_Sequence.NT'.format(region)],translation_frame)# Seq(query_results['{0}_Sequence.NT'.format(region)][translation_frame:],generic_dna).translate().tostring() #amino acid of gapped dna sequence 
	
						new_seq_nt += query_results['{0}_Sequence.NT'.format(region)]
						startingNt+=len(query_results['{0}_Sequence.NT'.format(region)])-translation_frame
					  	
					if numDel>0: #at least one mutation was made to the sequence
						if cdr3start>=0:#a cdr3 was found, so lets update the frame information for the CDR3 and FR4 which occur after V GENE
							cdr3Frame = (cdr3start-numDel)%3+1					   
							newFrame[5] = str(cdr3Frame)						
						if fr4start>=0:
							fr4Frame  = (fr4start-numDel)%3+1
							newFrame[6] = str(fr4Frame)
						query_results['Full_Length_Sequence.NT'] = new_seq_nt+query_results['Strand_Corrected_Sequence'][last_base:end_of_ab+1] #new_seq_nt is only the new v-gene sequence. we have to add remaining bases when translatign
						query_results['Full_Length_Sequence.AA'] = TranslateSeq(query_results['Full_Length_Sequence.NT'],adjustedFrame)# Seq(query_results['Full_Length_Sequence.NT'][adjustedFrame:],generic_dna).translate().tostring()
						mutation_notes = mutation_notes[:-1]+';'
						notes+= 'Insertions have been removed from: {0}'.format(mutation_notes)
				else: #there were no indels detected in alignment
					for i,region in annotations_present: #simply translate the provided nucleotide sequences. no gap correction is required 
						adjustedFrame = query_translating_codon - int(query_results['VGENE: Query_{0}_Start::End'.format(region)].split('::')[0]) if region == query_results['5_Prime_Annotation'] else 0					
						query_results['{0}_Sequence.AA'.format(region)] = TranslateSeq(query_results['{0}_Sequence.NT'.format(region)],adjustedFrame)# Seq(query_results['{0}_Sequence.NT'.format(region)][adjustedFrame:],generic_dna).translate().tostring()										
						query_results['{0}_Sequence.NT.Gapped'.format(region)] = ''
						query_results['{0}_Sequence.AA.Gapped'.format(region)] = ''							   
			   
			   	
				mutations = int(query_results['VGENE: Total_Mismatches']) +  query_results['VREGION.NUM_GAPS'] - numDel
				query_results['VRegion.SHM.NT'] = str(mutations) #adjust the number of gaps to reflect the removed insertions (numDel)
				query_results['VRegion.SHM.Per_nt'] = round(100*mutations/float(algnLenV),2)
																	 
			
			else:									
				if not(cdr3_nt == ''):
					start_of_ab = cdr3start				 
					query_results['5_Prime_Annotation'] = "CDR3"
					query_results['Full_Length_Sequence.NT'] = query_results['Strand_Corrected_Sequence'][cdr3start:end_of_ab+1]
					adjustedFrame = 0
				query_results['Full_Length_Sequence.AA'] = TranslateSeq(query_results['Full_Length_Sequence.NT'],adjustedFrame)# Seq(query_results['Strand_Corrected_Sequence'][cdr3start:end_of_ab],generic_dna).translate().tostring()	   
	
			
			if jpresent:
				jmutations = int(query_results['JGENE: Total_Mismatches']) +  int(query_results['JGENE: Total_Mismatches'])
				query_results['JRegion.SHM.NT'] = str(jmutations) #adjust the number of gaps to reflect the removed insertions (numDel)
				query_results['JRegion.SHM.Per_nt'] = round(100*jmutations/float(algnLenJ),2)
			
			#attempt to guess the isotype 
			if (vpresent or jpresent) and identifyIsotype:
				#take all dna sequences after substring 				
				start_p_iso = end_of_ab-overlap_len
				if start_p_iso<0:
					start_p_iso = 0 
				substring = algn_seq[start_p_iso:] 
				isotypes_results = isotype_aligner.AlignToSeq(substring)				
												
				if isotypes_results:
					query_results["Isotype"] = ','.join(isotypes_results['Isotype'])
					query_results["Isotype mismatches"] =  ','.join([str(s) for s in isotypes_results['Mismatches']])
					query_results['Isotype percent similarity'] = ','.join([str(s) for s in isotypes_results['Percent similarity']])								
					query_results['Isotype barcode direction'] = ','.join([str(s) for s in isotypes_results['Direction']])								
								
						
			query_results['Command'] = commandString				
			query_results['Stop_Codon'] = '*' in query_results['Full_Length_Sequence.AA']
			query_results['Full_Length'] = query_results['5_Prime_Annotation']=='FR1' and query_results['3_Prime_Annotation'] == 'FR4'
			
			query_results['CDR3_Junction_In_Frame'] = query_results['CDR3_Sequence.AA'] != '' and query_results['CDR3_Sequence.AA'] in query_results['Full_Length_Sequence.AA']
			query_results['Reading_Frames: FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4'] = ','.join(newFrame)
			if query_results['Stop_Codon']:
				query_results['Productive'] = 'NO'
			else:
			   if query_results['CDR3_Sequence.AA']=='':
				  notes+='There are no stop codons in the sequence, but the CDR3 sequence could not be found;'
				  query_results['Productive'] = 'MAYBE'
			   else:
				  if query_results['CDR3_Junction_In_Frame']:
				  	if query_results['5_Prime_Annotation'] == 'FR1':					
				  		query_results['Productive'] = 'YES'
					else:
						query_results['Productive'] = 'MAYBE'
					  
				  else:
					  notes+='The CDR3 is out-of-frame with respect to the translated sequence;'
					  query_results['Productive'] = 'NO'
			
			query_results["Notes"] = notes + query_results['Notes']
							
			if write_format == "TAB":								
				Write_Seq_TAB(query_results,tab_header_var,output_files['annotated'])
			else:																									   				
				Write_Seq_JSON(query_results,tab_header_var,output_files['annotated'])
				
		except Exception as e:
			total_parsing_errors+=1					
			exc_type, exc_obj, exc_tb = sys.exc_info()		
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		
			print("An error occurred when analyzing the output for sequence # "+str(count)+". This error has been reported in the error log file. Continuing with analysis")
			error_to_file = "Error Number: "+str(total_parsing_errors)
			error_to_file += "\nSequence count: "+str(count)
			error_to_file += "\nSequence header: "+seqHeader
			error_to_file += "\nSequence: "+seq
			error_to_file += "\n\nError Message: "+str(e)
			error_to_file += "\nLine number: "+str(exc_tb.tb_lineno)
			error_to_file += "\nError type: "+str(exc_type)
			error_to_file += "\nFunction name: "+ str(fname)
			error_to_file += "\n\nAlignment results for this Sequence:\n"
			error_to_file += json.dumps(query_results,sort_keys=True,indent=4, separators=(',', ': '))
			error_to_file += "\n**End of Error message**\n\n"
			ferrorlog.write(error_to_file)	
			error_dic = query_results			
			error_dic["Notes"] =  "Error"
			error_dic["Errors"] = "This query generated an error (see error log file): "+str(e)			
			
			if write_format == "TAB":								
				Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
			else:																									   				
				Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])			
																		
			continue
							
	foutfile.close()
	ferrorlog.close()
	f1.IFclass.close()
	resulting_files = [outfile,annotatedFile]
	if total_parsing_errors>0:
		print("The analysis has completed. However, {0} of the {1} total sequences analyzed contained errors when parsing the file. Please refer to the error log file to debug possible errors".format(str(total_parsing_errors), str(numSeqs)))
		resulting_files.append(outfile_errors)
	else:
		os.system("rm '{0}'".format(outfile_errors))		
		resulting_files.append('')		
	
	if unknown_hits == 0: 
		os.system("rm '{0}'".format(outfile_unknown_rec))
		resulting_files.append('')
	else:
		resulting_files.append(outfile_unknown_rec)	
	return resulting_files


def MakeQueryList(query, fields_list):
	"""
		Just simply convert lists in query to mongodb format 
	"""
	mongo_query = {}
	for field in fields_list:
		if field in query:
			if isinstance(query[field], list):
				mongo_query[field] = {'$in': query[field]}				
			elif isinstance(query[field], basestring):
				mongo_query[field] = query[field]
			else:
				raise Exception('Improper format for {0} query. We only allow a list of strings or a single string. Fields provided by user: {1}'.format(field, json.dumps(query)))
	return mongo_query	
	

def download_germlines_from_db(query_settings):
	"""
		Function for downloading germline locus sequence files from our MongoDB database		
		Parameters
		----------
		query_settings : dict
			dict has the following allowable keys: 
				_id: LIST OF IDS REFERRING TO SPECIFIC GERMLINE SETS YOU WANT TO DOWNLOAD (if this is passed in all other fields are ignored)
				SPECIES: LIST OF SPEICES OR SINGLE STRING 
				PRODUCTIVITY: STRINGS DEFINING PRODUCTIVITY OF GERMLINES TO QUERY
				LOCUS: LIST OF LOCUS OR SINGLE STRING	
				MOLECULAR_COMPONENT: LIST OR SINGLE STRING
				SOURCE: STRING DEFINING SOURCE
				VERSION: STRING DEFINING THE VERSION
		
		.. important::
			IF _id is not provided , then SPECIES IS REQUIRED 									
	"""

	query = copy.deepcopy(query_settings)

	# Create a special folder for germlines downloaded from database only
	output_directory = os.path.join(databaseFolder, 'databasedownloads/')
	
	if not os.path.isdir(output_directory):
		os.makedirs(output_directory)
		
	# Capitalize everything 
	ids = query.pop('_id', None)
	query = {k.upper(): v for k, v in query.iteritems()}
	
	if ids:
		query['_id'] = ids	
	if '_id' in query and query['_id']:
		id_list = query['_id']		
		query.pop('_id')
	elif 'SPECIES' in query:
		query.pop('_id', None)
		id_list = []
		# we will always choose the following database source and person as default if not defined
		source = 'IMGT'
		uploadedby = 'immunogrep'	
		temp = MakeQueryList(query, ['SPECIES', 'LOCUS', 'MOLECULAR_COMPONENT', 'SOURCE', 'VERSION', 'UPLOADED-BY', 'GENETYPE'])		
		query.update(temp)		
		if 'SOURCE' not in query:
			query['SOURCE'] = source
		if 'UPLOADED-BY' not in query:
			query['UPLOADED-BY'] = uploadedby										
	else:
		raise Exception('Either a list of germline ids or a SPECIES must be provided as query input. Fields provided by user: {0}'.format(json.dumps(query)))
	
	# Add extra filters/functionality for query only genes within a germline set with provided productivity 
	# i.e. add filter to the query if one was specifieid ('F','[F]','ORF')
	functionality = query.pop('FUNCTIONALITY', None)
	if not functionality:
		productive_gene_filters = {'$nin': []}
	else:
		if not isinstance(functionality, list):
			productive_gene_filters = {'$in': [functionality]}
		else:
			productive_gene_filters = {'$in': functionality}

	# create a query for getting germlines from database 
	db_class_var = query_germlines.GermlineDB()						
	
	# this is just for documenting purproses/keeping track of what the query request was. storing settings of query basically 	
	query.pop('GENETYPE', None)

	unique_settings = db_class_var.QueryDistinctValsByID(id_list, extra_filters=copy.deepcopy(query), distinct_fields=['SPECIES', 'GENETYPE', 'MOLECULAR_COMPONENT', 'SOURCE', 'VERSION', 'LOCUS'])					
	unique_settings['DB-FILE-SOURCE'] = unique_settings.pop('SOURCE')	
	if functionality:
		unique_settings['FUNCTIONALITY'] = functionality
	# Annoying, but make subtypes lowercase
	selected_genes = unique_settings['GENETYPE']
	
	# figure out a name for the database file that makes it easy to recognize when referred to later on 

	germlines = {}	
	if not selected_genes:
		print ('WARNING: the following settings returned no results from database {0}'.format(json.dumps(query_settings)))
	for subtype in ['V', 'J']:		
		if subtype in selected_genes:
			basenamefile = ('_'.join(unique_settings['SPECIES']) + '_' + '_'.join(unique_settings['LOCUS']) + '_' + subtype + '_' + re.sub('[\:_\- ]', '', str(datetime.now())) + '.txt').lower().replace(' ', '')
			query['GENETYPE'] = subtype			
			# actually run the query and download germlines as a TAB file for igfft format
			db_class_var.QueryGenesByID(id_list, extra_filters=copy.deepcopy(query), gene_functionality_filter=productive_gene_filters).PrintFFTDBFormat(parent_folder=output_directory, filename=basenamefile) 		
			# ensure sequences downloaded correctly 	
			if os.path.isfile(os.path.join(output_directory, basenamefile)):				
				germlines[subtype.lower()] = [os.path.join(output_directory, basenamefile)]
			else:								
				print ('WARNING: the following settings returned no results from database {0}'.format(json.dumps(query_settings)))
	
	return [germlines, unique_settings]


def CDR3_search(cdr3_search_class, algn_seq, cdr3start, end_of_ab, locus=None, max_lmotif_len=0, max_rmotif_len=0):
	
	guess_starting = max(0, cdr3start - max_lmotif_len)
	guess_ending = min(end_of_ab + 21, len(algn_seq))
	algn_seq = algn_seq[:guess_ending]
	
	[bestmotifset, MaxP, cdr3start, cdr3end, cdr3_nt, cdr3_aa, bestchain, bestscoreset, allscores] = cdr3_search_class.FindCDR3(algn_seq, suggest_chain=locus, start_pos=guess_starting, strand='+')	
	if cdr3_nt:		
		fr4start = cdr3end + 1
		cdr3Frame = cdr3start % 3 + 1
		fr4Frame = fr4start % 3 + 1
		result_notes = ''				
		fr4_nt = algn_seq[fr4start:end_of_ab + 1]				
		fr4_aa = TranslateSeq(fr4_nt, 0)		
	else:		
		result_notes = 'CDR3 not found because motif probability score was below threshold'
		cdr3start = -1
		fr4start = -1
		cdr3Frame = 0
		fr4Frame = 0
		fr4_nt = ''
		fr4_aa = ''
	
	return [cdr3_nt, cdr3_aa, fr4_nt, fr4_aa, cdr3Frame, fr4Frame, cdr3start, fr4start, result_notes]

		
def TranslateSeq(ntstring, frame):
	ntstring = ntstring[frame:]
	truncatedLength = 3 * (len(ntstring) / 3)
	ntstring = ntstring[:truncatedLength]
	return str(Seq(ntstring, generic_dna).translate())


	
####################################################################################################	

##########SUPPLMENTARY FUNCTIONS FOR DEFINING VARIABLE NAMES, WRITING, AND
##########READING IGBLAST DATA TO TEXT FILE########################


#this fucntion will define the order of how fields will appear in a text file
def TABFileHeader():
	chain_independent_results_order = ['Header',
									idIdentifier,
									'Sequence',
									'Quality_Score',
									'Document_Header',
									'Strand_Corrected_Sequence',
									'Notes',
									'Errors',
									'Command',
									'Recombination_Type',
									'Percent_Identity',
									'Alignment_Length',
				'Direction',				
				'Locus',
				'Chain',
				'Codon_Start',
				'Query_Start',
				'Query_End',
				'Codon_Frame',
				'5_Prime_Annotation',
				'3_Prime_Annotation',
				'Stop_Codon',
				'CDR3_Junction_In_Frame',
				'Full_Length',
				'Productive',
				'Full_Length_Sequence.NT',
				'Full_Length_Sequence.AA',				
				"Top_V-Gene_Hits",				
				"V-Gene_Alignment_Scores", 
				"Top_J-Gene_Hits", 
				"J-Gene_Alignment_Scores", 
				"FR1_Sequence.NT", 
				"CDR1_Sequence.NT", 
				"FR2_Sequence.NT", 
				"CDR2_Sequence.NT", 
				"FR3_Sequence.NT", 
				"CDR3_Sequence.NT", 
				"FR4_Sequence.NT", 
				"FR1_Sequence.AA", 
				"CDR1_Sequence.AA", 
				"FR2_Sequence.AA", 
				"CDR2_Sequence.AA", 
				"FR3_Sequence.AA", 
				"CDR3_Sequence.AA", 
				"FR4_Sequence.AA", 
				'FR1_Sequence.NT.Gapped',
				'CDR1_Sequence.NT.Gapped',
				'FR2_Sequence.NT.Gapped',
				'CDR2_Sequence.NT.Gapped',
				'FR3_Sequence.NT.Gapped',																		
				
				'FR1_Sequence.AA.Gapped',
				'CDR1_Sequence.AA.Gapped',
				'FR2_Sequence.AA.Gapped',
				'CDR2_Sequence.AA.Gapped',
				'FR3_Sequence.AA.Gapped',
				'VRegion.SHM.NT',
				'VRegion.SHM.Per_nt',
				"Reading_Frames: FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4",
				"VGENE: Total_Matches", 
				"VGENE: Total_Mismatches", 
				"VGENE: Total_Indel", 
				'JRegion.SHM.NT',
				'JRegion.SHM.Per_nt',
				"JGENE: Total_Matches", 
				"JGENE: Total_Mismatches", 
				"JGENE: Total_Indel", 
				"VGENE_Matches: FR1,CDR1,FR2,CDR2,FR3,CDR3", 
				"VGENE_Mismatches: FR1,CDR1,FR2,CDR2,FR3,CDR3", 
				"VGENE_Indels: FR1,CDR1,FR2,CDR2,FR3,CDR3", 				
				"VGENE: Query_Start", 
				"VGENE: Query_End", 
				"VGENE: Query_FR1_Start::End", 
				"VGENE: Query_CDR1_Start::End", 
				"VGENE: Query_FR2_Start::End", 
				"VGENE: Query_CDR2_Start::End", 
				"VGENE: Query_FR3_Start::End", 
				"VGENE: Germline_Start", 
				"VGENE: Germline_End", 
				"JGENE: Query_Start", 
				"JGENE: Query_End", 
				"JGENE: Germline_Start", 
				"JGENE: Germline_End", 
				"VGENE: Alignment_Sequence_Query", 
				"VGENE: Alignment_Sequence_Germline", 
				"VGENE: Alignment_FR1_Start::End",
				"VGENE: Alignment_CDR1_Start::End",
				"VGENE: Alignment_FR2_Start::End",
				"VGENE: Alignment_CDR2_Start::End",
				"VGENE: Alignment_FR3_Start::End",
				"JGENE: Alignment_Sequence_Query", 
				"JGENE: Alignment_Sequence_Germline",
				"Isotype",
				"Isotype mismatches",
				"Isotype percent similarity",
				"Isotype barcode direction"
	]
	
	igDict = OrderedDict()
	for i, field in enumerate(chain_independent_results_order):
		igDict[field] = i	
		
	return igDict


def DatabaseTranslator(input_dictionary={}):

	"""
		We will need to update the database with the results from IgFFT.  In order
		to update teh database, we need a translator, so that we know what fields go where in the database
	"""

	key = translation_var
	
	translator = {					
		"ANNOTATION": "GEORGIOU_INHOUSE",  # NAME OF THE ANALYSIS 
		"RECOMBINATION_FIELD": {  # THIS TELLS THE PROGRAM HOW TO DETERMINE WHETHER AN ANALYSIS/QUERY RESULT (from file) IS VDJ OR VJ
			"FIELD_NAME": "Recombination_Type",  # name of the field in the file that will give information regarding the recombination type (VDJ OR VJ)
			"EXPLICIT": True,   # IF EXPLICIT, then the RECOMBINATION_TYPE is equal to the EXACT value in this field (i.e. it will either list VDJ or VJ), IF false, then VDJ and VJ are defined by values in list below
			"INEXPLICIT_DEFINITIONS": {
				"VDJ": [],  # if expliity is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VDJ type (i.e. if Locus is used to determine recombination type then it woudl be VDJ:['IGH']
				"VJ": [],  # if expliity is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VJ type (i.e. if Locus is used to determine recombination type then it woudl be VJ:['IGK','IGL']
			}
		},
		"FIELDS_UPDATE": {  
			# key = field name  in database 
			# value = field name in file
			# this will map all of the fields in the file to the proper location in the database. I.E. If I list VGENES as the column name/field name, then i want to map VREGION.VGENES:VGENES (because VREGION.VGENES is the name in the database)						
			idIdentifier: idIdentifier,
			"SEQUENCE": "Sequence",
			"COMMAND": "Command",
			"SEQUENCE_HEADER": "Header",	
			"QUALITY_SCORE": "Quality_Score",
			"NOTES": "Notes",
			"PREDICTED_AB_SEQ.NT": "Full_Length_Sequence.NT",
			"PREDICTED_AB_SEQ.AA": "Full_Length_Sequence.AA",
			"STRAND": "Direction",
			"PREDICTED_CHAIN_TYPE": "Chain",
			"PRODUCTIVE": "Productive",
			"LOCUS_NAME": "Locus",
			"FULL_LENGTH": "Full_Length",
			"STOP_CODONS": "Stop_Codon",
			"VREGION.SHM.NT": "VRegion.SHM.NT",
			'VREGION.VGENE_QUERY_START': 'VGENE: Query_Start',
			'VREGION.VGENE_QUERY_END': 'VGENE: Query_End',
			"VREGION.FR1.NT": "FR1_Sequence.NT",
			"VREGION.FR1.AA": "FR1_Sequence.AA",
			"VREGION.CDR1.NT": "CDR1_Sequence.NT",
			"VREGION.CDR1.AA": "CDR1_Sequence.AA",			
			"VREGION.FR2.NT": "FR2_Sequence.NT",
			"VREGION.FR2.AA": "FR2_Sequence.AA",
			"VREGION.CDR2.NT": "CDR2_Sequence.NT",
			"VREGION.CDR2.AA": "CDR2_Sequence.AA",
			"VREGION.FR3.NT": "FR3_Sequence.NT",
			"VREGION.FR3.AA": "FR3_Sequence.AA",
			"VREGION.VGENES": "Top_V-Gene_Hits",
			"VREGION.VGENE_SCORES": "V-Gene_Alignment_Scores",
			"CDR3.NT": "CDR3_Sequence.NT",
			"CDR3.AA": "CDR3_Sequence.AA",
			"DREGION.DGENES": "Top_D-Gene_Hits",
			"DREGION.DGENE_SCORES": "D-Gene_Alignment_Scores",
			"JREGION.FR4.NT": "FR4_Sequence.NT",		
			"JREGION.FR4.AA": "FR4_Sequence.AA",
			"JREGION.JGENES": "Top_J-Gene_Hits",
			"JREGION.JGENE_SCORES": "J-Gene_Alignment_Scores",
			'JREGION.JGENE_QUERY_START': 'JGENE: Query_Start',
			'JREGION.JGENE_QUERY_END': 'JGENE: Query_End',																
			"ISOTYPE.GENE": "Isotype",
			"ISOTYPE.MISMATCHES": "Isotype mismatches",
			"ISOTYPE.PER_ID": "Isotype percent similarity"		
		}
	}
	
	input_dictionary[key] = translator			
	return input_dictionary


def PrintErrorsToScreen(filename):
	'''
		This function will be used to print to screen the error file generated when parsing program
	'''
	eof = False
	line_string = ""
	with open(filename) as f1:
		while not(eof):
			line = f1.readline()
			if line == "":
				eof = True
			else:				
				line_string += line
				if line.strip() == "**End of Error message**":
					yield line_string
					line_string = ""	
	yield line_string


def EnsureDatabaseDirectoryExists():
	if not os.path.isdir(databaseFolder):
		os.makedirs(databaseFolder)
		os.makedirs(os.path.join(databaseFolder, 'database'))
		
	for subtype in ['V', 'D', 'J']:
		folder = os.path.join(databaseFolder, subtype)
		if not os.path.isdir(folder):		
			os.makedirs(folder)			


def GetDefaultParameters():
	"""
		Get default settings for the igfft program
	"""		
	try:
		# This parameter creates a file containing default parameters
		subprocess.call(igfft_location + ' --defaults', shell=True)
		dataparameters = open('defaultsettings_fftprogram.txt').read().split('\n')	
		default_settings = {}				
		for p in dataparameters:
			if not p.strip():
				continue
			p = p.strip().split('\t')
			default_settings[p[0][1:]] = {}									
			if p[1] != "none":
				v = p[1].split('.')
				v[0] = v[0].lstrip('0')
				if len(v) > 1:
					v[-1] = v[-1].rstrip('0')
				v = '.'.join(v)
				v = v[:-1] if v[-1] == '.' else v
				default_settings[p[0][1:]] = v
			else:
				default_settings[p[0][1:]]['default'] = None					
	except:		
		default_settings = {
			'num_hits': None,
			'gap_open_fft': None,
			'gap_extend_fft': None,
			'pep_len': None,
			'gap_extend_sw': None,
			'gap_open_sw': None,
			's_pep': None,
			'cluster_per_id': None,
			'group_clusters': None,
			'gap': None,
			'similar_clusters': None,
			'match_sw': None,
			'mismatch_sw': None,
			'min_per_id': None,
			'ratio_cutoff': None,
			'times_above_ratio': None,
			's_fft': None						
		}
	os.remove('defaultsettings_fftprogram.txt')

	return default_settings


def Write_Seq_TAB(seqDic, output_fields, foutfile):
	foutfile.write('\t'.join([str(seqDic[field]) for field in output_fields]) + '\n')	
	

def GrabAdditionalHeaderInfo(header):
	tmp = header.split(fasta_file_delimiter)
	
	if len(tmp) > 1:
		additional_info = json.loads(tmp[1])
	else:
		additional_info = {}
	additional_info.pop('document_header', None)  # pop out this from dictionary if it was carried along somehow
	additional_info['Document_Header'] = header
	additional_info['Header'] = tmp[0]
	return additional_info


