'''
	Module for generating a FASTA file to use as the reference databsae in our proteomics pipeline.
	Requires the following unix/cygwin tools: awk
	Requires the following binary: GG lab clonotyping tool (written in D)
	Written by: Sebastian Schaetzle
	Date: 10-19-2015
'''
import sys
sys.path.append("C:\Users\cc35463\OneDrive\Appsoma Database\GitHub\IGREP\common_tools")

import os
# for making folders in igrep
import immunogrep_file_system_tools as filesystem
# for reading files
import immunogrep_read_file as readfile
import copy
from collections import defaultdict
import datetime
import re
import subprocess

import pandas as pd
# import numpy as np
# import gc
# from itertools import izip
from immunogrep_global_variables import descriptor_symbol as comment_line

allowed_file_formats = ['TAB', 'CSV', 'JSON', 'IMGT']
# Order in which fields should appear in the FASTA file. Any fields not defined here will appear afterwards as optional fields
# important: the first element in default_order_of_fields will be treated as the 'sequence' in a FASTA file >header\nsequence!!!
default_order_of_fields = ['ABSEQ.AA', 'RECOMBINATIONTYPE', 'DBIDENTIFIER','READS', 'TAG', 'VREGION.CDR1.AA', 'VREGION.CDR2.AA', 'CDR3.AA', 'VREGION.VGENES', 'JREGION.JGENES', 'Clonotype']
rename_must_be_present = {'CDR1': 'VREGION.CDR1.AA', 'CDR2': 'VREGION.CDR2.AA', 'CDR3': 'CDR3.AA'}
# variable for determining whether the antibody is heavy or light chain
chain_call = {
	'VDJ': ['TRB', 'TB', 'IGH', 'TRD'],
	'VJ': ['IGK', 'TRA', 'TA', 'IGL', 'TRG'],
}

ms_version = 1

# Only these characters will be allowed in amino acid sequences
aa_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWZXY'
pattern = r"[^" + aa_alphabet + "]"


def default_field_names():
	'''
		Defining a default set of field names to use for the program. Assume default set is from the database.
		Therefore a default input file for this program would have the following field names defined in any order:
		'PREDICTED_AB_SEQ', 'VREGION.CDR1.AA', 'CDR3.AA', 'VREGION.SHM.NT_PER', 'JREGION.SHM.NT_PER',
		'VREGION.VGENES', 'JREGION.JGENES', 'RECOMBINATION_TYPE', 'ISOTYPE.GENE'
	'''
	default_fields = {
		'ABSEQ.AA': 'PREDICTED_AB_SEQ.AA',
		'VREGION.CDR1.AA': 'VREGION.CDR1.AA',
		'VREGION.CDR2.AA': 'VREGION.CDR2.AA',		
		'CDR3.AA': 'CDR3.AA',
		'VREGION.VGENES': 'VREGION.VGENES',
		'JREGION.JGENES': 'JREGION.JGENES'
	}
	return default_fields


def iffile_to_pandas(ifclass, chunks):
	'''
		Convert files to dataframes using our read class
	'''
	list_seqs = []
	chunker = 0
	for i in ifclass.read():
		list_seqs.append(i)
		chunker += 1
		if chunker == chunks:
			yield pd.DataFrame(list_seqs)
			list_seqs = []
			chunker = 0
	if not list_seqs:
		yield list_seqs
	ifclass.IFclass.close()


def cdr_in_seq(row, cdr_field):
	return row[cdr_field] in row['ABSEQ.AA']


def tab_to_fasta(tabfile, fastafilename, ofs='_', append=False):
	'''
		Converts tab file created by module/pd.to_csv into a fasta file
	'''
	method = ">" if append is False else ">>"
	h = 1 if append is False else 0
	
	output_path = fastafilename  # os.path.join(os.path.dirname(tabfile), fastafilename)
	awk_convert = ''' awk 'BEGIN{{FS="\t"; OFS="{ofs}"; header={h} }}
		{{seq=$1; $1=""}}
		int(header)==1&&FNR==1{{print "#"substr($0,2) {method} "{output}"}}
		FNR==1{{next}}
		{{print ">"substr($0,2) {method} "{output}"; print seq {method} "{output}" }}' {input}
	'''.format(ofs=ofs, output=output_path, input=tabfile, method=method, h=h)
	error = subprocess.call(awk_convert, shell=True)
	if error != 0:
		raise Exception('Error occurred when converting to FASTA file')


def pandas_read_chunks(input_file, filetype, fields, chunks=10000):
	'''
		Function for creating pandas dataframes using the provided input_files.
		We will read x lines from each file at a time where x = chunks
		If the input file is not a delimited file (not CSV, TAB), then we will have to read the file using our read class and create a dataframe from those results
	'''
	# FLIPPED FIELDS -> the parameter fields is as follows: key => field name we use in program, value => field name in provided file.
	# We will flip this structure to make it easier for parsing the input files
	flipped_fields = {value: key for key, value in fields.iteritems()}

	each_file = input_file
	if isinstance(each_file, list):
		print('Reading file: ', each_file[0])
	else:
		print('Reading file: ', each_file)
	# First Lets figure out how we will read the files. Will we be able to simply load the dataframe using a delimited file, or will we use our methods for reading the files
	if not filetype:
		# We need to guess the filetype
		temp = readfile.immunogrepFile(each_file)
		guessed_type = temp.getFiletype()
	else:
		guessed_type = filetype

	# Next lets figure out what fields from each file we want to load
	if guessed_type == 'IMGT':
		# we will need these fields from IMGT
		field_names = ['V-D-J-REGION_5', 'V-J-REGION_5', 'CDR1-IMGT_5', 'CDR2-IMGT_5', 'CDR3-IMGT_5', 'Functionality_1', 'V-GENE and allele_1', 'J-GENE and allele_1', 'V-REGION Nb of mutations_8']
	else:
		field_names = flipped_fields.keys()

	if guessed_type in ['TAB', 'CSV']:
		# this is easy to read into a pandas dataframe
		seps = {'TAB': '\t', 'CSV': ','}
		skip_lines = 0
		# We need to figure out if we have to skip any lines because we aren't using our immunogrep class reader
		for line in open(each_file):
			line = line.strip()
			if not line.startswith(comment_line):
				break
			skip_lines += 1
		tmp = readfile.immunogrepFile(each_file, 'TAB')
		header_names_in_file = tmp.getDescription()
		
		field_names = [f for f in field_names if f in header_names_in_file]
		tmp.IFclass.close()
		reader = pd.read_table(each_file, sep=seps[guessed_type], chunksize=chunks, usecols=field_names, dtype=object, skip_blank_lines=True, skiprows=skip_lines)
	else:
		# We need to use our class for reading files
		if guessed_type == 'IMGT':
			reader = iffile_to_pandas(readfile.immunogrepFile(each_file, 'IMGT', required_files=[1, 5, 8], field_names=field_names, chunk_size=chunks))
		else:
			reader = iffile_to_pandas(readfile.immunogrepFile(each_file, guessed_type, field_names=field_names, chunk_size=chunks))
	for table in reader:

		# Rename the columns names based on the translator provided (this df will have column names as we expect in program)
		table.rename(columns=flipped_fields, inplace=True)
		# append dataset name to column for input file

		if 'READS' not in table.columns:
			table['READS'] = 1
		else:
			table['READS'] = table['READS'].astype(int)
			table['READS'].fillna(1, inplace=True)
		if 'RECOMBINATIONTYPE' not in table.columns:
			table['RECOMBINATIONTYPE'] = ''
		table.fillna('', inplace=True)
		yield table


def mongodb_generate_msdb_file(expquery, productivity_filter=['PRODUCTIVE', 'PRODUCTIVE (SEE COMMENT)', 'MAYBE', 'YES'], additional_filter={}, dbidentifier=None, output_folder_path=None, cluster_id=0.9, must_be_present=['CDR1', 'CDR2', 'CDR3'], use_vl_sequences=False):
	'''
		Main function for generating a FASTA mass spec db file data from our GGLAB mongo db. User passes in a query to define which experiments from the database are desired.
	'''
	pass


def generate_msdb_file(input_files, filetype=None, output_folder_path=None, dbidentifier=None, dataset_tags=None, translate_fields={}, cluster_id=0.9, must_be_present=['CDR1', 'CDR2', 'CDR3'], use_vl_sequences=False, pc_file_location=None, vl_file_location=None):
	'''
		Main function for generating a FASTA mass spec db file using annotation. User passes in a list of input files generated from an annotation program,
		and this function will output a FASTA file that can be used as a reference db for identifying peptides from mass spec. It will also generate a summary file
		describing which filters were passed.

		Parameters
		----------
		input_files : list of strings
			corresponding to filepaths for each annotated file.
			.. note::IMGT filepath format
				If the input files come from IMGT analyses, then input_files can either be a list of lists (i.e. 11 files per experiment) or a single list of all experiment filenames
				(i.e. program will split files into proper experiments)
		filetype : string, default None
			Describes the filetype of the input files. file type can be either TAB, CSV, FASTA, FASTQ, IMGT
		output_folder_path : string
			Filepath for returning results. If not defined, will make a new folder
		dbidentifier : string, default None
			A string identifier to label the mass spec database file
		dataset_tags : list of strings; default None
			If defined, then this will be used as the identifier for each file/dataset. If Not defined (dataset_tags is None) then the experiment identifier will be equal to the filepath
			names defined by input_files.
			.. important::
				If defined, the length of the string must be equal to the length of the input files/each unique dataset provided
			.. note::
				The following characters are replaced from the tags: '_', ' ',':', '|' and ',' are replaced by '-'
		translate_fields : dict, default {}
			This defines which fields in the file corresponds to fields we need for the analysis.
			key = field we use in the program, value = field name in the provided file(s)
			If empty, then this variable will assume the field names are the exact same as the fields we use in this program (see default_fields)
		cluster_id : float, default 0.9
			The percent identity required for clonotyping
		must_be_present : list of strings, default ['CDR1','CDR2','CDR3']
			This will define which CDR fields MUST be present in the sequence to be considered for analysis.
			.. note::Fields
				Only CDR1, CDR2, and CDR3 can be defined in this list
			.. note::CDR3
				CDR3 will ALWAYS be required for this analysis
		use_vl_sequences : boolean, default False
			If True, then any VL sequences detected in the provided experiments will be appended to the database file.
			If False, then the program will append a default VL sequence list generated by Sebastian and stored in the database

		Outputs
		-------
		Path of the FASTA db file created
	'''
	# SETUP: Lets setup the function and make sure all provided settings are corect #####################################3333
	runtime = str(datetime.datetime.now()).replace(' ', '').replace('/', '').replace('-', '').replace(':', '').split('.')[0]

	# THE FOLLOWING IS USED TO CONFIRM INPUT SETTINGS AND FORMATS#
	if cluster_id <= 0 or cluster_id > 1:
		# The use provided an inaccurate number for clonotyping
		raise Exception('Cluster cutoff must be a floating number between 0 and 1')

	# User just passed in a single string for single file, then make it a list
	if not isinstance(input_files, list):
		input_files = [input_files]

	if isinstance(filetype, basestring):
		# user just passed in a string. copy this to all a list containing same string for all files
		annotated_file_formats = filetype.upper()
		if filetype not in allowed_file_formats:
			raise Exception('The provided file format, {0}, is not currently supported in this program'.format(annotated_file_formats))
		annotated_file_formats = [annotated_file_formats] * len(input_files)
	else:
		raise Exception('Only a single string can be used to define the filetypes for all input files; If this innacurate please notify administrator')

	# We dont HAVE to do this, but lets just do it...
	if filetype == 'IMGT':
		# MAKE sure the input files are in the proper IMGT groups of 11 files each
		imgt_groups = readfile.GroupIMGTFiles(input_files)
		input_files = [vals for group, vals in imgt_groups.iteritems()]

	# No dataset_tags were provided so make them equal to input files
	if not dataset_tags:
		dataset_tags = [os.path.basename(f[0]) for f in input_files] if filetype == 'IMGT' else [os.path.basename(f) for f in input_files]
	else:
		if len(dataset_tags) != len(input_files):
			raise Exception('If defining each dataset using tags, then you must define a tag for each provided input file')

	remove_from_tags = ['\_', '\ ', '\|', '\,']
	for dn in range(len(dataset_tags)):
		dataset_tags[dn] = re.sub(re.compile('|'.join(remove_from_tags)), '-', dataset_tags[dn])

	# Defining output file path
	if output_folder_path:
		# user provided ouput folder
		output_folder_path = output_folder_path.rstrip('/\\').rstrip(os.sep)
		if not os.path.isdir(output_folder_path):
			# Folder path provided does not exist
			# Make folder, but ONLY USE BASENAME. FORCE the folder to exist within IGREP project
			parent_dir = filesystem.ExperimentDirs().make_dir(os.path.basename(output_folder_path))
			print("WARNING, OUTPUT DIRECTORY NOT FOUND. CREATED A NEW DIRECTORY FOR OUTPUT CALLED: {0}".format(parent_dir))
		else:
			# this is the output folder
			parent_dir = output_folder_path
	else:
		# User did not pass in output folder, so make a new folder
		folder_name = 'MSDB_' + runtime
		# make a new folder
		parent_dir = filesystem.ExperimentDirs().make_dir(folder_name)
	# replace empty spaces in folder name (can be annoying..)
	if not dbidentifier:
		dbidentifier = 'MASSSPECDB' + runtime
	# Make sure dbidentifier is properly formated for FASTA output
	dbidentifier = re.sub(re.compile('|'.join(remove_from_tags)), '-', dbidentifier)
	parent_dir = parent_dir.replace(' ', '_')
	prefix = dbidentifier
	prefix = os.path.join(parent_dir, prefix)

	if not translate_fields:
		# User did not define the field names, so lets assume it matches the default file format
		fields = {key: value for key, value in default_field_names().iteritems()}
	else:
		# User has defined field names to use in the program
		fields = {key: value for key, value in translate_fields.iteritems()}

	# Check if user passed in a value other than CDR1, 2, OR 3 for must_be_present field
	if 'CDR3' not in must_be_present:
		must_be_present.append('CDR3')

	unexpected_values = set(must_be_present) - set(['CDR1', 'CDR2', 'CDR3'])

	if unexpected_values:
		raise Exception('User passed in an unexpected set of values for the parameter "must_be_present". Allowed values are "CDR1", "CDR2", and/or "CDR3". User included other fields: {0}'.format(','.join(list(unexpected_values))))

	must_be_present = [rename_must_be_present[c] for c in must_be_present]
	# THESE KEYS ARE REQUIRED TO RUN THE PROGRAM. ALL OTHER FIELDS ARE OPTIONAL
	required_fields_for_program = ['ABSEQ.AA', 'VREGION.VGENES', 'JREGION.JGENES'] + must_be_present

	fields_undefined = [fieldname for fieldname in required_fields_for_program if fieldname not in fields or not(fields[fieldname])]
	if fields_undefined:
		# User did not define all required fields
		raise Exception('You must define the fieldname in your input file for the following fields we need in the program: {0}'.format(','.join(fields_undefined)))
	# **************SETTINGS COMPLETED********************* #

	# ************PROCEED TO ACTUAL FUNCTION/STEPS ************* #
	# Step 1: parse the input fields
	print('Parsing provided files and generating list of unique sequences')
	[vdj_df, vj_df, stats] = parse_annotation_files(input_files, filetype, fields, must_be_present, dataset_tags, dbidentifier, use_vl_sequences)
	if vdj_df is None:
		print('Problem occurred: there is no VDJ/VH data, please double check input files are proper; Also check the summary file')
	else:
		print('Clonotyping heavy chains')
		clonotype_seqs(vdj_df, cluster_id, prefix + 'msDB.fasta', prefix + 'malGucken.txt', append=False)
	if not(vj_df is None):
		# append VJ data to files
		print('Clonotyping light chains')
		clonotype_seqs(vj_df, cluster_id, prefix + 'msDB.fasta', prefix + 'malGucken.txt', append=True)
	print('Writing summary file')
	write_summary_file(input_files, prefix + 'parsed_summary.txt', stats)

	# CONCATENATE DATABSE FILES

def ProcessGene(gene, noallele=True):
	"""
		Removes the allele calls from genes		
		Assumption: 			
			Multiple genes are seperated by ',' 			
			Alleles are seperated by '*' 			
			If a gene is seperated by multiple spaces, then the gene should be identified by a gene that contains either - or '*' 
			For example: 
				Imgt genes may be: 
					Homo sapiens IGHV1-3*01
					We only want to isolate the word IGHV1-3
	"""	
	if not gene:
		return ''
	# extract first genes in field 
	top_gene = gene.split(',')[0]	
	# split each gene by spaces. Go through each word
	split_words  = top_gene.split(' ')	
	if len(split_words)==1: #there is only one word 
		g = split_words[0]#.split('*')[0]
	else:
		g = ''
		for subv in split_words:
			#if we find a word with gene characters in it, it must be our gene 
			if '*' in subv or '-' in subv:						
				#extract everything before '*'
				g = subv
				
				break
	if noallele:
		top_gene = g.split('*')[0]
	else:
		top_gene = g
	return top_gene
		


def clean(df):
	'''
		Process a few fields from the data frame and clean them/remove annoying aspects
	'''
	df['ABSEQ.AA'] = df['ABSEQ.AA'].map(lambda x: x + 'ASTK')
	df['ABSEQ.AA'] = df['ABSEQ.AA'].str.replace('_', '')
	df['CDR3.AA'] = df['CDR3.AA'].str.replace('_', '')
	# => choose first gene only

	df['VREGION.VGENES'] = df['VREGION.VGENES'].map(lambda x: ProcessGene(x) )  # x.split(',')[0].split('*')[0])  # because germans dislike regular expressions
	df['JREGION.JGENES'] = df['JREGION.JGENES'].map(lambda x: ProcessGene(x, False) )  # x.split(',')[0].split('*')[0])
	if 'ISOTYPE.GENES' in df.columns:
		df['ISOTYPE.GENES'] = df['ISOTYPE.GENES'].map(lambda x: ProcessGene(x) )  # x.split(',')[0].split('*')[0])
	if 'DREGION.DGENES' in df.columns:
		df['DREGION.DGENES'] = df['DREGION.DGENES'].map(lambda x: ProcessGene(x) )  # x.split(',')[0].split('*')[0])
	return df


def determine_recomb(row):
	'''
		Fxn for determining if a sequence is VDJ (vh) or VJ (vl/k)
	'''
	recomb_type = ''
	if row['RECOMBINATIONTYPE']:
		recomb_type = row['RECOMBINATIONTYPE']
	else:
		# we need to determine recombination type using v gene or j gene
		vgene = row['VREGION.VGENES']
		jgene = row['JREGION.JGENES']
		if (vgene and not recomb_type) or (recomb_type != 'VDJ' and recomb_type != 'VJ'):
			# We need to predict recobmination using vgene
			if vgene[:3] in chain_call['VDJ']:
				recomb_type = 'VDJ'
			elif vgene[:3] in chain_call['VJ']:
				recomb_type = 'VJ'
		if (jgene and not recomb_type) or (recomb_type != 'VDJ' and recomb_type != 'VJ'):
			# We need to predict recobmination using jgene
			if jgene[:3] in chain_call['VDJ']:
				recomb_type = 'VDJ'
			elif jgene[:3] in chain_call['VJ']:
				recomb_type = 'VJ'
	return recomb_type


def filter_dataframe(df, filtered_events, must_be_present):
	'''
		Filter each dataframe based on the following rules
			1) must have an antibody sequence
			2) cannot have a stop codon
			3) cannot have unrecognizeable characters
			4) all cdrs defined by must_be_present must be found
	'''

	num_seqs = len(df)
	# Remove sequences that lack a full length antibody sequence
	df = df[df['ABSEQ.AA'] != ""]
	new_len = len(df)
	filtered_events['missing_seq'] += num_seqs - new_len
	num_seqs = new_len
	if num_seqs == 0:
		return None

	# Remove sequences that have '*' codons in sequence
	df = df[~df['ABSEQ.AA'].str.contains(r"\*")]
	new_len = len(df)
	filtered_events['has_stop_codons'] += num_seqs - new_len
	num_seqs = new_len
	if num_seqs == 0:
		return None

	df = clean(df)

	# Remove sequences that have unrecognizeable characters
	df = df[~df['ABSEQ.AA'].str.contains(pattern)]
	new_len = len(df)
	filtered_events['unk_chars'] += num_seqs - new_len
	num_seqs = new_len

	cdr_bool = df['CDR3.AA'] != ''
	for cdr in must_be_present:
		if cdr == 'CDR3.AA':
			# we already did this filter above
			continue
		cdr_bool = (cdr_bool) & (df[cdr] != '')

	# Remove rows that are missing required CDRs
	df = df[cdr_bool]
	new_len = len(df)
	filtered_events['missing_cdr'] += num_seqs - new_len
	num_seqs = new_len

	if num_seqs == 0:
		return None

	# Remove rows whose CDRs are out of frame
	for cdr in must_be_present:
		df = df[df.apply(lambda x: x[cdr] in x['ABSEQ.AA'], axis=1)]
	new_len = len(df)
	filtered_events['cdr_not_in_frame'] += num_seqs - new_len
	num_seqs = new_len
	if num_seqs == 0:
		return None

	# Make sure RECOMBINATION_TYPE is correct
	df['RECOMBINATIONTYPE'] = df.apply(determine_recomb, axis=1)

	# Remove recombination types that are not 'VDJ' or 'VJ'
	df = df[df['RECOMBINATIONTYPE'].isin(['VDJ', 'VJ'])]
	new_len = len(df)
	filtered_events['unk_recomb'] += num_seqs - new_len
	num_seqs = new_len
	if num_seqs == 0:
		return None
	return df


def merge_tag_fxn(unique_tag):
	IDs = []
	for i in unique_tag:
		if i not in IDs:
			IDs.append(i)
	return '|'.join(list(set(IDs)))


def clono(row, clonoDict):
	row['Clonotype'] = clonoDict[row['CDR3.AA']]
	return row


def merge_df(df):
	'''
		fxn for collapsing dataframes by full length antibody sequence
	'''
	groupby_cols = ['RECOMBINATIONTYPE', 'ABSEQ.AA']
	agg_commands = {col: 'first' for col in df.columns if col not in groupby_cols}  # we will only want to keep the first instance of a column during group by
	agg_commands['READS'] = 'sum'  # we will want to sum all the reads while aggregating

	"""
	#if large_df is not None:
	#	df = large_df.concat([df], axis=0)
		# dfShort = df
	# else:
		# dfMerged = pd.merge([large_df, df], axis=0)
	"""
	return df.groupby(groupby_cols).agg(agg_commands).reset_index()


def merge_all_data(dflist):
	'''
		fxn for concating all dataframes in list and then collapsing by unique antibody sequence,
		and choosing the unique tag/id for each unique sequence
	'''
	groupby_cols = ['RECOMBINATIONTYPE', 'ABSEQ.AA']
	dfMerged = pd.concat(dflist, axis=0).reset_index()
	agg_commands = {col: 'first' for col in dfMerged.columns if col not in groupby_cols}  # we will only want to keep the first instance of a column during group by
	agg_commands['READS'] = 'sum'  # we will want to sum all the reads while aggregating
	agg_commands['TAG'] = merge_tag_fxn

	dfShort = dfMerged.groupby(groupby_cols).agg(agg_commands)
	return dfShort


def parse_annotation_files(input_files, filetype, fields, must_be_present, dataset_tags, dbidentifier, use_vl_sequences):
	'''
		Method for parsing annotation files that are not in a delimited file format and filtering out sequences which do not match our requirements for mass spec file
		Using pandas we will:
		1) Perform the following filters:
			1) Remove sequences with stop codons
			2) Remove sequences that do not have CDRs defined by must_be_present
			3) Remove sequences whose CDRs are not in frame in the full length
			4) Only select VH sequences (VDJ)
		2) Will group sequences by unique amino acid in each input file (WILL NOT COLLAPSE BY ALL because sebastian does that)

		.. note::
			If use_vl_sequences is true, then a seperate file for VJ sequences will be made according to steps above

		.. note::
			This method should only be called by generate_msdb_file
	'''
	chunks = 100000
	filtered_events = defaultdict(int)
	df_list = []
	df_vdj = None
	df_vj = None
	df_alldata = None
	for file_num, each_file in enumerate(input_files):
		# Create a generator for reading input files as dataframe
		pandas_reader_generator = pandas_read_chunks(each_file, filetype, fields, chunks=chunks)

		# iterate..
		total_seqs = 0
		for df in pandas_reader_generator:
			num_seqs = len(df)
			total_seqs += num_seqs
			
			# Filter rows based on our rules for good seqs
			df = filter_dataframe(df, filtered_events, must_be_present)
			if df is None:
				# No results passed filtering
				continue
			filtered_events['passed_filter'] += len(df)
			r_counts = df['RECOMBINATIONTYPE'].value_counts()

			if 'VDJ' in r_counts:
				filtered_events['VDJ'] += r_counts['VDJ']
			if 'VJ' in r_counts:
				filtered_events['VJ'] += r_counts['VJ']

			if use_vl_sequences is False:
				df = df[df['RECOMBINATIONTYPE'] == 'VDJ']

			# Merge and group results by full length antibody sequence
			df['TAG'] = dataset_tags[file_num]
			df_list.append(merge_df(df))

			if (total_seqs) / 100000 > (total_seqs - chunks) / 100000:
				print('Processed {0} sequences'.format(total_seqs))
	print('Merging alldata')
	if not df_list:
		return [None, None, filtered_events]
	df_alldata = merge_all_data(df_list).reset_index()
	df_alldata['DBIDENTIFIER'] = dbidentifier
	print('Sortingdata')
	df_alldata = df_alldata.sort(['READS'], ascending=False)
	df_vdj = df_alldata[df_alldata['RECOMBINATIONTYPE'] == 'VDJ']
	if use_vl_sequences and filtered_events['VJ'] > 0:
		df_vj = df_alldata[df_alldata['RECOMBINATIONTYPE'] == 'VJ']

	df_alldata = None
	filtered_events['total_seqs'] = total_seqs
	return [df_vdj, df_vj, filtered_events]


def clonotype_seqs(dfShort, percentIdentity, dbfilename, tabfilename, append=False):
	'''
		This is our function for going through the results file and clonotyping sequences by their CDRH3.
		Assumptions:
		df => a dataframe that has already been collapsed by unique sequences and all the TAGS are properly formated and filters passed		
	'''

	'''		ok, for now we'll start with a dictionary (dfDict) of DFs and the keys are IDs, like "HD1_D273_VH"
	remember to leave/set numbers as integers/floats, not string...'''

	table_output = tabfilename
	fasta_output = dbfilename

	# Generate alist of unique CDR3s
	print('Generating unique CDR3 list')
	cdr3List = dfShort['CDR3.AA'].unique()
	# not used yet
	# read_count_by_cdr3 = dfShort.groupby(['CDR3.AA']).agg({'READS': 'sum'})
	print('Running clonotyping')
	# Call clonotyping
	# CDR3s = cdr3List
	sample_output_file = os.path.join(os.path.dirname(table_output), 'cdr3list.txt')
	with open(sample_output_file, 'w') as writecdr3:
		for cd in cdr3List:
			writecdr3.write(cd + '\n')
	output_file_cdr3_result = os.path.join(os.path.dirname(table_output), 'cdr3_clonotype_list.txt')
	run_clonotype_command = 'CDR3_Clonotyping --file {0} --thresh {1} --output {2}'.format(sample_output_file, str(1 - percentIdentity), output_file_cdr3_result)
	error = subprocess.call(run_clonotype_command, shell=True)
	if error != 0:
		raise Exception('Error occurred when running clonotyping')
	"""
	# Generate clonotyping dict
	# clonotyped CDR3s are in clonoDict
	#clonoDict = {}
	#with open(output_file_cdr3_result) as clono_file:
	#	for clonotypes_created in clono_file:
	#		if not clonotypes_created.strip():
	#			continue
	#		clonotypes_created = clonotypes_created.strip().split('\t')
	#		clonoDict[clonotypes_created[0]] = clonotypes_created[1]cd
	#clonoDict = pd.DataFrame(clonoDict)
	"""
	clonoDict = pd.read_csv(output_file_cdr3_result, sep='\t', skip_blank_lines=True, header=None, names=['CDR3.AA', 'Clonotype'])

	# Create a clonotype field in the dataframe
	# dfClono = dfShort.apply(lambda x: clono(x), axis=1)  # => instead of this, lets use a dataframe merge function
	print('Appending clono')
	dfClono = pd.merge(dfShort, clonoDict, on='CDR3.AA', how='left')
	dfClono.fillna('', inplace=True)

	# NOW lets figure out how we should reorder columns in dataframe based on the order we want them to appear in the FASTA file
	current_columns = list(dfClono.columns)
	output_fields = copy.deepcopy(default_order_of_fields)
	for col in current_columns:
		if col not in output_fields:
			# add a new column from dataframe to our output
			output_fields.append(col)
	# Now make sure all columns exist in dataframe
	for col in output_fields:
		if col not in current_columns:
			dfClono[col] = ''
	
	# Now reorder columns
	dfClono = dfClono[output_fields]
	# Finally RENAME COLUMNS AND ENSURE NO COLUMN NAMES HAVE A '_' IN IT; we do this because the FASTA file is delimited by '_'
	renamecol = {r: r.replace('_', '') for r in dfClono.columns}
	dfClono.rename(columns=renamecol, inplace=True)
	dfClono.set_index('ABSEQ.AA', inplace=True)
	dfClono.drop('index', axis=1, inplace=True)
	# new_field_order = list(dfClono.columns)

	# Export DF
	mode = 'a' if append is True else 'w'
	header = False if append is True else True
	print('Exporting MS DB as TAB')
	dfClono.to_csv(table_output, sep='\t', mode=mode, header=header)

	print('Exporting MS DB as FASTA')
	"""
	# For some reason its super super slow to export dataframe this way. is this because it is not using c?
	dfClono = dfClono.astype(str)
	with open(fasta_output, mode) as output:
		if mode == 'w':
			output.write('#' + '_'.join(new_field_order) + '\n')
		for row in dfClono.iterrows():
			header = row[1].str.cat(sep='_')  # '_'.join(row[1].astype(str))
			# output.write('>'+row[1]['ID']+'_'+row[1]['CDR1_Sequence.AA']+'\n'+row[0]+'\n')
			output.write('>{0}\n{1}\n'.format(header, row[0]))
	"""
	# Try doing it with AWK instead
	tab_to_fasta(table_output, fasta_output, "_", append=append)


def write_summary_file(input_files, summary_file, filtered_events):
	total_seqs = filtered_events['total_seqs']
	passed_filter = filtered_events['missing_seq']
	missing_seq = filtered_events['passed_filter']
	missing_cdr = filtered_events['missing_cdr']
	cdr_not_in_frame = filtered_events['cdr_not_in_frame']
	has_stop_codons = filtered_events['has_stop_codons']
	unk_chars = filtered_events['unk_chars']
	unk_recomb = filtered_events['unk_recomb']
	VH_count = filtered_events['VDJ']
	VL_count = filtered_events['VJ']

	# MAKE THIS SPECIAL IF ITS IMGT FILES (MERGE FILES INTO ONE GROUP NAME)	
	file_path_summary = ['\t' + f[0] + ',' if isinstance(f, list) else '\t' + f + ',' for f in input_files]

	summary = open(summary_file, 'w')
	# Write useful information in summary file
	summary.write('Data generated using version %s on: %s\n' % (ms_version, datetime.datetime.now()))
	summary.write("***************************************************\n")
	summary.write("***************   Summary Report  *****************\n")
	summary.write("***************************************************\n\n\n")

	summary.write('The following input files were used to make MS DB:\n')
	summary.write('\n'.join(file_path_summary) + '\n\n\n')

	summary.write('Parsed through {0} sequences\n'.format(str(total_seqs)))
	summary.write('{0} ({1}%) sequences did not pass filters:\n'.format(str(total_seqs - passed_filter), str(round(float((100 * (total_seqs - passed_filter)) / total_seqs), 1))) if total_seqs > 0 else '0')
	summary.write('\t{0} sequences did not have an antibody sequence\n'.format(str(missing_seq)))
	summary.write('\t{0} sequences did not have a required cdr\n'.format(str(missing_cdr)))
	summary.write('\t{0} sequences had stop codons\n'.format(str(has_stop_codons)))
	summary.write('\t{0} sequences had unrecognizeable characters\n'.format(str(unk_chars)))
	summary.write('\t{0} sequences had unrecognizeable recombination types\n'.format(str(unk_recomb)))
	summary.write('\t{0} sequences had an out of frame CDR\n\n'.format(str(cdr_not_in_frame)))

	summary.write('{0} ({1}%) sequences were VDJ recombination:\n\n'.format(str(VH_count), str(round(float(100 * (VH_count) / total_seqs), 1))) if total_seqs > 0 else '0')
	summary.write('{0} ({1}%) sequences were VJ recombination:\n\n'.format(str(VL_count), str(round(float(100 * (VL_count) / total_seqs), 1))) if total_seqs > 0 else '0')

	summary.close()


"""
# These functions are deprecated but can be useful in future (?)
def parse_annotation_files(input_files, filetype, outfile_prefix, fields, must_be_present, dataset_tags, hardcoded_fields, optional_fields, dbidentifier, use_vl_sequences):
	'''
		Method for parsing annotation files that are not in a delimited file format and filtering out sequences which do not match our requirements for mass spec file
		1) Will perform the following filters:
			1) Remove sequences with stop codons
			2) Remove sequences that do not have CDRs defined by must_be_present
			3) Remove sequences whose CDRs are not in frame in the full length
			4) Only select VH sequences (VDJ)
		2) Will group sequences by unique amino acid in each input file (WILL NOT COLLAPSE BY ALL because sebastian does that)

		.. note::
			If use_vl_sequences is true, then a seperate file for VJ sequences will be made according to steps above

		.. note::
			This method should only be called by generate_msdb_file
	'''
	has_stop_codons = 0
	missing_seq = 0
	missing_cdr = 0
	cdr_not_in_frame = 0
	VH_count = 0
	VL_count = 0
	total_seqs = 0
	unk_chars = 0
	unk_recomb = 0
	passed_filter = 0
	alphabet='ABCDEFGHIJKLMNOPQRSTUVWZXY'
	pattern = re.compile('[^'+alphabet+']')
	
	vh_file_name = outfile_prefix+'parsed_vh.txt'
	vl_file_name = outfile_prefix+'parsed_vl.txt'
	
	summary_file_name = outfile_prefix + 'parsed_summary.txt'	
	hardcode_copy = copy.deepcopy(hardcoded_fields)
	hardcode_copy.pop('READS')
	# hardcode_copy.pop('TAG')		
	# fields we will use for the output file
	output_header_row = hardcode_copy.keys() + optional_fields.keys()	
	# keys in files (harcoded fields have a key value format key = field name in new file, value = field in provided file)
	field_order = hardcode_copy.values() + optional_fields.values()		
	#vh_file.write('\t'.join(output_header_row) + '\n')
	#vl_file.write('\t'.join(output_header_row) + '\n')	
	# remove these two fields because we add them in later using the OrderedDicts		
	# Read through all of the input files
	all_seq_data = []
	
	vh_file = open(vh_file_name+'.temp1', 'w')
	vl_file = open(vl_file_name+'.temp1', 'w')
	vh_file2 = open(vh_file_name+'.temp2', 'w')
	vl_file2 = open(vl_file_name+'.temp2', 'w')
	
	for exp_num, each_file in enumerate(input_files):				
		if filetype != 'IMGT':
			iffile = readfile.immunogrepFile(each_file, filetype)
		else:
			iffile = readfile.immunogrepFile(each_file, 'IMGT', required_files=[1,5,8],field_names=['V-D-J-REGION_5', 'V-J-REGION_5', 'CDR1-IMGT_5', 'CDR2-IMGT_5', 'CDR3-IMGT_5', 'Functionality_1', 'V-GENE and allele_1', 'J-GENE and allele_1', 'V-REGION Nb of mutations_8'])		
				
		# keep track of unique amino acid sequences from VDJ	
		unique_aa_vdj = OrderedDict()
		# keep track of unique amino acid sequences from VJ
		unique_aa_vj = OrderedDict()
			
		# Read through all sequences in file
		for row in iffile.read():						
			if (total_seqs ) % 100000 == 0:
				print('Processed {0} sequences'.format(total_seqs))
			#if (total_seqs + 1) % 500000 == 0:
			#	break
			if not row:
				continue				
			row = defaultdict(str, row)
			total_seqs += 1
			if not row[fields['ABSEQ.AA']]:
				# no full length antibody sequence
				missing_seq += 1
				continue			
			row[fields['ABSEQ.AA']] = row[fields['ABSEQ.AA']].replace('_','')
			row[fields['CDR3.AA']] = row[fields['CDR3.AA']].replace('_','')
			full_seq = row[fields['ABSEQ.AA']]						
			if '*' in full_seq:
				# has a stopcodon
				has_stop_codons += 1
				continue			
			if pattern.search(full_seq):
				#has unknown characters
				unk_chars += 1
				continue						
			# Makes a list of elements to check if CDR is present and CDR is found in full length
			# Only CDRs present in row will be added to the list
			cdr_lists = [row[fields[cdr]] in full_seq for cdr in must_be_present if row[fields[cdr]]]			
			# If a CDR is missing, then it will not be added to the list so the length of the final list will be < the number of CDRs checked (must_be_present)
			if len(cdr_lists) < len(must_be_present):
				missing_cdr += 1
				continue
			# If all CDRS are found within full length then the sum of True, True, True (etc) will be = to the length of the list
			if sum(cdr_lists) < len(cdr_lists):
				cdr_not_in_frame += 1
				continue			

			recomb_type = row[fields['RECOMBINATIONTYPE']]
			vgene = ProcessGene(row[fields['VREGION.VGENES']])			
			jgene = ProcessGene(row[fields['JREGION.JGENES']])
			if (vgene and not recomb_type) or (recomb_type != 'VDJ' and recomb_type != 'VJ'):				
				# We need to predict recobmination using vgene
				if vgene[0][:3] in chain_call['VDJ']:
					recomb_type = 'VDJ'
				elif vgene[0][:3] in chain_call['VJ']:
					recomb_type = 'VJ'
			if (jgene and not recomb_type) or (recomb_type != 'VDJ' and recomb_type != 'VJ'):
				# We need to predict recobmination using jgene
				if jgene[0][:3] in chain_call['VDJ']:
					recomb_type = 'VDJ'
				elif jgene[0][:3] in chain_call['VJ']:
					recomb_type = 'VJ'
			if recomb_type != 'VDJ' and recomb_type != 'VJ':
				unk_recomb += 1
				continue
			row[fields['VREGION.VGENES']] = vgene
			row[fields['JREGION.JGENES']] = jgene
			passed_filter += 1
			# Now add in some fields we will use in program
			# row['READS'] = 1
			row['TAG'] = dataset_tags[exp_num]
			row['RECOMBINATIONTYPE'] = recomb_type
			row['DBIDENTIFIER'] = dbidentifier
			if 'ISOTYPE.GENES' in fields:
				row[fields['ISOTYPE.GENES']] = ProcessGene(row[fields['ISOTYPE.GENES']])						
			if recomb_type == 'VDJ':
				VH_count += 1												
				if full_seq not in unique_aa_vdj:
					# Store the unique amino acid sequence
					unique_aa_vdj[full_seq] = 0
					output_row = [ str(row[f1]) for f1 in field_order ]					
					vh_file.write('\t'.join(output_row) + '\n')				
				# Update count
				unique_aa_vdj[full_seq] += 1
			else:
				VL_count += 1
				if use_vl_sequences is False:
					# No need to store sequence info
					continue				
				if full_seq not in unique_aa_vj:
					unique_aa_vj[full_seq] = 0
					output_row = [ str(row[f1]) for f1 in field_order ]					
					vl_file.write('\t'.join(output_row) + '\n')				
				unique_aa_vj[full_seq] += 1
			
		for seq, val in unique_aa_vdj.iteritems():
			vh_file2.write(str(val)+'\n')
		
		for seq, val in unique_aa_vj.iteritems():
			vl_file2.write(str(val)+'\n')
				
	output_header_row.extend(['READS'])
	print('Reporting sequence counts')
	
	vh_file.close()
	vh_file2.close()
	vl_file.close()
	vl_file2.close()
	
	# Paste two files together	
	with open(vh_file_name, 'w') as vh_file, open(vh_file_name +'.temp1') as f1, open(vh_file_name +'.temp2') as f2:
		vh_file.write( '\t'.join(output_header_row) + '\n' )		
		for line1, line2 in izip(f1, f2):			
			vh_file.write("{}\t{}\n".format(line1.rstrip('\n'), line2.strip()))
		

	# Sort file by the number of reads (ignore first line in sort )
	#rint('Sorting file')
	#my_folder = os.path.dirname(vh_file_name)
	#ead_col = len(output_header_row)
	#rror = subprocess.call('''(head -n 1 "{0}" && tail -n+2 "{0}" | sort -T "{1}" -t '\t' -k{2}nr,{2} ) > "{0}.sorted"  '''.format(vh_file_name,my_folder, read_col), shell=True)
	#if error != 0:
	#	raise Exception('Error when processing subprocess terminal unix sort command')			
	#os.rename(vh_file_name+'.sorted', vh_file_name)			
	
	unique_aa_vdj = {}
	unique_aa_vj = {}
				
	if VL_count == 0 or use_vl_sequences == False:
		vl_file_name = None
	else:
		print('Saving VJ data')
			# Paste two files together	
		with open(vl_file_name, 'w') as vl_file, open(vl_file_name +'.temp1') as f1, open(vl_file_name +'.temp2') as f2:
			vl_file.write( '\t'.join(output_header_row) + '\n' )		
			for line1, line2 in izip(f1, f2):			
				vl_file.write("{} {}\n".format(line1.rstrip('\n'), line2.strip()))

		unique_aa_vj = {}		

	
	print('Parsing complete')
	os.remove(outfile_prefix+'parsed_vl.txt.temp1')		
	os.remove(vh_file_name+'.temp1')		
	write_summary_file(input_files, summary_file_name, total_seqs, passed_filter, missing_seq, missing_cdr, cdr_not_in_frame, has_stop_codons, unk_chars, unk_recomb, VH_count, VL_count)		
	gc.collect()
	gc.collect()
	return [vh_file_name, vl_file_name, summary_file_name]



def group_unique_seqs(filepath, seq_col, tag_col, read_col):
	'''
		This function will collapse an input file by unique sequences. 
		It will also keep track of all unique 'tags' in a unique sequence group.
		We assume that the input filepath was generated by the function 'parse_annotation_files'
		
		.. note::
			We could have probably improved speed by grouping unique sequences whil parsing the annotation files at the same time.
			But we instead wanted to focus on just not creating hash tables of unique sequences, and instead usng this function to 
			group sequences together
	'''
	# Now sort the input file by antibody sequence and by tag name
	print('Sorting by antibody sequence')
	my_folder = os.path.dirname(filepath)	
	subprocess.call('''tail -n+2 "{0}" | sort -T "{1}" -t '\t' -k{2},{2} -k{3},{3} > "{0}.sorted" '''.format(filepath,my_folder, seq_col + 1, tag_col + 1), shell=True)
	
	with open(filepath) as r:
		header = r.readline()
		
	# Now parse through the sorted file and group by unique amino acid sequence
	print('Grouping by AA seq')
	print tag_col
	with open(filepath, 'w') as out:		
		file_buffer = open(filepath+'.sorted')
		firstline = file_buffer.readline().split('\t')		
		ab_seq = firstline[seq_col]
		tag = [ firstline[tag_col] ]
		reads = int(firstline[read_col])
		dataline = firstline
		total_reads = 0		
		for line in file_buffer:
			line = line.split('\t')			
			ab_line = line[seq_col]
			tag_line = line[tag_col]
			read_line = int(line[read_col])
			if ab_seq != ab_line:
				# new sequence, output prevoius row
				dataline[read_col] = str(reads)
				dataline[tag_col] = '|'.join(tag)
				out.write('\t'.join(dataline) + '\n')
				total_reads += reads
				
				reads = read_line
				tag = [ tag_line ]
				ab_seq = ab_line
				dataline = line
			else:
				if tag_line not in tag:
					tag.append(tag_line)
				reads += read_line
		# report last row
		dataline[read_col] = str(reads)
		dataline[tag_col] = '|'.join(tag)
		out.write( '\t'.join(dataline) + '\n')
		total_reads += reads	
	
	print total_reads
	# Sort sequences by antibody counts		
	print('Sorting by Ab seq count')
	my_folder = os.path.dirname(filepath)	
	header_row_text = '\t'.join(header) + '\n'
	
	subprocess.call('''echo -e "{3}" | cat sort -T "{1}" -t '\t' "{0}" -k{2}nr,{2} > "{0}.sorted" '''.format(filepath, my_folder, read_col + 1, header_row_text), shell=True)
	os.rename(filepath+'.sorted', filepath)
"""
	

