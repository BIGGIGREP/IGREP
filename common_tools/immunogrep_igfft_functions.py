#VERSION 360-370 SOMEWHERE IN THIS VERSION WE CHANGED HOW RECOMBINATION IS DISPLAYED 

##################################### PARSE IGFFT RESULT FILE ####################################################################################
#THIS SCRIPT WILL PARSE THE ALIGNMENT FILE CREATED BY RUNNING IGFFT VERSION 0.6 and above
#The code is hardcoded and assumes that the output of the file is identical to the following template:
#this is sample output from the igblast program using the settings defined above
##############################################################################################

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import commands
import re

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import time
from time import gmtime, strftime

import immunogrep_useful_immunogrep_functions as useful

import json
import os
import traceback
from pprint import pprint

#import immunogrep_run_igblast_command as run_igblast
import immunogrep_immunogrepfile as readwrite

from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import fasta_file_delimiter
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict
import collections
import immunogrep_appsoma_cdr3tools as CDR3tools
from immunogrep_global_variables import descriptor_symbol #this is the symbol we will use to designate 'comment' lines in the file 
from immunogrep_global_variables import translation_var #this key will signify the translation/translator key 
from immunogrep_global_variables import filetype_var
import immunogrep_query_germline_functions as query_germlines #use this for querying germlines 
import cchrysostomou_nt_fft_align_tools as fftaligner
from cchrysostomou_isotype_fft import defaultBarcodes
import copy
#from cchrysostomou_readwritefiles import *



#This function will be used to print to screen the error file generated when parsing program
def PrintErrorsToScreen(filename):	
	eof = False
	line_string = ""
	with open(filename) as f1:
		while not(eof):
			line = f1.readline()
			if line == "":
				eof=True
			else:
				#line=line.strip('\r\n')
				line_string+=line
				if line.strip()=="**End of Error message**":
					yield line_string
					line_string = ""	
	yield line_string

databaseFolder = "scratch/fftannotation_db_files/"



def EnsureDatabaseDirectoryExists():
	if not os.path.isdir(databaseFolder):
		os.makedirs(databaseFolder)
		os.makedirs(databaseFolder+'/database')
		
	for subtype in ['V','D','J']:
		folder = databaseFolder+subtype
		if not os.path.isdir(folder):		
			os.makedirs(folder)			

def GetDefaultParameters():
	
	try:
		os.system('igfft --defaults')
		filename_created = 'defaultsettings_fftprogram.txt'
		with open(filename_created) as f:
			dataparameters = f.readlines()
		os.system('''rm '{0}' '''.format(filename_created))
		default_settings = {}
		
		for p in dataparameters:
			p = p.strip().split('\t')
			default_settings[p[0][1:]] = {}
			
			if p[1]!="none":
				v = p[1].split('.')
				v[0] = v[0].lstrip('0')
				if len(v)>1:
					v[-1] =v[-1].rstrip('0')
				v = '.'.join(v)
				v = v[:-1] if v[-1] == '.' else v
				default_settings[p[0][1:]]['default'] = v
			else:
				default_settings[p[0][1:]]['default'] = None
			
			default_settings[p[0][1:]]['description'] = p[2] if len(p) > 2 else ''
	except:		
		default_settings = {
			'num_hits':None,
			'gap_open_fft':None,
			'gap_extend_fft':None,
			'pep_len':None,
			'gap_extend_sw':None,
			'gap_open_sw':None,
			's_pep':None,
			'cluster_per_id':None,
			'group_clusters':None,
			'gap':None,
			'similar_clusters':None,
			'match_sw':None,
			'mismatch_sw':None,
			'min_per_id':None,
			'ratio_cutoff':None,
			'times_above_ratio':None,
			's_fft':None						
		}
		
	return default_settings

#just simply convert lists in query to mongodb format 
def MakeQueryList(query,fields_list):
	mongo_query = {}
	for field in fields_list:
		if field in query:
			if isinstance(query[field],list):
				mongo_query[field] = {'$in':query[field]}				
			elif isinstance(query[field],basestring):
				mongo_query[field] = query[field]
			else:
				raise Exception('Improper format for {0} query. We only allow a list of strings or a single string. Fields provided by user: {1}'.format(field,json.dumps(query)))
	return mongo_query	
	
#function for downloading germline files from our database 
#the structure of query_settings is as follows: 
#allowable keys: 
			#_id: LIST OF IDS REFERRING TO SPECIFIC GERMLINE SETS YOU WANT TO DOWNLOAD (if this is passed in all other fields are ignored)
			#IF _id is not provided , then SPECIES IS REQUIRED 									
			#SPECIES: LIST OF SPEICES OR SINGLE STRING 
			
			#the following are optional, if not provided, then SOURCE IS ASSUMED TO BE IMGT AND VERSION IS THE LATEST VERSION 
			#PRODUCTIVITY: STRINGS DEFINING PRODUCTIVITY OF GERMLINES TO QUERY 
			#LOCUS: LIST OF LOCUS OR SINGLE STRING 										
			#MOLECULAR_COMPONENT: LIST OR SINGLE STRING 
			#SOURCE: STRING DEFINING SOURCE 
			#VERSION: STRING DEFINING THE VERSION 
def download_germlines_from_db(query_settings):
	query = {}
	
	output_directory = databaseFolder+'databasedownloads/'
	
	if not os.path.isdir(output_directory):
		os.makedirs(output_directory)
		
	#capitalize everything 
	ids = query_settings.pop('_id',None)
	query_settings = {k.upper():v for k,v in query_settings.iteritems()}
	if ids:
		query_settings['_id'] = ids
	
	
	if '_id' in query_settings and query_settings['_id']:		
		#query['_id'] = {'$in':query_settings['_id']}
		id_list = query_settings['_id']		
		query_settings.pop('_id')
		
	elif 'SPECIES' in query_settings:
		query_settings.pop('_id',None)
		id_list = []
		#we will always choose the following database source and person as default if not defined
		source ='IMGT'
		uploadedby='immunogrep'
				
		query = MakeQueryList(query_settings,['SPECIES','LOCUS','MOLECULAR_COMPONENT','SOURCE','VERSION','UPLOADED-BY','GENETYPE'])
		if 'SOURCE' not in query:
			query['SOURCE'] = source
		if 'UPLOADED-BY' not in query:
			query['UPLOADED-BY'] = uploadedby								
		
	else:
		raise Exception('Either a list of germline ids or a SPECIES must be provided as query input. Fields provided by user: {0}'.format(json.dumps(query_settings)))
	
	#add extra filters/functionality for query only genes within a germline set with provided productivity 
		#i.e. add filter to the query if one was specifieid ('F','[F]','ORF')
	if 'PRODUCTIVITY' in query_settings:
		if not isinstance(query_settings['PRODUCTIVITY'],list):
			query_settings['PRODUCTIVITY'] = [query_settings['PRODUCTIVITY']]
		if query_settings['PRODUCTIVITY']:
			productive_gene_filters = {'$in':query_settings['PRODUCTIVITY']}
		else:
			productive_gene_filters={'$nin':[]}
	else:
		productive_gene_filters = {'$nin':[]}
					
	
	#create a query for getting germlines from database 
	db_class_var = query_germlines.GermlineDB()						
	
	#this is just for documenting purproses/keeping track of what teh query request was. storing settings of query basically 
	unique_settings = db_class_var.QueryDistinctValsByID(id_list,extra_filters=copy.deepcopy(query),distinct_fields=['SPECIES','GENETYPE','MOLECULAR_COMPONENT','SOURCE','VERSION','LOCUS'])					
	
	unique_settings['DB-FILE-SOURCE'] = unique_settings.pop('SOURCE')
	
	selected_genes = unique_settings['GENETYPE']	
	
	if 'PRODUCTIVITY' in query_settings:
		unique_settings['PRODUCTIVITY'] = query_settings['PRODUCTIVITY']	
	
	
	
	#figure out a name for the database file that makes it easy to recognize when referred to later on 
	name_data = defaultdict(list)
	for s in unique_settings['SPECIES']:
		for l in unique_settings['LOCUS']:
			name_data[s].append(l)	
	file_prefix = 'DatabaseDownload_'+'_'.join({s+'_'.join(v) for s,v in name_data.iteritems()}) #should create a file name of datbasedownload_species_all loci_speices_all loci.txt
	file_prefix=file_prefix.replace(' ','')			
		
	#actually run the query and download germlines as a TAB file for igfft format 
	db_class_var.QueryGenesByID(id_list,extra_filters=copy.deepcopy(query),gene_functionality_filter = productive_gene_filters).PrintFFTDBFormat(parent_folder=output_directory,filename=file_prefix) 
		
	#ensure sequences downloaded correctly 	
	selected_genes = unique_settings['GENETYPE']	
	germlines = {}	
	for subtype in ['V','D','J']:
		if subtype in selected_genes:					
			filepath = output_directory+file_prefix+'_'+subtype+'.txt'#the function QueryGenesById will create germline database files with these filesnames				
			print filepath
			germlines[subtype] = filepath if os.path.isfile(filepath) else None #make sure the file downloaded correctly							
	
	print germlines
	
	return [germlines,unique_settings]
	

		
		
	

#testSet -> location ofthe input file
#outfile -> desired location and name of the output file
#variable_parameters -> optional parameters for running igblast
#skipRunning -> dont run igblast, just generate a string command
#input_Filetype -> the file type of the input file (FASTA, TAB, CSV, ETC)
#header_field -> the field in the file that corresponds to the sequence header
#sequence_field -> the field in the file that corresponds to the seqeuence 

#germline parameter is as follows:
	#a two index tuple 
		#first element => string defining method for describing where germline came from: default, custom, db
			#default => use default files currently in folder 
			#custom => provide custom text files 
			#db => download germlines from database 
		#second element => defines the parameters for the germline file locations
			#if element 1 => default:
				#provide a two element tuple defining the species and the molecular componenet
					#allowed species: Mus musuculus, homo sapiens 
					#allowed molecular components: IG, TR 
			#if element 2=> provide a dictionary whose keys are {V,D, AND/OR J} and the values are the filepath for each germline file 
		#third element => download fields from database 
			#return a dictionary for querying database: 
				#allowable keys: 
					#_id: LIST OF IDS REFERRING TO SPECIFIC GERMLINE SETS YOU WANT TO DOWNLOAD (if this is passed in all other fields are ignored)
					#IF _id is not provided , then SPECIES IS REQUIRED 									
					#SPECIES: LIST OF SPEICES OR SINGLE STRING 
					
					#the following are optional, if not provided, then SOURCE IS ASSUMED TO BE IMGT AND VERSION IS THE LATEST VERSION 
					#PRODUCTIVITY: STRINGS DEFINING PRODUCTIVITY OF GERMLINES TO QUERY 
					#LOCUS: LIST OF LOCUS OR SINGLE STRING 										
					#MOLECULAR_COMPONENT: LIST OR SINGLE STRING 
					#SOURCE: STRING DEFINING SOURCE 
					#VERSION: STRING DEFINING THE VERSION 
def RunIgFFT(testSet, germlines=None,outfile=None, variable_parameters={}, input_filetype='FASTA',header_field='header',sequence_field='sequence',quality_field='phred'):	
	#first go through the variable germlines to figure out the germline files to use 
	default_germline_folder = 'scratch/fftannotation_db_files'
	
	if not germlines:
		#we assume that the user manually provided v and j parameters to variable parameters. therefore, this is a CUSTOM source 
		
		custom_details = {}
		v = variable_parameters.pop('v',None)
		j = variable_parameters.pop('j',None)
		if v:			
			custom_details['v']=v
		if j:			
			custom_details['j']=j
		if custom_details == {}:
			#ther user has not provided any information at all regarding a germline 
			raise Exception("You must provide at least one germline database. If you do not have a germline database then define the variable germlines to use default germline database in program or download germline from GG lab database. See documentation in function")
		else:
			germlines = ('custom',custom_details)
			
	
	#the possible ways of getting germlien files is from 'default' (use default files in folder) , 'custom' (provide your own files), or 'db' (database)
	if not isinstance(germlines,tuple):
		raise Exception('Germline definition must be a tuple of two elements. Element one is a string of either "custom","db", or "default". Element two defines details concerning the specific germlines desired')
	
	subtypes = ['v','j'] #IGNORE d FOR NOW  (igfft requires lowercalse subtypes)
	
	if germlines[0]=='custom':
		germline_settings = {'GERMLINE-SOURCE':'CUSTOM-FILE',
							 'PARAMS':{}
							 }							
		#the user is providing proper germline files for each 
		for keys,paths in germlines[1].iteritems():
			keys = keys.lower()
			if keys in subtypes and os.path.isfile(paths):				
				variable_parameters[keys.lower()] = paths				
				#germline_settings['PARAMS'][keys.lower()] = paths 
	elif germlines[0] == 'default':
		germline_settings = {'GERMLINE-SOURCE':'DEFAULT',
							 'PARAMS':{}
						 }	
		#the user wants to use the default germlines from the program 
		species = germlines[1][0].lower().replace(' ','')
		molcular_component = germlines[1][1].upper().replace(' ','')[:2]				
		if os.path.isdir(default_germline_folder+'/'+species+'/'+molcular_component):
			for s in subtypes:
				#annoying, but igfft requires lowercalse subtypes, but folder names are uppercased...
				s = s.upper()
				if os.path.isdir(default_germline_folder+'/'+species+'/'+molcular_component+'/'+s.upper()):
					#get teh file listed inside 										
					variable_parameters[s.lower()] =default_germline_folder+'/'+species+'/'+molcular_component+'/'+s.upper()+'/'+os.listdir(default_germline_folder+'/'+species+'/'+molcular_component+'/'+s.upper())[0]  #there shoudl only be one file 										
					#germline_settings['PARAMS'][s.lower()] = variable_parameters[s.lower()]
	elif germlines[0] == 'db':
		#the user wants to download germlines from database		
		germline_settings = {'GERMLINE-SOURCE':'DATABASE-DOWNLOAD'}		
		[germline_files, germline_settings['PARAMS']] = download_germlines_from_db(germlines[1])
		for s in subtypes:
			s = s.upper()
			if s in germline_files:
				#get teh file listed inside 				
				#annoying, but igfft requires lowercalse subtypes
				variable_parameters[s.lower()] = germline_files[s]	 										
	else:		
		raise Exception('Germline definition must be a tuple of two elements. Element one is a string of either "custom","db", or "default". Element two defines details concerning the specific germlines desired')						
		
	if type(variable_parameters) is not dict:
		raise Exception('Incorrect format of variable_parameters. The correct format of the function variable_parameters is a dictionary. keys correspond to the parameters used in the program and values correspond to the settings for that parameter')
		
	
	#make sure igfft parameters are good in variable_parameters dictionary  
	default_parameters = GetDefaultParameters()
	
	default_parameters = {v:default_parameters[v]['default'] for v in default_parameters}
	
	
	EnsureDatabaseDirectoryExists()		
	
	num_found = 0
	remove_germline = []
	for subset in ['v','d','j']: 
		#remove any germlines whose path are not defined or are incorrect
		if subset in variable_parameters:
			if not(variable_parameters[subset]) or not(os.path.isfile(variable_parameters[subset])):
				variable_parameters.pop(subset)
		else:
			num_found+=1
				
	if num_found == 0:
		raise Exception('Could not find any of the germline files passed in the parameters. Program terminating.')	
		
			
	if not outfile:
		outfile = testSet+'.igfft.alignment'
			
	variable_parameters['o'] = outfile
	ignore_seq_prefix = {var:var[1:] if var[0] in ['v','d','j'] and len(var)>1 else var for var in variable_parameters } #this will remove subtype preffix for parameters (for example, if parameter is num_gaps, but we are only modifying vgene, then parameter will be vnum_gaps. this will remove the v part of the string)	
	
	remove_vals_set_as_default = [var for var in variable_parameters  if ignore_seq_prefix[var] in default_parameters and default_parameters[ignore_seq_prefix[var]] and variable_parameters[var]==default_parameters[ignore_seq_prefix[var]] ] #remove any paramters that are already set as default in the program
	
	for var in remove_vals_set_as_default:
		variable_parameters.pop(var)
			
	
	deletetemp = False
	#FASTA/FASTQ FILE FORMAT IS ALREADY SUPPORTED BY IGFFT
	if input_filetype=='FASTA':
		tempFile = testSet
		filetype='FASTA'
	elif input_filetype=='FASTQ':
		tempFile = testSet
		filetype='FASTQ'
	else:	
		deletetemp = True
		#read in the input file and convert it to a TAB file so that that can be read by igfft
		inputIFFile = readwrite.immunogrepFile(filelocation=testSet,filetype=input_filetype)
		date = str(datetime.now())
		remove_chars = [':','-',' ','.']
		for char in remove_chars:
			date = date.replace(char,'')		
		tempFile=testSet+'_{0}'.format(date)#temporary TAB file we will make for running program 
		output_temp_file = open(tempFile,'w')
		
		output_temp_file.write('HEADER\tSEQUENCE\tSEQUENCE_QUALITY\n')
		for line in inputIFFile.read():		
			if line:			
				db_info = '{0}{1}'.format(fasta_file_delimiter,json.dumps({idIdentifier:line[idIdentifier]})) if idIdentifier in line else ''
				quality=line[quality_field] if quality_field in line else ''
					
				if sequence_field in line:													
					output_temp_file.write('{0}{2}\t{1}\t{3}\n'.format(line[header_field].strip(),line[sequence_field].strip(),db_info,quality))
				
		output_temp_file.close()
		inputIFFile.IFclass.close()
		filetype='TAB'
	
	variable_parameters['i'] = filetype
	
	
	print('\n\nRunning IgFFT using the following settings:\n')
	pprint(variable_parameters)		
	running_program_text = "Input file location: {0}\n".format(tempFile)
	running_program_text += "Output file will be saved as: {0}\n".format(outfile)
	running_program_text += "IgFFT Run has started at {0}\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()))			
	
	#NOW RUN THE BINARY 
	print(running_program_text)
	command_string = '''igfft "{0}" '''.format(tempFile) #important! use " " to accoutn for spaces in filename!. 
	print command_string
	for var in variable_parameters:
		if var in ['o','v','d','j']:
			#OUTPUT FILE ADD DOUBLE QUOTES
			command_string += '-{0} "{1}" '.format(var, variable_parameters[var])			
		else:
			command_string += '-{0} {1} '.format(var, variable_parameters[var])				
	os.system(command_string)	
	
	#progrrm complete 
	print("Analysis Completed at {0}\nAnalysis saved to: {1}\n\n\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()),outfile))
		
	if deletetemp:
		#dleete any temp files that were made 
		os.system("rm '{0}'".format(tempFile))
				
	command = variable_parameters
	command.pop('o',None) #dont store output path 
	command.pop('i',None) #dont store input path 
	if 'v' in command:
		command['v'] = os.path.basename(command['v'])
	if 'j' in command:
		command['j'] = os.path.basename(command['j'])
	if 'd' in command:
		command['d'] = os.path.basename(command['d'])
	command['germline-details'] = germline_settings	
	command['annotation']='IGFFT'
	return command


def IdentifyCDR3UsingRegExp(algn_seq,cdr3start,end_of_ab,j_start,fr3present,jpresent,cdr3_search_parameters,locus,var_type,Jgermline_alignments,maxleftmotif=None,maxrightmotif=None):
	result_notes = ''
	fr4start = 0
	if fr3present:
		if jpresent:
			cdr3_fr4_seq = algn_seq[cdr3start:end_of_ab+1]				
			#Determine whether we are searching for a motif or not
			motif = cdr3_search_parameters[var_type] if cdr3_search_parameters and var_type in cdr3_search_parameters and cdr3_search_parameters[var_type] else None								
			if motif:						
				[fr4start,result_notes] = Find_CDR3_Motif(motif,Jgermline_alignments,cdr3_fr4_seq,cdr3start)			
			else:					
				fr4start = cdr3start + 3 * int((j_start - cdr3start + 1) / 3)					
			
			if fr4start != -1 and fr4start > cdr3start:
				cdr3_nt = algn_seq[cdr3start:fr4start - 1]		
				cdr3_aa = TranslateSeq(cdr3_nt,0)# Seq(cdr3_nt,generic_dna).translate().tostring()
				
				fr4_nt = algn_seq[fr4start:end_of_ab+1]
				fr4_aa = TranslateSeq(fr4_nt,0)#Seq(fr4_nt,generic_dna).translate().tostring()
						
				fr4Frame = (fr4start) % 3 + 1		
			else:
				if j_start > cdr3start:
					cdr3_nt = algn_seq[cdr3start:j_start] 
					cdr3_aa = TranslateSeq(cdr3_nt,0)#Seq(cdr3_nt,generic_dna).translate().tostring()				
				else:
					cdr3_nt = ""
					cdr3_aa = ""								
				fr4start = -1
				fr4_nt = ""
				fr4_aa = ""
				fr4Frame = 0
		else:
			result_notes='No JGene alignment was identified to find CDR3'
			fr4Frame = 0			
			cdr3_nt=''
			fr4_nt=''
			cdr3_aa=''
			fr4_aa=''
			fr4start = -1
		cdr3Frame = cdr3start%3+1		
	else:
		cdr3_nt = ''
		cdr3_aa = ''		
		fr4_nt = ''
		fr4_aa = ''
		fr4Frame = 0
		cdr3Frame = 0
		cdr3start = -1
		fr4start = -1
		result_notes = 'No FR3 was identified to find CDR3'
	   
	return [cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,cdr3start,fr4start,result_notes]


#this is the function for finding the CDR3 using a motif.  There are two
#options for finding the motif.
#a) search through all of the possible J-GERMLINE hits for the motif of
#interest.  If we find any of the j-germline-gene alignments, then we will map
#the position back to the proper position in the sequence
#b) search the CDR3-FR4 subsequence (antibody after the FR3 more or less) for
#the motif
def Find_CDR3_Motif(motif,j_gene_algn_info,cdr3_fr4_subseq,cdr3start):
			
	result_notes = ''
	
	motif_found = False

	#first we loop through all possible j-germline alignments.  in each germline
	#alignment, we search for the motif int he germline
	for j_algn in j_gene_algn_info:
		j_gene_query = j_algn['query_seq'] #query alignment
		j_gene_germline = j_algn['germline_seq'] #germline alignemnt
		j_gene_start = j_algn['start']
		
		j_search = j_gene_query.replace('-','')
		j_search_germ = j_gene_germline.replace('-','')
			
		mtch = re.finditer(motif,j_search_germ,re.IGNORECASE) #first search germline for motif
		regcount = 0
		for loopreg in mtch:																							
			m = loopreg #store teh last occurrence of the motif
			regcount+=1
			
		if regcount > 0:
			motif_found = True
			
			germ_fr4start = m.start()
			fr4start = j_gene_start - 1
			
			i = 0
			counter_char = 0
			
			#The motif was found.  Now we want to map the motif found in the germline to
			#the proper position in the sequence
			while i != germ_fr4start:					
				if j_gene_query[counter_char] != '-':
					fr4start+=1	
				if j_gene_germline[counter_char] != '-':
					i+=1						
				counter_char+=1			
			fr4start+=1
			break #exit_loop
			
	if not(motif_found): #we couldnt find the motif in the germline, so we need to search the cdr3-fr4
					  #sequence
		mtch = re.finditer(motif,cdr3_fr4_subseq,re.IGNORECASE) #search for all instances of motif in the query sequence j-region
	
		regcount = 0		
		for loopreg in mtch:																							
			m = loopreg #store the last occurrence of the motif
			regcount+=1
		
		if regcount > 0: 
			motif_found = True
			fr4start = m.start() + cdr3start
			
			
	if not(motif_found):					
		fr4start = -1
		result_notes = 'RegExp not found;'
		
				
	return [fr4start,result_notes]
					

def IdentifyCDR3UsingGGJunctionalMotif(algn_seq,cdr3start,end_of_ab,cdr3_index,fr3present,jpresent,cdr3_search_parameters=None,locus=None,var_type=None,Jgermline_alignments=None,max_lmotif_len=0,max_rmotif_len=0):
			
	Motif=cdr3_search_parameters
	guess_starting =max(0,cdr3start-max_lmotif_len)
	guess_ending = min(end_of_ab+21,len(algn_seq))
	algn_seq=algn_seq[:guess_ending]
	
	 
	[bestmotifset,MaxP,cdr3start,cdr3end,cdr3_nt,cdr3_aa,bestchain,bestscoreset,allscores] = CDR3tools.FindCDR3(algn_seq,Motif,suggest_chain=locus,start_pos=guess_starting,strand='+',motif_type='AA')
	if cdr3_nt:		
		fr4start = cdr3end+1
		cdr3Frame = cdr3start%3+1
		fr4Frame = fr4start%3+1
		result_notes = ''
		fr4_nt = algn_seq[fr4start:end_of_ab+1]
		fr4_aa = TranslateSeq(fr4_nt,0)
	else:		
		result_notes = 'CDR3 not found because motif probability score was below threshold'
		cdr3start = -1
		fr4start = -1
		cdr3Frame = 0
		fr4Frame = 0
		fr4_nt = ''
		fr4_aa = ''
	
	return [cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,cdr3start,fr4start,result_notes]
	
#input_format = list of parametesr for the input file.  first value = filetype,
#second value = field name for header sequence, third value = field name for
#the sequence
def Parse_Alignment_File(annotatedFile,outfile=None,commandVal={},numSeqs=0,input_format=['FASTA','header','sequence',''],write_format='TAB',return_results_dict=False, annotation_settings={}):
		
	###JUST SETTING UP SETTINGS AND PARAMETERS####
	dna_alphabet = "ACTGUKMRYSWBVHDXN"
	bad_dna_char_pattern = re.compile('[^' + dna_alphabet + ']',re.IGNORECASE)
	
	print('Parsing and summarizing alignment file')	
	if not(outfile):
		outfile = useful.removeFileExtension(annotatedFile)+'.igfft.annotation'
		
	outfile_unknown_rec = outfile + ".unk_recombination"
	outfile_errors = outfile + ".igfft.error_log" #this file will be used to note any errors that occur while parsing the file
	
	if annotatedFile == outfile:
		os.system('''mv '{0}' '{0}.temp' '''.format(annotatedFile))
		annotatedFile+='.temp'
	
	input_file_type = input_format[0]
	header_name = input_format[1]
	seq_name = input_format[2]	
	quality_name = input_format[3] if len(input_format)>3 else ''
	
	#decides what min cutoff to use to consider a successful sequence match
	if 'min_cutoff' in annotation_settings:
		min_per_algn_cutoff = annotation_settings['min_cutoff']
	else:
		min_per_algn_cutoff = 0.4
	
	#minimum number of sequences that must align to germlines to be considered
	#antibody
	if 'min_len_cutoff' in annotation_settings:
		min_algn_len_cutoff = annotation_settings['min_len_cutoff']
	else:
		min_algn_len_cutoff = 100
	
	#remove insertions in sequence
	if 'remove_insertions' in annotation_settings:
		remove_insertions = annotation_settings['remove_insertions']
	else:
		remove_insertions = 0 #0 -> dont remove insertions, 1 -> remove insertions always, 2-> remove
						#insertions only if there is a stop codon
	
	
	#determine whether we will perform isotyping
	isotype_aligner = None
	if 'isotype' in annotation_settings:	
		if annotation_settings['isotype']['method'] == 'gg_inhouse':
			identifyIsotype = True
			if 'param' in annotation_settings['isotype']:
				iso_settings = annotation_settings['isotype']['param']
			else:
				iso_settings = {}
			if 'barcode-list' in iso_settings:			
				isotype_barcodes= iso_settings['barcode-list']
			else:				
				isotype_barcodes=copy.deepcopy(defaultBarcodes())												
			
			p_t = iso_settings['penalize_truncations'] if 'penalize_truncations' in iso_settings else True
			num_mismatch = iso_settings['maxmismatch'] if 'maxmismatch' in iso_settings else 2
			minimum_iso_alignment_length = iso_settings['min-len'] if 'min-len' in iso_settings else 15
			search_rc_isotype = iso_settings['iso_search_direction'] if 'iso_search_direction' in iso_settings else 0 #only consider forward direction of barcodes
			
			isotype_aligner = fftaligner.BarcodeAligner(isotype_barcodes,p_t,search_rc_isotype,num_mismatch,minimum_iso_alignment_length)
	else:
		identifyIsotype = False
		
	

		
	
	#determine which CDR3 function we will use to identify cdr3
	if 'cdr3_search' in annotation_settings:
		identifyCDR3 = True
		if annotation_settings['cdr3_search']['method'] == 'gg_inhouse':
			#format of querying for cdr3: 
				#list of tuples: 
					#element 1) => species 
					#element 2) => list of loci
				#(species, list of loci)
			
			#query database for proper motif:
			#generate a query for motis using the combination of species and loci using the 									
			
			germline_query = []
			for possible_combos in annotation_settings['cdr3_search']['param']:
				if not isinstance(possible_combos[1],list):					
					germline_query.append({'Species':possible_combos[0],'Locus':possible_combos[1]}) #using list(set to ensure only unique values in list 
				else:
					germline_query.append({'Species':possible_combos[0],'Locus':{'$in':list(set(possible_combos[1]))}})
			germline_query = {'$or':germline_query}
			
			#now run the query 
			cdr3_motif_class = query_germlines.CDR3MotifDB()
			annotation_settings['cdr3_search']['param'] = cdr3_motif_class.GetMotifForProgram(query=germline_query)
												
			CDR3_SEARCH_FUNCTION = IdentifyCDR3UsingGGJunctionalMotif	 #point the fucntion cdr3_seearch_function to the method defined for motifs
			cdr3_search_parameters = annotation_settings['cdr3_search']['param']
			unique_locus_sets = set([a[0].upper() for a in cdr3_search_parameters])
			max_l_motif_len = 4*max(sorted([int(k) for a in cdr3_search_parameters for k in a[1]])) #find the maximum length of the poossible motifs
			max_r_motif_len = 4*max(sorted([int(k) for a in cdr3_search_parameters for k in a[2]])) #find the maximum lenght of the possible right motifs
			cdr3_search_name = list(unique_locus_sets)
		elif annotation_settings['cdr3_search']['method'] == 'regexp':
			max_l_motif_len = None
			max_r_motif_len = None			
			unique_locus_sets = set(['VDJ','VJ'])
			CDR3_SEARCH_FUNCTION = IdentifyCDR3UsingRegExp #point the fucntion cdr3_search_function to the method defined for regular
												  #expression
			cdr3_search_parameters = annotation_settings['cdr3_search']['param']			
			cdr3_search_name = cdr3_search_parameters
		else:
			raise Exception('Only allowed values for cdr3_search method are: gg_inhouse or regexp')
	else:		
		annotation_settings['cdr3_search'] = {
				'method':'CDR3 analysis not selected',
				'param':None
			}

		identifyCDR3 = False
		cdr3_search_name = "None selected"
		print ('No parameters were provided for identifying the CDR3. CDR3 analysis will not be included')
				
  

	#CHECK IF IGBLAST FILE WAS ACTUALLY CREATED
	
	if not(os.path.isfile(annotatedFile)):			   
		raise Exception('ERROR: IGFFT FILE WAS NOT CREATED.PLEASE MODIFY SETTINGS AND/OR RE-RUN')		
	
	if numSeqs == 0:
		numSeqs = useful.file_line_count(annotatedFile)-1						   

	folderOutputLocation = os.path.dirname(annotatedFile)
			
	foutfile = open(outfile,'w')
	foutfile_unknown_recombination = open(outfile_unknown_rec,'w')
	ferrorlog = open(outfile_errors,'w')
	output_files = {'annotated':foutfile,					
					'UNK':foutfile_unknown_recombination,
					'ERROR':ferrorlog}
	
	
	
	
	
	
	translator = DatabaseTranslator()
	translator[filetype_var] = write_format 
	translator_comment = descriptor_symbol#textFileCommentTranslator
	translator_string = json.dumps(translator)
	foutfile.write(translator_comment+translator_string+'\n')#in the first line of the file output a comment line definining how to convert this file into a database file 
		
	#read fasta file
	filename = annotatedFile
	
	commandVal['Cdr3-Fr4Identification'] = json.dumps({'method':annotation_settings['cdr3_search']['method'],'param':cdr3_search_name})	
	commandVal['Min_Alignment_Idendity_Cutoff'] = str(min_per_algn_cutoff)
	commandVal['Min_Alignment_Length_Cutoff'] = str(min_algn_len_cutoff)
	if identifyIsotype:
		commandVal['Isotyping'] = {'Barcodes':isotype_barcodes,'mismatch_cutoff':num_mismatch,'penalize_truncations':p_t,'minimum_length_cutoff':minimum_iso_alignment_length}
				
	
	if remove_insertions == 0:
		commandVal['Fix Mutations'] = 'Never' 
	elif remove_insertions == 1:
		commandVal['Fix Mutations'] = 'Always'
	elif remove_insertions == 2:
		commandVal['Fix Mutations'] = 'WhenStopCodon'
		
	
	commandString = json.dumps(commandVal)
	
	if 'isotype' in annotation_settings:
		isostring = '\t\tPerforming isotyping on sequences using provided barcodes: '+'\n'#+json.dumps(isotype_barcodes)+'\n'
	else:
		isostring= ''
	
	parse_igfft_notification = '''
	The resulting annotation output will be parsed using the following settings:
	\tMinimum alignment length cutoff: {0},
	\tMinimum percent identity: {1},
	\tFix detected insertions: {2},
	\tMethod to identify CDR3 and start of fr4: {3}\n{4}
	\n\n
	'''.format(str(min_algn_len_cutoff), str(min_per_algn_cutoff), commandVal['Fix Mutations'], json.dumps(annotation_settings['cdr3_search']['method']),isostring)
	print(parse_igfft_notification)
		
	
	useDebug = False		
	#open the igblast alignment file that was just made
	#f1 = open(annotatedFile)
	
	###PARAMETERS SET####
	vdj_hits = 0
	vj_hits = 0
	unknown_hits = 0
	
	f1 = readwrite.immunogrepFile(filelocation=annotatedFile,filetype='TAB')
	
	count = 0
	
	startPer = 0
	
	chainDic = {
		'IGH':["heavy","VDJ","IGH"],
		'IGK':["light","VJ","IGK"],
		'IGL':["light","VJ","IGL"],
		'TRA':["alpha","VJ","TRA"],
		'TRB':["beta","VDJ","TRB"],
		'VH':["heavy","VDJ","IGH"],
		'VK':["light","VJ","IGK"],
		'VL':["light","VJ","IGL"],
		'TA':["alpha","VJ","TRA"],
		'TB':["beta","VDJ","TRB"]
	}
	
	chainTypes = {
		'heavy': 'VDJ',
		'light': 'VJ',
		'alpha': 'VJ',
		'beta': 'VDJ'
	}
	
	#readfastafile =
	#readwrite.immunogrepFile(filelocation=filename,filetype=input_file_type)
	
	debugMe = False
	useDebug = False
		
	igblast_read_line = True
	
	total_parsing_errors = 0

	loop_status_gen = useful.LoopStatusGen(numSeqs,10) 
	
	tab_header_var = TABFileHeader()
	
	overlap_len = 10 
	
	#create header row
	if write_format=='TAB':		
		output_files['annotated'].write('\t'.join([field for field in tab_header_var])+'\n')
	
	t0 = 0;
	t1 = 0;
	t2 = 0;
	t3 = 0; 
	t4 = 0; 
	
	#while count<numSeqs: #read through every sequence in the file
	for query_results in f1.read():
		if count==numSeqs:						
			continue
		
		count += 1
		
		#allow the user to monitor what percent of the sequences have been processed
		loop_status_gen.next()
				
		if not(query_results):
			continue	
		
		[seqHeader,additionalInfo] = readwrite.GetAdditionalInfo(query_results['Header'])				
		query_results['Document_Header'] = query_results['Header']
		query_results['Header'] = seqHeader
												
		#additionalInfo = GrabAdditionalHeaderInfo(seqHeader)
		#query_results['Header'] = additionalInfo['Header']
	   	
		additionalInfo.pop('Header',None)
		additionalInfo.pop('document_header',None)		
		query_results = defaultdict(str,dict(query_results.items() + additionalInfo.items()))
		
		query_results['Quality_Score'] = query_results['Sequence quality']
		
		error_dic =defaultdict(str, {'Sequence':query_results['Sequence'],
					 'Header':query_results['Header'],
					 'Document_Header':query_results['Document_Header']					 
					})
					
		if idIdentifier in query_results:
			error_dic[idIdentifier] = query_results[idIdentifier]
					
		if query_results['Sequence'] != "":
			seq = query_results['Sequence']			
		else:			 
			error_dic["Notes"] = "Sequence not found"
			error_dic["Errors"] = "Sequence not found"
			error_dic["Percent_Identity"] = None
			error_dic["Alignment_Length"] = None					
			
			#Write_Seq_JSON(error_dic,tab_header_var,output_files['ERROR'])			
			
			if write_format == "TAB":								
				Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
			else:																									   				
				Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])									
			continue						

		try:				
			tt = time.time()
			if bad_dna_char_pattern.search(seq): #CHECK TO SEE IF SEQUENCE HAS WEIRD CHARACTERS IN IT
				notes = "Sequence contains unknown characters"
				print "Seq # " + str(count) + " contains unusual characters and therefore we are ignoring this sequence: " + seqHeader				
				error_dic["Percent_Identity"] = None
				error_dic["Alignment_Length"] = None
				error_dic['Notes'] = notes
				error_dic['Errors'] = notes				
				
				#Write_Seq_JSON(error_dic,chain_ind_fields,chain_dep_fields,'',output_files['ERROR'])
				#Write_Seq_JSON(error_dic,tab_header_var,output_files['ERROR'])
				if write_format == "TAB":								
					Write_Seq_TAB(error_dic,tab_header_var,output_files['annotated'])
				else:																									   				
					Write_Seq_JSON(error_dic,tab_header_var,output_files['annotated'])			
				continue				
			t0+=time.time()-tt
			
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
				print 'This chain type was not recognized in the parsing script. Analysis information for this sequence will be placed in the file "{0}". Consider updating the variable chainTypes in the funcion "Parse_Alignment_File"'.format(outfile_unknown_rec)				
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
				tt = time.time()
				[cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,cdr3start,fr4start,cdr3notes] = CDR3_SEARCH_FUNCTION(algn_seq,cdr3start,end_of_ab,j_start,fr3present,jpresent,cdr3_search_parameters,locus,var_type,jGermlineInfo,max_l_motif_len,max_r_motif_len)																	
				t1+=time.time()-tt
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
			
			tt =time.time()
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
	
			t2+=time.time()-tt
	
	
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
								
			
			tt= time.time() 
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
			
			t3+=time.time()-tt
			
			query_results["Notes"] = notes + query_results['Notes']
							
			if write_format == "TAB":								
				Write_Seq_TAB(query_results,tab_header_var,output_files['annotated'])
			else:																									   				
				Write_Seq_JSON(query_results,tab_header_var,output_files['annotated'])
				
		except Exception as e:
			total_parsing_errors+=1					
			exc_type, exc_obj, exc_tb = sys.exc_info()		
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		
			print "An error occurred when analyzing the output for sequence # "+str(count)+". This error has been reported in the error log file. Continuing with analysis"
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
	
	
def TranslateSeq(ntstring,frame):
	ntstring = ntstring[frame:]
	truncatedLength = 3*(len(ntstring)/3)
	ntstring = ntstring[:truncatedLength]
	return str(Seq(ntstring,generic_dna).translate())


	
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
				'VRegion..NT',
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
	for i,field in enumerate(chain_independent_results_order):
		igDict[field] = i	
		
	
	return igDict#,chain_independent_results_order,chain_dependent_results_order]



#We will need to update the database with the results from IgFFT.  In order
#to update teh database, we need a translator, so that we know what fields go where in the database
def DatabaseTranslator(input_dictionary = {}):
	key = translation_var
	
	translator = {					
			"ANALYSIS_NAME": "GEORGIOU_INHOUSE", # NAME OF THE ANALYSIS 
			"RECOMBINATION_FIELD":{ #THIS TELLS THE PROGRAM HOW TO DETERMINE WHETHER AN ANALYSIS/QUERY RESULT (from file) IS VDJ OR VJ
					"FIELD_NAME": "Recombination_Type", #name of the field in the file that will give information regarding the recombination type (VDJ OR VJ)
					"EXPLICIT":True, #IF EXPLICIT, then the RECOMBINATION_TYPE is equal to the EXACT value in this field (i.e. it will either list VDJ or VJ), IF false, then VDJ and VJ are defined by values in list below
					"INEXPLICIT_DEFINITIONS":{
						"VDJ":[],#if expliity is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VDJ type (i.e. if Locus is used to determine recombination type then it woudl be VDJ:['IGH']
						"VJ":[],#if expliity is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VJ type (i.e. if Locus is used to determine recombination type then it woudl be VJ:['IGK','IGL']
	 				}
			},
			"FIELDS_UPDATE":{ 
				#key = field name  in database 
				#value = field name in file
				#this will map all of the fields in the file to the proper location in the database. I.E. If I list VGENES as the column name/field name, then i want to map VREGION.VGENES:VGENES (because VREGION.VGENES is the name in the database)			
				
				idIdentifier:idIdentifier,
				"SEQUENCE":"Sequence",
				"COMMAND":"Command",
				"SEQUENCE_HEADER":"Header",	
				"QUALITY_SCORE":"Quality_Score",
				"NOTES":"Notes",
				"PREDICTED_AB_SEQ.NT":"Full_Length_Sequence.NT",
				"PREDICTED_AB_SEQ.AA":"Full_Length_Sequence.AA",
				"STRAND":"Direction",
				"PREDICTED_CHAIN_TYPE":"Chain",
				"PRODUCTIVE":"Productive",
				"LOCUS_NAME":"Locus",
				"FULL_LENGTH":"Full_Length",
				"STOP_CODONS":"Stop_Codon",
				"VREGION.SHM.NT":"VRegion.SHM.NT",
				'VREGION.VGENE_QUERY_START':'VGENE: Query_Start',
				'VREGION.VGENE_QUERY_END':'VGENE: Query_End',
				"VREGION.FR1.NT":"FR1_Sequence.NT",
				"VREGION.FR1.AA":"FR1_Sequence.AA",
				"VREGION.CDR1.NT":"CDR1_Sequence.NT",
				"VREGION.CDR1.AA":"CDR1_Sequence.AA",			
				"VREGION.FR2.NT":"FR2_Sequence.NT",
				"VREGION.FR2.AA":"FR2_Sequence.AA",
				"VREGION.CDR2.NT":"CDR2_Sequence.NT",
				"VREGION.CDR2.AA":"CDR2_Sequence.AA",
				"VREGION.FR3.NT":"FR3_Sequence.NT",
				"VREGION.FR3.AA":"FR3_Sequence.AA",
				"VREGION.VGENES":"Top_V-Gene_Hits",
				"VREGION.VGENE_SCORES":"V-Gene_Alignment_Scores",
				"CDR3.NT":"CDR3_Sequence.NT",
				"CDR3.AA":"CDR3_Sequence.AA",
				"DREGION.DGENES":"Top_D-Gene_Hits",
				"DREGION.DGENE_SCORE":"D-Gene_Alignment_Scores",
				"JREGION.FR4.NT":"FR4_Sequence.NT",		
				"JREGION.FR4.AA":"FR4_Sequence.AA",
				"JREGION.JGENES":"Top_J-Gene_Hits",
				"JREGION.JGENE_SCORES":"J-Gene_Alignment_Scores",
				'JREGION.JGENE_QUERY_START':'JGENE: Query_Start',
				'JREGION.JGENE_QUERY_END':'JGENE: Query_End',																
				"ISOTYPE.GENE":"Isotype",
				"ISOTYPE.MISMATCHES":"Isotype mismatches",
				"ISOTYPE.PER_ID":"Isotype percent similarity"		
			}
	}
	
	input_dictionary[key] = translator		
	
	return input_dictionary

def Write_Seq_TAB(seqDic,output_fields,foutfile):
	foutfile.write('\t'.join([str(seqDic[field]) for field in output_fields])+'\n')	
	
#def Write_Seq_JSON(seqDic,chain_independent_fields,chain_dep_fields,recomb_type,foutfile):
def Write_Seq_JSON(seqDic,output_fields,foutfile):
	
	#convert dictionary to json-string
	#dict_to_write = {recomb_type+'.'+field:seqDic[field] for field in chain_dep_fields if field in seqDic}
	#ind_field = {field:seqDic[field] for field in chain_independent_fields if field in seqDic}
	dict_to_write = {field:seqDic[field] for field in output_fields if field in seqDic}
	#dict_to_write = dict(dict_to_write.items()+ind_field.items())	
	
	st = json.dumps(dict_to_write)		
	foutfile.write(st + '\n')
###############################################################################

def GrabAdditionalHeaderInfo(header):
	tmp = header.split(fasta_file_delimiter)
	
	if len(tmp)>1:
		additional_info = json.loads(tmp[1])
	else:
		additional_info = {}
	additional_info.pop('document_header',None) #pop out this from dictionary if it was carried along somehow
	additional_info['Document_Header'] = header
	additional_info['Header'] = tmp[0]
	return additional_info



###############################################################################


#THIS FUNCTION HAS BEEN DEPRECATED#

#HEADERLINE VARIABLE#
#headerLine is a dictionary that defines the headerrow for the text file
#the key of the dictionary refers to the column name used for the text file
#and the value of each key is the "column number" for that column name
#if you do not have a headerLine defined, then pass in the variable as {}

#INCLUDE DECORATORS VARIABLE#
#we put comments into our analysis files using special text identifies that
#refer to "comments"
#if INCLUDE DECORATORS is set to true, then when converting the text file, it
#will also copy these "comment" lines in the datafile.  If not, then it will
#skip them
def JSON_to_TXT(input_filename,output_filename,includeDecorators,headerLine): #converts a JSON file into a text tab delimited file
	
	try:
		sorted_fields = collections.OrderedDict(sorted(headerLine.items() ,key=lambda t: t[1])).keys()	
		temp_file_path = str(datetime.now()).replace(' ','_').replace(':','').replace('.','').replace('-','_')+'.temp'	
	
		if input_filename == output_filename:
			overwriteFile = True	
			newfilename = "scratch/temp_file_name_" +temp_file_path
			f = open(newfilename,'w')							
		else:
			overwriteFile = False
			f = open(output_filename,'w')
					
		i = open(input_filename,'r')	
		
		decoratorFound = True
		
		dictHeader = {}
		dictList = []
		currentRow = 1
		headerRow = []
		
		error = False
		commentLine = descriptor_symbol

		lenComment = len(commentLine)
		
		#print dictHeader
		i = i.close()
		i = open(input_filename)
	
		decoratorFound = True
		
		while(decoratorFound):
			line_one = i.readline().strip()
			if line_one[:lenComment] != commentLine:
				decoratorFound = False		
			else:
				if includeDecorators:
					f.write(line_one + '\n')
		
		
		header_row = ''
		for header in sorted_fields:
			header_row+=header+'\t'
		header_row = header_row[:-1]+'\n'
		f.write(header_row)
										
		line = line_one.strip()
		if line and line != "" and line[:lenComment] != commentLine: #check to see if this line in the text file is a "comment line", if not, then																	   #write tofile
			myDict = json.loads(line)			
			row = ''
			for keys in sorted_fields:							
				if keys in myDict and myDict[keys] is not None:					
					if type(myDict[keys]) is list:
						row+=','.join([str(d) for d in myDict[keys]])+'\t'
					elif type(myDict[keys]) is dict:
						row+=json.dumps(myDict[keys])+'\t'
					else:											
						row+=str(myDict[keys])+'\t'
				else:
					row+='\t'
			row = row[:-1]+'\n'										
			f.write(row)
		else:
			if includeDecorators:
				f.write(line_one + '\n')			
		for line in i:				
			line = line.strip()
			if line[:lenComment] != commentLine:
				myDict = json.loads(line)
				row = ''
				for keys in sorted_fields:			
					if keys in myDict and myDict[keys] is not None:						
						if type(myDict[keys]) is list:
							row+=','.join([str(d) for d in myDict[keys]])+'\t'
						elif type(myDict[keys]) is dict:
							row+=json.dumps(myDict[keys])+'\t'
						else:											
							row+=str(myDict[keys])+'\t'
					else:
						row+='\t'
				row = row[:-1]+'\n'										
				f.write(row)
			else:
				if includeDecorators:
					f.write(line + '\n')
		
		f.close()
		i.close()
		
		if overwriteFile:
			os.system("mv '{0}' '{1}'".format(newfilename,input_filename))# + "scratch/temp_file_name_" + a + " " + input_filename)	
			
		return error
		
	except Exception as e:
		error = True
		
		exc_type, exc_obj, exc_tb = sys.exc_info()
		
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		
		
		print("there was an error when writing to file: "+str(e)+".  Please Restart")
		print "Line Number: " + str(exc_tb.tb_lineno) + " type: " + str(exc_type) + " fname: " + str(fname)
		sys.exit("ERROR FOUND SEE ABOVE")
		return error



###########FUNCTIONS FOR TEXT FILES###########################################
#DEPRECATED FUNCTION#
def WriteTABFileHeader(headerDic, outfile, filetype):
	filename = open(outfile,'w')
	
	#first write the "translator" json line so that we can update the database
	#with igblast results
	translator = DatabaseTranslator()
	translator_string = json.dumps(translator)
	translator_comment = textFileCommentTranslator
	filename.write(translator_comment + translator_string + '\n')
	
	rowData = [None] * len(headerDic)
	
	for key in headerDic:
		rowData[headerDic[key] - 1] = key
	
	for i in range(len(rowData) - 1):
		filename.write(rowData[i] + '\t')
	
	filename.write(rowData[len(rowData) - 1] + '\n')
