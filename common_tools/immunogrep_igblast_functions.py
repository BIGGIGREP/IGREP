#VERSION 350-360 SOMEWHERE IN THIS VERSION WE CHANGED HOW RECOMBINATION IS DISPLAYED 

##################################### PARSE IGBLAST RESULT FILE ####################################################################################
#THIS SCRIPT WILL PARSE THE ALIGNMENT FILE CREATED BY RUNNING IGBLAST VERSION 2
#The code is hardcoded and assumes that the output of the file is identical to the following template:
#this is sample output from the igblast program using the settings defined above


	
##############SAMPLE IGBLAST OUTPUT TEMPLATE using igblast 2.0, outfmt -7 *************************************************************
## IGBLASTN 2.2.28+
## Query: F1G46TU02JQT2D length=176 xy=3879_0291 region=2 run=R_2009_08_27_18_48_19_
## Database: scratch/igblast_db_files/database/mouse_gl_V scratch/igblast_db_files/database/mouse_gl_D scratch/igblast_db_files/database/mouse_gl_J
## Domain classification requested: imgt
	
## Note that your query represents the minus strand of a V gene and has been converted to the plus strand. The sequence positions refer to the converted sequence. 
	
## V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.
#ai4	JK5	VK	Yes	In-frame	No	-
	
## V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
#TCCCC	G	CTCAC	
	
## Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)
#FR3-IMGT	1	109	109	106	2	1	97.2
#CDR3-IMGT (germline)	110	129	20	20	0	0	100
#Total	N/A	N/A	129	126	2	1	97.7
	
## Hit table (the first field indicates the chain type of the hit)
## Fields: subject id, query length, subject length, q. start, q. end, s. start, s. end, query seq, subject seq, evalue, alignment length, identical, mismatches, gaps, query frame, sbjct frame, subject strand
## 6 hits found
#V	ai4	176	291	1	129	160	287	AAACTGGCTTCTGGAGTCCCAGCTTCGCTTCAGTGGCAGTGGGTCTGGGACCTCTTGCTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCACCAGTATCATCGTTCCCC	AACCTGGCTTCTGGAGTCCCAGCT-CGCTTCAGTGGCAGTGGGTCTGGGACCTCTTACTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCACCAGTATCATCGTTCCCC	1e-49	129	126	2	1	1	1	plus
#V	aa4	176	285	1	129	154	281	AAACTGGCTTCTGGAGTCCCAGCTTCGCTTCAGTGGCAGTGGGTCTGGGACCTCTTGCTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCACCAGTATCATCGTTCCCC	AACCTGGCTTCTGGAGTCCCTGCT-CGCTTCAGTGGCAGTGGGTCTGGGACCTCTTACTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCAGCAGTATCATAGTTACCC	6e-46	129	122	6	1	1	1	plus
#V	ar4	176	285	1	129	154	281	AAACTGGCTTCTGGAGTCCCAGCTTCGCTTCAGTGGCAGTGGGTCTGGGACCTCTTGCTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCACCAGTATCATCGTTCCCC	AACCTGGCTTCTGGAGTCCCTGCT-CGCTTCAGTGGCAGTGGATCTGGGACCTCTTATTCTCTCACAATCAGCAGCATGGAGGCTGAAGATGCTGCCACTTATTACTGCCAGCAATATCATAGTTACCC	4e-43	129	119	9	1	1	1	plus
#J	JK5	176	39	131	153	1	23	CTCACGTTCGGTGCTGGGACCAA	CTCACGTTCGGTGCTGGGACCAA	9e-10	23	23	0	0	1	1	plus
#J	JK1	176	39	153	170	22	39	AAGCTGGAAATCAAACGT	AAGCTGGAAATCAAACGT	9e-07	18	18	0	0	1	1	plus
#J	JK2	176	39	153	170	22	39	AAGCTGGAAATCAAACGT	AAGCTGGAAATAAAACGT	2e-04	18	17	1	0	1	1	plus
	
## IGBLASTN 2.2.28+
## Query: F1G46TU02I7MEL length=224 xy=3660_0923 region=2 run=R_2009_08_27_18_48_19_
## Database: scratch/igblast_db_files/database/mouse_gl_V scratch/igblast_db_files/database/mouse_gl_D scratch/igblast_db_files/database/mouse_gl_J
## Domain classification requested: imgt
	
## V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.
#IgK9-128	N/A	VK	Yes	N/A	No	+
	
## V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself
#GGGCA	N/A	N/A	
	
## Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)
#FR1-IMGT	15	92	78	75	3	0	96.2
#CDR1-IMGT	93	110	18	18	0	0	100
#FR2-IMGT	111	160	51	50	0	1	98
#CDR2-IMGT	161	169	9	9	0	0	100
#FR3-IMGT	170	218	50	49	0	1	98
#Total	N/A	N/A	206	201	3	2	97.6
	
## Hit table (the first field indicates the chain type of the hit)
## Fields: subject id, query length, subject length, q. start, q. end, s. start, s. end, query seq, subject seq, evalue, alignment length, identical, mismatches, gaps, query frame, sbjct frame, subject strand
## 3 hits found
#V	IgK9-128	224	287	15	218	1	206	GACATTGTGATGACCCAGTCTCCATCCTCCATGTATGCATCGCTGGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAAAAGCTATTTAAGCTGGTACCAGCAG-AACCATGGAAATCTCCTAAGACCCTGATCTATTATGCAACAAGCTTGGCAGATGGGGTCCCATCAAGATTCAGT-GCAGTGGATCTGGGCA	GACATCAAGATGACCCAGTCTCCATCCTCCATGTATGCATCGCTGGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAAAAGCTATTTAAGCTGGTACCAGCAGAAACCATGGAAATCTCCTAAGACCCTGATCTATTATGCAACAAGCTTGGCAGATGGGGTCCCATCAAGATTCAGTGGCAGTGGATCTGGGCA	6e-82	206	201	3	2	1	1	plus
#V	br9	224	286	15	218	1	207	GACATTGTGATGACCCAG-TCTCCATCCTCCATGTATGCATCGCTGGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAAAAGCTATTTAAGCTGGTACCAGCAG-AACCATGGAAATCTCCTAAGACCCTGATCTATTATGCAACAAGCTTGGCAGATGGGGTCCCATCAAGATTCAGT-GCAGTGGATCTGGGCA	GACATCAAGATGACCCAGATCTCCATCCTCCATGTATGCATCGCTGGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAAAAGCTATTTAAGCTGGTACCAGCAGAAACCATGGAAATCTCCTAAGACCCTGATCTATTATGCAACAAGCTTGGCAGATGGGGTCCCATCAAGATTCAGTGGCAGTGGATCTGGGCA	1e-79	207	201	3	3	1	1	plus
#V	ba9	224	285	15	218	1	206	GACATTGTGATGACCCAGTCTCCATCCTCCATGTATGCATCGCTGGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAAAAGCTATTTAAGCTGGTACCAGCAG-AACCATGGAAATCTCCTAAGACCCTGATCTATTATGCAACAAGCTTGGCAGATGGGGTCCCATCAAGATTCAGT-GCAGTGGATCTGGGCA	GACATCAAGATGACCCAGTCTCCATCTTCCATGTATGCATCTCTAGGAGAGAGAGTCACTATCACTTGCAAGGCGAGTCAGGACATTAATAGCTATTTAAGCTGGTTCCAGCAGAAACCAGGGAAATCTCCTAAGACCCTGATCTATCGTGCAAACAGATTGGTAGATGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGCA	1e-69	206	188	16	2	1	1	plus
	
##################################################################################################

from Bio import SeqIO
import sys
import commands
import re
import appsoma_api
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import time
from time import gmtime, strftime

from collections import defaultdict
from collections import OrderedDict
import immunogrep_useful_immunogrep_functions as useful

import json
import os
import traceback
from pprint import pprint

#import immunogrep_run_igblast_command as run_igblast
import immunogrep_immunogrepfile as readwrite

from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import fasta_file_delimiter
from immunogrep_global_variables import descriptor_symbol #this is the symbol we will use to designate 'comment' lines in the file 
from immunogrep_global_variables import translation_var #this key will signify the translation/translator key 
from immunogrep_global_variables import filetype_var
from datetime import datetime

import immunogrep_appsoma_cdr3tools as CDR3tools
#from cchrysostomou_readwritefiles import *


#download bash script and run. this bash command actuall runs igblast 
import immunogrep_make_igblast_database as newigblastdb
import immunogrep_query_germline_functions as query_germlines #use this for querying germlines 

appsoma_api.resource_pull("https://www.appsoma.com/programs/get/immunogrep_igblast_run.bash","run_igblast_bash.bash",override_cache=True)

default_germline_folder = 'scratch/igblast_db_files/database'

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
def download_germlines_from_db(query_settings,db_type='nucl'):
	query = {}
	output_directory = default_germline_folder+'/databasedownloads/'
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
	elif 'SPECIES' in query_settings:
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
		productive_gene_filters = {'$in':query_settings['PRODUCTIVITY']}
	else:
		productive_gene_filters = {'$nin':[]}
					
	
	#create a query for getting germlines from database 
	db_class_var = query_germlines.GermlineDB()						
	#this is just for documenting purproses/keeping track of what teh query request was. storing settings of query basically 
	unique_settings = db_class_var.QueryDistinctValsByID(id_list,extra_filters=query,distinct_fields=['SPECIES','GENETYPE','MOLECULAR_COMPONENT','SOURCE','VERSION','LOCUS'])					
	unique_settings['DB-FILE-SOURCE'] = unique_settings.pop('SOURCE')
	selected_genes = unique_settings['GENETYPE']	
	if 'PRODUCTIVITY' in query_settings:
		unique_settings['PRODUCTIVITY'] = query_settings['PRODUCTIVITY']	

	
	#figure out a name for the database file that makes it easy to recognize when referred to later on 
	name_data = defaultdict(list)
	for s in unique_settings['SPECIES']:
		for l in unique_settings['LOCUS']:
			name_data[s].append(l)	
	file_prefix = 'DatabaseDownload_'+'_'.join({s+'_'.join(v) for s,v in name_data.iteritems()})+db_type #should create a file name of datbasedownload_species_all loci_speices_all loci.txt
	file_prefix=file_prefix.replace(' ','')			
	
	#actually run the query and download germlines as a FASTAfile for igblast format
	db_class_var.QueryGenesByID(id_list,extra_filters=query,gene_functionality_filter = productive_gene_filters).PrintIgBlastDBFormat(output_directory,filename=file_prefix) #download the necessary germline files from the database
		
	#ensure sequences downloaded correctly 	
	selected_genes = unique_settings['GENETYPE']	
	germlines = {}	
	for subtype in ['V','D','J']:
		if subtype in selected_genes:			
			filepath = output_directory+file_prefix+'_'+subtype+'.txt'#the function QueryGenesById will create germline database files with these filesnames				
			if os.path.isfile(filepath):													
				output_db_germ_name = useful.removeFileExtension(os.path.basename(filepath))																				
				#make an inglbast database file 
				newigblastdb.MakeIgBlastDB(fasta_filename=filepath,output_filename=default_germline_folder+'/'+output_db_germ_name,db_type=db_type) 						
				germlines[subtype] = output_db_germ_name		
	
	return [germlines,unique_settings]





#this python script will run an igblast command. 
#you pass in:
	#testSEt -> location ofthe input file
	#outfile -> desired location and name of the output file
	#variable_parameters -> optional parameters for running igblast
	#skipRunning -> dont run igblast, just generate a string command
	#input_Filetype -> the file type of the input file (FASTA, TAB, CSV, ETC)
	#header_field -> the field in the file that corresponds to the sequence header
	#sequence_field -> the field in the file that corresponds to the seqeuence 

#germline parameter is as follows:
	#a three (or four) index tuple  => see 'HACK' below to see when to use a 'four' element tuple 
		#first element => string defining method for describing where germline came from: default, custom, db
			#default => use default files currently in folder 
			#custom => provide custom text files 
			#db => download germlines from database 
			#make => pass in a FASTA file and we make the blast database for you 
		#second element => defines the parameters for the germline file locations
			#if element 1 => default:
				#provide a two element tuple defining the species and the molecular componenet
					#allowed species: Mus musuculus, homo sapiens 
					#allowed molecular components: IG, TR 
			#if element 2=> provide a dictionary whose keys are {V,D, AND/OR J} and the values are the filepath for each germline file 			
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
		#third element => whether it is a 'nucl' or 'aa' database 
		#fourth elmeent => whether we have already downloaded germliens from db: (see HACK below)
			#in this case, the second element must be {'FILES':{'v','d','j'}, and 'PARAMS'}
def RunIgBlast(testSet,germlines={},outfile=None,variable_parameters={}, skipRunning=False,input_filetype='FASTA', sequence_type = 'NT',header_field='header',sequence_field='sequence'):
	if variable_parameters==None:
		variable_parameters = {}
	if not(outfile):
		outfile=testSet+'.igblast_results'
		
	igblast_internal_files = 'scratch/igblast_db_files'
	
	
	if not germlines:
		#we assume that the user manually provided v and j parameters to variable parameters. therefore, this is a CUSTOM source 
		
		custom_details = {}
		v = variable_parameters.pop('germline_db_V',None)
		j = variable_parameters.pop('germline_db_J',None)
		d = variable_parameters.pop('germline_db_D',None)
		if v:			
			custom_details['germline_db_V']=v
		if j:			
			custom_details['germline_db_J']=j
		if d:
			custom_details['germline_db_D']=d
		if custom_details == {}:
			if variable_parameters['organism'] == 'human':
				custom_details['germline_db_V'] = 'human_gl_V'
				custom_details['germline_db_D'] = 'human_gl_D'
				custom_details['germline_db_J'] = 'human_gl_J'
			elif variable_parameters['organism'] == 'mouse':
				custom_details['germline_db_V'] = 'mouse_gl_V'
				custom_details['germline_db_D'] = 'mouse_gl_D'
				custom_details['germline_db_J'] = 'mouse_gl_J'
			else:
				#ther user has not provided any information at all regarding a germline 
				raise Exception("You must provide at least one germline database. If you do not have a germline database then define the variable germlines to use default germline database in program or download germline from GG lab database. See documentation in function")
		else:
			germlines = ('custom',custom_details)
	


	#the possible ways of getting germlien files is from 'default' (use default files in folder) , 'custom' (provide your own files), or 'db' (database)
	if not isinstance(germlines,tuple):
		raise Exception('Germline definition must be a tuple of two elements. Element one is a string of either "custom","db", or "default". Element two defines details concerning the specific germlines desired')
	
	subtypes = ['v','d','j'] 
	
	if germlines[0]=='custom': #we assume you have already created a blast database format 
		germline_settings = {'GERMLINE-SOURCE':'CUSTOM-FILE',
							 'PARAMS':{}
							 }							
		#the user is providing proper germline files for each 
		for keys,paths in germlines[1].iteritems():
			keys = keys.lower()
			if keys in subtypes and os.path.isfile(paths):				
				variable_parameters['germline_db_'+keys.upper()] = paths				
				
	elif germlines[0]=='make': #we assume that you want to make a blast database format for search. You pass in keys referring to database file names, and we will make the datatype
		#it shoudl be three elements ('make',filesdictionary,method/type of database (nucl/aa))
		if len(germlines)<=2: 
			db_type = 'nucl'
		else:
			db_type = germlines[2]
			
		germline_settings = {'GERMLINE-SOURCE':'CUSTOM-FILE',
							 'PARAMS':{}
							 }							
		
		
		#the user is providing proper FASTA germline files for each 
		for keys,paths in germlines[1].iteritems():
			keys = keys.lower()
			if keys in subtypes and os.path.isfile(paths):	
				input_db_germ_name = paths #PATHS MUST BE A FASTA FILE!!!
				output_db_germ_name = useful.removeFileExtension(os.path.basename(input_db_germ_name))				
				#make an inglbast database file 
				newigblastdb.MakeIgBlastDB(fasta_filename=input_db_germ_name,output_filename=default_germline_folder+'/'+output_db_germ_name,db_type=db_type) 						
				variable_parameters['germline_db_'+keys.upper()] = output_db_germ_name																											
		
	
	elif germlines[0] == 'default':
		germline_settings = {'GERMLINE-SOURCE':'DEFAULT',
							 'PARAMS':{}
						 }			
		#if organism nto defined, then assume it is human. 
		if 'organism' not in variable_parameters or ('organism' in variable_parameters and variable_parameters['organism'].lower() == 'human'):
			variable_parameters['germline_db_V']='human_gl_V', # name of v germline database files (currently assumes database files are located within scrath/igblast_db_files/database/{db_name*}
			variable_parameters['germline_db_D']='human_gl_D', # name of d germline database files
			variable_parameters['germline_db_J']='human_gl_J', # name of j germline database files			
			germline_settings['PARAMS']['organism']='human'
		elif 'organism' in variable_parameters and variable_parameters['organism'].lower() == 'mouse':
			variable_parameters['germline_db_V']='mouse_gl_V', # name of v germline database files (currently assumes database files are located within scrath/igblast_db_files/database/{db_name*}
			variable_parameters['germline_db_D']='mouse_gl_D', # name of d germline database files
			variable_parameters['germline_db_J']='mouse_gl_J', # name of j germline database files			
			germline_settings['PARAMS']['organism']='mouse'
		else:
			raise Exception('We only support default germlines for human or mouse. please provide a custom germline database for other speices')
												
	elif germlines[0] == 'db':
		if len(germlines)<=2: 
			db_type = 'nucl'
		else:
			db_type = germlines[2]
			
		#the user wants to download germlines from database		
		germline_settings = {'GERMLINE-SOURCE':'DATABASE-DOWNLOAD'}		
		
		#HACK:
			#this ifstatement is a hack for the following reason: We will want to be able to download germlines from the database. 
			#But if we also want to do multithreading, then all of the downloaded files will keep overwriting one another 
			#therefore, for situations when we do multithreading, we add an extra term to the tuple defining the database
				#this fourth element basically is TRUE/FALSE. IF true it means we did want to download germlines from the database, but the germlines were already downloaded before hand, so no need to 
				#run the funciton. Instead the value of the database will correspond to the blast format database 
		if(len(germlines)>3) and germlines[3]==True: #so we have already downloaded the data from the database 
			#IF THIS IS THE CASE THEN THE PROVIDED FORMAT MUST BE 'FILES':{'v','d','j'}, 'PARAMS':
			germline_files = germlines[1]['FILES']
			germline_settings['PARAMS']  = germlines[1]['PARAMS']
		else:
			#for all other situations that ARE NOT USING OUR MULTITHREADING FUCTNION, THE FORMAT IS THE STANDARD FORMAT DISCUSSED ABOVE 
			[germline_files, germline_settings['PARAMS']] = download_germlines_from_db(germlines[1],db_type)
		#end of hack 
		
		for s in subtypes:
			s = s.upper()
			if s in germline_files:
				#get teh file listed inside 				
				#annoying, but igfft requires lowercalse subtypes
				variable_parameters['germline_db_'+s] = germline_files[s]	 										
	else:		
		raise Exception('Germline definition must be a tuple of two elements. Element one is a string of either "make","custom","db", or "default". Element two defines details concerning the specific germlines desired')						
		
	

	#make sure igblast parameters are good in variable_parameters dictionary  
	default_igblast_parameters = {
		'organism':'human',
		'ig_seqtype':'Ig',
		'germline_db_V':'human_gl_V', # name of v germline database files (currently assumes database files are located within scrath/igblast_db_files/database/{db_name*}
		'germline_db_D':'human_gl_D', # name of d germline database files
		'germline_db_J':'human_gl_J', # name of j germline database files															
		'num_alignments_V':5, #desired number of v genes
		'num_alignments_D':5, #desired number of d genes
		'num_alignments_J':5,  #desired number of j genes
		'domain_system':'imgt'
	}
		 	
	
	
	for var in default_igblast_parameters: #if an igblast variable is not define,d then use the default variable
		if not (var in variable_parameters):
			variable_parameters[var] = default_igblast_parameters[var]
	
	var_keys = list(variable_parameters.keys())
	for var in var_keys: #remove any keys that are not pertinent to running igblast -> leave out for now
		if not (var in default_igblast_parameters):
			variable_parameters.pop(var,None)
	
	variable_parameters['organism'] = variable_parameters['organism'].lower()
	
	if variable_parameters['organism']=="human":
		variable_parameters['auxiliary_data'] = "{0}/optional_file/human_gl.aux".format(igblast_internal_files)
	elif variable_parameters['organism']=="mouse":
		variable_parameters['auxiliary_data'] = "{0}/optional_file/mouse_gl.aux".format(igblast_internal_files)
	
	variable_parameters["germline_db_V"] = igblast_internal_files+'/database/'+variable_parameters["germline_db_V"]
	variable_parameters["germline_db_D"] = igblast_internal_files+'/database/'+variable_parameters["germline_db_D"]
	variable_parameters["germline_db_J"] = igblast_internal_files+'/database/'+variable_parameters["germline_db_J"]
	##variables defined##
	
	
	####first, read in the input file and conver it to a fasta file that can be read by IgBlast
	numSeqs = 0;	
	inputIFFile = readwrite.immunogrepFile(filelocation=testSet,filetype=input_filetype)
	date = str(datetime.now())
	remove_chars = [':','-',' ','.']
	for char in remove_chars:
		date = date.replace(char,'')		
	tempFASTAFile=testSet+'_{0}'.format(date)#temporary FASTA file we will make for running IgBlast 
	#MAKE SURE THIS FILE IS NOT ALREADY PRESENT  (for multithreading)
	while True:
		if not os.path.isfile(tempFASTAFile):			
			output_temp_file = open(tempFASTAFile,'w')
			break
		tempFASTAFile+='1' #add some 1's 
		
	
	
	while not(inputIFFile.IFclass.eof):
		line = inputIFFile.IFclass.read()
		if line:
			numSeqs+=1
			#db_info = '{0}{1}'.format(fasta_file_delimiter,json.dumps({idIdentifier:line[idIdentifier]})) if idIdentifier in line else ''				
			if sequence_field in line and line[sequence_field] != '':							
				output_temp_file.write('>{0}\n{1}\n'.format(line[header_field].strip(),line[sequence_field].strip()))#,db_info))
	output_temp_file.close()
	inputIFFile.IFclass.close()
	####temp fasta file made ##
			
	
	############INPUT VARIABLES FOR THIS SCRIPT#############################			
			
	outfileIgBlastAlgn = outfile
	
	##generate igblast command 
	file_source_command = '-query {0} -out {1} '.format(tempFASTAFile,outfileIgBlastAlgn)
	igblast_command=""
	appsoma_api_text = '\n\nRunning IgBlast for {0} sequences using the following settings:\n'.format(numSeqs)
	appsoma_api_text += "\tInput file: {0}\n".format(tempFASTAFile)
	appsoma_api_text += "\tOutput file will be saved as: {0}\n".format(outfileIgBlastAlgn)
	
	for var in variable_parameters:
		if (var == 'auxiliary_data' or var=='num_alignments_D' or var=="num_alignments_J" or var == 'germline_db_D' or var=='germline_db_J') and sequence_type == 'AA':
			continue
		igblast_command+='-{0} {1} '.format(var,variable_parameters[var])
		appsoma_api_text+='\tParameters {0}: {1}\n'.format(var,variable_parameters[var])	
	#####bash line / string commmand prompt for runing igblast/ made
	
	
	###############RUN IGBLAST###################################################################			
	#run IgBlast in bash using igBlastCommand line
		
	if not(skipRunning):
		appsoma_api_text+= "IgBlast Run has started at {0}\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()))
		
		print(appsoma_api_text)
		os.system(''' bash run_igblast_bash.bash "{0}" {1} '''.format(file_source_command+igblast_command,sequence_type.upper()))	
		print("IgBlast Analysis Completed at {0}.\nAnalysis saved to: {1}.\n\n\n".format(strftime("%a, %d %b %Y %X +0000", gmtime()),outfileIgBlastAlgn))
						
	
	os.system("rm '{0}'".format(tempFASTAFile))
	
	#for documentation purposes in database 
	commandVal = variable_parameters
	commandVal.pop('auxiliary_data',None)
	vgenes = commandVal.pop('germline_db_V',None)
	jgenes = commandVal.pop('germline_db_J',None)
	dgenes = commandVal.pop('germline_db_D',None)	
	commandVal['annotation'] = 'IGBLAST'
	
	commandVal['germline_db_V']=os.path.basename(vgenes) if vgenes else vgenes
	commandVal['germline_db_D']=os.path.basename(dgenes) if dgenes else dgenes
	commandVal['germline_db_J']=os.path.basename(jgenes) if jgenes else jgenes
	commandVal['germline-details']=germline_settings
		
	return [commandVal,numSeqs]

	#############################################################################################	

def IdentifyCDR3UsingRegExp(algn_seq,cdr3start,end_of_ab,cdr3_index,cdr3_search_parameters,locus,var_type,Jgermline_alignments):
	
	cdr3_fr4_seq = algn_seq[cdr3start-1:end_of_ab]
			
	#Determine whether we are searching for a motif or not 
	motif = cdr3_search_parameters[var_type] if cdr3_search_parameters and var_type in cdr3_search_parameters and cdr3_search_parameters[var_type] else None
							
	if motif:						
		[fr4_start,result_notes] = Find_CDR3_Motif(motif,Jgermline_alignments,cdr3_fr4_seq,cdr3start)
		
	else:			
		result_notes = ''
		fr4_start = cdr3start+3*int( (cdr3_index - cdr3start + 1) / 3 )					
		
	if fr4_start!=-1 and fr4_start>cdr3start:
		cdr3_nt = algn_seq[cdr3start-1:fr4_start-1]		
		cdr3_aa = str(Seq(cdr3_nt,generic_dna).translate())#.tostring()
			
		fr4_nt = algn_seq[fr4_start-1:end_of_ab]
		fr4_aa = str(Seq(fr4_nt,generic_dna).translate())#.tostring()
		
		cdr3Frame = (cdr3start-1)%3+1
		fr4Frame = (fr4_start-1)%3+1		
	else:
		if cdr3_index>cdr3start:
			cdr3_nt = algn_seq[cdr3start-1:cdr3_index] 
			cdr3_aa = str(Seq(cdr3_nt,generic_dna).translate())#.tostring()
			cdr3Frame = (cdr3start-1)%3+1
		else:
			cdr3_nt = None
			cdr3_aa = None
			cdr3Frame = None
			
		fr4_nt= None
		fr4_aa = None
		fr4Frame = None
	
	return [cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,result_notes]


#this is the function for finding the CDR3 using a motif. There are two options for finding the motif. 
#a) search through all of the possible J-GERMLINE hits for the motif of interest.  If we find any of the j-germline-gene alignments, then we will map the position back to the proper position in the sequence
#b) search the CDR3-FR4 subsequence (antibody after the FR3 more or less) for the motif 

def Find_CDR3_Motif(motif,j_gene_algn_info,cdr3_fr4_subseq,cdr3start):
			
	result_notes=''
	
	motif_found = False

	#first we loop through all possible j-germline alignments. in each germline alignment, we search for the motif int he germline
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
			
		if regcount>0:
			motif_found = True
			
			germ_fr4start = m.start()
			fr4start = j_gene_start-1
			
			i = 0
			counter_char = 0
			
			#The motif was found. Now we want to map the motif found in the germline to the proper position in the sequence 
			while i != germ_fr4start:					
				if j_gene_query[counter_char] != '-':
					fr4start+=1	
				if j_gene_germline[counter_char] != '-':
					i+=1						
				counter_char+=1			
			fr4start+=1
			break #exit_loop 
			
	if not(motif_found): #we couldnt find the motif in the germline, so we need to search the cdr3-fr4 sequence			
		mtch=re.finditer(motif,cdr3_fr4_subseq,re.IGNORECASE) #search for all instances of motif in the query sequence j-region 
	
		regcount = 0		
		for loopreg in mtch:																							
			m = loopreg #store the last occurrence of the motif 
			regcount+=1
		
		if regcount > 0: 
			motif_found = True
			fr4start = m.start() + cdr3start
			
			
	if not(motif_found):					
		fr4start=-1
		result_notes = 'RegExp not found;'
		
				
	return [fr4start,result_notes]
					

def IdentifyCDR3UsingDefaultIgBlast(algn_seq,cdr3start,end_of_ab,cdr3_index,cdr3_search_parameters=None,locus=None,var_type=None,Jgermline_alignments=None):
	fr4_start = cdr3start+3*int( (cdr3_index - cdr3start + 1) / 3 )
	result_notes = ''
	if fr4_start!=-1 and fr4_start>cdr3start:
		cdr3_nt = algn_seq[cdr3start-1:fr4_start-1]		
		cdr3_aa = str(Seq(cdr3_nt,generic_dna).translate())#.tostring()
			
		fr4_nt = algn_seq[fr4_start-1:end_of_ab]
		fr4_aa = str(Seq(fr4_nt,generic_dna).translate())#.tostring()
		
		cdr3Frame = (cdr3start-1)%3+1
		fr4Frame = (fr4_start-1)%3+1		
	else:
		if cdr3_index>cdr3start:
			cdr3_nt = algn_seq[cdr3start-1:cdr3_index] 
			cdr3_aa = str(Seq(cdr3_nt,generic_dna).translate())#.tostring()
			cdr3Frame = (cdr3start-1)%3+1
		else:
			cdr3_nt = None
			cdr3_aa = None
			cdr3Frame = None
			
		fr4_nt= None
		fr4_aa = None
		fr4Frame = None
	
	return [cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,result_notes]
	

def IdentifyCDR3UsingGGJunctionalMotif(algn_seq,cdr3start,end_of_ab,cdr3_index,cdr3_search_parameters=None,locus=None,var_type=None,Jgermline_alignments=None):
# PWM() function will check a sequence for the starting position and the ending position of the best motif
# Required inputs: PWM(seq<>text,LMotif<>dict,frame<>[0,1,2,3,4,5],start_pos<>intDefault0,motif_type<>AA/NTDefaultAA)
# Note: the start_pos is referencing on the original sequence input
	cdr3start-=1
	sub_string = algn_seq[cdr3start:end_of_ab+1]
	results = [CDR3tools.PWM(sub_string,x[2],frame=[0,1,2]) if x[0] in locus or locus == [] else [-100]*4 for x in cdr3_search_parameters] 
	
	maxMotif = -1
	maxP = 1e-8
	for i,r in enumerate(results):
		if r[3]>maxP:
			maxMotif = i
			maxP = r[3]
			beststart = r[0]
			bestend = r[1]
			bestframe = r[2]
	
	if maxMotif>-1 and beststart>0:
		trim = cdr3_search_parameters[maxMotif][4] if len(cdr3_search_parameters[maxMotif])==5 else 1
		fr4_start = beststart+cdr3start+(trim-1)-1
		cdr3_nt = algn_seq[cdr3start:fr4_start]
		fr4_nt = algn_seq[fr4_start:end_of_ab]
		cdr3Frame = cdr3start%3+1
		fr4Frame = fr4_start%3+1
		cdr3_aa = TranslateSeq(cdr3_nt,0)
		fr4_aa = TranslateSeq(fr4_nt,0)
		result_notes = ""
	else:
		fr4_start = -1
		cdr3_nt = ""
		fr4_nt = ""
		fr4_aa = ""
		cdr3_aa = ""
		cdr3Frame = cdr3start%3+1
		fr4Frame = 0
		result_notes = "FR4 Not Found;"
	
	return [cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,result_notes]	

def TranslateSeq(ntstring,frame):
	ntstring = ntstring[frame:]
	truncatedLength = 3*(len(ntstring)/3)
	ntstring = ntstring[:truncatedLength]
	return str(Seq(ntstring,generic_dna).translate())
	
#we use this function belwo in Parse_Igblast_file
def ReadIgBlastQueryBlock(filename):	
	f1 = open(filename)		
	tempLine = f1.readline() 
	if tempLine[1]+" "+tempLine[2]=="BLAST processed":
		igblast_eof = True
	else:
		igblast_eof = False	
	igblast_query_block = ""
	igblast_query_seq_header = ""		
	while (not(igblast_eof)): 		
		line_info = f1.readline()				
		if not(line_info):
			igblast_eof = True			
			break		
		if line_info.strip() == "": #skip empty lines 
			continue							
		if line_info[:9].upper()== "# IGBLAST": #start of a new query 
			igblast_read_line = False
			####NOW WE SHOULD HAVE READ IN AN ENTIRE BLOCK FOR THE CURRENT IGBLAST query
			yield [igblast_query_block,igblast_query_seq_header]				
			igblast_query_block = "" #reset query_block to empty
			igblast_query_seq_header = ""#reset seq_header to empty		
		elif line_info[:17].upper() =="# BLAST PROCESSED": #at the complete end of the igblast files							
			igblast_eof = True						
		else:
			igblast_query_block += line_info					
		if line_info[:8].upper() == "# QUERY:":
			igblast_query_seq_header = line_info[9:].strip()					
	yield [igblast_query_block,igblast_query_seq_header]					
	f1.close()


#This function will ONLY be used for printing an igblast result to a file or to html using the PreviewFile class
def PrintToScreenIgBlastQueryBlock(filename):	
	f1 = open(filename)		
	tempLine = f1.readline() 
	if tempLine[1]+" "+tempLine[2]=="BLAST processed":
		igblast_eof = True
	else:
		igblast_eof = False	
	igblast_query_block = ""
	igblast_query_seq_header = ""		
	while (not(igblast_eof)): 		
		line_info = f1.readline()				
		if not(line_info):
			igblast_eof = True			
			break		
		if line_info.strip() == "": #skip empty lines 
			continue							
		if line_info[:9].upper()== "# IGBLAST": #start of a new query 
			igblast_read_line = False
			####NOW WE SHOULD HAVE READ IN AN ENTIRE BLOCK FOR THE CURRENT IGBLAST query
			yield igblast_query_block
			igblast_query_block = "" #reset query_block to empty
			igblast_query_seq_header = ""#reset seq_header to empty		
		elif line_info[:17].upper() =="# BLAST PROCESSED": #at the complete end of the igblast files							
			igblast_eof = True						
		else:
			igblast_query_block += line_info					
		if line_info[:8].upper() == "# QUERY:":
			igblast_query_seq_header = line_info[9:].strip()					
	yield igblast_query_block
	f1.close()
	
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
				
	
				
#input_format = list of parametesr for the input file.  first value = filetype, second value = field name for header sequence, third value = field name for the sequence	
def Parse_IgBlast_File(inputFile,igBlastFile,outfile=None,commandVal = {},numSeqs=0,seq_type = 'NT',input_format=['FASTA','header','sequence','phred'],write_format='JSON',return_results_dict=False, annotation_settings={}):
	
	seq_type = seq_type.upper()
	
	dna_alphabet = "ACDEFGHIJKLMNPQRSTVWXYZ*" if seq_type == 'AA' else "ACTGUKMRYSWBVHDXN"
	bad_dna_char_pattern = re.compile('[^'+dna_alphabet+']',re.IGNORECASE)
	
	print('\n\nParsing and summarizing alignment file\n')
	
	if not(outfile):
		outfile = useful.removeFileExtension(igBlastFile)+'.igfft.annotation'
		
	#outfile = useful.removeFileExtension(outfile)+".igblast.annotation"
	outfile_errors=outfile+".igblast.error_log" #this file will be used to note any errors that occur while parsing the file
	
	if igBlastFile == outfile:
		os.system('''mv '{0}' '{0}.temp' '''.format(igBlastFile))
		igBlastFile+='.temp'
	
	
	input_file_type = input_format[0]
	header_name = input_format[1]
	seq_name = input_format[2]
	quality_name = input_format[3] if len(input_format)>3 else ''
		
		
	#decides what min cutoff to use to consider a successful sequence match 
	if 'min_cutoff' in annotation_settings:
		min_per_algn_cutoff = annotation_settings['min_cutoff']
	else:
		min_per_algn_cutoff = 0.4
	
	#minimum number of sequences that must align to germlines to be considered antibody
	if 'min_len_cutoff' in annotation_settings:
		min_algn_len_cutoff = annotation_settings['min_len_cutoff']
	else:
		min_algn_len_cutoff = 100
	
	#remove insertions in sequence
	if 'remove_insertions' in annotation_settings and seq_type != 'AA':
		remove_insertions = annotation_settings['remove_insertions']
	else:
		remove_insertions = 0 #0 -> dont remove insertions, 1 -> remove insertions always, 2-> remove insertions only if there is a stop codon
	
	if seq_type=='AA':
		force_to_v_germ_beg = False
		force_to_j_germ_end = False
	else:
		if 'align-to-beg' in annotation_settings:
			force_to_v_germ_beg = annotation_settings['align-to-beg']
		else:
			force_to_v_germ_beg = True # if the sequence alignment starts at position that is not 1 of the germline (i.e. starts at position 15), but the sequence has enough bases,  then force the sequence to align to the front of ther germline
		
		if 'align-to-end' in annotation_settings:
			force_to_j_germ_end = annotation_settings['align-to-end']
		else:
			force_to_j_germ_end = True # if the sequence alignment starts at position that is not 1 of the germline (i.e. starts at position 15), but the sequence has enough bases,  then force the sequence to align to the front of ther germline	
		
	
	#determine which CDR3 function we will use to identify cdr3
	if 'cdr3_search' in annotation_settings:
		
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
					germline_query.append({'Species':possible_combos[0],'Locus':list(set(possible_combos[1]))}) #using list(set to ensure only unique values in list 
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
			CDR3_SEARCH_FUNCTION = IdentifyCDR3UsingRegExp #point the fucntion cdr3_search_function to the method defined for regular expression
			cdr3_search_parameters = annotation_settings['cdr3_search']['param']
			cdr3_search_name = cdr3_search_parameters
		
		elif annotation_settings['cdr3_search']['method'] == 'default':
			CDR3_SEARCH_FUNCTION = IdentifyCDR3UsingDefaultIgBlast 
			cdr3_search_parameters = annotation_settings['cdr3_search']['param']				
			cdr3_search_name = 'IgBlast Default CDR3 End'
		else:
			raise Exception('Only allowed values for cdr3_search method are: gg_inhouse,regexp,or default')
	else:
		annotation_settings['cdr3_search'] = {'method':None,'param':None}
		CDR3_SEARCH_FUNCTION = IdentifyCDR3UsingDefaultIgBlast 
		cdr3_search_parameters = annotation_settings['cdr3_search']['param']				
		cdr3_search_name='None selected'
				
	
	#CHECK IF IGBLAST FILE WAS ACTUALLY CREATED
	
	if not(os.path.isfile(igBlastFile)):   			
   		raise Exception("ERROR: IGBLAST FILE WAS NOT CREATED.  PLEASE MODIFY SETTINGS AND/OR RE-RUN")				
	
	folderOutputLocation = os.path.dirname(igBlastFile)
		
	foutfile = open(outfile,'w')
	ferrorlog = open(outfile_errors,'w')
	
	if write_format == "TAB":		
		runningHeader = TABFileHeader()
	
	
	translator = DatabaseTranslator()
	translator[filetype_var] = write_format 
	translator_string = json.dumps(translator)
	
	translator_comment = descriptor_symbol#textFileCommentTranslator
	#translator_comment = textFileCommentTranslator #from cchrysostomou_databaseschema file
	#filetype_comment = textFileType+"json" #if the file type is a text file, it will be converted afterwards
	foutfile.write(translator_comment+translator_string+'\n')#in the first line of the file output a comment line definining how to convert this file into a database file 
	#foutfile.write(filetype_comment+'\n')	
			
	#read fasta file
	filename = inputFile;
	
	commandVal['Cdr3-Fr4Identification'] = json.dumps({'method':annotation_settings['cdr3_search']['method'],'param':cdr3_search_name})	 
	
	commandVal['Min_Alignment_Idendity_Cutoff'] = str(min_per_algn_cutoff)
	commandVal['Min_Alignment_Length_Cutoff'] =  str(min_algn_len_cutoff)
	
	if remove_insertions == 0:
		commandVal['Fix Mutations'] = 'Never' 
	elif remove_insertions == 1:
		commandVal['Fix Mutations'] = 'Always'
	elif remove_insertions == 2:
		commandVal['Fix Mutations'] = 'WhenStopCodon'
		
	commandString = json.dumps(commandVal)
	
	parse_igblast_notification = '''
	The resulting IgBlast Output will be parsed using the following settings:
	Minimum alignment length cutoff: {0},
	Minimum percent identity: {1},
	Fix detected mutations: {2},
	Method to identify CDR3 and start of fr4: {3}\n\n
	'''.format(str(min_algn_len_cutoff), str(min_per_algn_cutoff), commandVal['Fix Mutations'], json.dumps(annotation_settings['cdr3_search']['method']))
	
	print(parse_igblast_notification)
	
	useDebug = False		
	
	num_errors_found = 0;
	count = 0
	
	#seq_record = list(SeqIO.parse(testSet,"fasta"))
	startPer = 0
	
	
	chainDic = {
		'VH':["heavy","VDJ","IGH"],
		'VK':["light","VJ","IGK"],
		'VL':["light","VJ","IGL"],
		'TA':["alpha","VJ","TRA"],
		'TB':["beta","VDJ","TRB"]
	}
	
	readfastafile = readwrite.immunogrepFile(filelocation=filename,filetype=input_file_type)
	
	debugMe = False
	useDebug = False
		
	igblast_seq_found = True
	
	total_parsing_errors = 0

	generator_for_reading_igblastfile = ReadIgBlastQueryBlock(igBlastFile)
	
	while count<numSeqs:	#read through every sequence in the file 
		count += 1
		#allow the user to monitor what percent of the sequences have been processed					
		startPer = useful.LoopStatus(count,numSeqs,10,startPer,div=None)
																						
		fasta_line = readfastafile.IFclass.read()
		
		if fasta_line == None: #if end of file but didnt know that from before/accident parsing file 
			continue
				
		seqHeader = fasta_line[header_name]
		
		if seq_name  in fasta_line and fasta_line[seq_name]!="":
			seq = fasta_line[seq_name]
			
		else: #there is no sequence, so there will be no igblast output, skip the sequence  
			[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)
			seqDic["NOTES"] = "Sequence not provided"
			seqDic["ERRORS"] = "Sequence not provided"
			seqDic["PERCENT_IDENTITY"] = None
			seqDic["ALIGNMENT_LENGTH"] = None
			#if write_format=="TAB":				
				#for checkNewVals in newVals:				
				#	if checkNewVals not in runningHeader:					
				#		newColumn = len(runningHeader)+1
				#		runningHeader[checkNewVals] = newColumn									
			Write_Seq_JSON(seqDic,foutfile)			
			print "Seq # "+str(count)+" did not include a sequence: "+seqHeader			
			continue						
		
		seq_found = False										
		
		
		if igblast_seq_found: #In the previous iteration, the esequence in the fasta file was found in the igblast file. So we can read in the next block of info from the igblast file
			#READ an entire result block from igblast				
			[igblast_query_block,igblast_query_seq_header] = generator_for_reading_igblastfile.next() 
				
		if seqHeader.replace(" ","") != igblast_query_seq_header.replace(" ",""): #sequence header does not match the current igblast query block, go to the next block then
			[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)#IgBlastDef(seq,seqHeader)	
			seqDic["NOTES"] = "Sequence header not found in igblast"
			seqDic["ERRORS"] = "Sequence header not found in igblast"
			
			#if write_format=="TAB":				
#				for checkNewVals in newVals:				
#					if checkNewVals not in runningHeader:					
#						newColumn = len(runningHeader)+1
#						runningHeader[checkNewVals] = newColumn								
			Write_Seq_JSON(seqDic,foutfile)				
			print "Seq # "+str(count)+" not present in IgBlast file: "+seqHeader		
			igblast_seq_found = False
			continue
																							
		igblast_seq_found = True #the current sequence was found in igblast, so when done parsing this sequence, continue reading igblast in the next iteration
		
		notes=""
		parsing_error = ""
		t = {}
		perIden = 0.0 
		algnLen = 1
		t['STRAND']='+'
		matchVals = 0
		
		cdr3=""
		cdr3start = -1
		jstart = -1
		vgenestart = -1
		vend = -1
		alignmentFound=False
		antibody_alignV={}
		
		igblast_hit_field = {}
		type_list = {}							
		
		#FOR EACH SESQUENCE, EXTRACT THE RELEVANT INFO FROM THE IGBLAST FILE/block of results
		try:
			#SEQUENCE HEADER FOUND IN IGBLAST files
			query_results = igblast_query_block.split('\n') #break up results line by line into list
			var_type = 'UNK'
			line_num = 0
			no_results = False
			while line_num < len(query_results):
				
				line_info = query_results[line_num].strip('\r\n')						
				
				if line_info == '# 0 hits found':
					print "ok"
					no_results = True
									
				elif len(line_info)>0 and line_info[0]=="#": #hash tags refer to new section 
					tempLine = line_info.split(' ')
										
					#LIST THE TOP V,D,AND J GENES
					#DETERMINE IF PRODUCTIVE
					if tempLine[1]+" "+tempLine[2]=="V-(D)-J rearrangement": #-> this line indicates that we are at the section called "V-(D)-j rearragenemtn summary..."										
						line_num += 1
						topHits = query_results[line_num].split('\t')
																
						if len(topHits)==7:					
							var_type = "VJ"															
							if topHits[2]!="N/A":
								t['CHAIN'] = topHits[2];										
							
							#if topHits[3] == "Yes": ---> we will look for stop codons our self, at the end
							#	notes=notes+"Contains Stop Codon; "									
							if topHits[4] != "N/A":
								notes=notes+"Frame: "+topHits[4]+"; "
							else:
								notes=notes+"Frame: Unknown; "						
							if topHits[5] != "N/A":
								t['PRODUCTIVE'] = topHits[5].strip()
							else:
								t['PRODUCTIVE'] = "Unknown"																
							t['STRAND'] = topHits[6].strip();																				
							
						elif len(topHits)==8:						
							var_type = "VDJ"											
							
							if topHits[3]!="N/A":
								t['CHAIN'] = topHits[3];					
																										
							if topHits[5] != "N/A":
								notes=notes+"Frame: "+topHits[5]+"; "
							else:
								notes=notes+"Frame: Unknown; "
								
							if topHits[6] != "N/A":
								t['PRODUCTIVE'] = topHits[6].strip()
							else:
								t['PRODUCTIVE'] = "Unknown"											
																							
							t['STRAND'] = topHits[7].strip(); 																														
												
						else:					
							print "problems"
						
						if t['STRAND'] == "-":
							algn_seq = str(Seq(seq,generic_dna).reverse_complement())
						else:
							algn_seq = seq # str(Seq(seq,generic_dna))
						
						
					#GET DETAILED INFO ABOUT THE CDR3/V(D)-J JUNCTION AREA	
					elif tempLine[1]+" "+tempLine[2]=="V-(D)-J junction": #-> this line indicates that we are at the section called "V-(D)-j junction details..."																	
					#	vdj_cdr3seq = f1.readline().strip().split('\t')
						
						line_num += 1
						vdj_cdr3seq = query_results[line_num].strip().split('\t')
						
						initial_cdr3_len = 0
						
						for cdr3_substring in vdj_cdr3seq[1:]:
							if cdr3_substring.strip() != "N/A":
								cdr3_substring = cdr3_substring.replace(')','')
								cdr3_substring = cdr3_substring.replace('(','')
								initial_cdr3_len+=len(cdr3_substring)
														
							
					#GET INFORMATION ABOUT THE ANTIBODY ->FR1,FR2,CDR1-3,ETC
					elif tempLine[1]+" "+tempLine[2]=="Alignment summary": #-> this line inddicates that we are at the section, "Alignment summary ..." It lists all the relevant region of antibody that are present					
						annotation_list = []					
						line_num += 1
						
						tempLine = query_results[line_num].strip().split('\t')
						
						#tempLine = f1.readline().strip().split('\t')				
						alignmentFound = True
											
												
						while tempLine[0]!="Total": 					
							
							region=tempLine[0].split('-')#region[0] can be either, "FR1","CDR1","FR2","CDR2","FR3","CDR3"								
							
							for temp_count in range(len(tempLine)):
								if tempLine[temp_count] == "N/A":
									tempLine[temp_count] = 0
													
							start = int(tempLine[1])
							end = int(tempLine[2])
													
							#store all of annotated regions in a list. i.e. [{'name':'VREGION.FR1','start':3,....},{'name':'VREGION.CDR1'....}]
							annotation_list.append( 
								{								
									'name': region[0],
									'start': start,
									'end': end,
									'matches':int(tempLine[4]), #---> found an occurence  where matches and length did nto match the alignment discussed below in hit talbe, so we will do this manually
									'mismatches':int(tempLine[5]),
									'length':int(tempLine[3]),
									'gaps':int(tempLine[6]),
									'percent_identity':float(tempLine[7])/float(100),
									'frame':1 if seq_type == 'AA' else (start-1)%3+1
								}						
							)
							
							
							#start storing the positions of the alignment for each region of the antibody 						
							if region [0] == "FR3":
								cdr3start = end+1 #just in-case CDR3 below is not listed , this ahppens if igblast detects no "germline-V-gene"  to the CDR3 
								vd_start = cdr3start
								
							if region[0]=="CDR3": #this will replace the above code if CDR3(germline) is listed 
								cdr3start = start
								vd_start = end+1
						
							else:						
								matchVals+=int(tempLine[4])
								algnLen+=int(tempLine[3])							
								
						#	tempLine = f1.readline().strip().split('\t')																			
							line_num +=1
							tempLine = query_results[line_num].strip().split('\t')
							
						#algnLen = int(tempLine[3])
						numGaps = int(tempLine[6]) # total gaps from the line "total"
						algnNumMut = int(tempLine[5])+int(tempLine[6])
						#perIden = float(tempLine[7])/100
						totalLenReport = int(tempLine[3]) #we will use this later on to make sure the length is the same as the alignment  
					
					#DETAILED ALIGNMENT DATA
					#GET INFORMATION FOR THE JGENE START HERE
					#FIND THE END OF THE CDR3 USING HIT TABLE
					elif tempLine[1]+" "+tempLine[2]=="Hit table":#-> this line indicates that we are at the section called "Hit Table". It lists all of the info defined by outfmt -7.  Also stores the actual alignment between the two sequences
						line_num +=1
						fields = query_results[line_num].strip().split(',')
						line_num +=1
						numhits = query_results[line_num].strip().split(' ')						
						
						numhits = int(numhits[1])
											
						align_hits = {}
						
						#parse the "fields" line to find out which columns contain our fields of interest
						if not(igblast_hit_field):
							#Igblast returns a section called: ' Hit table (the first field indicates the chain type of the hit) '
							## Fields: subject id, query length, subject length, q. start, q. end, s. start, s. end, query seq, subject seq, evalue, alignment length, identical, mismatches, gaps, query frame, sbjct frame, subject strand
							#this just parses the fields returned by igblast so we know which column/tab returns the parameter we want				
							desired_fields={ #keys are the field names listed by igblast, values are the fields names I use in function
								'subject id':['gene_name','str'],
								'q. start':['start','int'],
								'q. end':['end','int'],
								's. start':['germline_start','int'],
								's. end': ['germline_end','int'],
								'subject length':['germline_len','int'],
								'alignment length': ['algn_len','int'],
								'identical':['matches','int'],
								'mismatches':['mismatches','int'],
								'gaps':['gaps','int'],
								'query seq':['query_seq','str'],
								'subject seq': ['germline_seq','str']
							}
																						
							fields_list = fields
							fields_list[0] = fields_list[0].replace('# Fields: ','')
							for i in range(len(fields_list)):
								fields_list[i] = fields_list[i].strip()
								
								if fields_list[i] in desired_fields:
									igblast_hit_field[ desired_fields[fields_list[i]][0] ] = i+1
									type_list[ desired_fields[fields_list[i]][0] ] = desired_fields[fields_list[i]][1]
									
						#FIND THE NON TOP-HIT V,D,AND J GENES
						#FIND OUT WHERE THE V GENE ENDS AND J GENE STARTS 
						for i in range(numhits):				
							line_num += 1
							algnSum = query_results[line_num].strip().split('\t')
						#	algnSum = f1.readline().strip().split('\t')				
							abregion = algnSum[0] #abregion would be V, D, or J (its the first column of "hittable")
												 
							if not(abregion in align_hits): #if havent defined V,D, or J yet 												
								align_hits[abregion] = [] #start as empty list
							
							temp_dict = {}
							
							for vals in igblast_hit_field:
								temp_dict[vals] = algnSum[igblast_hit_field[vals]] # iglbast_hit_field[vals] corresonds to the column location of the field or "vals" we want
								if type_list[vals]=='int':
									temp_dict[vals] = int(temp_dict[vals])
							
																	  
							temp_dict['percent_identity'] = temp_dict['matches']/float(temp_dict['algn_len'])
							
							align_hits[abregion].append(temp_dict)
				
				line_num += 1
			
			if no_results:
				[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)
				seqDic["NOTES"] = "No matches were found to germlines"				
				seqDic["PERCENT_IDENTITY"] = None
				seqDic["ALIGNMENT_LENGTH"] = None
				#if write_format=="TAB":				
#					for checkNewVals in newVals:				
#						if checkNewVals not in runningHeader:					
#							newColumn = len(runningHeader)+1
#							runningHeader[checkNewVals] = newColumn									
				Write_Seq_JSON(seqDic,foutfile)							
				continue			
				
		except Exception as e:
			total_parsing_errors+=1					
			exc_type, exc_obj, exc_tb = sys.exc_info()    	
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		
			print("An error occurred when parsing the igblast file for sequence # "+str(count)+". This error has been reported in the error log file. Continuing with analysis")
			error_to_file = "Error Number: "+str(total_parsing_errors)
			error_to_file += "\nSequence count: "+str(count)
			error_to_file += "\nSequence header: "+seqHeader
			error_to_file += "\nSequence: "+seq
			error_to_file += "\n\nError Message: "+str(e)
			error_to_file += "\nLine number: "+str(exc_tb.tb_lineno)
			error_to_file += "\nError type: "+str(exc_type)
			error_to_file += "\nFunction name: "+ str(fname)
			error_to_file += "\n\nIgBlast Output for Sequence:\n"
			error_to_file += igblast_query_block
			error_to_file += "\n**End of Error message**\n\n"
			ferrorlog.write(error_to_file)
			
			[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)
			seqDic["NOTES"] = "Error"
			seqDic["ERRORS"] = "This query generated an error (see error log file): "+str(e)			
			
			#if write_format=="TAB":				
#				for checkNewVals in newVals:				
#					if checkNewVals not in runningHeader:					
#						newColumn = len(runningHeader)+1
#						runningHeader[checkNewVals] = newColumn									
			Write_Seq_JSON(seqDic,foutfile)			
			
			continue
								
		try:
			if seq_type == 'AA':
				algn_seq = seq
				t['CHAIN'] = ''
				t['PRODUCTIVE'] = ''
				t['STRAND'] = '+'
				var_type = 'VDJ'
			
			recomb_type = var_type
			var_type = '' #if we want to go back to the old system, then make var_type = var_type+'.'
				
			#request the default variables from the IgBlastDef function
			#makes a new dictionary with these values defined as None
			[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)#IgBlastDef(seq,seqHeader)	
			
			perIden = matchVals/float(algnLen)
			
			seqDic["PERCENT_IDENTITY"] = perIden
			seqDic["ALIGNMENT_LENGTH"] = algnLen
																					
			if bad_dna_char_pattern.search(seq): #CHECK TO SEE IF SEQUENCE HAS WEIRD CHARACTERS IN IT						
				notes = "Sequence contains unknown characters"
				print "Seq # "+str(count)+" contains unusual characters: "+seqHeader
				seqDic["PERCENT_IDENTITY"] = None
				seqDic["ALIGNMENT_LENGTH"] = None
				seqDic['NOTES'] = notes
				seqDic['ERRORS'] = notes
					
#				if write_format=="TAB":				
#					for checkNewVals in newVals:				
#						if checkNewVals not in runningHeader:					
#							newColumn = len(runningHeader)+1
#							runningHeader[checkNewVals] = newColumn						
				
				Write_Seq_JSON(seqDic,foutfile)
	
				continue							
			
			seqDic['RECOMBINATION_TYPE'] = recomb_type
			if 'CHAIN' in t:
				if t['CHAIN'] in chainDic:
					seqDic[var_type+'PREDICTED_CHAIN_TYPE'] = chainDic[t['CHAIN']][0]					
					locus = [chainDic[t['CHAIN']][2]]
					seqDic[var_type+'LOCUS_NAME'] = locus[0]
				else:
					seqDic[var_type+'PREDICTED_CHAIN_TYPE'] = 'UNK'
					seqDic[var_type+'LOCUS_NAME'] = 'UNK'
					locus=[]
					notes+= 'Warning: Unknown chain type;' 
			else:
				locus = []
			
			
			#this will parse each annotated region of the sequence and find the position of nucleotides, ins/del, and translate result  
			if perIden>=min_per_algn_cutoff and algnLen>=min_algn_len_cutoff:			
				
				#highly condition/rare occurrence for this if statement
				if totalLenReport != align_hits['V'][0]['algn_len']: #---> found an occurence  where matches and length did nto match the alignment discussed below in hit talbe, so we will do this manually			
					start = align_hits['V'][0]['start']
					count_char = 0
					
					for region in range(len(annotation_list)):				
						#annotation_list[region]['start'] = start
						end = annotation_list[region]['end']
						matches = 0
						mismatches = 0
						gaps = 0
						length = 0
					
						while start<end+1:
							length+=1
							
							if align_hits['V'][0]['germline_seq'][count_char]=='-':
								gaps+=1
								start+=1
							elif align_hits['V'][0]['query_seq'][count_char] == '-':
								gaps+=1
							elif align_hits['V'][0]['query_seq'][count_char] == align_hits['V'][0]['germline_seq'][count_char]:
								start+=1
								matches+=1
							else:
								start+=1
								mismatches +=1
							count_char+=1
						
						annotation_list[region]['length'] = length
						annotation_list[region]['matches'] = matches
						annotation_list[region]['mismatches'] = mismatches
						annotation_list[region]['gaps'] = gaps
						
						if length>0:
							annotation_list[region]['percent_identity']=matches/float(length)
						else:							
							annotation_list[region]['percent_identity'] = 0
							
						annotation_list[region]['frame']=(annotation_list[region]['start']-1)%3+1
				
				[start_of_ab,start_frame] = Find_Start_Pos(annotation_list,align_hits['V'][0],force_to_v_germ_beg,seqtype=seq_type) #find the beginnin of the antibody using top 'V gene' hits 			
				
				annotation_list[0]['length'] = annotation_list[0]['length']+(annotation_list[0]['start']-start_of_ab)
				annotation_list[0]['start'] = start_of_ab
				annotation_list[0]['frame'] = start_frame
											
				#find end of antibody
				if 'J' in align_hits:					
					end_of_ab = Find_End_Pos( align_hits['J'][0],len(seq),force_to_j_germ_end,seqtype=seq_type  ) #find the end of the antibody using the top 'J' hits
				elif 'D' in align_hits:			
					end_of_ab = Find_End_Pos( align_hits['D'][0],len(seq),False,seqtype=seq_type ) #no j gene found, so find the end of the antibody using the top 'D' hits
				elif 'V' in align_hits:			
					end_of_ab = Find_End_Pos( align_hits['V'][0],len(seq),False,seqtype=seq_type ) #no d or j gene found, so find the end of the antibody using the top 'V' hits
				
				querySeq = algn_seq[:align_hits['V'][0]['start']-1] + align_hits['V'][0]['query_seq']
				germSeq = 'n'*(align_hits['V'][0]['start']-1) + align_hits['V'][0]['germline_seq']
																
				ab_seq_nt = ''
				ab_seq_aa =''
				
				eachFrame = []		
				start_align = start_of_ab-1
				mutation_notes = ''
				for annotated_region in annotation_list: # loop through each of the antibodies 
					startPos = annotated_region['start']
					endPos = annotated_region['end']
					name = annotated_region['name']
																								
					if name !='CDR3': #for the non-cdr3 part 
						
						end_align = start_align+annotated_region['length']-1
						frame = 1 if seq_type=='AA' else (startPos-1)%3 + 1
						
						[nt,gap_corrected_nt,aa,gap_corrected_aa,gap_adjusted] = GetAlignment(start_align,end_align,startPos,remove_insertions,querySeq,germSeq,seqtype=seq_type)
						 
						if seq_type!='AA':
							ab_seq_nt+=nt
						else:
							ab_seq_aa+=aa
																
						if gap_adjusted:
							mutation_notes+=name+','
																		
						#store FR1,CDR1,FR2,CDR2,FR3 in dictionary						
						seqDic['{0}VREGION.{1}.NT'.format(var_type,name)] = nt	#i.e. VDJ.VREGION.FR1.NT
						seqDic['{0}VREGION.{1}.AA'.format(var_type,name)] = aa	
						
						seqDic['{0}VREGION.{1}.NT.GAPPED'.format(var_type,name)] = gap_corrected_nt	#i.e. VDJ.VREGION.FR1.NT
						seqDic['{0}VREGION.{1}.AA.GAPPED'.format(var_type,name)] = gap_corrected_aa	
						
						seqDic['{0}VREGION.{1}.FRAME'.format(var_type,name)] = frame	
						seqDic['{0}VREGION.{1}.NT.Per_Identity'.format(var_type,name)] = annotated_region['percent_identity']
						eachFrame.append(frame)				
						start_align = end_align+1
						remaining_bases = endPos+1
									
				if mutation_notes != '':
					tempNotes = list(mutation_notes)
					tempNotes[-1] = ';'
					mutations_notes = ''.join(tempNotes)					
					notes+="Insertions/Deletions to germline have been removed from: {0}; ".format(mutation_notes)						
				
				if seq_type == 'AA':
					ab_seq_aa += algn_seq[remaining_bases-1:]
					ab_seq_nt = ''
				else:
					ab_seq_nt +=  algn_seq[remaining_bases-1:end_of_ab]
					
												
				##CDR3 IDENTIFICATION###
				if 'J' in align_hits and cdr3start>-1: #must have a J gene identified to find a CDR3 					
					cdr3_index = vd_start-1
					cdr3_index += initial_cdr3_len	
					
					[cdr3_nt,cdr3_aa,fr4_nt,fr4_aa,cdr3Frame,fr4Frame,cdr3notes] = CDR3_SEARCH_FUNCTION(algn_seq,cdr3start,end_of_ab,cdr3_index,cdr3_search_parameters,locus,var_type,align_hits['J'])																	
					notes+=cdr3notes
					seqDic['{0}CDR3.NT'.format(var_type)] = cdr3_nt
					seqDic['{0}CDR3.AA'.format(var_type)] = cdr3_aa	
					seqDic['{0}JREGION.FR4.NT'.format(var_type)] = fr4_nt
					seqDic['{0}JREGION.FR4.AA'.format(var_type)] = fr4_aa	
												
					seqDic['{0}CDR3.FRAME'.format(var_type)] = cdr3Frame	
					seqDic['{0}JREGION.FR4.FRAME'.format(var_type)] = fr4Frame																
				#######
				
				###just store any peritinent remaining variables into the dictionary
				seqDic[var_type+"PREDICTED_AB_SEQ.NT"] = ab_seq_nt
				
				seqDic[var_type+"PREDICTED_AB_SEQ.AA"] = ab_seq_aa if seq_type == 'AA' else str(Seq(ab_seq_nt,generic_dna).translate())#.tostring()		
				
				fullseq_has_stop = seqDic[var_type+"PREDICTED_AB_SEQ.AA"].find('*')		
				
				if fullseq_has_stop>-1:
					seqDic[var_type+'STOP_CODONS'] = True
				else:
					seqDic[var_type+'STOP_CODONS'] = False		
								
				#store all info on each gene name detected and percent identity 
				for detected_germline in align_hits: #detected_germline will either be V, D, or J (see above)
					seqDic['{0}{1}REGION.{1}GENES'.format(var_type,detected_germline)] = [] #i.e. VDJ.VREGION.VGENES
					seqDic['{0}{1}REGION.{1}GENE_SCORE'.format(var_type, detected_germline)] = []#i.e. VDJ.VREGION.VGENE_SCORE
					for gene in align_hits[detected_germline]:
						seqDic['{0}{1}REGION.{1}GENES'.format(var_type,detected_germline)].append(gene['gene_name']) #append all gene names to this list
						seqDic['{0}{1}REGION.{1}GENE_SCORE'.format(var_type,detected_germline)].append(gene['percent_identity']) #append the mtuation profile/percent identity to list
						
					seqDic['{0}{1}GENE_ALIGNMENT_SEQUENCE'.format(var_type,detected_germline)] = align_hits[detected_germline][0]['query_seq']
					seqDic['{0}{1}GENE_ALIGNMENT_GERMLINE'.format(var_type,detected_germline)] = align_hits[detected_germline][0]['germline_seq']
																	
				if len(annotation_list)==6 and 'D' in align_hits and 'J' in align_hits:
					seqDic[var_type+"FULL_LENGTH"] = True
				else:
					seqDic[var_type+"FULL_LENGTH"] = False					
											
				if seq_type == 'AA':
					seqDic[var_type+'VREGION.SHM_AA'] = str(algnNumMut) if seq_type == 'AA' else ''
					seqDic[var_type+'VREGION.SHM_AA_PER'] = round(100*algnNumMut/float(algnLen),2)
				else:
					seqDic[var_type+'VREGION.SHM_NT'] =  str(algnNumMut) if seq_type != 'AA' else ''				
					seqDic[var_type+'VREGION.SHM_NT_PER'] = round(100*algnNumMut/float(algnLen),2)
				
				seqDic[var_type+'VREGION.NUM_GAPS'] = numGaps
				
				
				seqDic['COMMAND'] = commandString		
				seqDic['ORIENTED_SEQ'] = algn_seq
			else:			
				notes='Sequence results did not pass alignment filter settings declared in app;'+notes				
						
			if 'PRODUCTIVE' in t:
					seqDic[var_type+"PRODUCTIVE"]=t['PRODUCTIVE']		
			
			if 'STRAND' in t:
				seqDic[var_type+'STRAND'] = t['STRAND']			
			
			
	
			seqDic["NOTES"] = notes
			
				
			
		except Exception as e:
			total_parsing_errors+=1					
			exc_type, exc_obj, exc_tb = sys.exc_info()    	
			fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	
			print "An error occurred when analyzing the igblast output for sequence # "+str(count)+". This error has been reported in the error log file. Continuing with analysis"
			error_to_file = "Error Number: "+str(total_parsing_errors)
			error_to_file += "\nSequence count: "+str(count)
			error_to_file += "\nSequence header: "+seqHeader
			error_to_file += "\nSequence: "+seq
			error_to_file += "\n\nError Message: "+str(e)
			error_to_file += "\nLine number: "+str(exc_tb.tb_lineno)
			error_to_file += "\nError type: "+str(exc_type)
			error_to_file += "\nFunction name: "+ str(fname)
			error_to_file += "\n\nIgBlast Output for Sequence:\n"
			error_to_file += igblast_query_block
			error_to_file += "\n**End of Error message**\n\n"
			ferrorlog.write(error_to_file)	
		
		
			[seqDic,newVals] = IgBlastDef(fasta_line,seq_name,header_name,quality_name)
			seqDic["NOTES"] = "Error"
			seqDic["ERRORS"] = "This query generated an error (see error log file): "+str(e)			
		
			#if write_format=="TAB":				
#				for checkNewVals in newVals:				
#					if checkNewVals not in runningHeader:					
#						newColumn = len(runningHeader)+1
#					runningHeader[checkNewVals] = newColumn									
			Write_Seq_JSON(seqDic,foutfile)			
		
			print "Seq # "+str(count)+" did not include a sequence: "+seqHeader
															
			continue
		
		#if write_format=="TAB":				
#			for checkNewVals in newVals:				
#				if checkNewVals not in runningHeader:					
#					newColumn = len(runningHeader)+1
#					runningHeader[checkNewVals] = newColumn																					
		
		Write_Seq_JSON(seqDic,foutfile)		
							
	foutfile.close()
	ferrorlog.close()	
	
	if write_format=="TAB":
		includeDecorators=True		
		print("The IgBlast analysis has completed and the summarized file has been generated. The file will now be converted to a TAB delimited file...")
		writeErrorsFound=JSON_to_TXT(outfile,outfile,includeDecorators,runningHeader)	
		if writeErrorsFound:			
			print("Error: An Error occurred when trying to write the text tab delimited file.  But, the IgBlast Alignment File should still be present")
				
	if total_parsing_errors > 0:
		print("The analysis has completed. However, {0} of the {1} total sequences analyzed contained errors when parsing the file. Please refer to the error log file to debug possible errors".format(str(total_parsing_errors), str(numSeqs)))
		return [outfile,igBlastFile,outfile_errors] 
	else:
		os.system("rm '{0}'".format(outfile_errors))		
		return [outfile,igBlastFile]
	
	

def Check_Frame(seq):
	nt_length = len(seq)%3
	if nt_length==1:
		seq = "NN"+seq
	elif nt_length==2:
		seq = "N"+seq
	return seq

#find the beginning of the antibody sequence
def Find_Start_Pos(annotation_list,v_gene_info,adjust_start,seqtype='NT',deb=False):
	#assumes the next annotated region is in frame.  i.e. if there is FR1, CDR1, then assumes taht CDR1 is in frame relative to the antibody sequence, but FR1 may not be.  
	seq = v_gene_info['query_seq']
	
	seq_start = v_gene_info['start'] 
	
	#frame = annotation_list[0]['frame']
			
	germline_seq = v_gene_info['germline_seq']
	
	algnLen = annotation_list[0]['length'] #length of this alignment
	
	substr_seq = seq[:algnLen]
	substr_germ = germline_seq[:algnLen]
	
	numGaps_seq = substr_seq.count('-')
	numGaps_germ = substr_germ.count('-')
	
	framework_start = v_gene_info['germline_start']		
	
	framework_end = framework_start+algnLen-1
	framework_end = framework_end - numGaps_germ
	
	if adjust_start: #adjust_start means the user wants to force the alignment of the seqence to the front of the germline if there are nucleotides upstream
		seq_start = seq_start-(framework_start-1)
		framework_start = 1
	
	frame = 1 if seqtype=='AA' else (framework_start-1)%3+1
	
	next_region = framework_end+1 #start of the next annoatedregion. we assume this region is in frame

	germline_frame = 1 if seqtype=='AA' else (next_region-1)%3+1 #this should be the actual frame of the germline sequence based on the alignment  (ususally, but not always, germline should be frame 1)

	shift_start_pos = frame-germline_frame #this tells us how to shift the sequence 
	
	
	seq_start = seq_start - shift_start_pos

	#if adjust_start: 
	#	seq_start = seq_start - (framework_start - germline_frame)

	
	while seq_start<1: #make sure we are in the poisitive positions
		seq_start += 1 if seqtype=='AA' else 3
	frame = 1 if seqtype=='AA' else (seq_start-1)%3+1		
		
	return seq_start,frame
	
def Find_End_Pos( algn_info,  seq_len, adjust_end,seqtype='NT'):		
				
	if adjust_end: #the user wants to force the alignment to the end of the predicted j-germline
		extraPos = algn_info['germline_len'] - algn_info['germline_end'] 		
	else:
		extraPos = 0
	
	endPos = algn_info['end'] + extraPos
		
	while endPos>seq_len:
		endPos-= 1 if seqtype=='AA' else 3
	
	return endPos	
	


############TRANSLATE ALIGNMENT OF ANTIBDOY TO GERMLINE INTO AMINO ACID AND CORRECT FOR STOPCODON FORMATION (if requested)############
#remove_insertions: 0 -> dont remove insertions, 1 -> remove insertions always, 2-> remove insertions only if there is a stop codon
def GetAlignment(start_align,end_align,nt_query_pos,remove_insertions,querySeq,germSeq,deb=False,seqtype='NT'):
	
	query_substring = str(querySeq[start_align:end_align+1])
	germ_substring = str(germSeq[start_align:end_align+1])

	num_mutated = 0
	if deb:
		print 'seq_info'
		print query_substring
		print germ_substring
		
	ins_event = []
	del_event = []
	
	ins_event_substring = []
	del_event_substring = []
	gap_adjusted = False			 
	
	nt_seq_gappless = query_substring
	
	if ('-' in query_substring) or ('-' in germ_substring):

		#go through tthe alignment looking for gaps	
		for char_pos in range(len(germ_substring)):		
			if germ_substring[char_pos] == '-':  #there is a gap int he germline sequence, must be due to insertion
				ins_event.append(nt_query_pos)
				ins_event_substring.append(char_pos)
				nt_query_pos+=1			
			elif query_substring[char_pos] == '-': #there is a gap in the query sequence, must be due to a deletion event 
				del_event.append(nt_query_pos)
				del_event_substring.append(char_pos)
			else:
				nt_query_pos+=1		
		
		if deb:
			print 'mutations'
			print ins_event_substring
			print del_event_substring
			
		if remove_insertions == 0:			
			nt_seq_gappless = nt_seq_gappless.replace('-','')	
			nt_seq_gap = query_substring.replace('-','n')
		elif remove_insertions == 1:
			nt_seq_gappless = list(query_substring)
			if deb:
				print nt_seq_gappless
			for j in ins_event_substring:
				nt_seq_gappless[j] = ''
				num_mutated -= 1
			
			if deb:
				print nt_seq_gappless
			
			nt_seq_gappless = ''.join(nt_seq_gappless)			
									
			nt_seq_gap = nt_seq_gappless.replace('-','n')
			nt_seq_gappless = nt_seq_gappless.replace('-','')
			gap_adjusted = True
		
		if seqtype == 'AA':
			aa_seq_gappless = nt_seq_gappless
			nt_seq_gappless = ''
		else:
			aa_seq_gappless = str(Seq(nt_seq_gappless,generic_dna).translate())#.tostring()
			
		if remove_insertions == 2 and '*' in aa_seq_gappless:
			nt_seq_gappless = list(nt_seq_gappless)
			for j in ins_event_substring:
				nt_seq_gappless[j] = ''
				num_mutated -= 1
			
			nt_seq_gappless = ''.join(nt_seq_gappless)			
			
			nt_seq_gap = nt_seq_gappless.replace('-','n')
			nt_seq_gappless = nt_seq_gappless.replace('-','')
			
			aa_seq_gappless = str(Seq(nt_seq_gappless,generic_dna).translate())#.tostring()
			gap_adjusted = True
		
		if seqtype=='AA':
			aa_seq_gap = nt_seq_gap
			nt_seq_gap = ''
		else:
			aa_seq_gap = str(Seq(nt_seq_gap,generic_dna).translate())#.tostring()								
	
	else: #no gaps detected in this region. 			
		
		if seqtype=='AA':
			nt_seq_gap = None
			aa_seq_gap = None
			aa_seq_gappless = nt_seq_gappless
			nt_seq_gappless = ''
		else:
			aa_seq_gappless = str(Seq(nt_seq_gappless,generic_dna).translate())#.tostring()
			nt_seq_gap = None
			aa_seq_gap = None
		
	return nt_seq_gappless, nt_seq_gap,aa_seq_gappless, aa_seq_gap,gap_adjusted



	
####################################################################################################	

##########SUPPLMENTARY FUNCTIONS FOR DEFINING VARIABLE NAMES,  WRITING, AND READING IGBLAST DATA TO TEXT FILE########################

###########FUNCTIONS FOR TEXT FILES###########################################
#write the header row if using text tabe delimited format###################

def WriteTABFileHeader(headerDic, outfile, filetype):
	filename = open(outfile,'w')
	
	#first write the "translator" json line so that we can update the database with igblast results
	translator = DatabaseTranslator()
	translator_string = json.dumps(translator);
	translator_comment = descriptor_symbol#textFileCommentTranslator
	filename.write(translator_comment+translator_string+'\n')
	
	rowData = [None]*len(headerDic)
	
	for key in headerDic:
		rowData[headerDic[key]-1] = key
	
	for i in range(len(rowData)-1):
		filename.write(rowData[i]+'\t')
	
	filename.write(rowData[len(rowData)-1]+'\n')

#We will need to update the database with the results from IgBlast.  In order to update teh database, we need a translator
#So that we know what fields go where in the database
def DatabaseTranslator(input_dictionary = {}):
	key = translation_var
	translator = {				
			"ANALYSIS_NAME": "IGBLAST", # NAME OF THE ANALYSIS 
			"RECOMBINATION_FIELD":{ #THIS TELLS THE PROGRAM HOW TO DETERMINE WHETHER AN ANALYSIS/QUERY RESULT (from file) IS VDJ OR VJ
					"FIELD_NAME": "RECOMBINATION_TYPE", #name of the field in the file that will give information regarding the recombination type (VDJ OR VJ)
					"EXPLICIT":True, #IF EXPLICIT, then the RECOMBINATION_TYPE is equal to the EXACT value in this field (i.e. it will either list VDJ or VJ), IF false, then VDJ and VJ are defined by values in list below
					"INEXPLICIT_DEFINITIONS":{
						"VDJ":[],#if expliit is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VDJ type (i.e. if Locus is used to determine recombination type then it woudl be VDJ:['IGH']
						"VJ":[],#if expliit is FALSE, then this list MUST NOT be empty. It must list all values from this field that will correspond to a VJ type (i.e. if Locus is used to determine recombination type then it woudl be VJ:['IGK','IGL']
	 				},
			},
			
			
			"FIELDS_UPDATE":{ 
				#this will map all of the fields in the file to the proper location in the database. 
				#KEY = field name in database
				#VALUE = field name in file 
				#I.E. If I list VGENES as the column name/field name, then i want to map VGENES:VREGION.VGENES (because VREGION.VGENES is the name in the database)							
				#"COMMAND":"COMMAND",								
				idIdentifier:idIdentifier,
				"SEQUENCE":"FULL_SEQ",
				"COMMAND":"COMMAND",
				"SEQUENCE_HEADER":"SEQ_HEADER",				
				"QUALITY_SCORE":"QUALITY_SCORE",
				"FULL_LENGTH":"FULL_LENGTH",
				"STOP_CODONS":"STOP_CODONS",
				"PREDICTED_CHAIN_TYPE":"PREDICTED_CHAIN_TYPE",
				"PRODUCTIVE":"PRODUCTIVE",
				"NOTES":"NOTES",							
				"PREDICTED_AB_SEQ.NT":"PREDICTED_AB_SEQ.NT",
				"PREDICTED_AB_SEQ.AA":"PREDICTED_AB_SEQ.AA",
				"STRAND":"STRAND",
				"LOCUS_NAME":"LOCUS_NAME",
				"VREGION.SHM.NT":"VREGION.SHM_NT",
				"VREGION.FR1.NT":"VREGION.FR1.NT",
				"VREGION.FR1.AA":"VREGION.FR1.AA",
				"VREGION.CDR1.NT":"VREGION.CDR1.NT",			
				"VREGION.CDR1.AA":"VREGION.CDR1.AA",
				"VREGION.FR2.NT":"VREGION.FR2.NT",
				"VREGION.FR2.AA":"VREGION.FR2.AA",
				"VREGION.CDR2.NT":"VREGION.CDR2.NT",
				"VREGION.CDR2.AA":"VREGION.CDR2.AA",
				"VREGION.FR3.NT":"VREGION.FR3.NT",
				"VREGION.FR3.AA":"VREGION.FR3.AA",
				"VREGION.VGENES":"VREGION.VGENES",
				"VREGION.VGENE_SCORES":"VREGION.VGENE_SCORE",
				"CDR3.NT":"CDR3.NT",
				"CDR3.AA":"CDR3.AA",
				"DREGION.DGENES":"DREGION.DGENES",				
				"DREGION.DGENE_SCORES":"DREGION.DGENE_SCORE",
				"JREGION.FR4.NT":"JREGION.FR4.NT",
				"JREGION.FR4.AA":"JREGION.FR4.AA",
				"JREGION.JGENES":"JREGION.JGENES",
				"JREGION.JGENE_SCORES":"JREGION.JGENE_SCORE",								
		}
	}
	
	input_dictionary[key] = translator
	return input_dictionary


#KEYS REPRESENT THE NAME FO THE COLUMN, VALUES REPRESENT THE COLUMN NUMBER###
def TABFileHeader():
	igBlDic = {							
			"SEQ_HEADER":1,
			"DOCUMENTHEADER":2,
			"FULL_SEQ":3,
			"QUALITY_SCORE":4,
			"COMMAND":5,	
			"RECOMBINATION_TYPE":6,
			
			"PERCENT_IDENTITY":7,
			"ALIGNMENT_LENGTH":8,
						
			"PREDICTED_CHAIN_TYPE": 9,
			"PRODUCTIVE": 10,
			"NOTES": 11, 
			
			"LOCUS_NAME":12,
			"PREDICTED_AB_SEQ.NT":13,
			"PREDICTED_AB_SEQ.AA":14,
			"STRAND":15,			
			
			"VREGION.VGENES":16,
			"VREGION.VGENE_SCORE":17,
			"JREGION.JGENES":18,
			"JREGION.JGENE_SCORE":19,			
			"DREGION.DGENES":20,
			"DREGION.DGENE_SCORE":21,
			"VREGION.FR1.NT":22,
			"VREGION.FR1.AA":23,
			"VREGION.CDR1.NT":24,
			"VREGION.CDR1.AA":25,
			"VREGION.FR2.NT":26,
			"VREGION.FR2.AA":27,
			"VREGION.CDR2.NT":28,
			"VREGION.CDR2.AA":29,
			"VREGION.FR3.NT":30,
			"VREGION.FR3.AA":31,						
			
			
			"CDR3.NT":32,
			"CDR3.AA":33,
			"JREGION.FR4.NT":34,
			"JREGION.FR4.AA":35,			
			
			
			"VREGION.FR1.NT.GAPPED":36,
			"VREGION.FR1.AA.GAPPED":37,
			"VREGION.CDR1.NT.GAPPED":38,
			"VREGION.CDR1.AA.GAPPED":39,
			"VREGION.FR2.NT.GAPPED":40,
			"VREGION.FR2.AA.GAPPED":41,
			"VREGION.CDR2.NT.GAPPED":42,
			"VREGION.CDR2.AA.GAPPED":43,
			"VREGION.FR3.NT.GAPPED":44,
			"VREGION.FR3.AA.GAPPED":45,						
												
			"VREGION.SHM_NT":46,
			"VREGION.SHM_NT_PER":47,
			"VREGION.NUM_GAPS":48,
			
			"VREGION.FR1.NT.Per_Identity":49,			
			"VREGION.CDR1.NT.Per_Identity":50,			
			"VREGION.FR2.NT.Per_Identity":51,			
			"VREGION.CDR2.NT.Per_Identity":52,			
			"VREGION.FR3.NT.Per_Identity":53,			
			
			"VGENE_ALIGNMENT_SEQUENCE":54,
			"DGENE_ALIGNMENT_SEQUENCE":55,
			"JGENE_ALIGNMENT_SEQUENCE":56,
			"VGENE_ALIGNMENT_GERMLINE":57,
			"DGENE_ALIGNMENT_GERMLINE":58,
			"JGENE_ALIGNMENT_GERMLINE":59,
			"FULL_LENGTH":60,
			"STOP_CODONS":61,							
					
			'VREGION.FR1.FRAME': 62,
			'VREGION.CDR1.FRAME': 63,
			'VREGION.FR2.FRAME': 64,
			'VREGION.CDR2.FRAME': 65,
			'VREGION.FR3.FRAME': 66,
			'CDR3.FRAME': 67,
			'JREGION.FR4.FRAME': 68,
			'ORIENTED_SEQ':69,			
			'ERRORS': 70
			
			
			
	};
	
	numCols = len(igBlDic)
	
	#headerList = [None]*numCols
	
	#fill list with column anmes to write to
	#for headerName in igBlDic:		
	#	colNum = igBlDic[headerName]
	#	headerList[colNum-1] = headerName
			
	#for i in range(numCols-1):
	#	output_file.write(headerList[i]+'\t')
	
	#output_file.write(headerList[numCols-1]+'\n')
	
	return igBlDic
	
def Write_Seq_JSON(seqDic,foutfile):
	#convert dictionary to json-string
	st = json.dumps(seqDic,sort_keys=True)
	#output string in json format
	foutfile.write(st+'\n')
###############################################################################

##################DEFAULT DICTIONARY FOR IGBLAST FILES#######################
def DefaultIgBlastDic():

	igBlDic = {
			"FULL_SEQ":None,
			"SEQ_HEADER":None,
			"QUALITY_SCORE":None,
			"DOCUMENTHEADER":None,
			"COMMAND":None,
			"RECOMBINATION_TYPE":None,
			
			"PERCENT_IDENTITY":None,
			"ALIGNMENT_LENGTH":None,
			
			"PREDICTED_CHAIN_TYPE":None,
			"PREDICTED_CHAIN_TYPE":None,
			
			"NOTES":None,
			"PRODUCTIVE":None,			
			
			#heavy chain terms
			"VREGION.VGENES":None,
			"VREGION.VGENE_SCORE":None,
			"JREGION.JGENE_SCORE":None,
			"DREGION.DGENE_SCORE":None,			
			"VREGION.FR1.NT":None,
			"VREGION.FR1.AA":None,
			"VREGION.CDR1.NT":None,
			"VREGION.CDR1.AA":None,
			"VREGION.CDR2.NT":None,
			"VREGION.CDR2.AA":None,
			"VREGION.FR2.NT":None,
			"VREGION.FR2.AA":None,			
			"CDR3.NT":None,
			"CDR3.AA":None,			
			"VREGION.FR3.NT":None,
			"VREGION.FR3.AA":None,			
			"DREGION.DGENES":None,
			"JREGION.JGENES":None,			
			"JREGION.FR4.NT":None,
			"JREGION.FR4.AA":None,			
			"STRAND":None,						
			"PREDICTED_AB_SEQ.NT":None,
			"PREDICTED_AB_SEQ.AA":None,			
			"VREGION.SHM_NT":None,
			"VREGION.SHM_NT_PER":None,			
			"VREGION.NUM_GAPS":None,
			"VGENE_ALIGNMENT_SEQUENCE":None,
			"DGENE_ALIGNMENT_SEQUENCE":None,
			"JGENE_ALIGNMENT_SEQUENCE":None,			
			"VGENE_ALIGNMENT_GERMLINE":None,
			"DGENE_ALIGNMENT_GERMLINE":None,
			"JGENE_ALIGNMENT_GERMLINE":None,			
			"FULL_LENGTH":None,
			"STOP_CODONS":None,												
			"LOCUS_NAME":None,						
			"ERRORS": None
	};
	return igBlDic

#fasta_line => dictionary of all fields detected in fasta file
	#=> seqheader
	#=> sequence
	#=> and all fields detected after the fasta_file_delimiter as a json (see read function in immunogrepfile)
def IgBlastDef(fasta_line,seq_field,header_field,quality_field):
		
	
	igBlDic = DefaultIgBlastDic()
	
	igBlDic["FULL_SEQ"] = fasta_line.pop(seq_field,None)
	
	#SEQ header only returns the sequence header and ignores json fields after the fasta_file_delimiter. these fields are included as keys in fasta_line variable
	igBlDic["SEQ_HEADER"] = fasta_line.pop(header_field,None) 	
	
	#quality score if its a fastqfile (yes fasta_line is a bad varialbe name)
	igBlDic["QUALITY_SCORE"] = fasta_line.pop(quality_field,None)
	
	#documentheader returns the entire header line where the JSON string is included.
	igBlDic['DOCUMENTHEADER'] = fasta_line.pop('document_header') if 'document_header' in fasta_line else igBlDic['SEQ_HEADER']		
	
	if idIdentifier in fasta_line:
		igBlDic[idIdentifier] = fasta_line.pop(idIdentifier)
	
	#igBlDict = dict(igBlDic.items()+fasta_line.items()) #essential..this ensures that if seq_id is in the original file, it will be passed on TO THE OUTPUT/ANNOTATION FILE!! =L> commented out 

	return igBlDic, fasta_line
###############################################################################

#HEADERLINE VARIABLE#
#headerLine is a dictionary that defines the headerrow for the text file
#the key of the dictionary refers to the column name used for the text file
#and the value of each key is the "column number" for that column name
#if you do not have a headerLine defined, then pass in the variable as {}

#INCLUDE DECORATORS VARIABLE#
#we put comments into our analysis files using special text identifies that refer to "comments"
#if INCLUDE DECORATORS is set to true, then when converting the text file, it will also copy these "comment" lines in the datafile. If not, then it will skip them
def JSON_to_TXT(input_filename,output_filename,includeDecorators,headerLine): #converts a JSON file into a text tab delimited file
	try:
		a = str(datetime.now()).replace(' ','_')
		a = a.replace(':','-')
	
		if input_filename==output_filename:
			overwriteFile=True						
			f = open("scratch/temp_file_name_"+a,'w')							
		else:
			overwriteFile=False
			f = open(output_filename,'w')
			
		
		i = open(input_filename)	
		
		decoratorFound = True
		
		dictHeader = {}
		dictList = []
		currentRow = 1;
		headerRow = []
		
		error = False
		commentLine = descriptor_symbol
		lenComment = len(commentLine)
		
		if headerLine == {}:				
			
			while(decoratorFound):
				line_one = i.readline().strip()
				if line_one[:lenComment] != commentLine:
					decoratorFound = False		
							
			myDict = json.loads(line_one)					
			for keys in myDict:
				if not keys in dictHeader:
					dictHeader[keys] = currentRow
					headerRow.append(keys)
					currentRow+=1
			
			for line in i:
				line = line.strip()
			
				myDict = json.loads(line)					
				for keys in myDict:
					if not keys in dictHeader:
						dictHeader[keys] = currentRow
						headerRow.append(keys)
						currentRow+=1
			
		else:
			dictHeader = headerLine
			headerRow = [None]*len(dictHeader)		
			for keys in dictHeader:			
				#print keys
				#print dictHeader[keys]
				#print len(dictHeader)
				headerRow[dictHeader[keys]-1] = keys
		
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
		
		
		
		for header in headerRow:
			if header:
				f.write(header+'\t')
			else:
				f.write(''+'\t')
				
		f.write('\n')
		
		rowList = [None]*len(headerRow)
		line = line_one.strip()
		if line and line!="" and line[:lenComment]!=commentLine: #check to see if this line in the text file is a "comment line", if not, then write tofile
			myDict = json.loads(line)
		
			
			for keys in myDict:			
				if keys in dictHeader:
					rowList[dictHeader[keys]-1] = myDict[keys]
	
			for rowVal in rowList:
				if rowVal is not None:				
					if type(rowVal) is list:						
						f.write(','.join([str(r) for r in rowVal])+'\t')
					elif type(rowVal) is dict:
						f.write(json.dumps(rowVal)+'\t')
					else:											
						f.write(str(rowVal)+'\t')
				else:
					f.write('\t')
				
			f.write('\n')
		else:
			if includeDecorators:
				f.write(line_one+'\n')
			
		for line in i:	
			rowList = [None]*len(headerRow)
			line = line.strip()
			if line[:lenComment]!=commentLine:
				myDict = json.loads(line)
			
				for keys in myDict:					
					if keys in dictHeader:
						rowList[dictHeader[keys]-1] = myDict[keys]
						
				for rowVal in rowList:
					if rowVal is not None:				
						if type(rowVal) is list:
							f.write(','.join([str(r) for r in rowVal])+'\t')
						elif type(rowVal) is dict:
							f.write(json.dumps(rowVal)+'\t')
						else:											
							f.write(str(rowVal)+'\t')
					else:
						f.write('\t')
					
				f.write('\n')
			else:
				if includeDecorators:
					f.write(line+'\n')
		
		f.close()
		i.close()
		
		if overwriteFile:
			os.system("mv "+"scratch/temp_file_name_"+a+" "+input_filename)	
			
		return error
		
	except Exception as e:
		error = True
		
		exc_type, exc_obj, exc_tb = sys.exc_info()
    	
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
				
		print("there was an error when writing to file: "+str(e)+". Please Restart")		
		print("Line Number: " + str(exc_tb.tb_lineno) + " type: " + str(exc_type) + " fname: " + str(fname))
		
		return error
		
