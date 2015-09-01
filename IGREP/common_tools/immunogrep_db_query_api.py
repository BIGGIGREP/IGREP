#6/25/2015
#new version for querying: 
	#allows conversions to CSV, TAB, JSON, IGREP, 
	#allows requesting to search across anotation methods 
	#allows requesting to  reprot sequence adn sequence header for any seqs collection query 
import inspect
import sys
import json
from functools import wraps
from datetime import datetime
import random

import bson
from bson.json_util import dumps as bson_dumps 
from bson.json_util import loads as bson_loads 
from bson.objectid import ObjectId
import pymongo
import copy

import os
import time
import immunogrep_useful_functions as useful_functions
from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import expIdentifier
from immunogrep_global_variables import seqRawData
from immunogrep_global_variables import fasta_file_delimiter
import itertools
import immunogrep_database_schema as schema
from immunogrep_database_schema import convert_to_objectid
from immunogrep_database_schema import Exps_Collection
from immunogrep_database_schema import convert_text_to_index_field_text
from collections import defaultdict
import re
from collections import OrderedDict
import copy
import sets

from collections import namedtuple

import immunogrep_proxy_tools
import shutil


fields_above_data_key = ['_id','ANALYSIS_NAME','EXP_ID','RECOMBINATION_TYPE','SEQ_ID','DATE_UPDATED','SETTINGS']
data_fields_for_raw_data = ['SEQUENCE_HEADER','SEQUENCE','QUALITY_SCORE','FILENAME']
query_data_key_name = 'QUERY_DATA'

#processeses the output query faster by using multithreading
allow_multithreading=True
max_multithreading_processes = 6
chunk_size = 1

oid_type=type(ObjectId())
test2=type({})
#try:
#	import appsoma_api
#	appsoma_api.resource_pull("https://biotseq.ut.appsoma.com:5001/kamhonhoi_public/Programs/immunogrep_cython_db_tools.so","immunogrep_cython_db_tools.so")
#except:
#	pass

try:
	#cython code defining useful_functions	
	from immunogrep_cython_db_tools import flatten_dictionary
	from immunogrep_cython_db_tools import RemoveObjId				
	#print('aaahhh im commented out!!!! visual studio jazz!!!')
except:
	#it wont work if there is no cython module (currently immunogrep_cython_db_tools). when dit doesnt work, import non cyhton version
	from immunogrep_useful_functions import flatten_dictionary 
	from immunogrep_useful_functions import RemoveObjId
	print('Not using cython for flatten dictionary')


global_time =0
gt1=0
gt2=0
gt3=0
gt4=0
# use simplejson if available
# else use json
#try: import simplejson as json
#except ImportError: import json



#-----------------------------------------------------------------------------------------------------------------------------------------
# MongoDB connection
#-----------------------------------------------------------------------------------------------------------------------------------------
# opens the database connection.
#DELETE FUNCTION WHEN DONE!
def connectToIgDatabase(dbpath,username,password):
	print "aaaargconnecting7"					
	connection = pymongo.MongoClient( host=dbpath, port=27017 )	
	db = connection.appsoma			
	db.authenticate(username,password)		
	return db,connection



#NOTE1: PYMONO AGGREGATION SORT REQUIRES BSON SON OR ORDEREDDICT!!!!http://api.mongodb.org/python/current/examples/aggregation.html

#use this value to delimit "filename information" and "real data"
file_def_seperator = "}{" #filename}{>somerealdatacomes after

#a default variable controlling the preferred order of fields output to delimited files  
#also shows the default schema as of 6/26/2015
default_sorting_order = [	
	idIdentifier,
	'SEQUENCE_HEADER',
	'SEQUENCE',
	'PREDICTED_AB_SEQ.AA',
	'PREDICTED_AB_SEQ.NT',

	'STRAND',
	'QUALITY_SCORE',
	'PRODUCTIVE',
	'RECOMBINATION_TYPE',
	'PREDICTED_CHAIN_TYPE',
	'LOCUS_NAME',

	'VREGION.VGENES',
	'DREGION.DGENES',
	'JREGION.JGENES',	
	'FULL_LENGTH',
	'STOP_CODONS',
	
	'VREGION.FR1.AA',
	'VREGION.CDR1.AA',
	'VREGION.FR2.AA',
	'VREGION.CDR2.AA',
	'VREGION.FR3.AA',	
	'CDR3.AA',
	'JREGION.FR4.AA',
	
	'VREGION.FR1.NT',	
	'VREGION.CDR1.NT',
	'VREGION.FR2.NT',
	'VREGION.CDR2.NT',
	'VREGION.FR3.NT',	
	'CDR3.NT',
	'JREGION.FR4.NT',
	
	'VREGION.CDR1.AA_LENGTH',
	'VREGION.CDR1.NT_LENGTH',
	'VREGION.CDR2.AA_LENGTH',
	'VREGION.CDR2.NT_LENGTH',
	'CDR3.AA_LENGTH',
	'CDR3.NT_LENGTH',
	'VREGION.VGENE_SCORES',
	'JREGION.JGENE_SCORES',
	'DREGION.DGENE_SCORES',
	'VREGION.SHM.NT',
	'VREGION.SHM.AA',
	
	'GAPPED.VREGION.FR1.AA',
	'GAPPED.VREGION.CDR1.AA',
	'GAPPED.VREGION.FR2.AA',
	'GAPPED.VREGION.CDR2.AA',
	'GAPPED.VREGION.FR3.AA',	
	'GAPPED.CDR3.AA',
	'GAPPED.JREGION.FR4.AA',
	
	'GAPPED.VREGION.FR1.NT',
	'GAPPED.VREGION.CDR1.NT',
	'GAPPED.VREGION.FR2.NT',
	'GAPPED.VREGION.CDR2.NT',
	'GAPPED.VREGION.FR3.NT',	
	'GAPPED.CDR3.NT',
	'GAPPED.JREGION.FR4.NT',
	
	'VREGION.VGENE_QUERY_START',
	'VREGION.VGENE_QUERY_END',
	'JREGION.JGENE_QUERY_START',
	'JREGION.JGENE_QUERY_END',
	'NOTES',	
	 expIdentifier,
	'ANALYSIS_NAME',
	'DATE_UPDATED',
	'FILENAME',
	'SETINGS'
]

keys_with_object_id = ['_id','SEQ_ID','EXP_ID']

#current filetypes we account for:
#TAB => TAB DELIMITED FILE. ALWAYS FLATTEN RESULTS FROM DATABASE QUERY. CONVERT ALL RESULTS TO STRING FOR PROPER FILE WRITING
#CSV => COMMA SPERATED FILE. ALWAYS FLATTEN RESULTS FROM DATABASE QUERY. CONVERT ALL RESULTS TO STRING FOR PROPER FILE WRITING
#IGREP => IGREP FILE FORMAT => DATA STORED AS JSON FORMAT, BUT RESULTS/KEYS IN DICTIONARY ARE ALWAYS FLATTENED AND VALUES ARE ALWAYS CONVERTED TO STRING
#JSON => JSON FILE FORMAT. NEVER FLATTEN DICTIONARY OR CONVER RESULTS TO STR (JUST USE JSON.DUMPS)
#FASTA => FASTA FILE FORMAT. DATA ALWAYS FLATTENED. ANY RESULTS THAT ARE NOT DEFINED IN SEQUENCE HEADER ARE REPORTED AS JSON FORMAT >SEQHEADER >{JSONDATA}
#FASTQ => FASTQ FILE FORMAT. DATA ALWAYS FLATTENED. ANY RESULTS THAT ARE NOT DEFINED IN SEQUENCE HEADER ARE REPORTED AS JSON FORMAT @SEQHEADER >{JSONDATA}
allowed_file_types = ['TAB','JSON','IGREP','CSV','FASTA','FASTQ']

#ALL results from proxy were returned to the intermediate file named to_temp_filename
#SEARCH through each line and move sequences to their proper filename defined at the front of the line before file_def_seperator				
#if default_filename is blank, then a coyp of the input file is made and teh default_filename defaults to inputfile name
def split_file_to_multiple_files(inputfile,default_filename = '',ext='query',max_lines_per_file=0):
	dir_path = os.path.dirname(inputfile)
	if not default_filename:		
		tempname = inputfile+'.temp'
		#os.system("cp '{0}' '{1}'".format(inputfile,tempname))
		shutil.copyfile(inputfile,tempname)
		default_filename = os.path.basename(inputfile)
		inputfile = tempname
	time_now = str(datetime.now()).replace(' ','').replace(':','')
	location_of_files_created = "query_files_created_"+time_now+str(random.randint(1, 100))+".txt"
	
	awk_command = '''awk 'BEGIN{{FS="{0}";OFS="";limitseq=({5}>0);num_header=0;check_for_header=1}}						   
					   check_for_header && $1=="HeaderDescription"{{header_lines[num_header]=$2;num_header+=1;next}}					   
					   check_for_header{{check_for_header=0}}
					   {{											
							if (NF>1 && $1!=""){{
								#split($1,files,",");																												
								files=$1;										
								$1="";																														
							}}
							else{{			
								$1=$1;#AWK REQUIRES A FIELD TO BE CHANGED FOR OFS								
								files="{3}"
							}}
							counter[files]+=1
							
						}}							
						limitseq&&counter[files]>{5}{{files=files"."int((counter[files]-1)/{5}); }} #if we are limiting by max_seq_per_file, then run this pattern/line check. this will add an extra suffix to file base don its divisor
						{{							
							files=files".{6}"							
							filenames[files]+=1
							f = "{1}/"files;							
							if (filenames[files]==1){{
								for(i = 0; i< num_header; i++)
									print header_lines[i] > f							
							}}
							print $0 > f														
						}}
						END{{for (files in filenames){{
									print files"\t"filenames[files]>"{1}/{4}"
									}}
							}}' '{2}' '''.format(file_def_seperator,dir_path,inputfile,default_filename,location_of_files_created,max_lines_per_file,ext)			

	os.system(awk_command)
	#os.system("rm '{0}'".format(inputfile))
	os.remove(inputfile)
	
	if os.path.isfile(dir_path+'/'+location_of_files_created):				
		with open(dir_path+'/'+location_of_files_created) as f:		
			files_created = []
			for file in f.readlines():				
				file = file.split('\t')
				if file[0]:
					files_created.append({'filename':dir_path+'/'+file[0],'linecount':int(file[1])})					
		#os.system("rm '{0}'".format(dir_path+'/'+location_of_files_created))
		os.remove(dir_path+'/'+location_of_files_created)
		return files_created
		
	else:		
		return []
		

#chunk size => number of lines to read in the file before returnign results to user
#results will always be returned as a list of dictionaries 
def pass_query_results(filepath,chunk_size=1000):			
	with open(filepath,'r') as reader:
		counter = 0 
		list = []		
		for lines in reader:
			line = lines.strip()
			if not line:
				continue
			line = json.loads(line)
			counter+=1
			list.append(line)
			if counter%chunk_size == 0:
				yield list
				list = []		
	#os.system("rm '{0}'".format(filepath))
	os.remove(filepath)
	if list:				
		yield list		
	if counter == 0:
		#no results 		
		yield []
				
			
def post_to_proxy(query_class=None, path="http://biotseq.ut.appsoma.com:5998",chunk_size=1,filename_suffix='query',max_doc_per_file=0):				
	
	query_to_file = query_class.to_file			
	temp_file_name_path = 'IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')			
			
	if query_to_file:	
		if query_class.file_prefix == None:#user did not pass in path, so we will maek a temp name 
			query_class.file_prefix = temp_file_name_path								
				
		file_prefix = os.path.basename(query_class.file_prefix) if query_class.file_prefix!='' else ''				
		
		if len(query_class.file_prefix.split('/'))==1:			
			#there is no directory defined
			if file_prefix == '':
				file_path_if_empty = 'scratch/'+temp_file_name_path
			else:
				file_path_if_empty = 'scratch/'+file_prefix
		else:
			if file_prefix == '':
				file_path_if_empty = query_class.file_prefix+temp_file_name_path
			else:
				file_path_if_empty = query_class.file_prefix
		
		to_temp_filename = file_path_if_empty+'.intermediate'
		
	else:
		file_prefix = None
		to_temp_filename='IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')+'.temp'		
	
	class_attibutes = {
		'to_file':query_to_file,		
		'modify_query_values':query_class.modify_query_values_to_follow_db_schema,
		'redirect_query_fields':query_class.redirect_fields_for_improved_queries,
		'file_prefix':file_prefix,
	}
			
	print('FIX THE AUTHENTICATION METHODS!!')
	#data to send to http command 
	data = {
		'db_action': 'query',
		'command': query_class.query_command_list, #list of functions to run on proxy
		'authkey':'',#appsoma_api.environment_get_authkey(), #identification of username running in appsoma
		'query_object_id':query_class.name,
		'connection_type':query_class.db_method, #how are we conneting to database (read access or write access, or better yet -> username to mongo connect)		
	}			
	data.update(class_attibutes)
			
	print 'trying here'
	query_string =  immunogrep_proxy_tools.http(
		path,
		data=bson_dumps(data),
		action="POST",
		headers={ "Content-Type":"application/json" },
		toFilename=to_temp_filename# always output data to file ,
		#progressCallback=progress		
	)	
	print 'endinghere'
	
	if query_to_file:
		#just save the results from the file. run split_file_to_multiple_files function so that it will seperate documents in the file use a simple awk command
		#if query_to_file, then return a list of filenames created		
		query_result = split_file_to_multiple_files(to_temp_filename,file_path_if_empty,filename_suffix,max_doc_per_file)		
	else: 		
		#if to memeory then return a generator that can be used to read the file 
		#query_result = json.loads(query_string) #bson_loads(query_string)# json.loads(query_string)		
		query_result = pass_query_results(to_temp_filename,chunk_size)	
	query_class.query_command_list = []
	return query_result


#FOR REFERENCE:
#http://docs.mongodb.org/manual/reference/operator/query/

#currently unused variables##
top_level_operators = ['$and','$or','$nor']	
field_operators = ['$in','$nin']
single_field_value_operators = ['$eq','$ne','$gt','$gte','$lt','$lte']
array_field_value_operators = ['$all']
more_operators = ['$elemMatch']
recursive_field_operators = ['$all']
#####


#do not change the values to any keys in a query in the following list
#these values are standard numbers/booleans based on a mongodb query and not the field it self (for example fields that are strings still will use True/False for $exists)
operators_do_not_modify = ['$size','$exists','$mod','$options','$type','$meta','$language','$geoWithin','$geoIntersects','$near','$minDistance','$maxDistance','$geometry','$nearSphere']
#the following should be top_level_operators that do not 
#1) we do not use a 'text' index in our schema currently, so we are not handling how to parse this query
	#without a 'text index' on a specific field then any query that defines this field will raise an error in mongodb
#2) the '$where' command uses javascript functions to perform a query. again we will not parse these, and leave these special cases to the user
	#not sure when using '$where' would be useful 
#3) comment mongo operator is just for adding comments to a query 
operators_do_not_parse = ['$text','$where','$comment'] 


#now we will go through each of the field defined by the DOT NOTATION, and see if any are NUMBERS
#basically what we are doing is removing any NUMBERS from the field_name. 
#for example the following query: {'VREGION.VGENES.1':'APPLE'} ASKS WHETHER THE SECOND ELEMENT IN THE ARRAY VREGION.VGENES IS APple.
	#but the field name we care about is VREGION.VGENES
#this function is for handling MULTIKEY STRUCTURES IN MONGODB
def Remove_Numbers_From_Field_Name(db_field_name):
	#first split the db_field_name by '.'		
	db_sub_fields = db_field_name.split('.')	
	
	outer_field = []
	reached_integer = False
	db_field_name = []
	
	for f in db_sub_fields:
		try:
			#test if number
			int(f)
			reached_integer=True
		except:
			#if its not a number, then it must be a string so , this is part of the field name 
			if not(reached_integer):
				#keep track of any field before a number is reached
				outer_field.append(str(f))
			#keep track of string where numbers are removed 
			db_field_name.append(str(f))
					
	db_field_name = '.'.join(db_field_name)
	outer_field = '.'.join(outer_field)
	#outer_fields => fields name before the first number was encountered
	#db_field_name = > fieldname wihtout number
	return [db_field_name,outer_field]


###THE FOLLOWING FUNCTIONS/VARIABLES HELP PROCESS QUERY FUNCTIONS TO RETURN ACCURATE RESULTS##
#SUMMARY/FLOW OF CODE: Parse_Mongo_Query_Expression -> uses redirection_fields (redirect_seq_collection_fields,redirect_exp_collection_fields), and calls -> Process_Field_Value -> uses dict_defining_value_transformation (fields_for_queries_seqs_collection, or fields_for_queries_exps_collection)
#first function 
#Parse_Mongo_Query_Expression -> PASS IN A MONGO QUERY DICTIONARY
#	USER DEFINES WHETHER TO redirect_fields_for_improved_queries (for example, instead of query by EXPERIMENT_NAME, query by DUPLICATED_FIELDS.EXPERIMENT_NAME (default = True)
#	USER DEFINES WHETHER TO modify_query_values_to_follow_db_schema (i.e. lower case or uppercase values based on field name) (default = True)
#	->Parse_Mongo_Query_Expression -> will recursively loop through the dictionary and subdictionaries. 
#variables used in Parse function: redirect_seq_collection_fields or redirect_exp_collection_fields
#	->Parse_Mongo_Query_Expression function uses the variable redirection_fields which is either set to redirect_seq_collection_fields or redirect_exp_collection_fields define
#		-> If redirect_fields_for_improved_queries is set to True, then each key in the dictionary/subdictionaries that is found to not be a mongo operatory ('$some-operator') is evaluated
#			The variables redirect_seq_collection_fields and redirect_exp_collection_fields define which query field names should be renamed to other fields 
#		-> If modify_query_values_to_follow_db_schema is set to True
#Next function
#			-> Parse_Mongo_Query_Expression will call the function: Process_Field_Value. Parse_Mongo_Query_Expression will pass in the database field name, the value requested by the user, and either
#			-> the varaible fields_for_queries_seqs_collection or fields_for_queries_exps_collection
#IN Process_Field_Value:
#	-> if the variable db_field_name is found within the dictionary dict_defining_value_transformation, then the variable values is modified based on the function pointed to in the dict_defining_value_transformation
#		function


#DETAILS
#these fields are duplicated in the exp collection 
#the duplicated values to the fields are transformed to imporve queries 
#for example, the word "hello there Sam" becomes "hellotheresam"
metadata_duplicated_fields = schema.Exps_Collection()['duplicated_fields']

#If desired/used, then this dictionary will say what fields should be renamed to before query.
#these variables will be used by the function 'Parse_Mongo_Query_Expression'
#for example, a user may request a query on EXPERIMENT_NAME:value, and if in the function Parse_Mongo_Query_Expression, the variable 'redirect_fields' is True, then the query will 
#be changed to 'DUPLICATED_FIELDS.EXPERIMENT_NAME':value. This is because the values in this field are better made for general queries (i.e. case insensitivity)
redirect_exp_collection_fields = {field_name: lambda x:'DUPLICATED_FIELDS.'+x for field_name in metadata_duplicated_fields}

#now make a similar variable for fields in the seq collection that we specially make for imporved quieres. Currently we only do special treatments fo gene/allele names
redirect_seq_collection_fields = {
	'DATA.VREGION.VGENES':lambda x: 'QUERY_'+x+'.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.'+x+'.PARSED_ALLELES',#for example, VREGION.VGENES.0 => QUERY_DATA.VREGION.VGENES.0.PARSED_ALLELES
	'DATA.DREGION.DGENES':lambda x: 'QUERY_'+x+'.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.'+x+'.PARSED_ALLELES',#for example, VREGION.VGENES.0 => QUERY_DATA.VREGION.VGENES.0.PARSED_ALLELES
	'DATA.JREGION.JGENES':lambda x: 'QUERY_'+x+'.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.'+x+'.PARSED_ALLELES',#for example, VREGION.VGENES.0 => QUERY_DATA.VREGION.VGENES.0.PARSED_ALLELES
}


#this variable will tell the server how to properly transform query requests on the seqs collection sent by user 
#for example, if a user queries by CDR3.NT = 'aactg' , this variable will make the query uppercased -> CDR3.NT='AACTG'
#if a key is not defined in variable below, then by default the value will be treated by the function default_fields_data_types in immunogrep_database_schema
fields_for_queries_seqs_collection = defaultdict(lambda:schema.default_fields_data_types,{		
	#THE FOLLOWING FIELDS ARE AUTOMATICALLY HANDLED AND INSERTED AS DOCUMENT FIELDS INTO THE DATABASE BY THE SERVER 
	#NO LONGER CONVERT DATATYPE TO STR, ASSUME DATATYPE IS A STRING OR USER HAS DECIDED TO QUERY BY DIFFERNET DATATYPE
	"ANALYSIS_NAME":lambda x:x.upper(), #lambda x:str(x).upper(), #uppercase analysisname
	"RECOMBINATION_TYPE":lambda x:x.upper(), #str(x).upper(), #uppercase RECOMBINATION_TYPE
	expIdentifier:convert_to_objectid, #should always be object_id
	'DATE_UPDATED':lambda x:x.upper(), #str(x).upper(), #uppercase DATE. date exists in DB as a STR of day/month/year
	'_id':convert_to_objectid, #should always be object_id
	idIdentifier:convert_to_objectid, #should always be object_id
	#'SETTINGS':lambda x:int(x), #must always be an integer ==> assume user knows to make it int 				
	
	'QUERY_DATA.VREGION.VGENES.PARSED_ALLELES':lambda x: x.upper(),#str(x).upper(), #will always uppercase the query for gene names
	'QUERY_DATA.DREGION.DGENES.PARSED_ALLELES':lambda x: x.upper(),#str(x).upper(), #will always uppercase the query for gene names
	'QUERY_DATA.JREGION.JGENES.PARSED_ALLELES':lambda x: x.upper(),#str(x).upper(), #will always uppercase the query for gene names
	
	
	#THE FOLLOWING FIELDS ARE DEFINED BY RESULTS IN THE INSERT/UPDATE DATABASE FILES (THAT IS WHEN RUNNING THE INSERT/UPDATE FUNCTION)
	#they are pased in by the user
	#!!EACH OF THESE FIELDS WILL GO UNDER 'DATA'. that is, they appear as a subdocument under DATA!!
	###
	#THESE FIELDS ARE ONLY RELEVENT FOR DOCUEMNTS WITH @SEQ ANALYSIS_NAME
	#"DATA.SEQUENCE_HEADER":None, #sequence header of the sequence, datatype = string ==> modified as a 'default value' parameter  => see immunogrep_database_schema.default_fields_data_types() function 
	#"DATA.SEQUENCE":None	#actual raw nucleotide sequence , datatype = string	==> modified as a 'default value' parameter  => see immunogrep_database_schema.default_fields_data_types() function 
	#"DATA.QUALITY_SCORE":None, #quality score of the sequence ==> modified as a 'default value' parameter => see immunogrep_database_schema.default_fields_data_types() function 
	'DATA.FILENAME':lambda x: x, #str(x),#do not change casing 
	#
	
	#THESE FIELDS ARE RELEVENT FOR DOCUMENTS IN ANY OTHER ANALYSIS_NAME
	"DATA.COMMAND": lambda x:x.upper(), #str(x).upper(),	
	"DATA.NOTES":lambda x: x.upper(),#str(x).upper(),	
	"DATA.PREDICTED_AB_SEQ.NT":lambda x: x.upper(),#str(x).upper(),	
	"DATA.PREDICTED_AB_SEQ.AA":lambda x: x.upper(),#str(x).upper(),	
	"DATA.STRAND":lambda x:x.upper(),# str(x).upper(),
	"DATA.PREDICTED_CHAIN_TYPE":lambda x:x.upper(),# str(x).upper(),
	"DATA.PRODUCTIVE":lambda x: x.upper(),#str(x).upper(),
	"DATA.LOCUS":lambda x: x.upper(),#str(x).upper(),
	"DATA.VREGION.SHM.NT":lambda x:float(x),
	"DATA.VREGION.SHM.AA":lambda x:float(x),	
	"DATA.VREGION.SHM.NT_PER":lambda x:float(x),
	"DATA.VREGION.SHM.AA_PER":lambda x:float(x),	
	"DATA.JREGION.SHM.NT":lambda x:float(x),
	"DATA.JREGION.SHM.AA":lambda x:float(x),	
	"DATA.JREGION.SHM.NT_PER":lambda x:float(x),
	"DATA.JREGION.SHM.AA_PER":lambda x:float(x),	
	'DATA.VREGION.VGENE_QUERY_START':lambda x:int(x),
	'DATA.VREGION.VGENE_QUERY_END':lambda x:int(x),
	"DATA.VREGION.FR1.NT":lambda x:x.upper(), #str(x).upper(),
	"DATA.VREGION.FR1.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.CDR1.NT":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.CDR1.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.FR2.NT":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.FR2.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.CDR2.NT":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.CDR2.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.FR3.NT":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.FR3.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.VREGION.VGENES":lambda x:x.upper(),#str(x),#.upper(),
	"DATA.VREGION.VGENE_SCORES": lambda x:float(x),	
	"DATA.CDR3.NT":lambda x: x.upper(),# str(x).upper(),
	"DATA.CDR3.AA":lambda x: x.upper(),#:str(x).upper(),
	"DATA.DREGION.DGENES":lambda x:x.upper(),#str(x),#.upper(),
	"DATA.DREGION.DGENE_SCORES": lambda x:float(x),	
	"DATA.JREGION.FR4.NT":lambda x:x.upper(),#str(x).upper(),
	"DATA.JREGION.FR4.AA":lambda x:x.upper(),#str(x).upper(),
	"DATA.JREGION.JGENES":lambda x:x.upper(),#str(x),#.upper(),
	"DATA.JREGION.JGENE_SCORES":lambda x:float(x),	
	'DATA.JREGION.JGENE_QUERY_START':lambda x:int(x),
	'DATA.JREGION.JGENE_QUERY_END':lambda x:int(x),
	'DATA.CDR3.AA_LENGTH':lambda x:int(x),
	'DATA.CDR3.NT_LENGTH':lambda x:int(x),
	'DATA.VREGION.CDR1.AA_LENGTH':lambda x:int(x),
	'DATA.VREGION.CDR1.NT_LENGTH':lambda x:int(x),
	'DATA.VREGION.CDR2.AA_LENGTH':lambda x:int(x),
	'DATA.VREGION.CDR2.NT_LENGTH':lambda x:int(x),	
	'DATA.ISOTYPE.GENE':lambda x:x.upper(),
	'DATA.ISOTYPE.MISMATCHES':lambda x:int(x),
	'DATA.ISOTYPE.SCORES':lambda x:int(x),
	'DATA.ISOTYPE.PER_ID':lambda x:float(x)
})

#this variable will tell the server how to properly transform query requests on the EXPS collection sent by user 
#for example, if a user queries by DUPLICATED_FIELDS.EXPERIMENT_NAME = 'howdy there' , this variable will lowercase the query, remove empty spaces, and change it to => DUPLICATED_FIELDS.EXPERIMENT_NAME='howdythere'
#if key does not exist, then default function/treatment of the value/query will be defined by the function default_metadata_fields in schema file 
fields_for_queries_exps_collection = defaultdict(lambda:schema.default_metadata_fields,
	dict(
		{		
			#THE FOLLOWING FIELDS ARE AUTOMATICALLY HANDLED AND INSERTED AS DOCUMENT FIELDS INTO THE DATABASE BY THE SERVER 				
			'_id':convert_to_objectid, #should always be object_id	
			'PROJECT_NAME':lambda x:str(x),
			'EXPERIMENT_NAME':lambda x:str(x),
			'SEQUENCING_PLATFORM':lambda x:x.lower(),# str(x).lower(),
			'SPECIES': lambda x:' '.join([word.title() if i==0 else word.lower() for i,word in enumerate(x.split(' '))]), #REMOVED STR(X).SPLIT because, if it doesnt work, will revert to orogina setting by user. #capitalize first word only 
			'LAB': lambda x:x.lower(),#str(x).lower(),#lowercase
			'OWNERS_OF_EXPERIMENT': lambda x:x.lower(),#str(x).lower(),#lowercase
			'READ_ACCESS': lambda x: x.lower(),# str(x).lower(),#lowercase
			'CHAIN_TYPES_SEQUENCED': lambda x: str(x).lower(),#lowercase
			'CELL_TYPES_SEQUENCED': lambda x: str(x), #leave alone 
			'ISOTYPES_SEQUENCED': lambda x: str(x), #leave alone 
			'DESCRIPTION': lambda x: str(x), #leave alone
			'MID_TAG': lambda x: x.replace(' ','').lower(),#lowercase adn remove spaces
			'CELL_MARKERS_USED': lambda x:x.lower(), #str(x).lower(),#lowercase
			'PAIRING_TECHNIQUE':lambda x: x, #str(x),#leave alone
			'PUBLICATIONS': lambda x: x, #str(x),#leave alone
			'CELL_NUMBER':lambda x:int(x),#integer
			'TARGET_READS':lambda x:int(x),#integer
			'CONTAINS_RNA_SEQ_DATA': schema.AttemptToConvertToBool, #for boolean operations 
			'VH:VL_PAIRED': schema.AttemptToConvertToBool, #for boolean operations 
			'LIST_OF_POLYMERASES_USED':lambda x: str(x).lower(),#lowercase
			'PRIMER_SET_NAME': lambda x:x.lower(),# str(x).lower(),#lowercase
			"POST_SEQUENCING_PROCESSING:PHI_X_FILTER":schema.AttemptToConvertToBool,
			"POST_SEQUENCING_PROCESSING:QUALITY_FILTER":lambda x:x, #str(x),
			"POST_SEQUENCING_PROCESSING:PROCESS_R1_R2_FILE": lambda x:x,# str(x),
			"PERSON_WHO_PREPARED_LIBRARY":lambda x: str(x).title(),#capitalize every word 
			"REVERSE_PRIMER_USED_IN_RT_STEP":lambda(x):str(x).lower() #lowercase 
		}.items()+
		{'DUPLICATED_FIELDS.'+field_name:convert_text_to_index_field_text for field_name in metadata_duplicated_fields}.items()#now any of the fields taht we duplicate for indexing/better queries, we add them here to the variable	
	)
)


def default_fields_to_file(x,splitter=','):
	if type(x) is list:
		return (splitter).join([str(e) for e in x])
	else:
		return str(x)	
		
def default_fields_to_file_no_commas(x):	
	if type(x) is list:
		return ('|').join([str(e) for e in x])
	else:
		return str(x)		

#simply convert to str
def ReturnStr(x):
	return str(x)

#convert list into string. seperate elements by d char
def ReturnListStr(x,d):
	return (d).join(x)

#convert list containing NON STRINGS into string. seperate elements by d char
def ReturnListNotStr(x,d):
	v = str(x[0])
	for j in range(1,len(x)):
		v+=d+str(x[j])
	return v

#schema_fields_to_file => for each key, transforms the value from the database to the defined after getting a query result 
#if schema_fields_to_file =
	#0 => run default convert to str function
	#1 => convert to str
	#2 => convert from list, delim by commas 
	#3 => DO nothing, its already a 
schema_fields_to_file = defaultdict(lambda: default_fields_to_file, {
	"COMMAND":lambda x: str(x),				
	"NOTES":lambda x: str(x),	
	"PREDICTED_AB_SEQ.NT":lambda x: str(x),
	"PREDICTED_AB_SEQ.AA":lambda x: str(x),
	"STRAND":lambda x: str(x),	
	"PREDICTED_CHAIN_TYPE":lambda x: str(x),
	"PRODUCTIVE":lambda x: str(x),
	"LOCUS":lambda x: str(x),
	"VREGION.SHM.NT":lambda x: str(x),
	"VREGION.SHM.AA":lambda x: str(x),	
	"VREGION.SHM.NT_PER":lambda x: str(x),
	"VREGION.SHM.AA_PER":lambda x: str(x),	
	"JREGION.SHM.NT":lambda x: str(x),
	"JREGION.SHM.AA":lambda x: str(x),	
	"JREGION.SHM.NT_PER":lambda x: str(x),
	"JREGION.SHM.AA_PER":lambda x: str(x),	
	'VREGION.VGENE_QUERY_START':lambda x: str(x),	
	'VREGION.VGENE_QUERY_END':lambda x: str(x),
	"VREGION.FR1.NT":lambda x: str(x),
	"VREGION.FR1.AA":lambda x: str(x),	
	"VREGION.CDR1.NT":lambda x: str(x),	
	"VREGION.CDR1.AA":lambda x: str(x),	
	"VREGION.FR2.NT":lambda x: str(x),
	"VREGION.FR2.AA":lambda x: str(x),	
	"VREGION.CDR2.NT":lambda x: str(x),	
	"VREGION.CDR2.AA":lambda x: str(x),	
	"VREGION.FR3.NT":lambda x: str(x),	
	"VREGION.FR3.AA":lambda x: str(x),	
	"VREGION.VGENES":lambda x: ReturnListStr(x,','),	
	"VREGION.VGENE_SCORES": lambda x: ReturnListNotStr(x,','),
	"CDR3.NT":lambda x: str(x),
	"CDR3.AA":lambda x: str(x),
	"DREGION.DGENES":lambda x: ReturnListStr(x,','),
	"DREGION.DGENE_SCORES":lambda x: ReturnListNotStr(x,','),	
	"JREGION.FR4.NT":lambda x: str(x),	
	"JREGION.FR4.AA":lambda x: str(x),	
	"JREGION.JGENES":lambda x: ReturnListStr(x,','),	
	"JREGION.JGENE_SCORES":lambda x: ReturnListNotStr(x,','),	
	'JREGION.JGENE_QUERY_START':lambda x: str(x),	
	'JREGION.JGENE_QUERY_END':lambda x: str(x),	
	'CDR3.AA_LENGTH':lambda x: str(x),	
	'CDR3.NT_LENGTH':lambda x: str(x),	
	'VREGION.CDR1.AA_LENGTH':lambda x: str(x),	
	'VREGION.CDR1.NT_LENGTH':lambda x: str(x),	
	'VREGION.CDR2.AA_LENGTH':lambda x: str(x),	
	'VREGION.CDR2.NT_LENGTH':lambda x: str(x),	
	'FILENAME':lambda x: str(x),
	'GAPPED.VREGION.FR1.NT':lambda x: str(x),
	'GAPPED.VREGION.FR2.NT':lambda x: str(x),
	'GAPPED.VREGION.FR3.NT':lambda x: str(x),
	'GAPPED.VREGION.CDR1.NT':lambda x: str(x),
	'GAPPED.VREGION.CDR2.NT':lambda x: str(x),	
	'GAPPED.VREGION.FR1.AA':lambda x: str(x),
	'GAPPED.VREGION.FR2.AA':lambda x: str(x),
	'GAPPED.VREGION.FR3.AA':lambda x: str(x),
	'GAPPED.VREGION.CDR1.AA':lambda x: str(x),
	'GAPPED.VREGION.CDR2.AA':lambda x: str(x),
	'GAPPED.JREGION.FR4.NT':lambda x: str(x),
	'GAPPED.JREGION.FR4.AA':lambda x: str(x),
	'GAPPED.CDR3.NT':lambda x: str(x),
	'GAPPED.CDR3.AA':lambda x: str(x),
	'SETTINGS':lambda x: str(x),
	'ANALYSIS_NAME':lambda x: str(x),
	'RECOMBINATION_TYPE':lambda x: str(x),
	'SEQ_ID':lambda x: str(x),
	'DATE_UPDATED':lambda x: str(x),
	'EXP_ID':lambda x: str(x),
	'_id':lambda x: str(x),
	'ISOTYPE.GENE':lambda x: ReturnListStr(x,','),
	'ISOTYPE.MISMATCHES':lambda x: ReturnListNotStr(x,','),
	'ISOTYPE.SCORES':lambda x: ReturnListNotStr(x,','),
	'ISOTYPE.PER_ID':lambda x: ReturnListNotStr(x,',')
})



#schema_fields_to_file => for each key, transforms the value from the database to the defined after getting a query result 
#if schema_fields_to_file =
	#0 => run default convert to str function
	#1 => convert to str
	#2 => convert from list, delim by commas 
	#3 => DO nothing, its already a 
schema_fields_to_file_avoid_commas = defaultdict(lambda: default_fields_to_file_no_commas)
for s_keys,s_values in schema_fields_to_file.iteritems():
	schema_fields_to_file_avoid_commas[s_keys] = s_values
schema_fields_to_file_avoid_commas["VREGION.VGENES"] =  lambda x: ReturnListStr(x,'|')
schema_fields_to_file_avoid_commas["VREGION.VGENE_SCORES"] =  lambda x: ReturnListNotStr(x,'|')
schema_fields_to_file_avoid_commas["JREGION.JGENES"] =  lambda x: ReturnListStr(x,'|')	
schema_fields_to_file_avoid_commas["JREGION.JGENE_SCORES"] =  lambda x: ReturnListNotStr(x,'|')
schema_fields_to_file_avoid_commas["DREGION.DGENES"] =  lambda x: ReturnListStr(x,'|')	
schema_fields_to_file_avoid_commas["DREGION.DGENE_SCORES"] =  lambda x: ReturnListNotStr(x,'|')
schema_fields_to_file_avoid_commas["ISOTYPE.GENE"] =  lambda x: ReturnListStr(x,'|')	
schema_fields_to_file_avoid_commas["ISOTYPE.MISMATCHES"] =  lambda x: ReturnListNotStr(x,'|')
schema_fields_to_file_avoid_commas["ISOTYPE.SCORES"] =  lambda x: ReturnListNotStr(x,'|')
schema_fields_to_file_avoid_commas["ISOTYPE.PER_ID"] =  lambda x: ReturnListNotStr(x,'|')




#USING either the fields_for_queries_seqs_collection or the fields_for_queries_exps_collection variables (that is, dict_defining_value_transformation is set to one of these variables depending on query function)
#this program will properly convert the passed in values for a specicific db_field_name to match the format in the database
#db_field_name will correspond to fields in the database using '.' notation 
#the values, will the be request query values for that field name 
#the values shoudl only be of type list or other individual elemsnts such as str, int, float. they shoudlnt be dict for example 
def Process_Field_Value(db_field_name,values,dict_defining_value_transformation):	
	
	[db_field_name,outer_field] = Remove_Numbers_From_Field_Name(db_field_name)	
	if type(values) is list:
		temp = []		
		for list_element in values:
			try: #use a try/except, because we assume that teh user may have a better knowledge of the query. so if it should be an int, but now its a string, then default to what user passed in. worst case scenario, the query will return no results
				if isinstance(list_element, type(re.compile(''))) or isinstance(list_element, type(bson_loads(bson_dumps(re.compile(''))))) or isinstance(list_element, type(bson_loads(bson_dumps({'$regex':''})))):#always ignore types that are REGULAR EXPRESSIONS. Cannot compile those to match datatype  										
					#this means that the value is a regular expression. 
					#in order to convert a re format from python object to a dictionary, do the following: 
					re_dict_format = json.loads(bson_dumps(list_element)) # first dump the object to a string using bson, then load it as a dic using json
					#now we can actually modify the string that has been made into a regular expression (it is found in '$regex key'
					re_dict_format['$regex'] = dict_defining_value_transformation[db_field_name](re_dict_format['$regex'])
					#now that we have formated the variable correctly, lets update v, but convert  this dict back into python re format
					v = bson_loads(bson_dumps(re_dict_format))
				else:
					v = dict_defining_value_transformation[db_field_name](list_element)
			except Exception as e:
				print "Error in compiling fields values: {0}. Error {1}".format(str(values),str(e))
				v = list_element
			temp.append(v)		
	else:
		try:#use a try/except, because we assume that teh user may have a better knowledge of the query. so if it should be an int, but now its a string, then default to what user passed in. worst case scenario, the query will return no results
			if isinstance(values, type(re.compile(''))) or isinstance(values, type(bson_loads(bson_dumps(re.compile(''))))) or isinstance(values, type(bson_loads(bson_dumps({'$regex':''})))):#always ignore types that are REGULAR EXPRESSIONS. Cannot compile those to match datatype  									
					#this means that the value is a regular expression. 
					#in order to convert a re format from python object to a dictionary, do the following: 
					re_dict_format = json.loads(bson_dumps(values)) # first dump the object to a string using bson, then load it as a dic using json
					#now we can actually modify the string that has been made into a regular expression (it is found in '$regex key'
					re_dict_format['$regex'] = dict_defining_value_transformation[db_field_name](re_dict_format['$regex'])
					
					#now that we have formated the variable correctly, lets update v, but convert  this dict back into python re format
					temp = bson_loads(bson_dumps(re_dict_format))					
			else:
				temp = dict_defining_value_transformation[db_field_name](values)
		except Exception as e:
			print "Error in compiling fields values: {0}. Error {1}".format(str(values),str(e))
			temp = values		
	return temp 

def Modify_Key_Name_To_Include_Data(key):
	if key not in fields_above_data_key:
		first_field = key.split('.')[0]
		if first_field!=query_data_key_name and first_field!='DATA':
			key = 'DATA.'+key
	return key
	
	
#when performing a query, this will attempt to ensure that the values for each field name are properly formated
#for example. if we always uppercase the values in a field, then make sure that this query also has only uppercased letters
#also, if we have a better field for query, such as DUPLICATED_FIELDS in the exp_collection or QUERY_DATA in seqs collection, then again
#redirect the field name to include this field 
#redirect_fields_for_improved_queries => this variable defines whether queiry fields shoudl be renamed,
#for example, if a user wants to query by 'DATA.VREGION.VGENES', it woudl be much for flexible if the query was on the field 'QUERY_DATA.VREGION.VGENES.PARSED_ALLELES' => this field int he database allows user to query by allele/genename in any case-senistiive
#so when redirect_fields_for_improved_queries = True, select queries on fields such as 'DATA.VREGION.VGENES' will be renamed to the field designed for improved queires
#modify_query_values_to_follow_db_schmea => this variable defines whther queries on specific fields should be modified to better match its value in the database. for example, most string fields in SEQS collection are capitalized. 
#so when this variable is true, the function will look for any fields in the query that should be modified (i.e. whose values should be captialized)
#when both modify_query_values_to_follow_db_schema=False AND redirect_fields_for_improved_queries=False, then NO MODIFICATION OF THE query occurs 
def Parse_Mongo_Query_Expression(query,field_name='',dict_defining_value_transformation={},redirection_fields={},modify_query_values_to_follow_db_schema=True, redirect_fields_for_improved_queries=True):
	
	def RunParser(query,field_name,dict_defining_value_transformation,redirection_fields,modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries):
		ignore_mongo_operators = operators_do_not_modify+operators_do_not_parse	
		if type(query) is dict: 
			copied_query = {f:v for f,v in query.iteritems()}# in copy.deepcopy(query) => copy.deepcopy does not work on re.compile 
			#looop through all keys in dictionary 
			for key,values in copied_query.iteritems():
				#rename keys that should be redirected
				#this will only work if the field is explicitly written in dot notation. For example, it will not work if its CDR3:{'AA':'APPL'}. we cannot redirect this field to  QUERY_DATA.CDR3.AA
				if redirect_fields_for_improved_queries and key[0]!='$':
					
					#key = Modify_Key_Name_To_Include_Data(key) ===> was attempting to add 'DATA' key to queries, but this would make confusing problems with exp collection vs seqs collection
					
					[key_no_int,outer_key] = Remove_Numbers_From_Field_Name(key)
					#this key should be redirected 
					
					if key_no_int in redirection_fields:
						
						popped_value = query.pop(key)
						key = redirection_fields[key_no_int](key)
						query[key] = popped_value
			
				
				if key in ignore_mongo_operators: #operators_do_not_modify:
					#do not change the values passed in for this field
					#i.e. {'field_name_a':{'$exists':True}} => even if field name is a string, obviously the value to exists is only True/False, so do not modify the value of this field			
					query[key] = values #just reutrn the original query 			
					continue
				elif key[0] == '$':#key is a specific mongo operator			
					sub_field_name = field_name #sub_field_name stays the same, because a specific mongo operator does not add to a field name, just the mongo function
				else:
					#its not a mongo operator, so it must be a field name, add this as a subdocument in the query command using '.' notation 
					#for exampple if we are querying: {CDR3:{AA:'ALPHA','NT':'BETA'}} then we will want to LOOK AT THE FOLLOWING FIELD NAMES: CDR3.AA, CDR3.NT
					sub_field_name = field_name+'.'+key if field_name != '' else key
				
				
				if type(values) is list:
					for i,each_value in enumerate(values):
						#if its a list, then recursively parse through each value of list 
						query[key][i] = RunParser(each_value,sub_field_name,dict_defining_value_transformation,redirection_fields,modify_query_values_to_follow_db_schema,redirect_fields_for_improved_queries)								
				elif type(values) is dict:
					#now that we always load dictionary using bson in above/parent fucntion, then the following if/elif statements should
					#not occur. there should no longer be keys that say {'$oid'] or {'$regex'}. instead, they exist as ObjectId and re.copmile 
					#hardcoding object id conversions...
					#no need to check recursifely this dictionary. it is an object id 
					if values.keys()==['$oid']:						
						query[key] = Process_Field_Value(sub_field_name,values,dict_defining_value_transformation)																				
					elif '$regex' in values:
						#hardocindg regular expressions
						#again the user requested using regular expressions , so no need to recursively check this dictionary.
						#first process the regular epxression using process_field_value (only process '$regex' key
						#once complete, use bson_loads to convert it into a python re variable 
						values['$regex'] = Process_Field_Value(sub_field_name,values['$regex'],dict_defining_value_transformation)
						#first uson json.dump to dump dictionary, then RELOAD IT using bson
						query[key] = bson_loads(json.dumps(values))														
					else:				
						#if its a dict, then just recursively pass in the dictionary rather than each index  in a list as above 			
						query[key] = RunParser(values,sub_field_name,dict_defining_value_transformation,redirection_fields,modify_query_values_to_follow_db_schema,redirect_fields_for_improved_queries)								
				else: #assume anythign else is a specific query parameter/value. therefore, prcoess its field value 				
					if modify_query_values_to_follow_db_schema:#modify the value for this query 
						query[key] = Process_Field_Value(sub_field_name,values,dict_defining_value_transformation)			
					
		else: #if the query was not actually a dictionary then it is probably from values in a list, or just singular values, so just process each field name 
			values = query				
			if modify_query_values_to_follow_db_schema:#modify teh value for this query 
				query = Process_Field_Value(field_name,values,dict_defining_value_transformation)
	
		return query
	#first take the provided query and ENSURE that it will be in BSON format:
	#dump query as a string using bson dumps => handles ObjectIds, and regular expression values 
	#once dumped, reload is as a bson to ensure that any previous dictionaries which should be objects (such as object id and regex) are represtend properly
	query = bson_loads(bson_dumps(query))
	
	#now that we have an ensured BSON format dictoinary, run the parsing function if user wants to modify field names or values
	if modify_query_values_to_follow_db_schema or redirect_fields_for_improved_queries:
		query = RunParser(query,field_name,dict_defining_value_transformation,redirection_fields,modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries)
	
	return query
###END OF Parse_Mongo_Query_Expression FUNCTIONS/DOCUMUENTATION##

#THIS FUNCTION WILL ONLY BE RELEVANT FOR QUERIES ON '_id' in exps collection and 'EXP_ID' in seqs collection#
#this is a simple function for getting the exp_id intersection of possible ojbect IDS
#1) first, it takes the intersection of all values in list_of_allowed_ids
#2) then once this id list has been filtered down, it will analyze the id_query 
#3) The id_query will allow the following values:
	#Single Value expression that can be made into an ObjectId
	#A dictionary containing the following keys only: '$eq','$in','$ne', '$nin'
def GetIdIntersection(list_of_allowed_ids, id_query=None):	
	#as an input, list_of_allowed_ids can be either 1) just an objectid, 2) a list of objectids, or 3)list of lists of objectids, in this situation all elements in 'top' list, must be lists.	
	#but, in the end, we want list_of_allowed_ids to b lists of lists
	
	#single element passe din 
	if not isinstance(list_of_allowed_ids,list):
		#convert to lists of lists...
		list_of_allowed_ids =[[convert_to_objectid(list_of_allowed_ids)]]
	else:	#list passed in 	
		#make sure if the user wants a LIST OF LISTS, then every elment in TOP list is LIST 
		#do not allow a 'MIXTURE' of some lists, some not lists 
		are_elements_lists = [isinstance(input,list) for input in list_of_allowed_ids]
	
		if sum(are_elements_lists) == 0:
			#no element wa a list, so convert it into a list of lists
			list_of_allowed_ids = [list_of_allowed_ids]
		elif sum(are_elements_lists) != len(list_of_allowed_ids):
			#at least one elmeent WAS A LIST, but NOT ALL elements 
			raise Exception("The list variable 'list_of_allowed_ids' is improperly formatted. You may only pass in 1) a single value, 2) a list of possible values, or 3) a list of LISTS of possitble values. You passed in the following "+str(list_of_allowed_ids))		
			
	#obtain intersection of all ids defined in list_of_allowed_ids
	combined_ids = set([convert_to_objectid(id) for id in list_of_allowed_ids[0]])
	
	for i in range(1,len(list_of_allowed_ids)):
		combined_ids = combined_ids&set([convert_to_objectid(id) for id in list_of_allowed_ids[i]])	
	
	if id_query:				
		if isinstance(id_query,list):
			combined_ids = combined_ids&set([convert_to_objectid(exp) for exp in id_query])
		elif isinstance(id_query,dict):
			#first we need to ensure that the id_query has been properly loaded into ObjectId's using bson_loads.
			#to ensure this, we will first use dumps on the dictionary, and then reload using loads
			#so if an object id is {'oid:""} , this function will make it into a proper object => ObjectId("")
			id_query = bson_loads(bson_dumps(id_query))
			
			operators_supported = ['$in','$eq','$ne','$nin']
			for each_key,query_val in id_query.iteritems():
				#perform set intersection
				if each_key == '$in':
					if not isinstance(query_val,list):
						query_val = [query_val]
					combined_ids = combined_ids&set([convert_to_objectid(exp) for exp in query_val])
				elif each_key == '$eq':
					#for '_id' and 'EXP_ID' the value can never be EQUAL to a list 
					if isinstance(query_val,list):
						raise Exception("The ObjectId can never be equal ('$eq') to a list. If you want to allow a range of posible object id's, use '$in' operator\nPassed in query: {0}".format(bson_dumps(id_query)))
					combined_ids = combined_ids&set([convert_to_objectid(query_val)])
				#perform set difference
				elif each_key == '$ne':					
					#for '_id' and 'EXP_ID' the value can never be EQUAL to a list 
					if isinstance(query_val,list):
						raise Exception("The ObjectId can never be not-equal ('$ne') to a list. If you want to allow a range of posible object id's, use '$in' operator\nPassed in query: {0}".format(bson_dumps(id_query)))
					combined_ids = combined_ids-set([convert_to_objectid(query_val)])
				elif each_key == '$nin':		
					if not isinstance(query_val,list):
						query_val = [query_val]
					combined_ids = combined_ids-set([convert_to_objectid(exp) for exp in query_val])
				else:
					raise Exception("You have requested the following query on experiment id: {0}. Currently we only support the following mongo operators on the exp id: {1}".format(bson_dumps(id_query),','.join(operators_supported)))
		else:			
			combined_ids = combined_ids&set([convert_to_objectid(id_query)])
		
	combined_ids = list(combined_ids)
	
	#return a proprely formatted mongo query 
	return combined_ids[0] if len(combined_ids)==1 else {'$in':combined_ids}

#when using simple PROJECT (that is when nto using aggregation framework), this will parse the project query and predict which fields will be output from query 
def Get_Schema_Details(possible_metadata_fields,projected_fields,allTrue):
	analyses_types = []
	recombination_types = []
	
	overall_schema = {f:1 for f in fields_above_data_key}
	
	for each_metadata in possible_metadata_fields:
		if 'ANALYSIS_SCHEMA' in each_metadata:			
			schema_keys = flatten_dictionary(each_metadata['ANALYSIS_SCHEMA']).keys()
			for field in schema_keys:
				f = 'DATA.'+'.'.join(field.split('.')[:-1])
				overall_schema[f] = 1
		
		if 'ANALYSES_COUNT' in each_metadata:
			for each_analysis in each_metadata['ANALYSES_COUNT']:
				analyses_types.append(each_analysis)
				recombination_types.extend([recomb for recomb,count in each_metadata['ANALYSES_COUNT'][each_analysis].iteritems() if count>0])	#only choose recombination types whose count > 0
	
	projected_fields = copy.deepcopy(projected_fields)
	if projected_fields == None:
		projected_fields = {}
	not_integer_projections = {}
	for projections in projected_fields.keys():
		values = projected_fields[projections]
		if isinstance(values,dict):
			projected_fields.pop(projections)
			not_integer_projections[projections] = values
	
	if projected_fields == {}:
		for i in overall_schema:
			overall_schema[i] = 1
	elif allTrue:
		for i in overall_schema:
			overall_schema[i] = 0			
	else:
		for i in overall_schema:
			overall_schema[i] = 1
	
			
	
	for projections, values in projected_fields.iteritems():
		if values == 0:
			for fields in overall_schema:
				if fields == projections:
					overall_schema[fields] = 0
				elif fields.startswith(projections+'.'):
					overall_schema[fields] = 0
		else:
			for fields in overall_schema:
				if fields == projections:
					overall_schema[fields] = 1
				elif fields.startswith(projections+'.'):
					overall_schema[fields] = 1
	
	#any field in the projection that was not a number was probably  acomplex projection operator (such as slice operator)
	#so for all of these fields, change their schema values to 1 
	for projections in not_integer_projections:
		for fields in overall_schema:
			if fields == projections:
				overall_schema[fields] = 1
			elif fields.startswith(projections+'.'):
				overall_schema[fields] = 1
			
	
	#remove fields that have the DATA prefix
	possible_fields = [field[5:] if field.startswith('DATA.') else field for field,value in overall_schema.iteritems() if value == 1]
	
	
				
	return [sorted(list(set(possible_fields))),sorted(list(set(analyses_types))),sorted(list(set(recombination_types)))]


def get_new_field_name(new_field_name,original_field,overall_schema):
	new_field = new_field_name.split('.$')[0]
	if original_field[0] =='$':					
		found_a_value=False
		if original_field.startswith('$$CURRENT.'):
			old_field = original_field[10:]
		else: #id field has to start with $ mongo operator
			old_field = original_field[1:] 													
		for each_key in overall_schema.keys():
			#so we found all values that contain the old key name 
			if each_key == old_field:							
				found_a_value=True
				overall_schema[new_field] = 1 #new field name 
			elif each_key.startswith(old_field+'.'):
				found_a_value=True
				new_key = new_field+'.'+( old_field+'.').join(each_key.split(old_field+'.')[1:])							
				overall_schema[new_key] = 1		
		if found_a_value==False:
			overall_schema[new_field] = 1
	else:
		overall_schema[new_field] = 1
	
	
#projects from aggreagtion will be slightly more complicated than simple projects from queries. 
#fields can be renamed, or projected, or new fields can be added. 
#we need to keep track of these fields so that we can adequately predict the field names for outputing to TAB/CSV files
def Get_Schema_Details_Aggregation_Pipeline(possible_metadata_fields,aggregation_pipeline_query):
	
	possible_fields = fields_above_data_key
	overall_schema = {}
	for i in possible_fields:
		overall_schema[i] = 1
	analyses_types = []
	recombination_types = []
	for each_metadata in possible_metadata_fields:
		if 'ANALYSIS_SCHEMA' in each_metadata:
			schema = flatten_dictionary(each_metadata['ANALYSIS_SCHEMA']).keys()
			for i in range(len(schema)):
				schema[i] = '.'.join(schema[i].split('.')[:-1])			
				overall_schema['DATA.'+schema[i]] = 1			
			
		if 'ANALYSES_COUNT' in each_metadata:
			for each_analysis in each_metadata['ANALYSES_COUNT']:
				analyses_types.append(each_analysis)
				recombination_types.extend([recomb for recomb,count in each_metadata['ANALYSES_COUNT'][each_analysis].iteritems() if count>0])	#only choose recombination types whose count > 0
	
	
	for pipe in aggregation_pipeline_query:	
		if pipe.keys()[0]=='$project':
			sub_project_query = flatten_dictionary(pipe['$project'])
			#turn off all key fields (by default nothing is projected, unless explicitiley defined)
			for each_key in overall_schema:
				overall_schema[each_key] = 0						
			#by defaault, '_id' is always shown
			overall_schema['_id'] = 1						
			for field,operator in sub_project_query.iteritems():				
				if '.$literal.' in field:
					overall_schema[field.replace('.$literal.','.')] = 1					
				elif operator == 0:
					for each_key in overall_schema:
						if each_key==field or each_key.startswith(field+'.'):
							overall_schema[each_key] = 0
				elif operator == 1:
					#include these fields from overall schema 
					for each_key in overall_schema:
						if each_key==field or each_key.startswith(field+'.'):
							overall_schema[each_key] = 1
				elif isinstance(operator,basestring) and operator[0] == '$':
					get_new_field_name(field,operator,overall_schema)																									
				else:
					#it is most likely using a mongo operator to project a NEW FIELD, so add the NEW field to overall_schema
					#ONly consider the field name up until MONGO operator '$'
					new_field = field.split('.$')[0]
					overall_schema[new_field] = 1																			
		elif pipe.keys()[0] == '$group':
			#turn off all key fields (by default nothing is projected, unless explicitiley defined)
			for each_key in overall_schema:
				overall_schema[each_key] = 0			
				
			operation_keys = copy.deepcopy(pipe['$group'])			
			
			#traverse the ID operator. ID determiens how to group sequences together
			id_field = operation_keys.pop('_id',None)			
			if not id_field: #user sayd group:{_id:null}
				overall_schema['_id']=1 #only an id field will be made			
			elif isinstance(id_field,basestring): #user is only grouping by specific field. this field will be called 'ID'
				overall_schema['_id']=1 #id must appear in data now 																												
			else: #id is a complex dictionary group values by multiple fields 
				overall_schema['_id'] = 0 
				id_field = flatten_dictionary(id_field)
				for id_keys,id_values in id_field.iteritems():
					if isinstance(id_values,basestring):
						get_new_field_name('_id.'+id_keys,id_values,overall_schema)																	
					else:
						#no idea, so must rename the fields based on the new field name 
						new_field = id_keys.split('.$')[0]
						overall_schema['_id.'+new_field] = 1
								
			operation_keys = flatten_dictionary(operation_keys)
			for keys,values in operation_keys.iteritems():
				if isinstance(values,basestring):
					get_new_field_name(keys,values,overall_schema)																	
				else:
					#no idea, so must rename the fields based on the new field name 
					new_field = keys.split('.$')[0]
					overall_schema[new_field] = 1				
				
	#remove fields that have the DATA prefix
	possible_fields = [field[5:] if field.startswith('DATA.') else field for field,value in overall_schema.iteritems() if value == 1]
			
	return [sorted(list(set(possible_fields))),sorted(list(set(analyses_types))),sorted(list(set(recombination_types)))]
	


#redirect the field name to include this field 
	#redirect_fields_for_improved_queries => this variable defines whether queiry fields shoudl be renamed,
	#for example, if a user wants to query by 'DATA.VREGION.VGENES', it woudl be much for flexible if the query was on the field 'QUERY_DATA.VREGION.VGENES.PARSED_ALLELES' => this field int he database allows user to query by allele/genename in any case-senistiive
	#so when redirect_fields_for_improved_queries = True, select queries on fields such as 'DATA.VREGION.VGENES' will be renamed to the field designed for improved queires
	#modify_query_values_to_follow_db_schmea => this variable defines whther queries on specific fields should be modified to better match its value in the database. for example, most string fields in SEQS collection are capitalized. 
	#so when this variable is true, the function will look for any fields in the query that should be modified (i.e. whose values should be captialized)
	#when both modify_query_values_to_follow_db_schema=False AND redirect_fields_for_improved_queries=False, then NO MODIFICATION OF THE query occurs 
	#by default, just suppress DUPLICATED_FIELDS





#ONLY USE THIS FUNCTION FOR THE FOLLOWING:
	#DUMPING RESULTS IN NON MODIFIED JSON FORMAT (NOT IGREP FORMAT)	
#Process cursor for output:
	#Remove 'DATA' key from documents
	#SIMPLY USE JSON.DUMPS to maintain nested dictionary AND data formats (i.e. lists remain lists)
def Simple_Process_Output(document):
	
	#IF DOCUMENT IS {'_id':ObjectId(),'DATA':{'VREGION':{'VGENES':[a]},'DREGION':,'CDR3':}}, this line of code will do the following:
	#1) remove data key from document (moves everything underneath/nested-within the data key up 1 level)
	#2) flatten results from dictionary so it becomes unnested: {'_id':ObjectId(),'VREGION.VGENES':[a],'DREGION':,'CDR3'}			
		#convert ojbectId values to string
		#for k,v in document.iteritems():
		#	if isinstance(v,test):
		#		document[k]=str(v)
	
		#only convert certain keys to object-id	
	#for key in keys_with_object_id:		
	#	if key in document and isinstance(document[key],oid_type):			
	#		document[key] = str(document[key])
	RemoveObjId(document)
	
	data_info=document.pop('DATA',None)
	if data_info:
		document.update(data_info)
			
	#dump data 
	try:		
		return json.dumps(document)
	except:
		#try to use bson dump 
		return bson_dumps(document)

#ONLY USE THIS FUNCTION FOR THE FOLLOWING:
	#CONVERTING TO FASTA,FASTQ,CSV,TAB,IGREP

#Process cursor for output: 		
	#Remove 'DATA' key from documents
	#flatten all documents so they are not nested -> while flattening, convert any objectids to str()
	#FOR each key, convert values to 'STRING' output (for example split lists by ',')
def Process_Cursor_For_Output_File(document,schema_output_fnc = schema_fields_to_file):
	t=time.time()
	#if not isinstance(document,dict):
	#	return [document,all_keys_dict]		
			
	#print type(document)
	#IF DOCUMENT IS {'_id':ObjectId(),'DATA':{'VREGION':{'VGENES':[a]},'DREGION':,'CDR3':}}, this line of code will do the following:
	#1) remove data key from document (moves everything underneath/nested-within the data key up 1 level)
	#2) flatten results from dictionary so it becomes unnested: {'_id':ObjectId(),'VREGION.VGENES':[a],'DREGION':,'CDR3'}
		#1 million lines => 1 second to perform both step
	data_info=document.pop('DATA',None)
			
	if data_info:
		document.update(data_info)		
	
	tr = time.time()
	
	#1 million lines => takes 30 seconds to perform step using CYTHON, takes 1 minute to perform step using python code
	document = flatten_dictionary(document)			
	tf = time.time()
	
	#result_document = OrderedDict()
	#for keys in all_keys_dict:
	#	result_document[keys] = ""	
	
	#force values from docs to be strings
	for field,value in document.iteritems():					
		if not isinstance(value,basestring): 
			#its not a string, we need to modify value to string using schema_fields_to_file
			document[field] = schema_output_fnc[field](value)											
			
	tb = time.time()	
	global global_time
	global gt1
	global gt2
	global gt3
	global gt4
	global_time+=tb-t
	gt1+=tr-t
	gt2+=tf-tr
	gt3+=tb-tf
	
	
	return document


#Creating an instance of the class
#db_connect_data=>tuple defining authentication for database => four element tuple: database path, database user, database password, immunogrep user name
#modify_query_values_to_follow_db_schema=>modify the query passed in by user. it will properly format query values to match those in original schema 
#redirect_fields_for_improved_queries=>modify the query to change which fields are queried.
#to_file => whether or not results will be saved to a file a the end (as opposed to saving results as json documents and automatically loading them in memory at the end of query)
#file_prefix => if to_file = True, then all generated files will have this prefix
#proxy_path 
class RunQuery(object):
	def __init__(self,db_connect_data=None,modify_query_values_to_follow_db_schema=True,redirect_fields_for_improved_queries=True,to_file=True,file_prefix=None,proxy_path = ''):									
		self.to_file = to_file
		db_user_name='' #we use this variable to authenticate the person coming into the proxy 
		
		if self.to_file == False:						
			self.file_prefix = ''						
		else:
			if file_prefix==None:
				self.file_prefix = 'IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')															
			else:
				self.file_prefix = file_prefix		
		
		#most of the following parameters will not matter if proxy_path == None
		self.db_method = 'reader'										
		self.modify_query_values_to_follow_db_schema = modify_query_values_to_follow_db_schema
		self.redirect_fields_for_improved_queries = redirect_fields_for_improved_queries
		self.delim_file_headers = []
		#self.to_file=to_file
#		self.filename_suffix=filename_suffix
#		self.filename=filename
		self.query_command_list = [] #list of functions to execute on proxy side 		
		self.process_cursor_fnc = None				
		self.query_results = None					
							
		self.debug_file ='query_debug_log.txt'		
		self.exp_metadata = {}
		self.allowed_experiment_ids = []		
											
		#END of potentially unncecessary variables
				
		#if proxy_path is None then we want to actually run database queries and have passed in correct database information using the db_connect_data 
		if proxy_path == None:				
			if(len(db_connect_data)>4 or len(db_connect_data)<3):
				raise Exception('db_connect_data must be a tuple of either 3 or 4 elements')
			if len(db_connect_data)==4:		
				db_user_name = db_connect_data[3]		
			else:
				db_user_name = ''		
			self.proxy_path = None			
			self.name=None # we use self.name as a method for 'reconnecting' to the same varialbe via proxy. this is not implemented however
			#connect to database
			[self.db_path,self.mongo_connect]= connectToIgDatabase(db_connect_data[0],db_connect_data[1],db_connect_data[2]) #igdbconnect						
			#get user permissions
			self.user = self.db_path.users.find_one({'user':db_user_name})
			if not self.user:
				self.user = defaultdict(str)
				self.no_user_found = True
			else:
				self.user = defaultdict(str,self.user) # copy.deepcopy(user_info)							
				self.no_user_found = False
			
			#createa  log file
			if not os.path.isdir(db_user_name):
				os.mkdir(db_user_name)
			self.logged_output = open(db_user_name+'/logfile.txt','a')	
			self.logged_output.write('####################\nUser {0} has logged on to proxy at time: {1}.\n'.format(db_user_name,str(datetime.now())))
			
			self.location = 'proxy'			
			#get all metadata user has access to 
			self.exp_metadata = self.query_exp_docs_with_read_access()
			#get a list of experiments user has access to 
			self.allowed_experiment_ids = [ObjectId(id) for id in self.exp_metadata]			
			
			self.name=None
		
		#IN THE FOLLOWING CASES: WE WANT TO RUN THE FUNCTIONS VIA A PROXY THAT HAS ACCESS TO A DATABASE FROM ANOTHER ADDRESS/SOURCE
		elif proxy_path == '':			
			self.proxy_path = "http://biotseq.ut.appsoma.com:5998"									
			self.location = 'appsoma'			
			self.name = str(datetime.now()) # we use self.name as a method for 'reconnecting' to the same varialbe via proxy. this is not implemented however
			for bad_char in [':','-',' ','.']:
				self.name = self.name.replace(bad_char,'')												
			#we want to run these functions using the proxy defined by proxy path. 
			#so run the MODIFY_FUNCTION (see def modify_fucntion below) on all the function names in this file 
			#rather than running any functions defined for the class, it will instead generate a list of functions and parameters to run on the proxy 
			class_methods_to_modify = [obj for name, obj in inspect.getmembers(RunQuery) if inspect.ismethod(obj) and not (name.startswith('_') or name.startswith('wraps'))] 		
			for method in class_methods_to_modify:
				setattr(RunQuery, method.__name__, _modify_function(method))										
		else:
			self.proxy_path = proxy_path
			self.location = 'appsoma'			
			self.name = str(datetime.now()) # we use self.name as a method for 'reconnecting' to the same varialbe via proxy. this is not implemented however
			for bad_char in [':','-',' ','.']:
				self.name = self.name.replace(bad_char,'')												
			#we want to run these functions using the proxy defined by proxy path. 
			#so run the MODIFY_FUNCTION (see def modify_fucntion below) on all the function names in this file 
			#rather than running any functions defined for the class, it will instead generate a list of functions and parameters to run on the proxy 
			class_methods_to_modify = [obj for name, obj in inspect.getmembers(RunQuery) if inspect.ismethod(obj) and not (name.startswith('_') or name.startswith('wraps'))] 		
			for method in class_methods_to_modify:
				setattr(RunQuery, method.__name__, _modify_function(method))																			
				
	def __delete_query_object__(self): #delete object on proxy
		if self.location == 'proxy': 
			query_objects.pop(self.name,None)				
		
	def __iter__(self):		
		self.query_results = post_to_proxy(self) #OK now lets call proxy and send over the list of commands we will run
		return iter(self.query_results)				

	
	#THIS FUNCTION IS ONLY TO BE USED WHEN USING PROXY (WHEN PROXY_PATH !=None) WITHIN APPSOMA 		
	#filename_suffix=>how to end the filename. fileextension	
	#WHEN to_file = True then the following parameters are used:
		#max_doc_per_file=>maximum number of documents in each output file
	#WHEN to_file = False, then the following parameters are used 
		#chunk_size=> defines how many documents to return at one time to memory	
	def _return(self,chunk_size=10000,to_file=True,filename=None,max_doc_per_file=0,filename_suffix='query'):					
		if self.proxy_path == None:
			raise Exception("Error: You cannot use the function _return when proxy_path is set to None. Either set proxy_path to an address referring to the database OR use teh _return_results function")
			
		self.to_file = to_file		
		self.file_prefix = filename				
		return post_to_proxy(self,self.proxy_path,chunk_size,filename_suffix,max_doc_per_file)			
			
		
	#wfile can be either a filebuffer set in write mode OR the http proxy variable wfile
	#to_file => output query result to a file 
	def _return_results(self,wfile):														
		#no results found		
		
		if not self.query_results:			
			wfile.write('[]')
		else:		
			#if self.query_to_file: #if we are writing the results to a file, then we will iterate through query_results, and, after each iteration, add a new_line character (delim character)
	#				start_char = ''
	#				end_char = ''
	#				delim_char = '\n' #for files we seperate results using new line
	#			else: #but if we are returning the results to memory, then instead we will iterate through query_reulsts, and return all rsults a list of documents/results
	#				start_char = '[' #for conversion to list in memory
	#				end_char = ']'
	#				delim_char = ',' #for memory, we seperate results using ',' so that can be read in as list
			
			#wfile.write(start_char) 
			if isinstance(self.query_results,dict): #usually happens with a find_one command (it doesnt return a cursor that is)			
				wfile.write(Simple_Process_Output(self.query_results)+'\n')
				#If returned to memory it will always be returned as a list, if its to file, then we will not encapuslate data in list
				#[document,self.unique_field_keys] = Process_Cursor_For_Output(self.query_results,self.unique_field_keys,self.make_output_strings)
				#wfile.write(start_char+json.dumps(document)+end_char) #then just dump dictionary and write results.  
								
			elif isinstance(self.query_results,basestring): #if the results are simply string/text info
				wfile.write(self.query_results+'\n')
				#If returned to memory it will always be returned as a list, if its to file, then we will not encapuslate data in list
				#wfile.write(start_char+self.query_results+end_char) #its just a string, so just return that
						
			elif isinstance(self.query_results,list): #if its a list, then just dump the contents of the list to a file or to memory
				self.query_results = [Simple_Process_Output(list_result) if isinstance(list_result,dict) else list_result for list_result in self.query_results]
				wfile.write(json.dumps(self.query_results)+'\n')
				#wfile.write(json.dumps(self.query_results)) #no need for start_char and end_char in this situation because it is encapsulated by the list already
			
			elif isinstance(self.query_results,int) or isinstance(self.query_results,float) or isinstance(self.query_results,long) or isinstance(self.query_results,complex):				
				wfile.write(json.dumps(self.query_results)+'\n')
				#wfile.write(json.dumps(self.query_results))
			
			else: # if its a cursor or a generator, then we can iterate through it with a "next()" function 									
				
				#diry = []
				t1 = 0
				aa = time.time()
				total_s = ''
				for count,result_document in enumerate(self.query_results):					
					if result_document==None:
						continue
					#if this document/result is a string, then no need for dumps
					if isinstance(result_document,basestring): 
						wfile.write(result_document)																		
						
						#total_s+=result_document
						#if count%chunk_size==0:
							#wfile.write(total_s)
							#total_s=''
					elif isinstance(result_document,dict):							
						#if its a dictionary, then we assume that the user did not run through any CONVERT_TO_FASTA/FASTQ/TAB/IGREP GENERATOR FUNCTIONS (THESE GENERATORS RETURN STRINGS)
						#THEREFORE, WE JUST NEED TO DUMPS THE RESULTS AS IS 						
						wfile.write(Simple_Process_Output(result_document)+'\n')																							
						#doc=result_document
						#seq_doc = Process_Cursor_For_Output_File(result_document)
						#wfile.write(json.dumps(seq_doc)+'\n')						
						#document=result_document
					else:
						#A non-dict document was return, so just use json.dumps to format it correctly
						try:
							wfile.write(json.dumps(result_document)+'\n')																	
						except:
							wfile.write(bson_dumps(result_document)+'\n')																	
				#t1 = time.time()		
				#for i in range(len(diry)):
				#	Simple_Process_Output(diry[i])				
				#t2 = time.time()					
				t1+=time.time()-aa
				print t1
				#print global_time
				print gt1
				print gt2
				print gt3
		self.logged_output.write('Successfully executed function at time: {0}.\n#################################\n'.format(str(datetime.now())))
		self.logged_output.close()
	def WriteToDebugFile(self,text):
		with open(self.debug_file,'a') as w:
			w.write(str(text)+'\n')
		
					
	##############QUERY FUNCTIONS##########################
	
	#mongo cursor functions
	def count(self,with_limit_and_skip=False):
		try:
			self.query_results = self.query_results.count(with_limit_and_skip)
		except Exception as e:
			raise Exception("count function can only be used on a mongodb cursor! Error message: "+str(e))
		return self
	
	#set a limit for number of results 
	def limit(self,limit):
		try:
			self.query_results = self.query_results.limit(limit)
		except Exception as e:
			raise Exception("Limit function can only be used on a mongodb cursor! Error message: "+str(e))
		return self
	
	#set a skip command for general queries
	def skip(self,num_skip):
		try:
			self.query_results = self.query_results.skip(num_skip)
		except Exception as e:
			raise Exception("Skip function can only be used on a mongodb cursor! Error message: "+str(e))
		return self
	
	#set up a sort command for general queries
	def sort(self,sort_tuples):
		try:
			self.query_results = self.query_results.sort(sort_tuples)
		except Exception as e:
			raise Exception("Sort function can only be used on a mongodb cursor! Error message: "+str(e))
		return self 
		
	#for mongo queries debugging. returns a dictionary explaining the query performance
	def explain(self):
		try:
			self.query_results = self.query_results.explain()
		except Exception as e:
			raise Exception("Explain function can only be used on a mongodb cursor! Error message: "+str(e))
		return self
	
	def distinct(self,field_name):
		try:
			self.query_results = self.query_results.distinct(field_name)
		except Exception as e:
			raise Exception("Distinct function can only be used on a mongodb cursor! Error message: "+str(e))		
		return self
	
	#QUERIES ON USERS COLLECTION#
	#returns list of all users currrently registered with the database
	def get_db_user_list_info(self):
		db = self.db_path.users 				
		
		if self.user['administrator']:
			self.query_results =[u for u in db.find()]
		else:
			#if not administrator then just report username, full name, and lab
			self.query_results = [u for u in db.find({},{'_id':0,'user':1,'name':1,'lab':1})]
			self.delim_file_headers =[]		
		self.delim_file_headers = ['username','name','administrator','curator','write_access','lab','email']		
		self.list_to_gen()
		return self
		
	#returns information about current user 
	def get_user_access_info(self):								
		self.query_results = copy.deepcopy(self.user)						
		if self.no_user_found:
			self.user['user'] = ''
		self.delim_file_headers = ['username','name','administrator','curator','write_access','lab','email']		
		return self
	
	#return object IDs converted to string 
	def get_accessible_exp_ids(self):				
		self.query_results = [str(id) for id in self.allowed_experiment_ids]
		self.delim_file_headers = ['EXPERIMENT_IDS']
		self.list_to_gen()
		return self
		
	#return object IDs converted to string 
	def get_exp_ids_owned_by_current_user(self):				
		exp_data = copy.deepcopy(self.exp_metadata)
		write_access_exps = [str(id) for id in exp_data if self.user['user'].lower() in exp_data[id]['OWNERS_OF_EXPERIMENT']]
		self.query_results = copy.deepcopy(write_access_exps)		
		self.delim_file_headers = ['EXPERIMENT_IDS']	
		self.list_to_gen()
		return self
		
	#return object IDs converted to string
	def get_write_accessible_exp_ids(self):
		if self.user['administrator']: #return all allowed experiments
			self.query_results = [str(id) for id in self.allowed_experiment_ids]
		else: #only return experiments with write access 
			self.get_exp_ids_owned_by_current_user()
		self.delim_file_headers = ['EXPERIMENT_IDS']		
		self.list_to_gen()
		return self
	
	#QUERIES ON metadat/exps COLLECTION#
	
	#return experiment metadata of experiments with read access
	#write_access = only return documents with write access	
	#administrators will recieve all experiment metadata
	#object_id_list => filter experiment metadata to return by experiment id
	def get_exp_docs(self,object_id_list = [], write_access_only = False):				
		
		if object_id_list:			
			object_id_list = list(set(self.allowed_experiment_ids)&set([convert_to_objectid(o) for o in object_id_list]))
		else:
			object_id_list = self.allowed_experiment_ids			
		if self.user['administrator']:
			self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata if ObjectId(ids) in object_id_list]
		else:
			if write_access_only:
				self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata if ObjectId(ids) in object_id_list and self.user['user'].lower() in self.exp_metadata[ids]['OWNERS_OF_EXPERIMENT']]
			else:
				self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata if ObjectId(ids) in object_id_list]
		self.delim_file_headers = []
		for each_doc in self.query_results:
			self.delim_file_headers.extend(each_doc.keys())
		self.delim_file_headers=list(set(self.delim_file_headers))		
		
		self.list_to_gen()

		return self
	
	#return experiment metadata of experiments a user included in the field 'OWNERS_OF_EXPERIMENT'
	#similar to previous function with write_access_only = True, except will also only return write_acces_exps for administrators
	def get_exp_docs_with_ownership(self,object_id_list=[]):
		if object_id_list:
			object_id_list = list(set(self.allowed_experiment_ids)&set([convert_to_objectid(o) for o in object_id_list]))
		else:
			object_id_list = self.allowed_experiment_ids			
		self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata if ObjectId(ids) in object_id_list and self.user['user'].lower() in self.exp_metadata[ids]['OWNERS_OF_EXPERIMENT']]
		
		self.delim_file_headers = []
		for each_doc in self.query_results:
			self.delim_file_headers.extend(each_doc.keys())
		self.delim_file_headers=list(set(self.delim_file_headers))		
		self.list_to_gen()
		return self
	
	
	#if analyses is empty, then returns only documents with sequences in experiment
	#if analyses is not empty, then returns only documents with squences containing annotation info from provided experiments 
	def get_exp_docs_containing_sequences(self,write_access_only=False,analyses = []):
		db = self.db_path	
		object_id_list = self.allowed_experiment_ids
		query = {'_id':{'$in':object_id_list},'SEQ_COUNT':{'$gt':0}}
		if write_access_only and self.user['administrator']==False:
			query['OWNERS_OF_EXPERIMENT'] = self.user['user'].lower()
		if analyses and analyses is not list:			
			analyses = [analyses]
		for a in analyses:			
			query['ANALYSES_COUNT.{0}'.format(a.upper())] = {'$gt':0}
		self.query_results = [exp for exp in db.exps.find(query,{'DUPLICATED_FIELDS':0})]
	
		self.delim_file_headers = []
		for each_doc in self.query_results:
			self.delim_file_headers.extend(each_doc.keys())
		self.delim_file_headers=list(set(self.delim_file_headers))		
		self.list_to_gen()
		return self
		
	#QUERY experiment collections using standard pymongo .find function
	#exp_id = > filter results by experiemnts with ObjectId listed in exp_id
	#fields for query are defined by 'q'
	#p => what fields should be projected	
	def query_exps_collection(self, exp_id = None, q={},p={'DUPLICATED_FIELDS':0},include_duplicated_fields = False):		
		if exp_id == None:
			#default var which means undefined, so pass in allowed_experiment_ids
			#basically assume that user does not knwo what to search,so wants to search all avaialabe experiments 
			exp_id = self.allowed_experiment_ids
		
		if not isinstance(exp_id,list):
			exp_id = [exp_id]
			
		#logged_output.write('The following query in function "Query_Exps_Collection" was desired: {0}\n'.format(str(q)))
		db = self.db_path
		
		#modify the query to improve success of a match 
		#redirect select fieldname to fields created by server to imporve queries
		#ensure that the format of the values provided to the field matches the format of the field put in database			
		
		q = Parse_Mongo_Query_Expression(q, dict_defining_value_transformation=fields_for_queries_exps_collection,redirection_fields = redirect_exp_collection_fields,modify_query_values_to_follow_db_schema=self.modify_query_values_to_follow_db_schema,redirect_fields_for_improved_queries=self.redirect_fields_for_improved_queries)
	
	
		
		#we will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)		
		q['_id'] = GetIdIntersection([self.allowed_experiment_ids,exp_id],q.pop('_id',None))
		
		self.logged_output.write('Query was modified to: {0}\n'.format(str(q)))							
		project = copy.deepcopy(p)		
				
		if project:
			allTrue = sum(project.values())
			if allTrue==len(project): #the user projected fields using '1'. they listed select fields to project
				project['_id'] = 1 #also always project '_id' field 
				if include_duplicated_fields:
					project['DUPLICATED_FIELDS'] = 1
			elif allTrue==0:#the user suppressed feilds using '0'
				if not(include_duplicated_fields):
					project['DUPLICATED_FIELDS'] = 0
				else:
					#just incase user explicityl said to include_duplicated_fields
					project.pop('DUPLICATED_FIELDS',None)
		else:
			project = None
		
		self.logged_output.write('The following values will be projected: {0}\n'.format(str(project)))
		
		if type(q) is not dict:
			raise Exception("ERROR: query_command must be a dictionary for query experiment collection")																
		#if q:
		query = {f:v for f,v in q.iteritems()}# copy.deepcopy(q) ==> just in case re.expression is in there? 
		self.query_results =db.exps.find(query,project)
		
		
		return self
		
		#just incase the query was innaccurate, go through the metadat results and enusre '_id' returned is within allowed_experiment_ids
		
		#return self.query_results
		
		#else:
			
			#no query passed, so simply return metadata for all experiments 
		#	self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata]
		
	#RETURNS distinct values from fields within experiment collection
	#ONLY returns values from CURATED experiment documents
	#DOES not return unique values from all experiments if requested field (unique_fields) is not found within allow_unique_metadata	
	#unique_fields => list of fields the user wants to show distinct values for 
	def get_distinct_exp_collection_values(self,unique_fields=[]):
		db = self.db_path	
		allow_unique_metadata = ['PAIRING_TECHNIQUE','CELL_MARKERS_USED','LIST_OF_POLYMERASES_USED','POST_SEQUENCING_PROCESSING: PHI_X_FILTER', 
						'POST_SEQUENCING_PROCESSING: QUALITY_FILTER', 'POST_SEQUENCING_PROCESSING: PROCESS_R1_R2_FILE',
						'SPECIES','SEQUENCING_PLATFORM','CHAIN_TYPES_SEQUENCED','CELL_TYPES_SEQUENCED','ISOTYPES_SEQUENCED',
						'TEMPLATE_TYPE','REVERSE_PRIMER_USED_IN_RT_STEP']
		all_experiments_metadata = db.exps.find({'CURATED':True}) 
		#self.query_exp_collection({'_id':{'$in':self.allowed_experiment_ids},'CURATED':True})# db.expreal.find({ self.get_exp_ids_for_current_user()
		allowed_experiments_metadata = db.exps.find({'_id':{'$in':self.allowed_experiment_ids},'CURATED':True})
		if unique_fields:
			if type(unique_fields) is not list:
				unique_fields = [unique_fields]			
			
			if self.user['administrator']:
				allow_unique_metadata = unique_fields
								
			self.query_results = {field:all_experiments_metadata.distinct(field) if field in allow_unique_metadata else allowed_experiments_metadata.distinct(field) for field in unique_fields}			
		else:
			self.query_results = None
			
		self.delim_file_headers = unique_fields
		return self
	
	#PRIVATE FUNCTION..Cannot be used by user
	def query_exp_docs_with_read_access(self):	#search in appsoma this function 
		self.query_results = None
		db = self.db_path				
		if self.user['administrator']:
			exp_docs = {str(exp['_id']):exp for exp in db.exps.find({})}				
		else:
			read_access_names = ['all'] + [self.user['user'].lower()] + ['lab_'+lab.lower() for lab in self.user['lab'] if lab]					
			exp_docs = {str(exp['_id']):exp for exp in db.exps.find({'READ_ACCESS':{'$in':read_access_names}})}				
		return exp_docs
	
	#PRIVATE FUNCTION..Cannot be used by user
	def	query_exp_docs_with_write_access(self):	#search in appsoma this function 	
		self.query_results = None
		db = self.db_path		
		if self.user['administrator']:
			exp_docs = {str(exp['_id']):exp for exp in db.exps.find({})}				
		else:
			write_access_name = self.user['user'].lower()
			exp_docs = {str(exp['_id']):exp for exp in db.exps.find({'OWNERS_OF_EXPERIMENT':write_access_name})}		
		return exp_docs
	
	########SEQUENCE COLLECTION QUERIES#####							
	
	#GENERAL QUERY FUNCTION#
	#QUERY seq collections using standard pymongo find function
	#exp_id = > filter results by experiemnts with ObjectId listed in exp_id
	#fields for query are defined by 'query'
	#project => what fields should be projected
	
	#analysis_name => a list of analysis types to filter by. note -> more complex queries using analysis_name can be defined in the query field 
	#recombination_type => a list of recombination_types to filter by. note -> more complex queries using recombination_type can be defined in the query field 
	#metadata_query => in addition to explicity passing in a list of ObjectId's, a user can further filter experiemnts by querying teh exps collection/metadata
		#for example, if metadata_query = {'PROJECT_NAME':'Demo'}, then first the this metadata will be queried to return exp_ids whose PROJECT_NAME contain 'Demo' 	
	def query_seqs_collection(self,exp_id=None, analysis_name = [], recombination_type = [], query={},project={'QUERY_DATA':0},metadata_query={},include_original_ngs_seq=False,limit=0):
		
		self.query_results = None
		self.logged_output.write('The following query in function "Query_Seqs_Collection" was desired: {0}\n'.format(str({'experiments':exp_id,'query':query,'metadata_query':metadata_query})))		
		db = self.db_path#db_reader						
		
		q = {f:v for f,v in query.iteritems()} #make a copy of query 
		
		fields_to_project = copy.deepcopy(project)		
		if fields_to_project:
			allTrue = 0
			counts_in = 0
			for p,v in fields_to_project.iteritems():
				if isinstance(v,dict):
					pass
				elif v == 0:
					counts_in+=1
				else:
					counts_in+=1
					allTrue+=1		
			
			#allTrue = sum(fields_to_project.values())
			if allTrue==counts_in: #the user projected fields using '1'. they listed select fields to project
				allTrue=True
				fields_to_project['_id'] = 1 #also always project '_id' field 
				fields_to_project[idIdentifier] = 1 #also always project 'SEQ_ID' field 
				fields_to_project[expIdentifier] = 1 #also always project 'EXP_ID' field 				
				#fields_to_project['SETTINGS'] = 1 #also always project 'EXP_ID' field 				
				fields_to_project['ANALYSIS_NAME'] = 1 #also always project 'EXP_ID' field 				
				fields_to_project['RECOMBINATION_TYPE'] = 1 #also always project 'EXP_ID' field 				
				#fields_to_project['DATE_UPDATED'] = 1 #also always project 'EXP_ID' field 				
			elif allTrue==0:#the user suppressed feilds using '0'
				allTrue=False
				fields_to_project['QUERY_DATA'] = 0
		else:
			allTrue=True
			fields_to_project = None
			
		if exp_id == None: #if exp_id is passed in as empty, then make it self.allowed_experiment_ids
			exp_id = self.allowed_experiment_ids
		elif not isinstance(exp_id,list):
			exp_id = [exp_id]
						
		if metadata_query:		
			#query exp_collection
			#result of query will be stored in class name 'self.query_results'			
			self.query_exps_collection(exp_id,q={f:v for f,v in metadata_query.iteritems()})				
			if self.query_results:
				self.exp_metadata = {str(exp['_id']):exp for exp in self.query_results if exp['_id'] in self.allowed_experiment_ids}
				#erase metadata query 							
				exp_id = [vals['_id'] for meta,vals in self.exp_metadata.iteritems()]
				self.query_results = None
			else:
				self.exp_metadata = {}
				#no experiments to searchby. no experiments were found !
				exp_id = []			
						
		if analysis_name:
			if not isinstance(analysis_name,list):
				q['ANALYSIS_NAME']=analysis_name
			elif len(analysis_name)==1:
				q['ANALYSIS_NAME']=analysis_name[0]
			else:
				q['ANALYSIS_NAME']={'$in':analysis_name}
		
		if recombination_type:
			if not isinstance(recombination_type,list):
				q['RECOMBINATION_TYPE']=recombination_type
			elif len(recombination_type)==1:
				q['RECOMBINATION_TYPE'] = recombination_type[0]			
			else:
				q['RECOMBINATION_TYPE']={'$in':recombination_type}
				
		
		#modify the query to improve success of a match 
		#redirect select fieldname to fields created by server to imporve queries
		#ensure that the format of the values provided to the field matches the format of the field put in database			
		q = Parse_Mongo_Query_Expression(q, dict_defining_value_transformation=fields_for_queries_seqs_collection,redirection_fields = redirect_seq_collection_fields,modify_query_values_to_follow_db_schema=self.modify_query_values_to_follow_db_schema,redirect_fields_for_improved_queries=self.redirect_fields_for_improved_queries)
		
		#we will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)
		#there are exp_ids passed in by users OR found by the metadata_query above
		#therefore, find the intersection between results and allowed id's		
		#THIS WILL ALSO CONSIDER ANY TOP_LEVEL EXP_ID requests in the query field 
		q[expIdentifier] = GetIdIntersection([self.allowed_experiment_ids,exp_id],q.pop(expIdentifier,None))
		
		
		#based on the actual experiment IDs used in query, figure out which fields will be projected by using analysis schema in metadata
		if isinstance(q[expIdentifier],dict): #the query is dictionary using '$in' command
			possible_ids = q[expIdentifier]['$in']
		else:
			possible_ids = [q[expIdentifier]]
		possible_metadata_fields = [self.exp_metadata[str(id)] for id in possible_ids]
		[schema_fields,possible_analyses,possible_recombination_types] = Get_Schema_Details(possible_metadata_fields,fields_to_project,allTrue)
		
		if 'ANALYSIS_NAME' not in q and possible_analyses!=[]:
			if len(possible_analyses)==1:
				q['ANALYSIS_NAME']= possible_analyses[0]
			else:
				q['ANALYSIS_NAME']={'$in': possible_analyses}
		
		if 'RECOMBINATION_TYPE' not in q and possible_recombination_types!=[]:
			if len(possible_recombination_types)==1:
				q['RECOMBINATION_TYPE']=possible_recombination_types[0]
			else:
				q['RECOMBINATION_TYPE']={'$in': possible_recombination_types}
												
		self.logged_output.write('Query was modified to: {0}\n'.format(str(q)))							
		self.logged_output.write('The following values will be projected: {0}\n'.format(str(fields_to_project)))
	
		self.delim_file_headers = schema_fields
		self.query_results =db.seqs.find(q,fields_to_project).limit(limit)
		
		if include_original_ngs_seq:
			self.delim_file_headers.extend(data_fields_for_raw_data)
			#for every result in the query, also get the raw sequence data associated with result 
			self.get_rawseq_info()		

		return self
	
	#GENERAL AGGREGATION PIPELINE FUNCTION#
	#QUERY seq collections using standard pymongo find function
	#exp_id = > filter results by experiemnts with ObjectId listed in exp_id	
	#aggregation_fxn => a list of mongo aggregation commands to run pipeline 
		#IMPORTANT=> FIRST ELEMENT IN AGGREGATION FUNCTION MUST BE '$MATCH'. IF THE USER DOES NTO PASS IT IN, THEN THIS FUNCTION WILL ADD A $MATCH FILTER USE EXP_ID
	
	#initial_limit => this value tells the program to insert a '$limit' command right after the first filter/$match command in the aggregation function
	
	#analysis_name => a list of analysis types to filter by. note -> more complex queries using analysis_name can be defined in the query field 
	#recombination_type => a list of recombination_types to filter by. note -> more complex queries using recombination_type can be defined in the query field 
	#metadata_query => in addition to explicity passing in a list of ObjectId's, a user can further filter experiemnts by querying teh exps collection/metadata
		#for example, if metadata_query = {'PROJECT_NAME':'Demo'}, then first the this metadata will be queried to return exp_ids whose PROJECT_NAME contain 'Demo' 	
	
	def aggregate_seqs_collection(self,exp_id=None, analysis_name = [], recombination_type = [],aggregation_fxn =[],initial_limit=0,metadata_query={},allowDiskUse=True):	
		self.query_results = None
		self.logged_output.write('The following aggregation function in "Aggregate_Seqs_Collection" was desired: {0}\n'.format(str({'experiments':exp_id,'aggregation':aggregation_fxn,'metadata_query':metadata_query})))		
		db = self.db_path#db_reader						
		
		aggregate = [pipe for pipe in aggregation_fxn] #make a copy of aggregation
		
		if exp_id == None: #if exp_id is passed in as empty, then make it self.allowed_experiment_ids
			exp_id = self.allowed_experiment_ids
		elif not isinstance(exp_id,list):
			exp_id = [exp_id]
		
		#FIRST CHECK WHETHER THE FIRST ELEMENT OF AGGREGATION FUNCTION IS A MATCH COMMAND, IF NOT, ADD IN A MATCH COMMAND 
		if aggregate[0].keys()!=['$match']:			
			aggregate.insert(0,{'$match':{}}) #prepend and empty filter command to top of list, we will add exp_ids to this later
									
		if metadata_query:		
			#query exp_collection
			#result of query will be stored in class name 'self.query_results'			
			self.query_exps_collection(exp_id,q={f:v for f,v in metadata_query.iteritems()})				
			if self.query_results:
				self.exp_metadata = {str(exp['_id']):exp for exp in self.query_results if exp['_id'] in self.allowed_experiment_ids}
				#erase metadata query 							
				exp_id = [vals['_id'] for meta,vals in self.exp_metadata.iteritems()]
				self.query_results = None
			else:
				self.exp_metadata = {}
				#no experiments to searchby. no experiments were found !
				exp_id = []			
						
		if analysis_name:
			#add to first fitler in aggregation function
			aggregate[0]['$match']['ANALYSIS_NAME'] = analysis_name if not isinstance(analysis_name,list) else {'$in':analysis_name} #if analysis_name is a dict or a string, then we do not modify query. if its a list, we add an '$in' command 
									
		if recombination_type:
			#add to first fitler in aggregation function
			aggregate[0]['$match']['RECOMBINATION_TYPE'] = recombination_type if not isinstance(recombination_type,list) else {'$in':recombination_type}
					
		
		#for aggregation functions, we only modify fields and field names in the first MATCH command. 
		#We trust the user to handle the remaining match commands. this is mainly assuming that the user will transform document sin aggregation
		#So we will not know how to process the filter commands 		
		#CONVERT query to proper format/objects using bson loads 
		#modify the query to improve success of a match 
		#redirect select fieldname to fields created by server to imporve queries
		#ensure that the format of the values provided to the field matches the format of the field put in database			
		aggregate[0]['$match'] = Parse_Mongo_Query_Expression(aggregate[0]['$match'], dict_defining_value_transformation=fields_for_queries_seqs_collection,redirection_fields = redirect_seq_collection_fields,modify_query_values_to_follow_db_schema=self.modify_query_values_to_follow_db_schema,redirect_fields_for_improved_queries=self.redirect_fields_for_improved_queries)
	
		#THIS LINE SHOULD ALWAYS COME AFTER Parse_Mongo_Query_Expression because sometimes object id will be {'oid':""} rather than ObjectId(). previousf ucntion fixes it
		
		#we will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)
		#there are exp_ids passed in by users OR found by the metadata_query above
		#therefore, find the intersection between results and allowed id's		
		#THIS WILL ALSO CONSIDER ANY TOP_LEVEL EXP_ID requests in the query field 
		aggregate[0]['$match'][expIdentifier] = GetIdIntersection([self.allowed_experiment_ids,exp_id],aggregate[0]['$match'].pop(expIdentifier,None))
		
		#based on the actual experiment IDs used in query, figure out which fields will be projected by using analysis schema in metadata
		if isinstance(aggregate[0]['$match'][expIdentifier],dict): #the query is dictionary using '$in' command
			possible_ids = aggregate[0]['$match'][expIdentifier]['$in']
		else:
			possible_ids = [aggregate[0]['$match'][expIdentifier]]
		possible_metadata_fields = [self.exp_metadata[str(id)] for id in possible_ids]
		
		[schema_fields,possible_analyses,possible_recombination_types] = Get_Schema_Details_Aggregation_Pipeline(possible_metadata_fields,copy.deepcopy(aggregate))
		
			
		if 'ANALYSIS_NAME' not in aggregate[0]['$match'] and possible_analyses!=[]:
			if len(possible_analyses)==1:
				aggregate[0]['$match']['ANALYSIS_NAME']= possible_analyses[0]
			else:
				aggregate[0]['$match']['ANALYSIS_NAME']={'$in': possible_analyses}
		
		if 'RECOMBINATION_TYPE' not in aggregate[0]['$match'] and possible_recombination_types!=[]:
			if len(possible_recombination_types)==1:
				aggregate[0]['$match']['RECOMBINATION_TYPE']=possible_recombination_types[0]
			else:
				aggregate[0]['$match']['RECOMBINATION_TYPE']={'$in': possible_recombination_types}
		
		
		self.delim_file_headers = schema_fields
		
		
		if initial_limit>0:
			#insert directly after first match 
			aggregate.insert(1,{'$limit':initial_limit})
		
		self.logged_output.write('AggregationPipeline was modified to: {0}\n'.format(json.dumps([str(pipe) for pipe in aggregate],indent=4)))
	
	
		self.query_results =db.seqs.aggregate(aggregate,allowDiskUse=allowDiskUse)			
		
	
		return self
	
	#SPECIALIZED QUERY FUNCTIONS#
	#function name handles most of the queries, should not be a requirement for mongo-db knowledge#
	#this function will return NGS reads for given experiments defined by exp_id
	#By Default, this function will only query NGS reads for each experiment defined in exp_id query
	#IF all annotated_data for a specific NGS read is desired, then include_all_data_by_seq_id = True
    	#=> IF TRUE, THEN THE FIELDS INCLUDE_ANALYSES AND INCLUDE_RECOMBINATION FIELDS ARE IGNORED
	#IF a select subset of annotated data for a specific NGS read is desired, then include a list of string names for
    	#a) include_analyses (analysis_types) and 
    	#b) include_recombination_type (recombniation types)
	#IF RESULTS from this output include annotated data, then results from annotated data will be sorted by SEQ_ID
	
	#PASS in a list of object ids corresponding to a experiments 
	#OUTPUTS sequences from each defined experiment and fields defined in 'project' variable
	#metadata_query => in addition to explicity passing in a list of ObjectId's, a user can further filter experiemnts by querying teh exps collection/metadata
		#for example, if metadata_query = {'PROJECT_NAME':'Demo'}, then first the this metadata will be queried to return exp_ids whose PROJECT_NAME contain 'Demo' 	
	def get_sequences_from_experiment(self,exp_id=None,project={'_id':1,'DATE_UPDATED':1,idIdentifier:1,expIdentifier:1,'DATA':1,'ANALYSIS_NAME':1,'RECOMBINATION_TYPE':1}, metadata_query = {} ,include_all_data_by_seq_id=False, include_analyses = None, include_recombination_type = None,limit=0, same_document=False):	
		self.query_results = None
		self.logged_output.write('The following query in function "Get_Sequences_From_Experiment" was desired: {0}\n'.format(bson_dumps({'experiments':exp_id,'additional_analyses':include_analyses,'include_recombination':include_recombination_type,'include_all_data':include_all_data_by_seq_id})))
		db = self.db_path#db_reader						
		
		if include_analyses!=None and not isinstance(include_analyses,list):
			include_analyses = [include_analyses]
		if include_recombination_type!=None and not isinstance(include_recombination_type,list):
			include_recombination_type = [include_recombination_type]
		
		fields_to_project = copy.deepcopy(project)		
		
		if fields_to_project:
			numTrue = 0
			counts_in = 0
			for p,v in fields_to_project.iteritems():
				if isinstance(v,dict):
					pass
				elif v == 0:
					counts_in+=1
				else:
					counts_in+=1
					numTrue+=1		
			
			#allTrue = sum(fields_to_project.values())
			if numTrue==counts_in: #the user projected fields using '1'. they listed select fields to project
				allTrue=True
				fields_to_project['_id'] = 1 #also always project '_id' field 
				fields_to_project[idIdentifier] = 1 #also always project 'SEQ_ID' field 
				fields_to_project[expIdentifier] = 1 #also always project 'EXP_ID' field 				
				fields_to_project['ANALYSIS_NAME'] = 1
				fields_to_project['RECOMBINATION_TYPE'] = 1
				fields_to_project['DATA.SEQUENCE'] = 1 #for this function, always project sequence
				fields_to_project['DATA.SEQUENCE_HEADER'] = 1 #for this function, always project sequence
				fields_to_project['DATA.QUALITY_SCORE'] = 1 #for this function, always project QUALITY_SCORE
			elif numTrue==0:#the user suppressed feilds using '0'
				allTrue=False
				fields_to_project['QUERY_DATA'] = 0
		else:
			allTrue=True
			fields_to_project = None
		
		
		
		if exp_id == None: #if exp_id is passed in as empty, then make it self.allowed_experiment_ids
			exp_id = self.allowed_experiment_ids
		elif not isinstance(exp_id,list):
			exp_id = [exp_id]				
		
		if metadata_query:		
			#query exp_collection
			#result of query will be stored in class name 'self.query_results'			
			self.query_exps_collection(exp_id,q={f:v for f,v in metadata_query.iteritems()})				
			if self.query_results:
				self.exp_metadata = {str(exp['_id']):exp for exp in self.query_results if exp['_id'] in self.allowed_experiment_ids}
				#erase metadata query 							
				exp_id = [vals['_id'] for meta,vals in self.exp_metadata.iteritems()]
				self.query_results = None
			else:
				self.exp_metadata = {}
				#no experiments to searchby. no experiments were found !
				exp_id = []			
				
		#we will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)
		#there are exp_ids passed in by users OR found by the metadata_query above
		#therefore, find the intersection between results and allowed id's
		exp_query = GetIdIntersection([self.allowed_experiment_ids,exp_id])
				
		self.logged_output.write('The following experiments were filtered: {0}\n'.format(str({expIdentifier:exp_query})))
		
		#based on the actual experiment IDs used in query, figure out which fields will be projected by using analysis schema in metadata
		if isinstance(exp_query,dict): #the query is dictionary using '$in' command
			possible_ids = exp_query['$in']
		else:
			possible_ids = [exp_query]
			
		
		#seqRawData analysis schema: 
		seq_analysis_schema = {'ANALYSIS_SCHEMA':{'SEQUENCE.ANALYSES':'RAW','SEQUENCE_HEADER.ANALYSES':'RAW','QUALITY_SCORE.ANALYSES':'RAW','FILENAME.ANALYSES':'RAW'}}
								
		if include_all_data_by_seq_id==False and include_analyses==None and include_recombination_type==None:
			#we only want sequences, do not want annotation data 
			possible_metadata_fields = [seq_analysis_schema]
			[schema_fields,possible_analyses,possible_recombination_types] = Get_Schema_Details(possible_metadata_fields,fields_to_project,allTrue)						
			
			#only query NGS reads where analysis_name = seqRawData variable
			#NO NEED TO SORT BY SEQID
			self.query_results = db.seqs.find({
				expIdentifier:exp_query,
				'ANALYSIS_NAME':seqRawData
			},fields_to_project).limit(limit)			
			
		elif include_all_data_by_seq_id==True:
			#we want all data associated with experiment 
			possible_metadata_fields = [self.exp_metadata[str(id)] for id in possible_ids]
			possible_metadata_fields.append(seq_analysis_schema)
			[schema_fields,possible_analyses,possible_recombination_types] = Get_Schema_Details(possible_metadata_fields,fields_to_project,allTrue)						
			total_docs_per_seq = len(possible_analyses)+1		
			#query ALL DATA from an experiment
			#sort by SEQ_ID
			#process output specially
			self.query_results = db.seqs.find({
				expIdentifier:exp_query
			},fields_to_project).sort(idIdentifier).hint('exp_seq_id_an_rt').limit(limit*total_docs_per_seq)						
			#create a generator to parse the cursor and group together documents by SEQ_ID
			self.group_documents_by_seq_id(limit,same_document)    
			
		elif include_all_data_by_seq_id == False:
			possible_metadata_fields = [self.exp_metadata[str(id)] for id in possible_ids]
			possible_metadata_fields.append(seq_analysis_schema)
			[schema_fields,possible_analyses,possible_recombination_types] = Get_Schema_Details(possible_metadata_fields,fields_to_project,allTrue)						
			total_docs_per_seq = len(possible_analyses)+1					
			#QUERY ONLY SEQUENCE DATA OR SEQUENCES WITH 
			#SPECIFIED ANALYSIS TYPES AND RECOMBINATION TYPES                        
			if include_recombination_type and possible_recombination_types!=include_recombination_type: #ADD A RECOMBINATION_TYPE FILTER , but if there is only one recombination type (i.e. include_recombination_type = VDJ and only VDJ is found in possible_recombination_types, then filtering by recombination type is unnecessary)
				or_statement = [
					{'ANALYSIS_NAME':seqRawData}, #seqRawData = '@SEQ' currently 
					{                
					 'RECOMBINATION_TYPE': {'$in':include_recombination_type}
					}
				]
				
				if include_analyses:#update to include analyses names				
					include_analyses = [a for a in include_analyses if a != seqRawData]
					or_statement[1]['ANALYSIS_NAME']= {'$in':include_analyses} #now or statment becomes {'ANALYSIS_NAME':'seq' or  [analysis_name:$in list, AND RECOMBINATION_TYPE: $in list ]
				
				self.query_results =  db.seqs.find({
					expIdentifier:exp_query,
					'$or':or_statement
				},fields_to_project).sort('SEQ_ID').hint('exp_seq_id_an_rt').limit(limit*total_docs_per_seq)
			else: #DO NOT FILTER BY RECOMBINATION_TYPE                       
				if seqRawData not in include_analyses:
					include_analyses.append(seqRawData)                    
				self.query_results = db.seqs.find({
						expIdentifier:exp_query,
						'ANALYSIS_NAME':{'$in':include_analyses}
				},fields_to_project).sort('SEQ_ID').hint('exp_seq_id_an_rt').limit(limit*total_docs_per_seq)
				
			#create a generator to parse the cursor and group together documents by SEQ_ID
			self.group_documents_by_seq_id(limit,same_document)                  
		
		self.delim_file_headers = schema_fields
		self.logged_output.write('The following values will be projected: {0}\n'.format(str(fields_to_project)))		
		
		return self
		
	#this function will query num_doc_search within each experiment defined for every analysis type possible in that experiment
	#field name => single string defining unique field to use 
	#distinct_fields => reports which fields we want distinct values from. a list of fields
	#if analyses is not blank, then this function will report distinct values from only those listed analyses
	#num_doc_search => the query will limit the results to this many results 
	#num_doc_skip => the query will skip the first X documents (useful for preforming back to back functions. i.e. getting first 100 then next) 
	#filter_by_query => a general MONGODB query dictionary that defines hoDw to first filter results from database before counting unique values 
	#exp_id => PASS in a list of object ids corresponding to a experiments 	
	#metadata_query => in addition to explicity passing in a list of ObjectId's, a user can further filter experiemnts by querying teh exps collection/metadata
		#for example, if metadata_query = {'PROJECT_NAME':'Demo'}, then first the this metadata will be queried to return exp_ids whose PROJECT_NAME contain 'Demo' 	
	def get_distinct_values_from_experiment(self,field_name,exp_id=None,analysis_name=[],recombination_type=[],filter_by_query={},metadata_query={},num_doc_search = 0,num_doc_skip=0,treat_analysis_types_seperate=False):
		self.query_results = None
		self.logged_output.write('The following function was desired "get_distinct_values_from_experiment". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'distinct field name':field_name, 'treat analyses seperately':treat_analysis_types_seperate,'metadata_query':metadata_query})))
		db = self.db_path#db_reader						
	
		#if treat_analysis_types_seperate:
		#	group_by_field = {'ANALYSIS_NAME':'$ANALYSIS_NAME',field_name:'$'+field_name}
		#else:
		group_by_field = {'distinct_field':'$'+field_name}

		compiled_schema_fields = []
		if not(analysis_name):			
			#in this special situation, we do not know which analyses to search by...
			#so we will obtain the list of possible exp_ids using exp_id and metadata_query (normally this is handled in function aggregate_seqs_collection, but we will do it here to get a list of possible experiments
			
			#if exp_id is passed in as empty, then make it self.allowed_experiment_ids
			if exp_id == None: 
				exp_id = self.allowed_experiment_ids
			elif not isinstance(exp_id,list):
				exp_id = [exp_id]
										
			if metadata_query:		
				#query exp_collection
				#result of query will be stored in class name 'self.query_results'			
				self.query_exps_collection(exp_id,q={f:v for f,v in metadata_query.iteritems()})				
				if self.query_results:
					self.exp_metadata = {str(exp['_id']):exp for exp in self.query_results if exp['_id'] in self.allowed_experiment_ids}
					#erase metadata query 							
					exp_id = [vals['_id'] for meta,vals in self.exp_metadata.iteritems()]
					self.query_results = None
				else:
					self.exp_metadata = {}
					#no experiments to searchby. no experiments were found !
					exp_id = []			
										
			#now we have  alist of possible EXP_ID...	
			if expIdentifier in filter_by_query:					
				possible_exp_id = GetIdIntersection([self.allowed_experiment_ids,exp_id],filter_by_query[expIdentifier])
			else:
				possible_exp_id = GetIdIntersection([self.allowed_experiment_ids,exp_id])
			
			#using list of possible EXP_ID's, query for analysis names 
			temp = db.exps.find({'_id':possible_exp_id},{'ANALYSES_COUNT':1})
			list_of_analyses = []
			for each_exp in temp:#go through each query result/exp metadata,
				if 'ANALYSES_COUNT' in each_exp:
					list_of_analyses.extend(each_exp['ANALYSES_COUNT'].keys())#each key in analyses_count corresponds to a new analyses type put in db
			analysis_name = list(set(list_of_analyses))
			if not analysis_name: #there are no analysis types for this experiment, so just exit function
				self.query_results = None
				return self

	

		#performs a distinct query search on EACH analysis type SEPERATELY 
		if treat_analysis_types_seperate:
			analysis_loops = analysis_name #will be used to loop through each analyses_name seperately for performing aggregation 
		else:				
			#instead of each analysis treated seperately, just group all analysises names in one search
			analysis_loops = [analysis_name] 
		
		#first QUERY THE seqs collection, using exp_id, filter_by_query command, and metadatQuery		
		#self.query_seqs_collection(exp_id,analysis_name,recombination_type,query=filter_by_query,metadata_query=metadata_query)
		
		#CREATE A MONGODB AGGREGATION FUNCTION RATHER THAN USE DISTINCT FUNCTION  => had problems using distinct function with skip and limit 
		#add cursor function to self.query_results 
		#second add a SKIP command to skip X documents
		#third LIMIT command on the query
		#fourth add a DISTINCT command at the end 
		chained_results = [] #this will turn into a cursor keeping track of each of the curosrs called to the database
		
		all_schema_fields = []
		for list_of_analyses in analysis_loops: 			
			#BY putting this in a loop, and query multiple times, we ensure that, rather than return the first X results, it will return the first X results for Y analysis names 
			aggregation_pipeline_query=[
				{'$match':filter_by_query},			
				{'$skip':num_doc_skip},
				{'$group':{
					'_id':group_by_field,
					'ANALYSES_NAMES':{'$addToSet':'$ANALYSIS_NAME'}
					}
				},
				{'$project':{
						'ANALYSES_NAMES':1,
						'_id':0,
						field_name:'$_id.distinct_field'
					}
				},
				{'$match':{field_name:{'$ne':None}}} #ignore null values 
			]
			
			if num_doc_search>0:
				aggregation_pipeline_query.insert(2,{'$limit':num_doc_search})																									
			
			#call aggregation function using this pipelin
			self.aggregate_seqs_collection(exp_id,list_of_analyses,aggregation_fxn=aggregation_pipeline_query,metadata_query=metadata_query)
			
			#keep running tally of all field names (should in theory be the same everytime)
			compiled_schema_fields.extend(self.delim_file_headers)
			
			#chain together the results for each query...			
			chained_results = itertools.chain(chained_results,self.query_results)
			#self.query_results = self.query_results.skip(num_doc_skip).limit(num_doc_search)
		
		#all queries set up, now reset se.fquery reulsts to the CHAIN_TYPES_SEQUENCED
		self.delim_file_headers = list(set(compiled_schema_fields))
		self.query_results = chained_results		
		return self

	#this field will group results from a query into distinct groups based on a specific field in the query result/filtered result  
	#group_by_field_name => a list of fields from the dtabase that you want to group by. analysis_name and recombination_type will ALWAYS be included in group function 
	#analysis_name => query only results form this analysis type 
	#recombination_type => query only results from this recombination type 
	#filter_by_query => a general MONGODB query dictionary that defines how to first filter results formdatabase before counting unique values 
	#include_fields => in addition to unique fields combinations (in gruop by fields) alo include the first occurrence of any fields defined in this list 
	#treat_experiments_seperately => DO NOT INCLUDE identical values from SEPERATE EXPERIMENTS as same. Seperate them into differnet groups
	#limit => only count distinct values from this limit . that is LIMIT THE RESULTS FO THE INITIAL FILTER/QUERY 
	#filter_counts => ONLY SHOW RESULTS WHOSE UNIQUE VALUES HAVE A COUNT > THIS VALUE . for example filter_counts = 1=> shows only sequences/unique values found more than one time in query 
	#sort_counts => SORT results by their respective counts
	def count_unique_values_of_field(self,group_by_field_name,exp_id=None, analysis_name=[],recombination_type=[],filter_by_query = {},metadata_query = {},include_fields=[],treat_experiments_seperately = True,limit=0,filter_counts=0,sort_counts=True,allowDiskUse=True):        		
		self.query_results = None
		self.logged_output.write('The following aggregation function was desired "count unique values of a field". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'group_by_these_fields':group_by_field_name, 'seperateByExp':treat_experiments_seperately,'filter_count_results':filter_counts,'metadata_query':metadata_query})))
		db = self.db_path#db_reader						
				
		if not(isinstance(group_by_field_name,list)):
			group_by_field_name = [group_by_field_name]
		
		if not (isinstance(include_fields,list)):
			include_fields = [include_fields]
			
		if not(isinstance(analysis_name,dict)):
			if not(analysis_name):
				analysis_name = {'$ne':seqRawData}
			elif not(isinstance(analysis_name,list)):			
				analysis_name=[analysis_name]
		
		#for i,f in enumerate(field_name):
		#	if not(f.startswith('DATA')):
		#		field_name[i]='DATA.'+f       
		
			
		#this means we want to treat identical values from different experiments as SEPERATE GROUPS 
		#IN OTHER WORDS: do not group together sequences from differnet experiments
		#ADD EXP_ID to group parameters
		if treat_experiments_seperately:
			group_by_field_name.append(expIdentifier)    
		
		field_name = list(set(group_by_field_name))
		
		#IMPORTANT, WE ARE GROUPING RESULTS BY SPECIFIC FIELD NAMES, AND THEN PROJECTING THE FIRST OCCURRENCE OF RESULTS (FOR EXAMPLE WE MIGHT WANT TO PROJECT UNIQUE CDR3AA BUT ONLY PROJECT THE FIRST OCCURRENCE OF ALL ANTIBODY FIELDS (DATA = 1). 
		#IN THIS CASE , THEN CDR3.AA WILL BE INCLUDED IN THE PROJECT OF DATA, SO WE NEED TO EXLCUDE IT FROM BEING PROJECTED SEPERATELY
		#IF WE DONT EXCLUDE THE SEQUENCE, THEN THE AGGREGATION FUNCTION GENERATES AN ERROR 
		included_grouped_field = {}
		for fields in group_by_field_name:			
			included_grouped_field[fields] = 1 #bydefault project the field 
			for projected_fields in include_fields:
				if projected_fields in ['$$ROOT','$$CURRENT']: #we are projecting the entire document (see mongo documentation), ===> CURRENTLY ROOT OR CURRENT DOES NOT SEEM TO WORK PROPERLY
					included_grouped_field[fields] = 0 	#so no need to project this field
				elif projected_fields == fields or fields.startswith(projected_fields+'.'): #we are projecting the node/parent that already includes this field
					included_grouped_field[fields] = 0 	#so no need to project this field
						
		
		#now build teh aggregation pipeline that we will use 
		#note => analysis names, RECOMBINATION_TYPE, metadata_query, and exp_id will be used in the next function aggregate_seqs_collection
		#this just getnerates the aggregation function that we pass into aggregate_seqs_collection
		
		#first define teh 'MATCH'/filter command of the pipeline
		aggregation_pipeline_query = [
			{'$match':filter_by_query}, #return only docs that match this query 
			#GROUP BY DEFINED FIELDNAME
			{'$group':dict(
					{
						#group results by 1) analysis name, 2) recombination type, 3) all other fields defined by field_name
						'_id':dict({
								'ANALYSIS_NAME':'$ANALYSIS_NAME',
								'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
								}.items()+
								{"FIELDNAME_{0}".format(str(i)):"${0}".format(f) for i,f in enumerate(field_name)}.items() #these are all the fields we use for grouping
							),
						'COUNTS':{'$sum':1},#create a new COUNTS field that will add one for every member added into group (basically count num unique)
					}.items()+
					{each_field.replace('.','____'):{'$first':'$'+each_field} for each_field in include_fields}.items() #project the FIRST OCCURRENCE of any field defined by include_fields			
				)
			},
			{'$project':dict({#RENAME THE RESULTING GROUPS, BASICALLY SEPERATE GROUPS FROM THE ID field 
				'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',#this is first key under '_id' formed from group fxn above
				'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',#this is the second key under '_id' fromed from group fxn above
				'COUNTS':1,#include 'COUNTS' field 
				'_id':0#exclude id field
				}.items()+
				#the custom dictionary of all other ekys in 'id'/group field. 
				{"{0}".format(f):"$_id.FIELDNAME_{0}".format(str(i)) for i,f in enumerate(field_name) if included_grouped_field[f] == 1}.items()+ #these are all the fields we project after grouping. take them outside of the 'ID' field, makes results prettier    
				{"{0}".format(f):"${0}".format(f.replace('.','____')) for f in include_fields}.items() #project the ORIGINAL FIELD NAME of the projected field taken from groouping
				)
			}			
		]
		
		if limit>0:
			#add a limit command to second index position in array 
			aggregation_pipeline_query.insert(1,{'$limit':limit})
			    								
		#add a FINAL filter step if the value, 'filter_counts' > 0 
		if filter_counts>0:			
			aggregation_pipeline_query.append({'$match':{'COUNTS':{'$gt':filter_counts}}})
		
		#add a sort stage at the end to sort by most occuring values
		if sort_counts:
			aggregation_pipeline_query.append({'$sort':{'COUNTS':-1}})
					
		#RUN THE QUERY using the custom aggregation function 
		#THIS FUNCION WILL use the exp_id value and the metadata_query value 
		self.aggregate_seqs_collection(exp_id,analysis_name,recombination_type,aggregation_pipeline_query,metadata_query=metadata_query,allowDiskUse=allowDiskUse)
	
		
	def v_allele_usage(self,exp_id=None,analysis_name=[],recombination_type=[],filter_by_query={},metadata_query={},only_consider_first_allele_element=True, treat_experiments_seperately=True,limit=0,sort_counts=True,allowDiskUse=True):
		self.query_results = None
		self.logged_output.write('The following aggregation function was desired "v_allele_usage". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'only_consider_first_allele_element':only_consider_first_allele_element, 'seperateByExp':treat_experiments_seperately,'metadata_query':metadata_query})))
		
		db = self.db_path#db_reader						
									
		if not(isinstance(analysis_name,dict)):
			if not(analysis_name):
				analysis_name = {'$ne':seqRawData}
			elif not(isinstance(analysis_name,list)):			
				analysis_name=[analysis_name]
		
		aggregation_pipeline_query=[
			{'$match':filter_by_query},#filter results by query 
			{'$unwind':'$DATA.VREGION.VGENES'},#unwind all vgene elements in list 			
		]
		
		if only_consider_first_allele_element: #so we want to ignore all other vgenes selected from the unwind command (cannot use slice, to only project first...yet?)
			aggregation_pipeline_query.append(
				{'$group':{#group all 'unwound elements' by '_id. all vgene from list in same documetn should have same '_id'
						'_id':'$_id',
						'V_ALLELE':{'$first':'$DATA.VREGION.VGENES'}, #store the FIRST occurrence of a V GENE 						
						'ANALYSIS_NAME':{'$first':'$ANALYSIS_NAME'},
						'RECOMBINATION_TYPE':{'$first':'$RECOMBINATION_TYPE'},
						 expIdentifier:{'$first':'$'+expIdentifier}
					}
				}
			)
		else:
			#instead of just taking the first element, just rename the v elements 
			aggregation_pipeline_query.append(
				{'$project':{					
						'V_ALLELE':'$DATA.VREGION.VGENES',
						'ANALYSIS_NAME':1,
						'RECOMBINATION_TYPE':1,
						 expIdentifier:1
					}
				}
			)
		
		#NO complette the piipeline by grouping together documents with identical V ALLELE 
		if treat_experiments_seperately:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'V_ALLELE':'$V_ALLELE',
							 expIdentifier:'$'+expIdentifier,
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique vallele/exp occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',
						 expIdentifier:'$_id.'+expIdentifier,
						'DATA.VREGION.VGENES':'$_id.V_ALLELE'						
					}					
				}
			])
		else:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'V_ALLELE':'$V_ALLELE',							
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique vallele occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',					
						'DATA.VREGION.VGENES':'$_id.V_ALLELE'						
					}					
				}
			])
		
		if limit>0:
			#add a limit command to second index position in array 
			aggregation_pipeline_query.insert(1,{'$limit':limit})
		
		#add a sort stage at the end to sort by most occuring values
		if sort_counts:
			aggregation_pipeline_query.append({'$sort':{'COUNTS':-1}})
					
		#RUN THE QUERY using the custom aggregation function 
		#THIS FUNCION WILL use the exp_id value and the metadata_query value 
		self.aggregate_seqs_collection(exp_id,analysis_name,recombination_type,aggregation_pipeline_query,metadata_query=metadata_query,allowDiskUse=allowDiskUse)
		
		return self
	
	def j_allele_usage(self,exp_id=None,analysis_name=[],recombination_type=[],filter_by_query={},metadata_query={},only_consider_first_allele_element=True, treat_experiments_seperately=True,limit=0,sort_counts=True,allowDiskUse=True):
		self.query_results = None
		self.logged_output.write('The following aggregation function was desired "j_allele_usage". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'only_consider_first_allele_element':only_consider_first_allele_element, 'seperateByExp':treat_experiments_seperately,'metadata_query':metadata_query})))
		
		db = self.db_path#db_reader						
									
		if not(isinstance(analysis_name,dict)):
			if not(analysis_name):
				analysis_name = {'$ne':seqRawData}
			elif not(isinstance(analysis_name,list)):			
				analysis_name=[analysis_name]
		
		aggregation_pipeline_query=[
			{'$match':filter_by_query},#filter results by query 
			{'$unwind':'$DATA.JREGION.JGENES'},#unwind all Jgene elements in list 			
		]
		
		if only_consider_first_allele_element: #so we want to ignore all other Jgenes selected from the unwind command (cannot use slice, to only project first...yet?)
			aggregation_pipeline_query.append(
				{'$group':{#group all 'unwound elements' by '_id. all Jgene from list in same documetn should have same '_id'
						'_id':'$_id',
						'J_ALLELE':{'$first':'$DATA.JREGION.JGENES'}, #store the FIRST occurrence of a J GENE 						
						'ANALYSIS_NAME':{'$first':'$ANALYSIS_NAME'},
						'RECOMBINATION_TYPE':{'$first':'$RECOMBINATION_TYPE'},
						 expIdentifier:{'$first':'$'+expIdentifier}
					}
				}
			)
		else:
			#instead of just taking the first element, just rename the J elements 
			aggregation_pipeline_query.append(
				{'$project':{					
						'J_ALLELE':'$DATA.JREGION.JGENES',
						'ANALYSIS_NAME':1,
						'RECOMBINATION_TYPE':1,
						 expIdentifier:1
					}
				}
			)
		
		#NO complette the piipeline by grouping together documents with identical J ALLELE 
		if treat_experiments_seperately:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'J_ALLELE':'$J_ALLELE',
							 expIdentifier:'$'+expIdentifier,
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique Jallele/exp occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',
						 expIdentifier:'$_id.'+expIdentifier,
						'DATA.JREGION.JGENES':'$_id.J_ALLELE'						
					}					
				}
			])
		else:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'J_ALLELE':'$J_ALLELE',							
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique Jallele occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',					
						'DATA.JREGION.JGENES':'$_id.J_ALLELE'						
					}					
				}
			])
		
		if limit>0:
			#add a limit command to second index position in array 
			aggregation_pipeline_query.insert(1,{'$limit':limit})
		
		#add a sort stage at the end to sort by most occuring values
		if sort_counts:
			aggregation_pipeline_query.append({'$sort':{'COUNTS':-1}})
					
		#RUN THE QUERY using the custom aggregation function 
		#THIS FUNCION WILL use the exp_id value and the metadata_query value 
		self.aggregate_seqs_collection(exp_id,analysis_name,recombination_type,aggregation_pipeline_query,metadata_query=metadata_query,allowDiskUse=allowDiskUse)
	
		return self
	
	def vj_allele_usage(self,exp_id=None,analysis_name=[],recombination_type=[],filter_by_query={},metadata_query={},only_consider_first_allele_element=True, treat_experiments_seperately=True,limit=0,sort_counts=True,allowDiskUse=True):
		self.query_results = None
		self.logged_output.write('The following aggregation function was desired "v_allele_usage". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'only_consider_first_allele_element':only_consider_first_allele_element, 'seperateByExp':treat_experiments_seperately,'metadata_query':metadata_query})))
		
		db = self.db_path#db_reader						
									
		if not(isinstance(analysis_name,dict)):
			if not(analysis_name):
				analysis_name = {'$ne':seqRawData}
			elif not(isinstance(analysis_name,list)):			
				analysis_name=[analysis_name]
		
		aggregation_pipeline_query=[
			{'$match':filter_by_query},#filter results by query 
			{'$unwind':'$DATA.VREGION.VGENES'},#unwind all vgene elements in list 			
			{'$unwind':'$DATA.JREGION.JGENES'},#unwind all Jgene elements in list 			
		]
		
		if only_consider_first_allele_element: #so we want to ignore all other vgenes selected from the unwind command (cannot use slice, to only project first...yet?)
			aggregation_pipeline_query.append(
				{'$group':{#group all 'unwound elements' by '_id. all vgene from list in same documetn should have same '_id'
						'_id':'$_id',
						'V_ALLELE':{'$first':'$DATA.VREGION.VGENES'}, #store the FIRST occurrence of a V GENE 						
						'J_ALLELE':{'$first':'$DATA.JREGION.JGENES'}, #store the FIRST occurrence of a J GENE 						
						'ANALYSIS_NAME':{'$first':'$ANALYSIS_NAME'},
						'RECOMBINATION_TYPE':{'$first':'$RECOMBINATION_TYPE'},
						 expIdentifier:{'$first':'$'+expIdentifier}
					}
				}
			)
		else:
			#instead of just taking the first element, just rename the v elements 
			aggregation_pipeline_query.append(
				{'$project':{					
						'V_ALLELE':'$DATA.VREGION.VGENES',
						'J_ALLELE':'$DATA.JREGION.JGENES',
						'ANALYSIS_NAME':1,
						'RECOMBINATION_TYPE':1,
						 expIdentifier:1
					}
				}
			)
		
		#NOW complette the piipeline by grouping together documents with identical V ALLELE 
		if treat_experiments_seperately:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'V_ALLELE':'$V_ALLELE',
							'J_ALLELE':'$J_ALLELE',
							 expIdentifier:'$'+expIdentifier,
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique vallele/exp occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',
						 expIdentifier:'$_id.'+expIdentifier,
						'DATA.VREGION.VGENES':'$_id.V_ALLELE',						
						'DATA.JREGION.JGENES':'$_id.J_ALLELE'						
					}					
				}
			])
		else:
			aggregation_pipeline_query.extend([
				{'$group':{			
						'_id':{
							'V_ALLELE':'$V_ALLELE',							
							'J_ALLELE':'$J_ALLELE',							
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
						},
						'COUNTS':{'$sum':1}#count each unique vallele occurrence 
					}
				},
				{'$project':{
						'_id':0,
						'COUNTS':1,
						'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',
						'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',					
						'DATA.VREGION.VGENES':'$_id.V_ALLELE',
						'DATA.JREGION.JGENES':'$_id.J_ALLELE'						
					}					
				}
			])
		
		if limit>0:
			#add a limit command to second index position in array 
			aggregation_pipeline_query.insert(1,{'$limit':limit})
		
		#add a sort stage at the end to sort by most occuring values
		if sort_counts:
			aggregation_pipeline_query.append({'$sort':{'COUNTS':-1}})
					
		#RUN THE QUERY using the custom aggregation function 
		#THIS FUNCION WILL use the exp_id value and the metadata_query value 
		self.aggregate_seqs_collection(exp_id,analysis_name,recombination_type,aggregation_pipeline_query,metadata_query=metadata_query,allowDiskUse=allowDiskUse)
	
		return self
	
		
	#this field will group results from a query into by the field CDR3_AA_LENGTH	
	#feature=> allowed values 'NT', or 'AA' => whether or not to perform length distribution on cdr3 aa or cdr3 nt 
	#analysis_name => query only results form this analysis type 
	#recombination_type => query only results from this recombination type 
	#filter_by_query => a general MONGODB query dictionary that defines how to first filter results formdatabase before counting unique values 	
	#treat_experiments_seperately => DO NOT INCLUDE identical values from SEPERATE EXPERIMENTS as same. Seperate them into differnet groups
	#limit => only count distinct values from this limit . that is LIMIT THE RESULTS FO THE INITIAL FILTER/QUERY 	
	#sort_counts => SORT CDR3_LENGTHS by their respective counts	
	def cdr3_length_distribution(self,exp_id=None,feature='AA',analysis_name = [],recombination_type=[],filter_by_query={},metadata_query = {},treat_experiments_seperately = True,limit=0,sort_counts=True,exclude_cdr3_lengths_not_multiple_three=False, exclude_empty_cdr3_strings=True,allowDiskUse=True):
		self.query_results = None
		feature = feature.upper()
		if feature=='NT':
			field_name = ['DATA.CDR3.NT_LENGTH']
			if exclude_empty_cdr3_strings:#only select cdr3 lengths > 2 (only length of 3 or more is releveant AA) 
				if exclude_cdr3_lengths_not_multiple_three:
					#also only choose values whose cdr3 nt length is a factor of 3 
					filter_by_query['DATA.CDR3.NT_LENGTH'] = {'$gt':2,'$mod':[3,0]}
				else:
					filter_by_query['DATA.CDR3.NT_LENGTH'] = {'$gt':2}
		elif feature=='AA':
			field_name = ['DATA.CDR3.AA_LENGTH']
			if exclude_empty_cdr3_strings:#only select cdr3 lengths > 0 
				filter_by_query['DATA.CDR3.AA_LENGTH'] = {'$gt':0} 
				if exclude_cdr3_lengths_not_multiple_three:				
					#also only choose values whose cdr3 nt length is a factor of 3 
					filter_by_query['DATA.CDR3.NT_LENGTH'] = {'$mod':[3,0]}
			
		else:
			raise Exception("The variable feature can only be a string with the value 'NT' or 'AA'")
		
		self.logged_output.write('The following aggregation function was desired "cdr3_length_distribution". The following settings were requested: {0}\n'.format(str({'experiments':exp_id,'cdr3_feature':feature, 'seperateByExp':treat_experiments_seperately,'metadata_query':metadata_query})))
		db = self.db_path#db_reader						
				
					
		if not(isinstance(analysis_name,dict)):
			if not(analysis_name):
				analysis_name = {'$ne':seqRawData}
			elif not(isinstance(analysis_name,list)):			
				analysis_name=[analysis_name]
		
							
		#this means we want to treat identical values from different experiments as SEPERATE GROUPS 
		#IN OTHER WORDS: do not group together sequences from differnet experiments
		#ADD EXP_ID to group parameters		
		if treat_experiments_seperately:
			field_name.append(expIdentifier)    
		
		#now build teh aggregation pipeline that we will use 
		#note => analysis names, RECOMBINATION_TYPE, metadata_query, and exp_id will be used in the next function aggregate_seqs_collection
		#this just getnerates the aggregation function that we pass into aggregate_seqs_collection
		
		#first define teh 'MATCH'/filter command of the pipeline
		aggregation_pipeline_query = [
			{'$match':filter_by_query}, #return only docs that match this query 
			#GROUP BY DEFINED FIELDNAME
			{'$group':{
					#group results by 1) analysis name, 2) recombination type, 3) all other fields defined by field_name
					'_id':dict({
							'ANALYSIS_NAME':'$ANALYSIS_NAME',
							'RECOMBINATION_TYPE':'$RECOMBINATION_TYPE'
							}.items()+
							{"FIELDNAME_{0}".format(str(i)):"${0}".format(f) for i,f in enumerate(field_name)}.items() #these are all the fields we use for grouping
						),
					'COUNTS':{'$sum':1},#create a new COUNTS field that will add one for every member added into group (basically count num unique)
				
				}
			},
			{'$project':dict({#RENAME THE RESULTING GROUPS, BASICALLY SEPERATE GROUPS FROM THE ID field 
					'ANALYSIS_NAME':'$_id.ANALYSIS_NAME',#this is first key under '_id' formed from group fxn above
					'RECOMBINATION_TYPE':'$_id.RECOMBINATION_TYPE',#this is the second key under '_id' fromed from group fxn above
					'COUNTS':1,#include 'COUNTS' field 
					'_id':0#exclude id field
					}.items()+
					#the custom dictionary of all other ekys in 'id'/group field. 
					{"{0}".format(f):"$_id.FIELDNAME_{0}".format(str(i)) for i,f in enumerate(field_name)}.items() #these are all the fields we project after grouping. take them outside of the 'ID' field, makes results prettier    				
				)
			}
		]
		
		if limit>0:
			#add a limit command to second index position in array 
			aggregation_pipeline_query.insert(1,{'$limit':limit})
		
		#add a sort stage at the end to sort by most occuring values
		if sort_counts:
			aggregation_pipeline_query.append({'$sort':{'COUNTS':-1}})
					
		#RUN THE QUERY using the custom aggregation function 
		#THIS FUNCION WILL use the exp_id value and the metadata_query value 
		self.aggregate_seqs_collection(exp_id,analysis_name,recombination_type,aggregation_pipeline_query,metadata_query=metadata_query,allowDiskUse=allowDiskUse)
	
		return self
			
	
	#################Cursor Modification Functions######################
	#these functions can be used for modifying the result from the database. All functions will be generator functions
	
	
	#convert the sequences to a delimited string format for output 
	#IF to_file => false, then we ignore this function since the user will be reading the results to memory as json variable any. 
	#BUT IF TO_FILE=> TRUE, then we need to create a file delimeted by the char defined in 'delimiter'
	#header_var => the user can pass in the preferred order/output for fields: 
		#if its a list => then the user would like to request the fields defined in the list 
		#if header_var is a dictionary => then user would like to request the fields DEFINED BY KEYS BUT rename the fields by the value 
			#i.e. {field_name_in_db:column_name_in_file}
		#if header var is a tuple => then user would like to request the fields DEFINED BY INDEX POSITION 0 BUT rename the fields by INDEX POSITION 1 
	#split_results_by => split documents by the query into multiple files using the fields defined here
	#keep_all_info => if TRUE, then save all fields from output doucment. 
					 #IF false then ONLY save values defined by header_var
	#save_by_filename => split results based on the filename inputed into experiments (filename field in database)
	#save_by_exp_name => split results/documents based on respective experiment name 
	#avoid_commas_in_string_output => if TRUE, then we wil not use ',' to split lists. 
	def convert_to_delim(self,delimiter,header_var = [], split_results_by=[],keep_all_info=False,save_by_filename=False,save_by_exp_name=False,avoid_commas_in_string_output=False):				
		delim_len = len(delimiter)
		def generate_delim(query_results,output_sorted_fields):			
			keys_to_be_reported = output_sorted_fields.keys()
			h_l = len(keys_to_be_reported)-1
			range_count = range(h_l)
			if avoid_commas_in_string_output:
				schema_schema_output_fnc_for_delim = schema_fields_to_file_avoid_commas
			else:
				schema_schema_output_fnc_for_delim = schema_fields_to_file
				
			if self.to_file:
				yield 'HeaderDescription'+file_def_seperator+(delimiter).join(output_sorted_fields.values())+'\n'				
			for seq_document in query_results:
				seq_document =	defaultdict(str,Process_Cursor_For_Output_File(seq_document,schema_schema_output_fnc_for_delim))				
				if self.to_file:	
					prefix_line = self.file_prefix+'.' if self.file_prefix else ''
					
					if save_by_exp_name and expIdentifier in seq_document:						
						prefix_line+=self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME']+'.'
					
					if save_by_filename and 'FILENAME' in seq_document:
						prefix_line+=seq_document['FILENAME']+'.'
											
					for f in split_results_by:						
						prefix_line+=seq_document[f]+'.' if f in seq_document else '.'						
					
					prefix_line = prefix_line[:-1]
															
					
					if prefix_line:
						#prefix_line+=file_def_seperator					
						yield prefix_line+file_def_seperator						
					yield (delimiter).join([seq_document[k] for k in keys_to_be_reported])+'\n'
				else:										
					yield json.dumps({output_field:seq_document[db_field] for db_field,output_field in output_sorted_fields.iteritems()})+'\n'																				
		
		output_sorted_fields = OrderedDict()					
		
		##LETS TRY AND GET THE ORDER OF THE COLUMNS CORRECT 
		#ACOUNT FOR USER REQUESTS, DEFAULT SORT ORDERS, AND FOUND FIELDS FROM QUERIES#
		user_requested_fields = True
		if header_var == [] or not(header_var):	
			#user does not care about order of results, so we  will use the sorting order , we will also force keep_all_info to be true 
			header_var = default_sorting_order
			keep_all_info = True							
			user_requested_fields = False
			header_type = 'list'
		elif isinstance(header_var,basestring): # make a list
			header_var = [header_var]	
			header_type = 'list'
		elif isinstance(header_var,dict):
			header_type = 'dict'
		elif isinstance(header_var,list):
			first_val = type(header_var[0])
			if isinstance(header_var[0],basestring):
				header_type='list'				
			elif isinstance(header_var[0],type(('db_field','new_field'))):
				header_type='tuple'
			else:
				raise Exception('The parameter, header_var, can only be a list of strings or tuples.')							
			for all_vals in header_var:
				if not(isinstance(all_vals,first_val)):
					raise Exception('The parameter, header_var, can only be a list of strings or list of tuples. It can not have multiple types within the list')
		else:			
			raise Exception('The parameter, header_var, can only be a string, list, tuples or dictionary of strings')		
		
													
		
		self.delim_file_headers = [s.upper() if s!='_id' else s for s in self.delim_file_headers]
		
		append_id = True if idIdentifier in self.delim_file_headers else False
			
		
		
		if header_type=='list':	
			header_var = [h.upper() for h in header_var]
			if append_id and idIdentifier not in header_var:
				header_var.insert(0,idIdentifier)
			#user passed in a list defining the order to output fields.
			#BECAUSE its a list, the NAMES of each field will not change 
			
			#THE USER EXPLICITELY SAID THEY WANTED THESE FIELDS, SO THIS IS HOW HTE COLUMNS SHALL BE 
			if user_requested_fields:
				#so all fields in header_var will be output to file 
				for fields_requested in header_var:
					output_sorted_fields[fields_requested] = fields_requested											
				
			else:
				#if the user does not explicitely define which fields they want, then we will only output fields 
				#found in the query (self.delim_file_headers)
				for fields_requested in header_var:
					if fields_requested in self.delim_file_headers:
						output_sorted_fields[fields_requested] = fields_requested											
				
		elif header_type=='dict':
			header_var = {h.upper():nv if h!='_id' else h  for h,nv in header_var.iteritems()}
			if idIdentifier not in header_var and append_id==True:
				header_var[idIdentifier] = idIdentifier
			
			#cannot control order of fields, so use the order defined in default_sorting_order
			sub_dict = {}
			dic_len = 0
			max_pos = 0
			new_requested_fields = []
			for each_requested_key in header_var:								
				#OK now we have a default sorting order of the keys 
				if each_requested_key in default_sorting_order:
					sub_dict[each_requested_key] = default_sorting_order.index(each_requested_key)										
				else:
					new_requested_fields.append(each_requested_key)
					
			max_pos = max(sub_dict.values())+1 if sub_dict else 0 
			#any keys that were not found in our defualt variable get appende to dictionary 
			for each_requested_key in new_requested_fields:				
				sub_dict[each_requested_key] = max_pos
				max_pos+=1
			
			#add in the keys sorted by their psostiion to the ordered dictionary 
			for key, value in sorted(sub_dict.iteritems(), key=lambda (k,v): (v,k)):
				output_sorted_fields[key.upper()] = header_var[key]
				#key = field name 
				#value = column name in output file 							
																													
		elif header_type=='tuple':
			header_var = [(row[0].upper(),row[1]) if row[0]!='_id' else row[0] for row in header_var]
			
			if append_id:
				found = False
				for i in header_var:
					if i[0] == idIdentifier:
						found = True
						break
				if found == False:
					header_var.insert(0,(idIdentifier,idIdentifier))
				
			#user passed in a tuple defining both the order to output fields AND the new field names.			
			for fields_requested in header_var:
				output_sorted_fields[fields_requested[0]] = fields_requested[1]			
				
		
			
			
		#NOW go through the fields produced by the query 
		#IF they are not currently found in output_sorted_fields, then add them to the fields to be output 
		if keep_all_info:
			sub_dict = {}
			new_requested_fields = []
			#loop through fields generated by query 
			for generated_fields in self.delim_file_headers:
				#this field is not currently defined in output_sorted_fields
				if generated_fields not in output_sorted_fields:
					if generated_fields in default_sorting_order:
						sub_dict[generated_fields] = default_sorting_order.index(generated_fields)
					else:
						new_requested_fields.append(generated_fields)
			
			max_pos = max(sub_dict.values())+1 if sub_dict else 0
			
			#any keys that were not found yet gets appended to dictionary 
			for each_requested_key in new_requested_fields:				
				sub_dict[each_requested_key] = max_pos
				max_pos+=1
			
			#add in the keys sorted by their psostiion to the ordered dictionary 
			for key, value in sorted(sub_dict.iteritems(), key=lambda (k,v): (v,k)):
				output_sorted_fields[key] = key		
								
		########ORDER OF OUTPUT FILES DETERMINED###
		
					
		if self.file_prefix == '' and split_results_by==[] and save_by_filename==False and save_by_exp_name==False:
			self.file_prefix = 'IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')																
		
		if isinstance(self.query_results,dict):
			self.query_results = [self.query_results]		
				
		self.query_results = generate_delim(self.query_results,output_sorted_fields)
	
	
	#convert the sequences to a TAB file format
	#IF to_file => false, then we ignore this function since the user will be reading the results to memory as json variable any. 
	#BUT IF TO_FILE=> TRUE, then we need to create a file delimeted by the char defined in 'delimiter'
	#header_var => the user can pass in the preferred order/output for fields: 
		#if its a list => then the user would like to request the fields defined in the list 
		#if header_var is a dictionary => then user would like to request the fields DEFINED BY KEYS BUT rename the fields by the value 
			#i.e. {field_name_in_db:column_name_in_file}
	#split_results_by => split documents by the query into multiple files using the fields defined here
	#keep_all_info => if TRUE, then save all fields from output doucment. 
					 #IF false then ONLY save values defined by header_var
	#save_by_filename => split results based on the filename inputed into experiments (filename field in database)
	#save_by_exp_name => split results/documents based on respective experiment name 
	def convert_to_tab(self,header_var = [], split_results_by=[],keep_all_info=False,save_by_filename=False,save_by_exp_name=False):
		self.convert_to_delim('\t',header_var, split_results_by,keep_all_info,save_by_filename,save_by_exp_name)
	
	
	#convert the sequences to a CSV file format
	#IF to_file => false, then we ignore this function since the user will be reading the results to memory as json variable any. 
	#BUT IF TO_FILE=> TRUE, then we need to create a file delimeted by the char defined in 'delimiter'
	#header_var => the user can pass in the preferred order/output for fields: 
		#if its a list => then the user would like to request the fields defined in the list 
		#if header_var is a dictionary => then user would like to request the fields DEFINED BY KEYS BUT rename the fields by the value 
			#i.e. {field_name_in_db:column_name_in_file}
	#split_results_by => split documents by the query into multiple files using the fields defined here
	#keep_all_info => if TRUE, then save all fields from output doucment. 
					 #IF false then ONLY save values defined by header_var
	#save_by_filename => split results based on the filename inputed into experiments (filename field in database)
	#save_by_exp_name => split results/documents based on respective experiment name 
	def convert_to_csv(self,header_var = [], split_results_by=[],keep_all_info=False,save_by_filename=False,save_by_exp_name=False):
		self.convert_to_delim(',',header_var, split_results_by,keep_all_info,save_by_filename,save_by_exp_name,avoid_commas_in_string_output=True)
			
	
	#convert query results into IGREP format
	#flatten all results from database 
	#convert all values in document to string format
	def convert_to_igrep(self,split_results_by=[],save_by_filename=False,save_by_exp_name=False):							
		def generate_igrep(results):			
			for seq_document in results:				
				seq_doc = Process_Cursor_For_Output_File(seq_document)				
				if self.to_file:										
					prefix_line = self.file_prefix if self.file_prefix else ''									
					
					if save_by_exp_name and expIdentifier in seq_doc:						
						prefix_line+=self.exp_metadata[seq_doc[expIdentifier]]['EXPERIMENT_NAME']+'.'
					
					if save_by_filename and 'FILENAME' in seq_doc:
						prefix_line+=seq_doc['FILENAME']+'.'
											
					for f in split_results_by:						
						prefix_line+=seq_doc[f]+'.' if f in seq_doc else '.'						
					
					prefix_line = prefix_line[:-1]
					
					if prefix_line:
						#prefix_line+=file_def_seperator					
						yield prefix_line+file_def_seperator						
				
				#else:
				#	prefix_line = ''												
				yield json.dumps(seq_doc)+'\n' #=> DO NOT YIELD PREFIX_LINE+JSON.DUMPS...THAT ADDS UNNECESSARY TIME!!					
				#yield ''
		self.query_results = generate_igrep(self.query_results)
	
	#convert query results in to FASTA format
	#IMPORTANT:
		#=> IF to_file defined above (see _retrun and _return_results) is True, the result is returned in fasta format >header\ndata\n because the user wants the result to file 
		#=> BUT if to_file = False, then the user will want to use the result in memeory. so this will return a dictionary of the following {'header':,'sequence':}
	#if keep_all_info: Then dump any remaining fields in seq document that are not accounted for in sequence header 		
	#save_by_exp_name => if TRUE, then the experiment name is added to the output filename  
	#save_by_filename => if TRUE AND IF Filename field is present in document, then the value for 'filename' is added to the output filename	
	#split_results_by => ONLY valid if self.to_file = TRUE. PASS in list of fields that you want to use to group together documents into multiple files. 
		#=> split_results_by = ['EXP_ID']  => all documents containing the same value for 'EXP_ID' will be grouped into a single same file
	def convert_to_fasta(self,sequence_key = 'SEQUENCE', header_var = ['SEQUENCE_HEADER'],keep_all_info=False,save_by_exp_name=False,save_by_filename=False,split_results_by = []):							
		def generate_fasta(results,to_file,split_results_by):						
			#if len(header_var)>1 and to_file ==True:
			#	delimHeader=True
			#if to_file:
			#	yield 'HeaderDescription'+file_def_seperator+'#'+'|'.join(header_var)+'\n'
			for seq_document in results:
				#process document:
					#flatten dictionary 
					#convert all fields to strings 
				seq_document =	Process_Cursor_For_Output_File(seq_document)
				if sequence_key in seq_document:				
					if to_file:	
						if self.file_prefix:
							prefix_line = [self.file_prefix]
						else:
							prefix_line = []
						
						
						if save_by_exp_name and expIdentifier in seq_document:						
							prefix_line.append(self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'])
						
						if save_by_filename and 'FILENAME' in seq_document:
							prefix_line.append(seq_document['FILENAME'])
							
						for f in split_results_by:
							if f in seq_document:
								prefix_line.append(seq_document[f])
							else:
								prefix_line.append('')
							
						prefix_line = '.'.join(prefix_line)
						
					else:
						prefix_line = ''
					
					
					if idIdentifier in seq_document:
						seq_id = seq_document.pop(idIdentifier)
						added_data = {idIdentifier:seq_id}
					else:					
						added_data = {}	
						
					
					#header = seq_document.pop(header_key) if header_key in seq_document else seq_id
					sequence = seq_document.pop(sequence_key)																		
					#generate sequence header					
					header_str = seq_document.pop(header_var[0]) if header_var[0] in seq_document else ''
					
					for ind in range(1,len(header_var)):
						header_str += '|'+seq_document.pop(header_var[ind]) if header_var[ind] in seq_document else '|'															
					
					if keep_all_info:							
						added_data.update(seq_document)
										
					descriptor_str=header_str+fasta_file_delimiter+json.dumps(added_data) if added_data else header_str
										
					if to_file:
						if prefix_line:
							prefix_line+=file_def_seperator
							yield '{0}>{1}\n{0}{2}\n'.format(prefix_line,descriptor_str,sequence)
						else:
							yield '>{1}\n{2}\n'.format(prefix_line,descriptor_str,sequence)
					else:
						if added_data:		
							added_data.update({'sequence':sequence,'header':header_str,'document_header':descriptor_str})
							yield json.dumps(added_data)+'\n'
						else:
							yield json.dumps({'sequence':sequence,'header':header_str,'document_header':descriptor_str})+'\n'
		
		if not(header_var):
			keep_all_info = True	
																				
		if not isinstance(header_var,list):
			header_var = [header_var]
		#header_var = [s.upper() for s in header_var]
		#make header_var unique 
		temp_header_var = header_var 
		header_var = []
		for each_field in temp_header_var:
			if each_field.startswith('DATA.'):
				each_field = each_field[5:]
			if each_field not in header_var:				
				header_var.append(each_field.upper())
		
		#always remvoe sequence key from header_var 
		if sequence_key in header_var:
			header_var.remove(sequence_key)
		
		
		if self.file_prefix == '' and split_results_by==[] and save_by_filename==False and save_by_exp_name==False:
			self.file_prefix = 'IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')																
		
		if isinstance(self.query_results,dict):
			self.query_results = [self.query_results]		
		

		self.query_results = generate_fasta(self.query_results,self.to_file,split_results_by)							

	
	
	#convert query results in to FASTA format
	#IMPORTANT:
		#=> IF to_file defined above (see _retrun and _return_results) is True, the result is returned in fasta format >header\ndata\n because the user wants the result to file 
		#=> BUT if to_file = False, then the user will want to use the result in memeory. so this will return a dictionary of the following {'header':,'sequence':}
	#if keep_all_info: Then dump any remaining fields in seq document that are not accounted for in sequence header 	
	#split_results_by => ONLY valid if self.to_file = TRUE. PASS in list of fields that you want to use to group together documents into multiple files. 
		#=> split_results_by = ['EXP_ID']  => all documents containing the same value for 'EXP_ID' will be grouped into a single same file
	def convert_to_fastq(self,sequence_key = 'SEQUENCE', header_var = ['SEQUENCE_HEADER'],keep_all_info=False,save_by_exp_name=False,save_by_filename=False,split_results_by = [],null_quality = 40):							
		min_valid_quality_scores = 33
		max_valid_quality_score = 126
		default_quality_score = 40
		if null_quality:
			if type(null_quality) is int:
				if null_quality<min_valid_quality_scores or null_quality>max_valid_quality_score:
					raise Exception('Default quality score must be between: {0} and {1} '.format(str(min_valid_quality_scores,max_valid_quality_score)))
				default_char = chr(null_quality)
			elif (type(null_quality) is str or type(null_quality) is unicode) and len(null_quality)==1:
				ascii = ord(null_quality)
				if ascii<min_valid_quality_scores or ascii>max_valid_quality_score:
					raise Exception('Default quality score must be between: {0} and {1}. Parameter passed was an ascii value of {2}'.format(str(min_valid_quality_scores),str(max_valid_quality_score),str(ascii)))
				default_char = null_quality
			else:
				raise Exception('Invalid ascii value passed to null_quality parameter')
		else:
			default_char = chr(default_quality_score)
		
						
		def generate_fastq(results,to_file,split_results_by):					
			default_char_str = ''.join([default_char]*500)
			char_str_len = len(default_char_str)
			#if len(header_var)>1 and to_file ==True:
			#	delimHeader=True
			#if to_file:
			#	yield 'HeaderDescription'+file_def_seperator+'#'+'|'.join(header_var)+'\n'
			for seq_document in results:
				#process document:
					#flatten dictionary 
					#convert all fields to strings 
				seq_document =	Process_Cursor_For_Output_File(seq_document)								
				if sequence_key in seq_document:				

					if to_file:	
						if self.file_prefix:
							prefix_line = [self.file_prefix]
						else:
							prefix_line = []										
						if save_by_exp_name and expIdentifier in seq_document:						
							prefix_line.append(self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'])					
						if save_by_filename and 'FILENAME' in seq_document:
							prefix_line.append(seq_document['FILENAME'])						
						for f in split_results_by:
							if f in seq_document:
								prefix_line.append(seq_document[f])
							else:
								prefix_line.append('')
							
						prefix_line = '.'.join(prefix_line)
						
					else:
						prefix_line = ''
					
					sequence = seq_document.pop(sequence_key)																		
					if 'QUALITY_SCORE' in seq_document:
						quality_score = seq_document['QUALITY_SCORE']
					else:					
						if len(sequence)>char_str_len:
							default_char_str = ''.join([default_char]*len(sequence))
							char_str_len = len(sequence)
							quality_score = default_char_str
						else:
							quality_score = default_char_str[:len(sequence)]
					
					if idIdentifier in seq_document:
						seq_id = seq_document.pop(idIdentifier)
						added_data = {idIdentifier:seq_id}
					else:					
						added_data = {}	
						
					
					
					
					#generate sequence header					
					header_str = seq_document.pop(header_var[0]) if header_var[0] in seq_document else ''
					
					for ind in range(1,len(header_var)):
						header_str += '|'+seq_document.pop(header_var[ind]) if header_var[ind] in seq_document else '|'															
					
					if keep_all_info:							
						added_data.update(seq_document)
										
					descriptor_str=header_str+fasta_file_delimiter+json.dumps(added_data) if added_data else header_str
										
					if to_file:
						if prefix_line:
							prefix_line+=file_def_seperator
							yield '{0}@{1}\n{0}{2}\n{0}+\n{0}{3}\n'.format(prefix_line,descriptor_str,sequence,quality_score)
						else:
							yield '@{1}\n{2}\n+\n{3}\n'.format(prefix_line,descriptor_str,sequence,quality_score)
					else:
						if added_data:		
							added_data.update({'sequence':sequence,'header':header_str,'document_header':descriptor_str,'phred':quality_score})
							yield json.dumps(added_data)+'\n'
						else:
							yield json.dumps({'sequence':sequence,'header':header_str,'document_header':descriptor_str,'phred':quality_score})+'\n'												
		
		if not header_var:
			keep_all_info = True
																				
		#header_var = [s.upper() for s in header_var]
		#make header_var unique 
		temp_header_var = header_var 
		header_var = []
		for each_field in temp_header_var:
			if each_field.startswith('DATA.'):
				each_field = each_field[5:]
			if each_field not in header_var:				
				header_var.append(each_field.upper())
		
		#always remvoe sequence key from header_var 
		if sequence_key in header_var:
			header_var.remove(sequence_key)
		
		#always remove quality score 
		if 'QUALITY_SCORE' in header_var:
			header_var.remove('QUALITY_SCORE')
		
		if self.to_file and self.file_prefix == '' and split_results_by==[] and save_by_filename==False and save_by_exp_name==False:
			self.file_prefix = 'IGREP_Query_'+str(datetime.now()).replace(':','').replace('_','').replace('-','').replace(' ','').replace('.','')																
		
		if isinstance(self.query_results,dict):
			self.query_results = [self.query_results]		
		

		self.query_results = generate_fastq(self.query_results,self.to_file,split_results_by)							

	#creates a generator for going through the cursor results and grouping results by seq_id
	#if same_document, then all results are merged into one document, if False, then only the information from @SEQ is appended to each dcoument
	def group_documents_by_seq_id(self,limit=None ,same_document=False):   		
		def group_by_id(results):				
			reported_docs = 0
			seq_id = ''
			output_doc = {}
			
			for doc_results in results:				
				if doc_results[idIdentifier]!=seq_id:					
					#new doc found
					if output_doc:
						if same_document: #output all results to the same result/document	/same row 						
							yield output_doc
						else: #split up results by annotation 
							sub_doc = output_doc.pop('ANALYSES',None)							
							if sub_doc: #NO annotation results found
								output_doc.pop('ANALYSIS_NAME',None)
								output_doc.pop('DATE_UPDATED',None)
								output_doc.pop('_id',None)
								for analysis,each_doc in sub_doc.iteritems():
									each_doc.update(output_doc)								
									yield each_doc
							else:
								yield output_doc						
						reported_docs+=1
					seq_id = doc_results[idIdentifier]				
					output_doc = {'ANALYSES':{}}#initialize outputdoc							
					
					if limit and reported_docs==limit:
						break
				data = doc_results.pop('DATA',None)	#remove data key	
				doc_results.update(data)
				if doc_results['ANALYSIS_NAME'] == seqRawData:					
					output_doc.update(doc_results) #ALL OF DOC INFORMATION FROM SEQUENCES/RAW DATA WILL APPEAR IN THIS DICTIONARY
				else:
					doc_results.update(data)
					output_doc['ANALYSES'][doc_results['ANALYSIS_NAME']] = doc_results
					
			
			#OUTPUT the last result if there is still information:
			if output_doc and output_doc!={'ANALYSES':{}}:
				if same_document: #output all results to the same result/document							
					yield output_doc
				else: #split up results by annotation 
					sub_doc = output_doc.pop('ANALYSES',None)
					if sub_doc: #NO annotation results found 
						for analysis,each_doc in sub_doc.iteritems():
							each_doc.update(output_doc)								
							yield each_doc
					else:
						yield output_doc								

			
					#output_doc['ANALYSES'][doc_results['ANALYSIS_NAME']

			
		self.query_results = group_by_id(self.query_results)


	#this generator function will return teh RAW sequence data information for every doocument returned from the query 
	def get_rawseq_info(self):
		def get_data(query_data):
			for doc in query_data:
				if 'ANALYSIS_NAME' in doc and doc['ANALYSIS_NAME']!=seqRawData and idIdentifier in doc:
					raw_data = self.db_path.seqs.find({'_id':doc[idIdentifier]},{'DATA':1} ).next()#{idIdentifier:doc[idIdentifier],'ANALYSIS_NAME':seqRawData},{'DATA':1})
					if 'DATA' in raw_data:
						doc.update(raw_data['DATA'])
						yield doc		
					else:
						yield doc
				else:
					yield doc
		
		self.query_results = get_data(self.query_results)
		return self

	#this generator function will take a list of docs and yield them one by one 
	def list_to_gen(self):
		def make_gen(query_data):
			for doc in query_data:				
				yield doc		
		self.query_results = make_gen(self.query_results)		
		return self	
	
def _modify_function(f):
	"""A decorator function that replaces the function it wraps with a function that captures all information necessary to call the function again:
	the function name, its module, and the arguments it was called with. It then sends a POST request to the Immunogrep proxy server with
	the information, passed as a dictionary.
	"""
	@wraps(f)  # This decorator from functools takes care of keeping __name__ and similar attributes of the function f constant, despite being replaced by _post_args
	def _store_calls(*args, **kwargs):
		"""This is the function that the module functions are replaced with.
		"""
		command_dict = {}
		command_dict['command'] = f.__name__
		command_dict['args'] = args[1:]
		command_dict['kwargs'] = kwargs
		#print command_dict
		args[0].location = 'appsoma'
		args[0].query_command_list.append(command_dict)
		return args[0]
	return _store_calls
	
