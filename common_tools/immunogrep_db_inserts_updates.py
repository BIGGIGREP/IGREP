from bennigoetz_mongo_functions import DotAccessible
from bson.objectid import ObjectId
import pymongo
from collections import defaultdict,namedtuple, Counter, deque,OrderedDict
import time
from immunogrep_database_schema import schema_fields_and_data_types
import appsoma_api
from bson.json_util import dumps as bson_dumps
from immunogrep_global_variables import seqRawData
from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import expIdentifier
import immunogrep_database_schema as schema
import immunogrep_useful_functions as useful 
import copy

from cchrysostomou_igdbtools_proxy import connectToIgDatabase as connect_to_ig_database
from immunogrep_database_schema import convert_to_objectid
import subprocess
#for multiprocessing 
from multiprocessing import Process
import os
from pprint import pprint

multiprocessing_threads = 1
#to update:
	#convert settings to list rather than number
	#move all 'filenames' as a number. move actual filename into experiment metadata 
	
	#TO DO IN FUTURE
	#convert locus to list
	#chang name of query data to a shorter name 
	#rename PREDICTED_AB_SEQ to just AB_SEQ

#function for storing datatypes of fields in schema
#this is desired to ensure that all querable fields will be 
#accounted for and we will know HOW to query for the fields
def return_data_types(x):
	if isinstance(x,list):
		string_type='list_of_'
		x=x[0]
	else:
		string_type=''
	
	if isinstance(x,basestring):
		return string_type+'string'
	elif isinstance(x,(long,float,complex,int)):
		return string_type+'number'	
	elif isinstance(x,bool):
		return string_type+'boolean'
	elif isinstance(x,dict):
		return string_type+'dict'	
	else:
		return string_type+str(type(x))


#################WARNING FIX THE CONNECT TO IGDB CALLS IN EVERY FUNCTION WHNE TRANSFERING TO PROXY!!#######################

import json

#OLD INSERT/UPDATE COMMANDS ARE FOUND IN VERSION 493

def delete_entire_experiment(experiment_id_list):
	global db
	global username
	global temp_user_info
	
	[db,connection_data] = connect_to_ig_database()
	
	try:
		#import appsoma_api
		username = appsoma_api.environment_get_username()
		temp_user_info = defaultdict(str,db.users.find_one({'user':username}))
	except:
		pass

	if not type(experiment_id_list) is list:
		experiment_id_list = [convert_to_objectid(experiment_id_list)]
	else:
		experiment_id_list = [convert_to_objectid(exp) for exp in experiment_id_list]			
	if temp_user_info['administrator']!=True:
		exps_to_delete =[result['_id'] for result in db.exps.find({'OWNERS_OF_EXPERIMENT':temp_user_info['user'],'_id':{'$in':experiment_id_list}},{'_id':1})]	
	else:
		exps_to_delete = experiment_id_list		
	if exps_to_delete == []:
		raise Exception("User {0} does not have access to the listed experiments: {1}".format(temp_user_info['user'],str(experiment_id_list)))	
	operation_result = db.seqs.remove({expIdentifier:{'$in':exps_to_delete}})	
	db.exps.remove({'_id':{'$in':exps_to_delete}})
	num_deleted = operation_result['n']
	print "Deleted {0} sequences".format(str(num_deleted))
	connection_data.close()

def delete_sequences_from_experiment(experiment_id_list):
	global db
	global username
	global temp_user_info
	
	
	[db,connection_data] = connect_to_ig_database()
	
	try:		
		username = appsoma_api.environment_get_username()
		temp_user_info = defaultdict(str,db.users.find_one({'user':username}))
	except:
		pass


	if not type(experiment_id_list) is list:
		experiment_id_list = [convert_to_objectid(experiment_id_list)]
	else:
		experiment_id_list = [convert_to_objectid(exp) for exp in experiment_id_list]	
	if temp_user_info['administrator']!=True:
		exps_to_delete =[result['_id'] for result in db.exps.find({'OWNERS_OF_EXPERIMENT':temp_user_info['user'],'_id':{'$in':experiment_id_list}},{'_id':1})]
	else:
		exps_to_delete = experiment_id_list			
	if exps_to_delete == []:
		raise Exception("User {0} does not have access to the listed experiments: {1}".format(temp_user_info['user'],str(experiment_id_list)))
	operation_result = db.seqs.remove({expIdentifier:{'$in':exps_to_delete}})
	db.exps.update({'_id':{'$in':exps_to_delete}},{'$set':{'SEQ_COUNT':0},'$unset':{'ANALYSES_SETTINGS':"",'ANALYSIS_SCHEMA':"",'ANALYSES_COUNT':""}},multi=True)

	num_deleted = operation_result['n']
	print "Deleted {0} sequences".format(str(num_deleted))
	connection_data.close()
	
def delete_analyses_from_experiment(experiment_id_list,analysis_types=[],recombination_types=[]):
	global db
	global username
	global temp_user_info
	
	
	[db,connection_data] = connect_to_ig_database()
	
	try:	
		username = appsoma_api.environment_get_username()
		temp_user_info = defaultdict(str,db.users.find_one({'user':username}))
	except:
		pass


	if not type(experiment_id_list) is list:
		experiment_id_list = [convert_to_objectid(experiment_id_list)]
	else:
		experiment_id_list = [convert_to_objectid(exp) for exp in experiment_id_list]			
	if temp_user_info['administrator']!=True:
		exps_to_delete =[result['_id'] for result in db.exps.find({'OWNERS_OF_EXPERIMENT':temp_user_info['user'],'_id':{'$in':experiment_id_list}},{'_id':1})]
	else:
		exps_to_delete = experiment_id_list	
	if exps_to_delete == []:
		raise Exception("User {0} does not have access to the listed experiments: {1}".format(temp_user_info['user'],str(experiment_id_list)))		
	
	if not analysis_types:
		analysis_query = {'$ne':seqRawData}		
	else:
		if not type(analysis_types) is list:
			analysis_types = [analysis_types]
		#user can never delete 'raw data' sequences. use previous function to do that
		analysis_types = [a for a in analysis_types if analysis_types != seqRawData] 
		analysis_query = {'$in':analysis_types}		
	
	remove_these_analyses_counts = { 'ANALYSES_COUNT.'+analyses_deleted:"" for analyses_deleted in analysis_types}
		
	if recombination_types:
		if not type(recombination_types) is list:
			recombination_types = [recombination_types]						
		operation_result =  db.seqs.remove({expIdentifier:{'$in':exps_to_delete},'ANALYSIS_NAME':analysis_query,'RECOMBINATION_TYPE':{'$in':recombination_types}})						
	else:
		operation_result = db.seqs.remove({expIdentifier:{'$in':exps_to_delete},'ANALYSIS_NAME':analysis_query})
	
	#update counts for each of the analysis types
	exp_metadata = db.exps.find({'_id':{'$in':exps_to_delete}},{'ANALYSES_COUNT':1,'_id':1})	 #get the analyses counts BEFORE the update
	for exp_results in exp_metadata:		
		updated_analyses_counts = {}
		exp_id = exp_results['_id'] 
		#for each experiment in the query, go through each possible analysis
		if 'ANALYSES_COUNT' in exp_results:
			for each_analysis_type in exp_results['ANALYSES_COUNT']:
				#RE-COUNT the number of sequences with this analysis type after 
				num_counts = db.seqs.find({expIdentifier:exp_id,'ANALYSIS_NAME':each_analysis_type}).count()
				if num_counts>0:
					updated_analyses_counts[each_analysis_type] = num_counts
		if updated_analyses_counts:
			db.exps.update({'_id':exp_id},{'$set':{'ANALYSES_COUNT':updated_analyses_counts}})
		else:
			db.exps.update({'_id':exp_id},{'$unset':{'ANALYSES_COUNT':""}})
			
	num_deleted = operation_result['n']
	print "Deleted {0} sequences".format(str(num_deleted))
	
	connection_data.close()


#assume niput file is a TAB delimited file 
#the first column is a JSON dumped dictionary 
#all other columns are used for sorting/grouping together information for the same sequence
#run an awk command to sort sequences 
#merge together any lines, whose sorting fields are the same, the first column from each identical 'sorting field' line is merged into one line seperated by tabs
#once complete, split the file into smaller files based on the number of threads that will be required 
def merge_and_split_db_files(filename):
	working_directory = '/'.join(filename.split('/')[:-1])+'/'
	renamed_input_file = working_directory+os.path.basename(filename).replace(' ','')+'.temp' 
	os.rename(filename,renamed_input_file)
	
	
	awk_init = '''awk -v output_file="{0}"'''.format(filename)
	
	awk = awk_init+''' 'BEGIN{FS="\t";OFS=","}	
		{sort_field="";for(i=2;i<=NF;i++)sort_field=sort_field""OFS""$i}
		FNR==1{out_line=$1;unique_sort=sort_field;next}
		sort_field!=unique_sort{print "["out_line"]" > output_file;out_line=$1;unique_sort=sort_field;next}
		{out_line=out_line""OFS""$1}
		END{print "["out_line"]">output_file}' '''
	
	#sort all sequences sent to file by columsn 2-> end 
	#then run the awk script which will merge all of the same/identical sequences into a single row 
	subprocess.call("sort -t '\t' '{0}' -k2 | {1}".format(renamed_input_file,awk), shell=True)
			
	#now split all of the files into multiple little files so we can use multithreading 
	created_files = useful.split_files_by_seq(filename,multiprocessing_threads,1,False)
	os.remove(filename)
	os.remove(renamed_input_file)
	return created_files




#function for INSERTING new sequence data to the database 
#inserting is different from updating as updating assumes sequence information is already present in database, and annotation is being 'updated' for that sequence
#User can only update a single experiment at a time with the prvoided file => all sequences within the file will be inserted to the experiment defined by  "exp_id"
# append_seqs = True => do not delete sequence data that is currently in the database
#Pseudo code for function 
#1) user passes in a single TAB delimited file that contains information for inserting new sequences to the database 
	#=> The preferred method for generating this TAB file is to use the translator functions found in cchrysostomou_translator_functions. i.e. using the "ConvertFilesForDB" functions. 
	#=> The tab file may be a multi-column file. 
	#=> The tab file  MAY NOTT have a sequence header
	#=> At a minimum the file must have a single column.Column 1 of the tab file must ALWAYS be a json dumped dictionary containing the following keys: 
		#=> SEQUENCE, SEQUENCE_HEADER
		#=> other keys such as 'QUALITY_SCORE', 'FILENAME' are optional (for example we can support fasta files which lack quality score or fastq files which provide a quality score)
		#=> The Key, 'ANALYSIS' is also optional. 
			#IF ANALYSIS IS PROVIDED, THEN ANALYSIS MUST BE A LIST OF SUB DICTIONARIES. Each element of list must have keys and values defining the following: 
				#ANALYSIS_NAME: string
				#RECOMBINATION_TYPE:string 
				#DATA: dictionary of all annotation fields to update the database with 
#2) Elements in the provided TAB delimited file will be sorted by any fields defined after column1 in the input file. Therefore, columns 2=> the end of line are treated as a sorting field 
#3) lines whose sorting fields are identical are merged into a single line where the first column of lines containing the same sorting field are seperated by a tab character 
	#using teh function "ConvertFilesForDB" will always use the sequence_header field as the sorting field 
#4) Files are then split into multiple files to support multithreading 
#5) each file is passed into the function insert_seqs_from_file_real
	#for each line of the generated file, it will read each provided JSON dict by splitting the line by tabs and using JSON.loads for each column 
	#Any JSON dicts in the same line whose SEQUENCE_HEADER AND SEQUENCE DO NOT MATCH will be treated as SEPERATE SEQUENCES
	#BUT Any JSON dicts in the same line whose SEQUENCE_HEADER AND SEQUENCE/REVERSE_COMPLEMENT SEQUENCE MATCH will be treated as the SAME sequence their data will be merged as a new elementin the ANALYSIS_KEY
#6) Each unique sequence will be parsed for special database fields and will be converted into the proper format required for the database and queries. A list of resulting sequences to add to the database are returned
	#at the end. All elements are then inserted into the database using a multiinsert function. 
def benni_insert_seqs_from_file(filename,exp_id='',append_seqs = False,bulk_seqs =50000):	
	global db
	global username
	global temp_user_info
			
	created_files = merge_and_split_db_files(filename)		
	
	if append_seqs == False:
		delete_sequences_from_experiment(exp_id)
	
	#count the number of lines in filename
#	num_lines = useful.file_line_count(filename)
#	#determine the number of lines to split file by
#	split = int(num_lines/multiprocessing_threads)+1
		
#	parent_path = '/'.join(filename.split('/')[:-1])+'/'
#	new_folder = parent_path+'folder_'+str(time.time()).replace('.','')
#	#create a subfolder in the currentfodler, then split the input file by the number of lines determined in split var
#	os.system("mkdir '{0}';split -l {2} '{1}' '{0}/tempfiles' ".format(new_folder,filename,str(split)))
	
#	#get a list of the created files
#	created_files = [new_folder+'/'+f for f in os.listdir(new_folder)]	
	
	inputs=created_files
	
	d = time.time()			
	#setup multithreading
	proc = []
	for fn in inputs:
		#for each input file run the update function
		p = Process(target=benni_insert_seqs_from_file_real,args=(fn,exp_id,append_seqs,bulk_seqs,))
		p.start()
		proc.append(p)
  	for p in proc:
		p.join()
	print time.time()-d
	
	#delete temp database files 
	for fn in inputs:
		os.remove(fn)
	#os.system("rm -r '{0}'".format(new_folder))
	
	

#this should be treated as a private function. Users should be using the public function benni_insert_seqs_from_file for inserting sequence data 
#the input file, filename, should always be created by the merge_and_split_db_files function deefined above 
#assumes the input file is a text file. each line is a JSON DUMP DICT of LISTS.  
#each line SHOULD represent a unique sequence to add to the database, but this function will ensure that the sequence and sequence_header for each line is indeed the same 
#each element in the list of each line is a  json dump dictionary containing different annotation/sequence information
#any line that contains information for different sequences/sequence headers will be split into multiple lines 
#this will append all column data under analysis into the same variable and then proceed to call the insert_seqs function
def benni_insert_seqs_from_file_real(filename,exp_id='',append_seqs = False,bulk_seqs =50000):	
	
	[db,connection_data] = connect_to_ig_database()
		
	try:
		import appsoma_api
		username = appsoma_api.environment_get_username()
		temp_user_info = defaultdict(str,db.users.find_one({'user':username}))
	except:
		pass
	
	
	tr = time.time()	
	count = 0 
	with open(filename) as f:
		seqnum = 1
		seq_list = []
		for seqdata in f:			
			count+=1			
			seqdata = seqdata.strip()			
			unique_sequence = OrderedDict()
			if not seqdata:								
				continue
			#each row of a file can have multiple sequence documents to insert for the same sequence. These dicts are seperated by tabs. 
			all_annotations =json.loads(seqdata.strip()) #[json.loads(s.strip()) for s in seqdata.strip().split('\t') if s.strip()]
											
			for num,each_annotation_field in enumerate(all_annotations):
				#ensure that each sequence is indeed unique: use the sequence and sequence_header field to ensure that we are updating the information for the same exact sequence and not a new sequence
				if 'SEQUENCE' not in each_annotation_field or 'SEQUENCE_HEADER' not in each_annotation_field:		
					print each_annotation_field
					print 'Error code insertion: 0'
					continue
				unique_key = each_annotation_field['SEQUENCE_HEADER']+'.'+each_annotation_field['SEQUENCE'].upper()
				if not unique_sequence:
					unique_sequence[unique_key] = each_annotation_field
					if 'ANALYSIS' not in unique_sequence[unique_key]:
						unique_sequence[unique_key]['ANALYSIS'] = []
				elif unique_key in unique_sequence:
					if 'ANALYSIS' in each_annotation_field:
						unique_sequence[unique_key]['ANALYSIS'].extend(each_annotation_field['ANALYSIS'])
				else:				
					#this sequence/sequence_header combination has not been observed yet 
					#try looking for reverse complement					
					unique_key_2 = each_annotation_field['SEQUENCE_HEADER']+'.'+useful.Reverse_Complement(each_annotation_field['SEQUENCE'].upper())
					if unique_key_2 in unique_sequence:						
						#the reverse complement of the sequence was found 
						if 'ANALYSIS' in each_annotation_field:
							unique_sequence[unique_key_2]['ANALYSIS'].extend(each_annotation_field['ANALYSIS'])															
					else:
						#OK so we need to add a new set of sequences because this specific sequence/header combination has not been seen 
						unique_sequence[unique_key] = each_annotation_field
						if 'ANALYSIS' not in each_annotation_field:
							unique_sequence[unique_key]['ANALYSIS'] = []														
			seq_list.extend(unique_sequence.values())	
																			
			#main_annotation = all_annotations[0]
			#analyses = main_annotation.pop('ANALYSIS',None)
			#if not analyses:
			#	analyses = []
			#for other_annotations_ind in range(1,len(all_annotations)):
			#	other_annotations = all_annotations[other_annotations_ind]
			#	if 'ANALYSIS' in other_annotations:
			#		analyses.extend(other_annotations['ANALYSIS'])
			#if analyses:
			#	main_annotation['ANALYSIS'] = copy.deepcopy(analyses)														
			#seq_list.append(main_annotation)				
			seqnum+=1
			if seqnum%bulk_seqs == 0:																
				benni_insert_sequences(seq_list,exp_id)	
				seq_list = []
			if seqnum%100000==0:
				print "Inserted {0} sequences".format(str(seqnum))
			
		if seq_list:
			benni_insert_sequences(seq_list,exp_id)	
		print "Inserted {0} sequences".format(str(seqnum))
	print 'myfinaltime: '+str(time.time()-tr)
	connection_data.close()		

def benni_insert_sequences(seq_data, experiment_obj_id, user_info=None):
	"""
	Input: a JSON string of sequence info, or list of JSON documents for sequences, and the Object ID of the experiment the sequence belongs to. 
	(This experiment Object ID should be checked before the function is called.)
	Output: ?
	Insert single record into sequence database
	"""		
	experiment_obj_id = convert_to_objectid(experiment_obj_id)
	
	if temp_user_info['administrator']:
		experiment = db.exps.find_one({'_id': experiment_obj_id})		
		if not(experiment):
			raise Exception("Can't find experiment with Object ID: {}".format(str(experiment_obj_id)))
	else:
		experiment = db.exps.find_one({'_id': experiment_obj_id,'OWNERS_OF_EXPERIMENT':temp_user_info['user']})	
		if not(experiment):
			raise Exception("Either the experiment does not exist or the User {0} does not have write access to the following experiment: {1}".format(temp_user_info['user'],str(experiment_obj_id)))	
	
	#if we have already stored a schema for this experiment, then store the previous schema in this varaible
	#we use flatten_dictionary function to conver the dictionary into a flat dicionatry of DOT notation (subdocuments are seperated by '.')
	#convert current_schema from lists into sets 
	current_schema = {field:set(value) for field,value in useful.flatten_dictionary(experiment['ANALYSIS_SCHEMA']).iteritems()} if 'ANALYSIS_SCHEMA' in experiment else {}				
	current_schema = defaultdict(set,current_schema)
	
	#if we have already stored a list of settings for the experiment from previous sequences then store these analyses into our variable, if not then use empty list 	
	current_analyses_commands_list = experiment['ANALYSES_SETTINGS'] if 'ANALYSES_SETTINGS' in experiment else []
	current_analyses_commands_dict = defaultdict(lambda:-1,{command:index for index,command in enumerate(current_analyses_commands_list)})
	
	#if we have already stored a list of filenames for sequences in expermient from previous sequences then store these filenames into our variable, if not then use empty list 	
	current_filenames_list = experiment['FILENAMES'] if 'FILENAMES' in experiment else []
	current_filenames_dict = defaultdict(lambda:-1,{fn:index for index,fn in enumerate(current_filenames_list)})
	
	command_list_len = len(current_analyses_commands_list)	
	initial_command_list_len = command_list_len
	
	initial_filenames_list_len = len(current_filenames_list)
	
	dateupdated = time.strftime('%D')
					
	# add exp obj_id to seq data, insert sequence, update exp count
	if not isinstance(seq_data, list):
		seq_data = [seq_data]

	seq_analysis_insert_list = []
	seq_analysis_counter = Counter()

	#print seq_data
	#t = time.time()
	for seq in seq_data:	# casting to list in case seq_data is a single dictionary
		##the following fields are added to raw seq doc by server
		seq_doc = {}
		seq_doc[expIdentifier] = experiment_obj_id #currently expIdentifier = 'EXP_ID'
		seq_key = ObjectId()  #for sequences we want _id and SEQ_KEY to be identical, so we generate them ahead of insertion
		seq_doc['_id'] = seq_key
		seq_doc[idIdentifier] = seq_key #currently idIdentifier = 'SEQ_ID'
		seq_doc['ANALYSIS_NAME'] = seqRawData #currently this string is '@SEQ'
		seq_doc['DATE_UPDATED'] = dateupdated
		##		
		
		if not set(seq.keys()) >= set(schema.Seqs_Collection()['required_fields'].keys()) - set([expIdentifier, idIdentifier, 'ANALYSIS_NAME']):
			missing_fields = set(schema.Seqs_Collection()['required_fields'].keys()) - set(seq_data.keys())
			#### DO WE WANT TO RAISE AN EXCEPTION OR WRITE TO A LOG FILE??
			# raise Exception('The experiment data is missing these fields required by the database schema: '.format(', '.join(missing_fields)))
			print 'The sequence data is missing these fields required by the database schema: '.format(', '.join(missing_fields))
			continue#do nothing for now
		
		analysis_list = seq.pop('ANALYSIS', None)
								
		
		mod_seq_data = {}
		for key,value in seq.iteritems():
			key = key.upper().replace(' ','_')				
			value = schema_fields_and_data_types[key]['to_db'](value)# if key in schema_fields_and_data_types else schema.default_fields_data_types(value)
			mod_seq_data[key] = value		
			
		
		seq_doc['DATA'] =mod_seq_data #(append all user defined data about the sequence to DATA subdocument)
		seq_analysis_insert_list.append(useful.removeEmptyVals(seq_doc))		
		if analysis_list:
			seq_analysis_insert_list.extend(generate_analyses_list(analysis_list, seq_key, experiment_obj_id,dateupdated))	
	
	seq_analysis_counter = Counter()
	#print 'step1: '+str(time.time()-t)
	#go through each possible analysis 
	#t = time.time()
	for i,analysis in enumerate(seq_analysis_insert_list):
		
		if analysis['ANALYSIS_NAME']!=seqRawData:
			#if there is a setting for how this sequence was edited, then remove it from DATA field
			#Next, check if setting currently exists in our list of analyses_settings
			#if not append to list 
			analysis_command = analysis['DATA'].pop('COMMAND',None)
			
			#this is not rawd ata, it is annotated data
			seq_analysis_counter['ANALYSES_COUNT.'+analysis['ANALYSIS_NAME']+'.'+analysis['RECOMBINATION_TYPE']]+=1
			
			flat_data = analysis['DATA'] #useful.flatten_dictionary(analysis['DATA'])
			
			if analysis_command:
				analysis_command_num = []
				for each_analysis_command in analysis_command:
					index_position = current_analyses_commands_dict[each_analysis_command]
					if index_position==-1: #command/string not found yet
						#append to list of commands
						current_analyses_commands_list.append(each_analysis_command)
						index_position = command_list_len
						current_analyses_commands_dict[each_analysis_command] = index_position #sotre new index position
						command_list_len+=1#add to total length of list
					analysis_command_num.append(index_position)
				analysis['SETTINGS'] = list(set(analysis_command_num))
			else:
				analysis_command = None
															
			for field,value in flat_data.iteritems():
				#go through each schema/annotated field
				#add analysis type to the current schema for that field			
								
				current_schema[field+'.ANALYSES'].add(analysis['ANALYSIS_NAME'])
				#add the datatype for the value of this field for that current schema
				current_schema[field+'.DATATYPE'].add(return_data_types(value))
		else:
		
			#for raw data, just increase sequence count
			seq_analysis_counter['SEQ_COUNT']+=1
			if 'FILENAME' in analysis['DATA']:
				if analysis['DATA']['FILENAME'] not in current_filenames_dict:
					current_filenames_dict[analysis['DATA']['FILENAME']] = len(current_filenames_list)
					current_filenames_list.append(analysis['DATA']['FILENAME'])
				analysis['DATA']['FILENAME'] = current_filenames_dict[analysis['DATA']['FILENAME']]
					
				
		
		#fields under DATA and QUERY_DATA are currently in dot notation 
		#before inserting to database, make sure to unflatten the fields under DATA and QUERY_DATA using DotAccessible
		if 'QUERY_DATA' in analysis:
			seq_analysis_insert_list[i]['QUERY_DATA'] = dict(DotAccessible(analysis['QUERY_DATA']))
		if 'DATA' in analysis:
			seq_analysis_insert_list[i]['DATA'] = dict(DotAccessible(analysis['DATA']))
	
	#with open('scratch/testingdbstuff.txt','w') as fout:
	#	for line in seq_analysis_insert_list:
	#		fout.write(bson_dumps(line,indent=4)+'\n')
	
	#print 'step 2: '+str(time.time()-t)
		
	#r1 = time.time()			
	db.seqs.insert_many(seq_analysis_insert_list,ordered=False)
	#print 'totaltime: '+str(time.time()-r1)	
	update_command = {'$inc':seq_analysis_counter}
	
					
	set_fields = {'ANALYSIS_SCHEMA.'+field:list(value) for field,value in current_schema.iteritems()} #add any fields detected in current schema 
	if command_list_len != initial_command_list_len: #add any fields detected in analyses 
		set_fields['ANALYSES_SETTINGS'] = current_analyses_commands_list	
		
	if initial_filenames_list_len!=len(current_filenames_list):#add new filenames to the experiment metadata
		set_fields['FILENAMES'] = current_filenames_list
	
	if set_fields: #there are fields that need to be modified in document 
		update_command['$set'] = set_fields
						
	db.exps.update({'_id':experiment_obj_id},update_command)
	
	#db.experiments.update({'_id': ObjectId(experiment_obj_id)}, {'$inc': exp_seq_analyses_counter})
	#print seq_data
	#return inserted_seq_ids
	
def generate_analyses_list(analysis_data, sequence_obj_id, experiment_obj_id, dateupdated=None, user_info=None):
	# Assuming that checking that the caller has permissions has already been done.
	if not dateupdated:
		dateupdated = time.strftime('%D')
	if not isinstance(analysis_data, list):
		analysis_data = [analysis_data]				
		
	# check there isn't a repeated analysis/recomb combo in the list of analysis types
	# (this is only necessary for inserts, for updates repeated types will overwrite each other)
	analysis_recomb_combos = {} #(analysis['ANALYSIS_NAME'], analysis['RECOMBINATION_TYPE']) for analysis in analysis_data]
	#if len(analysis_recomb_combos) != len(set(analysis_recomb_combos)):
	#	raise Exception('One ANALYSIS_NAME/RECOMB_TYPE combination appears more than once in the analysis data.')
	
	analysis_data_modified = []
	
	# check that analysis required fields are in analysis_data
	for i,analysis in enumerate(analysis_data):
		
		modified_data = {}			
		#print analysis['DATA']
		for key,value in analysis['DATA'].iteritems():
			key = key.upper().replace(' ','_')				
			#defaultdict in the variable (see immunogrep_database_schema.schema_fields_and_data_types) ensures that any key not accoutned for in varaible is given a default 
			#value for to_db and to_file settings
			value = schema_fields_and_data_types[key]['to_db'](value)# if key in schema_fields_and_data_types else schema.default_fields_data_types(value)
			modified_data[key] = value
		#analysis_data[i]['DATA'] = modified_data.copy()
	
		
		#this function in database schema module defines a few finalized treatments of the data
		#for example, it adds CDR lengths to the data and creates a sepreate fields called query_data_fields for improved queries of certain fields such as V-GENES 
		[modified_data,query_data_fields] = schema.AddedSeqCollectionFields(modified_data)
		#these fields contain special treatment of cretain fields for improved queries. 
		#for example we transform VGENES, such that both genes and alleles are stored as a list 
		#SEE AddedSeqCollectionFelds function in database schema for details
		if query_data_fields:
			analysis_data[i]['QUERY_DATA'] = query_data_fields# dict(DotAccessible(query_data_fields))
				
		if not set(analysis.keys()) >= set(['ANALYSIS_NAME', 'RECOMBINATION_TYPE', 'DATA']): # these required fields should probably be defined in schema
			missing_fields = set(['ANALYSIS_NAME', 'RECOMBINATION_TYPE', 'DATA']) - set(analysis.keys())
			raise Exception('The analysis data is missing these fields required by the database schema: '.format(', '.join(missing_fields)))
		analysis_data[i][idIdentifier] = sequence_obj_id  # Note ObjectId is an idempotent function
		analysis_data[i][expIdentifier] = experiment_obj_id
		analysis_data[i]['DATE_UPDATED'] = dateupdated
		#analysis['DATA'] = dict(DotAccessible(analysis['DATA']))
		analysis_data[i]['DATA'] = modified_data# dict(DotAccessible(modified_data))
				
		# check there isn't a repeated analysis/recomb combo in the list of analysis types
		# (this is only necessary for inserts, for updates repeated types will overwrite each other)
		if analysis['ANALYSIS_NAME']+'.'+analysis['RECOMBINATION_TYPE'] not in analysis_recomb_combos:
			a_info = analysis_data[i]
			a_info['DATA'] = useful.removeEmptyVals(a_info['DATA'])
			if 'COMMAND' in a_info['DATA']:
				a_info['DATA']['COMMAND'] = [a_info['DATA']['COMMAND']]							
			analysis_recomb_combos[analysis['ANALYSIS_NAME']+'.'+analysis['RECOMBINATION_TYPE']]=len(analysis_recomb_combos)
			analysis_data_modified.append(a_info) #remove NONE VALUES, EMPTY LISTS, EMPTY STRINGS 
		#	pprint(analysis_data_modified[-1])
		else:
			#now we have to modify the analysis name/recombinatio type with new information
			more_info = useful.removeEmptyVals(analysis_data[i]['DATA'])
			ind_pos = analysis_recomb_combos[analysis['ANALYSIS_NAME']+'.'+analysis['RECOMBINATION_TYPE']]
			#pprint(analysis_data_modified[ind_pos])
			if 'COMMAND' in more_info:
				if 'COMMAND' in analysis_data_modified[ind_pos]['DATA']:
					analysis_data_modified[ind_pos]['DATA']['COMMAND'].append(more_info.pop('COMMAND'))
				else:
					analysis_data_modified[ind_pos]['DATA']['COMMAND'] = [more_info.pop('COMMAND')]
			analysis_data_modified[ind_pos]['DATA'].update(more_info) 
			
		
#	print analysis_data
#	print analysis_counter
	return analysis_data_modified

def benni_update_seqs_from_file(filename,bulk_seqs = 50000,update_replace=False):	
		
	'''
	#count the number of lines in filename
	num_lines = useful.file_line_count(filename)
	#determine the number of lines to split file by
	split = int(num_lines/multiprocessing_threads)+1
		
	parent_path = '/'.join(filename.split('/')[:-1])+'/'
	new_folder = parent_path+'folder_'+str(time.time()).replace('.','')
	#create a subfolder in the currentfodler, then split the input file by the number of lines determined in split var
	os.system("mkdir '{0}';split -l {2} '{1}' '{0}/tempfiles' ".format(new_folder,filename,str(split)))
	
	#get a list of the created files
	created_files = [new_folder+'/'+f for f in os.listdir(new_folder)]		
	'''
	created_files = merge_and_split_db_files(filename)
	
	inputs=created_files
	
	d = time.time()			
	#setup multithreading
	proc = []
	for fn in inputs:
		#for each input file run the update function
		p = Process(target=benni_update_seqs_from_file_real,args=(fn,bulk_seqs,update_replace,))
		p.start()
		proc.append(p)
  	for p in proc:
		p.join()
	print time.time()-d
	
	#delete new subfolder
	#os.system("rm -r '{0}'".format(new_folder))
	
	
	

def benni_update_seqs_from_file_real(filename,bulk_seqs = 50000,update_replace=False):	
	global db
	global username
	global temp_user_info
	print 'starting'
	
	[db,connection_data] = connect_to_ig_database()
	
	try:		
		username = appsoma_api.environment_get_username()
		temp_user_info = defaultdict(str,db.users.find_one({'user':username}))
	except:
		pass


	t1 = time.time()
	with open(filename) as f:
		seqnum = 0
		seq_list = []
		for seqdata in f:
			if not seqdata.strip():
				continue
			#each row of a file can have multiple sequence documents to update, each document is an element in the json list 
			all_annotations =json.loads(seqdata.strip()) #[json.loads(s.strip()) for s in seqdata.strip().split('\t') if s.strip()]
			main_annotation = all_annotations[0]
			analyses = main_annotation.pop('ANALYSIS',None)
			if not analyses:
				analyses = []
			for other_annotations_ind in range(1,len(all_annotations)):
				other_annotations = all_annotations[other_annotations_ind]
				if 'ANALYSIS' in other_annotations:
					analyses.extend(other_annotations['ANALYSIS'])
			if analyses:
				main_annotation['ANALYSIS'] = copy.deepcopy(analyses)														
			seq_list.append(main_annotation)								
		
			seqnum+=1
			if seqnum%bulk_seqs == 0:				
				benni_update_analyses(seq_list,update_replace)
				seq_list = []
				print 'UPDATED {0} SEQUENCES'.format(str(seqnum))
		if seq_list:
			benni_update_analyses(seq_list,update_replace)
	print 'total time: '+str(time.time()-t1)

	connection_data.close()

def benni_update_analyses(analysis_data, update_replace, user_info=None):
	dateupdated = time.strftime('%D')
	#get list of allowed experiments 
	allowed_exps =[result['_id'] for result in db.exps.find({'OWNERS_OF_EXPERIMENT':temp_user_info['user']},{'_id':1})]	
	
	if not isinstance(analysis_data, list):
		analysis_data = [analysis_data]

	experiment_upsert_counts = defaultdict(Counter)
	
	#store schema for all experiments currently being updated 
	#for each experiment/key in variable store:
		#sets of the following: each annotated field from database for that experiment
			#for each field/key, store: analysis name which contains it and the datatype for that field 
		#example: current_schema_by_exp = {'EXP':{'VREGION.VGENES.ANALYSES':['IMGT','IGBLAST'],
			#'VREGION.VGENES.DATATYPES':[list_of_string]}}	
			
	current_exp_data = {}
	current_schema_by_exp = {}
	time_counter = 0
	time_counter2 = 0
	schema_timer=0
	dict_timer=0
	modify_timer=0
	command_timer=0
	datatype_timer=0
			
	issues = 0
	
	updates = []
	#out_updates = []
	tc = time.time()
	for analysis in analysis_data:
		# Get seq and exp ids, and check that they already exist in the database.
		if idIdentifier not in analysis.keys():
			#no seq identifier found 
			print "USER DID NOT PASS SEQID"
			issues+=1
			continue
		
		analysis_list = analysis.pop('ANALYSIS', None)
		if not analysis_list:
			#user did not pass in 'ANALYSIS' as key 
			print "USER DID NOT PASS ANALYSIS KEY"
			issues+=1
			continue
				
		sequence_obj_id = convert_to_objectid(analysis[idIdentifier])
				
		#query to ensure that the sequence object id exists 
		#'_id' of RAW nucleotide sequences will ALWAYS be equal to the sequence_obj_id. 
		#in addition, RAW nucleotide sequences will ALWAYS have an analysis_name corresponding to seqRawData variable in global variables script
		
		
		if temp_user_info['administrator']:
			#if administrator then dont need to search by exp_id
			seq_query  = db.seqs.find({'_id': sequence_obj_id}, {expIdentifier: 1}).limit(1).hint('_id_').next()
		else:
			#for non-administrators, only allow updates on experiments where write access was granted
			seq_query  = db.seqs.find({'_id': sequence_obj_id,expIdentifier:{'$in':allowed_exps}}, {expIdentifier: 1}).limit(1).hint('_id_').next()					
				
		if seq_query:
			#sequence_ID found and user allowed to edit sequence 
			experiment_obj_id = seq_query[expIdentifier]
		else:			
			#sequence_ID not found or do not have WRITE PERMISSION to this sequence 
			print(sequence_obj_id)
			issues+=1
			print("SEQ_ID_NOT_FOUND!!")
			continue		
			#raise Exception("User {0} does not have write access to the following experiment: {1}".format(temp_user_info['user'],str(experiment_obj_id)))
		
				
		#this is the first time this experiment has been updated
		if str(experiment_obj_id) not in current_schema_by_exp:
			new_experiment = db.exps.find_one({'_id': experiment_obj_id})	
			#get current schema for experiment
			current_analyses_command_list = new_experiment['ANALYSES_SETTINGS'] if 'ANALYSES_SETTINGS' in new_experiment else []		#get current analysis commands 
			current_schema_by_exp[str(experiment_obj_id)] = {
				'schema': defaultdict(set,{field:set(value) for field,value in useful.flatten_dictionary(new_experiment['ANALYSIS_SCHEMA']).iteritems()} if 'ANALYSIS_SCHEMA' in new_experiment else {}),
				'analyses_commands': current_analyses_command_list,
				'analyses_commands_dict':defaultdict(lambda:-1,{command:index for index,command in enumerate(current_analyses_command_list)}),#get current anlayses commands as a dict
				'analyses_commands_dict_len':len(current_analyses_command_list),#get length of analyses list
				'analyses_commands_initial_dict_len':len(current_analyses_command_list)#store initial len for later
			}
	
		
		#first go through each of the analyses and group them by ANALYSIS_NAME and RECOMBINATION_TYPE
		grouped_analysis_list = []
		analysis_command_list = []
		combos = {}
		for each_analysis in analysis_list:			
			
			# Check necessary fields are present
			if not set(each_analysis.keys()) >= set(['ANALYSIS_NAME', 'RECOMBINATION_TYPE', 'DATA']): # these required fields should probably be defined in schema
				missing_fields = set(['ANALYSIS_NAME', 'RECOMBINATION_TYPE', 'DATA']) - set(each_analysis.keys())
				issues+=1
				print("Error code 0")
				continue
				#raise Exception('The analysis data is missing these fields required by the database schema: '.format(', '.join(missing_fields)))
			
			if each_analysis['ANALYSIS_NAME'] == seqRawData:
				#this is not allowed. no user can update an analysis type using the analysis type that corresponds to the raw sequence 
				issues+=1
				print("Error code 1")
				continue			
			
			modified_data = {}	
			#modify the keys and values for each analysis 
			for key,value in each_analysis['DATA'].iteritems():				
				key = key.upper().replace(' ','_')
				if value!=None:				
					#defaultdict in the variable (see immunogrep_database_schema.schema_fields_and_data_types) ensures that any key not accoutned for in varaible is given a default 
					#value for to_db and to_file settings
					value = schema_fields_and_data_types[key]['to_db'](value)# if key in schema_fields_and_data_types else schema.default_fields_data_types(value)
					modified_data[key] = value
			
			analysis_command = modified_data.pop('COMMAND',None)		
						
			#now use the analysis_name and recombination_type to group together results for the same doucment together 
			if each_analysis['ANALYSIS_NAME']+'.'+each_analysis['RECOMBINATION_TYPE'] in combos:
				index_pos = combos[each_analysis['ANALYSIS_NAME']+'.'+each_analysis['RECOMBINATION_TYPE']]								
				grouped_analysis_list[index_pos]['DATA'].update(modified_data)								
				if analysis_command:					
					grouped_analysis_list[index_pos]['SETTINGS'].append(analysis_command)																						
			else:
				each_analysis['DATA'] = modified_data
				combos[each_analysis['ANALYSIS_NAME']+'.'+each_analysis['RECOMBINATION_TYPE']] = len(grouped_analysis_list)
				grouped_analysis_list.append(each_analysis)					
				if analysis_command:
					grouped_analysis_list[-1]['SETTINGS'] = [analysis_command]																
				else:
					grouped_analysis_list[-1]['SETTINGS'] = []
			
		for ind,each_analysis in enumerate(grouped_analysis_list):
			query_dict = {idIdentifier: sequence_obj_id,expIdentifier: experiment_obj_id, 'ANALYSIS_NAME': each_analysis['ANALYSIS_NAME'], 'RECOMBINATION_TYPE': each_analysis['RECOMBINATION_TYPE']}												
			#this function in database schema module defines a few finalized treatments of the data
			#for example, it adds CDR lengths to the data and creates a sepreate fields called query_data_fields for improved queries of certain fields such as V-GENES 
			modified_data = each_analysis['DATA']
			[modified_data,query_data_fields] = schema.AddedSeqCollectionFields(modified_data)			
			
			#update EXPERIMENT SCHEMA with nonempty fields. 
			for field,value in modified_data.iteritems():# useful.removeEmptyVals(modified_data).iteritems():
#				#go through each schema/annotated field
#				#add analysis type to the current schema for that field			
				if value or value==False:
					current_schema_by_exp[str(experiment_obj_id)]['schema'][field+'.ANALYSES'].add(each_analysis['ANALYSIS_NAME'])
	#				#add the datatype for the value of this field for that current schema
					current_schema_by_exp[str(experiment_obj_id)]['schema'][field+'.DATATYPE'].add(return_data_types(value))
			
			analysis_command = each_analysis.pop('SETTINGS',None)
			if analysis_command:	
				analysis_number_command = []
				for each_command in analysis_command:
					index_position = current_schema_by_exp[str(experiment_obj_id)]['analyses_commands_dict'][each_command]
					if index_position==-1: #command/string not found yet
						#append to list of commands
						current_schema_by_exp[str(experiment_obj_id)]['analyses_commands'].append(each_command)
						index_position = current_schema_by_exp[str(experiment_obj_id)]['analyses_commands_dict_len']
						current_schema_by_exp[str(experiment_obj_id)]['analyses_commands_dict'][each_command] = index_position #sotre new index position
						current_schema_by_exp[str(experiment_obj_id)]['analyses_commands_dict_len']+=1#add to total length of list
					analysis_number_command.append(index_position)
				analysis_number_command = list(set(analysis_number_command))
			else:
				analysis_number_command = []
												
			if update_replace: #do not use set command. instead replace everything under 'DATA' and 'QUERY_DATA'																			
				modified_data = dict(DotAccessible(modified_data)) #remove NONE VALUES, EMPTY LISTS, EMPTY STRINGS 																						
				
				update_query = {'DATE_UPDATED':dateupdated,'DATA':modified_data}				
				if analysis_number_command:
					update_query['SETTINGS'] = analysis_number_command# list(set(analysis_number_command))
				
				#these fields contain special treatment of cretain fields for improved queries. 
				#for example we transform VGENES, such that both genes and alleles are stored as a list 
				#SEE AddedSeqCollectionFelds function in database schema for details
				query_data_fields = dict(DotAccessible(query_data_fields))
												
				update_query = useful.removeEmptyVals(update_query)
				if query_data_fields:
					#response = db.seqs.update(query_dict, {'$set':update_query}, upsert=True)
					update_query['QUERY_DATA'] = query_data_fields
					updates.append(pymongo.UpdateOne(query_dict,{'$set':update_query}, upsert=True))
				else:
					#response = db.seqs.update(query_dict, {'$set':update_query}, upsert=True)
					updates.append(pymongo.UpdateOne(query_dict,{'$set':update_query,'$unset':{'QUERY_DATA':""}}, upsert=True))
				#out_updates.append(query_dict)
				#out_updates.append(update_query)
			else:	#use set command and therefore only change individual fields 									
				modified_data = {'DATA.'+field:value for field,value in modified_data.iteritems()}
				modified_data['DATE_UPDATED'] = dateupdated
												
				add_commands = {'SETTINGS':{'$each':list(set(analysis_number_command))}}
												
				#these fields contain special treatment of cretain fields for improved queries. 
				#for example we transform VGENES, such that both genes and alleles are stored as a list 
				#SEE AddedSeqCollectionFelds function in database schema for details
				for query_field,value in query_data_fields.iteritems():
					modified_data['QUERY_DATA.'+query_field] = value
				
				[set_non_empty_fields,unset_empty_fields] = useful.divideEmptyAndNonEmptyVals(modified_data) #go through dictionary and seperate the fields based on non-empty fields and empty-fields 				
				
				if unset_empty_fields: #we need to remove empty fields					
					updates.append(pymongo.UpdateOne(query_dict, {'$set': set_non_empty_fields,'$unset':unset_empty_fields,'$addToSet':add_commands}, upsert=True))
					#out_updates.append(query_dict)
					#out_updates.append(set_non_empty_fields)
					#out_updates.append(unset_empty_fields)
					#out_updates.append(add_commands)
					#response = db.seqs.update(query_dict, {'$set': set_non_empty_fields,'$unset':unset_empty_fields}, upsert=True)
				else: #no need to remove empty fields
					#response = db.seqs.update(query_dict, {'$set': set_non_empty_fields}, upsert=True)
					updates.append(pymongo.UpdateOne(query_dict, {'$set': set_non_empty_fields,'$addToSet':add_commands}, upsert=True))
					#out_updates.append(query_dict)
					#out_updates.append(set_non_empty_fields)
					#out_updates.append(add_commands)
			
			#keep track of all analysis types and recombination types added in this experiment 
			experiment_upsert_counts[str(experiment_obj_id)][each_analysis['ANALYSIS_NAME']+'.'+each_analysis['RECOMBINATION_TYPE']]+=1
			# Check if the update was an insert, and if so, add 1 to the count dictionary for that analysis type						
			#if (response['n'] == 1) and not response['updatedExisting']:				
			#	experiment_upsert_counts[str(experiment_obj_id)][each_analysis['ANALYSIS_NAME']+'.'+each_analysis['RECOMBINATION_TYPE']] += 1
	print('and inserting')
	ttest=time.time()-tc
	
	#with open('scratch/testingdbstuff.txt','w') as fout:
	#	for line in out_updates:
	#		fout.write(bson_dumps(line,indent=4)+'\n')
	
	
	db.seqs.bulk_write(updates,ordered=False)
#	dbtime=time.time()-ttest
#	#print 'total: '+str(time.time()-tc)
	print 'dbtime: '+str(ttest)
##	print 'querytime: '+str(time_counter)
##	print 'settingup_updates: '+str(time_counter2)
##	print 'dict timer: '+str(dict_timer)
##	print 'modify_timer: '+str(modify_timer)
##	print 'command timer: '+str(command_timer)
##	print 'datatype_timer: '+str(datatype_timer)
##	print 'schema_timer: '+str(schema_timer)
		
	print issues
	#print experiment_upsert_counts
	
	# Update counts and experiment schema in the affected experiments
	for experiment_obj_id_str, update_settings in current_schema_by_exp.iteritems():
		exp_schema = update_settings['schema']
		
		#for each experiment, figure out the number of documents with a specfic recombination type and analysis name 
		analyses_counter = {}
		count_these_analysis_types  = experiment_upsert_counts[experiment_obj_id_str]# in experiment_upsert_counts.iteritems():				
		for analysis,v in count_these_analysis_types.iteritems():
			analysis_recomb = analysis.split('.')
			analyses_counter[analysis] = db.seqs.find({expIdentifier:ObjectId(experiment_obj_id_str),'ANALYSIS_NAME':analysis_recomb[0],'RECOMBINATION_TYPE':analysis_recomb[1]}).count()	
		
		set_fields = {'ANALYSIS_SCHEMA.'+field:list(value) for field,value in exp_schema.iteritems()}	#add the schema information 			
		if update_settings['analyses_commands_initial_dict_len']!=update_settings['analyses_commands_dict_len']: #check to see if we need to add new analyses commands
			set_fields['ANALYSES_SETTINGS'] = update_settings['analyses_commands']
		
		update_command = {}
		if analyses_counter:
			#update_command['$inc'] = analyses_counter # do we need to update analyses counts
			for all_analyses,counts in analyses_counter.iteritems():
				set_fields['ANALYSES_COUNT.'+all_analyses] = counts
		if set_fields:			
			update_command['$set'] = set_fields
		
		
		
		if update_command:#do we need to update analyses
			db.exps.update({'_id':ObjectId(experiment_obj_id_str)},update_command)
	
		
