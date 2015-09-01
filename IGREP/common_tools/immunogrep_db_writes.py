from collections import defaultdict, namedtuple, Counter, deque
import immunogrep_database_schema as schema
from immunogrep_database_schema import convert_text_to_index_field_text
import immunogrep_proxy_tools as proxy_tools
#from immunogrep_database_schema import *
from immunogrep_useful_functions import DotAccessible
#try:
#	import cchrysostomou_appsoma_ab_members as members
#except:
#	import members
import inspect
import sys
import os
from functools import wraps
import json
# modules used for md5 calculation
from functools import partial
import hashlib
# module used for running function in new processing
# from multiprocessing import Process
import subprocess, shlex, signal
from itertools import islice
from bson.objectid import ObjectId

#-- v -- This import should be removed in production, connecting to the database should be handled with some other function.
from pymongo import MongoClient

import time #CC
import datetime #CC 
import immunogrep_useful_functions as useful #CC 

db = None


def allele_name_to_list(allele):
	gene_name = []
	if ' ' in allele:
		gene_name.append(allele)
	for word in allele.split():
		gene_name.append(word)
		if '*' in word:
			gene_name += word.split('*')
	return gene_name



#-------------------------------------
# Inserting experiments and sequences
#-------------------------------------

#NEW insert_experiment function (previous function found in version 52#
#TO DO: sample preperation date -> enforce it as a date variable of month/year??
#TO DO: IF CURATED = FALSE --> SEND EMAIL TO ADMINISTRATOR 
def insert_experiment(exp_data, user_info=None):
	allowed_data_types = [float,int,str,unicode,bool]
	"""
		Input: Database collection information, dictionary describing the data for an experiment.
		Output: _id of the new experiment.
		Insert the new experiment given in the data. Check that the EXPERIMENT_ID doesn't already exist; if it does, throw error. Check that the user inserting the experiment is listed in WRITE_ACCESS.
	"""
	# check if user is in the Immunogrep database and has write access
	user_info = defaultdict(str,user_info)
	if not user_info or not user_info['user'] or user_info['write_access'] == False:
		raise Exception('{0}({1}) is not allowed to create immunogrep experiments.'.format(user_info['user'],user_info['name']))	
		
	#check if dict passed in
	if type(exp_data) is not dict:
		raise Exception('A dictionary must be passed in containing experiment metadata')
		
	#convert all keys to upercase and remove spaces
	#remove any doublequotes found in the values of every field 
	exp_data_modified = {field.strip().upper().replace(' ','_'):exp_data[field].replace('"','') if isinstance(exp_data[field],basestring) else exp_data[field] for field in exp_data}
	
	#REMOVE any keys that client adds in/has control over values
	for field in schema.Exps_Collection()['protected_fields']:
		exp_data_modified.pop(field,None)	
	
	required_fields = schema.Exps_Collection()['required_fields']	
	#remove any fields that are empty lists, empty dicts, strings, or None
	exp_data_modified = useful.removeEmptyVals(exp_data_modified)				
	# check if keys in exp_data are a superset of required field keys	
	if not set(exp_data_modified.keys()) >= set(required_fields):
		missing_fields = set(required_fields) - set(exp_data_modified.keys())
		raise Exception('The experiment data is missing these fields required by the database schema: {0} '.format(','.join(missing_fields)))
			
	#first, go through each of the fields, and if we encounter a field whose data type we control, then modify it	
	#next, ensure that all datatypes provided are allowable
	for field in exp_data_modified:
		#the variable exp_fields_and_data_types defines how we will treat specific fields in the database. codes a function to process field values
		#if the key is not present in the vaiable exp_fields_and_data_types, then the defaultdict will ensure that we treat novel fields as defined by the defaultdict function. see the variable immunogrep_database_schema.exp_fields_and_data_types variable for details
		exp_data_modified[field] = schema.exp_fields_and_data_types[field](exp_data_modified[field]) #if field in schema.exp_fields_and_data_types else exp_data_modified[field] 
		if type(exp_data_modified[field]) is list:
			for each_member in exp_data_modified[field]:
				if type(each_member) not in allowed_data_types:
					raise Exception('The following field, {0}, contains an invalid value: {1} datatype is {2}. Only the following datatypes are allowed: {3}'.format(field,str(each_member),str(type(each_member)), ','.join([str(dt) for dt in allowed_data_types])))
		elif type(exp_data_modified[field]) not in allowed_data_types:
			raise Exception('The following field, {0}, contains an invalid value: {1} datatype is {2}. Only the following datatypes are allowed: {3}'.format(field,str(exp_data_modified[field]),str(type(exp_data_modified[field])), ','.join([str(dt) for dt in allowed_data_types])))
			
	#AGAIN remove any fields that are empty lists, strings, or None
	#just incase the processing functions above (defined in exp_fields_and_data_types) made empty values
	exp_data_modified = useful.removeEmptyVals(exp_data_modified)					
	
	
	# check if keys in exp_data are a superset of required field keys
	# DO this again at the end to ensure that no required fields are removed in the process above 
	if not set(exp_data_modified.keys()) >= set(required_fields):
		missing_fields = set(required_fields) - set(exp_data_modified.keys())
		raise Exception('The experiment data is missing these fields required by the database schema: {0} '.format(','.join(missing_fields)))
	
	exp_data_modified['EXPERIMENT_NAME'] = exp_data_modified['EXPERIMENT_NAME'].replace(',','').replace(';','').replace('/','').replace('\\','')
	exp_data_modified['PROJECT_NAME'] = exp_data_modified['PROJECT_NAME'].replace(',','').replace(';','').replace('/','').replace('\\','')
		
	#ADD duplicated fields to fields whose case-sensitivity we will not maintain and field that may have spaces. We will search these duplicated fields during queries
	exp_data_modified['DUPLICATED_FIELDS'] ={}
	for indexed_field in set(schema.Exps_Collection()['duplicated_fields']) & set(exp_data_modified.keys()):  # set intersection, make sure we only check fields that do exist in exp_data
		#exp_data_modified[indexed_field + '_INDEX'] = convert_text_to_index_field_text(exp_data_modified[indexed_field])	
		exp_data_modified['DUPLICATED_FIELDS'][indexed_field] = convert_text_to_index_field_text(exp_data_modified[indexed_field])	
	
	#ensure that all users granted write access actually exist in the database
	list_of_db_users = [db_user['user'] for db_user in db.users.find({},{'user':1,'_id':0})]
	admin_emails = ','.join([db_user['name']+' ('+db_user['email']+')' for db_user in db.users.find({'curator':True}) if 'name' in db_user and 'email' in db_user])
	if not set(list_of_db_users) >= set(exp_data_modified['OWNERS_OF_EXPERIMENT']):
		undefined_users = set(exp_data_modified['OWNERS_OF_EXPERIMENT'])-set(list_of_db_users)
		raise Exception('The experiment was not added to the database. The following users granted write access to experiment are not defined in the database: {0}. Please only list allowed database users or request the administrator to add a new user. The following are administrators you can contact: {1}'.format(','.join(undefined_users),admin_emails))
	
	#make sure users with write access also have read access
	if 'READ_ACCESS' in exp_data_modified:
		exp_data_modified['READ_ACCESS'] = list(set(exp_data_modified['READ_ACCESS']+exp_data_modified['OWNERS_OF_EXPERIMENT']))
	else:
		exp_data_modified['READ_ACCESS'] = [o for o in exp_data_modified['OWNERS_OF_EXPERIMENT']]
	
	#ensure that no 'owners of experiments' can own experiments with the same experiment name and project name 
	owner_check = {'OWNERS_OF_EXPERIMENT':{'$in':exp_data_modified['OWNERS_OF_EXPERIMENT']}}
	project_name_check = {'PROJECT_NAME_INDEX':exp_data_modified['PROJECT_NAME_INDEX']} if 'PROJECT_NAME_INDEX' in schema.Exps_Collection() else {'DUPLICATED_FIELDS.PROJECT_NAME':exp_data_modified['DUPLICATED_FIELDS']['PROJECT_NAME']}
	exp_name_check = {'EXPERIMENT_NAME_INDEX':exp_data_modified['EXPERIMENT_NAME_INDEX']} if 'EXPERIMENT_NAME_INDEX' in schema.Exps_Collection() else {'DUPLICATED_FIELDS.EXPERIMENT_NAME':exp_data_modified['DUPLICATED_FIELDS']['EXPERIMENT_NAME']}
	check_project_exp_name = dict(owner_check.items()+project_name_check.items()+exp_name_check.items())			
	if db.exps.find(check_project_exp_name).count()>0:
		raise Exception('We do not allow for owners of experiments to insert experiments that have the same name as other experiments in that project. Please change the name of experiment or project name. Experiment name: {0}, Project name: {1}'.format(exp_data_modified['EXPERIMENT_NAME'],exp_data_modified['PROJECT_NAME']))
	
	#IF pairing technique is defined, then just make 'VH:VL_PAIRED' field true
	if 'PAIRING_TECHINQUE' in exp_data_modified:
		exp_data_modified['VH:VL_PAIRED'] = True
					
	#keep track of who uploaded experiment
	exp_data_modified['UPLOADED_BY'] = user_info['user']
	#keep track of when experiment was added to database
	exp_data_modified['EXPERIMENT_CREATION_DATE'] = time.strftime('%D') #=> use string method because JSON cannot dump datetime data...but bson can...could be problematic though if we always use json loads...maybe not..datetime.datetime.now()#mongo requires time with dates, cannot pass in just date -> datetime.date.today()#time.strftime('%D')
	
	#IF a non-administrator has inserted an experiment, make curated field 'FALSE'
	exp_data_modified['CURATED'] = True if 'curator' in user_info and user_info['curator']==True else False
	#always default seq count to 0
	exp_data_modified['SEQ_COUNT'] = 0				
	#db.exps.insert(exp_data_modified)
	return db.exps.insert(exp_data_modified)


