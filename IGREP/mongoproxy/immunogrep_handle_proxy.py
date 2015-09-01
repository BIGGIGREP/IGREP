import urllib2
import json
import pymongo
from datetime import datetime
import time
from collections import namedtuple
import collections
import __builtin__
from bson.objectid import ObjectId
import json
from itertools import imap
from functools import partial
import inspect
import sys
import socket
import os
import traceback
import ast
#import immunogrep_proxy_tools as proxy_tools
from bson.json_util import dumps as bson_dumps 
from pprint import pprint
#import bennigoetz_appsoma_ab_members as members
import immunogrep_authenticate_ab_members as members
#import bennigoetz_proxytest_query_functions 
#import proxytest_query_functions 	
import immunogrep_db_query_api
import immunogrep_db_writes as bennigoetz_writes
#import cchrysostomou_proxy_test_query_functions 
import  immunogrep_igdbtools as igdbconnect
	
# import bennigoetz_writes

def handle_GET( password, path ):
	return { "command":"GET" }

def handle_POST(password_dict, path, data, wfile):		
	try:
		#TRY TO RUNFUNCTIONS ON PROXY....IF ANY ERRORS ARE ENCOUNTERED, ERRORS ARE REPORTED IN ERROR LOG BELOW (SEE EXCEPT). ERROR LOG IS CALLED "mongoproxy_error_log.txt"
		if 'authkey' not in data:			
			data['authkey'] = ''
		if data['authkey'] == '':
			print('get rid of this authkey hack!!')
		ap_user = members.getFullName(data['authkey'])
		
		db_connect_data = (password_dict['dbpath'],'reader',password_dict['reader'],ap_user) 		
		
		#print datetime.now()				
				
		if not os.path.isdir(ap_user):
			os.mkdir(ap_user)
		#[db_writer,_ig_]= igdbconnect.connectToIgDatabase('writer',password_dict['writer'])
		[db_reader,_ig_]= igdbconnect.connectToIgDatabase('reader',password_dict['reader'])						
		user_dict = db_reader.users.find_one({'user':ap_user},{'_id':0})
				
		if not user_dict: #user not found in database		
			user_dict = {
				'user':ap_user,
				'name':'',
				'email':'',
				'administrator':False,
				'lab':[],
				'write_access': False
			}
			
		if data.get('db_action') == 'updates': #for insert and updating to databse
			print "debug: updates"
			bennigoetz_writes.db = db_writer
			data.pop('db_action', None)
			#[db_writer,_ig_] = igdbconnect.connectToIgDatabase('writer',password_dict['writer'])
			print data['module']
			data['module'] = 'bennigoetz_writes'
			print data['module']
			print data['command']
#			print data['args']
			print data['kwargs']
			print globals()['bennigoetz_writes']
			data['kwargs']['user_info'] = user_dict
	#		getattr(globals()['bennigoetz_writes'], data['command'])(db=db_writer, user_info=user_dict, *data['args'], **data['kwargs'])
			return_value = getattr(globals()['bennigoetz_writes'], data['command'])(*data['args'], **data['kwargs'])
			print "return value: {}".format(return_value)
			wfile.write(return_value)
	#		db_writer.test_exps.insert({'animal': 'walrus'})
	
		elif data.get('db_action') == 'query': #this is what actually processes a request sent to proxy 
			
			#query_id = data['query_object_id']
	
			#if query_id not in query_objects:
			#	query_objects[query_id] = proxytest_query_functions.RunQuery(query_id,password_dict,db_method=data['connection_type'])
									
			#RunQuery Parameters => db_connect_data,db_user_name,modify_query_values_to_follow_db_schema=True,redirect_fields_for_improved_queries=True,to_file=True,file_prefix=None,proxy_path = ''			
			#see _return function to see how the data dictionary is populated
			#important => set proxy_path to None because we are already on the proxy 
			class_args = [db_connect_data]			
			class_karg = {'proxy_path':None}
			if 'to_file' in data:
				class_karg['to_file'] = data['to_file']
			if 'modify_query_values' in data:
				class_karg['modify_query_values_to_follow_db_schema'] = data['modify_query_values']
			if 'redirect_query_fields' in data:
				class_karg['redirect_fields_for_improved_queries'] = data['redirect_query_fields']
			if 'file_prefix' in data:
				class_karg['file_prefix'] = data['file_prefix']						
									
			new_query_object=immunogrep_db_query_api.RunQuery(*class_args,**class_karg)
			
			#new_query_object= immunogrep_db_query_api.RunQuery(db_connect_data,ap_user,proxy_path=None,modify_query_values_to_follow_db_schema=data['modify_query_values'],redirect_fields_for_improved_queries=data['redirect_query_fields'],to_file=data['query_to_file'],file_prefix=data['file_prefix'])
									
			for function_to_run in data['command']: #data['command'] contains a list of functions to run and all of the accompanying parameters passed in by teh user on appsoma side
				#print function_to_run['args']
				#print dumps(function_to_run['args'])
				#print function_to_run
				#print 'args'
				#print str(function_to_run['args'])
				#print function_to_run['kwargs']
				#print 'keyargs'
				#print bson_dumps(function_to_run['kwargs'])
				#print str(function_to_run['kwargs'])
				
				#if we use bson_loads, then parameters such as ObjectId and re.compile are not passed in correcto to the funciton using exec 
				#however, if we jsut use bson_dumps or json_dumps, then variables such as True are converted to true which are not correct 
				#so instead, we => load as json by first dumps using bson and then reload using json, then convert to string
				exec('new_query_object.{0}(*{1},**{2})'.format(function_to_run['command'],str(json.loads(bson_dumps(function_to_run['args']))),str(json.loads(bson_dumps(function_to_run['kwargs'])))))
				
				#getattr(globals()['__builtins__']['new_query_object'],function_to_run['command'])(*function_to_run['args'], **function_to_run['kwargs'])		
								
			#return_query = query_objects.pop(query_id)
			new_query_object._return_results(wfile) #ok all commands have been run, lets return the results via proxy			
			del new_query_object				
			#OK, all done.... now, lets delete me...
			#query_objects.pop(query_id,None)				
		
	except Exception as e: #CATCH any errors						
		exc_type, exc_obj, exc_tb = sys.exc_info()
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]    
		tb_error = traceback.format_exc()		
		
		with open(ap_user+'/'+str(ap_user)+"_mongoproxy_error_log.txt","a") as error_log: #WRITE errors to file 							
			#tb_error = tb_error.replace('\n',';')						
			#tb_error = tb_error.replace('\t',',')					
			error_log.write('###ERROR AT: '+str(datetime.now())+'\n')
			error_log.write('Username: '+str(ap_user)+'\n')
			error_log.write('Proxy Request: '+bson_dumps(data,indent=4, separators=(',', ': '))+'\n')
			error_log.write('ERROR Message: '+tb_error)
			error_log.write('########################################\n\n\n')			
			#error_log.write(str(datetime.now())+'\t'+str(ap_user)+'\t'+json.dumps(data)+'\t'+tb_error+'\n')
		
		raise Exception(tb_error)
	
	
	
	