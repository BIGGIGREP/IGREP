import pymongo
import os
import re
import sys
import copy
import sets
import datetime
from array import *
from collections import namedtuple

# use simplejson if available
# else use json
import json
#try: import simplejson as json
#except ImportError: import json
from bson.json_util import dumps as bson_dumps 
#-----------------------------------------------------------------------------------------------------------------------------------------
# ig_database_tools.py
#
# 	Contains the utility functions for
#		- connecting to MongoDB at TACC
#		- adding FASTA file data to MongoDB
#		- adding data from IMGT files to MongoDB (configurable to allow changing which analysis data is used in the database)
#-----------------------------------------------------------------------------------------------------------------------------------------


def ig_DB_PATH(whichserver='biotseq.icmb.utexas.edu',port=27017):

	#return [None, None]
	
	connection = pymongo.MongoClient(host=whichserver,port=port)	
	db = connection.appsoma	

	return [db,connection]


#-----------------------------------------------------------------------------------------------------------------------------------------
# MongoDB connection
#-----------------------------------------------------------------------------------------------------------------------------------------
# opens the database connection.
def connectToIgDatabase(version='writer',Ig_DB_passwd='rag1rag2',path=None):
	print "aaaargconnecting4"
	#whichserver=int(raw_input("Which database 0=TACC 1=BIGG 2=BIOTSEQ:"))
		
	connection_data_template = {
	'connection': None,
	'db': None,
	'db_seqs': None,
	'db_info': None,
	'exps':None,
	'speciesid':None,
	'motifs':None
	}
	

	return_values = connection_data_template.copy()
	if path:
		[db,connection] = ig_DB_PATH(path)
	else:
		[db,connection] = ig_DB_PATH()
		
		
	
	db.authenticate(version,Ig_DB_passwd)
	# TODO: these two variables don't make accessing collections any more convenient -
	# these should be removed and references to them replaced with db.{collection_name}
	#db_sequences = db.seqs
	
	
	return_values['connection'] = connection
	return_values['db'] = db
	return_values['seqs'] = db.seqs
	return_values['exps'] = db.exps
	return_values['exps_test'] = db.test_exps
	return_values['seqs_test'] = db.test_seqs
	#return_values['speciesid'] = db.SpeciesID
	return_values['motifs'] = db.Motifs
	return_values['users'] = db.users
	#return_values['seqsarray'] = db.seqs_array_v2
	return_values['seqsarray'] = db.test_seqs
	#return_values['exp_old'] = db.exps_old
	#print "aheraek"

	db_path = return_values['db']
	db_sequences = return_values['seqs']
				
	#to use namedtuple function, must import collections namedtuple (from collections import namedtuple)
	#initialize the tuple
	db_info = namedtuple('db_info','seqsarray,seqs,exps,motifs,users')
			
	#create the tuple using dbTuple
	#db = db_info(seqsarray=return_values['seqsarray'],sequences=return_values['seqs'],experiments=return_values['exps'],motifs=return_values['motifs'],species = return_values['speciesid'],seqreal=return_values['seqs_real'],expreal=return_values['exps_real'],users=return_values['users'],exp_old=return_values['exp_old']);# {'seqs':connectiondata['seqs'],'exps':connectiondata['exps']}
	db = db_info(seqsarray=return_values['seqsarray'],seqs=return_values['seqs'],exps=return_values['exps'],motifs=return_values['motifs'],users=return_values['users']);# {'seqs':connectiondata['seqs'],'exps':connectiondata['exps']}
	#print (str(return_values['seqs']))
	# returns the connection and the database.
	return db,return_values['connection']

