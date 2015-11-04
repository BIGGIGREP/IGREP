"""
Query API made for interacting with the igrep mongodb database

10-4-2015:
	Seperated the api into two parts:

		Main part - Class name = RunQuery (this file/module): a method for interacting with the immunogrep mongodb database. This is the main class and will assume you have credentials for accessing
		the exps and seqs collection.

		Part two - Class name = QueryProxy (new file/module immunogrep_query_dbproxy_api): a proxy class for calling the main class above (Main part). This is for users you do not know the password for access to each collection,
		but have experiments in the database or want to query publicly available data. They just would need to know the address of the proxy to the database
		and optionally their login credentials (if they have their own experiments on that database)

	Current design for a user who has access to the seqs and exps collection (usually an administrator):

		User -> RunQuery Class (API for optimizing/standardizing queries to seqs and exps collection) -> Pymongo -> Immunogrep DB using MongoDB
		RunQuery Class will require:

			1) MongoDB database address
			2) MongoDB Username/password authentication to the seqs/exps collections

	Current design for a user who wants to query database without access (non-administrator):
		User -> QueryProxy Class -> Proxy to controlling access to database -> RunQuery Class -> Pymongo -> Immunogrep DB using MongoDB

6/25/2015:

	new version for querying:

		allows conversions to CSV, TAB, JSON, IGREP,
		allows queries across multiple anotation methods
		allows user to query the NGS sequence and sequence header for any seqs collection query

"""

import inspect
import sys
import json
import random
from functools import wraps
from datetime import datetime
from bson.json_util import dumps as bson_dumps
from bson.json_util import loads as bson_loads
from bson.objectid import ObjectId
import copy
import os
import time
import itertools
from collections import defaultdict
import re
from collections import OrderedDict

# immunogrep modules:
# variables/characters we will use

from immunogrep_global_variables import idIdentifier  # currently idIdentifier = 'SEQ_ID', but check global_variables for changes
from immunogrep_global_variables import expIdentifier  # currently expIdentifier = 'EXP_ID' , but check global_variables for changes
from immunogrep_global_variables import seqRawData  # currently seqRawData = '@SEQ'
# currently fasta_file_delimiter = '< '. We use this for seperating a fasta header, from our internal data. For example, here is a possible FASTQ/FASTA header with JSON data:
# >sequence_1< {'SEQ_ID':ObjectId(),'PROJECT_NAME':'demo'}; or @sequence_1< {'SEQ_ID':ObjectId(),'PROJECT_NAME':'demo'}
from immunogrep_global_variables import fasta_file_delimiter

import immunogrep_database_schema as schema
import immunogrep_useful_functions as useful

from mongo_igdb_tools import connectToIgDatabase
from mongo_igdb_tools import convert_to_objectid
from mongo_igdb_tools import convert_text_to_index_field_text
from mongo_igdb_tools import default_fields_data_types
from mongo_igdb_tools import AttemptToConvertToBool
from mongo_igdb_tools import default_metadata_fields

try:
	# cython-ized version of the functions from immunogrep_useful functions. These versions are about 2X faster
	from immunogrep_cython_db_tools import flatten_dictionary
	from immunogrep_cython_db_tools import RemoveObjId
except:
	# it wont work if there is no cython module (currently immunogrep_cython_db_tools). when it doesnt work, import non cyhton version
	from immunogrep_useful_functions import flatten_dictionary
	from immunogrep_useful_functions import RemoveObjId
	print('Not using cython')

# NOTE1 :PYMONGO AGGREGATION SORT REQUIRES BSON OR ORDEREDDICT! -> http://api.mongodb.org/python/current/examples/aggregation.html

# ######GLOBAL VARIABLE DEFINITIONS###############
'''
	Defining global variables used in module
'''
oid_type = type(ObjectId())

# these variables assume a default structure for our SEQS collection
# We assume that any AB data with regard to a specific NGS read (SEQ_ID) is included as a subdocument under the key called 'DATA' (DATA.CDR3.NT).
# The following list defines fields that are NOT fields within DATA
fields_above_data_key = ['_id', 'ANALYSIS_NAME', 'EXP_ID', 'RECOMBINATION_TYPE', 'SEQ_ID', 'DATE_UPDATED', 'SETTINGS']

# The following lists defines fields within the DATA key that exist for raw NGS data (not AB annotated data; analysis_name = "@SEQ")
data_fields_for_raw_data = ['SEQUENCE_HEADER', 'SEQUENCE', 'QUALITY_SCORE', 'FILENAME']

# This key name is a hardcode key that exists in the SEQS collection. It is used to store information in a document that IS NOT wanted by the user,
# but is added in by the us when inserting sequence data to the database. It is meant to only be used for improving queries on specific fields.
# It should barely be used in projections. It should also rarely be used as a field in queries, because this API is meant to re-route any fields
# automatically. For example, any queries on DATA.VREGION.VGENES will be re-routed to QUERY_DATA.VREGION.VGENES.PARSED_ALLELES (see redirect_seq_collection_fields)
query_data_key_name = 'QUERY_DATA'

# These next three variables are currently untested and not implemented
# Initial tests on Biotseq showed this method was not faster than running single processes
# When true, the results in a cursor from a specific query are loaded into lists and processed in parallel.
allow_multithreading = False
max_multithreading_processes = 6
chunk_size = 1


# a harcoded variable controlling the preferred order of fields output to delimited files (that is, if the user does not define an order, then these fields will always appear first)
# also shows the default schema as of 6/26/2015
default_sorting_order = [
	idIdentifier,
	'SEQUENCE_HEADER', 'SEQUENCE', 'QUALITY_SCORE',
	'PREDICTED_AB_SEQ.AA', 'PREDICTED_AB_SEQ.NT', 'STRAND', 'PRODUCTIVE',
	'RECOMBINATION_TYPE', 'PREDICTED_CHAIN_TYPE', 'LOCUS_NAME',
	'VREGION.VGENES', 'DREGION.DGENES', 'JREGION.JGENES',
	'FULL_LENGTH', 'STOP_CODONS',
	'VREGION.FR1.AA', 'VREGION.CDR1.AA', 'VREGION.FR2.AA', 'VREGION.CDR2.AA', 'VREGION.FR3.AA', 'CDR3.AA', 'JREGION.FR4.AA',
	'VREGION.FR1.NT', 'VREGION.CDR1.NT', 'VREGION.FR2.NT', 'VREGION.CDR2.NT', 'VREGION.FR3.NT', 'CDR3.NT', 'JREGION.FR4.NT',
	'VREGION.CDR1.AA_LENGTH', 'VREGION.CDR1.NT_LENGTH', 'VREGION.CDR2.AA_LENGTH', 'VREGION.CDR2.NT_LENGTH', 'CDR3.AA_LENGTH', 'CDR3.NT_LENGTH',
	'VREGION.VGENE_SCORES', 'JREGION.JGENE_SCORES', 'DREGION.DGENE_SCORES',
	'VREGION.SHM.NT', 'VREGION.SHM.AA',
	'GAPPED.VREGION.FR1.AA', 'GAPPED.VREGION.CDR1.AA', 'GAPPED.VREGION.FR2.AA', 'GAPPED.VREGION.CDR2.AA', 'GAPPED.VREGION.FR3.AA', 'GAPPED.CDR3.AA', 'GAPPED.JREGION.FR4.AA',
	'GAPPED.VREGION.FR1.NT', 'GAPPED.VREGION.CDR1.NT', 'GAPPED.VREGION.FR2.NT', 'GAPPED.VREGION.CDR2.NT', 'GAPPED.VREGION.FR3.NT', 'GAPPED.CDR3.NT', 'GAPPED.JREGION.FR4.NT',
	'VREGION.VGENE_QUERY_START', 'VREGION.VGENE_QUERY_END',
	'JREGION.JGENE_QUERY_START', 'JREGION.JGENE_QUERY_END',
	'NOTES',
	expIdentifier, 'ANALYSIS_NAME', 'DATE_UPDATED', 'FILENAME', 'SETINGS'
]


# ###########END OF GLOBAL VARIABLES DEFINITIONS ########################

# GLOBAL FUNCTIONS USED BY THE QUERY CLASS #################3

def get_allowed_file_types():
	'''
		Results from the query can be saved in many different file formats. This summarizes each of the file formats we currently support.

		.. note::

			The word 'flattened' refers to how we handle documents returned from the database. Documents from cursors are returned as dictionaries.
			But most documents are nested using JSON format. For example DATA.JGENES would be a nested dictionary of:
				{'DATA':
					'JGENES':[]
				}
			So when we refer to flattened dictionaries as those whose dictionaries have only one level of keys. All sub-documents are seperated using
			'.' as a string in the key name. So the example above becomes:
				{'DATA.JGENES':[]}

		..TAB::

			TAB delimited file.
			Results from a database query will always be flattened prior to saving to file.
			All field value types will be converted to strings. Lists are converted to strings using ','.join(x)

		..CSV::

			Comma delimited file.
			Results from a database query will always be flattened prior to saving to file.
			All field value types will be converted to strings. Lists are converted to strings using '|'.join(x)

		..FLATJSON::

			Our custom JSON format where ever line in a file represents a 'flattened' document in JSON format.

		..JSON::

			Our custom JSON forma where each line in a file represents a document in JSON format.
			Does not modify the results from a query. Only json.dumps is used to save results to file.
			Field value types are maintained (i.e. lists remain lists)

		..FASTA::

			FASTA file format.
			Results from a database query will always be flattened prior to saving to file.
			All field value types will be converted to strings.
			Any results that are not defined to be part of the sequence header will be reported in json format using the 'fasta_file_delimiter' field
			(i.e. >SEQHEADER <{'SEQ_ID':ObjectId()}

		..FASTQ::

			FASTQ file format.
			Results from a database query will always be flattened prior to saving to file.
			All field value types will be converted to strings.
			For results that do not have a QUALITY_SCORE, we will output a default quality value for each character
			Any results that are not defined to be part of the sequence header will be reported in json format using the 'fasta_file_delimiter' field
			(i.e. @SEQHEADER <{'SEQ_ID':ObjectId()}
	'''
	allowed_file_types = ['TAB', 'JSON', 'IGREP', 'CSV', 'FASTA', 'FASTQ']
	return allowed_file_types


# ###### FUNCTIONS AND VARIABLES USED FOR MODIFYING QUERIES ###############
'''
	This section defines a series of variables and functions that will be used when modifying user queries.

	Some variables are used to redefine field names and query data types. It is used alongside the function process_field_value
	These variables define what are mongo-operators as compared to user defined field names, how to rename field names, and how to change the datatypes of select queires

	Flow of functions will go as follows:
	Key function, Parse_Mongo_Query_Expression -> uses the parameter redirection_fields and calls the function Process_Field_Value
		The parameter, redirection_fields will either be set to:
			(a) redirect_seq_collection_fields OR
			(b) redirect_exp_collection_fields
	Key function Process_Field_Value -> uses the parameter dict_defining_value_transformation
		The parameter dict_defining_value_transformation will either be set to:
			(a) fields_for_queries_seqs_collection (defined in modifying_seqs_collection_queries) OR
			(b) fields_for_queries_exps_collection

'''

# FIRST LETS DEFINE HOW WE WANT TO ADDRESS/HANDLE QUERIES ON SPECIFIC FIELD NAMES

# Any time a user wants to perform a query on the seqs collection using any of these fieldnames defined as keys, we will change the field name to the key's value
# For example, DATA.VREGION.VGENES will be convereted to QUERY_DATA.VREGION.VGENES.PARSED_ALLELES
# USAGE: new_field = redirect_seq_collection_fields['DATA.VREGION.VGENES'] ==> new_field is now equal to: 'QUERY_DATA.VREGION.VGENES.PARSED_ALLELES
redirect_seq_collection_fields = {
	'DATA.VREGION.VGENES': lambda x: 'QUERY_' + x + '.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.' + x + '.PARSED_ALLELES',
	'DATA.DREGION.DGENES': lambda x: 'QUERY_' + x + '.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.' + x + '.PARSED_ALLELES',
	'DATA.JREGION.JGENES': lambda x: 'QUERY_' + x + '.PARSED_ALLELES' if x.startswith('DATA') else 'QUERY_DATA.' + x + '.PARSED_ALLELES',
}

# lets say a user wants to query for an experiment called "MyFirstExp". Unfortunately, the user actually called their experiment "My first exp".
# We want to make sure that the user can still find his/her proper experiment without having to remember to use 'regular expressions' or case insensitivity.
# This is where this variable is useful. When creating experiments, we create duplicated fields of select fields in the experiments collection.
# Any time a user wants to perform a query on the exps collection using any of the fields we know are duplicated (defined as keys in this variable), then we will change the field name to the key's value
# For example, a user may request a query on EXPERIMENT_NAME : 'MyFirstExp'. This query will be converted to 'DUPLICATED_FIELDS.EXPERIMENT_NAME':'MyFirstExp'. This will ensure that his/her query is found

metadata_duplicated_fields = schema.Exps_Collection()['duplicated_fields']
redirect_exp_collection_fields = {field_name: 'DUPLICATED_FIELDS.' + field_name for field_name in metadata_duplicated_fields}


def Process_Field_Value(db_field_name, values, dict_defining_value_transformation):
	'''
		This function converts the queries provided for specific database field names into the format defined by the dictionary dict_defining_value_transformation

		Parameters
		----------
		db_field_name : The field name of the database query. The field name should be in '.' notation (i.e. DATA.VREGION.VGENES)
		values : The specific query requested for that field name. Allowed formats: are str, list, int, float or lists of those formats
		dict_defining_value_transformation : A key-value pair defining how to transform values for each db_field_name. Unless modified by user, this variable will be equal to either fields_for_queries_seqs_collection or fields_for_queries_exps_collection (see functions below)

	'''

	[db_field_name, outer_field] = Remove_Numbers_From_Field_Name(db_field_name)

	if type(values) is list:
		temp = []
		for list_element in values:
			try:  # use a try/except, because we assume that teh user may have a better knowledge of the query. so if it should be an int, but now its a string, then default to what user passed in. worst case scenario, the query will return no results
				# always ignore types that are REGULAR EXPRESSIONS. Cannot compile those to match datatype
				if isinstance(list_element, type(re.compile(''))) or isinstance(list_element, type(bson_loads(bson_dumps(re.compile(''))))) or isinstance(list_element, type(bson_loads(bson_dumps({'$regex': ''})))):
					# this means that the value is a regular expression.
					# in order to convert a re format from python object to a dictionary, do the following:
					re_dict_format = json.loads(bson_dumps(list_element))  # first dump the object to a string using bson, then load it as a dic using json
					# now we can actually modify the string that has been made into a regular expression (it is found in '$regex key')
					re_dict_format['$regex'] = dict_defining_value_transformation[db_field_name](re_dict_format['$regex'])
					# now that we have formated the variable correctly, lets update v, but convert  this dict back into python re format
					v = bson_loads(bson_dumps(re_dict_format))
				else:
					v = dict_defining_value_transformation[db_field_name](list_element)
			except Exception as e:
				print "Error in compiling fields values: {0}. Error {1}".format(str(values), str(e))
				v = list_element
			temp.append(v)
	else:
		try:  # use a try/except, because we assume that the user may have a better knowledge of the query. so if it should be an int, but now its a string, then default to what user passed in. worst case scenario, the query will return no results
			# always ignore types that are REGULAR EXPRESSIONS. Cannot compile those to match datatype
			if isinstance(values, type(re.compile(''))) or isinstance(values, type(bson_loads(bson_dumps(re.compile(''))))) or isinstance(values, type(bson_loads(bson_dumps({'$regex': ''})))):
					# this means that the value is a regular expression.
					# in order to convert a re format from python object to a dictionary, do the following:
					re_dict_format = json.loads(bson_dumps(values))  # first dump the object to a string using bson, then load it as a dic using json
					# now we can actually modify the string that has been made into a regular expression (it is found in '$regex key'
					re_dict_format['$regex'] = dict_defining_value_transformation[db_field_name](re_dict_format['$regex'])

					# now that we have formated the variable correctly, lets update v, but convert  this dict back into python re format
					temp = bson_loads(bson_dumps(re_dict_format))
			else:
				temp = dict_defining_value_transformation[db_field_name](values)
		except Exception as e:
			print "Error in compiling fields values: {0}. Error {1}".format(str(values), str(e))
			temp = values
	return temp


def Modify_Key_Name_To_Include_Data(key):
	'''
		Adds a 'DATA' string before a field name. i.e. 'VREGION.VGENES' becomes 'DATA.VREGION.VGENES'
	'''
	if key not in fields_above_data_key:
		first_field = key.split('.')[0]
		if first_field != query_data_key_name and first_field != 'DATA':
			key = 'DATA.' + key
	return key


def Parse_Mongo_Query_Expression(query, field_name='', dict_defining_value_transformation={}, redirection_fields={}, modify_query_values_to_follow_db_schema=True, redirect_fields_for_improved_queries=True):
	'''
		This function will attempt to ensure that the values for each field name are properly formatted.
		It will redirect field names to make sure they use indexed fields. It will ensure that the values for each query are in the proper format.

		For example: It will always make sure queries on CDR3.AA are uppercased, and that any queries involving DATA.VREGION.VGENES to QUERY_DATA.VREGION.VGENES.PARSED_ALLELES

		Parameters
		----------
		query : dict
			This should be a dictionary used as a MongoDB query
		field_name : string
			name of the current 'parent' document. Should start as empty string when first calling function.
		dict_defining_value_transformation : dict
			A key-value pair defining how to transform values for each db_field_name. Unless modified by user, this variable will be equal to either fields_for_queries_seqs_collection or fields_for_queries_exps_collection (see functions below)
		redirection_fields : dict
			key-value pair where fields/keys defined in this dictionary will be redirected to their values. If not defined by user then it would be equal to redirect_seq_collection_fields or redirect_exp_collection_fields
		modify_query_values_to_follow_db_schema : boolean
			If True, then it will modify the values in a database query to match those defined by dict_defining_value_transformation
		redirect_fields_for_improved_queries : boolean
			If True, then change/rename the field names of the query using the values defined by redirection_fields

	'''

	# do not change the values to any keys in a query in the following list. these values are standard numbers/booleans based on a mongodb query and not the field it self (for example fields that are strings still will use True/False for $exists)
	operators_do_not_modify = ['$size', '$exists', '$mod', '$options', '$type', '$meta', '$language', '$geoWithin', '$geoIntersects', '$near', '$minDistance', '$maxDistance', '$geometry', '$nearSphere']

	# 1) we do not use a 'text' index in our schema currently, so we are not handling how to parse this query; without a 'text index' on a specific field then any query that defines this field will raise an error in mongodb
	# 2) the '$where' command uses javascript functions to perform a query. again we will not parse these, and leave these special cases to the user
	# 3) comment mongo operator is just for adding comments to a query
	operators_do_not_parse = ['$text', '$where', '$comment']

	def RunParser(query, field_name, dict_defining_value_transformation, redirection_fields, modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries):
		ignore_mongo_operators = operators_do_not_modify + operators_do_not_parse
		if type(query) is dict:
			copied_query = {f: v for f, v in query.iteritems()}
			# loop through all keys in dictionary
			for key, values in copied_query.iteritems():
				# rename keys that should be redirected. this will only work if the field is explicitly written in dot notation. For example, it will not work if its CDR3:{'AA':'APPL'}. we cannot redirect this field to  QUERY_DATA.CDR3.AA
				if redirect_fields_for_improved_queries and key[0] != '$':
					# key = Modify_Key_Name_To_Include_Data(key) ===> was attempting to add 'DATA' key to queries, but this would make confusing problems with exp collection vs seqs collection

					[key_no_int, outer_key] = Remove_Numbers_From_Field_Name(key)
					if key_no_int in redirection_fields:
						popped_value = query.pop(key)
						key = redirection_fields[key_no_int](key)
						query[key] = popped_value

				if key in ignore_mongo_operators:  # operators_do_not_modify:
					# do not change the values passed in for this field
					# i.e. {'field_name_a':{'$exists':True}} => even if field name is a string, obviously the value to exists is only True/False, so do not modify the value of this field
					query[key] = values  # just reutrn the original query
					continue
				elif key[0] == '$':  # key is a specific mongo operator
					sub_field_name = field_name  # sub_field_name stays the same, because a specific mongo operator does not add to a field name, just the mongo function
				else:
					# its not a mongo operator, so it must be a field name, add this as a subdocument in the query command using '.' notation
					# for exampple if we are querying: {CDR3:{AA:'ALPHA','NT':'BETA'}} then we will want to LOOK AT THE FOLLOWING FIELD NAMES: CDR3.AA, CDR3.NT
					sub_field_name = field_name + '.' + key if field_name != '' else key

				if type(values) is list:
					for i, each_value in enumerate(values):
						# if its a list, then recursively parse through each value of list
						query[key][i] = RunParser(each_value, sub_field_name, dict_defining_value_transformation, redirection_fields, modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries)
				elif type(values) is dict:
					# now that we always load dictionary using bson in above/parent fucntion, then the following if/elif statements should
					# not occur. there should no longer be keys that say {'$oid'] or {'$regex'}. instead, they exist as ObjectId and re.copmile
					# hardcoding object id conversions...
					# no need to check recursifely this dictionary. it is an object id
					if values.keys() == ['$oid']:
						query[key] = Process_Field_Value(sub_field_name, values, dict_defining_value_transformation)
					elif '$regex' in values:
						# hardocindg regular expressions
						# again the user requested using regular expressions , so no need to recursively check this dictionary.
						# first process the regular epxression using process_field_value (only process '$regex' key
						# once complete, use bson_loads to convert it into a python re variable
						values['$regex'] = Process_Field_Value(sub_field_name, values['$regex'], dict_defining_value_transformation)
						# first uson json.dump to dump dictionary, then RELOAD IT using bson
						query[key] = bson_loads(json.dumps(values))
					else:
						# if its a dict, then just recursively pass in the dictionary rather than each index  in a list as above
						query[key] = RunParser(values, sub_field_name, dict_defining_value_transformation, redirection_fields, modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries)
				else:  # assume anythign else is a specific query parameter/value. therefore, prcoess its field value
					if modify_query_values_to_follow_db_schema:  # modify the value for this query
						query[key] = Process_Field_Value(sub_field_name, values, dict_defining_value_transformation)

		else:  # if the query was not actually a dictionary then it is probably from values in a list, or just singular values, so just process each field name
			values = query
			if modify_query_values_to_follow_db_schema:  # modify teh value for this query
				query = Process_Field_Value(field_name, values, dict_defining_value_transformation)
		return query

	# first take the provided query and ENSURE that it will be in BSON format:
	# dump query as a string using bson dumps => handles ObjectIds, and regular expression values
	# once dumped, reload is as a bson to ensure that any previous dictionaries which should be objects (such as object id and regex) are represtend properly
	query = bson_loads(bson_dumps(query))

	# now that we have an ensured BSON format dictoinary, run the parsing function if user wants to modify field names or values
	if modify_query_values_to_follow_db_schema or redirect_fields_for_improved_queries:
		query = RunParser(query, field_name, dict_defining_value_transformation, redirection_fields, modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries)

	return query


def Remove_Numbers_From_Field_Name(db_field_name):
	'''
		This function is for handling MULTIKEY STRUCTURES IN MONGODB
		This will parse a field name provided in a query. If the field contains numbers, then it will remove the number.
		For example: 'VREGION.VGENES.1':'IGHV3-2' will be changed to 'VREGION.VGENES':'IGHV3-2'
		We remove the number because when modifying/redirecting a query to QUERY_DATA, we are not interested in the specific array element,
		but instead, we are interested in the main field name 'VREGION.VGENES'

		Parameters
		----------
		db_field_name : string
			This variable corresponds to the field name in the query

		Outputs
		------
		db_field_name : string
			Returns fieldname without number (i.e. VREGION.VGENE.1.AB returns VREGION.VGENE.AB)
		outer_fields : string
			Returns fieldname before the first number was encountered (i.e. VREGION.VGENE.1.AB returns VREGION.VGENE )

	'''
	# first split the db_field_name by '.'
	db_sub_fields = db_field_name.split('.')

	outer_field = []
	reached_integer = False
	db_field_name = []

	for f in db_sub_fields:
		try:
			# test if integer
			int(f)
			if int(f) == float(f):
				reached_integer = True
		except:
			# if its not a number, then it must be a string so , this is part of the field name we care about
			if not(reached_integer):
				# keep track of any field before a number is reached
				outer_field.append(str(f))
			# keep track of string where numbers are removed
			db_field_name.append(str(f))
	db_field_name = '.'.join(db_field_name)
	outer_field = '.'.join(outer_field)
	# outer_fields => fields name before the first number was encountered
	# db_field_name = > fieldname wihtout number
	return [db_field_name, outer_field]


def modifying_seqs_collection_queries():
	'''
		Outputs
		-------
		This function just returns a global variable that is used for handling user queries on the SEQS collection.

		We need to ensure that the queries on the seqs collection are in the proper data types for select fields.
		For example, we need to ensure that queries on VREGION.VGENE_SCORES are numbers, whereas queries on VREGION.VGENES are strings that are uppercased.

		This variable is meant to greatly help the user in getting the proper query result returned. It will define how to modify the value of any field
		quieried by the user. For each defined field in a query, the datatype of the query will be transformed to match the datatype of the requested field.

		Examples:
			This variable can be used to handle the following query:
				{DATA.VREGION.VGENES:'3','DATA.CDR3.NT':'aacg'}
			will be converted to:
				{DATA.VREGION.VGENES:3,'DATA.CDR3.NT':'AACG'}

		The structure of this variable is as follows: dict
			key: the field name whose datatype must be modified
			value: a lambda function defining how to modify the value in the query

			If a key is not defined in this variable, then those values will be handled by the function: default_fields_data_types

	'''

	# these keys from mongo-db will usually be reported as an ObjectId.
	keys_with_object_id = ['_id', 'SEQ_ID', 'EXP_ID']
	# this will modify '_id', 'EXP_ID', idIdentifier (anything that needs to be an object id)
	fields_that_should_query_using_objectids = {oid_field: convert_to_objectid for oid_field in keys_with_object_id}

	# this will modify any queries on the following fields by making them uppercase
	fields_with_strings = [
		'DATE_UPDATED', 'ANALYSIS_NAME', 'RECOMBINATION_TYPE',
		'QUERY_DATA.VREGION.VGENES.PARSED_ALLELES', 'QUERY_DATA.DREGION.DGENES.PARSED_ALLELES', 'QUERY_DATA.JREGION.JGENES.PARSED_ALLELES',
		'DATA.COMMAND', 'DATA.NOTES',
		'DATA.PREDICTED_AB_SEQ.NT', 'DATA.PREDICTED_AB_SEQ.AA', 'DATA.STRAND', 'DATA.PREDICTED_CHAIN_TYPE', 'DATA.PRODUCTIVE', 'DATA.LOCUS',
		'DATA.VREGION.FR1.NT', 'DATA.VREGION.FR1.AA', 'DATA.VREGION.CDR1.NT', 'DATA.VREGION.CDR1.AA', 'DATA.VREGION.FR2.NT', 'DATA.VREGION.FR2.AA',
		'DATA.VREGION.CDR2.NT', 'DATA.VREGION.CDR2.AA', 'DATA.VREGION.FR3.NT', 'DATA.VREGION.FR3.AA', 'DATA.CDR3.NT', 'DATA.CDR3.AA', 'DATA.JREGION.FR4.NT', 'DATA.JREGION.FR4.AA',
		'DATA.VREGION.VGENES', 'DATA.JREGION.JGENES', 'DATA.DREGION.DGENES', 'DATA.ISOTYPE.GENE'

	]
	fields_that_should_query_using_uppercase_strings = {field: lambda x: x.upper() for field in fields_with_strings}

	# this will modify any queries on the following fields by making them into floats
	fields_with_floats = [
		'DATA.VREGION.SHM.NT', 'DATA.VREGION.SHM.AA', 'DATA.VREGION.SHM.NT_PER', 'DATA.VREGION.SHM.AA_PER',
		'DATA.JREGION.SHM.NT', 'DATA.JREGION.SHM.AA', 'DATA.JREGION.SHM.NT_PER', 'DATA.JREGION.SHM.AA_PER',
		'DATA.VREGION.VGENE_SCORES', 'DATA.JREGION.JGENE_SCORES', 'DATA.DREGION.DGENE_SCORES',
		'DATA.ISOTYPE.SCORES', 'DATA.ISOTYPE.PER_ID'
	]
	fields_that_should_query_using_floats = {field: lambda x: float(x) for field in fields_with_floats}	

	# this will modify any queries on the following fields by making them into int
	fields_with_ints = [
		'DATA.VREGION.VGENE_QUERY_START', 'DATA.VREGION.VGENE_QUERY_END',
		'DATA.JREGION.JGENE_QUERY_START', 'DATA.JREGION.JGENE_QUERY_END',
		'DATA.CDR3.NT_LENGTH', 'DATA.VREGION.CDR1.NT_LENGTH', 'DATA.VREGION.CDR2.NT_LENGTH',
		'DATA.CDR3.AA_LENGTH', 'DATA.VREGION.CDR1.AA_LENGTH', 'DATA.VREGION.CDR2.AA_LENGTH',
		'DATA.ISOTYPE.MISMATCHES'
	]
	fields_that_should_query_using_int = {field: lambda x: int(float(x)) for field in fields_with_ints}

	# combine all defined fields into a global variable; any fields not defined in these lists will be treated by the defaultdict function: schema.default_fields_data_types
	fields_for_queries_seqs_collection = defaultdict(lambda: default_fields_data_types)
	fields_for_queries_seqs_collection.update(fields_that_should_query_using_objectids)
	fields_for_queries_seqs_collection.update(fields_that_should_query_using_uppercase_strings)
	fields_for_queries_seqs_collection.update(fields_that_should_query_using_floats)
	fields_for_queries_seqs_collection.update(fields_that_should_query_using_int)

	return fields_for_queries_seqs_collection


def modifying_exps_collection_queries():
	'''
		outputs
		-------
		This function just returns a global variable that is used for handling user queries on the EXPS collection.

		We need to ensure that the queries on the exps collection are in the proper data types for select fields.

		This variable is meant to greatly help the user in getting the proper query result returned. It will define how to modify the value of any field
		quieried by the user. For each defined field in a query, the datatype of the query will be transformed to match the datatype of the requested field.

		Examples:
			This variable can be used to handle the following query:
				{'DUPLICATED_FIELDS.EXPERIMENT_NAME':'hOWDY THere'}
			will be converted to:
				{'DUPLICATED_FIELDS.EXPERIMENT_NAME':'howdythere'}

		The structure of this variable is as follows: dict
			key: the field name whose datatype must be modified
			value: a lambda function defining how to modify the value in the query

			If a key is not defined in this variable, then those values will be handled by the function: default_fields_data_types
	'''

	# these keys from mongo-db will usually be reported as an ObjectId.
	keys_with_object_id = ['_id']
	# this will modify '_id', 'EXP_ID', idIdentifier (anything that needs to be an object id)
	fields_that_should_query_using_objectids = {oid_field: convert_to_objectid for oid_field in keys_with_object_id}

	# this will modify any queries on the following fields by making them lowercase
	fields_to_lower = [
		'SEQUENCING_PLATFORM', 'LAB', 'OWNERS_OF_EXPERIMENT', 'READ_ACCESS',
		'CHAIN_TYPES_SEQUENCED', 'CELL_MARKERS_USED', 'LIST_OF_POLYMERASES_USED', 'PRIMER_SET_NAME',
	]
	fields_that_should_query_using_lowercase_strings = {field: lambda x: x.lower() for field in fields_to_lower}

	# this will modify any queries on the following fields by removing spaces and then lowercasign them
	fields_to_lower_no_spaces = [
		'MID_TAG'
	]
	fields_that_should_query_using_lowercase_strings_no_spaces = {field: lambda x: x.replace(' ', '').lower() for field in fields_to_lower_no_spaces}

	# this will handle how to modify a species name query
	species_query = {'SPECIES': lambda x: ' '.join([word.title() if i == 0 else word.lower() for i, word in enumerate(x.split(' '))])}  # seperate value by spaces, then only capitalize the first word and lowercase all other words

	# this will moidfy any queries on the following fields to booleans
	fields_to_boolean = [
		'CONTAINS_RNA_SEQ_DATA', 'VH:VL_PAIRED', 'POST_SEQUENCING_PROCESSING:PHI_X_FILTER', 'REVERSE_PRIMER_USED_IN_RT_STEP'
	]
	fields_that_should_query_using_boolean = {field: AttemptToConvertToBool for field in fields_to_boolean}

	# this will handle how a name should be queried
	names = {'PERSON_WHO_PREPARED_LIBRARY': lambda x: str(x).title()}

	# this will modify any queries on the following fields by making them into int
	fields_with_ints = [
		'CELL_NUMBER', 'TARGET_READS'
	]
	# should we consider using locale???
	fields_that_should_query_using_int = {field: lambda x: int(float(str(x).replace(',', ''))) for field in fields_with_ints}

	# this will modify any queries on the fields that are present within the DUPLICATED_FIELDS subdocument of the exps metadata
	fields_that_are_duplicated = {'DUPLICATED_FIELDS.' + field_name: convert_text_to_index_field_text for field_name in schema.Exps_Collection()['duplicated_fields']}

	fields_for_queries_exps_collection = defaultdict(lambda: default_metadata_fields)
	fields_for_queries_exps_collection.update(fields_that_are_duplicated)
	fields_for_queries_exps_collection.update(fields_that_should_query_using_objectids)
	fields_for_queries_exps_collection.update(fields_that_should_query_using_int)
	fields_for_queries_exps_collection.update(fields_that_should_query_using_boolean)
	fields_for_queries_exps_collection.update(names)
	fields_for_queries_exps_collection.update(species_query)
	fields_for_queries_exps_collection.update(fields_that_should_query_using_lowercase_strings_no_spaces)
	fields_for_queries_exps_collection.update(fields_that_should_query_using_lowercase_strings)

	return fields_for_queries_exps_collection


# see function above to see how these variables work
fields_for_queries_exps_collection = copy.deepcopy(modifying_exps_collection_queries())
fields_for_queries_seqs_collection = copy.deepcopy(modifying_seqs_collection_queries())


# #############END OF VARIABLES AND FUNCTION USED TO MODIFY QUERIES #########################

# ##################FUNCTIONS USED WHEN SAVING QUERY RESULTS TO FILE ####################
'''

	The following functions help the query class determine how to handle different fields when saving results to file.
	For example, TAB files must always have strings variables saved to file, so these function will assist in determining how to convert fields to strings

'''


def default_fields_to_file(x, splitter=','):
	'''
		Converts all text in x into strings

		Parameters
		----------
		x : dict/mongo document
		splitter : char, the character defines what will be used as the delimiter
	'''
	if x is None:
		return ''
	if type(x) is list:
		return (splitter).join(map('{}'.format, x))
	else:
		return str(x)


def default_fields_to_file_no_commas(x):
	'''
		Converts all text in x into strings, but wants to avoid a comma delimiter (used for CSV files). 
		The delimiter will always be '|'

		Parameters
		----------
		x : dict/mongo document

	'''
	if x is None:
		return ''
	if type(x) is list:
		return ('|').join(map('{}'.format, x))
	else:
		return str(x)


def ReturnStr(x):
	'''
		Convert x to str
	'''
	if x is None:
		return ''
	return str(x)


def ReturnListStr(x, d):
	'''
		Converts a list into a delimited string
		Parameters
		----------
		x : list of string characters
		d : char, the character defines what will be used as the delimiter
	'''
	if x is None:
		return ''
	return (d).join(x)


def ReturnListNotStr(x, d):
	'''
		Converts a list of non-strings into a delimited string
		Parameters
		----------
		x : list of non string values
		d : char, the character defines what will be used as the delimiter
	'''
	# v = str(x[0])
	# for j in range(1,len(x)):
	# v+=d+str(x[j])
	# return v
	if x is None:
		return ''
	return (d).join(map('{}'.format, x))


def define_saving_fields_to_file():
	'''
		This function will define variables we will use to determine how to save non-string fields to a Non-JSON file.
		For example, VREGION.VGENES exists as a list of strings whereas VREGION.VGENE_SCORES exists as a list of numbers.
		When saving as a tab file we would want VREGION.VGENES to appear as a ',' or '|' (if using CSV format) seperated strings
		In this function we will create two variables:

			1) schema_fields_to_file -> A key-value dictionary
				keys: represent fields from the database that we will save to the file
				values: a lambda function defining how to save the values from they field names

			2) schema_fields_to_file_avoid_commas -> this is identical to the variable above (schema_fields_to_file) except it has special
				operations on fields that may have commas. So rather than using a ',' delimiter it uses a | delimiter for lists. This is useful
				when saving results as a CSV

		Usage
		-----
		Assuming the query returns the following document q={'VREGION.VGENES:['A','B,'C']}
		Then the following transformation will convert the query document:

			a = schema_fields_to_file['VREGION.VGENES'](q['VREGION.VGENES'])
			a is now equal to a string 'A','B','C'
	'''
	# schema_fields_to_file => for each key, transforms the value from the database to the defined after getting a query result
	# default functionality : use default_fields_to_file function

	# these fields in the database should already be strings
	should_be_strings = [
		'COMMAND', 'NOTES', 'DATE_UPDATED',
		'PREDICTED_AB_SEQ.NT', 'PREDICTED_AB_SEQ.AA', 'STRAND', 'PREDICTED_CHAIN_TYPE', 'PRODUCTIVE', 'ANALYSIS_NAME', 'RECOMBINATION_TYPE',
		'VREGION.FR1.NT', 'VREGION.FR1.AA', 'VREGION.CDR1.NT', 'VREGION.CDR1.AA', 'VREGION.FR2.NT', 'VREGION.FR2.AA',
		'VREGION.CDR2.NT', 'VREGION.CDR2.AA', 'VREGION.FR3.NT', 'VREGION.FR3.AA',
		'CDR3.NT', 'CDR3.AA',
		'JREGION.FR4.NT', 'JREGION.FR4.AA',
		'GAPPED.VREGION.FR1.NT', 'GAPPED.VREGION.FR1.AA', 'GAPPED.VREGION.CDR1.NT', 'GAPPED.VREGION.CDR1.AA', 'GAPPED.VREGION.FR2.NT', 'GAPPED.VREGION.FR2.AA',
		'GAPPED.VREGION.CDR2.NT', 'GAPPED.VREGION.CDR2.AA', 'GAPPED.VREGION.FR3.NT', 'GAPPED.VREGION.FR3.AA',
		'GAPPED.JREGION.FR4.NT', 'GAPPED.JREGION.FR4.AA',
		'GAPPED.CDR3.NT', 'GAPPED.CDR3.AA'
	]

	# these fields in the database are not in string format (they are usually numbers)
	convert_to_strings = [
		'SETTINGS', 'FILENAME',  # these are exceptions to the format. they are strings in the exps collection, but numbers in this colleciton. If its done right then they will actually be converted to their values in the exps collection by this API when requested
		'SEQ_ID', 'EXP_ID', '_id',
		'VREGION.SHM.NT', 'VREGION.SHM.AA', 'VREGION.SHM.NT_PER', 'VERGION.SHM.AA_PER',
		'JREGION.SHM.NT', 'JREGION.SHM.AA', 'JREGION.SHM.NT_PER', 'JERGION.SHM.AA_PER',
		'VREGION.VGENE_QUERY_START', 'VREGION.VGENE_QUERY_END',
		'JREGION.JGENE_QUERY_START', 'JREGION.JGENE_QUERY_END'
		'CDR3.AA_LENGTH', 'CDR3.NT_LENGTH',
		'VREGION.CDR1.AA_LENGTH', 'VREGION.CDR1.NT_LENGTH', 'VREGION.CDR2.AA_LENGTH', 'VREGION.CDR2.NT_LENGTH',
	]

	# these fields in the database are lists of strings
	convert_liststring_to_strings = [
		'VREGION.VGENES', 'JREGION.JGENES', 'DREGION.DGENES', 'ISOTYPE.GENE'
	]

	convert_nonliststring_to_strings = [
		'VREGION.VGENE_SCORES', 'JREGION.JGENE_SCORES', 'DREGION.DGENE_SCORES',
		'ISOTYPE.MISMATCHES', 'ISOTYPE.SCORES', 'ISOTYPE.PER_ID'
	]

	# Initialize the main variable that will define how to transform different fields from the database
	schema_fields_to_file = defaultdict(lambda: default_fields_to_file)

	# WE WILL BE CHANGING THE FORMAT OF LOCUS!!!
	schema_fields_to_file['LOCUS'] = lambda x: '' if x is None else str(x)

	# If we think this is a waste of computation since they are already strings, then we can remove this line
	schema_fields_to_file.update({field: lambda x: '' if x is None else str(x) for field in should_be_strings})

	# Now set the lambda function for making non-string fields
	schema_fields_to_file.update({field: lambda x: str(x) for field in convert_to_strings})

	# Now set the lambda function for making lists of strings into strings
	schema_fields_to_file.update({field: lambda x: ReturnListStr(x, ',') for field in convert_liststring_to_strings})

	# Now set the lambda function for making lists of numbers into strings
	schema_fields_to_file.update({field: lambda x: ReturnListNotStr(x, ',') for field in convert_nonliststring_to_strings})

	#
	# Now we need to set a special variable that is identical to the schema_fields_to_file, except is to be used when saving as CSV format
	schema_fields_to_file_avoid_commas = defaultdict(lambda: default_fields_to_file_no_commas)

	# Copy the variable schema_fields_to_file into this variable
	for s_keys, s_values in schema_fields_to_file.iteritems():
		schema_fields_to_file_avoid_commas[s_keys] = s_values

	# Now set the lambda function for making list of strings into strings BUT dont use comma delimiter, use '|' instead
	schema_fields_to_file_avoid_commas.update({field: lambda x: ReturnListStr(x, '|') for field in convert_liststring_to_strings})

	# Now set the lambda function for making lists of numbers into strings BUT dont use comma delimiter, use '|' instead
	schema_fields_to_file_avoid_commas.update({field: lambda x: ReturnListNotStr(x, '|') for field in convert_nonliststring_to_strings})

	return schema_fields_to_file, schema_fields_to_file_avoid_commas

# See the function Process_Cursor_For_Output_File to see how these two variables are used
[schema_fields_to_file, schema_fields_to_file_avoid_commas] = define_saving_fields_to_file()


def Get_Schema_Details(possible_metadata_fields, projected_fields, allTrue):
	'''
		This function is used to predict what fields will be returned from the database based on the project query.
		This function assumes the user is using a simple PROJECT command and is not using the aggregation framework methods. See function below for handling aggregation pipelines (Get_Schema_Details_Aggregation_Pipeline).
		This function is useful because when saving results as TAB or CSV format, then we can predict what fields we should see given the schema and query.
		So saving to a TAB or CSV file will not create unnecessary columns unless it has to.

		Parameters
		----------
		possible_metadata_fields : A list of documents from the exps collection. Each document should have an 'anlysis_schema' field if it was annotated.
		projected_fields : A dictionary that corresponds to a mongo-db projection (see http://docs.mongodb.org/manual/reference/method/db.collection.find/)
		allTrue	: boolean varaible defining whether all projected_fields have a '1'

		Outputs
		-------
		possible_fields : A list of fields that should appear from the query and therefore should appear in the file
		analyses_types 	: A list of analyses_types (i.e. IMGT) that will also appear based on the projection
		recombination_types : A list of recombination types that will appear based on the projection
	'''

	analyses_types = []
	recombination_types = []

	# During database insertion, fields that are stored above 'DATA' field (i.e. analysis_name) will not appear in the experiment field 'analysis_schema'
	# But because we know about these fields since these fields are added in by our insert function, then we will explicitly state them here
	overall_schema = {f: 1 for f in fields_above_data_key}

	# Traverse through the list of possible fields (based on the exp ids selected for queries (see function GetIdIntersection))
	# For each experiment, go through the analysis schema, and add to the list of possible fields that can be reported
	for each_metadata in possible_metadata_fields:
		if 'ANALYSIS_SCHEMA' in each_metadata:
			schema_keys = flatten_dictionary(each_metadata['ANALYSIS_SCHEMA']).keys()
			for field in schema_keys:
				f = 'DATA.' + '.'.join(field.split('.')[:-1])
				# Add this field to the overall_schema variable. We expect this field to appear in the output at least once
				overall_schema[f] = 1

		if 'ANALYSES_COUNT' in each_metadata:
			for each_analysis in each_metadata['ANALYSES_COUNT']:
				analyses_types.append(each_analysis)
				# Only choose recombination types whose count > 0
				recombination_types.extend([recomb for recomb, count in each_metadata['ANALYSES_COUNT'][each_analysis].iteritems() if count > 0])

	projected_fields = copy.deepcopy(projected_fields)
	if projected_fields is None:
		# User did not define any fields to project
		projected_fields = {}
	not_integer_projections = {}
	for projections in projected_fields.keys():
		values = projected_fields[projections]
		if isinstance(values, dict):
			projected_fields.pop(projections)
			# Keep track of fields in the projection that are dictionaries instead of integers (again refer to MongoDB find format)
			not_integer_projections[projections] = values

	if projected_fields == {}:
		# If undefined, then mongo will return all possible fields for a document, so we know that all fields that appear in overall_schema should appear
		for i in overall_schema:
			overall_schema[i] = 1
	elif allTrue:
		# All of the fields in the projection are listed to appear, so initialize overall_schema as 0 so that we can convert back fields in projection to 1
		for i in overall_schema:
			overall_schema[i] = 0
	else:
		# All of the fields in the project are listed to not appear, so initialize overall_schema to 1 so that we can convert back fields to 0 based on projections
		for i in overall_schema:
			overall_schema[i] = 1

	for projections, values in projected_fields.iteritems():
		if values == 0:
			# Suppress these fields/these will not appear in the file
			for fields in overall_schema:
				if fields == projections:
					overall_schema[fields] = 0
				elif fields.startswith(projections + '.'):
					overall_schema[fields] = 0
		else:
			# These fields will appear in the file
			for fields in overall_schema:
				if fields == projections:
					overall_schema[fields] = 1
				elif fields.startswith(projections + '.'):
					overall_schema[fields] = 1

	# Any field in the projection that was not a number was probably  a complex projection operator (such as slice operator)
	# So for all of these fields, change their schema values to 1 because they will probably APPEAR in the file
	for projections in not_integer_projections:
		for fields in overall_schema:
			if fields == projections:
				overall_schema[fields] = 1
			elif fields.startswith(projections + '.'):
				overall_schema[fields] = 1

	# When processing results from the database, we alwasy want to remove 'DATA' from the results because it just looks ugly, so remove fields that have the DATA prefix.
	possible_fields = [field[5:] if field.startswith('DATA.') else field for field, value in overall_schema.iteritems() if value == 1]

	return [sorted(list(set(possible_fields))), sorted(list(set(analyses_types))), sorted(list(set(recombination_types)))]


def Get_Schema_Details_Aggregation_Pipeline(possible_metadata_fields, aggregation_pipeline_query):
	'''
		This function is used to predict what fields will be returned from the database based on a mongodb aggregation function
		This function is useful because when saving results as TAB or CSV format, then we can predict what fields we should see given the schema and query.
		So saving to a TAB or CSV file will not create unnecessary columns unless it has to.

		Projections from an aggregation will be slightly more complicated than simple projection from a find command. Therefore we needed a special function to handle aggregation functions.
		For example, fields in an aggregation function can be renamed, or projected, or new fields can be added.
		We need to keep track of all possible fields so that we can predict proper field names that appear in TAB/CSV file.
		It would be annoying to a have a file with a different number of tabs in each line; this function helps prevent that by saying before hand what fields should appear.

		Understanding this function requires a good undestanding of mongodb aggreagation methods.
		Reference: http://docs.mongodb.org/manual/aggregation/

		Parameters
		----------
		possible_metadata_fields : A list of documents from the exps collection. Each document should have an 'anlysis_schema' field if it was annotated.
		projected_fields : A dictionary that corresponds to a mongo-db projection (see http://docs.mongodb.org/manual/reference/method/db.collection.find/)
		allTrue	: boolean varaible defining whether all projected_fields have a '1'

		Outputs
		-------
		possible_fields : A list of fields that should appear from the query and therefore should appear in the file
		analyses_types 	: A list of analyses_types (i.e. IMGT) that will also appear based on the projection
		recombination_types : A list of recombination types that will appear based on the projection
	'''

	def get_new_field_name(new_field_name, original_field, overall_schema):
		'''
			Parses a projection given by new_field_name which contains a  mixture of a field to project and mongod-db operators on how to project the field
			Once successfully parsed and the actual field name being projected is determined, updates overall_schema variable 
		'''
		new_field = new_field_name.split('.$')[0]
		if original_field[0] == '$':
			found_a_value = False
			if original_field.startswith('$$CURRENT.'):
				old_field = original_field[10:]
			# id field has to start with $ mongo operator
			else:
				old_field = original_field[1:]
			for each_key in overall_schema.keys():
				# so we found all values that contain the old key name
				if each_key == old_field:
					found_a_value = True
					overall_schema[new_field] = 1  # new field name
				elif each_key.startswith(old_field + '.'):
					found_a_value = True
					new_key = new_field + '.' + (old_field + '.').join(each_key.split(old_field + '.')[1:])
					overall_schema[new_key] = 1
			if found_a_value is False:
				overall_schema[new_field] = 1
		else:
			overall_schema[new_field] = 1

	possible_fields = fields_above_data_key
	overall_schema = {}
	for i in possible_fields:
		overall_schema[i] = 1
	analyses_types = []
	recombination_types = []

	# Traverse through the list of possible fields (based on the exp ids selected for queries (see function GetIdIntersection))
	# For each experiment, go through the analysis schema, and add to the list of possible fields that can be reported
	for each_metadata in possible_metadata_fields:
		if 'ANALYSIS_SCHEMA' in each_metadata:
			schema = flatten_dictionary(each_metadata['ANALYSIS_SCHEMA']).keys()
			for i in range(len(schema)):
				schema[i] = '.'.join(schema[i].split('.')[:-1])
				overall_schema['DATA.' + schema[i]] = 1

		if 'ANALYSES_COUNT' in each_metadata:
			for each_analysis in each_metadata['ANALYSES_COUNT']:
				analyses_types.append(each_analysis)
				# Only choose recombination types whose count > 0
				recombination_types.extend([recomb for recomb, count in each_metadata['ANALYSES_COUNT'][each_analysis].iteritems() if count > 0])

	# Traverse through the users aggregation pipeline query
	for pipe in aggregation_pipeline_query:
		# Fields are only projected when the pipeline query = $project. so when this occurs, keep track of what fields are being reported.
		if pipe.keys()[0] == '$project':
			# Flatten the projection
			sub_project_query = flatten_dictionary(pipe['$project'])
			# Turn off all key fields by setting value to 0 (by default nothing is projected, unless explicitiley defined)
			for each_key in overall_schema:
				overall_schema[each_key] = 0
			# By defaault, '_id' is always shown
			overall_schema['_id'] = 1
			# Traverse the key/values of the projection. The keys or fields SHOULD either be a fieldname to project or a mongodb project operator
			# We need to extract just the fieldname parts of the projection and ignore mongo operators
			for field, operator in sub_project_query.iteritems():
				if '.$literal.' in field:
					overall_schema[field.replace('.$literal.', '.')] = 1
				elif operator == 0:
					for each_key in overall_schema:
						if each_key == field or each_key.startswith(field + '.'):
							overall_schema[each_key] = 0
				elif operator == 1:
					# Include these fields from overall schema
					for each_key in overall_schema:
						if each_key == field or each_key.startswith(field + '.'):
							overall_schema[each_key] = 1
				elif isinstance(operator, basestring) and operator[0] == '$':
					# Get the field name from this operator
					get_new_field_name(field, operator, overall_schema)	
				else:
					# It is most likely using a mongo operator to project a NEW FIELD, so add the NEW field to overall_schema
					# Only consider the field name up until MONGO operator '$'
					new_field = field.split('.$')[0]
					overall_schema[new_field] = 1
		# Fields can be renamed in the $group stage of the aggregation pipeline, so try to predict any time that happens
		elif pipe.keys()[0] == '$group':
			# Turn off all key fields (by default nothing is projected, unless explicitiley defined)
			for each_key in overall_schema:
				overall_schema[each_key] = 0
			operation_keys = copy.deepcopy(pipe['$group'])
			# Traverse the ID operator. ID determines how to group sequences together
			id_field = operation_keys.pop('_id', None)

			# Analyze the field names under '_id' for grouping
			if not id_field:
				# User said in the group stage that:{_id:null}.
				# Only an id field will be made
				overall_schema['_id'] = 1
			elif isinstance(id_field, basestring):
				# User is only grouping by specific field. this field will be called 'ID'
				# Id must appear in data now
				overall_schema['_id'] = 1
			else:
				# Id is a complex dictionary that groups values by multiple fields
				overall_schema['_id'] = 0
				id_field = flatten_dictionary(id_field)
				for id_keys, id_values in id_field.iteritems():
					if isinstance(id_values, basestring):
						get_new_field_name('_id.' + id_keys, id_values, overall_schema)
					else:
						# No idea, so must rename the fields based on the new field name 
						new_field = id_keys.split('.$')[0]
						overall_schema['_id.' + new_field] = 1

			# Analyze all other field names under group
			operation_keys = flatten_dictionary(operation_keys)
			for keys, values in operation_keys.iteritems():
				if isinstance(values, basestring):
					get_new_field_name(keys, values, overall_schema)
				else:
					# No idea, so must rename the fields based on the new field name
					new_field = keys.split('.$')[0]
					overall_schema[new_field] = 1

	# Remove fields that have the DATA prefix
	possible_fields = [field[5:] if field.startswith('DATA.') else field for field, value in overall_schema.iteritems() if value == 1]

	return [sorted(list(set(possible_fields))), sorted(list(set(analyses_types))), sorted(list(set(recombination_types)))]


'''
	These next two function define how to process documents in a database cursor. We will pass each document from a cursor
	into one of these functions before saving results to a file. The cursor will modify the document based on the variables above
'''


def Simple_Process_Output(document, dump_to_string=False):
	'''
		This is the simplest method that we use to process documents resulting from a database query.
		In this function, the following operations are performed on a query document:
		1) Remove the 'DATA' key from documents. So every field under 'DATA' are now raised one level -> 'DATA.VREGION.VGENES' becomes 'VREGION.VGENES'
		2) Convert an object id from a document into a string
		3) Concatenate the SeqId and ExpId field into a single string

		Parameters
		----------
		document : dict
			A query document
		dump_to_string : boolean
			If true, then a string is returned

	'''
	RemoveObjId(document)
	data_info = document.pop('DATA', None)
	if data_info:
		document.update(data_info)
	# Concatenate the exp id to the seq id
	if idIdentifier in document and expIdentifier in document:
		document[idIdentifier] = str(document[idIdentifier]) + '::' + str(document[expIdentifier])
	if dump_to_string is False:
		return document
	# Dump data
	try:
		return json.dumps(document)
	except:
		# Try to use bson dump
		return bson_dumps(document)


def Process_Cursor_For_Output_File(document, schema_output_fnc):
	'''
		This function is another method used for processing mongo documents resuling from a database query.
		It is meant to be used when saving a document as a flattened dictionary in a specific file format only: Either FASTQ, FASTA, TAB, CSV, or IGREP files.
		It is not used when saving results as a JSON format.

		In this function, the following operation are performed on a query document:

		1) Remove the 'DATA' key from documents. So every field under 'DATA' are now raised one level -> 'DATA.VREGION.VGENES' becomes 'VREGION.VGENES'

		2) Flatten all documents so that they are no longer nested; while flattenign documents, convert any objectids to strings

		3) Use the variable 'schema_output_fnc' to define exactly how to convert specific fields in the document to string format

		4) Concatenate the SeqId and ExpId field into a single string

		Parameters
		----------
		document : A document or dictionary from a mongo db query
		schema_output_fnc : A dictionary whose keys refer to a specific database field and whose value is a lamda function defining how to convert that field to a string
		See the global variables schema_fields_to_file and schema_fields_to_file_avoid_commas
	'''

	if schema_output_fnc is None:
		schema_output_fnc = schema_fields_to_file

	data_info = document.pop('DATA', None)
	if data_info:
		document.update(data_info)
	# 1 million lines => takes 30 seconds to perform step using CYTHON, takes 1 minute to perform step using python code
	document = flatten_dictionary(document)

	# Force values from docs to be strings
	for field, value in document.iteritems():
		if not isinstance(value, basestring):
			# Its not a string, we need to modify value to string using schema_fields_to_file
			document[field] = schema_output_fnc[field](value)

	# Concatenate the exp id to the seq id
	if idIdentifier in document and expIdentifier in document:
		document[idIdentifier] = str(document[idIdentifier]) + '::' + str(document[expIdentifier])

	return document
# ### END OF FUNCTIONS USED FOR HANDLING METHODS FOR SAVING RESULTS TO A FILE ######

# ####FUNCTIONS FOR FILTERING QUERIES BASED ON USER PERMISSIONS #####
'''
	This following section are function(s) we use to filter out results from users based on their permissions.
'''


def GetIdIntersection(list_of_allowed_ids, id_query=None):
	'''
		This function will take the intersection of all possible experiments the user wants to search AND has access to
		So it will make sure that the user only runs a query on the seqs collection using exp_ids defined by the intersection of all exp_ids requested

		Parameters
		----------
		list_of_allowed_ids : It corresponds to all experiments the user has access to
								The following values are allowed for list_of_allowed_ids: (a) a single object ID, (b) a list of ObjectIds, (c) a list of lists of ObjectIds
		id_query : A mongo query on '_id' (exps collection) or 'EXP_ID' (seqs collection)
					The following values are allowed for id_query: (a) single value that can be converted to an ObjectId, (b) dictionary containing the following keys only: '$eq','$in','$ne', '$nin',

		Summary
		-------
			a) First, it takes the intersection of all values in list_of_allowed_ids
			b) Then once this id list has been filtered down, it will analyze the input 'id_query' using the rules discussed above
	'''

	# single element passed in for list_of_allowed_ids
	if not isinstance(list_of_allowed_ids, list):
		# convert to lists of lists...
		list_of_allowed_ids = [[convert_to_objectid(list_of_allowed_ids)]]
	# list passed in for list_of_allowed_ids
	else:
		# make sure if the user wants a LIST OF LISTS, then every elment in TOP list is LIST; DO NOT allow a 'MIXTURE' of some lists, some not lists
		are_elements_lists = [isinstance(input, list) for input in list_of_allowed_ids]

		if sum(are_elements_lists) == 0:
			# no element was a list, so convert it into a list of lists
			list_of_allowed_ids = [list_of_allowed_ids]
		elif sum(are_elements_lists) != len(list_of_allowed_ids):
			# at least one elmeent was A LIST, but NOT ALL elements
			raise Exception("The list variable 'list_of_allowed_ids' is improperly formatted. You may only pass in 1) a single value, 2) a list of possible values, or 3) a list of LISTS of possitble values. You passed in the following ", str(list_of_allowed_ids))

	# List_of_allowed_ids has been properly formatted as list of lists of object_ids

	# First convert all values into a unique set of ObjectId format
	combined_ids = set([convert_to_objectid(id) for id in list_of_allowed_ids[0]])

	# Now obtain intersection of all ids defined in list_of_allowed_ids
	for i in range(1, len(list_of_allowed_ids)):
		combined_ids = combined_ids & set([convert_to_objectid(id) for id in list_of_allowed_ids[i]])

	# Now if the user also passed in an id_query, calculate the intersection of that value
	if id_query:
		if isinstance(id_query, list):
			combined_ids = combined_ids & set([convert_to_objectid(exp) for exp in id_query])
		elif isinstance(id_query, dict):
			# First we need to ensure that the id_query has been properly loaded into ObjectId's using bson_loads.
			# To ensure this, we will first use dumps on the dictionary, and then reload using loads (So if an object id is {'oid:""} , this function will make it into a proper object => ObjectId(""))
			id_query = bson_loads(bson_dumps(id_query))
			operators_supported = ['$in', '$eq', '$ne', '$nin']  # We only alow dictionaries with these keys
			for each_key, query_val in id_query.iteritems():
				# Perform set intersection
				if each_key == '$in':
					if not isinstance(query_val, list):
						query_val = [query_val]
					combined_ids = combined_ids & set([convert_to_objectid(exp) for exp in query_val])
				elif each_key == '$eq':
					# For '_id' and 'EXP_ID' the value can never be EQUAL to a list
					if isinstance(query_val, list):
						raise Exception("The ObjectId can never be equal ('$eq') to a list. If you want to allow a range of posible object id's, use '$in' operator\nPassed in query: {0}".format(bson_dumps(id_query)))
					combined_ids = combined_ids & set([convert_to_objectid(query_val)])
				elif each_key == '$ne':
					# For '_id' and 'EXP_ID' the value can never be EQUAL to a list
					if isinstance(query_val, list):
						raise Exception("The ObjectId can never be not-equal ('$ne') to a list. If you want to allow a range of posible object id's, use '$in' operator\nPassed in query: {0}".format(bson_dumps(id_query)))
					# Perform  the set difference between object ids (because we are saying != or $ne)
					combined_ids = combined_ids - set([convert_to_objectid(query_val)])
				elif each_key == '$nin':
					if not isinstance(query_val, list):
						query_val = [query_val]
					combined_ids = combined_ids - set([convert_to_objectid(exp) for exp in query_val])
				else:
					raise Exception("You have requested the following query on experiment id: {0}. Currently we only support the following mongo operators on the exp id: {1}".format(bson_dumps(id_query), ','.join(operators_supported)))
		else:
			combined_ids = combined_ids & set([convert_to_objectid(id_query)])

	combined_ids = list(combined_ids)

	# Return a proprely formatted mongo query: use $in if we have a range of ids to search, or just use a single value if we have only one object id to search
	return combined_ids[0] if len(combined_ids) == 1 else {'$in': combined_ids}

# #####END OF FUNCTIONS FOR FILTER RESULTS BASED ON USER PERMISSIONS ####


# Creating an instance of the class
class RunQuery(object):
	'''
		A class for interacting with the BIGG Lab MongoDB database for NGS sequencing immunological repertoires.
		To run this query you need to define the path to the database with proper admin credentials,
		Optionally, you can define the name of the user requesting access to the database.
		If the name of the user is not provided, then it will automatically treat this person as an admin since the person accessing this class already knows admin privelages to the database.
		The name of the user requesting access to the database can be used to filter results from the database using the proxy query-class (See immunogrep_query_dbproxy_api.py)

		.. note::

			We assume that users directly using this class knows the admin username and password required to access the exps and seqs collections of the database.
		If this is not true, then users need to access the database via the proxy. Refer to module immunogrep_query_dbproxy_api.py

		Parameters
		----------
		db_connect_data : three element tuple. Mandatory.
			The tuple defines how to connect to the mongodb instance: (username with access to seqs and exps collections, username password, path to database)
		dbuser : string. Optional. default = None
			This string should refer to a specific database user. If this string is defined, then all requested queries will be filtered based on experiments this user has access to.
			This is how the proxy (see immunogrep_query_dbproxy_api) handles requests from different users.
			.. important::

				If it is equal to None then it will assume that the user is an administrator
			If it is equal to an empty string then it will assume that the user is a public user (basically has the minimium access to the database)
		modify_query_values_to_follow_db_schema : boolean. Optional. default = True
			When set to true, this will ensure that any queries on fields we define in our default schema are in the proper format. For example, if a user
			queries for 'VREGION.VGENE_SCORES': {$gt:'3'} where '3' is a string, when set to true, it will modify the query to ensure that it is a number => 'VREGION.VGENE_SCORES': {$gt:3}
		redirect_fields_for_improved_queries : boolean. Optional. default = True
			When set to true, this will modify some queries to use fields optimized to return better results. For example, if a user wants to query for experiments whose name is
			'Experiment 1 test', then they would write 'EXPERIMENT_NAME':'Experment 1 test'. However, if it actually exists in the database as Experiment_1_test or Experiment 1 Test, then their
			query will not work. When creating experiments, we add extra fields for common fields that may have idiosyncracies in spacing and case sensitivity; these duplicated fields are stored in a new field called 'DUPLICATED_FIELDS'
			Therefore, when set to true, the query above will become 'DUPLICATED_FIELDS.EXPERIMENT_NAME':'Experiment 1 test'. In this situation, the query would work.

		Attributes
		----------
		query_results : int, list, dict, db cursor, generator
			This attribute stores the result from the latest query run. Iterating through the class, saving results to file, or running get_results,
			will return values stored in this attribute.
			The format of query_results will vary depending on the type of query run by the user.
		db_path : mongodb path to collections
			Stores the current mongodb connection to the database. Access to the exps, seqs, and users collections can be provided from this attribute
			i.e. db_path.seqs, db_path.exps, etc..
		user : dict defining current user querying database
			Using the value from dbuser, this is a dict containing all information about that db_user. keys for this attribute are: 'user', 'administrator', 'write_access', 'curator', 'email'
		allowed_experiment_ids : list of strings
			A list containing the _id of all experiments the user has access to
		exp_metadata : dict
			A dictionary storing all experiment metadata the user has read access to. Key of dictionary = the _id field of each experiment in the exps collection
		delim_file_headers : list of strings
			This will store all fields from the most current query that can be reported as columns to a delimited file. We use this attribute when saving results to
			a delimited file
		default_filename : string
			A string representing a default filename to use for saving query results to a file when not defined by user
		to_file : boolean
			A boolean keeping track of whether a query was saved to a file or not (currently not used)
		files_created : Write_Files object
			This creates a default_dict for keeping track of all files currently created when saving results from queries

		Examples
		--------
		Query for experiments whose experiment name is Demo 1, return the project name for each experiment

		>>> myquery = RunQuery(('biotseq.icmb.utexas.edu','someadminuser','somepassword'),'cchrysostomou')
		>>> exp_results = [exp for exp in myquery.query_exps_collection(exps_query={'EXPERIMENT_NAME':'Demo 1'},project_fields={'PROJECT_NAME':1})]

		Query for sequences from experiments whose experiment name is Demo 1, annotated by IMGT, and whose VGENE is IGHV3-2; Convert the results into FASTQ format, and save the results to a file
		>>> myquery.query_seqs_collection(exps_query={'EXPERIMENT_NAME':'Demo 1'},analysis_name = ['IMGT'], seqs_query={'DATA.VREGION.VGENES':'IGHV3-2'}).save_as_fastq(prefix_file_path='tesetquery')

		Get the amino acid CDR3 length distribution for the experiment, Demo 1, of sequences annotated by IMGT and containing the VGENE IGHV3-2
		>>> myquery.cdr3_length_distribution(feature='AA', analysis_name='IMGT', metadata_query = {'EXPERIMENT_NAME':'Demo 1'}, filter_by_query={'DATA.VREGION.VGENES':'IGHV3-2'}).save_as_json(prefix_file_path='cdr3dist')

		See also
		--------
		immunogrep_query_dbproxy_api

	'''

	def __init__(self, db_connect_data=None, dbuser=None, modify_query_values_to_follow_db_schema=True, redirect_fields_for_improved_queries=True, to_file=True, file_prefix=None):
		# Store values user passed into as attributes
		self.modify_query_values_to_follow_db_schema = modify_query_values_to_follow_db_schema
		self.redirect_fields_for_improved_queries = redirect_fields_for_improved_queries

		# All query results will be stored in this db cursor
		self.query_results = None
		# This variable will be used to figure out what fields are present in delimited files (CSV/TAB).  See the function(s) Get_Schema_* above. They will fill in this variable
		self.delim_file_headers = []
		self.process_cursor_fnc = None

		self.debug_file = 'query_debug_log.txt'
		# A dictionary of metadata that the user has access to
		self.exp_metadata = {}
		# A list of experiment ids, the user has access to
		self.allowed_experiment_ids = []

		# Connect to database using path, username, password
		self.db_path = connectToIgDatabase(db_connect_data[0], db_connect_data[1], db_connect_data[2])
		self.default_filename = 'IGREP_Query_' + re.sub('[\:_\- ]', '', str(datetime.now()))
		self.files_created = Write_Files()  # This is basically a default dict for writing results to files (see the class definition below)
		self.to_file = False
		# Get user permissions
		if dbuser is None:
			# if it was not explicitly said, then assume its an administrator
			self.user = {'user': '', 'administrator': True}
			dbuser = 'admin'
		else:
			self.user = self.db_path.users.find_one({'user': dbuser})

		if not self.user:
			self.user = defaultdict(str)
			self.no_user_found = True
		else:
			self.user = defaultdict(str, self.user)
			self.no_user_found = False

		# Create a log file
		if not os.path.isdir(dbuser):
			os.mkdir(dbuser)
		self.logged_output = open(os.path.join(dbuser, 'logfile.txt'), 'a')
		self.logged_output.write('{1}: User {0} has accessed database\n'.format(dbuser, str(datetime.now())))

		# Get all metadata user has access to
		self.exp_metadata = self._query_exp_docs_with_read_access()
		# Get a list of experiments user has access to
		self.allowed_experiment_ids = [ObjectId(id) for id in self.exp_metadata]

		self.name = None

	def __iter__(self):
		'''
			Creates an iterator for the class. When run, it will iterate throug the query results
			attribute, self.query_results

			yields
			------
			See function iterresults() below
		'''
		return self.iterresults()

	def get_results(self, chunk_size=None, flatten_documents=False):
		'''
			Method for return results from a query. This will return results, from the query results
			attribute self.query_results, to memory and not to a file.
			Results can be returned in chunks, or once at a time

			.. important::
				If you would like to save results to a file please see the 'save_as_*' functions below

			.. note::
				The DATA key will always be removed from the query results (DATA.VREGION becomes just VREGION)

			.. note::
				The database cursor from a query will be exhausted when running this function because we do not use iter tools tee

			Parameters
			----------
			chunk_size : int, default=None
				This defines how to return results back
				When chunk_size is None, then results are returned based on their format in self.query_results
					i.e. If a dict, returns a dict; if a cursor returns a cursor, if a list returns list
				When chunk_size is an integer, then results are loaded into chunks inside a list.
					When chunk_size is 1, each iterator of generator will return [result1]
					When chunk_size is 5, each iterator of generator will return [result1,result2,result3,result4,result5]
					When chunk_size is 0, then all results will be returned one-by-one but NOT inside of a list
					.. warning::
						a chunk_size of 0 will NOT yield results in a list
			flatten_documents : boolean, default=False
				When True, all documents from a query will be flattened so that they are not nested:
				i.e. 'VREGION':{CDR1:{AA:'ACT'}} becomes 'VREGION.CDR1.AA': 'ACT'

			Returns
			-------
			int, long, float, complex
				Returns a single number when self.query_results is a number (i.e. when using count function) and chunk_size is None
			dict
				Returns a dictionary (flattened if flatten_documents is true) when self.query_results is a dict and chunk_size is None
			generator
				1) Returns a generator when self.query_results is a db cursor or a generator and chunk_size is none
				2) Returns a generator for iterating through self.query_results when chunk_size is 0 
				3) Returns a generator for any format of self.query_results when chunk_size > 0.

					.. note:: List of results
					In this case (chunk_size > 0), all results from this generator will be yielded inside a list

			Examples
			--------
			>>> import immunogrep_db_query_api as query
			>>> myquery = query.RunQuery((user,password,'biotseq.icmb.utexas.edu'))
			>>> for results in myquery.get_user_access_info().get_results():
			... 	print(results)
		'''
		if chunk_size is None:
			# Return results precisely as we see them
			if isinstance(self.query_results, (list, basestring, int, long, float, complex)):
				return self.query_results
			elif isinstance(self.query_results, dict):
				if flatten_documents:
					return flatten_dictionary(Simple_Process_Output(self.query_results))
				else:
					return Simple_Process_Output(self.query_results)
			else:
				return self._yield_results(flatten_documents)
		elif chunk_size == 0:
			self.iterresults(flatten_documents)
		else:
			return self._yield_chunks(chunk_size, flatten_documents)

	def iterresults(self, flatten_documents=False):
		'''
			Iterator function for parsing results from the query results attribute: self.query_results
			This function returns a generator that can be iterated in a for loop

			.. note::
				The DATA key will always be removed from the query results (DATA.VREGION becomes just VREGION)

			.. note::
				The database cursor from a query will be exhausted when running this function because we do not use iter tools tee

			Parameters
			----------
			flatten_documents : boolean, default=False
				When True, all documents from a query will be flattened so that they are not nested:
				i.e. 'VREGION':{CDR1:{AA:'ACT'}} becomes 'VREGION.CDR1.AA': 'ACT'

			Yields
			------
			dict
				1) If self.query_results is equal to a dict (i.e. from a find_one command)
				2) self.query_results is equal to a db cursor or generator (i.e. from query_seqs_collection or query_exps_collection)
			int, float, complex, long
				If self.query_results is equal to a number (i.e. from a count command)
			list
				If self.query_results is equal to a list (i.e. output from allowed_experiment_ids)
		'''
		if flatten_documents:
			if isinstance(self.query_results, (list, basestring, int, long, float, complex)):
				yield self.query_results
			elif isinstance(self.query_results, dict):
				yield flatten_dictionary(Simple_Process_Output(self.query_results))
			else:
				for doc in self.query_results:
					yield flatten_dictionary(Simple_Process_Output(doc))
		else:
			if isinstance(self.query_results, (list, basestring, int, long, float, complex)):
				yield self.query_results
			elif isinstance(self.query_results, dict):
				yield Simple_Process_Output(self.query_results)
			else:
				for doc in self.query_results:
					yield Simple_Process_Output(doc)

	def get_created_files(self):
		'''
			Returns the filenames of all queries saved to a file

			Returns
			-------
			list of dict
				Each element in the list represents a file that was created
				[{'filename': path_of_files_created, 'doc_count': documents per file}]

		'''
		return self.files_created.get_files()

	def _yield_chunks(self, chunk_size, flatten_documents):
		if isinstance(self.query_results, (list, basestring, int, long, float, complex)):
			yield [self.query_results]
		elif isinstance(self.query_results, dict):
			if flatten_documents:
				yield [flatten_dictionary(Simple_Process_Output(self.query_results))]
			else:
				yield [Simple_Process_Output(self.query_results)]
		else:
			current_chunk = 0
			return_chunks = []
			for doc in self.query_results:
				doc = Simple_Process_Output(doc)
				if flatten_documents:
					doc = flatten_dictionary(doc)
				return_chunks.append(doc)
				current_chunk += 1
				if current_chunk >= chunk_size:
					yield return_chunks
					return_chunks = []
					current_chunk = 0
			if return_chunks:
				yield return_chunks

	def _yield_results(self, flatten_documents):
		if flatten_documents:
			for doc in self.query_results:
				yield flatten_dictionary(Simple_Process_Output(self.query_results))
		else:
			for doc in self.query_results:
				yield Simple_Process_Output(doc)

	# #############QUERY CURSOR FUNCTIONS##########################

	'''
		The following methods add in some mongo cursor functions to a set of query commands
		i.e. You can do query_seqs_collection().count() or .limit()
	'''

	def count(self, with_limit_and_skip=False):
		'''
			Count the number of results in a query

			Parameters
			----------
			with_limit_and_skip : boolean, default = False
				When performing the count, consider limit and skip settings.
				For example, if a limit of 100 was set on the query, and you perform a count where with_limit_and_skip is True,
				then this function will return 100; If False it would ignore any limit or skip settings set and count the entire query.

			Returns
			-------
			self
		'''
		try:
			self.query_results = self.query_results.count(with_limit_and_skip)
		except Exception as e:
			raise Exception("The count function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	def limit(self, limit):
		'''
			Sets a limit for the number of documents to return.

			Parameters
			----------
			limit : integer
				.. note::
					A limit of 0 is the same as not setting a limit

			Returns
			-------
			self
		'''
		try:
			self.query_results = self.query_results.limit(limit)
		except Exception as e:
			raise Exception("Limit function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	def skip(self, num_skip):
		'''
			Skips the first set of documents in a query

			Parameters
			----------
			num_skip : integer

			Returns
			-------
			self
		'''
		try:
			self.query_results = self.query_results.skip(num_skip)
		except Exception as e:
			raise Exception("Skip function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	# set up a sort command for general queries
	def sort(self, sort_tuples):
		'''
			Sorts documents in a query by a set of tuples

			Parameters
			----------
			sort_tuples : list of tuples
				A list of (key, direction) pairs specifying the sort order for this query.
				Direction is 1 for ascending or -1 for descending (or pymongo.ascending/descending)

			Returns
			-------
			self
		'''
		try:
			self.query_results = self.query_results.sort(sort_tuples)
		except Exception as e:
			raise Exception("Sort function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	# For debugging mongo queries debugging
	def explain(self):
		'''
			Gets a dict explaining the query performance

			.. note::
			Refer to mongodb explain function

			Returns
			-------
			self
		'''
		try:
			self.query_results = self.query_results.explain()
		except Exception as e:
			raise Exception("Explain function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	def distinct(self, field_name):
		'''
			Adds a distinct function to the query. This will collapse the results into a list of distinct values based on a field_name

			Parameters
			----------
			field_name : a string defining what field we want to find distinct values for

			Returns
			-------
			self

		'''
		try:
			self.query_results = self.query_results.distinct(field_name)
		except Exception as e:
			raise Exception("Distinct function can only be used on a mongodb cursor! Error message: ", str(e))
		return self

	# QUERIES ON USERS COLLECTION# ########################################################################
	'''
		The following methods run queries on the users collection in the database
	'''

	def get_db_user_list_info(self):
		'''
			Queries users collection for a list of documents containing information of all users currently registered with the database

			Datatype returned to self.query_results : list of dicts
				This function returns a list of dictionaries to the query results attribute, self.query_results
				At a minimum, each dictionary will contain the following keys: _id, user, name, lab

			Returns
			-------
			self
		'''
		db = self.db_path.users

		if self.user['administrator']:
			self.query_results = [u for u in db.find()]
		else:
			# If not administrator then just report username, full name, and lab
			self.query_results = [u for u in db.find({}, {'_id': 0, 'user': 1, 'name': 1, 'lab': 1})]
			self.delim_file_headers = []
		self.delim_file_headers = ['username', 'name', 'administrator', 'curator', 'write_access', 'lab', 'email']
		self.list_to_gen()
		return self

	def get_user_access_info(self):
		'''
			Queries users collection for a dictionary containing information about the current user
			i.e. returns the following in key-value of user: {username:, name:, administrator:, curator:, write_access:, lab:, email:, }

			Datatype returned to self.query_results : dict
				This function returns a single dictionary to the query results attribute, self.query_results
				The dictionary will contain the following keys: username, name, administrator, curator, write_access, lab, email

			Returns
			-------
			self
		'''
		self.query_results = copy.deepcopy(self.user)
		if self.no_user_found:
			self.user['user'] = ''
		self.delim_file_headers = ['username', 'name', 'administrator', 'curator', 'write_access', 'lab', 'email']
		return self

	# QUERIES ON metadata/exps COLLECTION# ###############################################################
	'''
		The following methods correspond to queries performed on the exps collection
	'''

	def get_accessible_exp_ids(self):
		'''
			Queries for a list of all experiment ids where the current user has read access

			.. note::
				ObjectIds are converted into strings

			Datatype returned to self.query_results : list of strings
				This function returns a list of object ids converted to strings to the query results attribute, self.query_results

			Returns
			-------
			self
		'''
		self.query_results = [str(id) for id in self.allowed_experiment_ids]
		self.delim_file_headers = ['EXPERIMENT_IDS']
		self.list_to_gen()
		return self

	def query_exp_ids_owned_by_current_user(self):
		'''
			Return a list of all experiment ids where the user is explicity defined as an 'owner of the experiment'

			.. note::
				ObjectIds are converted into strings

			Returns
			-------
			self

		'''
		exp_data = copy.deepcopy(self.exp_metadata)
		write_access_exps = [str(id) for id in exp_data if self.user['user'].lower() in exp_data[id]['OWNERS_OF_EXPERIMENT']]
		self.query_results = copy.deepcopy(write_access_exps)
		self.delim_file_headers = ['EXPERIMENT_IDS']
		self.list_to_gen()
		return self

	def get_write_accessible_exp_ids(self):
		'''
			Return a list of all experiment ids the user has write access to

			.. note::
				ObjectIds are converted into strings

			Returns
			-------
			self
		'''
		if self.user['administrator']:
			# Return all allowed experiments
			self.query_results = [str(id) for id in self.allowed_experiment_ids]
		else:
			# Only return experiments with write access
			self.query_exp_ids_owned_by_current_user()
		self.delim_file_headers = ['EXPERIMENT_IDS']
		self.list_to_gen()
		return self

	def query_exp_docs(self, object_id_list=[], write_access_only=False):
		'''
			Queries for all metadata from exps collection that the user has read-access to

			Parameters
			----------
			object_id_list : list of strings
				Filters results such that only experiments whose '_id' is listed in this list is considered
			write_access_only : boolean
				When true, filters results by only considering experiments where user is listed as one of the 'OWNERS_OF_EXPERIMENT'

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results

			Returns
			-------
			self
		'''
		if object_id_list:
			object_id_list = list(set(self.allowed_experiment_ids) & set([convert_to_objectid(o) for o in object_id_list]))
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
		self.delim_file_headers = list(set(self.delim_file_headers))
		self.list_to_gen()
		return self

	def query_exp_docs_with_ownership(self, object_id_list=[]):
		'''
			Queries for all metadata from exps collection that the user has write access to/is listed as one of the 'OWNERS_OF_EXPERIMENT'

			.. note::
			This function is identical to the preceeding function query_exp_docs where write_access_only is True

			Parameters
			----------
			object_id_list : list of strings
				Filters results such that only experiments whose '_id' is listed in this list is considered

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results

			Returns
			-------
			self
		'''
		if object_id_list:
			object_id_list = list(set(self.allowed_experiment_ids) & set([convert_to_objectid(o) for o in object_id_list]))
		else:
			object_id_list = self.allowed_experiment_ids
		self.query_results = [self.exp_metadata[ids] for ids in self.exp_metadata if ObjectId(ids) in object_id_list and self.user['user'].lower() in self.exp_metadata[ids]['OWNERS_OF_EXPERIMENT']]
		self.delim_file_headers = []
		for each_doc in self.query_results:
			self.delim_file_headers.extend(each_doc.keys())
		self.delim_file_headers = list(set(self.delim_file_headers))
		self.list_to_gen()
		return self

	def query_exp_docs_containing_sequences(self, write_access_only=False, analyses=[]):
		'''
			Queries for all metadata from all experiments in exps collection that contain sequences in the seqs collection
			(Queries for experiments where SEQ_COUNT is > 0)

			Parameters
			----------
			write_access_only : boolean
				When true, filters results by only considering experiments where user is listed as one of the 'OWNERS_OF_EXPERIMENT'
			analyses : list of strings
				When not empty, then this parameter will filter results by only returning documents that contain annotated information listed in the list
				i.e. when analyses = ['IMGT','IBLAST'], will only return experiments that contain annotation information from either IMGT or IGBLAST

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results

			Returns
			-------
			self
		'''
		db = self.db_path
		object_id_list = self.allowed_experiment_ids
		query = {'_id': {'$in': object_id_list}, 'SEQ_COUNT': {'$gt': 0}}
		if write_access_only and self.user['administrator'] is False:
			query['OWNERS_OF_EXPERIMENT'] = self.user['user'].lower()
		if analyses and analyses is not list:
			analyses = [analyses]
		for a in analyses:
			query['ANALYSES_COUNT.{0}'.format(a.upper())] = {'$gt': 0}
		self.query_results = [exp for exp in db.exps.find(query, {'DUPLICATED_FIELDS': 0})]

		self.delim_file_headers = []
		for each_doc in self.query_results:
			self.delim_file_headers.extend(each_doc.keys())
		self.delim_file_headers = list(set(self.delim_file_headers))
		self.list_to_gen()
		return self

	def query_exps_collection(self, exp_id=None, exps_query={}, project_fields={'DUPLICATED_FIELDS': 0}, include_duplicated_fields=False):
		'''
			Runs a general query on the exps collection using a standard mongo query .find() function

			Parameter
			----------
			exp_id : list of strings, defaut = None
				Filters results such that only experiments whose '_id' is listed in this list is considered
			exps_query : dict, default = {}
				A mongo query using mongo operators
			project_fields : dict, default = {'DUPLICATED_FIELDS': 0}
				A mongo projection query. Fields to project are represented as keys, and values are 0 or 1 depending on whether to suppress or include field
			include_duplicated_fields : boolean, default = False
				If true, then will include the added in field DUPLICATED_FIELDS to the projection query

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results

			Returns
			-------
			self
		'''
		if exp_id is None:
			# Default var which means undefined, so pass in allowed_experiment_ids
			# Basically assume that user does not knwo what to search,so wants to search all avaialabe experiments
			exp_id = self.allowed_experiment_ids
		if not isinstance(exp_id, list):
			exp_id = [exp_id]
		db = self.db_path
		# Modify the query to improve success of a match
		# 1) Redirect select fieldname to fields created by server to imporve queries
		# 2) Ensure that the format of the values provided to the field matches the format of the field put in database
		q = Parse_Mongo_Query_Expression(exps_query, dict_defining_value_transformation=fields_for_queries_exps_collection, redirection_fields=redirect_exp_collection_fields, modify_query_values_to_follow_db_schema=self.modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries=self.redirect_fields_for_improved_queries)
		# We will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)
		q['_id'] = GetIdIntersection([self.allowed_experiment_ids, exp_id], q.pop('_id', None))

		self.logged_output.write('Query was modified to: {0}\n'.format(str(q)))
		project = copy.deepcopy(project_fields)

		if project:
			allTrue = sum(project.values())
			if allTrue == len(project):  # The user projected fields using '1'. they listed select fields to project
				project['_id'] = 1  # Also always project '_id' field
				if include_duplicated_fields:
					project['DUPLICATED_FIELDS'] = 1
			elif allTrue == 0:  # The user suppressed feilds using '0'
				if not(include_duplicated_fields):
					project['DUPLICATED_FIELDS'] = 0
				else:
					# Just incase user explicityl said to include_duplicated_fields
					project.pop('DUPLICATED_FIELDS', None)
		else:
			project = None
		self.logged_output.write('The following values will be projected: {0}\n'.format(str(project)))
		if type(q) is not dict:
			raise Exception("ERROR: query_command must be a dictionary for query experiment collection")

		query = {f: v for f, v in q.iteritems()}
		self.query_results = db.exps.find(query, project)
		return self

	def query_distinct_exp_collection_values(self, unique_fields=[]):
		'''
			Returns distinct values for select fields from CURATED experiments in the exps collection

			.. note::Filtered results
			If the current user is not an administrator, then it only returns distinct values from fields found within experiments user has read access to.

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results


			Parameters
			----------
			unique_fields : dict
				The keys of the dict corresponds to fields defined by unique_fields, and the values of each key is a list of unique values for that field in curated experiments

			Returns
			-------
			self
		'''
		db = self.db_path
		allow_unique_metadata = [
			'PAIRING_TECHNIQUE', 'CELL_MARKERS_USED', 'LIST_OF_POLYMERASES_USED', 'POST_SEQUENCING_PROCESSING: PHI_X_FILTER',
			'POST_SEQUENCING_PROCESSING: QUALITY_FILTER', 'POST_SEQUENCING_PROCESSING: PROCESS_R1_R2_FILE',
			'SPECIES', 'SEQUENCING_PLATFORM', 'CHAIN_TYPES_SEQUENCED', 'CELL_TYPES_SEQUENCED', 'ISOTYPES_SEQUENCED',
			'TEMPLATE_TYPE', 'REVERSE_PRIMER_USED_IN_RT_STEP'
		]

		# Get a list of curated experiments
		all_experiments_metadata = db.exps.find({'CURATED': True})
		# Determine which curated experiments current user has access to
		allowed_experiments_metadata = db.exps.find({'_id': {'$in': self.allowed_experiment_ids}, 'CURATED': True})
		if unique_fields:
			if type(unique_fields) is not list:
				unique_fields = [unique_fields]
			if self.user['administrator']:
				allow_unique_metadata = unique_fields
			self.query_results = {field: all_experiments_metadata.distinct(field) if field in allow_unique_metadata else allowed_experiments_metadata.distinct(field) for field in unique_fields}
		else:
			self.query_results = None
		self.delim_file_headers = unique_fields
		return self

	# PRIVATE FUNCTION...should not be used by user
	def _query_exp_docs_with_read_access(self):
		'''
			Get all experiments user has read access to
		'''
		self.query_results = None
		db = self.db_path
		if self.user['administrator']:
			exp_docs = {str(exp['_id']): exp for exp in db.exps.find({})}
		else:
			read_access_names = ['all'] + [self.user['user'].lower()] + ['lab_' + lab.lower() for lab in self.user['lab'] if lab]
			exp_docs = {str(exp['_id']): exp for exp in db.exps.find({'READ_ACCESS': {'$in': read_access_names}})}
		return exp_docs

	# PRIVATE FUNCTION...should not be used by user
	def _query_exp_docs_with_write_access(self):
		'''
			Get all experiments user has write access to
		'''
		self.query_results = None
		db = self.db_path
		if self.user['administrator']:
			exp_docs = {str(exp['_id']): exp for exp in db.exps.find({})}
		else:
			write_access_name = self.user['user'].lower()
			exp_docs = {str(exp['_id']): exp for exp in db.exps.find({'OWNERS_OF_EXPERIMENT': write_access_name})}
		return exp_docs

	# QUERIES ON seqs COLLECTION# ###############################################################
	'''
		The following methods correspond to queries performed on the seqs collection. It will query for NGS and annotated data.
	'''

	def query_seqs_collection(self, exp_id=None, analysis_name=[], recombination_type=[], seqs_query={}, project_fields={'QUERY_DATA': 0}, exps_query={}, include_original_ngs_seq=False, limit=0):
		'''
			Runs a general query on the seqs collection using a standard mongo query .find() function

<<<<<<< HEAD
			Parameter
=======
			Parameters
>>>>>>> 670fa3e5030646debf2cfc5f20d277782898924f
			----------
			exp_id : list of strings, defaut = None
				Filters results such that only seqs results whose 'EXP_ID' is listed in this list is considered
			analysis_name : list of strings, default = []
				When not empty, filters results to only include documents whose ANALYSIS_NAME is within this list
			recombination_type : list of strings, default = []
				When not empty, filters results to only include documents whose RECOMBINATION_TYPE is within this list
			seqs_query : dict, default = {}
				A mongo query on the seqs collection using mongo operators
			project_fields : dict, default = {'QUERY_DATA': 0}
				A mongo projection query. Fields to project are represented as keys, and values are 0 or 1 depending on whether to suppress or include field
			exps_query : dict, default = {}
				A mongo query on the exps collection using mongo operators
			include_original_ngs_seq : boolean,
				When true, will also include the RAW NGS read and sequence header document for each queried document
			limit : integer
				Sets a limit on the number of documents to return. A limit of 0 is the same as not setting a limit.

			Datatype returned to self.query_results : list of dicts
				This function returns a list of documents to the query results attribute, self.query_results

			Returns
			-------
			self

			Examples
			--------
			Return the CDR3 field from of 100 IMGT documents from experiments in project Demo that contain the VGENE IGHV1-3
<<<<<<< HEAD
=======

>>>>>>> 670fa3e5030646debf2cfc5f20d277782898924f
			>>> query_seqs_collection(analysis_name['IMGT'], seqs_query={'DATA.VREGION.VGENES':'IGHV1-3'}, exps_query={'PROJECT_NAME':'Demo'}, project_fields={'DATA.CDR3':1}, limit=100)
		'''
		self.query_results = None
		self.logged_output.write('The following query in function "Query_Seqs_Collection" was desired: {0}\n'.format(str({'experiments': exp_id, 'seqs_query': seqs_query, 'exps_query': exps_query})))
		db = self.db_path
		q = {f: v for f, v in seqs_query.iteritems()}  # make a copy of query
		fields_to_project = copy.deepcopy(project_fields)
		if fields_to_project:
			allTrue = 0
			counts_in = 0
			for p, v in fields_to_project.iteritems():
				if isinstance(v, dict):
					pass
				elif v == 0:
					counts_in += 1
				else:
					counts_in += 1
					allTrue += 1

			if allTrue == counts_in:
				# The user projected fields using '1'. they listed select fields to project
				allTrue = True
				fields_to_project['_id'] = 1  # Also always project '_id' field
				fields_to_project[idIdentifier] = 1  # Also always project 'SEQ_ID' field
				fields_to_project[expIdentifier] = 1  # Also always project 'EXP_ID' field
				fields_to_project['ANALYSIS_NAME'] = 1  # Also always project 'EXP_ID' field
				fields_to_project['RECOMBINATION_TYPE'] = 1  # Also always project 'EXP_ID' field
			elif allTrue == 0:
				# The user suppressed feilds using '0'
				allTrue = False
				fields_to_project['QUERY_DATA'] = 0
		else:
			allTrue = True
			fields_to_project = None

		if exp_id is None:
			# If exp_id is passed in as empty, then make it self.allowed_experiment_ids
			exp_id = self.allowed_experiment_ids
		elif not isinstance(exp_id, list):
			exp_id = [exp_id]
		if exps_query:
			# Query exp_collection
			# Result of query will be stored in class name 'self.query_results'
			self.query_exps_collection(exp_id, q={f: v for f, v in exps_query.iteritems()})
			if self.query_results:
				self.exp_metadata = {str(exp['_id']): exp for exp in self.query_results if exp['_id'] in self.allowed_experiment_ids}
				# Erase metadata query
				exp_id = [vals['_id'] for meta, vals in self.exp_metadata.iteritems()]
				self.query_results = None
			else:
				self.exp_metadata = {}
				# No experiments to searchby. no experiments were found !
				exp_id = []
		if analysis_name:
			if not isinstance(analysis_name, list):
				q['ANALYSIS_NAME'] = analysis_name
			elif len(analysis_name) == 1:
				q['ANALYSIS_NAME'] = analysis_name[0]
			else:
				q['ANALYSIS_NAME'] = {'$in': analysis_name}
		if recombination_type:
			if not isinstance(recombination_type, list):
				q['RECOMBINATION_TYPE'] = recombination_type
			elif len(recombination_type) == 1:
				q['RECOMBINATION_TYPE'] = recombination_type[0]
			else:
				q['RECOMBINATION_TYPE'] = {'$in': recombination_type}

		# Modify the query to improve success of a match
		# 1) redirect select fieldname to fields created by server to imporve queries
		# 2) ensure that the format of the values provided to the field matches the format of the field put in database
		q = Parse_Mongo_Query_Expression(q, dict_defining_value_transformation=fields_for_queries_seqs_collection, redirection_fields=redirect_seq_collection_fields, modify_query_values_to_follow_db_schema=self.modify_query_values_to_follow_db_schema, redirect_fields_for_improved_queries=self.redirect_fields_for_improved_queries)

		# We will always add a range query on '_id' for each query (this ensures that it only searches allowed experiments)
		# There are exp_ids passed in by users OR found by the metadata_query above
		# Therefore, find the intersection between results and allowed id's
		# THIS WILL ALSO CONSIDER ANY TOP_LEVEL EXP_ID requests in the query field
		q[expIdentifier] = GetIdIntersection([self.allowed_experiment_ids, exp_id], q.pop(expIdentifier, None))

		# Based on the actual experiment IDs used in query, figure out which fields will be projected by using analysis schema in metadata
		if isinstance(q[expIdentifier], dict):  # The query is dictionary using '$in' command
			possible_ids = q[expIdentifier]['$in']
		else:
			possible_ids = [q[expIdentifier]]
		possible_metadata_fields = [self.exp_metadata[str(id)] for id in possible_ids]
		[schema_fields, possible_analyses, possible_recombination_types] = Get_Schema_Details(possible_metadata_fields, fields_to_project, allTrue)

		if 'ANALYSIS_NAME' not in q and possible_analyses != []:
			if len(possible_analyses) == 1:
				q['ANALYSIS_NAME'] = possible_analyses[0]
			else:
				q['ANALYSIS_NAME'] = {'$in': possible_analyses}
		if 'RECOMBINATION_TYPE' not in q and possible_recombination_types != []:
			if len(possible_recombination_types) == 1:
				q['RECOMBINATION_TYPE'] = possible_recombination_types[0]
			else:
				q['RECOMBINATION_TYPE'] = {'$in': possible_recombination_types}

		self.logged_output.write('Query was modified to: {0}\n'.format(str(q)))
		self.logged_output.write('The following values will be projected: {0}\n'.format(str(fields_to_project)))
		self.delim_file_headers = schema_fields
		self.query_results = db.seqs.find(q, fields_to_project).limit(limit)

		if include_original_ngs_seq:
			self.delim_file_headers.extend(data_fields_for_raw_data)
			# For every result in the query, also get the raw sequence data associated with result
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

	# ################Cursor Modification Functions######################
	'''
		The following functions are used to actually modify the format of the document coming out of the database. For example, it could
		convert documents be renaming keys, or it could convert documents into strings (i.e. json.dumps), or it could

	'''
	# Creates a generator for going through the cursor results and grouping results by seq_id
	# If same_document, then all results are merged into one document, if False, then only the information from @SEQ is appended to each dcoument
	def group_documents_by_seq_id(self, limit=None, same_document=False):
		def group_by_id(results):
			reported_docs = 0
			seq_id = ''
			output_doc = {}
			for doc_results in results:
				if doc_results[idIdentifier] != seq_id:
					# New doc found
					if output_doc:
						# Output all results to the same result/document/same row
						if same_document:
							yield output_doc
						# Split up results by annotation
						else:
							sub_doc = output_doc.pop('ANALYSES', None)
							if sub_doc:
								# No annotation results found
								output_doc.pop('ANALYSIS_NAME', None)
								output_doc.pop('DATE_UPDATED', None)
								output_doc.pop('_id', None)
								for analysis, each_doc in sub_doc.iteritems():
									each_doc.update(output_doc)
									yield each_doc
							else:
								yield output_doc
						reported_docs += 1
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
				if 'ANALYSIS_NAME' in doc and doc['ANALYSIS_NAME'] != seqRawData and idIdentifier in doc:
					raw_data = self.db_path.seqs.find({'_id': doc[idIdentifier]}, {'DATA': 1}).next()  # {idIdentifier:doc[idIdentifier],'ANALYSIS_NAME':seqRawData},{'DATA':1})
					if 'DATA' in raw_data:
						doc.update(raw_data['DATA'])
						yield doc
					else:
						yield doc
				else:
					yield doc

		self.query_results = get_data(self.query_results)
		return self

	# This generator function will take a list of docs and yield them one by one
	def list_to_gen(self):
		def make_gen(query_data):
			for doc in query_data:
				yield doc
		self.query_results = make_gen(self.query_results)
		return self

	# #################RETURNING QUERIES AS A FILE FUNCTIONS ####################
	'''
		The following functions are used when you want to save a query result to a file
		We have functions for saving results as :
			TAB, CSV, JSON, FLATJSON, FASTA, and FASTQ
			See the function allowed_file_types to see a description of each

		When running these functions the documents will be converted to strings and the cursor will be extinguished. Decided not to use itertools tee because it might have been an unnecessary use of memory.
	'''

	def save_as_delim(self, delimiter, prefix_file_path=None, file_ext='txt', header_var=[], keep_all_info=False, split_results_by=[], save_by_filename=False, save_by_exp_name=False, avoid_commas_in_string_output=False):
		'''
			Exports the results of a query to a delimited file. Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			delimiter : character
				This defines what character should be used to delimit file
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			header_var : list of strings, dict, or list of tuples default = []
				This variable allows user to manually define exactly how they want fields to appear in the delimited file.
				When format is a list of strings: Fields will appear in order based on the list of strings. i.e. ['SEQUENCE','SEQUENCE_HEADER']
				When format is a dict: Fields will be renamed based on key value pairs where key is field name and value is new field name. i.e. {'SEQUENCE':'raw read','SEQUENCE_HEADER':'miseqread'}
				When format is a list of tuples: Fields will appear in order based on list, and also be renamed based on (fieldname, new field name). i.e. [('SEQUENCE_HEADER','miseqread'),('SEQUENCE','raw read')]
			keep_all_info : boolean, default = False
				This variable defines whether to report all results from a query/projection or only report those fields defined in the variable 'header_var'.
				If False, then only those fields explicitly reported by 'header_var' will appear in the file. If True, then all fields not included in 'header_var' will
				appear afterwards in the order defined by 'default_sorting_order'
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
			avoid_commas_in_string_output : boolean, default = False
				When set to True, this will try to avoid using commas to convert fields to strings. For examples, lists will not be saved as strings using ','.join. Instead it will use '|'.join.
		'''
		self.to_file = True

		# This variable will define the order fields should appear in the delimited file
		output_sorted_fields = OrderedDict()

		# LETS TRY AND GET THE ORDER OF THE COLUMNS CORRECT
		# ACOUNT FOR USER REQUESTS, DEFAULT SORT ORDERS, AND FOUND FIELDS FROM QUERIES#
		user_requested_fields = True
		# First lets figure out what type of variable the user passed in as the header_var
		if header_var == [] or not(header_var):
			# User does not care about order of results, so we  will use the sorting order , we will also force keep_all_info to be true
			header_var = default_sorting_order
			keep_all_info = True
			user_requested_fields = False
			header_type = 'list'
		elif isinstance(header_var, basestring):
			# make it a list
			header_var = [header_var]
			header_type = 'list'
		elif isinstance(header_var, dict):
			# user will be renaming the field names based on key-value pairs
			header_type = 'dict'
		elif isinstance(header_var, list):
			# It is either a list or a list of tuples
			first_val = type(header_var[0])
			if isinstance(header_var[0], basestring):
				header_type = 'list'
			elif isinstance(header_var[0], type(('db_field', 'new_field'))):
				header_type = 'tuple'
			else:
				raise Exception('The parameter, header_var, can only be a list of strings or tuples.')
			for all_vals in header_var:
				if not(isinstance(all_vals, first_val)):
					raise Exception('The parameter, header_var, can only be a list of strings or list of tuples. It can not have multiple types within the list')
		else:
			raise Exception('The parameter, header_var, can only be a string, list, tuples or dictionary of strings')

		# Uppercase all field names in the in the outputfiles except for '_id'
		self.delim_file_headers = [s.upper() if s != '_id' else s for s in self.delim_file_headers]

		append_id = True if idIdentifier in self.delim_file_headers else False

		# Based on the header_type (list , dict or tuples), determine how the order of fields should appear
		if header_type == 'list':
			# user passed in a list defining the order to output fields. Because its a list, the NAMES of each field will not change
			header_var = [h.upper() for h in header_var]
			if append_id and idIdentifier not in header_var:
				header_var.insert(0, idIdentifier)
			if user_requested_fields:
				# So all fields in header_var will be output to file
				for fields_requested in header_var:
					output_sorted_fields[fields_requested] = fields_requested
			else:
				# If the user does not explicitely define which fields they want, then we will only output fields we know exist/were part of query
				for fields_requested in header_var:
					if fields_requested in self.delim_file_headers:
						output_sorted_fields[fields_requested] = fields_requested

		elif header_type == 'dict':
			# User passed in a dict so they do not want to specify the order of the fields, but they do want to rename fieldnames
			header_var = {h.upper(): nv if h != '_id' else h for h, nv in header_var.iteritems()}
			if idIdentifier not in header_var and append_id is True:
				header_var[idIdentifier] = idIdentifier

			# Cannot control order of fields, so use the order defined in default_sorting_order
			sub_dict = {}
			max_pos = 0
			new_requested_fields = []
			for each_requested_key in header_var:
				if each_requested_key in default_sorting_order:
					sub_dict[each_requested_key] = default_sorting_order.index(each_requested_key)
				else:
					new_requested_fields.append(each_requested_key)
			max_pos = max(sub_dict.values()) + 1 if sub_dict else 0
			# Any keys that were not found in our defualt variable get appended to dictionary
			for each_requested_key in new_requested_fields:
				sub_dict[each_requested_key] = max_pos
				max_pos += 1
			# Add in the keys sorted by their psostiion to the ordered dictionary
			for key, value in sorted(sub_dict.iteritems(), key=lambda (k, v): (v, k)):
				output_sorted_fields[key.upper()] = header_var[key]

		elif header_type == 'tuple':
			# User passed in a tuple so we can specify the order of the fields and rename fields based on user response
			# for each tuple, element 0 => name of field in document, element 1 => what to rename the field to
			header_var = [(row[0].upper(), row[1]) if row[0] != '_id' else row[0] for row in header_var]
			if append_id:
				found = False
				for i in header_var:
					if i[0] == idIdentifier:
						found = True
						break
				if found is False:
					header_var.insert(0, (idIdentifier, idIdentifier))
			for fields_requested in header_var:
				output_sorted_fields[fields_requested[0]] = fields_requested[1]

		# Field names desired by user have been defined; if keep_all_info is true then we need to go through the fields produced by the query
		# If they are not currently found in output_sorted_fields, then add them to the fields to be output
		if keep_all_info:
			sub_dict = {}
			new_requested_fields = []
			# loop through fields generated by query
			for generated_fields in self.delim_file_headers:
				# This field is not currently defined in output_sorted_fields
				if generated_fields not in output_sorted_fields:
					if generated_fields in default_sorting_order:
						sub_dict[generated_fields] = default_sorting_order.index(generated_fields)
					else:
						new_requested_fields.append(generated_fields)
			max_pos = max(sub_dict.values()) + 1 if sub_dict else 0

			# Any keys that were not found yet gets appended to dictionary
			for each_requested_key in new_requested_fields:
				sub_dict[each_requested_key] = max_pos
				max_pos += 1
			# Add in the keys sorted by their psostiion to the ordered dictionary
			for key, value in sorted(sub_dict.iteritems(), key=lambda (k, v): (v, k)):
				output_sorted_fields[key] = key
		# The order of all fields to output to the file has now been determined

		# The user did not define a prefix path
		if prefix_file_path is None:
			prefix_file_path = self.default_filename
		prefix_file_path = os.path.normpath(prefix_file_path)

		# If the query is just a dict (i.e. comes from a find_one command)
		if isinstance(self.query_results, dict):
			self.query_results = [self.query_results]

		# NOW define a function to control how we will save results to a file
		def generate_delim(query_results, output_sorted_fields):
			'''
				A function for actually parsing the cursor and saving as text.
				query_results represents a cursor from the mongodb query.
			'''
			# First set the default header line for files in this section
			self.files_created.add_header_row(output_sorted_fields, delimiter)
			# Next if there are any files in the class, open them up
			self.files_created.open_files()
			current_counts = 0
			num_line_to_write_sim = 1
			keys_to_be_reported = output_sorted_fields.keys()
			if avoid_commas_in_string_output:
				schema_schema_output_fnc_for_delim = schema_fields_to_file_avoid_commas
			else:
				schema_schema_output_fnc_for_delim = schema_fields_to_file

			# Lets figure out how to split results into multiple files using fields from documents
			split_results = []
			if save_by_filename:
				split_results.append('FILENAME')
			for fields in split_results_by:
				if fields not in split_results:
					split_results.append(fields)

			for seq_document in query_results:
				current_counts += 1
				prefix_line = prefix_file_path + '.'
				if not isinstance(seq_document, dict):
					# If its not a dict, then we just need to write results to a file
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += str(seq_document) + '\n'
				else:
					# If it is a dict, then it should be a db document
					# Process document: (1) Flatten dictionary, (2) cconvert all fields to strings
					seq_document = defaultdict(str, Process_Cursor_For_Output_File(seq_document, schema_schema_output_fnc_for_delim))
					added_prefix = ''
					if save_by_exp_name and expIdentifier in seq_document:
						added_prefix += self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'] + '.'
					for f in split_results:
						added_prefix += seq_document[f] + '.' if f in seq_document else '.'
					if added_prefix:
						prefix_line += re.sub("[ ,\\/]", '_', added_prefix)
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += (delimiter).join([seq_document[k] for k in keys_to_be_reported]) + '\n'
				if current_counts % num_line_to_write_sim == 0:
					self.files_created.write_data()
			self.files_created.write_data()
			self.files_created.close_files()
			self.query_results = None

		# First we need to see if query_results is an iterable
		if isinstance(self.query_results, (dict, basestring, int, long, float, complex)):
			# Just convert any single dicts/strings/numbers into lists
			self.query_results = [self.query_results]
		# Finally, save results to a file
		generate_delim(self.query_results, output_sorted_fields)

	def save_as_tab(self, prefix_file_path=None, file_ext='txt', header_var=[], keep_all_info=False, split_results_by=[], save_by_filename=False, save_by_exp_name=False):
		'''
			Exports the results of a query to a TAB delimited file. Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			header_var : list of strings, dict, or list of tuples default = []
				This variable allows user to manually define exactly how they want fields to appear in the delimited file.
				When format is a list of strings: Fields will appear in order based on the list of strings. i.e. ['SEQUENCE','SEQUENCE_HEADER']
				When format is a dict: Fields will be renamed based on key value pairs where key is field name and value is new field name. i.e. {'SEQUENCE':'raw read','SEQUENCE_HEADER':'miseqread'}
				When format is a list of tuples: Fields will appear in order based on list, and also be renamed based on (fieldname, new field name). i.e. [('SEQUENCE_HEADER','miseqread'),('SEQUENCE','raw read')]
			keep_all_info : boolean, default = False
				This variable defines whether to report all results from a query/projection or only report those fields defined in the variable 'header_var'.
				If False, then only those fields explicitly reported by 'header_var' will appear in the file. If True, then all fields not included in 'header_var' will
				appear afterwards in the order defined by 'default_sorting_order'
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
		'''
		self.save_as_delim('\t', prefix_file_path, file_ext, header_var, keep_all_info, split_results_by, save_by_filename, save_by_exp_name)

	def save_as_csv(self, prefix_file_path=None, file_ext='txt', header_var=[], keep_all_info=False, split_results_by=[], save_by_filename=False, save_by_exp_name=False):
		'''
			Exports the results of a query to a CSV file. Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			header_var : list of strings, dict, or list of tuples default = []
				This variable allows user to manually define exactly how they want fields to appear in the delimited file.
				When format is a list of strings: Fields will appear in order based on the list of strings. i.e. ['SEQUENCE','SEQUENCE_HEADER']
				When format is a dict: Fields will be renamed based on key value pairs where key is field name and value is new field name. i.e. {'SEQUENCE':'raw read','SEQUENCE_HEADER':'miseqread'}
				When format is a list of tuples: Fields will appear in order based on list, and also be renamed based on (fieldname, new field name). i.e. [('SEQUENCE_HEADER','miseqread'),('SEQUENCE','raw read')]
			keep_all_info : boolean, default = False
				This variable defines whether to report all results from a query/projection or only report those fields defined in the variable 'header_var'.
				If False, then only those fields explicitly reported by 'header_var' will appear in the file. If True, then all fields not included in 'header_var' will
				appear afterwards in the order defined by 'default_sorting_order'
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
		'''
		self.save_as_delim(',', prefix_file_path, file_ext, header_var, keep_all_info, split_results_by, save_by_filename, save_by_exp_name, avoid_commas_in_string_output=True)

	def save_as_flatjson(self, prefix_file_path=None, file_ext='json', convert_values_to_strings=False, split_results_by=[], save_by_filename=False, save_by_exp_name=False):
		'''
			Exports the results of a query to a FLAT JSON document
			Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			convert_values_to_strings : boolean, default=False
				When True, all values in from a document query will be converted to strings. For example fields contaning lists (i.e. VREGION.VGENES) will be come strings using ','.join(list)
				When False, datatypes are maintained
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
		'''
		self.to_file = True
		# The user did not define a prefix path
		if prefix_file_path is None:
			prefix_file_path = self.default_filename
		prefix_file_path = os.path.normpath(prefix_file_path)
		# JSON files do not have header rows
		self.files_created.add_header_row('')

		# Lets figure out how to split results into multiple files using fields from documents
		split_results = []
		if save_by_filename:
			split_results.append('FILENAME')
		for fields in split_results_by:
			if fields not in split_results:
				split_results.append(fields)

		def generate_flat(query_results):
			'''
				A function for parsing the cursor and saving as text.
				query_results represents a cursor from the mongodb query.
			'''
			num_line_to_write_sim = 1
			current_counts = 0
			for seq_document in query_results:
				current_counts += 1
				prefix_line = prefix_file_path + '.'
				if not isinstance(seq_document, dict):
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					try:
						self.files_created[prefix_line]['lines'] += json.dumps(seq_document) + '\n'
					except:
						self.files_created[prefix_line]['lines'] += json.dumps(str(seq_document)) + '\n'
				else:
					if convert_values_to_strings:
						# Process document: (1) Flatten dictionary, (2) cconvert all fields to strings
						seq_document = Process_Cursor_For_Output_File(seq_document, schema_fields_to_file)
					else:
						# Flatten dictionary only, convert object id's to strings
						seq_document = flatten_dictionary(Simple_Process_Output(seq_document))
					added_prefix = ''
					if save_by_exp_name and expIdentifier in seq_document:
						added_prefix += self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'] + '.'
					for f in split_results:
						added_prefix += seq_document[f] + '.' if f in seq_document else '.'
					if added_prefix:
						prefix_line += re.sub("[ ,\\/]", '_', added_prefix)
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += json.dumps(seq_document) + '\n'
				if current_counts % num_line_to_write_sim == 0:
					self.files_created.write_data()

			self.files_created.write_data()
			self.files_created.close_files()
			self.query_results = None

		# First we need to see if query_results is an iterable
		if isinstance(self.query_results, (dict, basestring, int, long, float, complex, list)):
			# Just convert any single dicts/strings/numbers into lists
			self.query_results = [self.query_results]
		# Finally, save results to a file
		generate_flat(self.query_results)

	def save_as_json(self, prefix_file_path=None, file_ext='json', split_results_by=[], save_by_filename=False, save_by_exp_name=False):
		'''
			Exports the results of a query as JSON. JSON dumps on all documents.
			Once this function is called, the query cursor will be exhausted.

			.. note::
				Because dumps using json is faster than dumps using bson, all ObjectId values will be converted to a string

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
		'''
		self.to_file = True
		# The user did not define a prefix path
		if prefix_file_path is None:
			prefix_file_path = self.default_filename
		prefix_file_path = os.path.normpath(prefix_file_path)
		# JSON files do not have header rows
		self.files_created.add_header_row('')

		# Lets figure out how to split results into multiple files using fields from documents
		split_results = []
		if save_by_filename:
			split_results.append('FILENAME')
		for fields in split_results_by:
			if fields not in split_results:
				split_results.append(fields)

		def generate_jsons(query_results):
			'''
				A function for parsing the cursor and saving as text.
				query_results represents a cursor from the mongodb query.
			'''
			num_line_to_write_sim = 1
			current_counts = 0
			for seq_document in query_results:
				current_counts += 1
				prefix_line = prefix_file_path + '.'
				if not isinstance(seq_document, dict):
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					try:
						self.files_created[prefix_line]['lines'] += json.dumps(seq_document) + '\n'
					except:
						self.files_created[prefix_line]['lines'] += json.dumps(str(seq_document)) + '\n'
				else:
					seq_document = Simple_Process_Output(seq_document, False)
					added_prefix = ''
					if save_by_exp_name and expIdentifier in seq_document:
						added_prefix += self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'] + '.'
					for f in split_results:
						added_prefix += seq_document[f] + '.' if f in seq_document else '.'
					if added_prefix:
						prefix_line += re.sub("[ ,\\/]", '_', added_prefix)
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += json.dumps(seq_document) + '\n'
				if current_counts % num_line_to_write_sim == 0:
					self.files_created.write_data()

			self.files_created.write_data()
			self.files_created.close_files()
			self.query_results = None

		# First we need to see if query_results is an iterable
		if isinstance(self.query_results, (dict, basestring, int, long, float, complex, list)):
			# Just convert any single dicts/strings/numbers into lists
			self.query_results = [self.query_results]
		# Finally, save results to a file
		generate_jsons(self.query_results)

	def save_as_fasta(self, prefix_file_path=None, file_ext='fasta', sequence_key='SEQUENCE', header_var=['SEQUENCE_HEADER'], keep_all_info=False, split_results_by=[], save_by_exp_name=False, save_by_filename=False, include_header_row=True):
		'''
			Exports the results of a query to a FASTA file
			Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			sequence_key : string, default = 'SEQUENCE'
				This field will be used as the sequence for a fasta file (i.e. >header\nsequence)
			header_var : list of strings, default = ['SEQUENCE_HEADER']
				These fields will be used to form the sequence header for a faster file. All fields in this list will be seperated by | (i.e. >header field 1 | header field 2 | header field 3\nsequence)
			keep_all_info : boolean, default = False
				When True then all other fields not explicity defined in 'header_var' will be exported in json format to the sequence header (i.e. >header field 1 <{field 2: value, field 3: value}\nsequence)
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
			include_header_row : boolean, default = True
				When true, a header row, prepended by '#', is output to the top of each created file. All field names defined by header_var are reported in this line

			.. note::
				The field corresponding to the sequencey_key will always be removed from the header_var variable
		'''
		if not(header_var):
			# Just use a JSON document to dump results to sequence header; it won't look pretty but user didnt define a value to use a sequence header and changed default value
			keep_all_info = True
		elif not isinstance(header_var, list):
			header_var = [header_var]

		temp_header_var = header_var
		header_var = []
		for each_field in temp_header_var:
			# Always make sure fields do not start with 'DATA.' (we are removing in them during the function Process_Cursor...)
			if each_field.startswith('DATA.'):
				each_field = each_field[5:]
			if each_field not in header_var:
				header_var.append(each_field.upper())
		# Always remove sequence key from header_var. It shouldnt appear in header.
		if sequence_key in header_var:
			header_var.remove(sequence_key)
		self.to_file = True
		# The user did not define a prefix path
		if prefix_file_path is None:
			prefix_file_path = self.default_filename
		prefix_file_path = os.path.normpath(prefix_file_path)
		# JSON files do not have header rows
		if include_header_row:
			show_header_var = header_var
			show_header_var[0] = '#' + show_header_var[0]
			self.files_created.add_header_row(show_header_var, '|')

		# Lets figure out how to split results into multiple files using fields from documents
		split_results = []
		if save_by_filename:
			split_results.append('FILENAME')
		for fields in split_results_by:
			if fields not in split_results:
				split_results.append(fields)

		def generate_fasta(results):
			'''
				A function for parsing the cursor and saving as FASTA text.
				query_results represents a cursor from the mongodb query.
			'''
			current_counts = 0
			num_line_to_write_sim = 1
			for seq_document in results:
				current_counts += 1
				prefix_line = prefix_file_path + '.'
				if not isinstance(seq_document, dict):
					# We really wouldnt like to export these fields as fasta, but just for edge edge case scenarios
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += '>Unexpected FASTA format\n{0}\n'.format(str(seq_document))
				else:
					# Process document: (1) Flatten dictionary, (2) cconvert all fields to strings
					seq_document = Process_Cursor_For_Output_File(seq_document, schema_fields_to_file)
					added_prefix = ''
					if save_by_exp_name and expIdentifier in seq_document:
						added_prefix += self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'] + '.'
					for f in split_results:
						added_prefix += seq_document[f] + '.' if f in seq_document else '.'
					if added_prefix:
						prefix_line += re.sub("[ ,\\/]", '_', added_prefix)
					prefix_line += file_ext
					if idIdentifier in seq_document:
						# We always want exp_id and seq_id to appear in the header of fasta file, but they will appear as part of json doc
						seq_id = seq_document.pop(idIdentifier)
						added_data = {idIdentifier: seq_id}
					else:
						added_data = {}
					# header = seq_document.pop(header_key) if header_key in seq_document else seq_id
					seq_val = seq_document.pop(sequence_key) if sequence_key in seq_document else ''  # check to see whether this document has a sequence field defined by sequence_key
					# Generate sequence header use fields defined in header_var
					header_str = seq_document.pop(header_var[0]) if header_var[0] in seq_document else ''
					for ind in range(1, len(header_var)):
						header_str += '|' + seq_document.pop(header_var[ind]) if header_var[ind] in seq_document else '|'
					if keep_all_info:
						# When True, any fields not explicitly defined in header_var will be shown as a json string in header
						added_data.update(seq_document)
					descriptor_str = header_str + fasta_file_delimiter + json.dumps(added_data) if added_data else header_str
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += '>{0}\n{1}\n'.format(descriptor_str, seq_val)
				if current_counts % num_line_to_write_sim == 0:
					self.files_created.write_data()
			self.files_created.write_data()
			self.files_created.close_files()
			self.query_results = None

		# Make sure we can iterate through entire document results (that is if any of the queries are NOT cursors, then make them lists)
		if isinstance(self.query_results, (dict, basestring, int, long, float, complex, list)):
			self.query_results = [self.query_results]
		generate_fasta(self.query_results)

	def save_as_fastq(self, prefix_file_path=None, file_ext='fastq', sequence_key='SEQUENCE', header_var=['SEQUENCE_HEADER'], keep_all_info=False, split_results_by=[], save_by_exp_name=False, save_by_filename=False, include_header_row=False, null_quality=30):
		'''
			Exports the results of a query to a FASTQ file
			Once this function is called, the query cursor will be exhausted.

			Parameters
			----------
			prefix_file_path : string , default=None
				All exported files will start with this path. The prefix should also include parent folder paths (i.e scratch/query/myquery)
				If None, then function will default prefix as, 'current folder/'+default_query_name
			file_ext : string, default='txt'
				Save files using this extension
			sequence_key : string, default = 'SEQUENCE'
				This field will be used as the sequence for a fasta file (i.e. >header\nsequence)
			header_var : list of strings, default = ['SEQUENCE_HEADER']
				These fields will be used to form the sequence header for a faster file. All fields in this list will be seperated by | (i.e. >header field 1 | header field 2 | header field 3\nsequence)
			keep_all_info : boolean, default = False
				When True then all other fields not explicity defined in 'header_var' will be exported in json format to the sequence header (i.e. >header field 1 <{field 2: value, field 3: value}\nsequence)
			split_results_by : list of strings, default = []
				This variable defines how to seperate results into multiple files by specific field names. For example, if you wanted each file to only include sequences with
				a unique CDR3 AA sequence, you would set split_results_by = ['CDR3.AA']. In this sitation, each new CDR3 instance would create a new file with the path: filepath+'.'+cdr3_name
			save_by_filename : boolean, default = False
				When true, results are seperated into different files based on the original filename of the sequences uploaded to the database. For example, if R1/R2 files from MIseq were uploaded, this
				function will seperate sequences based on either R1/R2 file location.
			save_by_exp_name : boolean, default = False
				When true, results are seperated into different files based on their respective experiment name
			include_header_row : boolean, default = False
				When true, a header row, prepended by '#', is output to the top of each created file. All field names defined by header_var are reported in this line
			null_quality : char or integer, default = 30
				This is used for documents that do not have fields for the sequence quality score.
				For example if sequences from a FASTA file were inserted into the database, then they have no quality information. So when saving these sequences as a FASTQ file,
				we will use this null_quality value as the default value of quality for each base. i.e. when null_quality is 40, the quality information wil be ')' for each sequence. ('@ACTGG' will be have quality ')))))')
				.. note::
					Assumes Sanger phred scoring system

			.. note::
				The field corresponding to the sequencey_key and the 'QUALITY_SCORE' field will always be removed from the header_var variable
		'''
		min_valid_quality_scores = 0
		max_valid_quality_score = 93
		default_quality_score = 30
		# Sanger shift
		shift = 33
		if null_quality:
			if type(null_quality) is int:
				if null_quality < min_valid_quality_scores or null_quality > max_valid_quality_score:
					raise Exception('Default quality score must be between: {0} and {1} '.format(str(min_valid_quality_scores, max_valid_quality_score)))
				default_char = chr(null_quality + shift)
			elif (type(null_quality) is str or type(null_quality) is unicode) and len(null_quality) == 1:
				ascii = ord(null_quality)
				phred = ascii - shift
				if phred < min_valid_quality_scores or phred > max_valid_quality_score:
					raise Exception('Default quality score must be between: {0} and {1}. Parameter passed was an ascii value of {2}'.format(str(min_valid_quality_scores), str(max_valid_quality_score), str(ascii)))
				default_char = null_quality
			else:
				raise Exception('Invalid ascii value passed to null_quality parameter')
		else:
			default_char = chr(default_quality_score + shift)

		if not(header_var):
			# Just use a JSON document to dump results to sequence header; it won't look pretty but user didnt define a value to use a sequence header and changed default value
			keep_all_info = True
		elif not isinstance(header_var, list):
			header_var = [header_var]
		temp_header_var = header_var
		header_var = []
		for each_field in temp_header_var:
			# Always make sure fields do not start with 'DATA.' (we are removing in them during the function Process_Cursor...)
			if each_field.startswith('DATA.'):
				each_field = each_field[5:]
			if each_field not in header_var:
				header_var.append(each_field.upper())
		# Always remove sequence key from header_var. It shouldnt appear in header.
		if sequence_key in header_var:
			header_var.remove(sequence_key)
		# Always remove quality score
		if 'QUALITY_SCORE' in header_var:
			header_var.remove('QUALITY_SCORE')
		self.to_file = True
		# The user did not define a prefix path
		if prefix_file_path is None:
			prefix_file_path = self.default_filename
		prefix_file_path = os.path.normpath(prefix_file_path)
		# JSON files do not have header rows
		if include_header_row:
			show_header_var = header_var
			show_header_var[0] = '#' + show_header_var[0]
			self.files_created.add_header_row(show_header_var, '|')

		# Lets figure out how to split results into multiple files using fields from documents
		split_results = []
		if save_by_filename:
			split_results.append('FILENAME')
		for fields in split_results_by:
			if fields not in split_results:
				split_results.append(fields)

		def generate_fastq(results, to_file, split_results_by):
			'''
				A function for parsing the cursor and saving as FASTQ text.
				query_results represents a cursor from the mongodb query.
			'''
			default_char_str = ''.join([default_char] * 500)
			char_str_len = len(default_char_str)
			current_counts = 0
			num_line_to_write_sim = 1
			for seq_document in results:
				current_counts += 1
				prefix_line = prefix_file_path + '.'
				if not isinstance(seq_document, dict):
					# We really wouldnt like to export these fields as fasta, but just for edge edge case scenarios
					prefix_line += file_ext
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += '@Unexpected FASTQ format\n{0}\n+\n{1}\n'.format(str(seq_document), ''.join([default_char] * len(str(seq_document))))
				else:
					# Process document: (1) Flatten dictionary, (2) cconvert all fields to strings
					seq_document = Process_Cursor_For_Output_File(seq_document, schema_fields_to_file)
					added_prefix = ''
					if save_by_exp_name and expIdentifier in seq_document:
						added_prefix += self.exp_metadata[seq_document[expIdentifier]]['EXPERIMENT_NAME'] + '.'
					for f in split_results:
						added_prefix += seq_document[f] + '.' if f in seq_document else '.'
					if added_prefix:
						prefix_line += re.sub("[ ,\\/]", '_', added_prefix)
					prefix_line += file_ext
					if idIdentifier in seq_document:
						# We always want exp_id and seq_id to appear in the header of fasta file, but they will appear as part of json doc
						seq_id = seq_document.pop(idIdentifier)
						added_data = {idIdentifier: seq_id}
					else:
						added_data = {}
					seq_val = seq_document.pop(sequence_key) if sequence_key in seq_document else ''  # check to see whether this document has a sequence field defined by sequence_key
					if 'QUALITY_SCORE' in seq_document:
						quality_score = seq_document.pop('QUALITY_SCORE')
					else:
						if len(seq_val) > char_str_len:
							default_char_str = ''.join([default_char] * len(seq_val))
							char_str_len = len(seq_val)
							quality_score = default_char_str
						else:
							quality_score = default_char_str[:len(seq_val)]
					# Generate sequence header use fields defined in header_var
					header_str = seq_document.pop(header_var[0]) if header_var[0] in seq_document else ''
					for ind in range(1, len(header_var)):
						header_str += '|' + seq_document.pop(header_var[ind]) if header_var[ind] in seq_document else '|'
					if keep_all_info:
						# When True, any fields not explicitly defined in header_var will be shown as a json string in header
						added_data.update(seq_document)
					descriptor_str = header_str + fasta_file_delimiter + json.dumps(added_data) if added_data else header_str
					self.files_created[prefix_line]['doc_count'] += 1
					self.files_created[prefix_line]['lines'] += '@{0}\n{1}\n+\n{2}\n'.format(descriptor_str, seq_val, quality_score)
				if current_counts % num_line_to_write_sim == 0:
					self.files_created.write_data()
			self.files_created.write_data()
			self.files_created.close_files()
			self.query_results = None

		# Make sure we can iterate through entire document results (that is if any of the queries are NOT cursors, then make them lists)
		if isinstance(self.query_results, (dict, basestring, int, long, float, complex, list)):
			self.query_results = [self.query_results]
		generate_fastq(self.query_results, self.to_file, split_results_by)


class Write_Files(defaultdict):
	'''
		Make a special class for handling dicts of file buffers.
		This will define a default dict to allow opening files not yet created
	'''
	def add_header_row(self, header_row, delim=''):
		'''
			Define what should be the header_row for a file
		'''
		if isinstance(header_row, list):
			self.header_row = delim.join(header_row)
		elif isinstance(header_row, basestring):
			self.header_row = header_row
		elif isinstance(header_row, dict):
			self.header_row = delim.join(header_row.values())
		else:
			raise Exception(str(header_row) + ' has an unknown structure')

	def __missing__(self, key):
		'''
			Create a default dict definition. If a string is not present in the dictionary, then add it to dict and create a new file
		'''
		temp = open(key, 'w')
		try:
			if self.header_row:
				temp.write(self.header_row + '\n')
		except:
			self.header_row = ''
		self[key] = {'buffer': temp, 'doc_count': 0, 'lines': ''}
		return self[key]

	def get_files(self):
		return [{'filepath': key, 'doc_count': self[key]} for key in self]

	def write_data(self):
		'''
			Write all data to file
		'''
		for key in self:
			self[key]['buffer'].write(self[key]['lines'])
			self[key]['lines'] = ''

	def close_files(self):
		'''
			Close all open files
		'''
		for key in self:
			self[key]['buffer'].close()

	def open_files(self):
		'''
			Open all files. Use append to avoid overwriting
		'''
		for key in self:
			self[key]['buffer'] = open(key, 'a')
