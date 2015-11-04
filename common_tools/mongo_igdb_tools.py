'''
	Module contains all tools required for interacting with mongo implementation of our ig-db
'''
import pymongo
from bson.json_util import loads as bson_loads
import json
from bson.objectid import ObjectId


def connectToIgDatabase(username, password, dbpath='biotseq.icmb.utexas.edu', db_instance='appsoma', db_port=27017):
	'''
		This function is used to connect to a mongo database instance

		Input
		-----
		username : string
			String defining the mongodb user that has READ access to the following collections: users, seqs, exps
		password : string
			String defining the password for that username
		dbpath : string, default = 'biotseq.icmb.utexas.edu'
			String defining the url location of the mongodb
		db_instance : string, default = 'appsoma'
			The name of the database containing the ig-db
		db_port : int, default = 27017
			Port that the mongodb instance listens on

		.. important::
			This function should be used only by administrators who have access to the database. It will return a variable that grants
			full access to the database and collections the user has access to.

	'''
	connection = pymongo.MongoClient(host=dbpath, port=db_port)
	db = connection[db_instance]
	db.authenticate(username, password)
	return db


def convert_to_objectid(id):
	'''
		This function will convert a string or properly formatted dictionary into a MongoDB object ID.

		There could be a few ways that an object ID was mal-transformed while querying the database:
		Examples:
			1) str(ObjectId) will convert into a string. In order to convert it back, we simply do ObjectID(str)
			2) bson_dumps(ObjectId) will create a bson dumped string variable of an object ID. As a string it looks like '{$oid: id }'
				But, if we did json.loads on that value, then we dont get an objectId. Instead we get a dictionary {'$oid': id}

		So this function will attempt to convert such examples back into an ObjectId variable

		Inputs
		------
		id: string, dict, or ObjectId variable

		Outputs
		-------
		id converted into an ObjectId variable
	'''

	try:
		# the ObjecId was dumped using bson dumps, BUT JSON was used to load. so reload it using bson loads
		if isinstance(id, dict):
			# the wrong funciton was used to load this variable:
			# json.loads of an object id results in a dictionary where keys are {'oid':}
			# therefore, we need to first json.dumps back to a string, and THEN loads using bson loads			
			return bson_loads(json.dumps(id))

		elif isinstance(id, basestring):
			if '$oid' in id: 
				# the ObjecId was dumped using bson dumps, so reload it using bson loads
				# if bson_dumps is used,but it is not turned into dict using josn.loads. 
				# instead it is just a string. So we need to convert a bson_dumps string to object id, then we need to slightly modify it..hack for special situations			
				return bson_loads(id.strip())
			else:
				return ObjectId(id.strip())
		else:
			return ObjectId(id)
	except Exception as e:
		raise Exception('There is something wrong with the provided id: {0}\nThe following error was reported: {1}'.format(str(id), str(e)))


def convert_text_to_index_field_text(value):
	"""
		Input
		-----
		value: String or list of strings

		Output
		------
		String or list of strings with the strings converted as such:
			eliminate any whitespaces, eliminate any [',',';','/',':','-'], and convert the remaining characters to lowercase
	"""

	def _text_conversion_to_index_text(text):
		replace_these_chars = ['_', ',', ';', '-', '/', ':', ' ']
		val = "".join(text.split()).lower()  # changed to lowercase
		for r in replace_these_chars:
			val = val.replace(r, '')
		return val
	return [_text_conversion_to_index_text(str(item)) for item in value] if isinstance(value, list) else _text_conversion_to_index_text(str(value))


def default_fields_data_types(x):
	'''
		This function will create a default method for transforming fields in SEQS collection we dont account for in our schema.
		It checks to see whether the variable is an allowed data type (see allowed_data_types in function)

		If its a list, then all values in x are uppercased if they are strings
		All other non-string values are not modified

		Input
		-----
		x: A variable that is used in a database query

		Output:
		x variable converted to a default format
	'''
	allowed_data_types = [float, int, str, unicode, bool, long]
	if type(x) is list:
		x = []
		for e in x:
			if isinstance(e, basestring):
				x.append(str(e).strip().upper())
			elif type(e) in allowed_data_types:
				x.append(e)
		if not x:
			return None
		else:
			return x
	else:
		if isinstance(x, basestring):
			return str(x).strip().upper()
		elif type(x) not in allowed_data_types:
			return None
		else:
			return x


def AttemptToConvertToBool(x):
	'''
		This function will attempt to convert a variable into a boolean value.
		The following strings are converted to True: 'yes', 'true'
		The following strings are converted to False: 'no', 'false'

		Input
		-----
		x: A variable to convert to boolean

		Output
		------
		x: A successfully convereted boolean; if unsuccessful then None
	'''

	if type(x) is bool:
		return x
	elif type(x) is int:
		return bool(x)
	elif str(x).strip().lower() == 'false' or str(x).strip().lower() == 'no':
		return False
	elif str(x).strip().lower() == 'true' or str(x).strip().lower() == 'yes':
		return True
	else:
		return None


def default_metadata_fields(x):
	'''
		This function will create a default method for transforming fields in EXPS collection we dont account for in our schema.

		Input
		-----
		x: A variable that is used in a database query

		Output:
		x variable converted to a default format
	'''

	# for now, we do not modify any metadata fields that are not defined in exp_fields_and_data_types below
	return x
