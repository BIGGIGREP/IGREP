
#import immunogrep_igdbtools as igdbtools
from datetime import datetime
import json
import immunogrep_useful_functions as useful_functions
from immunogrep_global_variables import fasta_file_delimiter
from immunogrep_global_variables import idIdentifier
from bson.objectid import ObjectId
import re 
from bson.json_util import loads as bson_loads

database_schema_warehouse_file = "https://biotseq.ut.appsoma.com:5001/public_folder/gg_form_fields.txt"

allowed_chain_types = ['heavy','light','alpha','beta']

from collections import defaultdict


def convert_to_objectid(id):
	try:			
		#the ObjecId was dumped using bson dumps, BUT JSON was used to load. so reload it using bson loads
		if isinstance(id,dict):# and '$oid' in id:			
			#the wrong funciton was used to load this variable:
				#json.loads of an object id results in a dictionary where keys are {'oid':}
				#therefore, we need to first json.dumps back to a string, and THEN loads using bson loads			
			return bson_loads(json.dumps(id))  #ObjectId(id['$oid'])
		
		elif isinstance(id,basestring):					
			if '$oid' in id: 
			#the ObjecId was dumped using bson dumps, so reload it using bson loads	
				#if bson_dumps is used,but it is not turned into dict using josn.loads. instead it is just a string. So we need to convert a bson_dumps string to object id, then we need to slightly modify it..hack for special situations
				#id = id.replace("{","").replace("}","").replace('"','')
				#id = id.split(':')[1]
				#print id
				return bson_loads(id.strip())				
			else:
				return ObjectId(id.strip())
		else:			
			return ObjectId(id)
	except Exception as e:
		raise Exception('There is something wrong with the provided id: {0}\nThe following error was reported: {1}'.format(str(id),str(e)))


#this function will convert any metadata field that is 'indexed' into a format that is easier for queries 
#this is relevant for all fields within the 'DUPLICATED_FIELDS' document in the exps collection 
def convert_text_to_index_field_text(value):	
	"""
	Input: string or list of strings.
	Output: string or list of strings with the strings converted to the _INDEX field format.
	Eliminate any whitespace in data[field], eliminate any '_', and convert the remaining characters to lowercase. ==> originally uppercase
	"""
	def _text_conversion_to_index_text(text):		
		replace_these_chars = ['_',',',';','-','/',':',' ']
		val = "".join(text.split()).lower()#changed to lowercase
		for r in replace_these_chars:
			val = val.replace(r, '')
		return val
	return [_text_conversion_to_index_text(str(item)) for item in value] if isinstance(value, list) else _text_conversion_to_index_text(str(value))



#describes the structure of the mongo database

#default experiments collection
def Exps_Collection():
	
	OptionalFields = [
		"CELL_NUMBER",#:None, #number of cells sequenced, datatype=integer		
		"PAIRING_TECHNIQUE",#:None, #method of preparing library for pairing (usually RT step), datatype = string				
		"ANALYSIS_METHODS",#:None, #lists the types of analysis done on this experiment so far (i.e. IMGT IGBLAST, ETC), dataype=list of strings		
		"PUBLICATIONS",#:None, #list of references that use this dataset, datype = list of strings		
		"TARGET_READS",#:None, #requested targeted number of reads for sequencing		
		"PERSON_WHO_PREPARED_LIBRARY",#:None,
		"CELL_MARKERS_USED",#:None,
		"CONTAINS_RNA_SEQ_DATA",#:None,
		"PRIMER_SET_NAME",#:None,
		"LIST_OF_POLYMERASES_USED",#:None,
		"LAB_NOTEBOOK_SOURCE",#:None,
		"POST_SEQUENCING_PROCESSING:PHI_X_FILTER",#:None,
		"POST_SEQUENCING_PROCESSING:QUALITY_FILTER",#:None,
		"POST_SEQUENCING_PROCESSING:PROCESS_R1_R2_FILE",#:None,
		"MID_TAG",#:None
		"READ_ACCESS"
	]	
	
	RequiredFields =[# {
		"PROJECT_NAME",#: None, #defines the project (one project can have many experiments), datatype = string  // Formerly "PROJECT_ID"
		"SPECIES",#: None, #species sequenced, datatype=list of integers 
		"LAB",#: None, #describes teh laboratory containing the data
		"DESCRIPTION",#:None, #describes the experiment, datatype = list of strings  // Formerly "KEYWORDS"		
		"EXPERIMENT_NAME",#:None, #defines the experiment (the run from the sequencer), datatype = string  // Formerly "EXPERIMENT_ID"
		#"UNIQUE_IRODS_ID":None, #unique key for the experiment stored in IRODs, datatype = string (generated using timestamp)
		"SAMPLE_PREPARATION_DATE",#:None, #Date experiment performed, datatype = datetime  // Formerly "DATE"
		"SEQUENCING_PLATFORM",#:None, #platform used for sequencing, datatype = string
		#"READ_ACCESS",#this will no longer be required. this is because read acces by default becomes identical to owners of experiment  #:None, #list of users OR LABS allowed to read results from this experiment, datatype =list of strings
		"OWNERS_OF_EXPERIMENT",#:None,	#list of users allowed to write or update this experiment in the database with results, datatpe = list of string/user names  // Formerly "WRITE_ACCESS"
		"CHAIN_TYPES_SEQUENCED",#:None, #chain types sequenced, datatype = currently only a list of the following is allowed: [heavy, light,alpha, beta]
		"CELL_TYPES_SEQUENCED",#:None, #lists the types of immune cells sequenced, datatype = list of strings
		"ISOTYPES_SEQUENCED",#:None, #isotypes sequenced, datatype = list of strings
		"VH:VL_PAIRED",#:None, #whether or not sequences were paired, datatype = boolean  // Formerly "PAIRED"
		"TEMPLATE_TYPE",#:None,
		"REVERSE_PRIMER_USED_IN_RT_STEP",#:None		
	]
	#}	
	
	#USERS are not allowed to pass these values in to experiment metadata. 
	#IF passed in, client will remove values from dictionary
	FieldsAddedByClient = [
		"UPLOADED_BY",#:None,
		"EXPERIMENT_CREATION_DATE",#:None, #Date experiment was created in the database
		"CURATED",#determines whether an experiment has been curated by the curator or not 
		"SEQ_COUNT",#stores counts of all raw sequences in experiment
		'DUPLICATED_FIELDS',#fields stored in here have been dupolicated for imporved queries. these metadata have been formated to allow case-insensitive queries 
		"ANALYSES_COUNT",#stores the counts of all analyses types in experiment 
		#'ANALYSIS_TYPES', #stores a list of all analysis types/methods for a given experiment and the number of sequences with given analysis type 
		'ANALYSIS_SCHEMA', #stores internal data with regard to the analysis types stored in the experiment and the combined schema of all sequences within the experiment 
		"ANALYSES_SETTINGS" #stores internal data with regard to analyses settings used to annotate the sequences in the database. 
	]
	
	#we will need to rewrite a few fields so they can be used efficiently for indexing and queries downstream
	DuplicatedFieldsForIndexing = ['PROJECT_NAME','EXPERIMENT_NAME', 'PAIRING_TECHNIQUE','ISOTYPES_SEQUENCED','CHAIN_TYPES_SEQUENCED','CELL_TYPES_SEQUENCED','SPECIES','PUBLICATIONS']
	
	return {'required_fields': RequiredFields, 'optional_fields': OptionalFields, 'duplicated_fields': DuplicatedFieldsForIndexing,'protected_fields':FieldsAddedByClient}

def ModifyStrInList(x, modifyCasing=True,delim=','):	
	if modifyCasing:		
		return [value.strip().lower() for value in str(x).split(delim) if value.strip()] if type(x) is not list else [str(value).strip().lower() for value in x if str(value).strip()]
	else:
		return [value.strip() for value in str(x).split(delim) if value.strip()] if type(x) is not list else [str(value).strip() for value in x if str(value).strip()]	

def CapitalizeWordsInList(x, delim=','):	
	return [value.strip().title() for value in str(x).split(delim) if value.strip()] if type(x) is not list else [str(value).strip().lower() for value in x if str(value).strip()]
	
def RemoveLab(x):
	labs = ModifyStrInList(x,True,',')
	for i,l in enumerate(labs): #remove any references to 'lab'..dont use 'replace', just incase some name has 'lab' somewhere in it
		labs[i] = ' '.join([word for word in l.split(' ') if word.strip()!='lab'])
	return labs

def AttemptToConvertToBool(x):
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
	#for now, we do not modify any metadata fields that are not defined in exp_fields_and_data_types below 
	return x

#THE SERVER WIL MAINTAIN THE DATATYPSE FOR THE FOLLOWING FIELDS PASSED INTO AS METADATA TO EXPS COLLECTION
#ALL METADATA FIELDS WITH EACH KEY NAME ARE TREATED BY THE LAMBDA FUNCTIONS BELOW
#IF A KEY IS NOT PRESENT IN THE VARIABLE BELOW, THEN THE DATA WILL BE TRANSFORMED BY TEH FUNCTION default_metadata_fields ABOVE=> defaultdict 
exp_fields_and_data_types = defaultdict(lambda:default_metadata_fields,{
	'PROJECT_NAME':lambda x:str(x),
	'EXPERIMENT_NAME':lambda x:str(x),
	'SEQUENCING_PLATFORM':lambda x:str(x).lower(),
	'SPECIES': lambda x:' '.join([word.title() if i==0 else word.lower() for i,word in enumerate(str(x).split(' '))]), #capitalize first word only 
	'LAB': RemoveLab,# LowercaseStrInListByComma, #lower case all labs and store as array 
	'OWNERS_OF_EXPERIMENT': lambda x: ModifyStrInList(x,True,','),# LowercaseStrInListByComma, #lower case all owners and store as array 
	'READ_ACCESS': lambda x: ModifyStrInList(x,True,','),# LowercaseStrInListByComma, #lower case all read users and store as array 
	'CHAIN_TYPES_SEQUENCED': lambda x: ModifyStrInList(x,True,','),# LowercaseStrInListByComma, #lower case all chains types  and store as array 
	'CELL_TYPES_SEQUENCED': lambda x: ModifyStrInList(x,False,','), #lambda x: [cell.strip() for cell in str(x).split(',')] if type(x) is not list else [str(cell).strip() for cell in x], #dont modify casing of cell types 
	'ISOTYPES_SEQUENCED': lambda x: ModifyStrInList(x,False,','), # lambda x: [isotype.strip() for isotype in str(x).split(',')] if type(x) is not list else [str(isotype).strip() for isotype in x], #dont modify casing of isotypes 
	'DESCRIPTION': lambda x: ModifyStrInList(x,False,';'), #lambda x: [descriptor.strip() for descriptor in str(x).split(';')] if type(x) is not list else [str(descriptor).strip() for descriptor in x], #dont modify casing of descriptions
	'MID_TAG': lambda x: [bc.replace(' ','') for bc in ModifyStrInList(x,True,',')],# LowercaseStrInListByComma,#lowercase all DNA barcodes used in this xperiment, ALSO REMOVE ANY EMPTY SPACES
	'CELL_MARKERS_USED': lambda x: ModifyStrInList(x,True,','), # LowercaseStrInListByComma,#lowercase all markers in this xperiment 
	'PAIRING_TECHNIQUE':lambda x: str(x),#don't change anythign about pariing technique
	'PUBLICATIONS': lambda x: ModifyStrInList(x,False,';'),#lambda x:[pub.strip() for pub in str(x).split(';')] if type(x) is not list else [str(pub).strip() for pub in x], #dont modify casing of publications
	'CELL_NUMBER':lambda x:int(x),
	'TARGET_READS':lambda x:int(x),
	'CONTAINS_RNA_SEQ_DATA': AttemptToConvertToBool, #for boolean operations 
	'VH:VL_PAIRED': AttemptToConvertToBool, #for boolean operations 
	'LIST_OF_POLYMERASES_USED':lambda x: ModifyStrInList(x,True,','),# LowercaseStrInListByComma,#lowercase all polymerases 	
	'PRIMER_SET_NAME': lambda x: ModifyStrInList(x,True,','),#LowercaseStrInListByComma,
	"POST_SEQUENCING_PROCESSING:PHI_X_FILTER":AttemptToConvertToBool,
	"POST_SEQUENCING_PROCESSING:QUALITY_FILTER":lambda x: str(x),
	"POST_SEQUENCING_PROCESSING:PROCESS_R1_R2_FILE": lambda x: str(x),
	"PERSON_WHO_PREPARED_LIBRARY":lambda x: CapitalizeWordsInList(x,','),#Captialize each name/word ben smith => Ben Smith
	"REVERSE_PRIMER_USED_IN_RT_STEP":lambda(x):str(x).lower()	
})


#default sequences collection
def Seqs_Collection():
	OptionalFields={
		"QUALITY_SCORE":None, #quality score of the sequence
		"ANALYSIS": None, # subdocuments for each analysis type
		'FILENAME':None #the basename of the original file 
		#"PAIRED_ID":None
	}
	
	RequiredFields={
		"EXP_ID":None, #key or id corresponding to an experiment defined in the experiments collection, datattype = _id from exps collection (or string?)
		"SEQUENCE_HEADER":None, #sequence header of the sequence, datatype = string
		"SEQUENCE":None	#actual raw nucleotide sequence , datatype = string
	}
		
	
	IndexedFields = indexed_fields_in_analysis()
	
	return {'required_fields': RequiredFields, 'optional_fields': OptionalFields, 'indexed_fields': IndexedFields}

def default_fields_data_types(x):
	allowed_data_types = [float,int,str,unicode,bool,long]
	if type(x) is list:		
		x = []
		for e in x:
			if isinstance(e,basestring):
				x.append(str(e).strip().upper())
			elif type(e) in allowed_data_types:
				x.append(e)		
		if not x:			
			return None
		else:
			return x
	else:
		if isinstance(x,basestring):
			return str(x).strip().upper()
		elif type(x) not in allowed_data_types:
			return None			
		else:
			return x

def ConvertTo(x,convtype):	
	if x or x==False or x == 0 or x == 0.0:
		try:			
			if convtype == str or convtype == unicode or convtype == basestring:
				return str(x).strip().upper() #convert to str 
			else:				
				return convtype(x)#convert x to type defined by convtype	
			#return x	
		except Exception as e:
			#raise Exception(e)
			return None
	else:
		return x
def ConvertToPercent(x):
	if x == 0 or x == 0.0:
		return 0.0
	if x:
		try:
			return round(float(x),4)
		except:
			return None
	else:
		return None

def ConvertToList(x,convtype,delim):
	#here is a hack to catch outputs from CSV files 
	backup_delim = '|' 
	
	if x or x==False:
		try:
			if type(x) is not list:				
				if delim==',' and delim not in x: 	#here is a hack to catch outputs from CSV files 
					delim = backup_delim
				if convtype == str or convtype == unicode or convtype == basestring:				
					return [convtype(e).strip().upper() for e in str(x).split(delim) if convtype(e).strip()]
				else:
					return [convtype(e) for e in str(x).split(delim)]
			else:
				if convtype == str or convtype == unicode or convtype == basestring:					
					return [str(e).strip().upper() for e in x if str(e).strip()]
				else:					
					return [convtype(e) for e in x]
		except Exception as e:
			#raise Exception(e)
			return None				
	else:
		return []
		
def ConvertFromList(list_x,delim):
	if not type(list_x) is list:
		return ''
	
	str_results = [str(x) for x in list_x]
	return delim.join(str_results)	
	
def default_fields_to_file(x):
	if type(x) is list:
		return ','.join([str(e) for e in x])
	else:
		return str(x)
	
#create a list for allowing queries of V,D,AND J GENES and ALLEES
#input the current list of genes/alleles
#output a new list where each element is split by ' ', and '*', and '(' and ')'
def CreateGeneQueryList(genes):
	if type(genes) is not list:
		genes = [genes]			
	return genes_for_queries
	
def Allele_Name_To_List(allele_list):
	#take a list of alleles and create an array where each element has been split using characters below
	#therefore each split value will be stored in database
		#most important to split by => spaces and '*'. the '*' will seperate gene from allele IGHV1-3*01 => gene = IGHV1-3, allele = IGHV1-3*01
	
	#these characters may sometimes be present in the gene name, but they do not differentiate gene vs allele. we want to separate these possible delimiters into words
	delims = ' |;|,|\(|\)|\[|\]' 	
	gene_name = []
	for allele in allele_list:	
		allele = allele.upper()
		gene_name.append(allele)#store original gene name uppercased
		words = [x.strip().upper() for x in re.split(delims,allele) if x.strip()]
		possible_alleles = [w.split('*')[0] for w in words]  #if any of the words have an '*' then store text before the '*'...this is usually a gene name  
		gene_name.extend(possible_alleles+words) #store all 'words' and also any text before the '*'
	return list(set(gene_name))

	
#DEFAULT ANALYSIS SCHEMA
#maintain the proper format of the following fields in any schema
#'to_db': for each key, transforms the value using the defined function before inserting to database
#'to_file': for each key, transforms the value from the database to the defined after getting a query result 
#ANY KEYS NOT IN THE FOLLOWING VARIABLE WILL USE THE FUNCTIONS DEFINED BY TO_DB AND TO_FILE IN DEFAULTDICT EXPRESSION. 
schema_fields_and_data_types = defaultdict(lambda:{'to_db':default_fields_data_types,'to_file':default_fields_to_file},
{
	"COMMAND":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),				
	},
	"NOTES":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"PREDICTED_AB_SEQ.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"PREDICTED_AB_SEQ.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"STRAND":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"PREDICTED_CHAIN_TYPE":#lambda x:ConvertTo(x,str),#or list??
	{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"PRODUCTIVE":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"LOCUS":#lambda x:ConvertTo(x,str),#or list??
	{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.SHM.NT":{
		'to_db':lambda x:ConvertTo(x,float),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"VREGION.SHM.NT_PER":{
		'to_db':lambda x:ConvertToPercent(x),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"VREGION.SHM.AA":{
		'to_db':lambda x:ConvertTo(x,float),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"VREGION.SHM.AA_PER":{
		'to_db':lambda x:ConvertToPercent(x),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"JREGION.SHM.AA":{
		'to_db':lambda x:ConvertTo(x,float),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"JREGION.SHM.AA_PER":{
		'to_db':lambda x:ConvertToPercent(x),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"JREGION.SHM.NT":{
		'to_db':lambda x:ConvertTo(x,float),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	"JREGION.SHM.NT_PER":{
		'to_db':lambda x:ConvertToPercent(x),
		'to_file':lambda x:ConvertTo(x,str),
	},#,lambda x:ConvertTo(x,float),
	'VREGION.VGENE_QUERY_START':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str),
	},#lambda x:ConvertTo(x,int),
	'VREGION.VGENE_QUERY_END':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str),
	},#lambda x:ConvertTo(x,int),
	"VREGION.FR1.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.FR1.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.CDR1.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.CDR1.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.FR2.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.FR2.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.CDR2.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.CDR2.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.FR3.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.FR3.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"VREGION.VGENES":{
		'to_db':lambda x: ConvertToList(x,str,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	"VREGION.VGENE_SCORES": {
		'to_db':lambda x: ConvertToList(x,float,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	"CDR3.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"CDR3.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},	
	"DREGION.DGENES":{
		'to_db':lambda x: ConvertToList(x,str,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	"DREGION.DGENE_SCORES": {
		'to_db':lambda x: ConvertToList(x,float,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	"JREGION.FR4.NT":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"JREGION.FR4.AA":{
		'to_db':lambda x:ConvertTo(x,str),
		'to_file':lambda x:ConvertTo(x,str),
	},
	"JREGION.JGENES":{
		'to_db':lambda x: ConvertToList(x,str,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	"JREGION.JGENE_SCORES":{
		'to_db':lambda x: ConvertToList(x,float,','),
		'to_file':lambda x: ConvertFromList(x,',')
	},
	'JREGION.JGENE_QUERY_START':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str),
	},
	'JREGION.JGENE_QUERY_END':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str),
	},
	'CDR3.AA_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'CDR3.NT_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'VREGION.CDR1.AA_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'VREGION.CDR1.NT_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'VREGION.CDR2.AA_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'VREGION.CDR2.NT_LENGTH':{
		'to_db':lambda x:ConvertTo(x,int),
		'to_file':lambda x:ConvertTo(x,str)
	},
	'ISOTYPE.GENE':{
		'to_db':lambda x: ConvertToList(x,str,','),
	},
	'ISOTYPE.MISMATCHES':{
		'to_db':lambda x:ConvertToList(x,int,',')
	},
	'ISOTYPE.SCORES':{
		'to_db':lambda x:ConvertToList(x,int,',')
	},
	'ISOTYPE.PER_ID':{
		'to_db':lambda x:ConvertToList(x,float,',')
	},
	'PAIRING.CONFIDENCE':{
		'to_db':lambda x:ConvertTo(x,float)
	},
	'PAIRING.DOMINANCE':{
		'to_db':lambda x:ConvertTo(x,float)
	},
	
	'FILENAME':{
		'to_db':lambda x:str(x),#just dont do anything to field
		'to_file':lambda x:str(x)
	}
})

	
#THE SERVER WILL CREATE THE FOLLOWIN FIELDS DEFINED IN THIS FUNCTION 
#we want to maintain values for these fields in the database 
#CREATE  list of queriable fields for V, D, and J genes 
#ENSURE that there will always CDR1,2, and 3 lengths stored in database if field is present
def AddedSeqCollectionFields(current_fields):
	if 'CDR3.AA' in current_fields:
		if current_fields['CDR3.AA']:
			current_fields['CDR3.AA_LENGTH'] = len(current_fields['CDR3.AA'])
		else:
			current_fields['CDR3.AA_LENGTH'] = None  #it is important to set it to None so that we can use $unset to remove this field if present in database during update command
	
	if 'CDR3.NT' in current_fields:
		if current_fields['CDR3.NT']:	
			current_fields['CDR3.NT_LENGTH'] = len(current_fields['CDR3.NT'])
		else:
			current_fields['CDR3.NT_LENGTH'] = None
				
	if 'VREGION.CDR1.AA' in current_fields:
		if current_fields['VREGION.CDR1.AA']:
			current_fields['VREGION.CDR1.AA_LENGTH'] = len(current_fields['VREGION.CDR1.AA'])
		else:
			current_fields['VREGION.CDR1.AA_LENGTH'] = None
						
	if 'VREGION.CDR1.NT' in current_fields:
		if current_fields['VREGION.CDR1.NT']:
			current_fields['VREGION.CDR1.NT_LENGTH'] =  len(current_fields['VREGION.CDR1.NT'])
		else:
			current_fields['VREGION.CDR1.NT_LENGTH'] = None
			
	if 'VREGION.CDR2.AA' in current_fields:
		if current_fields['VREGION.CDR2.AA']:
			current_fields['VREGION.CDR2.AA_LENGTH'] = len(current_fields['VREGION.CDR2.AA'])
		else:
			current_fields['VREGION.CDR2.AA_LENGTH'] = None
			
	if 'VREGION.CDR2.NT' in current_fields:
		if current_fields['VREGION.CDR2.NT']:
			current_fields['VREGION.CDR2.NT_LENGTH'] =  len(current_fields['VREGION.CDR2.NT'])
		else:
			current_fields['VREGION.CDR2.NT_LENGTH'] = None
	
		
	#for create_modified_fields in modified_fields_for_queries
	#this dictionary will be added to a new subdocument in seq collection called "QUERY_DATA"
	#this includes creating a novel gene list that includes boht alleles and genenames
	#create_fields_for_query = {field:modified_fields_for_queries[field](value) for field,value in current_fields.iteritems() if field in modified_fields_for_queries}
	
	#THE FOLLOWING FIELDS WILL BE ADDED TO QUERY_DATA IN THE DOCUMENT 
	#GENE DATA WILL BE STORED IN A SEPERATE SUBDOCUMENT IN THE DOCUMENT. THIS SUBDOCUMETN WILL BE CALLED QUERY_DATA
	#THE FORMAT OF GENES STORED IN QUERY_DATA WILL BE AS SUCH: 
	#QUERY_DATA = {'VREGION.VGENES':[
		#					{'PARSED_ALLELES': [IGHV5, IGHV5*01, IGHV5*02] },
		#					{'PARSED_ALLELES': [IGHV5,IGHV5*03,IGHV6,IGHV6*02 ] },
		#			]
		#		}
	#BASICALLY FOR EACH GENE PROVIDED, WE WILL CREATE A SPECIAL SUBDOCUMETN FOR IMPROVED QUERYING. YOU CAN QUERY THE TOP GENES BY EITHER ALLELE NAME OR GENE NAME . IN ADDTIIONAL, ALL OTHER GENES PROVIDED THAT ARE NOT 'TOP GENES' WILL BE LISTED IN 
	#THE SECOND INDEX OF THE ARRAY. AGAIN THIS INDEX WILL CONTAIN ALL VGENES AND VALLELES IN THAT INDEX
	#THIS STRUCTURE WILL BE TEH SAME FOR JREGION.JGENES AND DREGION.DGENES
	create_fields_for_query = {}
	if 'VREGION.VGENES' in current_fields and current_fields['VREGION.VGENES']:
		top_genes = []
		additional_genes = []		
		#if the annotation data contains the field 'VGENE_SCORES', then we will only isolate all genes whose score is top-ranked (equal to the first element in vGENE scores)
		if "VREGION.VGENE_SCORES" in current_fields and current_fields["VREGION.VGENE_SCORES"]:
			num_scores = len(current_fields['VREGION.VGENE_SCORES'])
			top_score = max(current_fields['VREGION.VGENE_SCORES']) #this will be the top vgene score 
			for i,gene in enumerate(current_fields['VREGION.VGENES']):#_SCORES']):								
				if i>=num_scores:
					additional_genes.append(gene)
					continue
				score = current_fields['VREGION.VGENE_SCORES'][i]
				if score==top_score:
					top_genes.append(gene)#append the gene correspodning to this score to the data 
				else:
					additional_genes.append(gene)#append any other gene to this element 
		else:
			top_genes = current_fields['VREGION.VGENES'] #if no scores were provided, then we assume that all the genes are equally weighted/do not split them into top scores 
		
		create_fields_for_query['VREGION.VGENES'] = []		
		create_fields_for_query['VREGION.VGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(top_genes)})
		if additional_genes:
			create_fields_for_query['VREGION.VGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(additional_genes)})		
	
	#REPEAT FOR D GENES
	if 'DREGION.DGENES' in current_fields and current_fields['DREGION.DGENES']:
		top_genes = []
		additional_genes = []
		#if the annotation data contains the field 'VGENE_SCORES', then we will only isolate all genes whose score is top-ranked (equal to the first element in vGENE scores)
		if "DREGION.DGENE_SCORES" in current_fields and current_fields["DREGION.DGENE_SCORES"]:
			num_scores = len(current_fields['DREGION.DGENE_SCORES'])
			top_score = max(current_fields['DREGION.DGENE_SCORES']) #this will be the top vgene score 
			for i,gene in enumerate(current_fields['DREGION.DGENES']):#_SCORES']):								
				if i>=num_scores:
					additional_genes.append(gene)
					continue
				score = current_fields['DREGION.DGENE_SCORES'][i]
				if score==top_score:
					top_genes.append(gene)#append the gene correspodning to this score to the data 
				else:
					additional_genes.append(gene)#append any other gene to this element
		else:
			top_genes = current_fields['DREGION.DGENES'] #if no scores were provided, then we assume that all the genes are equally weighted/do not split them into top scores 
		
		create_fields_for_query['DREGION.DGENES'] = []		
		create_fields_for_query['DREGION.DGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(top_genes)})
		if additional_genes:
			create_fields_for_query['DREGION.DGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(additional_genes)})		
		
	#REPEAT FOR J genes	
	if 'JREGION.JGENES' in current_fields and current_fields['JREGION.JGENES']:
		top_genes = []
		additional_genes = []
		#if the annotation data contains the field 'VGENE_SCORES', then we will only isolate all genes whose score is top-ranked (equal to the first element in vGENE scores)
		if "JREGION.JGENE_SCORES" in current_fields and current_fields["JREGION.JGENE_SCORES"]:
			num_scores = len(current_fields['JREGION.JGENE_SCORES'])
			top_score = max(current_fields['JREGION.JGENE_SCORES']) #this will be the top vgene score 
			for i,gene in enumerate(current_fields['JREGION.JGENES']):#_SCORES']):								
				if i>=num_scores:
					additional_genes.append(gene)
					continue
				score = current_fields['JREGION.JGENE_SCORES'][i]
				if score==top_score:
					top_genes.append(gene)#append the gene correspodning to this score to the data 
				else:
					additional_genes.append(gene)#append any other gene to this element
		else:
			top_genes = current_fields['JREGION.JGENES'] #if no scores were provided, then we assume that all the genes are equally weighted/do not split them into top scores 
		
		create_fields_for_query['JREGION.JGENES'] = []		
		create_fields_for_query['JREGION.JGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(top_genes)})
		if additional_genes:
			create_fields_for_query['JREGION.JGENES'].append({'PARSED_ALLELES':Allele_Name_To_List(additional_genes)})		
	
	
	
	return [current_fields,create_fields_for_query]


def indexed_fields_in_analysis():
	# The "index" key is not part of the database schema. It's included here to keep track of which fields are to be indexed.
	indexed_fields = {
		'CDR1': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.CDR1.AA':1, 'ANALYSIS.INDEXED_FIELDS.CDR1.NT':1}},
		'CDR2': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.CDR2.AA':1, 'ANALYSIS.INDEXED_FIELDS.CDR2.NT':1}},
		'CDR3': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.CDR3.AA':1, 'ANALYSIS.INDEXED_FIELDS.CDR3.NT':1}},
		'FR1': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.FR1.AA':1, 'ANALYSIS.INDEXED_FIELDS.FR1.NT':1}},
		'FR2': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.FR2.AA':1, 'ANALYSIS.INDEXED_FIELDS.FR2.NT':1}},
		'FR3': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.FR3.AA':1, 'ANALYSIS.INDEXED_FIELDS.FR3.NT':1}},
		'FR4': {'ANALYSIS': None, 'NT': None, 'AA': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.FR4.AA':1, 'ANALYSIS.INDEXED_FIELDS.FR4.NT':1}},
		'VGENES': {'ANALYSIS': None, 'GENE_NAME': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.VGENES.GENE_NAME':1}},
		'DGENES': {'ANALYSIS': None, 'GENE_NAME': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.DGENES.GENE_NAME':1}},
		'JGENES': {'ANALYSIS': None, 'GENE_NAME': None, 'RECOMB': None, 'index':{'ANALYSIS.INDEXED_FIELDS.JGENES.GENE_NAME':1}}
	}
	
	indexed_field_value_options = {
		'RECOMB': {'analysis_fields': ['VDJ', 'VJ'], 'transform': 'convert_text_to_index_field_text'},
		'AA': {'analysis_fields': ['AA'], 'transform': 'convert_text_to_index_field_text'},
		'NT': {'analysis_fields': ['NT'], 'transform': 'convert_text_to_index_field_text'},
		'GENE_NAME': {'analysis_fields': ['VGENES', 'DGENES', 'JGENES'], 'transform': 'allele_name_to_list'}
	}
	
	return indexed_fields

#if we download fasta file from mongo database, then we will append the doument identifier and other ifnromation to the seqheader
#this parses the file and reports the results
def extract_DB_INFO(seq_header):	
	info = {}
	parsed_seqheader = seq_header.split(fasta_file_delimiter);
	if len(parsed_seqheader)==1:
		info["SEQUENCE_HEADER"] = parsed_seqheader[0].replace('>','')
	elif len(parsed_seqheader)==2:		
		mongo_info = json.loads(parsed_seqheader[1].strip())		
		id_found = False 
		
		info["SEQUENCE_HEADER"] = parsed_seqheader[0].replace('>','')
		for keyval in mongo_info:
			if keyval == idIdentifier:
				if type(mongo_info[keyval]) is dict:
					info[idIdentifier] = mongo_info[keyval]["$oid"]
					id_found = True
				else:
					info[idIdentifier] = mongo_info[keyval]
					id_found = True
			else:
				info[keyval] = mongo_info[keyval]			
		if id_found == False:
			info = {}										
	else:
		info = {}		
	return info
	
