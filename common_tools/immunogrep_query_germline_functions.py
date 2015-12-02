import pymongo
import pandas as pd
import time
from pprint import pprint
import json
from operator import itemgetter
import math
import os
from immunogrep_database_schema import convert_to_objectid
from collections import defaultdict

# To signify additional info in sequence header 
fasta_file_delimiter = ' <' 


def connectToIgDatabase():
	"""
		creates a database connection
	"""
	connection_data_template = {
		'connection': None,	
		'germline': None,
		'motifs': None
	}
	connection = pymongo.MongoClient(host="geordbas01.ccbb.utexas.edu", port=27017)
	db = connection.germlines
	db.authenticate('germlinereader', 'germlinereader')	

	connection_data_template['connection'] = connection
	connection_data_template['germline'] = db.Germlines
	connection_data_template['motifs'] = db.Motifs	
	return connection_data_template

# will run a query for cdr3 motif given settings
# settings is a list of tuples. Each tuple in list consists of (Speices name, locus)
# locus may be a list 
def ReturnCDR3Motif(settings):
	germline_query = []
	for possible_combos in settings:
		if not isinstance(possible_combos[1], list):					
			germline_query.append({'Species': possible_combos[0],'Locus':possible_combos[1]}) #using list(set to ensure only unique values in list 
		else:
			germline_query.append({'Species':possible_combos[0],'Locus':{'$in':list(set(possible_combos[1]))}})
	germline_query = {'$or':germline_query}
	
	#now run the query 
	cdr3_motif_class = CDR3MotifDB()
	selectedmotif = cdr3_motif_class.GetMotifForProgram(query=germline_query)
	return selectedmotif

##class for querying##
class CDR3MotifDB():
	def __init__(self):
		db_dict = connectToIgDatabase()
		self.db = db_dict['motifs']		
	def formatMotif(self,motiflist):
		outputlist=[]
		for i in motiflist:
			lmotif=defaultdict(dict)
			rmotif=defaultdict(dict)
			for j in i[1]:
				tmp=j.split('_')
				lmotif[tmp[1]][tmp[0]]=float(tmp[2])	
			for j in i[2]:
				tmp=j.split('_')
				rmotif[tmp[1]][tmp[0]]=float(tmp[2])
			chain=str(i[0])
			ltrim = int(i[3])
			rtrim = int(i[4])
			settemp=(chain,lmotif,rmotif,ltrim,rtrim)
			outputlist.append(settemp)
		return outputlist
	
	def GetAllDocs(self):
		cursor = self.db.find()
		results = []
		for i in cursor:			
			results.append(i)
		return results
	def GetUniqSpecies(self):
		return self.db.distinct("Species")
	def GetMotifForProgram(self,id_list=[],query={},project={'Locus':1,'Lmotif':1,'Rmotif':1,'Ltrim':1,'Rtrim':1}):
		if id_list:
			id_list = [convert_to_objectid(id) for id in id_list]
			query['_id'] = {'$in':id_list}
		cursor = self.db.find(query,project)
		#return results as tuple whose indices are as follows: ['Locus','Lmotif','Rmotif','Ltrim','Rtrim']]	
		return self.formatMotif([(x['Locus'],x['Lmotif'],x['Rmotif'],x['Ltrim'],x['Rtrim']) for x in cursor])

class GermlineDB:
	def __init__(self):		
		db_dict = connectToIgDatabase()
		self.db = db_dict['germline']
		self.cursor = None
		self.results = None
		
	def GetAllDocs(self):
		self.cursor = self.db.find()
		self.results = []
		for i in self.cursor:			
			self.results.append(i)
		return self.results
		
	def GetAllDocInfo(self,id_list=[],extra_filters = {}):#return top-level results of all documents but leave out information regarding genes within a document				
		if id_list != []:
			id_list = [convert_to_objectid(id) for id in id_list]
			extra_filters['_id'] = {'$in': id_list}
		self.cursor = self.db.find(extra_filters,{'GERMLINE_GENES':0})
		self.results = []
		for i in self.cursor:			
			self.results.append(i)
		return self.results
		
	def QueryAllDocInfoByID(self,id_list=[],extra_filters={},gene_functionality_filter={'$nin':[]}):
		if id_list:
			id_list = [convert_to_objectid(id) for id in id_list]
			extra_filters['_id'] = {'$in':id_list}		
		
		if type(gene_functionality_filter) != 'dict':		
			raise Exception('Gene functionality filter must be a dict')
		
		self.cursor = self.db.aggregate([
			{'$match':extra_filters},
			{'$unwind':"$GERMLINE_GENES"},
			{'$match': {"GERMLINE_GENES.FUNCTIONALITY": gene_functionality_filter} }			
		])
		self.results = [r.pop('GERMLINE_GENES') for r in self.cursor]  #self.cursor['result']
		#for i,r in enumerate(self.results):
		#	self.results[i].pop('GERMLINE_GENES')			
		return self
	
	def QueryDistinctValsByID(self, id_list=[], extra_filters={}, distinct_fields=[]):
		if id_list:
			id_list = [convert_to_objectid(id) for id in id_list]
			extra_filters['_id'] = {'$in':id_list}		
		results = {field: self.db.find(extra_filters).distinct(field) for field in distinct_fields}
		return results
		
	
	#function will return all information for all documents that match the input list of IDs
	def QueryAllDocsByID(self,id_list=[],extra_filters={},gene_functionality_filter={'$nin':[]}):
		if id_list:
			id_list = [convert_to_objectid(id) for id in id_list]
			extra_filters['_id'] = {'$in':id_list}				
		if type(gene_functionality_filter) != 'dict':		
			raise Exception('Gene functionality filter must be a dict')
	
		self.cursor = self.db.aggregate([
			{'$match':extra_filters},
			{'$unwind':"$GERMLINE_GENES"},
			{'$match': {"GERMLINE_GENES.FUNCTIONALITY": gene_functionality_filter} },
			{'$sort':{"GERMLINE_GENES.ALLELE_NAME":1}}
		])
		return self
			
	#function for query germline genes based on 1) list of ids, 2) any additional filters from top level keys in documents, 3) gene functionality
	#function will return results as a dictionary where each key in dictionary is either 'V', 'D', or 'J'. each value of the key will be an array of all genes returned by query
	#include_in_group = > a list of additional fields you want to "join" with the germline gene results (i.e. fields from top level results such as 'CHAIN' and 'SPECIES'
	def QueryGenesByID(self,id_list=[],extra_filters={},gene_functionality_filter={'$nin':[]},include_in_group = ["CHAIN","SPECIES","LOCUS"]):		
		if '$in' in gene_functionality_filter and gene_functionality_filter['$in']==[]:
			gene_functionality_filter = {'$nin':[]}
		
		if id_list:
			id_list = [convert_to_objectid(id) for id in id_list]
			extra_filters['_id'] = {'$in':id_list}					
		
		include_in_group = {field:'$'+field for field in include_in_group}			
		if type(gene_functionality_filter) != dict:		
			raise Exception('Gene functionality filter must be a dict')						
				
		#self.cursor = self.db.find({'_id':{'$in':id_list}})				
		
		
		self.cursor = self.db.aggregate([
					{'$match':extra_filters}, #search for documents which match variables passed in this dictionary	
					{'$unwind':"$GERMLINE_GENES"}, #unwind results such that each element is now a unique germline gene (germline genes function was originally list of all genes in set)
					{'$match':{"GERMLINE_GENES.FUNCTIONALITY": gene_functionality_filter } }, #filter genes by their productivity 					
					{'$sort':{"GERMLINE_GENES.ALLELE_NAME":1}}, #sort results by allele name
					{'$group':{'_id':"$GENETYPE",'genes':{'$push':"$GERMLINE_GENES"},'additional_info':{'$push':include_in_group } } } #and finally group together all genes by their gene type (V,D,AND J)
				])				
		#self.results = {subset['_id']: [dict(gene_dict.items()+subset['additional_info'][row].items()) for row,gene_dict in enumerate(subset['genes'])] for subset in self.cursor['result']} #summarize teh results as a dictionary where each key corresponds to genetype(V,D,J,ETC) and each value is a list of dictionaries showing document results, it merges results from "genes" and "additional_info" dictionaires into one								
		self.results = {subset['_id']: [dict(gene_dict.items()+subset['additional_info'][row].items()) for row,gene_dict in enumerate(subset['genes'])] for subset in self.cursor} #summarize teh results as a dictionary where each key corresponds to genetype(V,D,J,ETC) and each value is a list of dictionaries showing document results, it merges results from "genes" and "additional_info" dictionaires into one								
		
		return self
		
	def ReturnData(): #reutnr the result as a list of dictionaries
		if self.cursor:
			self.results = []
			for i in self.cursor:
				self.results.append(i)
		return self.results
	
	def ReturnRawData(): #dont do any fancy processing, just return the results
		return self.results
	
	def ReturnRawCursor(): #dont do any fancy processing, just return the cursor
		return self.cursor
	
	def PrintToFasta(self,filename,list_of_values,keys_to_print=['ALLELE_NAME','SEQUENCE']): #assumes a list of dictionaries							
		with open(filename,'w') as outfile:
			for gene_of_interest in list_of_values:
				new_dict = {key:str(gene_of_interest[key]) for key in keys_to_print if key in gene_of_interest } if keys_to_print else {key:str(gene_of_interest[key]) for key in gene_of_interest}
				sequence = new_dict['SEQUENCE']
				header = new_dict['ALLELE_NAME']
				new_dict.pop('SEQUENCE')
				new_dict.pop('ALLELE_NAME')				
				sub_header = fasta_file_delimiter+json.dumps(new_dict,sort_keys=True) if len(new_dict)>0 else ''
				outfile.write('>{0}{1}\n'.format(header,sub_header))
				outfile.write('{0}\n'.format(sequence))
				
	def PrintToTAB(self,filename,list_of_values,keys_to_print=['ALLELE_NAME','SEQUENCE']):		
		with open(filename,'w') as outfile:			
			outfile.write('\t'.join(keys_to_print)+'\n')
			for gene_of_interest in list_of_values:								
				vals_to_save = []
				for key in keys_to_print:
					if key in gene_of_interest:
						if type(gene_of_interest[key]) is float and math.isnan(gene_of_interest[key]): #this shoudl catch edge cases where the values in the database are NAN and we have nto fixed those erros yet
							vals_to_save.append('')
							print("Please notify administrator that the germline database contains values that are NAN")														
						else:
							vals_to_save.append(str(gene_of_interest[key]))						
						
				#str_result = '\t'.join([gene_of_interest[key] if key in gene_of_interest else '' for key in keys_to_print])																				
				str_result = '\t'.join(vals_to_save)
				outfile.write(str_result+'\n')
		
	def PrintIgBlastDBFormat(self,parent_folder="scratch/",filename='Germline'):
		filename = filename.split('.')
		filename = '.'.join(filename[:-1]) if len(filename)>1 else filename[0]
		parent_folder = parent_folder + '/' if parent_folder[-1] != '/' else parent_folder
		suffix = '_' if filename != '' else ''	
		
		for germlines in self.results:
			genes_to_write = self.results[germlines]
			#genes_to_write = sorted(self.results[germlines], key=itemgetter('ALLELE_NAME')) #sort list of dictionary
			self.PrintToFasta(parent_folder + filename + suffix + germlines + '.txt', genes_to_write, ['ALLELE_NAME', 'SEQUENCE'])
	
	def PrintFFTDBFormat(self, parent_folder=None, filename='Germline'):
		if parent_folder is None:
			os.path.dirname('.')
		
		output_path = os.path.join(parent_folder, filename)
		fields_to_remove = ['GENENAME', 'SEQUENCE_WITH_IMGT_GAPS', 'GENENAME_INDEX']				
		for germlines in self.results:	
			unique_fields = ['GENE', 'SEQUENCE', 'LOCUS', 'FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'CHAIN', 'FUNCTIONALITY', 'SPECIES'] if germlines == 'V' else ['GENE','SEQUENCE','LOCUS','FUNCTIONALITY','CHAIN','SPECIES']
			genes_to_write = []
			for gene in self.results[germlines]:
				for f in fields_to_remove:
					if f in gene:
						gene.pop(f)
				gene['GENE'] = gene.pop('ALLELE_NAME')								
				genes_to_write.append(gene)			
			# genes_to_write = sorted(genes_to_write, key=itemgetter('GENE')) #sort list of dictionary			
			self.PrintToTAB(output_path, genes_to_write, unique_fields)
			
	def UniqueSpecies(self):
		results = self.db.distinct('SPECIES')
		return results
	
	def UniqueLoci(self):
		results = self.db.distinct('LOCUS')
		return results
	
	def DBSources(self):
		results = self.db.distinct('SOURCE')		
		return results
	
	def MolecularComponent(self):		
		results = self.db.distinct('MOLECULAR_COMPONENT')
		return results
	
	def QueryDoc(self):
		pass
	
###end of class###

