import immunogrep_makeappfunctions as makeapp
import appsoma_api
import pymongo
import pandas as pd
import time
from pprint import pprint

#-----------------------------------------------------------------------------------------------------------------------------------------
# MongoDB connection
#-----------------------------------------------------------------------------------------------------------------------------------------
# opens the database connection.
def connectToIgDatabase():
	connection_data_template = {
	'connection': None,	
	'germline':None,
	'motifs':None
	}
	
	connection = pymongo.Connection ( host="biotseq.icmb.utexas.edu", port=27017 )
	db = connection.appsoma
	db.authenticate('reader','cdrom')
	
	connection_data_template['connection'] = connection
	connection_data_template['germline'] = db.Germlines
	connection_data_template['motifs'] = db.Motifs
	
	return connection_data_template
###end connection

##class for querying##
class CDR3MotifDB:
	def __init__(self):
		db_dict = connectToIgDatabase()
		self.db = db_dict['motifs']		
	def GetAllDocs(self):
		cursor = self.db.find()
		results = []
		for i in cursor:			
			results.append(i)
		return results

class GermlineDB:
	def __init__(self):
		db_dict = connectToIgDatabase()
		self.db = db_dcit['germline']
	
	def GetAllDocs(self):
		cursor = self.db.find()
		results = []
		for i in cursor:			
			results.append(i)
		return results
		
	def GetAllDocInfo(self):#return results of all documents but leave out information regarding sequences
		cursor = self.db.find({'GERMLINE_GENES':0})
	
	def UniqueSpecies(self);
		pass
	
	def UniqueLoci(self):
		pass
	
	def DBSources(self):
		pass
	
	def MolecularComponent(self):
		pass
	
	def QueryDoc(self):
		pass
	
###end of class###

