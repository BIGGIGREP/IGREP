#~~~~Updated on 6/16/2015 to support reading IMGT files (previous version is 213)~~~#

#~~~~Updated on 5/6/2015 to allow for field names to be custom passed into initialization (previous version is 207)~~~#


#~~~~beta version updated on 06/06/2014~~~~#
#~~~~beta version where FASTA and FASTQ reading changed updated on 3/19/2015 (previous beta version is 186) ~~~~~
import sys
import os
import json
import re

try:
	from Bio import SeqIO
except:
	print "BIO NOT INSTALLED!!!READING FASTA AND FASTQ FILES WILL BE AFFECTED"
	
from collections import defaultdict


#import immunogrep_useful_functions as useful#import removeNoneVals
#from immunogrep_useful_functions import fieldsForAnnotatingAb
#from immunogrep_useful_functions import removeFileExtension

# +++++++++++++++++++++++++++++
# +Importing Global variables +
# +++++++++++++++++++++++++++++ 
from immunogrep_global_variables import descriptor_symbol
from immunogrep_global_variables import fasta_file_delimiter
from immunogrep_global_variables import description_var
from immunogrep_global_variables import translation_var
from immunogrep_global_variables import annotationExt
from immunogrep_global_variables import analysisExt
from immunogrep_global_variables import queryExt
from immunogrep_global_variables import IFExt
from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import listofextension

from immunogrep_global_variables import imgt_files_info




# +++++++++++++++++++
# + File variable   +
# +++++++++++++++++++
current_filetypes={
	'FASTA':['fna','fasta'],
	'TAB':['txt','tab'],
	'CSV':['csv'],
	'FASTQ':['fastq'],
	'JSON':['json','annotation','analysis','query','iffile','IFFile'],
	'DELIM': ['delim'],	
	'IMGT':['imgt']
}

#LAMBDA FUNCTION FOR RECOMBINATION FIELD IS A HACK TO ENSURE THAT ONLY FILE ONE IS REALLY REQUIRED 
#FILE 1 WILL ALWAYS BE REQUIRED 
imgt_db_translator = {
	"ANALYSIS_NAME":'IMGT',
	"RECOMBINATION_FIELD":{		
		"LAMBDA_FUNCTION":"lambda x: 'VDJ' if x['V-D-J-REGION_3'] else ('VJ' if x['V-J-REGION_3'] else ('' if not(x['V-GENE and allele_1']) else ('VDJ' if ('IGHV' in x['V-GENE and allele_1'] or 'TRBV' in x['V-GENE and allele_1']) else 'VJ')))",
	},
	"FIELDS_UPDATE":{
		idIdentifier:idIdentifier,
		'COMMAND':'COMMAND',
		'SEQUENCE_HEADER':'Sequence ID_1',		
		'SEQUENCE':'Sequence_1',
		'STRAND':'Orientation_1',
		'PRODUCTIVE':'Functionality_1',
		'VREGION.VGENE_SCORES':'V-REGION score_1',		
		'JREGION.JGENE_SCORES':'J-REGION score_1',		
		'NOTES':'Functionality comment_1',
		
		'VREGION.CDR1.NT':'CDR1-IMGT_3',
		'GAPPED.VREGION.CDR1.NT':'CDR1-IMGT_2',
		
		'VREGION.CDR1.AA':'CDR1-IMGT_5',
		'GAPPED.VREGION.CDR1.AA':'CDR1-IMGT_4',
		
		'VREGION.CDR2.NT':'CDR2-IMGT_3',
		'GAPPED.VREGION.CDR2.NT':'CDR2-IMGT_2',
		
		'VREGION.CDR2.AA':'CDR2-IMGT_5',
		'GAPPED.VREGION.CDR2.AA':'CDR2-IMGT_4',
		
		'CDR3.NT':'CDR3-IMGT_3',		
		'GAPPED.CDR3.NT':'CDR3-IMGT_2',				
		
		'CDR3.AA':'CDR3-IMGT_5',
		'GAPPED.CDR3.AA':'CDR3-IMGT_4',		
		
		'VREGION.FR1.NT':'FR1-IMGT_3',
		'GAPPED.VREGION.FR1.NT':'FR1-IMGT_2',
		
		'VREGION.FR1.AA':'FR1-IMGT_5',
		'GAPPED.VREGION.FR1.AA':'FR1-IMGT_4',
		
		'VREGION.FR2.NT':'FR2-IMGT_3',
		'GAPPED.VREGION.FR2.NT':'FR2-IMGT_2',
		
		'VREGION.FR2.AA':'FR2-IMGT_5',
		'GAPPED.VREGION.FR2.AA':'FR2-IMGT_4',
				
		'VREGION.FR3.NT':'FR3-IMGT_3',
		'GAPPED.VREGION.FR3.NT':'FR3-IMGT_2',
		
		'VREGION.FR3.AA':'FR3-IMGT_5',	
		'GAPPED.VREGION.FR3.AA':'FR3-IMGT_4',
		
		'JREGION.FR4.NT':'FR4-IMGT_3',		
		'GAPPED.JREGION.FR4.NT':'FR4-IMGT_2',
		
		'JREGION.FR4.AA':'FR4-IMGT_5',								
		'GAPPED.JREGION.FR4.AA':'FR4-IMGT_4',								
		
		'VREGION.VGENE_QUERY_START':'V-REGION start_3',
		'VREGION.VGENE_QUERY_END':'V-REGION end_3',
		'JREGION.JGENE_QUERY_START':'J-REGION start_3',
		'JREGION.JGENE_QUERY_END':'J-REGION end_3',
		'VREGION.VGENES':'V-GENE and allele_1',
		'JREGION.JGENES':'J-GENE and allele_1',
		'DREGION.DGENES':'D-GENE and allele_1',
		
		'VREGION.SHM.NT':'VREGION.SHM.NT', #created in filereader class (uses the field from fithe file summary V_REGION IDENTITY)
		'VREGION.SHM.NT_PER':'VREGION.SHM.NT_PER', #created in filereader class 
		'VREGION.SHM.AA':'VREGION.SHM.AA', #created in filereader class
		'VREGION.SHM.AA_PER':'VREGION.SHM.AA_PER', #created in filereader class
		
		'JREGION.SHM.NT':'JREGION.SHM.NT', #created in filereader class (uses the field from fithe file summary J_REGION IDENTITY)
		'JREGION.SHM.NT_PER':'JREGION.SHM.NT_PER', #created in filereader class 
		
		"PREDICTED_AB_SEQ.NT":"PREDICTED_AB_SEQ.NT",#created in filereader class
		"PREDICTED_AB_SEQ.AA":"PREDICTED_AB_SEQ.AA"#created in filereader class
				
	}
}


# ++++++++++++++++++++++++++++++++
# ++++Default Decorator values++++
# ++++++++++++++++++++++++++++++++
FASTAdecoratorinfo={description_var:['header','sequence']}
FASTQdecoratorinfo={description_var:['header','sequence','phred']}

# ++++++++++++++++++++++++++++++++++++
# ++++for extract_header function ++++
# ++++++++++++++++++++++++++++++++++++
num_line_search = 10

#sometimes additional info such as sequence id is passed along in the sequence header
#this fucntion will check if any immunogrep data is in sequence header and return it
def GetAdditionalInfo(header_line):
	additional_info = {'document_header':header_line}
	if fasta_file_delimiter in header_line:
		tmp=header_line.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter
		header=tmp[0]									
		try:
			additional_info = dict(additional_info.items()+json.loads(tmp[1]).items())
		except:
			#occurs if json string in sequence header was truncated
			pass		
	else:
		header=header_line		
	
	return [header,additional_info]

def find_imgt_file_type_index(filename,important_headers_only=True):
		filename = os.path.basename(filename)
		imgt_file_type = filename.strip().split('_')
		imgt_file_type = '_'.join(imgt_file_type[0:2])	
		imgt_file = [i+1 for i,a in enumerate(imgt_files_info) if a['Name']==imgt_file_type]
		if len(imgt_file)>0:
			return [imgt_file[0],imgt_files_info[imgt_file[0]-1]['Important_Headers']] if important_headers_only else [imgt_file[0],imgt_files_info[imgt_file[0]-1]['Headers']]
		else:
			return [None,None]
	

#WE ASSUME IMGT FILES ARE LABELED USING IMGT RULES:
#=>FILENUMBER_IMGTFILENAME_USERFILENAME.txt
def GroupIMGTFiles(imgt_files):
	#imgt_files = sorted([str(im_f) for im_f in imgt_files],key=str.lower)	
	grouped_files = defaultdict(list)
	for each_file in imgt_files:
		filenum=find_imgt_file_type_index(each_file)
		if each_file.lower().endswith('.txt') and len(filenum)>0 and filenum[0]!=None:			
			#take basename of file , split by '_'
			split_file = os.path.basename(each_file).split('_')
			#only consider the filename after 1_imgtfilename_submittedfilename.txt
			if len(split_file)>=2:
				submitted_filename = '_'.join(split_file[2:])
			else:
				submitted_filename = ''
			grouped_files[submitted_filename].append(each_file)	
	#grouped_files = {group:files for group,files in grouped_files.iteritems() if files}
	
	return grouped_files



#++++++++++++++++++++++++++++++++++++
#Simple function for reading a block of
#data from an igblast file output##
#++++++++++++++++++
def ReadIgBlastQueryBlock(f1,igblast_eof = False):

	igblast_query_block = ""
	igblast_query_seq_header = ""
	igblast_read_line = True
	
	while (igblast_read_line and not(igblast_eof)): 
		line_info = f1.readline()		
		
		if not(line_info):
			igblast_eof = True
			break
		
		if line_info.strip() == "": #skip empty lines 
			continue
					
		# IGBLASTN
		if line_info[:9].upper()== "# IGBLAST": #start of a new query 
			igblast_read_line = False
			break
		elif line_info[:17].upper() =="# BLAST PROCESSED": #at the complete end of the igblast files				
			igblast_eof = True
			break
		else:
			igblast_query_block += line_info
			
		
		if line_info[:8].upper() == "# QUERY:":
			igblast_query_seq_header = line_info[9:].strip()				
			
	####NOW WE SHOULD HAVE READ IN AN ENTIRE BLOCK FOR THE CURRENT IGBLAST query
	return [igblast_query_block,igblast_query_seq_header]



# +++++++++++++++++++++++
# +immunogrepFile Class +
# +++++++++++++++++++++++ 
class immunogrepFile():
	#the following parameters are ONLY used for IMGT : important_headers_only =True,include_parameters_file=False,required_files=0
	#the following parameters are ONLY used for DELIM files/TAB/CSV: delimiter=None, contains_header=True
	def __init__(self,filelocation,filetype=None,decoratorinfo=None,delimiter=None,contains_header=True,mode='r',field_names=[],chunk_size=1,important_headers_only =True,include_parameters_file=False,required_files=0):
		 
		self.filelocation=filelocation
		self.decoratorinfo=decoratorinfo
		self.delimiter=delimiter
		self.contains_header=contains_header
		self.mode=mode
		self.field_names=field_names
		self.chunk_size=chunk_size
		
		if filetype:
			if filetype.lower() in current_filetypes['FASTA']:
				self.filetype='FASTA'
			elif filetype.lower() in current_filetypes['TAB']:			
				self.filetype='TAB'
				self.delimiter="\t"
			elif filetype.lower() in current_filetypes['CSV']:
				self.filetype='CSV'
				self.delimiter=','
			elif filetype.lower() in current_filetypes['FASTQ']:
				self.filetype='FASTQ'
			elif filetype.lower() in current_filetypes['JSON']:
				self.filetype='JSON'
			elif filetype.lower() in current_filetypes['DELIM']:
				if delimiter:
					self.filetype='DELIM'
					self.delimiter=delimiter
				else:
					sys.exit(" 'delim' filetype requires user-defined delimiter !")
			elif filetype.lower() in current_filetypes['IMGT']:
				self.filetype='IMGT'
			#else:
			#	raise Exception('The following filetype, {0}, is invalid'.format(filetype))
			else:
				self.filetype=self.guessFiletype()
		else:			
			self.filetype=self.guessFiletype()
	
		# Checking file exisitence, if not, it will be createdf
		#IMGT FILES MUST BE LIST OF FILES, ALL OTHER FILES MUST BE BASESTRING
		if mode == 'r':
			if filetype=='IMGT':
				if not isinstance(self.filelocation,list):
					self.filelocation = [self.filelocation]
					for each_file in self.filelocation:
						if not os.path.isfile(each_file):
							raise Exception("The File {0} does not exist".format(each_file))
			else:
				if not isinstance(self.filelocation,basestring):
					raise Exception('Only immunogrep files whose filetypes are defined as "IMGT" may be passed in as a list of files. all other file types must be strings')
				if not os.path.isfile(self.filelocation):							
					raise Exception("The file {0} does not exist".format(self.filelocation))
					#sys.exit("The file does not exist")
			

		if self.filetype:
			self.filetypedict={'IMGT':immunogrepIMGT,'FASTA':immunogrepFASTA,'FASTQ':immunogrepFASTQ,'TAB':immunogrepDELIM,'CSV':immunogrepDELIM,'JSON':immunogrepJSON,'DELIM':immunogrepDELIM}						
			
			if self.filetype == 'DELIM' or self.filetype=='CSV' or self.filetype=='TAB':
				self.IFclass=self.filetypedict[self.filetype](filelocation=self.filelocation,filetype=self.filetype,decoratorinfo=self.decoratorinfo,delimiter=self.delimiter,contains_header=self.contains_header,mode=self.mode,field_names=self.field_names,chunk_size=self.chunk_size)			
			elif self.filetype=='IMGT':			
				self.IFclass=immunogrepIMGT(imgtfiles=self.filelocation,field_names = field_names, chunk_size = chunk_size,important_headers_only=important_headers_only,include_parameters_file=include_parameters_file,required_files=required_files)
			else:
				self.IFclass=self.filetypedict[self.filetype](filelocation=self.filelocation,filetype=self.filetype,decoratorinfo=self.decoratorinfo,mode=self.mode,field_names=self.field_names,chunk_size=self.chunk_size)
		else:			
			self.IFclass=None			
		
	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return self.filetype

	def guessFiletype(self):
		check = []
		fileformat=['FASTA','FASTQ','JSON','TAB','CSV']
	
		check.append(self.isFASTA())
	
		check.append(self.isFASTQ())
	
		check.append(self.isJSON())

		check.append(self.isTAB())

		check.append(self.isCSV())
		
		checkcnt=check.count(True)
		
		if check[2] == True:			
			#if JSON is listed as true, then irespect to other checks and report as JSON 
			self.filetype=fileformat[2]			
		elif checkcnt==1:
			self.filetype=fileformat[check.index(True)]
			if self.filetype=='TAB':
				self.delimiter="\t"
			if self.filetype=='CSV':
				self.delimiter=","					
		else:			
			self.filetype=None
		
		return self.filetype

	def getDecoratorinfo(self):
		if self.isimmunogrepFile() and self.decoratorinfo==None:
			with open(self.filelocation,"r") as f:
				line=f.readline()
				if descriptor_symbol in line:
					try:
						self.decoratorinfo=json.loads(line[len(descriptor_symbol):].strip().replace("'",'"'))
					except:
						self.decoratorinfo={}
			return self.decoratorinfo
		else:
			return self.decoratorinfo

	def getDescription(self):
		if self.IFclass:
			if self.getDecoratorinfo()!=None:
				if self.getDecoratorinfo()[description_var]:
					d=self.getDecoratorinfo()[description_var]
				else:				
					d=self.IFclass.Extract_Header() 		
			else:			
				d=self.IFclass.Extract_Header()
			return d
		else:		
			return None
		
	def getTranslation(self):
		if self.getDecoratorinfo()!=None:
			if self.getDecoratorinfo()[translation_var]:
				t=self.getDecoratorinfo()[translation_var]
			else:
				t=None
		else:
			t=None
		return t

	# Note: create method will create/overwrite file specificed in filelocation and add decorator to the first line
	# This function is a residual function; might be obsolete
	def create(self):
		with open(self.filelocation+'.'+IFExt,"w") as f:
			f.write(descriptor_symbol+str(self.decoratorinfo)+'\n')
		self.filelocation=self.filelocation+'.'+IFExt
		self.filetype='JSON'

	def isimmunogrepFile(self):
		check=False
		with open(self.filelocation,"r") as f:
			line=f.readline().strip()
			if descriptor_symbol in line:
				try:
					json.loads(line[len(descriptor_symbol):].replace("'",'"'))
					check=True
				except:
					check=False
			else:
				check=False
			line=f.readline().strip()
			try:
				json.loads(line.replace("'",'"'))
				check=True
			except:
				check=False
		return check and IFExt in self.filelocation

	def isFASTA(self):
		check=False
		try:
			preview_lines = 100
			temp_file = self.filelocation+'__temp__'
			bash_command = 'head -{0} "{1}" > "{2}" '.format(preview_lines,self.filelocation,temp_file)
			os.system(bash_command)
			with open(temp_file,"r") as f:			
				line=f.readline()
				
				if line.startswith(descriptor_symbol):
					line=f.readline()				
				else:				
					f.seek(0,0)			
					
				g=SeqIO.parse(f,"fasta")
									
				try:
					g.next()
					check=True
				except:	
					check=False
			os.system('rm "{0}"'.format(temp_file))
			return check
		except:
			return False
		
	def isFASTQ(self):
		check=False
		try:
			preview_lines = 100
			temp_file = self.filelocation+'__temp__'
			bash_command = 'head -{0} "{1}" > "{2}" '.format(preview_lines,self.filelocation,temp_file)
			os.system(bash_command)
			with open(temp_file,"r") as f:
				line=f.readline()
				if line.startswith(descriptor_symbol):
					line=f.readline()
				else:
					f.seek(0,0)
				g=SeqIO.parse(f,"fastq")
				try:
					g.next()
					check=True
				except:
					check=False
			os.system('rm "{0}"'.format(temp_file))
			return check
		except:
			return False

	# Note: Must be a block delimited file format; jaggered delimited file format is considered invalid
	def isDelim(self,delimiter):
		check=False
		chkcnt=[]
		
		try:
			with open(self.filelocation,"r") as f:
				line=f.readline().strip('\r\n')
				if descriptor_symbol in line:
					line=f.readline().strip('\r\n')
				else:
					f.seek(0,0)
				for i in range(num_line_search):
					if line=='':						
						break									
					chkcnt.append(line.count(delimiter))
					line=f.readline().strip('\r\n')
	
			unique_chckcnt = set(chkcnt)
			if len(unique_chckcnt) == 1 and list(unique_chckcnt)[0]!=0:
				check=True
			else:
				check=False		
			return check
		except:
			return False

	def isTAB(self):
		check=self.isDelim("\t")
		return check
	
	def isCSV(self):
		check=self.isDelim(',')
		return check

	def isJSON(self):
		check=True
		chkcnt=0
		try:
			with open(self.filelocation,"r") as f:
				line=f.readline().strip()
				if descriptor_symbol in line:
					line=f.readline().strip()
				else:
					f.seek(0,0)
				for i in range(num_line_search):
					line=f.readline().strip()
					if line=='':
						break
					try:
						json.loads(line)
					except:
						check=False				
			return check
		except:
			return False

	def convert2FASTA(self,headervar,sequencevar,keepinfo=False):
		if self.filetype=='JSON':
			with open(self.filelocation+'.fna',"w") as fout:
				with open(self.filelocation,"r") as f:
					for line in f:
						line=line.strip()
						if descriptor_symbol not in line:
							body=json.loads(line)
							addbody=dict()
							if keepinfo:
								for k in body:
									if k not in [headervar,sequencevar]:
										addbody[k]=body[k]
								fout.write(">%s\n%s\n"%(body[headervar]+fasta_file_delimiter+json.dumps(addbody),body[sequencevar]))
							else:
								fout.write(">%s\n%s\n"%(body[headervar],body[sequencevar]))
						else:
							fout.write(line+"\n")
		self.filelocation=self.filelocation+'.fna'
		self.filetype='FASTA'

	def convert2FASTQ(self,headervar,sequencevar,qualvar,keepinfo=False):
		if self.filetype=='JSON':
			with open(self.filelocation+'.fastq',"w") as fout:
				with open(self.filelocation,"r") as f:
					for line in f:
						line=line.strip()
						if descriptor_symbol not in line:
							body=json.loads(line)
							addbody=dict()
							if keepinfo:
								for k in body:
									if k not in [headervar,sequencevar,qualvar]:
										addbody[k]=body[k]
								fout.write("@%s\n%s\n+\n%s\n"%(body[headervar]+fasta_file_delimiter+json.dumps(addbody),body[sequencevar],body[qualvar]))
							else:
								fout.write("@%s\n%s\n+\n%s\n"%(body[headervar],body[sequencevar],body[qualvar]))
						else:
							fout.write(line+"\n")
		self.filelocation=self.filelocation+'.fastq'
		self.filetype='FASTQ'

	def convert2Delim(self,descriptionlist,delimiter):
		if self.filetype=='JSON':
			with open(self.filelocation+'.txt',"w") as fout:
				with open(self.filelocation,"r") as f:
					for i in descriptionlist:
						print >>fout,"%s\t"%i,
					print >>fout,""
					for line in f:
						line=line.strip()
						if descriptor_symbol not in line:
							body=json.loads(line)
							bodytxt=None
							for i in descriptionlist:
								if body[i]==None:
									body[i]='N/A'
								if not bodytxt:
									bodytxt=str(body[i])
								else:
									bodytxt=bodytxt+delimiter+str(body[i])
							fout.write("%s\n"%bodytxt)
						else:
							fout.write(line+"\n")
							label=None
							for i in descriptionlist:
								if not label:
									label=i
								else:
									label=label+delimiter+i
							fout.write("%s\n"%label)
		self.filelocation=self.filelocation+'.txt'
		self.filetype='DELIM'	

	def convert2TAB(self,descriptionlist):
		if self.filetype=='JSON':
			with open(self.filelocation+'.txt',"w") as fout:
				with open(self.filelocation,"r") as f:
					for i in descriptionlist:
						print >>fout,"%s\t"%i,
					print >>fout,""
					for line in f:
						line=line.strip()
						if descriptor_symbol not in line:
							body=json.loads(line)
							bodytxt=None
							for i in descriptionlist:
								if body[i] == None:
									body[i]='N/A'
								if not bodytxt:
									bodytxt=str(body[i])
								else:
									bodytxt=bodytxt+"\t"+str(body[i])
							fout.write("%s\n"%bodytxt)
						else:
							fout.write(line+"\n")
							label=None
							for i in descriptionlist:
								if not label:
									label=i
								else:
									label=label+"\t"+i
							fout.write("%s\n"%label)
		self.filelocation=self.filelocation+'.txt'
		self.filetype='TAB'
		
	def convert2CSV(self,descriptionlist):
		if self.filetype=='JSON':
			with open(self.filelocation+'.csv',"w") as fout:
				with open(self.filelocation,"r") as f:
					for i in descriptionlist:
						print >>fout,"%s\t"%i,
					print >>fout,""
					for line in f:
						line=line.strip()
						if descriptor_symbol not in line:
							body=json.loads(line)
							bodytxt=None
							for i in descriptionlist:
								if body[i]==None:
									body[i]='N/A'
								if not bodytxt:
									bodytxt=str(body[i])
								else:
									bodytxt=bodytxt+","+str(body[i])
							fout.write("%s\n"%bodytxt)
						else:
							fout.write(line+"\n")
							label=None
							for i in descriptionlist:
								if not label:
									label=i
								else:
									label=label+","+i
							fout.write("%s\n"%label)
		self.filelocation=self.filelocation+'.csv'
		self.filetype='CSV'

	#generator constructor for simple file reading
	#dummy_var=>I had to call immunogrepfile.read() in a function , but that function assumed there had to be a parameter passed into read function. So I set a dummy_var that never gets Used
	def read(self,dummy_var=None):
		if self.IFclass == None:
			return
		while not(self.IFclass.eof):		
			line_data = self.IFclass.read()				
			if line_data:
				yield line_data

	def __str__(self):
		return "Please refer to http://bigg.icmb.utexas.edu:5000/immunogrepFile/ for Manual."


class immunogrepFASTA():
	def __init__(self,filelocation,filetype,decoratorinfo=FASTAdecoratorinfo,contains_header=False,mode='r',field_names=[],chunk_size=1):
		self.filelocation=filelocation
		fasta_file_comment = "#" #standard comment for all fastafiles
		self.filetype=filetype
		self.field_names=field_names
		self.chunk_size=chunk_size
		self.contains_header = False
		if decoratorinfo==None:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip()
				if descriptor_symbol in tmpline:
					try:
						self.decoratorinfo=json.loads(tmpline[len(descriptor_symbol):].replace("'",'"'))
					except:
						self.decoratorinfo={}
				else:
					self.decoratorinfo=FASTAdecoratorinfo
		else:
			self.decoratorinfo=decoratorinfo
		self.mode=mode
		self.filehandle=open(self.filelocation,self.mode)
		self.SeqIOfilehandle=SeqIO.parse(self.filehandle,'fasta')
		#By-passing the decorator line if it exists
		self.temphandle=open(self.filelocation,"r")
		if descriptor_symbol in self.temphandle.readline():
			self.SeqIOfilehandle=SeqIO.parse(self.temphandle,'fasta')
		self.eof=False
		
	
	def Extract_Header(self):
		
		linenum = 0
		readunit=dict()
		endchk=False
		f = open(self.filelocation,'r')
		
		line = f.readline()		
		if descriptor_symbol in line:
			line = f.readline()
			if self.contains_header:
				line = f.readline()
		
		numSeq = 0
		fields = []
		totallines = 0
		fail = False
		try:
			while line and numSeq<num_line_search and not(fail):
				
				line = line.strip()			
				if len(line)>=1 and '>' in line[0]:					
					[header,more_fields] = GetAdditionalInfo(line[1:])					
					fields+=more_fields.keys()
					
					#if fasta_file_delimiter in line:
					#	tmp=line.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter													
					#	additional_info = json.loads(tmp[1])						
					#	fields = fields+additional_info.keys()
					numSeq+=1
				line = f.readline()
				totallines+=1
				if totallines == num_line_search*2 and numSeq == 0:
					fail = True				
			f.close()
		
			#unique_fields = sorted(set(fields))
			fields =['header','sequence','document_header']+sorted(set(fields))
			fields = sorted(set(fields))
		
			if fail:
				return None
			else:				
				return fields
		except Exception as e:
			print e
			return None
	
	def convert2IF(self):
		if self.filetype=='FASTA':
			with open(self.filelocation+'.'+IFExt,"w") as fout:
				fout.write(descriptor_symbol+json.dumps(self.decoratorinfo)+'\n')
				decor=self.decoratorinfo
				body=dict()
				for k in decor[description_var]:
					body[k]=''
				with open(self.filelocation,"r") as f:
					for line in f:
						if '>' in line:
							if body['header']!='' and body['sequence']!='':
								fout.write(json.dumps(body)+'\n')
							if fasta_file_delimiter in line:
								tmp=line.split(fasta_file_delimiter)
								addbody=json.loads(tmp[1])
								if idIdentifier in addbody:
									addbody[idIdentifier]=addbody[idIdentifier]
								for k in addbody:
									if k not in self.decoratorinfo[description_var]:
										self.decoratorinfo[description_var].append(k)
									body[k]=addbody[k]
								body['header']=tmp[0][1:]
							else:
								body['header']=line[1:].strip()
							body['sequence']=''
						else:
							body['sequence']=body['sequence']+line.strip()
					if body['header']!='' and body['sequence']!='':
						fout.write(json.dumps(body)+'\n')
			self.filelocation=self.filelocation+'.'+IFExt
			self.filetype='JSON'
			if len(self.decoratorinfo[description_var])>2:
				os.system("sed -i '1 c\%s' %s"%(descriptor_symbol+json.dumps(self.decoratorinfo),self.filelocation))


	def read(self):
		readunit=dict()
		
		try:
			line=self.SeqIOfilehandle.next()
		except:
			self.eof=True
			readunit=None
			
		if not self.eof:						
			readunit['sequence'] = str(line.seq)
			header_line = line.description
			readunit['document_header'] = header_line
			if fasta_file_delimiter in header_line:
				tmp=header_line.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter
				readunit['header']=tmp[0]									
				additional_info = json.loads(tmp[1])						
		
				for keyval in additional_info:
					if keyval == idIdentifier:
						if type(additional_info[keyval]) is dict:							
							readunit[idIdentifier] = additional_info[keyval]['$oid']
						else:							
							readunit[idIdentifier] = additional_info[keyval]								
					else:						
						readunit[keyval] = additional_info[keyval]
			else:
		 		readunit['header']=header_line
		return readunit

	def create(self):
		self.filehandle.write("%s\n"%(descriptor_symbol+json.dumps(self.getDecoratorinfo())))

	def write(self,content,additionalinfo=None):
		try:
			json.dumps(content)
		except:
			sys.exit("Write Error: Content is not conforming to the JSON format!")
		self.filehandle.write(">%s\n%s\n"%(content['header']+fasta_file_delimiter+additionalinfo,content['sequence']))

	def close(self):
		self.filehandle.close()
		self.SeqIOfilehandle.close()
		self.temphandle.close()

	def getDecoratorinfo(self):
		return self.decoratorinfo

	def getDescription(self):
		return self.decoratorinfo[description_var]

	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return self.filetype

	def getTranslation(self):
		if self.getDecoratorinfo()[translation_var]:
			t=self.getDecoratorinfo()[translation_var]
		else:
			t=None
		return t
		

class immunogrepFASTQ():
	#working version 211...
	def __init__(self,filelocation,filetype,decoratorinfo=FASTQdecoratorinfo,contains_header=False,mode='r',field_names=[],chunk_size=1):
		self.filelocation=filelocation
		self.filetype=filetype
		self.field_names=field_names
		self.chunk_size=chunk_size
		self.contains_header = False
		if decoratorinfo==None:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip()
				if descriptor_symbol in tmpline:
					try:
						self.decoratorinfo=json.loads(tmpline[len(descriptor_symbol):].replace("'",'"'))
					except:
						self.decoratorinfo={}
				else:
					self.decoratorinfo=FASTQdecoratorinfo
		else:
			self.decoratorinfo=decoratorinfo
		self.mode=mode
		self.filehandle=open(self.filelocation,self.mode)
		#self.SeqIOfilehandle=SeqIO.parse(self.filehandle,'fastq')
		#By-passing the decorator line if it exists
		
		self.temphandle=open(self.filelocation,"r")
		skip_lines = 0
		
		while True:			
			new_line = self.temphandle.readline().strip()
			if not(new_line):
				break
			if new_line.startswith('@'):
				break
			else:
				skip_lines+=1
			if new_line.startswith(descriptor_symbol) or new_line.startswith('#'):
				continue
			else:
				raise Exception('This FASTQ file, '+self.filelocation+', has unusual header lines: '+new_line)
		
		self.SeqIOfilehandle = open(self.filelocation,'r')
		
		for _ in range(skip_lines):
			self.SeqIOfilehandle.readline()
					
		
	#	if descriptor_symbol in self.temphandle.readline():
#			#self.SeqIOfilehandle=SeqIO.parse(self.temphandle,'fastq')
#			self.SeqIOfilehandle=open(self.filelocation,'r') #SeqIO.parse(self.filelocation,self.mode)
#		else:
#			#self.SeqIOfilehandle=SeqIO.parse(self.filehandle,'fastq')
#			self.SeqIOfilehandle=open(self.filelocation,'r') #SeqIO.parse(self.filelocation,self.mode)
		
		self.eof=False

	def convert2IF(self):
		if self.filetype=='FASTQ':
			with open(self.filelocation+'.'+IFExt,"w") as fout:
				fout.write(descriptor_symbol+json.dumps(self.decoratorinfo)+'\n')
				decor=self.decoratorinfo
				body=dict()
				for k in decor[translation_var]:
					body[decor[translation_var][k]]=''
				with open(self.filelocation,"r") as f:
					for line in f:
						if '@' in line[0]:
							if body['header']!='' and body['sequence']!='' and body['phred']!='':
								fout.write(json.dumps(body)+'\n')
							body['header']=line[1:].strip()
							body['sequence']=''
							body['phred']=''
							plus=False
						elif not plus and '+' not in line[0]:
							body['sequence']=body['sequence']+line.strip()
						elif '+' in line[0]:
							plus=True
							continue
						elif plus:
							body['phred']=body[decor[translation_var]['QUAL_SCORE']]+line.strip()
					if body['header']!='' and body['sequence']!='' and body['phred']!='':
						fout.write(json.dumps(body)+'\n')
			self.filelocation=self.filelocation+'.'+IFExt
			self.filetype='JSON'
	
	def Extract_Header(self):
		try:
			linenum = 0
			readunit=dict()
			endchk=False
			f = open(self.filelocation,'r')
			
			line = f.readline()		
			if descriptor_symbol in line:
				line = f.readline()
				if self.contains_header:
					line = f.readline()
			
			numSeq = 0
			fields = []
			fail = False
			totallines=0
			while line and numSeq<num_line_search and not(fail):
				line = line.strip()			
				if '@' in line[0]:
					[header,more_fields] = GetAdditionalInfo(line[1:])					
					fields+=more_fields.keys()
					#if fasta_file_delimiter in line:
					#	tmp=line.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter													
					#	additional_info = json.loads(tmp[1])						
					#	fields = fields + additional_info.keys()
					numSeq+=1
				line = f.readline()
				totallines+=1
				if totallines == num_line_search*4 and numSeq == 0:
					fail = True
		
			
			f.close()
			#unique_fields = sorted(set(fields))
			fields = sorted(set(['header','sequence','phred','document_header']+sorted(set(fields))))
		
			if fail:
				return None
			else:
				return fields
		except:
			return None		

	def read(self):
		readunit=dict()				
				
		try:
			line1 = self.SeqIOfilehandle.readline().strip()
		except:
			self.eof=True
			readunit=None
		
		if not line1:
			self.eof=True
			readunit = None
		else:				
			if line1[0] != '@':
				raise Exception('This does not appear to be a FASTQ sequence as the sequence header does not start with "@": {0}'.format(line1))			
			header = line1[1:]
			seq=''
			next_line = self.SeqIOfilehandle.readline().strip()
			while next_line!='+':
				seq+=next_line
				next_line = self.SeqIOfilehandle.readline().strip()
			len_seq = len(seq)
			next_line = self.SeqIOfilehandle.readline().strip()
			qual=next_line
			qual_len = len(qual)			
			while qual_len<len_seq:
				qual+=self.SeqIOfilehandle.readline().strip()
				qual_len = len(qual)
			if qual_len!=len_seq:
				raise Exception('The sequence length does not match the length of the quality score for the following sequence: {0}'.format(header))
			#line = self.SeqIOfilehandle.next()
		
			
		if not self.eof:						
			readunit['sequence'] = seq # str(line.seq)
			#qualnum=line.letter_annotations["phred_quality"]
			#readunit['phred']=''
			#for i in qualnum:
			#	readunit['phred']=readunit['phred']+str(chr(i+33))
			readunit['phred'] = qual
			header_line = header #line.description
			readunit['document_header'] = header_line
			if fasta_file_delimiter in header_line:
				tmp=header_line.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter
				readunit['header']=tmp[0]									
				additional_info = json.loads(tmp[1])						
		
				for keyval in additional_info:
					if keyval == idIdentifier:
						if type(additional_info[keyval]) is dict:							
							readunit[idIdentifier] = additional_info[keyval]['$oid']
						else:							
							readunit[idIdentifier] = additional_info[keyval]								
					else:						
						readunit[keyval] = additional_info[keyval]
			else:
		 		readunit['header']=header_line
		return readunit

	def create(self):
		self.filehandle.write("%s\n"%(descriptor_symbol+json.dumps(self.getDecoratorinfo())))

	def write(self,content,additionalinfo=None):
		try:
			json.dumps(content)
		except:
			sys.exit("Write Error: Content is not conforming to the JSON format!")
		self.filehandle.write("@%s\n%s\n+\n%s\n"%(content['header']+fasta_file_delimiter+additionalinfo,content['sequence'],content['phred']))

	def close(self):
		self.filehandle.close()
		self.SeqIOfilehandle.close()
		self.temphandle.close()

	def getDecoratorinfo(self):
		return self.decoratorinfo

	def getDescription(self):
		return self.decoratorinfo[description_var]

	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return self.filetype

	def getTranslation(self):
		if self.getDecoratorinfo()[translation_var]:
			t=self.getDecoratorinfo()[translation_var]
		else:
			t=None
		return t

class immunogrepDELIM():
	def __init__(self,filelocation,filetype,decoratorinfo,delimiter,contains_header=True,mode='r',field_names=[],chunk_size=1):
		self.filelocation=filelocation
		self.filetype=filetype
		self.delimiter=delimiter
		self.contains_header=contains_header
		descriptorPresent = False
		self.field_names=field_names
		self.chunk_size=chunk_size		
		if decoratorinfo==None:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip('\r\n')				
				if descriptor_symbol in tmpline:
					descriptorPresent = True
					try:
						self.decoratorinfo=json.loads(tmpline[len(descriptor_symbol):].replace("'",'"'))
					except:
						self.decoratorinfo = {}
				else:
					self.decoratorinfo={}			
		else:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip('\r\n')				
				if descriptor_symbol in tmpline:
					descriptorPresent = True				
			self.decoratorinfo=decoratorinfo
					
		if contains_header:
			self.header_skipped = False
			with open(self.filelocation,"r") as f:
				line=f.readline().strip('\r\n')
				if descriptor_symbol in line:
					line=f.readline().strip('\r\n')
				tmp=line.split(self.delimiter)
				self.decoratorinfo[description_var] = [f for f in tmp]#filed names detected in file 
				self.header_names_map={f:col for col,f in enumerate(tmp) if f}#dictionary mapping header name to column number
				self.maxfields = len(tmp)					
		else:
			self.header_skipped = True
			self.maxfields=0
			with open(self.filelocation,'r') as f:
				for line in range(num_line_search):					
					line = f.readline().strip('\r\n')
					if descriptor_symbol in line:
						continue
					tmp=line.split(delimiter)
					if len(tmp)>self.maxfields:
						self.maxfields=len(tmp)
			descriptor = ['Column '+str(i+1) for i in range(self.maxfields)]#field names detected in file 
			self.header_names_map = {'Column '+str(i+1):i for i in range(self.maxfields)}#dictionary mapping header name to column number
			self.decoratorinfo[description_var]=descriptor		

		self.mode=mode
		self.filehandle=open(self.filelocation,self.mode)
		
		#if field names were provided, then filter the dictionary header map by field names
		if self.field_names:
			self.header_names_map = {f:v for f,v in self.header_names_map.iteritems() if f in self.field_names}
		
		if descriptorPresent:
			self.filehandle.readline()
		self.eof=False

	def Extract_Header(self):		
		return self.decoratorinfo[description_var]

	def convert2IF(self):
		if self.filetype=='TAB' or self.filetype=='CSV' or self.filetype=='DELIM':
			with open(self.filelocation+'.'+IFExt,"w") as fout:
				fout.write(descriptor_symbol+json.dumps(self.decoratorinfo)+'\n')
				decor=self.decoratorinfo
				body=dict()
				for k in decor[description_var]:
					body[k]=''
				with open(self.filelocation,"r") as f:
					tmpline=f.readline().strip('\r\n')
					if descriptor_symbol in tmpline:
						f.readline()
					if self.contains_header:
						f.readline()
					for line in f:
						if line=='':
							continue
						line=line.strip('\r\n')
						tmp=line.split(self.delimiter)
						for i in range(len(tmp)):
							body[decor[description_var][i]]=tmp[i]
						fout.write(json.dumps(body)+'\n')
			self.filelocation=self.filelocation+'.'+IFExt
			self.filetype='JSON'

	def read(self):
	
		if not(self.header_skipped):
			self.filehandle.readline()
			self.header_skipped = True
							
		line=self.filehandle.readline().strip('\r\n')
		
		if not line:
			self.eof=True
			return None
				
		tmp=line.split(self.delimiter)
		max_len = len(tmp)
		readunit = {field:tmp[col_num] if col_num < max_len else '' for field,col_num in self.header_names_map.iteritems() }								
		
		return readunit

	def create(self):
		self.filehandle.write("%s\n"%(descriptor_symbol+json.dumps(self.getDecoratorinfo())))

	def write(self,content,additionalinfo=None):
		try:
			json.dumps(content)
		except:
			sys.exit("Write Error: Content is not conforming to the JSON format!")
		tmp=self.getDecoratorinfo()[description_var]
		
		for k in tmp:
			if not body:
				body=content[k]
			else:
				body=body+self.delimiter+content[k]
		body=body+"\n"
		self.filehandle.write(body)

	def close(self):
		self.filehandle.close()

	def getDecoratorinfo(self):
		return self.decoratorinfo

	def getDescription(self):
		return self.decoratorinfo[description_var]

	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return self.filetype

	def getTranslation(self):
		if self.getDecoratorinfo()[translation_var]:
			t=self.getDecoratorinfo()[translation_var]
		else:
			t=None
		return t

class immunogrepJSON():
	def __init__(self,filelocation,filetype,decoratorinfo,contains_header=False,mode='r',field_names=[],chunk_size=1):
		self.filelocation=filelocation
		self.filetype=filetype
		descriptorPresent = False
		self.field_names=field_names
		self.chunk_size=chunk_size
		if decoratorinfo==None:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip()
				if descriptor_symbol in tmpline:
					descriptorPresent = True
					try:
						self.decoratorinfo=json.loads(tmpline[len(descriptor_symbol):].replace("'",'"'))
					except:
						self.decoratorinfo={}						
				else:
					self.decoratorinfo={}
			if description_var not in self.decoratorinfo:
				self.decoratorinfo[description_var] = self.Extract_Header() 
		else:
			with open(self.filelocation,"r") as f:
				tmpline=f.readline().strip()				
				if descriptor_symbol in tmpline:
					descriptorPresent = True		
			self.decoratorinfo=decoratorinfo

		self.mode=mode
		self.filehandle=open(self.filelocation,self.mode)
		if descriptorPresent:
			self.filehandle.readline()
		self.eof=False
	
	
	def Extract_Header(self):
		try:
			linenum = 0
			readunit=dict()
			endchk=False
			f = open(self.filelocation,'r')
			line = f.readline()		
			if descriptor_symbol in line:
				line = f.readline()			
			numSeq = 0
			fields = []			
			while line and numSeq<num_line_search:
				line = line.strip()			
				
				additional_info = json.loads(line)
				fields= fields+additional_info.keys()
				numSeq+=1
				line = f.readline()	
			f.close()
			
			unique_fields = sorted(set(fields))
			return unique_fields
		except:	
			return None
		
	def read(self):
		line=self.filehandle.readline()
		while descriptor_symbol in line:
			line=self.filehandle.readline()
		if not line:
			self.eof=True
			return None
		tmp=json.loads(line)
		return tmp

	def create(self):
		self.filehandle.write("%s\n"%(descriptor_symbol+json.dumps(self.getDecoratorinfo())))

	def write(self,content):
		try:
			json.dumps(content)
		except:
			sys.exit("Write Error: Content is not conforming to the JSON format!")
		self.filehandle.write(json.dumps(content)+"\n")

	def close(self):
		self.filehandle.close()

	def getDecoratorinfo(self):
		return self.decoratorinfo

	def getDescription(self):
		return self.decoratorinfo[description_var]

	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return self.filetype

	def getTranslation(self):
		if self.getDecoratorinfo()[translation_var]:
			t=self.getDecoratorinfo()[translation_var]
		else:
			t=None
		return t


class immunogrepIMGT():
	def __init__(self, imgtfiles,important_headers_only =True,field_names=[],required_files=0,include_parameters_file=False,chunk_size=1):
		self.filelocation=imgtfiles
		self.filetype='IMGT'				
		self.chunk_size=chunk_size
		self.mode='r'
		if not isinstance(required_files,list):
			self.required_files = [required_files]
		else:
			self.required_files = required_files
		self.include_parameters_file = include_parameters_file
		self.important_headers_only=important_headers_only
		self.field_names = field_names
		self.pat = re.compile(r"\(([0-9]+)\)") #search for any numbers within '(' and ')' use \( to make string literal
						
		#these are redundant columns in multiple IMGT files 
		self.identical_fields = [
		#	'Sequence ID',	
		#	'V-GENE and allele',
		#	'J-GENE and allele',
		#	'D-GENE and allele'
		]

		self.special_fields = {
			'V-GENE and allele_1':lambda x:x.replace(', or',',').replace(' or ',',').replace('(see comment)',''),
			'D-GENE and allele_1':lambda x:x.replace(', or',',').replace(' or ',',').replace('(see comment)',''),		
			'J-GENE and allele_1':lambda x:x.replace(', or',',').replace(' or ',',').replace('(see comment)',''),	
			'CDR1-IMGT_3':lambda x:x.replace('.',''),
			'CDR1-IMGT_5':lambda x:x.replace('.',''),
			'CDR2-IMGT_3':lambda x:x.replace('.',''),
			'CDR2-IMGT_5':lambda x:x.replace('.',''),
			'CDR3-IMGT_3':lambda x:x.replace('.',''),
			'CDR3-IMGT_5':lambda x:x.replace('.',''),
			'FR1-IMGT_3':lambda x:x.replace('.',''),
			'FR1-IMGT_5':lambda x:x.replace('.',''),
			'FR2-IMGT_3':lambda x:x.replace('.',''),
			'FR2-IMGT_5':lambda x:x.replace('.',''),
			'FR3-IMGT_3':lambda x:x.replace('.',''),
			'FR3-IMGT_5':lambda x:x.replace('.',''),		
			'FR4-IMGT_3':lambda x:x.replace('.',''),
			'FR4-IMGT_5':lambda x:x.replace('.',''),
			'V-REGION Nb of mutations_8':lambda x: self.pat.findall(x)[0] if '(' in x else x.strip(),
			'V-REGION Nb of AA changes_9':lambda x: self.pat.findall(x)[0] if '(' in x else x.strip()		
		}

		
		#first seperate the provided IMGT files by their file name. 
		#group all files from the same original filename together 
		self.grouped_imgt_files = GroupIMGTFiles(imgtfiles)
		all_field_names = []
		#go through each IMGT group. Make sure the required files are present
		for each_group,found_files in self.grouped_imgt_files.iteritems():				
			data_to_read = [{}]*11						
			for file in found_files:# imgtfiles:			
				#go through each of the found files in that group 
				#determine each files  file number 
				#also make sure to only read in fields requested by user 
				[filenumber,fields] = find_imgt_file_type_index(file,self.important_headers_only)				
				if not filenumber:				
					continue
				if filenumber==11:				
					if self.include_parameters_file == False:
						read_fields = []
					else:
						read_fields = fields				
				elif self.field_names:											
					fields = set([f+'_'+str(filenumber) for f in fields])
					#only read the intersection of fields  
					read_fields = list( fields&set(self.field_names) )										
				else:
					fields = set([f+'_'+str(filenumber) for f in fields])
					read_fields = fields	
				all_field_names.extend(read_fields)
				
				if read_fields:
					data_to_read[filenumber-1] = True
																		
			#ENSURE THAT THE REQUIRED IMGT FILES ARE PRESENT		
			for each_required_file in self.required_files:
				if each_required_file>0 and each_required_file<=11:
					#user requested that specific IMGT file is required for reading 
					if data_to_read[each_required_file-1]=={}:
						raise Exception('IMGT File {0} from {1} is required'.format(str(each_required_file),each_group))				
				
		
		#ITS IMGT so we already know what to expect 
		self.decoratorinfo = {
			description_var:sorted(list(set(all_field_names))),
			translation_var:imgt_db_translator
		}
					
		self.eof=False
		self.my_gen = self.IMGT_Reader_Gen()
	
	def CheckForDBID(self,sequence_header):	
		additional_info = {}
		if fasta_file_delimiter in sequence_header:		
			tmp=sequence_header.split(fasta_file_delimiter) #if extra info is included in the header, then it is passed using the the fasta_file_delimiter
			sequence_header=tmp[0]					
			try:
				additional_info = json.loads(tmp[1])
			except:
				#usually occurs if sequence header is truncated 
				additional_info = {}
			
		return additional_info	#if the ID idnetifier was included form th atabase then it will be in this dict

	
	def treat_these_fields_specially(self,key,value):	
		if key in self.special_fields:
			return self.special_fields[key](value)
		else:
			return value

	def treat_these_fields_specially2(self,key,value):
		if key=='V-GENE and allele_1' or key == 'D-GENE and allele_1' or key == 'J-GENE and allele_1':	
			return value.replace(', or',',').replace(' or ',',').replace('(see comment)',''),
		elif key in ['CDR1-IMGT_3','CDR2-IMGT_3','CDR3-IMGT_3','FR1-IMGT_3','FR2-IMGT_3','FR3-IMGT_3','CDR1-IMGT_5','CDR2-IMGT_5','CDR3-IMGT_5','FR1-IMGT_5','FR2-IMGT_5','FR3-IMGT_5']:
			return value.replace('.','')
		elif key in ['V-REGION Nb of mutations_8','V-REGION Nb of AA changes_9']:
			return self.pat.findall(value)[0] if '(' in value else value.strip()
		else:
			return value
	
	
	def Extract_Header(self):
		return self.decoratorinfo[description_var]	
		
	def read(self):
				
		try:
			readunit = self.my_gen.next()
		except Exception as e:
			if str(e).strip():				
				raise Exception('Unexpected error reading imgt file: '+str(e).strip())
			readunit = None
			self.eof=True
		return readunit
			
		
	#def close(self):		
	#	self.filehandle.close()

	def getDecoratorinfo(self):
		return self.decoratorinfo

	def getDescription(self):
		return self.decoratorinfo[description_var]

	def getFilelocation(self):
		return self.filelocation

	def getFiletype(self):
		return 'IMGT'

	def getTranslation(self):
		if self.getDecoratorinfo()[translation_var]:
			t=self.getDecoratorinfo()[translation_var]
		else:
			t=None
		return t
		
	def IMGT_Reader_Gen(self):		
		include_fields = self.field_names
		 
		#first seperate the provided IMGT files by their file name. 
		#group all files from the same original filename together 
		grouped_imgt_files = GroupIMGTFiles(self.filelocation)
		
		#next make sure all the required fieldnames are present 
		add_fields = []
		if 'PREDICTED_AB_SEQ.NT' in include_fields:
			add_fields.extend(['V-D-J-REGION_3','V-J-REGION_3']) #these fields from IGMT files will be required to make PREDICTED_AB_SEQ.NT
		if 'PREDICTED_AB_SEQ.AA' in include_fields:
			add_fields.extend(['V-D-J-REGION_5','V-J-REGION_5'])
		
		if 'VREGION.SHM.NT_PER' in include_fields:
			add_fields.extend(['V-REGION identity %_1'])
		if 'JREGION.SHM.NT_PER' in include_fields:
			add_fields.extend(['J-REGION identity %_1'])
		
		if 'VREGION.SHM.NT' in include_fields:
			add_fields.extend(['V-REGION identity nt_1'])
		
		if 'JREGION.SHM.NT' in include_fields:
			add_fields.extend(['J-REGION identity nt_1'])
		
		if 'VREGION.SHM.AA' in include_fields or 'VREGION.SHM.AA_PER' in include_fields:
			add_fields.extend(['CDR3-IMGT Nb of AA_9','CDR3-IMGT Nb of AA changes_9','V-REGION Nb of AA_9','V-REGION Nb of AA changes_9'])
		
		for new_fields in add_fields:
			if new_fields not in include_fields:
				include_fields.append(new_fields)
		
		#read results from each IMGT 'group'
		for each_group,found_files in grouped_imgt_files.iteritems():			
		
			data_to_read = [{}]*11				
			
			for file in found_files:# imgtfiles:
			
				#go through each of the found files in that group 
				#determine each files  file number 
				#also make sure to only read in fields requested by user 
				[filenumber,fields] = find_imgt_file_type_index(file,self.important_headers_only)				
				if not filenumber:				
					continue
				if filenumber==11:				
					if self.include_parameters_file == False:
						read_fields = []
					else:
						read_fields = fields				
				elif include_fields:												
					fields = set([f+'_'+str(filenumber) for f in fields])
					#only read the intersection of fields  
					read_fields = list( fields&set(include_fields) )
					
					#now that we have the intersection, remove the file number that was appended to fields 
					for i,rf in enumerate(read_fields):
						read_fields[i] = '_'.join(rf.split('_')[:-1])
				else:
					read_fields = fields
					
					
				#pprint(read_fields)
				#if we are reading fields from this file, then add to the list of files to read 
				if read_fields:
					data_to_read[filenumber-1] = {
						'filenum': filenumber,
						'fields':read_fields,
						'filepath':file,				
						'file_reader':immunogrepFile(filelocation = file,filetype='TAB',field_names=read_fields)				
					}		
						
					
			#ENSURE THAT THE REQUIRED IMGT FILES ARE PRESENT			
			for each_required_file in self.required_files:
				if each_required_file>0 and each_required_file<=11:
					#user requested that specific IMGT file is required for reading 
					if data_to_read[each_required_file-1]=={}:
						raise Exception('IMGT File {0} from {1} is required'.format(str(each_required_file),each_group))
			
			command_settings_keys = {
				'IMGT/V-QUEST programme version':'Version',
				'IMGT/V-QUEST reference directory release':'Database release',
				'Species':'Species',
				'Receptor type or locus':'Receptor or locus',
				'IMGT/V-QUEST reference directory set':'Database set',
				'Search for insertions and deletions':'Search for insertions and deletions'
				}
			#uppercase 
			command_settings_keys = {c.upper():val.upper() for c,val in command_settings_keys.iteritems()}
			
			#check if file 11 exists
			command=''
			#print command_settings_keys
			if data_to_read[10]:
				with open(data_to_read[10]['filepath'],'r') as e:
					command_lines = e.readlines()
				for each_command in command_lines:
					settings = each_command.split('\t')
					#column1 = > command/name/field/setting
					#column2 = > value			
					if settings[0].strip(' \r\t\n:').upper() in command_settings_keys and settings[1].strip():				
						#if we find the command we we want to store in the file and the second column is not empty, then add to command string 
						command+=command_settings_keys[settings[0].strip(' \r\t\n:').upper()]+':'+settings[1].strip().upper()+';'
														
			command = command.strip(';')
			
			at_least_one_file_open = True
			while at_least_one_file_open:
				row_data = defaultdict(str,{})
				
				at_least_one_file_open = False#if all files are eof, then thsi will be false at t thend of loop
				for files_dict in data_to_read[0:10]:		
					
					if not files_dict:				
						continue
					
					filenumber = files_dict['filenum']
					if filenumber<10 and not files_dict['file_reader'].IFclass.eof:								
						at_least_one_file_open = True
						current_file_data = files_dict['file_reader'].IFclass.read()				
						if current_file_data:					
							for each_field,each_value in current_file_data.iteritems():
								if each_field in files_dict['fields']:
									if each_field in self.identical_fields:								
										row_data['{0}'.format(each_field,str(filenumber))] = self.treat_these_fields_specially(each_field,each_value)							
										#row_data['{0}'.format(each_field,str(filenumber))] = each_value#treat_these_fields_specially(each_field,each_value)							
									else:								
										#row_data['{0}_{1}'.format(each_field,str(filenumber))] = each_value#treat_these_fields_specially(each_field,each_value)							
										row_data['{0}_{1}'.format(each_field,str(filenumber))] = self.treat_these_fields_specially(each_field+'_'+str(filenumber),each_value)							
				
				if at_least_one_file_open==False:					
					continue
				
				if 'V-D-J-REGION_3' in row_data:
					row_data['PREDICTED_AB_SEQ.NT'] = row_data['V-D-J-REGION_3'] if row_data['V-D-J-REGION_3'] else row_data['V-J-REGION_3']			
				if 'V-D-J-REGION_5' in row_data:
					row_data['PREDICTED_AB_SEQ.AA'] = row_data['V-D-J-REGION_5'] if row_data['V-D-J-REGION_5'] else row_data['V-J-REGION_5']		
				if 'Sequence ID' in row_data: #extract any other information included in sequence header including database key information
					more_fields = self.CheckForDBID(row_data['Sequence ID'])
					row_data = dict(row_data.items()+more_fields.items())
				
												
				
				row_data['VREGION.SHM.NT_PER'] = ''
				if 'V-REGION identity %_1' in row_data and row_data['V-REGION identity %_1']:												
					try:
						row_data['VREGION.SHM.NT_PER'] = 100-float(row_data['V-REGION identity %_1'])					
					except:					
						row_data['VREGION.SHM.NT_PER']=''
				
				row_data['VREGION.SHM.NT'] = ''
				
				if 'V-REGION identity nt_1' in row_data and row_data['V-REGION identity nt_1']:
					try:
						match_vals = row_data['V-REGION identity nt_1'].split('/')
						match_vals[1] = match_vals[1].split(' ')[0]
						row_data['VREGION.SHM.NT'] = int(match_vals[1])-int(match_vals[0])
					except:					
						row_data['VREGION.SHM.NT'] = ''
				
	
				row_data['JREGION.SHM.NT'] = ''
				
				if 'J-REGION identity nt_1' in row_data and row_data['J-REGION identity nt_1']:
					try:
						match_vals = row_data['J-REGION identity nt_1'].split('/')
						match_vals[1] = match_vals[1].split(' ')[0]
						row_data['JREGION.SHM.NT'] = int(match_vals[1])-int(match_vals[0])
					except:
						row_data['JREGION.SHM.NT'] = ''
										
	
				row_data['JREGION.SHM.NT_PER'] = ''
				if 'J-REGION identity %_1' in row_data and row_data['J-REGION identity %_1']:								
					try:
						row_data['JREGION.SHM.NT_PER']=100-float(row_data['J-REGION identity %_1'])
					except:
						row_data['JREGION.SHM.NT_PER']=''					
				
							
				try:
					
					if row_data['CDR3-IMGT Nb of AA_9'].strip() == '-':					
						row_data['CDR3-IMGT Nb of AA_9'] = '0'
					if row_data['CDR3-IMGT Nb of AA changes_9'].strip() == '-':					
						row_data['CDR3-IMGT Nb of AA changes_9'] = '0'
									
					num_aa_v = float(row_data['V-REGION Nb of AA_9'].split('(')[0].strip()) - float(row_data['CDR3-IMGT Nb of AA_9'].split('(')[0].strip())
					num_aa_mismatch =  float(row_data['V-REGION Nb of AA changes_9']) - float(row_data['CDR3-IMGT Nb of AA changes_9'].split('(')[0].strip())
					row_data['VREGION.SHM.AA'] = num_aa_mismatch
					row_data['VREGION.SHM.AA_PER']=round(100*num_aa_mismatch/num_aa_v,4)
				except Exception as e:				
					row_data['VREGION.SHM.AA_PER']=''
					row_data['VREGION.SHM.AA']=''
				
				
				#add recombination_type to fields 
				row_data['R_TYPE']  = 'VDJ' if row_data['V-D-J-REGION_3'] else ('VJ' if row_data['V-J-REGION_3'] else ('' if not(row_data['V-GENE and allele_1']) else ('VDJ' if ('IGHV' in row_data['V-GENE and allele_1'] or 'TRBV' in row_data['V-GENE and allele_1']) else 'VJ')))
				
				
				if command:
					row_data['COMMAND'] = command
				#if 'Sequence ID_3' in row_data:
				#	more_fields = CheckForDBID(row_data['Sequence ID_3'])
				#	row_data = dict(row_data.items()+more_fields.items())
				#if 'Sequence ID_8' in row_data:
				#	more_fields = CheckForDBID(row_data['Sequence ID_8'])
				#	row_data = dict(row_data.items()+more_fields.items())				
				if row_data:
					yield row_data
	
		
		
