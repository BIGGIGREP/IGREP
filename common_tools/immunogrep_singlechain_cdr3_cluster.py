import csv
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import json
import datetime
import numpy as np
from operator import itemgetter
import os
import gc
import time
from collections import deque
import sys 
import math

#for making folders in igrep
import immunogrep_file_system_tools as filesystem
#for reading files
import immunogrep_read_file as readfile 

#default dictionaries
from collections import defaultdict
#for printing pretty
from pprint import pprint

from immunogrep_global_variables import translation_var #this key will signify the translation/translator key 
from immunogrep_global_variables import descriptor_symbol #this key will signify the translation/translator key 
from immunogrep_global_variables import idIdentifier #variable defining what SEQ_ID variable is in database 
from immunogrep_global_variables import fasta_file_delimiter
#import immunogrep_useful_immunogrep_functions as useful
import immunogrep_useful_functions as useful

from random import randint

import subprocess

#VERSION 4.0 applied on Version 7
pairing_version = 'Pairing_v4.1'
#updates => allow isotyping 
	#=> store SHM as mean SHM and variance 
	#=> allow diversity calculations 
#updates 4.0 => produce a pairing 'annotation' file 
	#the annotation file will store the 'paired-id' of successful pairs , and their corresponding 'paired cluster id'
#updates 4.1 => 
	#use 'MODE' for VGENE usage call:
		#ONLY CONSIDER VGENES WHOSE SCORES ARE EQUAL TO THE TOP SCORE!!!
		#RULE/ASSUMPTION: IF THERE IS NOT VGENE SCORE FIELD OR THE VGENE SCORE FIELD LENGTH != THE LENGHT OF VGENES THEN ONLY CONSIDER ALL VGENES PROVIDED
	#we will no longer filter out sequences whose full length nt or aa field are empty 
	#add a new productivity rule => productive if there is no stop codon in the CDR3 sequence ONLY. This will become the default productivity rule  
	

def get_ab_field_loc():
	#structure of array in codeDict and collapsedict 
	#every R1-R2 read (that passes filters) will be stored as an array in the variable codeDict
	#ONLY the first occurrence of H or L data from a pair will be stored in codeDict. 
	#the location of each data stored for either the heavy or light chain will be mapped by the following variable
	#for example, lets if the Heavy chain is encountered first. It will be stored in codeDict as an array based on the indexes below. 
				 #the moment its Light chain pair is encountered, then the heavy chain data will be removed from the codeDict var and 
				 #both heavy and light chain data will be stored in collapsed_dict as a 2-D array where index 0 is heavy chain and index 1 is light chain . The third element of collapsed dict stores counts
				 #-> collapsedict[fullcode] = [heavy chain data array/list, light chain data array/list, COUNTS!!]			 

	ab_field_loc = {
		'CDR3_SEQ':0, #store the cdr3 sequence in the first element of the array 
		'VGENE':1,
		'DGENE':2,
		'JGENE':3,
		'LOCUS':4,
		'MUT':5,
		'CDR3_LEN':6,
		'AB_SEQ':7,
		'CHAIN_CALL':8,
		'ISOTYPE':9,
		'VSCORES':10,
		'JSCORES':11,
		'DSCORES':12
	}
	return ab_field_loc

def get_required_field_names():
	required_field_names = {
			'vgene':'',
			'jgene':'',
			'dgene':'',
			'vgene_scores':'',
			'jgene_scores':'',
			'dgene_scores':'',
			'raw_seq_nt':'',
			'cdr3_nt':'',
			'seq_header':'',
			'shm':'',
			'full_len_ab_aa':'',
			'full_len_ab_nt':'',
			'functionality':'',
			'recomb':''
		}		
	return required_field_names

#the array for each heavy/light chain will always be this long
num_elem_stored = max([v for i,v in get_ab_field_loc().iteritems()])+1  
allowed_file_formats = ['TAB','CSV','JSON','IMGT']

#We will need to update the database with the results from IgFFT.  In order
#to update teh database, we need a translator, so that we know what fields go where in the database
def DatabaseTranslator(input_dictionary = {}):
	key = translation_var
	
	translator = {								
			"RECOMBINATION_FIELD":{ #THIS TELLS THE PROGRAM HOW TO DETERMINE WHETHER AN ANALYSIS/QUERY RESULT (from file) IS VDJ OR VJ
					"FIELD_NAME": "Recombination Type", #name of the field in the file that will give information regarding the recombination type (VDJ OR VJ)					
					"EXPLICIT":True,
			},
			"FIELDS_UPDATE":{ 
				#key = field name  in database 
				#value = field name in file
				#this will map all of the fields in the file to the proper location in the database. I.E. If I list VGENES as the column name/field name, then i want to map VREGION.VGENES:VGENES (because VREGION.VGENES is the name in the database)							
				idIdentifier:idIdentifier,
				"SEQUENCE":"Sequence",				
				"SEQUENCE_HEADER":"Header",								
				"PAIRING.PAIRED_ID":'PairedID',
				"PAIRING.PAIRED_CLUSTER":'PairedClusterId',
				"PAIRING.CONFIDENCE":'Confidence',
				"PAIRING.DOMINANCE":'Dominance',
				'COMMAND':'Command'
			}
	}
	
	input_dictionary[key] = translator		
	
	return input_dictionary

def FloatRange(minv,maxv,step,num_dec=None):
	if step==0:
		step =0.01
	minv = max(0,float(minv))
	maxv = min(1,float(maxv))
	
	x = minv
	if maxv<minv:
		min2=minv
		minv=maxv
		maxv=min2
	if minv==maxv:
		yield minv
		x+=1
		minv+=1
		#yield min+1
	else:
		while x<=maxv:
			if num_dec:
				yield round(x,num_dec)
			else:
				yield x
			x+=float(step)
	yield maxv

#use output from IMGT program to determine productivity
def IMGTProductivity(field_row):
	global functionality_field
	if field_row[functionality_field] == 'unproductive' or field_row[functionality_field] == 'No results' or not(field_row[functionality_field]):
		return False
	else:
		return True

#use output from IGBLAST program to determine productivity
#(have not handled this function yet)
def IGBLASTProductivity(field_row):
	return True 

#use output from IGFFT program to determine productivity
def IGFFTProductivity(field_row):
	global functionality_field
	if field_row[functionality_field].upper() == 'NO' or not(field_row[functionality_field]):
		return False
	else:		
		return True

#use output from MIXCR program to determine productivity
def MIXCRProductivity(field_row):	
	global full_aa_field
	if '*' in field_row[full_aa_field] or '_' in field_row[full_aa_field]:
		return False
	else:
		return True


#general function for determining if a sequence is productive or not 
def GeneralProductivity(field_row):
	"""
		General function for determining if an entire Antibody amino acid sequence is productive or not 
	
		Function will return True or False based on rule 
		Input variable field_row: this will represent all antibody annotation field for a specific sequence the provided annotation file
		
		General steps:
			1) Extract the field referring to the amino acid sequence of the antibody 			
			4) If there is a stop codon (*) in sequence 
				=> return false 
			5) return true
	
	"""
	global full_aa_field
	if '*' in field_row[full_aa_field]:
		return False
	else:
		return True

def CDR3Productivity(field_row):
	"""
	
		Rule for CDR3 productivity. 		
		
		Function will return True or False based on rule 
		Input variable field_row: this will represent all antibody annotation field for a specific sequence the provided annotation file
		
		General steps:
			1) Extract the CDR3 nucleotide field 
			2) If length of CDR3 is not a multiple of 3 and length of cdr3 < 9
				=> return false 
			3) Translate CDR3 to amino acid 
			4) If there is a stop codon (*) in sequence
				=> return false 
			5) return true
		
	"""
	global cdr3_field
	cdr3_nt = field_row[cdr3_field]
	if len(cdr3_nt)%3!=0 or len(cdr3_nt)<9:
		return False		
	cdr3_aa = str(Seq(cdr3_nt,generic_dna).translate())
	if '*' in cdr3_aa:
		return False
	else:
		return True

	
def populate_clusters(filename,delimiter='\t'):
	d = {}
	if not(os.path.isfile(filename)):
		return d
	with open(filename) as f_in:
		for line_text in f_in:			
			line_text = line_text.strip()
			line = line_text.split(delimiter)
		
			SorH=line[0] #Seed or Hit
			number=line[1] #seed number
			
			if ":" in line[8]:
				info=line[8]
			else:
				info=line[9]
			counts=int(info.split(":")[0])
			
			#identify seeds (S) vs hits (H)
			if SorH=='S':
				d[number]=[counts,]
				d[number].append(info)
			if SorH=='H':
				d[number][0]+=counts
				d[number].append(info)		
	return d 

def gene_hist(d,filename,geneCall):
	
	summary = open(filename,'a')
	
	hist=dict()
	
	if geneCall==1:
		Title='VH gene analysis'
		genePos=5
	if geneCall==2:
		Title='DH gene analysis'
		genePos=6
	if geneCall==3:
		Title='JH gene analysis'
		genePos=7
	if geneCall==4:
		Title='VL gene analysis'
		genePos=8
	if geneCall==5:
		Title='JL gene analysis'
		genePos=9
	
	total=0
	for key in d:
		gene=d[key][1].split(":")[genePos].split(',')[0].split("*")[0]
		if gene in hist:
			hist[gene]+=1
		else:
			hist[gene]=1
		total+=1	

	summary.write('\n\n******  %s  ******\n\n' %(Title))
	for key in sorted(hist, key=lambda i: int(hist[i]), reverse=True):
		percent=round(100.0*hist[key]/total,1)
		summary.write('{:<12} {:<12} {:<12}\n'.format(key,hist[key],percent))	

	
#variable for determining whether the antibody is heavy or light chain 
chain_call = {
	'VDJ':['TRB','TB','IGH','TRD'],
	'VJ':['IGK','TRA','TA','IGL','TRG'],	
}


#this is a simple variable for mapping the fields we require to run the analysis to the field names in different annotation files 
#so when reading files we know that the field name 'Functionality_1' in IMGT files correponds to our functionality variable we use in the analysis 
pairing_settings = {
	'IMGT':{
		'fields_for_analysis':{
			'functionality':'Functionality_1',
			'raw_seq_nt':'Sequence ID_1',
			'vgene':'V-GENE and allele_1',
			'jgene':'J-GENE and allele_1',
			'dgene':'D-GENE and allele_1',
			'raw_seq_nt':'Sequence_1',
			'cdr3_nt':'CDR3-IMGT_3',
			'seq_header':'Sequence ID_1',
			'shm':'VREGION.SHM.NT_PER',			
			'full_len_ab_aa':'PREDICTED_AB_SEQ.AA',#THIS FIELD IS NOT PRESENT IN THE IMGT FILE, BUT IT DOES GET ADDED TO THE INFORMATION WHEN READING IMGT FILES USING OUR READER CLASS			
			'full_len_ab_nt':'PREDICTED_AB_SEQ.NT',#THIS FIELD IS NOT PRESENT IN THE IMGT FILE, BUT IT DOES GET ADDED TO THE INFORMATION WHEN READING IMGT FILES USING OUR READER CLASS			
			'vgene_scores':'V-REGION score_1',
			'dgene_scores':'D-REGION score_1',
			'jgene_scores':'J-REGION score_1',
			'recomb':'R_TYPE'
			
		},
		'productivity_function':IMGTProductivity #when handling IMGT anaysis, use the function described above (IMGTProductivity) to determine whether a sequence is productive or not 
	},
	'IGBLAST':{
		'fields_for_analysis':{
			'functionality':'PRODUCTIVE',
			'raw_seq_nt':"FULL_SEQ",
			'vgene':'VREGION.VGENES',
			'jgene':'JREGION.JGENES',
			'dgene':'DREGION.DGENES',
			'raw_seq_nt':'FULL_SEQ',
			'cdr3_nt':'CDR3.NT',
			'seq_header':'DOCUMENTHEADER',
			'shm':"VREGION.SHM_NT_PER",			
			'full_len_ab_aa':'PREDICTED_AB_SEQ.AA',#optional field but preferred
			'full_len_ab_nt':'PREDICTED_AB_SEQ.NT',			
			'vgene_scores':"VREGION.VGENE_SCORE",
			'jgene_scores':"JREGION.VGENE_SCORE",
			'dgene_scores':"DREGION.VGENE_SCORE",
			'recomb':'RECOMBINATION_TYPE'
		},                    
		'productivity_function':CDR3Productivity
	},
	'IGFFT':{
		'fields_for_analysis':{
			'functionality':'Productive',
			'vgene':'Top_V-Gene_Hits',
			'raw_seq_nt':'Sequence',
			'jgene':'Top_J-Gene_Hits',
			'dgene':'',#Top_D-Gene_Hits
			'raw_seq_nt':'Sequence',
			'cdr3_nt':'CDR3_Sequence.NT',
			'seq_header':'Header',
			'shm':'VRegion.SHM.Per_nt',#'VRegion.SHM.NT',
			'full_len_ab_aa':'Full_Length_Sequence.AA',#optional field but preferred
			'full_len_ab_nt':'Full_Length_Sequence.NT',
			'isotype':'Isotype',
			'vgene_scores':"V-Gene_Alignment_Scores",
			'jgene_scores':"J-Gene_Alignment_Scores",
			'dgene_scores':'',
			'recomb':'Recombination_Type'
		},
		'productivity_function':CDR3Productivity
	},
	'MIXCR':{
		'fields_for_analysis':{
			'functionality':'Productivity',
			'vgene':"All V hits",
			'raw_seq_nt':'Sequence',
			'jgene':"All J hits",
			'dgene':"All D hits",
			'raw_seq_nt':'Sequence',
			'cdr3_nt':'N. Seq. CDR3',
			'seq_header':'Seqheader',
			'shm': 'VGENE: Shm.per',#'VGENE: Shm.nt',
			'full_len_ab_aa':'Full AA',#optional field but preferred
			'full_len_ab_nt':'Full NT',
			'isotype':'All C hits',
			'vgene_scores':"All V scores",
			'jgene_scores':"All J scores",
			'dgene_scores':"All D scores",
			'recomb':'Recombination Type'
		},
		'productivity_function':CDR3Productivity
	},
	
	'DATABASE':{
		'fields_for_analysis':{
			'functionality':'PRODUCTIVE',
			'raw_seq_nt':'SEQUENCE',
			'vgene':'VREGION.VGENES',
			'jgene':'JREGION.JGENES',
			'dgene':'DREGION.DGENES',
			'raw_seq_nt':'PREDICTED_AB_SEQ.NT',
			'cdr3_nt':'CDR3.NT',
			'seq_header':'SEQUENCE_HEADER',
			'shm':'VREGION.SHM.NT_PER',
			'full_len_ab_aa':'PREDICTED_AB_SEQ.AA',#optional field but preferred
			'full_len_ab_nt':'PREDICTED_AB_SEQ.NT',
			'isotype':'ISOTYPE.GENE',
			'vgene_scores':"VREGION.VGENE_SCORES",
			'jgene_scores':"JREGION.JGENE_SCORES",
			'dgene_scores':"DREGION.DGENE_SCORES",
			'recomb':'RECOMBINATION_TYPE'
		},
		'productivity_function':CDR3Productivity
	}
}


supported_analyses = pairing_settings.keys()

def GetHeaderInfo(file_data,header_var):
	h = file_data[header_var]
	
	#check if ididnetifier is in the file first 
	if idIdentifier in file_data and file_data[idIdentifier]:
		id = file_data[idIdentifier]
	else:
		#its not in the file, so it might in header 
		id = ''
		if fasta_file_delimiter in h:
			#check header 
			tmp = h.split(fasta_file_delimiter)
			try:
				additional_data = json.loads(tmp[-1])
			except:
				additional_data={}
			if idIdentifier in additional_data:
				id = additional_data[idIdentifier]
	return [h,id]

def initialize_input_files(analysis_method,list_of_filetypes,list_of_files):
	global annotation_headers
	#load all of the filenames and open the files using our class for reading different files 
	if list_of_filetypes == None:
		list_of_file_reading = [readfile.immunogrepFile(file) for i,file in enumerate(list_of_files)]
	elif list_of_filetypes[0] =='IMGT':
		list_of_file_reading = [readfile.immunogrepFile(file,'IMGT',required_files=[1,3,8],field_names=['Sequence_1','V-D-J-REGION_3','V-D-J-REGION_5','V-J-REGION_3','V-J-REGION_5','CDR3-IMGT_5','Functionality_1','V-GENE and allele_1','D-GENE and allele_1','J-GENE and allele_1','CDR3-IMGT_3','Sequence ID_1','V-REGION Nb of mutations_8']) for i,file in enumerate(list_of_files)]
	else:
		list_of_file_reading = [readfile.immunogrepFile(file,list_of_filetypes[i]) for i,file in enumerate(list_of_files)]
		
	return list_of_file_reading


def ProcessGene(gene):
	"""
		Removes the allele calls from genes
		
		Assumption: 
			
			Multiple genes are seperated by ',' 
			
			Alleles are seperated by '*' 
			
			If a gene is seperated by multiple spaces, then the gene should be identified by a gene that contains either - or '*' 
			For example: 
				Imgt genes may be: 
					Homo sapiens IGHV1-3*01
					We only want to isolate the word IGHV1-3
	"""
	
	if not gene:
		return ''	
	#extract all genes in field 
	gene_array = gene.split(',')
	for gene_num,each_gene in enumerate(gene_array):
		#split each gene by spaces. Go through each word
		split_words  = each_gene.split(' ')
		if len(split_words)==1: #there is only one word 
			gene_array[gene_num] = split_words[0].split('*')[0]
		else:
			g = ''
			for subv in split_words:
				#if we find a word with gene characters in it, it must be our gene 
				if '*' in subv or '-' in subv:						
					#extract everything before '*'
					g = subv.split('*')[0]
					break
			gene_array[gene_num] = g	
	return ','.join(gene_array)
	
def get_top_genes(gene_string,gene_scores):
	"""

		Return a list of 'top' unique genes
		Go through the reported list of genes and only report genes whose scores are equal to the top score
		
		Assumptions:	
			We err on the side of consider 'more genes' than less
			1) If there is no score provided for a field, then consider all genes equally 
			2) If only a single gene score is provided, then again, assume all reported genes have equal alignment (i.e. IMGT reporting)
			3) For all other cases where there is more than one gene score provided, then traverse the scores and only report genes whose score = top score
	
	"""
	if not gene_string:
		return []
	#gene_string should only be genes not alleles. this should have been handled in the add_to_dict_memory_safe function 
	gene_array = gene_string.split(',')
	
	#no scores provided, so we assume they are all equal 
	if not gene_scores:
		return list(set(gene_array))
	else:
		try:			
			gene_numbers = [float(g) for g in gene_scores.split(',')]
		except:
			#gene scores are not numbers
			return list(set(gene_array))

		if len(gene_numbers) == 1:
			#IMGT will report only one gene score for multiple genes. we assume these genes are idetnically scored then 
			return list(set(gene_array))
		else:
			top_score = gene_numbers[0]
			top_gene_list = []
			for n,s in enumerate(gene_numbers):
				if s<top_score:
					break
				top_gene_list.append(gene_array[n])
			return list(set(top_gene_list))
						
def parse_sorted_paired_file(pairing_temp_file,suffix):
	"""
	
		Function for collapsing by CDR3 sequence
		We assume that the input file (pairing_temp_file) is a tab delimited file whose columns match the variable ab_field_loc		
		
	
	"""
	ab_field_loc = get_ab_field_loc()
	global codeDict_output_dict
	global dict_summary
	global collapsed_output_file
	global usearch_cluster_file
	
	#parent_folder = '/'.join(dict_summary.split('/')[:-1])+'/'
	#temp_barcode_file = parent_folder+'cdrh3_l3_barcodes.txt'
	
	my_folder = os.path.dirname(pairing_temp_file)# '/'.join(pairing_temp_file.split('/')[:-1])+'/'

	print('Will now sort file by CDR3 sequence')
	#fullcode_column = num_elem_stored+1#last column number in file should be equal to the number of elements stored (in ab_field_loc) +1
	#cdr3_column = ab_field_loc['CDR3']+1
	#sort all sequences sent to file by the full code field (the last column in file) 		
	#subprocess.call('''sort -T "{2}" -t '\t' "{0}" -k{1}>"{0}.sorted" '''.format(pairing_temp_file,str(cdr3_column),my_folder), shell=True)		
	#os.remove(pairing_temp_file)
	#os.rename(pairing_temp_file+'.sorted',pairing_temp_file)
		
	mapping_dict = {}	
	r_count = 0
	index_loc = {'VDJ':0,'VJ':1}	

	dict_summary+='.'+suffix+'.txt'
	collapsed_output_file+='.'+suffix+'.txt'
	usearch_cluster_file+='.'+suffix+'.txt'
	#codeDict_save=open(codeDict_output_dict,'w')
	write_dict_summary=open(dict_summary,'w')
	write_dict_summary_unsorted=open(dict_summary+'.unsorted.txt','w')
	#collapsed_cdr3_barcode_file = open(temp_barcode_file,'w')	
	cdrh3_l3_barcode_dict = defaultdict(int)
	#codeDict_save.write('{')#we will save the variable codeDict as a JSON var to file. JSON files start with { and end with }. ',' will be used to save each R1-R2 pair
	
	#read through current filtered file and store counts of unique CDRH3. We dont
	#also create the 'line-by-line' UNSORTED and line-by-line SORTED files 
	with open(pairing_temp_file) as unsorted_file:
		for line in unsorted_file:
			line=line.strip('\r\n').split('\t')
			fullcode=line[-1]			
			if not fullcode:
				#shouldnt happen, but check anyway...
				continue				
			temp_array = line			
			rtype = temp_array[ab_field_loc['CHAIN_CALL']]			
			pos = index_loc[rtype]
			array = [[],[],'']
			array[pos] = {field:temp_array[index_val] for field,index_val in ab_field_loc.iteritems()}
			barcode_cdr3 =temp_array[ab_field_loc['CDR3_SEQ']]
			write_line = json.dumps(array)
			write_dict_summary.write(barcode_cdr3+'\t'+write_line+'\n')						 
			write_dict_summary_unsorted.write(write_line+'\n')						
			cdrh3_l3_barcode_dict[barcode_cdr3]+=1							
												
	write_dict_summary.close()
	write_dict_summary_unsorted.close()
	
	#Now sort the write_dict_summary file by CDRH3-L3 barcode sequences
	print('Sorting paired sequences by unique CDR3 sequences')
	subprocess.call('''sort -T "{1}" -t '\t' "{0}" -k1,1 > "{0}.sorted" '''.format(os.path.abspath(dict_summary),os.path.abspath(my_folder)), shell=True)
	#Now open the file again and keep track of unique CDRH3-CDRL3 sequence data. Keep and updated list of SHM to calculate the averate and standard deviation. 
	#Also keep track of which line has the MOST data. i.e. has all vgenes and isotype information. By default we will store this information as the 'unique' sequence. 
	print('Collapsing sequences by unique CDR3. Finding average SHM and most occurring "Gene Mode" for unique sequences')
	group_barcode = ''
	barcode_count = 0
	summed_counts = 0 
	collapsed_output=open(collapsed_output_file+'.presorted','w')
	
	old_summary = open(dict_summary+'.sorted')
	new_dict_summary = open(dict_summary,'w')

	l = 0
	while True:
		try:
			each_seq_line = old_summary.readline()
		except:
			break
		if not each_seq_line:
			break
		l+=1
		each_seq_line=each_seq_line.strip('\r\n')
		each_seq_line = each_seq_line.split('\t')
		barcode_cdr_h3_l3 = each_seq_line[0]
		group_barcode = barcode_cdr_h3_l3
		num_seqs = cdrh3_l3_barcode_dict[barcode_cdr_h3_l3]
		current_seq = 1
		remaining_data = '\t'.join(each_seq_line[1:])
		new_dict_summary.write(str(num_seqs)+'\t'+remaining_data+'\n')							
		
		seed_barcode_data = json.loads(remaining_data)
		vdj_data = seed_barcode_data[index_loc['VDJ']]
		vj_data = seed_barcode_data[index_loc['VJ']]
		
		if not vdj_data and not vj_data:
			continue
		elif vdj_data:
			data = vdj_data
		elif vj_data:
			data = vj_data
		vdj_data = None
		vj_data = None

		seed_barcode_data = data

		num_not_empty = 0
		for t in data:
			if t and t!="none" and t!="N/A":
				num_not_empty+=1
						
		most_fields = num_not_empty
		shm_data = [0,0,0]
		
		mut_found=False
		
		if data['MUT']!='none':
			shm_val = float(data['MUT'])
			mut_found=True
			shm_data = [1,shm_val,pow(shm_val,2)]		
		
		vgene_counts = defaultdict(int)
		jgene_counts = defaultdict(int)
		dgene_counts = defaultdict(int)
								
		#for all of the V(D)J genes in the pair, determine which genes are the top gene hits. return these genes as a list and update the current gene counts to the variable 
		#at the end of the group, we will determine the mode/most likely gene. for now, we only chose the first instance of a top count 
		for top_genes in get_top_genes(data['VGENE'],data['VSCORES']):			
			vgene_counts[top_genes]+=1
		for top_genes in get_top_genes(data['DGENE'],data['DSCORES']):		
			dgene_counts[top_genes]+=1
		for top_genes in get_top_genes(data['JGENE'],data['JSCORES']):
			jgene_counts[top_genes]+=1
				
		#go through all remaining sequences in this file that contains the same CDRH3-CDRL3 PAIR and take the sum of the SHM and SHM^2 												
		while current_seq!=num_seqs:		
			current_seq+=1			
			each_seq_line = old_summary.readline()
			each_seq_line=each_seq_line.strip('\r\n')
			each_seq_line = each_seq_line.split('\t')
			barcode_cdr_h3_l3 = each_seq_line[0]
			if barcode_cdr_h3_l3!=group_barcode:
				print 'error matching barocdes'
				print group_barcode
				print barcode_cdr_h3_l3
				raise Exception('problem collapsing by unique cdrh3-cdrl3')																
			
			remaining_data = '\t'.join(each_seq_line[1:])
			new_dict_summary.write(str(num_seqs)+'\t'+remaining_data+'\n')							
			
			barcode_data = json.loads(remaining_data)
			vdj_data = barcode_data[index_loc['VDJ']]
			vj_data = barcode_data[index_loc['VJ']]
			
		
			if not vdj_data and not vj_data:
				continue
			elif vdj_data:
				data = vdj_data
			elif vj_data:
				data = vj_data
			vdj_data = None
			vj_data = None


			num_not_empty = 0
			for t in data:
				if t and t!="none" and t!="N/A":
					num_not_empty+=1
			
			if num_not_empty>most_fields:
				#if you find more non empty fields in this sequence, then it should be treated as the 'seed'
				seed_barcode_data = data
				most_fields = num_not_empty
			
			if data['MUT']!='none':
				shm_val = float(data['MUT'])
				mut_found=True					
				shm_data[0]+=1
				shm_data[1]+=shm_val
				shm_data[2]+=pow(shm_val,2)
			
			#for all of the V(D)J genes in the pair, determine which genes are the top gene hits. return these genes as a list and update the current gene counts to the variable 
			#at the end of the group, we will determine the mode/most likely gene. for now, we only chose the first instance of a top count 			
			for top_genes in get_top_genes(data['VGENE'],data['VSCORES']):			
				vgene_counts[top_genes]+=1
			for top_genes in get_top_genes(data['DGENE'],data['DSCORES']):		
				dgene_counts[top_genes]+=1
			for top_genes in get_top_genes(data['JGENE'],data['JSCORES']):
				jgene_counts[top_genes]+=1
							
		if mut_found:
			mut_sum = shm_data[1]
			mut_sum_sq = round(shm_data[2],6)
			mut_count = shm_data[0]				
			sum_sq = round(pow(mut_sum,2)/mut_count,6)		
			#average shm
			mut_avg = round(mut_sum/mut_count,3)
			#shm variance -> sum of squares formulat => sum(vals^2)-(sum(vals)^2/counts)
			if mut_count>1:
				#mut_var = round(pow((mut_sum_sq-(pow(mut_sum,2)/mut_count))/(mut_count-1),0.5),3)
				mut_var = round(pow((mut_sum_sq-sum_sq)/(mut_count-1),0.5),3)
			else:
				mut_var = 'none'
		else:
			mut_avg = 'none'
			mut_var = 'none'
			mut_sum = 0
			mut_sum_sq = 0
			mut_count = 0					
		
		data = seed_barcode_data
		cdrlen = data['CDR3_LEN']

						
		summed_counts+=num_seqs
		
		##store the FIRST gene matching MAX MODE of the V(D)J genes 
		#vdj_data['VGENE'] = max(vgene_counts_vdj, key=vgene_counts_vdj.get)
		#vdj_data['DGENE'] = max(dgene_counts_vdj, key=dgene_counts_vdj.get)
		#vdj_data['JGENE'] = max(jgene_counts_vdj, key=jgene_counts_vdj.get)
		#vj_data['VGENE'] = max(vgene_counts_vj, key=vgene_counts_vj.get)		
		#vj_data['JGENE'] = max(jgene_counts_vj, key=jgene_counts_vj.get)

		##store the ALL genes matching MAX MODE of the V(D)J genes 
		if vgene_counts:
			max_item = max(vgene_counts.values())
			data['VGENE'] = ','.join([g for g,num in vgene_counts.iteritems() if num==max_item]) if max_item > 0 else ''
		else:
			data['VGENE'] = ''
		if dgene_counts:
			max_item = max(dgene_counts.values())
			data['DGENE'] = ','.join([g for g,num in dgene_counts.iteritems() if num==max_item]) if max_item > 0 else ''
		else:
			data['DGENE'] = ''
	
		if jgene_counts:
			max_item = max(jgene_counts.values())
			data['JGENE'] = ','.join([g for g,num in jgene_counts.iteritems() if num==max_item]) if max_item > 0 else ''
		else:
			data['JGENE'] = ''
		
		if data['CHAIN_CALL'].upper()=='VDJ':
			#the variable ab_field_loc stores the index position for each antibody region in the array 
			CDRH3=data['CDR3_SEQ']
			VHgene=data['VGENE']
			DHgene=data['DGENE']
			JHgene=data['JGENE']
			CDRL3=''
			VLgene=''
			JLgene=''
			IgH=data['LOCUS']
			IgL=''
			h_iso = data['ISOTYPE']
			l_iso = ''
			Hlen=cdrlen
			Llen=''		
			try: 
				CDRH3_trans=str(Seq(CDRH3,generic_dna).translate())
			except Exception as e: 
				CDRH3= 'translation error'
				print('Error translating: '+str(e))
			CDRL3_trans=''

			Hmut_avg = mut_avg
			Hmut_var = mut_var
			Lmut_avg = ''
			Lmut_var = ''
			Hmut_sum=mut_sum
			Hmut_sum_sq=mut_sum_sq
			Hmut_count=mut_count
			Lmut_sum=''
			Lmut_sum_sq=''
			Lmut_count = ''

		elif data['CHAIN_CALL'].upper()=='VJ':
			#the variable ab_field_loc stores the index position for each antibody region in the array 
			CDRH3=''
			VHgene=''
			DHgene=''
			JHgene=''
			CDRL3=data['CDR3_SEQ']
			VLgene=data['VGENE']
			JLgene=data['JGENE']
			IgH=''
			IgL=data['LOCUS']
			h_iso = ''
			l_iso = data['LOCUS']
			Hlen = ''
			Llen = cdrlen

			Hmut_avg = ''
			Hmut_var = ''
			Lmut_avg = mut_avg
			Lmut_var = mut_var
			Hmut_sum=''
			Hmut_sum_sq=''
			Hmut_count=''
			Lmut_sum=mut_sum
			Lmut_sum_sq=mut_sum
			Lmut_count = mut_count
						
			try: 
				CDRL3_trans=str(Seq(CDRL3,generic_dna).translate())
			except Exception as e: 
				CDRL3= 'translation error'
				print('Error translating: '+str(e))
			CDRH3_trans=''
			

		results = [num_seqs,CDRH3_trans,CDRL3_trans,CDRH3,CDRL3,VHgene,DHgene,JHgene,VLgene,JLgene,IgH,IgL,Hmut_avg,Hmut_var,Lmut_avg,Lmut_var,Hlen,Llen,h_iso,l_iso,Hmut_sum,Hmut_sum_sq,Hmut_count,Lmut_sum,Lmut_sum_sq,Lmut_count]
		results = [str(r) for r in results]
			
		collapsed_output.write("\t".join(results)+'\n')
			
	old_summary.close()
	new_dict_summary.close()	
	
	collapsed_output.close()
	
	print('Compiling collapsed reads file: sorting unique CDRH3-CDRL3 sequences by their respective counts and exporting to FASTA file')
	
	subprocess.call('''sort -T "{1}" -t '\t' "{0}.presorted" -k1,1nr >"{0}" '''.format(collapsed_output_file,my_folder), shell=True)
	
	#convert output to fasta file using AWK, ONLY print counts ($1) above 1 
	
	if suffix.upper()=='VDJ':
		fasta_seq = "$4" #$4 => nucleotide CDRH3
	elif suffix.upper()=='VJ':		
		fasta_seq = "$5" #$4 => nucleotide CDRL3
	else:
		raise Exception('Error providing suffix field to function: '+suffix)


	awk_command='gawk -v offile="'+usearch_cluster_file
	awk_command+='''" 'BEGIN{OFS=":";FS="\t"};
					(int($1)>1&&'''+fasta_seq+'''!=""){s=$1;for(i=2;i<=NF;i++)s=s""OFS""$i; print">"s>offile;print '''+fasta_seq+'''>offile;}' '''+collapsed_output_file
	
	subprocess.call(awk_command,shell=True)
	
	#if there are no VH-VL pairs found, then the awk command above will not create an empty file.
	#so to account for that we will do a hack and create an empty file 
	if not os.path.isfile(usearch_cluster_file):
		with open(usearch_cluster_file,'w') as w:
			pass
	
	os.remove(dict_summary+'.sorted')
	os.remove(collapsed_output_file+'.presorted')
	
	print('Creating a sorted line-by-line cluster file')
	os.rename(dict_summary,dict_summary+'.presorted')	
	final_sort = '''sort -T "{1}" -t '\t' "{0}.presorted" -k1,1nr | gawk 'BEGIN{{FS="\t"}};{{print $2 > "{0}"}}' '''.format(dict_summary,my_folder)	
	subprocess.call(final_sort,shell=True)	
	os.remove(dict_summary+'.presorted')
	
	num_unique_h_l = len(cdrh3_l3_barcode_dict)
	cdrh3_l3_barcode_dict = {}
	del  cdrh3_l3_barcode_dict
	gc.collect()
	
	
	return [num_unique_h_l]
	

#THIS IS another version of add_to_dict. Its functionality should be identical except it does not store everything inside of memory 
#The previous ADD_TO_DICT FUNCTION is in the main function after this one !!!
#THIS FUNCTION DIFFERS FROM ADD_TO_DICT BELOW IN THAT IT WILL TRY TO MINIMIZE THE AMOUNT OF MEMORY REQUIRED FOR RUNNING ANALYSIS (it uses about 1/3 to 1/4th of the memory required by  add_to_dict
#initialitests show that speed is not affected using this method. That is, its not slower than the previous add_to_dict method 
#THIS FUNCTION WILL READ THROUGH THE LIST OF FILES PROVIDED
#SO INSTEAD OF STORING DATA FROM READS IN A DICTIONARY, IT WILL INSTEAD STORE THE RESULTS IN a single tab delimited file 
#this file will be parsed by the function parse_sorted_paired_file using a mixture of awk and python code (see parse_sorted_paired_file) function 
def add_to_dict_memory_safe(list_of_files,list_of_filetypes,required_field_names,productivity_function_call,parent_dir,analysis_method,include_vl_analysis):	
	##DEFINING SOME VARIABLES##

	#get global files names for opening files to write/use in this function 
	global dict_summary
	global usearch_cluster_file
	
	#set global variables for field names we require to run analysis 
	global vgene_field
	global jgene_field
	global dgene_field
	global vgene_scores
	global jgene_scores
	global dgene_scores
	global recomb_field
	global seq_field
	global cdr3_field
	global header_field
	global mut_field
	global functionality_field
	global full_aa_field
	global full_nt_field
	
	required_field_names = defaultdict(str,required_field_names) #convert this to default dict so that not all fields have to be explicitly defined
	
	#keys to that refer to field names. these are the key values we use for reading files 
	vgene_field = required_field_names['vgene']
	jgene_field = required_field_names['jgene']
	dgene_field = required_field_names['dgene']
	seq_field = required_field_names['raw_seq_nt']
	cdr3_field = required_field_names['cdr3_nt']
	isotype_field = required_field_names['isotype']
	header_field = required_field_names['seq_header']
	mut_field = required_field_names['shm']
	functionality_field = required_field_names['functionality']
	full_aa_field = required_field_names['full_len_ab_aa']
	full_nt_field = required_field_names['full_len_ab_nt']
	vgene_scores = required_field_names['vgene_scores']
	jgene_scores = required_field_names['jgene_scores'] 
	dgene_scores = required_field_names['dgene_scores']
	recomb_field = required_field_names['recomb']
	
	num_vdj = 0 
	num_vj = 0 

	#seq_field = required_field_names['sequence']
	
	if productivity_function_call == None:
		productivity_function_call = CDR3Productivity	#GeneralProductivity
										
	#write_dict_summary = open(dict_summary,'w')
	#collapsed_save=open(collapsed_output_dict,'w')	
	#collapsed_output=open(collapsed_output_file,'wb')
	#error=open(error_file,'wb')	
	
	pairing_temp_file =os.path.join(os.path.dirname(dict_summary),'temp_annotation_files') #'/'.join(dict_summary.split('/')[:-1])+'/temp_annotation_files.txt'
	temp_output_files_vdj  = open(pairing_temp_file+'.vdj.txt','w')		
	temp_output_files_vj  = open(pairing_temp_file+'.vj.txt','w') if include_vl_analysis else None
	
	#temp_seq_data = open(pairing_temp_file,'w')	
		
	counter=0
		
	codeDictSeqs = {} #stores whether or not a sequence was successfully found (i.e. passes all filters) => 0 = not successful, 1 = successful 
	collapsed_dict = defaultdict(int) #store counts for unique cdrH3/CDRL3 counts 
	
	nonpaired_len = 0	
	no_cdr3_error = 0 #keep track of sequences which lack cdr3
	invalid_chain = 0 #keep track of sequences which do not have a valid chain call based on the chain_call variable above 
	not_productive = 0 #keep track of sequences which are determined to be 'unproductive'
	total_seqs = 0 #keep track of all sequences
	no_result = 0 #keep track of sequences that have no antibody information 	
	passed_filter = 0 #keep track of sequences that pass the filters described above (no cdr3, unproductive, etc)	
	
	ab_field_loc = get_ab_field_loc()
	
	###VARIABLES DEFINED, CONTINUING WITH FUNCTION###
	
	print 'Reading all files at once and collapsing identical CDRH3-CDRL3 pairs into collapsed dict' 
	
	at_least_one_file_open = True
	
	#we will generate annotation files as we read through the files for pairing 
	list_of_file_reading = initialize_input_files(analysis_method,list_of_filetypes,list_of_files)
	# we will read each of the files simultaneously (open all the files at once and read line by line)
	while at_least_one_file_open:		
	#while counter<200000:
		counter+=1
		if total_seqs%100000==0:
			print 'Processed '+str(total_seqs)+' sequences'#: '+str(counter)			
		
		#keep track of how many files have been completely read
		num_eof = 0 
		
		#read each file line by line
		for fnum,reader in enumerate(list_of_file_reading): 
			if reader.IFclass.eof:
				num_eof+=1
				continue						
				
			my_line = reader.IFclass.read()
							
			if not my_line:#probably end of file or just empty line			
				continue
			
			my_line = defaultdict(str,my_line)
			h = my_line[header_field]				
			if not h:								
				print 'Error type 0: no sequence header'
				continue
			
			[header,id] = GetHeaderInfo(my_line,header_field)						
			seq = my_line[seq_field]
			#get miseq code/read info 
			fullcode = ':'.join(header.replace(' ','_').split(':')[3:7]).split('_')[0]															
						
			total_seqs+=1												
			
			#no amino acid sequence found 
			if not my_line[full_aa_field]:			
				
				#we no longer require the variable full_aa_field to be 'non empty'
				if not my_line[full_nt_field]:					
					#no nucleotide sequence found either 
					#do not save this sequences results to file 
					#write to temperoary file annotation file, store that there is no full length sequence for this sequence, so no need to save its 
					#paired annotation information 
					#annotation_file_writing[fnum]['buffer'].write('\t'.join([header,seq,id,'','',''])+'\n')						
					#no_result+=1
					#continue
					pass
				else:
					#aminoa acid sequence was not proviided, but a nucleotide sequence was provided 
					#translate nucloetid to amino acid..this can result in some potential problems with PRODUCTIVITY determination IF using GENERAL PRODUCTIVTY function rule 
					try:
						end = (len(my_line[full_nt_field])/3)*3
						my_line[full_aa_field] = str(Seq(my_line[full_nt_field][:end],generic_dna).translate())
					except:
						print('Error type 2: problem translating amino acid => '+str(my_line[full_nt_field]))
						#annotation_file_writing[fnum]['buffer'].write('\t'.join([header,seq,id,'','',''])+'\n')						
						#no_result+=1
						#continue
						my_line[full_aa_field] = ''
						#no cdr3 found 
			if not(my_line[cdr3_field]):											
				#write to temperoary file 				
				no_cdr3_error+=1
				continue				
									
			#maintain all vgenes as a list 
			#Process genes and only extract the GENE from each of the fields (that is ignore the alleles '*01')
			Vgene = ProcessGene(my_line[vgene_field])#.split(',')#[0]
			Jgene = ProcessGene(my_line[jgene_field])#.split(',')#[0]			
			Dgene = ProcessGene(my_line[dgene_field])#.split(',')#[0]				
			#Unlike the above change for genes, only keep the first isotype in the list 
			isotype = my_line[isotype_field].strip().split(',')[0]			
			#figure out locus using vgene call 														
		
			if Vgene:
				#use vgene to determine locus 
				locus = Vgene.split(',')[0][:3].upper()						
			elif Jgene:
				#if no vgene is present, then use jgene
				locus = Jgene.split(',')[0][:3].upper()			
			else:
				#ignore using D gene
				locus=''			
			
			chain = 'N/A'			
		
			#the provided file defines recombination/VDJ or VJ chain 
			if my_line[recomb_field] and my_line[recomb_field] in ['VDJ','VJ']:
				chain = my_line[recomb_field]
			elif locus:				
				#shoudl return 'VDJ' or 'VJ' based on LOCUS var
				for possible_chains,values in chain_call.iteritems():			
					if locus in values:
						chain = possible_chains.upper()
						break		
			
															
			#chain could not be determind with provided locus 
			if chain == 'N/A':		
				this_chain=''
				invalid_chain+=1
				print 'Error type 1' #isotype is mislabeled
				print my_line			
				continue
			
			#this is the recombination type of the current sequence 
			rtype=chain.upper()
									
			#unproductive
			if productivity_function_call(my_line)==False:																									
				#write to temperoary file 
				not_productive+=1
				continue
								
			CDR3 = my_line[cdr3_field].upper()
			Len=len(CDR3)/3	
						
			#get SHM info 
			if my_line[mut_field]:
				try:
					Mut = round(float(my_line[mut_field]),3)
				except:
					if '(' in my_line[mut_field]:
						try:
							Mut = round(float(my_line[mut_field].split('(')[0].strip()),3)
						except:
							Mut = 'none'
					else:				
						Mut = 'none'		
			else:
				Mut = 'none'
			
			#get isotype 
			isotype = my_line[isotype_field]			
			seq=my_line[seq_field]						
			v_scores = my_line[vgene_scores]
			d_scores = my_line[dgene_scores]
			j_scores = my_line[jgene_scores]
			
			passed_filter+=1
						
			#This specific sequence passed all filters. Write this sequence to a file that can be sorted afterwards by its 'full code' afterwards
			temp_array = ['']*(num_elem_stored)			
			temp_array[ab_field_loc['CDR3_SEQ']] = CDR3
			temp_array[ab_field_loc['VGENE']]=Vgene
			temp_array[ab_field_loc['DGENE']]=Dgene
			temp_array[ab_field_loc['JGENE']]=Jgene
			temp_array[ab_field_loc['LOCUS']]=locus
			temp_array[ab_field_loc['MUT']]=Mut
			temp_array[ab_field_loc['CDR3_LEN']]=Len
			temp_array[ab_field_loc['AB_SEQ']]=my_line[full_nt_field]
			temp_array[ab_field_loc['CHAIN_CALL']]=chain 		
			temp_array[ab_field_loc['ISOTYPE']] = isotype
			temp_array[ab_field_loc['VSCORES']] = v_scores
			temp_array[ab_field_loc['DSCORES']] = d_scores
			temp_array[ab_field_loc['JSCORES']] = j_scores			
			temp_array.append(fullcode)			
			
			#now store all of the importat fields we want to use for pairing to the temp file defined by temp_seq_data
			if rtype=='VDJ':
				num_vdj +=1
				temp_output_files_vdj.write('\t'.join([str(t) for t in temp_array])+'\n')
			elif rtype=='VJ':
				num_vj+=1
				if  include_vl_analysis:
					temp_output_files_vj.write('\t'.join([str(t) for t in temp_array])+'\n')						
				
		if num_eof==len(list_of_files):#all files have been read through
			at_least_one_file_open=False 	
	
	#temp_seq_data.close()
	
	temp_output_files_vdj.close()
	if temp_output_files_vj:
		temp_output_files_vj.close()
	
			
	#ok the temp file was made. it only consits of sequences which passed the filters above. We now need to sort this file by the MISEQ id and group sequences into proper VH-VL pairs 
	[num_unique_cdr3] = parse_sorted_paired_file(pairing_temp_file+'.vdj.txt','vdj')
	if os.path.isfile(usearch_cluster_file):		
		num_unique_cdr3_pair_above1 = useful.file_line_count(usearch_cluster_file)/2
	else:
		num_unique_cdr3_pair_above1 = 0
	os.remove(pairing_temp_file+'.vdj.txt')
	#os.remove(pairing_temp_file)				
				
	print 'Summary: '
	print 'Parsed through {0} sequences'.format(str(total_seqs))
	print '{0} ({1}%) sequences did not pass filters: '.format(str(total_seqs-passed_filter),str( round(float(  (100*(total_seqs-passed_filter))/total_seqs),1)) if total_seqs>0 else '0' )
	#print '		{0} sequences did not have an antibody sequence'.format(str(no_result))
	print '		{0} sequences did not have a cdr3'.format(str(no_cdr3_error))
	print '		{0} sequences were not productive'.format(str(not_productive))
	print '		{0} sequences had an unidentifiable chain type'.format(str(invalid_chain))

	print '{0} VDJ sequences'.format(str(num_vdj))
	print '{0} VJ sequences'.format(str(num_vj))
		
	print 'CDRH3 sequences were collapsed into {0} sequences containing unique CDRH3 sequences: '.format(str(num_unique_cdr3))
	print '{0} ({1}%) of these unique pairs were observed more than once'.format(num_unique_cdr3_pair_above1,str(round(float(100*num_unique_cdr3_pair_above1/num_unique_cdr3),1)) if num_unique_cdr3>0 else '0')
		   
	return [total_seqs,passed_filter,no_result,no_cdr3_error,not_productive,invalid_chain,num_unique_cdr3,num_unique_cdr3_pair_above1,num_vdj,num_vj]
	
def Compile_Clusters(clustered_dict_raw,cluster_val):	
	"""
	
		Generates a compiled output file summarizing which CDRH3 clusters were made. 
		Reports the cluster size and extra features such as diversity and dominance to determine whether clustering makes sense for H-L pair 
	
	"""
	global clustered_final_output_file
	
	clustered_final_output=open(clustered_final_output_file+'.'+str(cluster_val)+'.txt','wb')		
	header_output = ['Total Counts','Indiv_Count','Dominance', 'Confidence', 'Diversity (Shannon)','Diversity (Simpsons)', 'Cluster size', 'CDRH3 AA', 'CDRL3 AA', 'CDRH3','CDRL3','VHgene','DHgene','JHgene','VLgene','JLgene','IgH','IgL','Hmut_avg','Hmut_std','L_mut_avg','Lmut_std','HLen','LLen','Hisotype','Lisotype']
	clustered_final_output.write("\t".join(header_output)+'\n')
	
	for key in sorted(clustered_dict_raw, key=clustered_dict_raw.get, reverse=True):	
		best_hit = clustered_dict_raw[key][1].split(':')
		#total_counts found from all CDR3s in a single cluster 
		total_counts=float(clustered_dict_raw[key][0])
				
		indiv_count=float(best_hit[0]) #total count for the BESt/highest populated CDR3 in the cluster/largest size
		CDRH3_trans=best_hit[1]
		CDRL3_trans=best_hit[2]
		CDRH3=best_hit[3]
		CDRL3=best_hit[4]
		VHgene=best_hit[5]
		DHgene=best_hit[6]
		JHgene=best_hit[7]		
		VLgene=best_hit[8]
		JLgene=best_hit[9]
		IgH=best_hit[10]
		IgL=best_hit[11]
						
		
		HLen=best_hit[16]
		LLen=best_hit[17]
		
		HIso = best_hit[18]
		LIso = best_hit[19]
		
		if best_hit[20]:
			Hmut_sum = float(best_hit[20])
		else:
			Hmut_sum = 0 
		if best_hit[21]:
			Hmut_sum_sq = float(best_hit[21])
		else:
			Hmut_sum_sq=0
		if best_hit[22]:
			Hmut_count = float(best_hit[22])
		else:
			Hmut_count = 0  
		
		if best_hit[23]:
			Lmut_sum = float(best_hit[23])
		else:
			Lmut_sum = 0
		if best_hit[24]:
			Lmut_sum_sq = float(best_hit[24])
		else:
			Lmut_sum_sq = 0
		
		if best_hit[25]:
			Lmut_count = float(best_hit[25])
		else:
			Lmut_count = 0 
						
		dominance=int(indiv_count*100/total_counts) #what percent of the TOP cdr3 contributes to the TOTAL count 
		
		freq = indiv_count/total_counts
		diversity = freq*math.log(freq)
		diversity_v2 = pow(freq,2)
		
		cluster_size = len(clustered_dict_raw[key])-1
		
		if cluster_size>1:
			indiv_count2=float(clustered_dict_raw[key][2].split(":")[0]) #this is the counts for the second best CDR3
			confidence=int(100-100*indiv_count2/indiv_count) #confidence => ratio of best cdr3 count/2nd best cdr3 count									
			for k in range(2,len(clustered_dict_raw[key])):
				other_hits = clustered_dict_raw[key][k].split(':')
				if other_hits[20]:
					Hmut_sum+=float(other_hits[20])
				if other_hits[21]:
					Hmut_sum_sq+=float(other_hits[21])
				if other_hits[22]:
					Hmut_count+=float(other_hits[22])
				if other_hits[23]:
					Lmut_sum+=float(other_hits[23])
				if other_hits[24]:
					Lmut_sum_sq+=float(other_hits[24])
				if other_hits[25]:
					Lmut_count+=float(other_hits[25])
				freq = float(other_hits[0])/total_counts
				diversity+=freq*math.log(freq)
				diversity_v2+=pow(freq,2)
			diversity = math.exp(-1*diversity)
			diversity_v2 = 1/diversity_v2			
			if Hmut_count>0:
				Hmut_avg = round(Hmut_sum/Hmut_count,3)
			else:
				Hmut_avg = 'none'			
			if Hmut_count>1:
				h1 = round(Hmut_sum_sq,6)
				h2 = round(pow(Hmut_sum,2)/Hmut_count,6)				
				Hmut_var =  round(pow( (h1-h2)/(Hmut_count-1),0.5),3)
			else:
				Hmut_var = 'none'			
			if Lmut_count>0:
				Lmut_avg = round(Lmut_sum/Lmut_count,3)			
			else:
				Lmut_avg = 'none'				
			#shm variance -> sum of squares formulat => sum(vals^2)-(sum(vals)^2)/counts)
			if Lmut_count>1:
				l1 = round(Lmut_sum_sq,6)
				l2 = round(pow(Lmut_sum,2)/Lmut_count,6)				
				Lmut_var =  round(pow( (l1-l2)/(Lmut_count-1),0.5),3)				
			else:
				Lmut_var = 'none'			
		else:				
			Hmut_avg = best_hit[12]
			Hmut_var = best_hit[13]
			Lmut_avg = best_hit[14]
			Lmut_var = best_hit[15]
			confidence=100
			diversity = 1
			diversity_v2 = 1
						
		try: CDRH3_trans=(Seq(CDRH3,generic_dna)).translate()
		except: CDRH3= 'translation error'
		try: CDRL3_trans=(Seq(CDRL3,generic_dna)).translate()
		except: CDRL3= 'translation error'
		
		results = [total_counts,indiv_count,dominance,confidence,diversity,diversity_v2,cluster_size,CDRH3_trans,CDRL3_trans,CDRH3,CDRL3,VHgene,DHgene,JHgene,VLgene,JLgene,IgH,IgL,Hmut_avg,Hmut_var,Lmut_avg,Lmut_var,HLen,LLen,HIso, LIso]
		results = [str(r) for r in results]		
		clustered_final_output.write("\t".join(results)+'\n')	

				
def WriteSummaryFile(annotated_file_paths,total_seqs,passed_filter,no_result,no_cdr3_error,not_productive,invalid_chain,num_unique_cdrh3_l3_pair,num_unique_cdrh3_l3_pair_above1,cluster,num_clustered,num_vdj,num_vj):
	#MAKE THIS SPECIAL IF ITS IMGT FILES (MERGE FILES INTO ONE GROUP NAME)	
	file_path_summary = ['\t'+f+',' for f in annotated_file_paths]

	global summary_file
	summary = open(	summary_file,'wb')
	#Write useful information in summary file	
	summary.write('Data generated using %s on: %s\n' %(pairing_version,datetime.datetime.now()))
	summary.write("***************************************************\n")
	summary.write("***************   Summary Report  *****************\n")
	summary.write("***************************************************\n\n\n")
	
	summary.write('The following input files were used to run pairing analysis:\n')
	summary.write('\n'.join(file_path_summary)+'\n\n\n')
	
	summary.write('Parsed through {0} sequences\n'.format(str(total_seqs)))
	summary.write('{0} ({1}%) sequences did not pass filters:\n'.format(str(total_seqs-passed_filter),str( round(float(  (100*(total_seqs-passed_filter))/total_seqs),1))) if total_seqs>0 else '0')
	summary.write('\t{0} sequences did not have an antibody sequence\n'.format(str(no_result)))
	summary.write('\t{0} sequences did not have a cdr3\n'.format(str(no_cdr3_error)))
	summary.write('\t{0} sequences were not productive\n'.format(str(not_productive)))
	summary.write('\t{0} sequences had an unidentifiable chain type\n\n'.format(str(invalid_chain)))

	summary.write('{0} Sequences had were identified as a VDJ recombination event\n'.format(str(num_vdj)))
	summary.write('{0} Sequences had were identified as a VJ recombination event\n\n'.format(str(num_vj)))
		
	summary.write('Sequence-pairs were collapsed into {0} sequences containing unique CDR3s\n'.format(str(num_unique_cdrh3_l3_pair)))
	summary.write('{0} ({1}%) of these unique CDRH3 were observed more than once\n\n\n'.format(num_unique_cdrh3_l3_pair_above1,str(round(float(100*num_unique_cdrh3_l3_pair_above1/num_unique_cdrh3_l3_pair),1))) if num_unique_cdrh3_l3_pair>0 else '0')
		
	summary.write("Number of unique CDRH3 sequences: %d\n" %(num_unique_cdrh3_l3_pair))
	summary.write("CDRH3 sequences were clustered using the following cluster cutoffs\n")
	summary.write("\tCluster Cutoff\tNumber of clustered sequences\n")
	
	for cpos,c in enumerate(cluster):
		summary.write('\t'+str(c*100)+'\t'+str(num_clustered[cpos])+'\n')
	
	#summary.write("\t Number of clustered sequences to %d%% CDRH3 identity: %d\n" %(cluster*100,num_clustered))
	summary.close()



def RunSingleChainClustering(annotated_file_paths,analysis_method,output_folder_path='',prefix_output_files='', annotated_file_formats=None,field_names=None,include_vl_analysis=False,cluster_cutoff = [0.96,0.96,0],files_from_igrep_database=False,productivity_function=None):
	"""
		
		Brief Description 
		
		This function serves as an addition to the current function RunPairing found in immunogrep_gglab_pairing.py file. 
	    This function preforms the same outlined steps and outputs the same information as the pairing program except it only considers the CDR3 of each chain seperately.
		Whereas RunPairing function will group together sequences by their respective R1-R2 reads (VH-VL pairs), this function will consider each chain independently. 
		Clustering therefore only considers VH only generated libraries and VL only generated libraries.
		Clustering of CDR3 sequences is again performed using Usearch program, and the output file should match the format of the RunPairing program except excludes VH-VL pairs.
		
		General algorithm:
		
		The single chain analysis program follows the following steps: 
			Step A: Go through the provided AB annotated files (IMGT, parsed IGBLAST files, etc), and filter out sequences that lack good results 
				1) Loop through all sequences provided in all files 
				2) Only select sequences that have a sequence header. 
				3) Filter out sequences that do not have a cdr3 sequence 
				4) Filter out sequences that are deemed 'unproductive' using our Productivity Rules, or a custom productivity rule provided by user 						
				5) Seperate VDJ annotated sequences from VJ annoated sequences into seperate files
				6) The user provided parameter (include_vl_analysis) which determine whether to cluster both the VDJ and VJ files. By default we only consider VDJ data.

			Step B: Group together CDRH3 sequences for each VDJ and/or VJ file
				1) We group together sequences with identical CDRH3 sequences from the filtered sequences in Step B
				2) For each CDRH3-CDRL3 group, we calculate the following as extra information for the group: 
					a) average SHM (if an SHM field is provided) and standard deviation SHM 
					b) mode of V,D, and J gene calls (if V, D, or J genes are provided). Mode of gene calls are based on top genes only (genes whose alignment scores are equal to the top alignment score)
					c) All observed unique isotypes (if isotypes are provided)
			
			Step D: Cluster CDRH3-CDRL3 pairs 
				1) Unique CDRH3-CDRL3 paired sequences are sorted by their total counts (number of sequences containing respective CDRH3-CDRL3 pair)
				2) The results from each CDRH3-CDRL3, pair whose counts are found > 1 time, are exported as a FASTA file. 
					a) The sequence header contains all of the data generated/information we are interested in for the H3-L3 pair. 
					b) The sequence value for the FASTA is the nucleotide sequence of the CDRH3 (heavy chain) sequence only 
				3) This FASTA file is used as an input file to the program USEARCH VERSION 7. 
				4) We run the Usearch program for each cluster defined in the cluster_cutoff parameter defined in this function
			
			Step E: Summarize the output from each USEARCH cluster file 
				1) For EACH cluster cutoff defined, USEARCH creates an output .uc file 
				2) We parse this file and generate a finalized output file based on the CDRH3 seed/cluster
				3) For each seed we calculate the following
					a) The average SHM for all sequences which belong in the seed 
					b) The V,D, and J call of the SEED sequence 
					c) The Isotype call of the SEED sequence
					d) A 'dominance' calculation => 
						100*Unique CDR3 Count of the seed sequence/Unique CDR3 Count of the next largest CDR3 group in the cluster
					e) A 'confidence' calculation => 
						100*Unique CDR3 Count of the seed sequence/The total counts of all sequences in cluster
					f) Diversity of all members in seed using 
						i. shannon entropy 
						ii. simpsons index 
								
		General usage:
		
		When calling this function, the most important feature for this program to run correctly is defining the type of file(s) you have submitted. While there is no uniform file type/format, we do need to know which field names in the file correspond to field names we use in the pairing program. 
		The variable, required_field_names, stores the fields we need to perform the analysis. From these fields, only the following fields are absolutely necessary (although providing all fields are ideal):
			1) seq_header
			2) cdr3_nt
					
		We already know which field names are required for files that were generated using our common annotation methods: IMGT, IGFFT (in house program), MIXCR, IGBLAST. Therefore, if you are pairing files from one of these programs AND you have not modified the output of the files created by our python wrapper functions, there is no need to define the fields. 
		However, if the file being provided is a custom file or an annotation file we do not know about, then you need to define custom field names in the input variable: field_names
										
		Example usage for pairing IMGT files: 			
			list_of_imgt_files = a list of strings providing the location of all IMGT generated files for pairing. You do not need to select a specific IMGT file or differente the 11 files by experiment, we handle that in the function
			list_of_imgt_files = [1_summary_file1.txt,5_summary_file1.txt,1_summary_file2.txt,5_summary_file2.txt] (etc  etc for the files to pair)
			RunPairing(annotated_file_paths = list_of_imgt_files,analysis_method = IMGT, output_folder_path='paired_seq_data.txt,prefix_output_files='MYPAIRINGEXP',annotated_file_formats='IMGT',field_names=None,cluster_cutoff=[0.85,0.96,0.01])
		
		Example use for pairing MIXCR files:
			list_of_mixcr_files = a list of strings providing the location of all MIXCR annotation files generated by the function ParseMIXCR
			list_of_mixcr_files = [file1.mixcr.annotation,file2.mixcr.annotation,file3.mixcr.annotation]
			RunPairing(annotated_file_paths = list_of_mixcr_files,analysis_method = MIXCR, output_folder_path='paired_seq_data.txt,prefix_output_files='MYPAIRINGEXP',annotated_file_formats='TAB',field_names=None,cluster_cutoff=[0.85,0.96,0.01])
		
		Example use for any custom file (a file type not generated by the database, parseMIXCR function, parseIgblast function, parseIGFFT function, or IMGT):
			list_of_files = a list of strings providing the location of some CSV file that we have not seen before 
			In this example, this CSV (not TAB) file only contains CDR3 nucleotide information. So every other field provided will be treated as blank for the program and not considered. In the provided file, the CDR3 field is labeled as CDR3 Nt and the sequence header is labeled as, Miseq header
			field_mappings = { #key = name of the field we need in the program, value = name of the field in the provided file 
				'cdr3_nt':	'CDR3 Nt',
				'seq_header':Miseq header'
			}
			list_of_files = [fie1.txt,file2.txt,file3.txt]
			
			RunPairing(annotated_file_paths = list_of_files,analysis_method = custom, output_folder_path='paired_seq_data.txt,prefix_output_files='MYCUSTOMPAIRINGEXP',annotated_file_formats='CSV',field_names=field_mappings,cluster_cutoff=[0.85,0.96,0.01])
			
			
		#REQUIRED FIELDS###
		#annotated_file_paths = a list of files that you want to run pairing with 
		#analysis_method:
			#A string defining which annotation program generated the results 
			#IF an annotation program is not defined in this script (variable above) yet, then use 'CUSTOM'
			#note: if analysis_method == 'CUSTOM' OR not in the predefined list, then the parameter,field_names, is required 
		
		#OPTIONAL FIELDS###
		#annotated_file_formats: 
			#accepts three possible formats: 
				#=> None => do not define file type, let program guess 
				#=> list of file types => the file type FOR EACH PROVIDED FILE 
				#=> SINGLE STRING => THIS MEANS ALL INPUTED FILES ARE OF THE SAME TYPE 
		
		#field_names: 
			#A dictionary defining which fields in the file correspond to the fields required for this program 
			#IMPORTANT: If analysis_method == CUSTOM then you must define field_names in the file. Program will raise exception otherwise
			#if provided, the structure of field names is as follows: 
				#keys of variable: 'vgene','jgene','dgene','raw_seq_nt','cdr3_nt','seq_header','shm','full_len_ab_aa'
				#values => for each key,provide the field name in the file that corresponds to the key 
				
		#output_folder_path => path of output folder. If empty, then a new folder is created. If folder path does not exist, then creates a new folder using provided folder path
		#prefix_output_files => prefix string to use for naming output files
		#cluster_cutoff => clustering cutoff for determining clustered pairs  (accepts a single number or a list of 3 numbers)
			If its a single number, i.e. cluster = 0.96, then only that clustering id cutoff is used 
			If a list is provided, then a range of clusters are attempted defined by the list [min cluster id, max cluster id, step increase of cluster id's]
		#files_from_igrep_database => indicates whether the provided files were downloaded from the database (this is important because files from the database will have the same field names regardless of the annotation type. Therefore, we need to use these files from database
		
	"""						
		
	required_field_names = get_required_field_names()	
			
	###########INITIALIZATION#################
	##########################################
	#JUST SET UP INPUT PARAMETERS. GO THROUGH PROVIDED PARAMETERS AND DETERMINE HOW TO RUN PROGRAM#		
	#cluster cutof must be between 0 and 1 
	if isinstance(cluster_cutoff,list):
		if len(cluster_cutoff)!=3:
			raise Exception('The cluster range must be defined as a 3 element list of :[min cluster,max cluster, step increase in cluster]')
		cluster_cutoff = [c for c in FloatRange(cluster_cutoff[0],cluster_cutoff[1],cluster_cutoff[2])]
	else:
		if cluster_cutoff<=0 or cluster_cutoff>1:
			raise Exception('Cluster cutoff must be a floating number between 0 and 1')			
		cluster_cutoff = [cluster_cutoff]		
	

	cluster_cutoff = [c for c in sorted(list(set(cluster_cutoff)))]			
	
	#if user just passe in a single string for single file, then make it a list 
	if not isinstance(annotated_file_paths,list):
		annotated_file_paths = [annotated_file_paths]			
	#user passed in the file type for EACH provided file 
	if isinstance(annotated_file_formats,list):
		if len(annotated_file_formats)!=len(annotated_file_paths):
			raise Exception('The user defined a file format for each file. But the number of file formats does not match the number of input files. If the user wants to use a single format for all files, then simply pass in a string not list defining file format')
		
		for i in range(len(annotated_file_formats)):
			annotated_file_formats[i] = annotated_file_formats[i].upper()
			if annotated_file_formats[i] not in allowed_file_formats:
				raise Exception('The provided file format, {0}, is not currently supported in this program'.format(annotated_file_formats))		
	#user just passed in a string. copy this to all a list contaiing same string for all files 
	elif isinstance(annotated_file_formats,basestring):		
		annotated_file_formats = annotated_file_formats.upper()		
		if annotated_file_formats not in allowed_file_formats:
			raise Exception('The provided file format, {0}, is not currently supported in this program'.format(annotated_file_formats))	
		annotated_file_formats = [annotated_file_formats]*len(annotated_file_paths)	
	#user provided ouput folder 
	if output_folder_path:		
		#folder path provided does not exist 
		if not os.path.isdir(output_folder_path):			
			#make folder, but ONLY USE BASENAME. FORCE the folder to exist within IGREP project 
			if output_folder_path[-1]=='/':
				output_folder_path=output_folder_path[:-1]
			parent_dir = filesystem.ExperimentDirs().make_dir(os.path.basename(output_folder_path))+'/'
			print "WARNING, OUTPUT DIRECTORY NOT FOUND. CREATED A NEW DIRECTORY FOR OUTPUT CALLED: {0}".format(parent_dir)
		else:
			#this is the output folder 
			parent_dir = output_folder_path+'/' if output_folder_path[-1]!='/' else output_folder_path									
	else:
		#user did not pass in output folder, so make a new folder
		folder_name = 'PAIRING_RESULT_'+str(datetime.datetime.now()).replace(' ','').replace('/','').replace('-','').replace(':','').split('.')[0]
		#make a new folder
		parent_dir = filesystem.ExperimentDirs().make_dir(folder_name)+'/'		
	
	#replace empty spaces in folder name (can be annoying..)
	parent_dir = parent_dir.replace(' ','_')
	prefix = parent_dir+prefix_output_files.replace(' ','_') if prefix_output_files else parent_dir+'VH_VL_PAIRING'	
	#make sure provided analysis settings are correct/make sure we have all the data we need for pairing 
	analysis_method=analysis_method.upper()
	if analysis_method == 'CUSTOM' or analysis_method not in supported_analyses:
		if not field_names:
			raise Excpetion('You have defined the following analysis method: {0}. However, we cannot use this file until we know the field names for the following variables: {1}. Please rerun function defining these values as key:value format (key = required field name, value = field name in file)'.format(analysis_method,','.join(required_field_names.keys())))
		for each_req_field,value in field_names.iteritems():
			if each_req_field in required_field_names:
				required_field_names[each_req_field] = value							
	else:
		#the user did not provide field_names, but thats ok because we already know the analysis type and have it hardcoded in variable above 
		if not field_names:
			for each_req_field,value in pairing_settings[analysis_method]['fields_for_analysis'].iteritems():
				if each_req_field in required_field_names:
					required_field_names[each_req_field] = value							
		#the user provided the field name values 
		else:
			for each_req_field,value in field_names.iteritems():
				if each_req_field in required_field_names:
					required_field_names[each_req_field] = value				
		if productivity_function == None:
			productivity_function = pairing_settings[analysis_method]['productivity_function']
		
	#at the very least, both the field for seq header and the field for cdr3 nt must be defined 
	if not required_field_names['seq_header'] or not required_field_names['cdr3_nt']:
		raise Exception('At the very least you must define what field name corresponds to the seq_header and what field name corresponds to the nucleotide cdr3. The following fields were defined by the variable field_names: '+str(field_names))
	
	#NOW tell the user whether they are will be missing any fields with the provided file for pairing (therefore they know whether to expect weird results)
	for each_req_field,value in required_field_names.iteritems():
		if not value:
			print "WARNING, THE FOLLOWING FIELD IS NOT DEFINED: {0}. PAIRING WILL NOT INCLUDE INFORMATION FOR THIS FIELD".format(each_req_field)			
	
	#WE PREDICT TWO WAYS TO PROVIDE IMGT FILES: 
		#METHOD 1=> PROVIDE THE ORIGINAL 11 IMGT FILES, IN THIS CASE, THEN ANLAYSIS_METHOD == IMGT AND ALSO ANNOTATED_FILE_FORMATS = IMGT OR [IMGT,IMGT,IMGT]
		#METHOD 2=> PROVIDE A SUMMARIZED FILE COMBINING RESULTS FROM IMGT. FOR EXAMPLE DOWNLOADING IMGT RESULTS FROM DATABASE. IN THIS CASE THEN ANALYSIS_METHOD == IMGT BUT!! THE ANNOTATED_FILE_FORMATS DOES  NOT EQUAL IMGT 		
	#SO IF METHOD 1, THEN WE WILL REORGANIZE THE PROVIDED FILES. WE DO THIS BECAUSE WE WANT TO READ MULTIPLE EXPERIMENTS simultaneously (READ MULTIPLE SETS OF IMGT FILES AT THE SAME TIME TO MAINTAIN THE SIZE OF CODEDICT AS SMALL)
	if analysis_method=='IMGT':
		if annotated_file_formats == None:			
			annotated_file_formats=['IMGT']*len(annotated_file_paths)		
			imgt_groups = readfile.GroupIMGTFiles(annotated_file_paths)		
			annotated_file_paths = [vals for group,vals in imgt_groups.iteritems()]			
		elif isinstance(annotated_file_formats,list):
			if len(annotated_file_formats)!=len(annotated_file_paths):
				raise Exception('The provided number of file formats does not match the provided files')
			if annotated_file_formats[0] == 'IMGT':
				annotated_file_formats=['IMGT']*len(annotated_file_paths)		
				imgt_groups = readfile.GroupIMGTFiles(annotated_file_paths)		
				annotated_file_paths = [vals for group,vals in imgt_groups.iteritems()]			
		elif annotated_file_formats == 'IMGT':
			annotated_file_formats=['IMGT']*len(annotated_file_paths)		
			imgt_groups = readfile.GroupIMGTFiles(annotated_file_paths)		
			annotated_file_paths = [vals for group,vals in imgt_groups.iteritems()]			
		else:
			annotated_file_formats = [annotated_file_formats]*len(annotated_file_paths)		
		
		#print annotated_file_paths
								
	#SETUP COMPLETE. PROCEED TO RUNNING FUNCTION#			
	############################################################################
		
	#the following files will be generated during the program
	global dict_summary
	global codeDict_output_dict
	#global collapsed_output_dict
	global error_file
	global usearch_cluster_file
	global collapsed_output_file
	global summary_file
	global annotation_path	
	global clustered_final_output_file		
		
	
	#stores a file of unique CDR3 pairs. Its a tab delimited file where the first column is the counts for that CDRH3 pair. Remaining columns contain AB information such as CDR3, reocombination type, V,D,J genes, isotypes, average SHM
	collapsed_output_file = prefix + "_CDR3_identical_nt_pairs"
	#stores a file of the parsed output from USEARCH clustering. It will store the seed CDR3 sequence and then store all other sequences above 1 read that were clustered into that seed. 
	#it is a TAB delimited file 
	clustered_output_file = prefix + "_CDR3_raw_clustered_nt_pairs_over1read"
	#stores a summary file of all CDR3 clusters created from usearch. It stores the seed sequence, then includes information about this size of the cluster and develops "confidence" and "dominance" terms for each cluster size
	clustered_final_output_file = prefix + "_CDR3_clustered_nt_pairs_over1read"
	#this is a JSON file that stores all of the sequences which passed filtering. It stores all the sequence information used during pairing. 
	#the key in this file represents the fullcode from a miseq header and the value represents VDJ OR VJ sequence data. 
	codeDict_output_dict= prefix + "_codeDict"		
	#this is a JSON file created after USEARCH clustering. For each selected cluster it will store the seed as the key and all of the AB information for the CDR3 sequence found in the sequence header submitted to usearch
	#consider removing this file....
	clustered_output_dict= prefix + "_clustered_dict"
	#this is a series of JSON strings 
	#each line in the file represents either a VH or VL result.
	#when using json.loads(line) on each line, it yeilds a 2 element list. the first element is a dictionary of VDJ data only. The second element is a dictionary of VJ data only.
	dict_summary= prefix + "_dict_summary_line_by_line"	
	#this will generate a final summary file describing how many sequences were parsed, how many were sucessful, and how many clusters were formed. it also includes gene usage for H/L
	summary_file=prefix+"_Summary.txt"		
	#the input FASTA file submitted to usearch 
	usearch_cluster_file = parent_dir+"clustered_input_file"
	#the output file generated by usearch program 
	usearch_output_file = parent_dir+'clustered_output_file'	
	
	#THIS IS THE ACTUAL CODE FOR RUNNING PAIRING (THIS CALLS ALL THE FUNCTIONS FOR PAIRING)	
	#first create the variables 'codeDict' and 'collapsedDict' which will store information regarding R1/R2 pairs from the annotation files 	
	[total_seqs,passed_filter,no_result,no_cdr3_error,not_productive,invalid_chain,num_unique_cdr3,num_unique_cdr3_pair_above1,num_vdj,num_vj] = add_to_dict_memory_safe(annotated_file_paths,annotated_file_formats,required_field_names,productivity_function,parent_dir,analysis_method,include_vl_analysis)		
				
	number_of_clusters = []
	for c in cluster_cutoff:
		usearch_file = usearch_output_file+'.'+str(c)+'.uc'
	
		#CODEICT CREATED, RUN USEARCH  for all clustered settings
		print 'Clustering CDRH3 reads with %d%% identity.\n' %(c*100)
		#run USEARCH
		os.system("usearch7 -cluster_smallmem {2} -minhsp 10 -minseqlength 10 -usersort -id {0} -centroids c.fa -uc {1}.{0}.uc".format(str(c),usearch_output_file,usearch_cluster_file))
		print 'USEARCH complete'	
		print 'Compiling clustered reads file.\n'
		#After populate_clusters, clustered_dict_raw[cluster number]=[total count, infoSeed, infoHit 1, infoHit 2, etc]		
		
		#Dict clusters based on %identity of heavy chain. light chain ignored. cluster variable defined at top of script.	
		clustered_dict_raw = populate_clusters(usearch_file,delimiter='\t')			
			
		clustered_output=open(clustered_output_file+'.'+str(c)+'.txt','wb')		
		clustered_save=open(clustered_output_dict+'.'+str(c)+'.txt','w')	
		for key in sorted(clustered_dict_raw, key=clustered_dict_raw.get, reverse=True):			
			CDRH3=clustered_dict_raw[key][1].split(":")[3]
			clustered_output.write("%d\t%s\n" %(clustered_dict_raw[key][0],CDRH3))
			
			for i in range(1,len(clustered_dict_raw[key])):
				clustered_output.write("%s\n" %clustered_dict_raw[key][i].replace(":","\t"))
			clustered_output.write("\n")	
			
		print "Saving clustered data.\n"
		len_clustered=len(clustered_dict_raw)
		json.dump(clustered_dict_raw, clustered_save)	
		
		clustered_save.close()
		clustered_output.close()	
		
		#Prepare clustered file with over 1 read and only the top light chain paired to each heavy chain
		#1st number is rank. 
		#2nd number is dominance factor ( LC reads / total number of reads * 100) Good=100, Bad=0
		#3rd number is confidence factor (100 - 100 * (2nd rank LC)/(1st rank LC)) Good=100, Bad=0
		#clustered_dict_raw[cluster number]=[total count, infoSeed, infoHit 1, infoHit 2, etc]	
		#info=(counts,CDRH3_trans,CDRL3_trans,CDRH3,CDRL3,VHgene,DHgene,JHgene,VLgene,JLgene,IgH,IgL,Hmut,Lmut,HLen,LLen)
		print('Compiling finalized clustered output file.\n')
		Compile_Clusters(clustered_dict_raw,c)		
	
	
		number_of_clusters.append(len(clustered_dict_raw))
	
	#Add in gene histograms to summary file
	gene_hist(clustered_dict_raw,summary_file,1)
	gene_hist(clustered_dict_raw,summary_file,2)
	gene_hist(clustered_dict_raw,summary_file,3)
	gene_hist(clustered_dict_raw,summary_file,4)
	gene_hist(clustered_dict_raw,summary_file,5)
	#generate a summary file 
	if (analysis_method=='IMGT' and annotated_file_formats[0]=='IMGT'):
		annotated_file_path_summary = []
		for each_group in annotated_file_paths:
			path = each_group[0].split('/')			
			filename = path[-1]
			path = os.path.dirname(each_group)#+'/'# '/'.join(path[:-1])+'/'
			filename = filename.split('_')
			if len(filename)>=3:
				annotated_file_path_summary.append(os.path.join(path,'_'.join(filename[2:])))# path+'_'.join(filename[2:]))
			else:
				annotated_file_path_summary.append(os.path.join(path,'_'.join(filename)))# path+'_'.join(filename))
	else:
		annotated_file_path_summary = annotated_file_paths
	WriteSummaryFile(annotated_file_path_summary,total_seqs,passed_filter,no_result,no_cdr3_error,not_productive,invalid_chain,num_unique_cdr3,num_unique_cdr3_pair_above1,cluster_cutoff,number_of_clusters,num_vdj,num_vj)
	
			
	analysis_files_to_report = [collapsed_output_file,summary_file]
	for c  in cluster_cutoff:
		analysis_files_to_report.extend([usearch_output_file+'.'+str(c)+'.uc',clustered_output_file+'.'+str(c)+'.txt',clustered_final_output_file+'.'+str(c)+'.txt'])
	
	
	print('All paired clustering complete')
	return {'annotation_files':[],'analysis_files':analysis_files_to_report}



