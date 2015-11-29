#if you are using the gapless alignment method in this funciton, then it requires the executive node to have installed the pyfftw module and fftw3.  
#instructions for installation on machine are provided in function -> nt_fft_align_tools

import json
from Bio import SeqIO
from Bio.Seq import Seq
import os 
import sys
import re
from datetime import datetime
import time

import copy 
from collections import OrderedDict
from collections import defaultdict

import numpy as np


from immunogrep_global_variables import idIdentifier
from immunogrep_global_variables import translation_var #this key will signify the translation/translator key 
from immunogrep_global_variables import descriptor_symbol #this key will signify the translation/translator key 
from immunogrep_global_variables import fasta_file_delimiter #this key will signify the translation/translator key 
import immunogrep_fft_align_tools as fft_tools
import immunogrep_useful_functions as useful
import immunogrep_read_file as readwrite


FileDelimFields = ['Header','Sequence','SEQ_ID','Sequence strand','Isotype','Mismatches', 'MaxScore','Isotype barcode orientation','Alignment start position','Percent similarity','Alignment length', 'Recombination type', 'Notes','Command']
translator = {"ANALYSIS_NAME": "GEORGIOU_INHOUSE", # NAME OF THE ANALYSIS 
			  "RECOMBINATION_FIELD":{ #THIS TELLS THE PROGRAM HOW TO DETERMINE WHETHER AN ANALYSIS/QUERY RESULT (from file) IS VDJ OR VJ
					"FIELD_NAME": "Recombination type", #name of the field in the file that will give information regarding the recombination type (VDJ OR VJ)					
					"EXPLICIT":True,
			 },
			"FIELDS_UPDATE":{ 
				#key = field name  in database 
				#value = field name in file
				#this will map all of the fields in the file to the proper location in the database. I.E. If I list VGENES as the column name/field name, then i want to map VREGION.VGENES:VGENES (because VREGION.VGENES is the name in the database)							
				idIdentifier:idIdentifier,
				"SEQUENCE":"Sequence",				
				"COMMAND":"Command",
				"SEQUENCE_HEADER":"Header",												
				"ISOTYPE.GENE":"Isotype",
				"ISOTYPE.MISMATCHES":"Mismatches",
				"ISOTYPE.PER_ID":"Percent similarity"
			}
	}
	

#open fasta file and put sequences into a list. first column is sequence, second column is header. 
def moveSeqToList(seqHandle):
	#find the maximum barcode length and read barcodes into list
	#barHandle = open(barcode_fasta_location,"rU")		
	seqList=[]
	for sequence in SeqIO.parse(seqHandle, "fasta"): 
		seqList.append([str(sequence.seq).upper(),sequence.id])	
	return seqList

def defaultDictionary():
	isotypeDict = {
		'Sequence':None,
		'Header':None,		
		'MaxScore':None,
		'Isotype':None,
		'AlgnPos':None,
		'PercentSimilarity':None,
		'DebugData':None
	};
	return isotypeDict




def defaultBarcodes():
	barcodeSeqList = [['CCTCCACCAAGGGCCCATCGCAG','IGHG'],
					  ['GGAGTGCATCCGCCCCAACC','IGHM'],
					  ['CATCCCCGACCAGCCCCAAGC','IGHA'],
					  ['GAACTGTGGCTGCACCATCT','IGK'],
					  ['GTCACTCTGTTCCCGCCCTC','IGL']]
	
	return barcodeSeqList
	

#this will attemp to use the barcode isotype name to guess the recombination type 
def GuessRecombType(isotype):
	mc = isotype[:2].upper() #molecular component
	if mc!='IG' and mc!='TR':
		return 'UNK' #cannot predict this isotype 
	#excpecting light chain/vj recombination types => IGK, IGL, TRA chains
	if isotype[:3].upper() in ['IGK','IGL','TRA']:
		return 'VJ'
	else:
		return 'VDJ'
		

#find the maximum barcode length and read barcodes into list
def readBarcodeFile(barcode_fasta_location):
	barHandle = open(barcode_fasta_location,"rU")				
	barcodeSeqList = moveSeqToList(barHandle) #-> in function fasta_file_functions
	
	return barcodeSeqList
	
def maxSeqLen(IFFile,seq_string_name):
	numSeq = 0
	maxLen = 0
	while not(IFFile.IFclass.eof):
		line_row = IFFile.IFclass.read()
		if line_row:
			if seq_string_name in line_row:
				lenT = len(line_row[seq_string_name])
				if lenT>maxLen:
					maxLen =lenT
				numSeq+=1
	IFFile.IFclass.close()
	return [maxLen,numSeq]


strand_orientation_list = ['+','-']

#seq_fasta_location,barcode_fasta_location="",output_file_location="scratch/isotype_output.json"
#penalize_truncations
#if this is true then we will not score against alignments that do not cover 100% over the barcode. so if we only align to the first 8 base of a barcode 10 bases long, then we only look at that barcode length
#IF false, then aligning only 8 bases to a barcode 10 bases long will result in 2 mismatch calls 		
def isotype_sequences(input_file,input_file_type,barcode_file='',output_file=None,output_format='TAB',seq_var='sequence',header_var='header',helper_fields = {},alignment_settings = {},analysis_name = None):		
	#####OVER HEAD FUNCTIONS
	
	help_1 = defaultdict(str,copy.deepcopy(helper_fields))
	recombination_var = help_1['recombination_var']
	strand_field = help_1['strand_field']
	end_of_ab_field = help_1['end_of_ab_field']
		
	
	al_1 = copy.deepcopy(alignment_settings)
	
	penalize_truncations = al_1['penalize_truncations'] if 'penalize_truncations' in al_1 else True
	
	minimum_alignment_length = al_1['minimum_alignment_length'] if 'minimum_alignment_length' in al_1 else 15
	
	#0=> only consider barcodes as provided
	#1=> only consider the reverse complmeent of barcodes provided 
	#2=> consider both strands 
	search_rc = al_1['search_rc'] if 'search_rc' in al_1 else 2
	
	allowed_mismatches_in_alignment = al_1['allowed_mismatches_in_alignment'] if 'allowed_mismatches_in_alignment' in al_1 else 2
	
	#the sequence filed provided is the sequence of the SENSE AB gene not the antisense
	#when False, will consider both the forward and reverse copmlmement of sequence 
	strand_corrected = al_1['strand_corrected'] if 'strand_corrected' in al_1 else False
		
				
	#file locations
	seq_fasta_location =input_file#  functionVars["folder_location"]+functionVars["input_file"] #location of input file
	
	translator_field = copy.deepcopy(translator)
	
	if analysis_name:
		translator_field['ANALYSIS_NAME'] = analysis_name.upper()
	
	
	translator_field = {translation_var:translator_field}
	if output_file == None or output_file==input_file:
		output_file = useful.removeFileExtension(input_file)+'.isotype.annotation'
	
	output_file_location = output_file
		
		
	output_file_format = output_format #functionVars['write_format']
	#seqHandle = open(seq_fasta_location,"rU")
		
	outHandle = open(output_file_location,'w')		
	outHandle.write(descriptor_symbol+json.dumps(translator_field)+'\n')#write a translator line to this file so that we know how to add results to database 
	if output_format == 'TAB' or output_format == 'CSV':
		outHandle.write('\t'.join(FileDelimFields)+'\n')
	
	if not barcode_file:# 'barcodefilename' in functionVars:
		#manually using these primers
		barcodeSeqList = defaultBarcodes()
	elif not(os.path.isfile(barcode_file)):
		print('Barcode file not found! Using default barcodes')		
		#manually using these primers
		barcodeSeqList = defaultBarcodes()
	else:
		barcodeSeqList = readBarcodeFile(barcode_file)
		
	command_string = json.dumps({'Barcodes':barcodeSeqList,'mismatch_cutoff':allowed_mismatches_in_alignment,'penalize_truncations':penalize_truncations,'minimum_length_cutoff':minimum_alignment_length})
	
	
	
	iffile = readwrite.immunogrepFile(filelocation=seq_fasta_location,filetype=input_file_type)
		
	#get maximum length of sequences in file 
	[maxLen,numSeq] = maxSeqLen(iffile,seq_var) 	
	
	#make a call to the generator for alinging sequences to isotypes 
	guessed_num_bases_after_jgene = 60
	isotype_predictor =fft_tools.BarcodeAligner(barcodeSeqList,penalize_truncations,search_rc,allowed_mismatches_in_alignment,minimum_alignment_length,nmax=maxLen,nmin=guessed_num_bases_after_jgene)		
					
	###END OF OVERHEAD FUNCTIONS
	
	
	#now lets read through sequences and start alignining
	algnLim = 10
	currentSeq = 0
	overlap_len = 10
	
	#seqHandle=open(seq_fasta_location,"rU")
	counter = 0
	startPer = 0
	
	num_isotype_found = {}
	total_isotype_found = 0
	total_found_score=0
	total_notfound_score=0
	
	print("Starting isotyping analysis for {0} sequences".format(numSeq))

	
	totaltime = 0
	a = int(round(time.time()))
	found = 0 
	
	iffile = readwrite.immunogrepFile(filelocation=seq_fasta_location,filetype=input_file_type);
	summary_data = {'found':0,'top_isotype':defaultdict(int),'average_mismatch':0,'average_num_isotype_found':0}
	
	for line_row in iffile.read():			
		jsonVar = {}
		if not line_row:
			continue
		
		if header_var in line_row:
			if idIdentifier in line_row:
				jsonVar[idIdentifier] = line_row[idIdentifier]
				jsonVar['Header'] = line_row[header_var]
			else:
				[header,id] = GrabAdditionalHeaderInfo(line_row[header_var])
				jsonVar[idIdentifier] = id			
				jsonVar['Header'] = header
			
		
		if seq_var not in line_row or line_row[seq_var]=='':		
			jsonVar['Sequence']=''					
			jsonVar['Notes'] = 'No sequence found'			
			writeSeqResult(outHandle,jsonVar,output_format)			
			continue
								
		#allow the user to monitor what percent of the sequences have been processed					
		startPer = useful.LoopStatus(counter,numSeq,10,startPer)
		
		bestScore = 0;
		bestBarcode = -1;
			
		jsonVar['Sequence'] = line_row[seq_var]
		jsonVar['Command'] = command_string
		counter+=1		
				
		seqFwd = jsonVar['Sequence']
		
		if strand_corrected:
			all_seqs = [seqFwd]
		else:
			all_seqs = [seqFwd,str(Seq(seqFwd).reverse_complement())]
		
		
		found_strand =''
		for pos,each_seq in enumerate(all_seqs):										
			#determine if we should take a substring of the sequence 
			#basically, only consider nucleotides AFTER the end of the ab field 
			if end_of_ab_field in line_row and line_row[end_of_ab_field]!='':
				try:
					end_of_ab = int(line_row[end_of_ab_field])							
				except:
					end_of_ab = 0
				#take substring
				if end_of_ab-overlap_len<len(each_seq) and end_of_ab-overlap_len>=0:
					each_seq = each_seq[end_of_ab:]																							
										
			isotypes_results = isotype_predictor.AlignToSeq(each_seq)
			if isotypes_results:
				found_strand = strand_orientation_list[pos]
				break
		
		
		if isotypes_results:
			found += 1 
			jsonVar = dict(jsonVar.items()+isotypes_results.items())
			
			jsonVar['Sequence strand'] = found_strand			
			
			
			if recombination_var in line_row and line_row[recombination_var]:
				#always trust the recombination type from input file IF provided
				jsonVar['Recombination type'] = line_row[recombination_var]
			else:
				#if there is no results then attemp to guess it our selves
				jsonVar['Recombination type'] = GuessRecombType(jsonVar['Isotype'][0])
			
			summary_data['top_isotype'][jsonVar['Isotype'][0]]+=1
			summary_data['average_num_isotype_found']+=len(jsonVar['Isotype'])
			summary_data['average_mismatch']+=jsonVar['Mismatches'][0]
		else:
			if recombination_var in line_row and line_row[recombination_var]:
				#always trust the recombination type from input file IF provided
				jsonVar['Recombination type'] = line_row[recombination_var]
		
		
			jsonVar['Isotype'] = ''
			jsonVar['Notes'] = 'Could not identify isotype with alignment score above threshold'
			summary_data['top_isotype']['NotFound']+=1
				
		writeSeqResult(outHandle,jsonVar,output_format)
				
		
	
	b = int(round(time.time()))
	
	summary_data['found'] = found
	if found:
		summary_data['average_mismatch'] = summary_data['average_mismatch']/float(found) 
		summary_data['average_num_isotype_found'] = summary_data['average_num_isotype_found']/float(found)
		
	totaltime=(b-a)			
	
	print "time: "
	print totaltime
	
	print "Summary of identified isotypes:"
	print summary_data
	
	#if total_isotype_found>0:
	#	print "\nAverage score for identified isotypes:"	
	#	print str(total_found_score/float(total_isotype_found))		
	
	#if numSeq-total_isotype_found>0:	
	#	print "\nAverage score for unidentified isotypes:"	
	#	print str(total_notfound_score/float(numSeq-total_isotype_found))
			
	outHandle.close()	
	#if output_file_format=="txt":
	#	JSON_to_TXT(output_file_location, output_file_location, True,{'Header':1,'Seq':2,'dir':3,'isotype':4,'algnPos':5,'maxscore':6,'bestscore':7})
	return output_file 
	
#looped through list, check to see if max sscore is above threshold hodl		
#if bestScore>0:#cutoff_threshold*maxBarcode:

def GrabAdditionalHeaderInfo(header):
	tmp = header.split(fasta_file_delimiter)
	
	if len(tmp)>1:
		additional_info = json.loads(tmp[-1])
	
	else:				
		additional_info = {}
	id = additional_info[idIdentifier] if idIdentifier in additional_info else ''
	
	
	return [tmp[0],id]
	


def writeSeqResult(outHandle,jsonVar,output_format):
	jsonVar = defaultdict(str,jsonVar)
	if output_format == 'CSV':
		jsonVar['Mismatches'] = '|'.join([str(s) for s in jsonVar['Mismatches']])
		jsonVar['MaxScore'] = '|'.join([str(s) for s in jsonVar['MaxScore']])
		jsonVar['Isotype'] = '|'.join(jsonVar['Isotype'])
		jsonVar['Alignment start position'] = '|'.join([str(s) for s in jsonVar['Alignment start position']])
		jsonVar['Percent similarity'] = '|'.join([str(s) for s in jsonVar['Percent similarity']])
		jsonVar['Alignment length'] = '|'.join([str(s) for s in jsonVar['Alignment length']])
		jsonVar['Isotype barcode orientation'] = '|'.join([str(strand_orientation_list[s]) for s in jsonVar['Direction']])
	else:
		jsonVar['Mismatches'] = ','.join([str(s) for s in jsonVar['Mismatches']])
		jsonVar['MaxScore'] = ','.join([str(s) for s in jsonVar['MaxScore']])
		jsonVar['Isotype'] = ','.join(jsonVar['Isotype'])
		jsonVar['Alignment start position'] = ','.join([str(s) for s in jsonVar['Alignment start position']])
		jsonVar['Percent similarity'] = ','.join([str(s) for s in jsonVar['Percent similarity']])
		jsonVar['Isotype barcode orientation'] = ','.join([str(strand_orientation_list[s]) for s in jsonVar['Direction']])
		jsonVar['Alignment length'] = ','.join([str(s) for s in jsonVar['Alignment length']])
	
	if output_format == 'TAB':
		outHandle.write('\t'.join([jsonVar[f] for f in FileDelimFields])+'\n')				
	elif output_format == 'CSV':
		outHandle.write(','.join([jsonVar[f] for f in FileDelimFields])+'\n')				
	elif output_format == 'JSON':
		outHandle.write(json.dumps(jsonVar)+'\n')		


                  				
#isotype_sequences('scratch/igrep/Exp_00051_test/Demo_3.igfft.annotation','TAB',seq_var='Strand_Corrected_Sequence',header_var='Header')
#isotype_sequences('scratch/igrep/Exp_00048_071715_095242/DonorV5.fastq','FASTQ')
#seq_var='sequence',header_var='header',analysis_name=None,functionVars={}):
