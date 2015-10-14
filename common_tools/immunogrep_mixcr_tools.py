from immunogrep_read_file import immunogrepFile
from Bio.Seq import Seq 
from Bio.Alphabet import generic_dna
import re
import os
import json
from immunogrep_global_variables import translation_var #this key will signify the translation/translator key 
from immunogrep_global_variables import descriptor_symbol #this key will signify the translation/translator key 
from immunogrep_global_variables import idIdentifier #variable defining what SEQ_ID variable is in database 
from immunogrep_global_variables import fasta_file_delimiter
#import immunogrep_useful_immunogrep_functions as useful
import immunogrep_useful_functions as useful
import subprocess

from collections import defaultdict
from collections import OrderedDict
import time 
#VERSION 5 VS 27...27 MAY BE BAD
#A FASTER PARSING FUNCTION AND AN UPDATE OF MIXCR TOOLS WAS APPLIED ON 8-20-2015. Previosu versions were at V22
#Version 21 was previous version

chainDic = {
		'IGH':["heavy","VDJ","IGH"],
		'IGK':["light","VJ","IGK"],
		'IGL':["light","VJ","IGL"],
		'TRA':["alpha","VJ","TRA"],
		'TRB':["beta","VDJ","TRB"]
	}
	
#We will need to update the database with the results from IgFFT.  In order
#to update teh database, we need a translator, so that we know what fields go where in the database
def DatabaseTranslator(input_dictionary = {}):
	key = translation_var
	
	translator = {					
			"ANALYSIS_NAME": "MIXCR", # NAME OF THE ANALYSIS 
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
				"COMMAND":"Command",
				"SEQUENCE_HEADER":"Seqheader",				
				"QUALITY_SCORE":"Read(s) sequence qualities",
				"NOTES":"Notes",
				"PREDICTED_AB_SEQ.NT":"Full NT",
				"PREDICTED_AB_SEQ.AA":"Full AA",
				"STRAND":"Orientation",
				"PREDICTED_CHAIN_TYPE":"Chain",
				"PRODUCTIVE":"Productivity",
				"LOCUS_NAME":"Locus",
				"VREGION.SHM.NT":"VGENE: Shm.nt",
				"VREGION.SHM.NT_PER":"VGENE: Shm.per",
				"JREGION.SHM.NT":"JGENE: Shm.nt",
				"JREGION.SHM.NT_PER":"JGENE: Shm.per",
				'VREGION.VGENE_QUERY_START':'VGENE: Query start',
				'VREGION.VGENE_QUERY_END':'VGENE: Query end',
				"VREGION.FR1.NT":"N. Seq. FR1",
				"VREGION.FR1.AA":"AA. seq. FR1",
				"VREGION.CDR1.NT":"N. Seq. CR1",
				"VREGION.CDR1.AA":"AA. seq. CDR1",			
				"VREGION.FR2.NT":"N. Seq. FR2",
				"VREGION.FR2.AA":"AA. seq. FR2",
				"VREGION.CDR2.NT":"N. Seq. CDR2",
				"VREGION.CDR2.AA":"AA. seq. CDR2",
				"VREGION.FR3.NT":"N. Seq. FR3",
				"VREGION.FR3.AA":"AA. seq. FR3",
				"VREGION.VGENES":"All V hits",
				"VREGION.VGENE_SCORES":"All V scores",
				"CDR3.NT":"N. Seq. CDR3",
				"CDR3.AA":"AA. seq. CDR3",
				"DREGION.DGENES":"All D hits",
				"DREGION.DGENE_SCORE":"All D scores",
				"JREGION.FR4.NT":"N. Seq. FR4",		
				"JREGION.FR4.AA":"AA. seq. FR4",
				"JREGION.JGENES":"All J hits",
				"JREGION.JGENE_SCORES":"All J scores",
				"ISOTYPE.GENE":"All C hits",
				"ISOTYPE.SCORES":"All C scores",
				'JREGION.JGENE_QUERY_START':'JGENE: Query start',
				'JREGION.JGENE_QUERY_END':'JGENE: Query end',
				'COMMAND':'Command'
			}
	}
	
	input_dictionary[key] = translator		
	
	return input_dictionary


def RunMixcr(inputname,outputlocation,filetype='FASTQ',loci = [],species='',exportPrettyAlignment=False):
	filetype=filetype.upper()
	renamed_file = False
	if filetype=='FASTA' and not(filetype.endswith('.fasta')):
		oldname = inputname
		inputname = inputname+'.fasta'
		renamed_file = True
		#IN ORDER FOR MIXCR TO RUN FASTA FILES, YOU MUST MAKE SURE IT ENDS IN A FASTA FORMAT 
		os.rename(oldname,inputname)
		
	inputname_run = inputname.replace(' ','\ ')
	return_output = outputlocation #python will not read files with \  as literal.
	outputlocation=outputlocation.replace(' ','\ ')
	cmdsS='java -jar -Xms4g -Xmx8g mixcr.jar '
	alignment_settings = ''
	if species:
		alignment_settings+='-OjParameters.parameters.mapperMaxSeedsDistance=5 -s '+species
	if loci:
		alignment_settings+=' -l '+','.join([l.strip().upper() for l in loci])
	
	#-OvParameters.geneFeatureToAlign=VRegion 
	cmds=cmdsS+'align '+alignment_settings+' --save-description -f '+inputname_run+' '+outputlocation+'.vdjca'+' -r '+outputlocation+'-alignment.log'+' -t 2 '
	
		
	subprocess.call(cmds,shell=True)
	subprocess.call(cmdsS+'exportAlignments -readId -descrR1 --preset full '+outputlocation+'.vdjca '+outputlocation,shell=True)
	if exportPrettyAlignment:
		subprocess.call(cmdsS+'exportAlignmentsPretty '+outputlocation+'.vdjca '+outputlocation+'.alignmentpretty',shell=True)	
	if renamed_file:
		os.rename(inputname,oldname)
	#os.system(cmds)
	#os.system(cmdsS+'exportAlignments '+outputlocation+'.vdjca '+outputlocation)
	command_val = {'mixcr_version':'1.3','loci':','.join(loci),'species':species}
	return [return_output,command_val]


# Note: If there is N in the original Input raw sequence file, MixCR will attempt to correct the sequence and report the result thus making 
# matching the sequences and to recover the header more challenging
# So the solution for now is to replace every N in the original raw sequence with [ATGC] and use Regular Expression to match the variations



def extractScores(genestr):
	tmp=genestr.split(',')
	geneusage=[]
	scores=[]
	for i in tmp:
		tmp2=i.split('(')
		geneusage.append(tmp2[0])
		scores.append(tmp2[1].replace(')',''))
	locus=geneusage[0][0:3]
	if locus in chainDic:
		recomb = chainDic[locus][1]
		chain = chainDic[locus][0]
	else:
		recomb = 'UNK'
		chain = ''
	
	return geneusage,scores,locus,chain,recomb

def ParseAlignment(alignment):
	#only consider the first alignment 
	alignment = alignment.split(';')[0] 
	#different fields of alignmetn are separated by '|'
	sub_strings = alignment.split('|')
	germ_start = sub_strings[0]
	germ_end = sub_strings[1]
	query_start = sub_strings[3]
	query_end = sub_strings[4]
	algn_len = sub_strings[2]
	alignment_string = sub_strings[5]
	num_mismatch = alignment_string.count('S')
	num_del = alignment_string.count('D')
	num_ins = alignment_string.count('I')
	shm = (num_mismatch+num_del+num_ins)/float(algn_len)
	return [query_start,query_end,germ_start,germ_end,algn_len,num_mismatch,num_ins,num_del,shm,alignment_string]


	
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

def guess_strand(nt_ab_seq,raw_seq,repeat=0):
	if nt_ab_seq in raw_seq:
		return '+'
	
	rc_seq = str(Seq(raw_seq,generic_dna).reverse_complement())		
	if nt_ab_seq in rc_seq:
		return '-'
	
	if 'N' in raw_seq:				
		forward=raw_seq.replace('N','[ATGCN]')
		mforward=re.match(r"%s"%forward,nt_ab_seq)				
		if mforward:
			return '+'
		reverse=rc_seq.replace('N','[ATGCN]')
		mreverse=re.match(r"%s"%reverse,nt_ab_seq)				
		if mreverse:
			return '-'
	else:
		if repeat>0:
			return ''
		else:
			return guess_strand(nt_ab_seq[10:-10],raw_seq,1)		


def match_sequence(input_seq,mixcr_seq):
	#try to match the input_seq to the sequence in the mixcr alignment file 
	if input_seq == mixcr_seq:
		#found the forward sequence
		return [True,'+']		
	#try comparing the reverse complement
	rc_seq = str(Seq(input_seq,generic_dna).reverse_complement())
	if rc_seq == mixcr_seq:
		#found the reverse complement
		return [True,'-']	
	if 'N' in input_seq:
		#this sequence has 'N's in it so mixcr coudl be auto-correcting the base		
		forwardcompseq=input_seq.replace('N','[ATGCN]')
		mforward=re.match(r"%s"%forwardcompseq,mixcr_seq)				
		if mforward:
			return [True,'+']		
		reversecompseq=rc_seq.replace('N','[ACTGN]')
		mreverse=re.match(r"%s"%reversecompseq,mixcr_seq)		
		if mreverse:
			return [True,'-']
		else:
			return [False,'']															
	else:
		#this sequence is not in the mixcr file 
		return [False,'']		

def GetFullAA(content,missing_fields):
	#our productiviety rules: 
		#if stop codon => no 
		#if there are no stop codons but we did not find all of the annotated regions or there is a potential indel, we say maybe 
		#all other situations, we say yes
	fields_to_concatenate = ['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']# ['AA. seq. FR1','AA. seq. CDR1','AA. seq. FR2','AA. seq. CDR2','AA. seq. FR3','AA. seq. CDR3','AA. seq. FR4']
	
	if not content['5_Prime_Annotation']:
		return ['','NO']
	
	starting_translation = fields_to_concatenate.index(content['5_Prime_Annotation'])
	aa = ''
	for f in range(starting_translation,len(fields_to_concatenate)):
		aa_region = content['AA. seq. '+fields_to_concatenate[f]]
		if not aa_region:			
			#once we hit a 'non-annotated' region of the antibody, stop adding to the amino acid sequence. For example if they return only FR1, FR2, FR4 then return FR1,FR2
			break		
		aa+=aa_region
	
	if content['5_Prime_Annotation'] == 'FR1' and content['3_Prime_Annotation'] == 'FR4':
		all_found = True
	else:
		all_found = False
	
	if '*' in aa:
		return [aa,'NO']
	if missing_fields or '_' in aa or all_found == False:
		return [aa,'MAYBE']	
	return [aa,'YES']	
	

def return_full_nt(content):
	all_regions = ['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']
	five_prime =''
	three_prime=''
	nt_seq = ''
	for region in all_regions:
		if content['N. Seq. '+region]:
			five_prime=region
			break
	
	if five_prime:
		for region in reversed(all_regions):#['FR4','CDR3','FR3','CDR2','FR2','CDR1','FR1']			
			if content['N. Seq. '+region]:
				three_prime=region
				break
			if region==five_prime:
				break
	missing_fields=False
	if five_prime and three_prime:
		index_pos_five = all_regions.index(five_prime)
		index_pos_three = all_regions.index(three_prime)		
		for inbetween in range(index_pos_five,index_pos_three+1):		
			region = content['N. Seq. '+all_regions[inbetween]]
			
			if not region:				
				missing_fields=True
				break
			nt_seq+=region
	else:
		missing_fields=True
	
	return [nt_seq,five_prime,three_prime,missing_fields]
			
		
		
	
	

#column names for output file 
presetlabels=['Seqheader',idIdentifier,	'Sequence',
	'Strand corrected sequence','Locus','Orientation',
	'FirstVgene','FirstDgene','FirstJgene',
	'Read(s) sequence','Read(s) sequence qualities','Recombination Type',
	'Chain','5_Prime_Annotation','3_Prime_Annotation','Productivity','Full length',
	'CDR3_Junction_In_Frame','Full AA',	'Full NT','All V hits','All V scores',
	'All D hits','All D scores','All J hits',
	'All J scores','All C hits','All C scores',
	'VGENE: Shm.nt','VGENE: Shm.per','JGENE: Shm.nt',
	'JGENE: Shm.per','N. Seq. FR1','Min. qual. FR1',
	'N. Seq. CDR1','Min. qual. CDR1','N. Seq. FR2',
	'Min. qual. FR2','N. Seq. CDR2','Min. qual. CDR2',
	'N. Seq. FR3','Min. qual. FR3','N. Seq. CDR3',
	'Min. qual. CDR3','N. Seq. FR4','Min. qual. FR4',
	'AA. seq. FR1','AA. seq. CDR1','AA. seq. FR2',
	'AA. seq. CDR2','AA. seq. FR3','AA. seq. CDR3',
	'AA. seq. FR4','AB start','AB end',
	'VGENE: Query start','VGENE: Query end','VGENE: Germline start',
	'VGENE: Germline end', 'VGENE: Alignment length', 'VGENE: Mismatch',
	'VGENE: Insertion', 'VGENE: Deletion','VGENE: Alignment',
	'JGENE: Query start','JGENE: Query end','JGENE: Germline start',
	'JGENE: Germline end', 'JGENE: Alignment length', 'JGENE: Mismatch',
	'JGENE: Insertion', 'JGENE: Deletion','JGENE: Alignment',
	'All D alignment','All C alignment','Command','Notes']			

def parseMIXCR(originalfileloc,resultfileloc,inputype,outfile=None,header_var='document_header',sequence_var='sequence',command_val = {}):
								
	command_string = json.dumps(command_val) if command_val else json.dumps({'MIXCR V1.3': 'Unknown settings'})
	if not outfile:
		outfile = "%s-parsed.annotation"%resultfileloc
	
	print('Parsing mixcr file')
	number_of_annotation_lines = useful.file_line_count(resultfileloc)	
	
	seqfile=immunogrepFile(originalfileloc,inputype) #the original file used as an input file for mixcr annotation  	
	iffile=immunogrepFile(resultfileloc,'TAB',None)#,"\t",True,"r") #the mixcr generated alignment file 
	
	parent_folder = '/'.join(resultfileloc.split('/')[:-1])+'/'
	error_file = open(resultfileloc+'.errorlog.txt','w')
	unfound_seqs = open(resultfileloc+'.notfound.txt','w') 
	notfound=0
	seq_num=0
	errors=0
	needcapture = True
	
	looper = useful.LoopStatusGen(number_of_annotation_lines,10)
	t1 = time.time()
	with open(outfile,"w") as f:
		f.write(descriptor_symbol+json.dumps(DatabaseTranslator())+'\n')#write a translator line to this file so that we know how to add results to database 
		f.write('\t'.join(presetlabels)+'\n')		
		
		#read each input sequence/file
		for fastseq in seqfile.read():			
			try:
				content={}
				if not fastseq:
					continue
				
				#read in the annotation information from mixcr
				seq = fastseq[sequence_var].upper()
				#extract sequence header and the SEQ_ID field from input file 
				[header,id] = GetHeaderInfo(fastseq,header_var)
				
				if needcapture:
					#we need to match this sequence to mixcr program output
					if iffile.IFclass.eof:
						mixcr_data = None
					else:
						mixcr_data = iffile.IFclass.read()
						#print percent status completed
						looper.next()
				
				#check whether mixcr data matches the current sequence 
				strand=''
				if mixcr_data:
					if 'Read id' in mixcr_data:
						if int(mixcr_data['Read id']) == seq_num:
							matched_seqs = True
						else:
							matched_seqs=False							
					elif 'Description R1' in mixcr_data:
						if mixcr_data['Description R1'].strip() == fastseq[header_var].strip():
							matched_seqs=True							
						else:
							matched_seqs = False							
					else:
						mixcr_data['Read(s) sequence'] = mixcr_data['Read(s) sequence'].upper()				
						[matched_seqs,strand] = match_sequence(seq,mixcr_data['Read(s) sequence'])					
					mixcr_seq = mixcr_data['Read(s) sequence']
				else:
					mixcr_seq = ''
					matched_seqs=False
					strand=''
					needcapture=True
								
				if matched_seqs==False:
					#these results did not match mixcr sequence, so this sequence probably did not yield any results
					#so we do not need to recapture a new miseq sequence. We will just stay with this one 
					needcapture = False
					content['Sequence']=seq
					content['Seqheader']=header
					content['Notes'] = 'Sequence not found in mixcr file;'
					content[idIdentifier] = id
					unfound_seqs.write('\t'.join([content['Seqheader'],content['Sequence'],mixcr_seq])+'\n')														
					content['Command'] = command_string
					content = defaultdict(str,content)
					output_line = [str(content[lab]) for lab in presetlabels]
					f.write('\t'.join(output_line)+'\n')					
					notfound+=1
					seq_num+=1		
					continue
				
				seq_num+=1		
				#in the next iteration of the code, we will need to get a fresh mixcr result
				needcapture=True						
				content = mixcr_data
				content['Notes'] = ''
				content[idIdentifier] = id
				content['Seqheader'] = header
				r_j = ''
				r_v = ''
				chain_v = ''
				
				content['Sequence']=seq
				content['Strand corrected sequence'] = content['Read(s) sequence']				
				[content['Full NT'],content['5_Prime_Annotation'],content['3_Prime_Annotation'],missing_fields]=return_full_nt(content)				
				if missing_fields:
					content['Notes']+='The sequence is missing features between the 5 prime and 3 prime region;'
					content['3_Prime_Annotation']=content['3_Prime_Annotation']+'*'				
					content['Full length'] = 'FALSE'
				else:
					if content['5_Prime_Annotation'] == 'FR1' and content['3_Prime_Annotation'] == 'FR4':						
						content['Full length'] = 'TRUE'
					else:
						content['Full length'] = 'FALSE'
					
				[content['Full AA'],content['Productivity']] =GetFullAA(content,missing_fields)
				if content['AA. seq. CDR3'] and content['AA. seq. CDR3'] in content['Full AA']:
					content['CDR3_Junction_In_Frame']= 'TRUE'
				else:
					content['CDR3_Junction_In_Frame']= 'FALSE'
					
				if content['All V hits']:
					[vgenelist,vscorelist,vlocus,chain_v,r_v]=extractScores(content['All V hits'])
					content['All V hits']=','.join(vgenelist)
					content['All V scores']=','.join(vscorelist)
					content['FirstVgene']=vgenelist[0]
					content['Locus']=vlocus
																	
				else:
					content['All V hits']=''
					content['All V scores']=''
					content['FirstVgene']=''
					content['Locus']=''
					
				if content['All D hits']:
					[dgenelist,dscorelist,dlocus,chain,recomb]=extractScores(content['All D hits'])
					content['All D hits']=','.join(dgenelist)
					content['All D scores']=','.join(dscorelist)
					content['FirstDgene']=dgenelist[0]
				else:
					content['All D hits']=''
					content['All D scores']=''
					content['FirstDgene']=''
					
					
				if content['All J hits']:
					[jgenelist,jscorelist,jlocus,chain,r_j]=extractScores(content['All J hits'])
					content['All J hits']=','.join(jgenelist)
					content['All J scores']=','.join(jscorelist)
					content['FirstJgene']=jgenelist[0]
				else:
					r_j = r_v
					content['All J hits']=''
					content['All J scores']=''
					content['FirstJgene']=''
					
				if content['All C hits']:
					[cgenelist,cscorelist,clocus,chain,recomb]=extractScores(content['All C hits'])
					content['All C hits']=','.join(cgenelist)
					content['All C scores']=','.join(cscorelist)
				else:
					content['All C hits']=''
					content['All C scores']=''
				if r_j == r_v:
					content['Recombination Type'] = r_v
					content['Chain'] = chain_v
				else:
					content['Recombination Type'] = ''
					content['Chain'] = ''
				
				if content['All V alignment']:
					[query_start,query_end,germ_start,germ_end,algn_len,num_mismatch,num_ins,num_del,shm,alignment_string] = ParseAlignment(content['All V alignment'])
					content['VGENE: Query start'] = query_start
					content['VGENE: Query end'] = query_end
					content['VGENE: Germline start'] = germ_start
					content['VGENE: Germline end'] = germ_end
					content['VGENE: Shm.nt'] = num_ins+num_del+num_mismatch
					content['VGENE: Mismatch'] = num_mismatch
					content['VGENE: Insertion'] = num_ins
					content['VGENE: Deletion'] = num_del
					content['VGENE: Alignment'] = alignment_string
					content['VGENE: Shm.per'] = round(100*shm,3)
					content['VGENE: Alignment length'] = algn_len						
					content['AB end'] = query_end
					content['AB start'] = query_start
					
				if content['All J alignment']:
					[query_start,query_end,germ_start,germ_end,algn_len,num_mismatch,num_ins,num_del,shm,alignment_string] = ParseAlignment(content['All J alignment'])
					content['JGENE: Query start'] = query_start
					content['JGENE: Query end'] = query_end
					content['JGENE: Germline start'] = germ_start
					content['JGENE: Germline end'] = germ_end
					content['JGENE: Shm.nt'] = num_ins+num_del+num_mismatch
					content['JGENE: Mismatch'] = num_mismatch
					content['JGENE: Insertion'] = num_ins
					content['JGENE: Deletion'] = num_del
					content['JGENE: Alignment'] = alignment_string
					content['JGENE: Shm.per'] = round(100*shm,3)
					content['JGENE: Alignment length'] = algn_len
					content['AB end'] = query_end
					if 'AB start' not in content:
						content['AB start'] = query_start
								
				content['Orientation'] = guess_strand(content['Full NT'],content['Sequence'])
				content['Command'] = command_string
				content = defaultdict(str,content)
				output_line = [str(content[lab]) for lab in presetlabels]
				f.write('\t'.join(output_line)+'\n')
																				
			except Exception as e:
				errors+=1
				print('There was an error in sequence: '+str(seq_num))
				print('Error: '+str(e))				
				error_file.write('****ERROR FOUND IN SEQUENCE:{0}  ****\n'.format(str(seq_num)))
				error_file.write(useful.print_error_string(e)+'\n')
				error_file.write('MIXCR DATA: \n')
				error_file.write(json.dumps(content,indent=4)+'\n')
				error_file.write('*************END OF ERROR*********\n')
														
					
	iffile.IFclass.close()
	seqfile.IFclass.close()
	error_file.close()
	unfound_seqs.close()
	
	if errors==0:
		os.remove(resultfileloc+'.errorlog.txt')
	
	if notfound==0:
		os.remove(resultfileloc+'.notfound.txt')
	t2 =time.time()
	print(str(t2-t1))
	return outfile
	
