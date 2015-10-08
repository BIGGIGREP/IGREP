#!/usr/bin/python
# This is the tools necessary for fetching motif positions and probability from MongoDB
# The MongoDB collection is Motifs
# Require Spee as TaxID (for the species) and Loc as Locus
# Return Left Flank position, Right Flank position, Left Flank probability table, Right Flank probability table,... 
# ...Common LMotif, Common RMotif
# Note: to be eligible to be common motif, the aa needs to be at least a prob of 0.3
#---------------------------------------------------------------------------------------------------------------------
# This tool will also provide ways of generating the MotifTables file (timestamp.motifT)
#
# 072814 - Addition of function --- PWM() where it processes on a per sequence basis and reports best score, 
#          position of best score, and frame [0,5]
# See below for the PWM() function requirements
#
# 072914 - Addition of function --- FindCDR3() where it will process on a per sequence basis based on a list of given Motifs,
#          to report which motif suits best and report best score, CDR3NT, CDR3AA, and report all scores from all motifs
# See below for FindCDR3() function requirements

from collections import defaultdict
#from immunogrep_taxonomy_tools import convert2TaxID
#from immunogrep_query_germline_functions import connectToIgDatabase
#from immunogrep_query_germline_functions import CDR3MotifDB
import json
import time
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import prod

# FindCDR3() function will use a given DNA sequence and list of possible LMotif and RMotifs to delimit the CDR3 regions
# Note: the format of the motif is changed because the JQuery is reporting the motif selection in its own format
# Required inputs: FindCDR3(seq<>DNA,Motif<>[(LMotif1,RMotif1,VDJ),(LMotif2,RMotif2,VJ)],suggest_chain<>VDJ/VJ/None,start_pos<>int,strand<>+/-,
#                  motif_type<>AA/NT,Ltrim<>the trim position for left motif,Rtrim<>the trim position for right motif)
# Note: the spec is that the input seq will be DNA !!!!
def FindCDR3(seq,Motif,suggest_chain=None,start_pos=0,strand=None,motif_type='AA'):
	with open('scratch/motiffiles.json','w') as w:
		w.write(json.dumps(Motif,indent=4))
	
	# initialization
	allscores=[]
	MaxP=1e-99*1e-99
	bestmotifset=0
	bestscoreset=()
	bestchain='NA'
	cdr3nt=''
	cdr3aa=''
	cdr3start=0
	cdr3end=0
	motifset=-1
	if suggest_chain !=None and not type(suggest_chain) is list:
		suggest_chain = [suggest_chain]
	for m in Motif:
		# map a motif set to variables
		motifset+=1
		chain=m[0]
		LMotif=m[1]
		RMotif=m[2]
		RMotifLength=int(sorted([int(k) for k in RMotif])[-1])
		Ltrim=m[3]
		Rtrim=m[4]
	
		if strand==None:
			frame=[0,1,2,3,4,5]
		elif strand=='+':
			frame=[0,1,2]
		elif strand=='-':
			frame=[3,4,5]

		# check for appropriate motif set, if suggest_chain provided then ignore those motif sets that aren't matching
		if suggest_chain!=None and suggest_chain !=[] and chain not in suggest_chain:
			scoreset=('N/A','-bypassed','by','user','!')
			allscores.append(scoreset)
			continue

		# scanning for best left flank
		[bestLstart,bestLend,Lframe,LMaxP]=PWM(seq,LMotif,frame,start_pos,motif_type)

		# scanning for best right flank
		[bestRstart,bestRend,Rframe,RMaxP]=PWM(seq,RMotif,frame,start_pos,motif_type)

		
		# storing scoreset
		scoreset=(LMaxP*RMaxP,bestLstart,bestRend,Lframe,Rframe,chain)
		allscores.append(scoreset)

		# checking for best possible CDR3 regions for reporting
		if LMaxP*RMaxP > 1e-16 and LMaxP*RMaxP > MaxP and Lframe==Rframe:
			MaxP=LMaxP*RMaxP
			bestmotifset=motifset
			bestscoreset=scoreset
			bestchain=chain
			if Lframe>=3:
				seq=str(Seq(seq,generic_dna).reverse_complement())
			# Perform trimming here to report the CDR3
			cdr3ntseq=seq[scoreset[1]:scoreset[2]]
			cdr3start=scoreset[1]+(Ltrim-1)*3
			cdr3nt=cdr3ntseq[(Ltrim-1)*3:-1*(RMotifLength-Rtrim+1)*3]
			cdr3end=cdr3start+len(cdr3nt)-1
			cdr3aa=str(Seq(cdr3nt,generic_dna).translate())
			
	return bestmotifset,MaxP,cdr3start,cdr3end,cdr3nt,cdr3aa,bestchain,bestscoreset,allscores



# PWM() function will check a sequence for the starting position and the ending position of the best motif
# Required inputs: PWM(seq<>text,LMotif<>dict,frame<>[0,1,2,3,4,5],start_pos<>intDefault0,motif_type<>AA/NTDefaultAA)
# Note: the start_pos is referencing on the original sequence input
def PWM(seq,LMotif,frame=[0,1,2,3,4,5],start_pos=0,motif_type='AA',end_pos=None):
	# if end_pos is not provided, set it as the length of the sequence
	if end_pos==None:
		end_pos=len(seq)
	else:
		end_pos=int(end_pos)
	# isolate only the section that needs to be checked
	seq=seq[start_pos:end_pos+1]
	# initialization
	LMaxP=1e-99 # initial highest probability
	beststart=0 # initial start of motif position
	bestend=0 # initial end of motif position
	bestframe=0 # initial best frame
	lwlen=int(sorted([int(k) for k in LMotif])[-1])
	lpos=[int(k) for k in LMotif]
	aaseq=defaultdict(str)
	# Checking sequence and motif compatibility and proceed with the pwm
	if validateSeq(seq,chartype='AA') and motif_type=='NT':
		print "!!!! Critical Error: Cannot compare Amino Acids sequence to a Nucleotide motif !!!!"
		return
	elif validateSeq(seq,chartype='NT') and motif_type=='AA':
		# Going through all supplied frames stored in the frame list
		for f in frame:
			# Original strand
			if not f>2:
				aaseq[f]=str(Seq(seq[f:],generic_dna).translate())
				for pos in range(len(aaseq[f])):
					lw=aaseq[f][pos:pos+lwlen]
					if (len(lw)<lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob=[LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP=prod(lwprob)
					if LP>LMaxP:
						LMaxP=LP
						beststart=start_pos+f+(pos+1)*3
						bestend=start_pos+f+pos*3+lwlen*3
						bestframe=f
			# Reverse complement strand
			else:
				rseq=str(Seq(seq,generic_dna).reverse_complement())
				aaseq[f]=str(Seq(rseq[(f-3):],generic_dna).translate())
				for pos in range(len(aaseq[f])):
					lw=aaseq[f][pos:pos+lwlen]
					if (len(lw)<lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob=[LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP=prod(lwprob)
					if LP>LMaxP:
						LMaxP=LP
						beststart=start_pos+f-3+(pos+1)*3
						bestend=start_pos+f-3+pos*3+lwlen*3
						bestframe=f
		
		# ---- This section is now obsolete ---- #
		# tempbeststart=defaultdict(int)
		# Note: the AA position is trying its best to give the original frame AA position
		# tempbeststart['AA']=beststart
		# tempbeststart['NT']=int(beststart*3)
		# tempbestend=defaultdict(int)
		# tempbestend['AA']=bestend
		# tempbestend['NT']=int(bestend*3)
		# beststart=tempbeststart
		# bestend=tempbestend
		# ---- The above section is now obsolete ---- #
		
		return beststart,bestend,bestframe,LMaxP
	else:
		# Going through all supplied frames stored in the frame list
		for f in frame:
			if not f>2:
				for pos in range(len(seq[f:])):
					lw=seq[f:][pos:pos+lwlen]
					if (len(lw)<lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob=[LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP=prod(lwprob)
					if LP>LmaxP:
						LMaxP=LP
						beststart=start_pos+f+pos
						bestend=start_pos+f+pos+lwlen
						bestframe=f
			else:
				rseq=str(Seq(seq,generic_dna).reverse_complement())
				for pos in range(len(rseq[(f-3):])):
					lw=rseq[(f-3):][pos:pos+lwlen]
					if (len(lw)<lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob=[LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP=prod(lwprob)
					if LP>LmaxP:
						LMaxP=LP
						beststart=start_pos+f-3+pos
						bestend=start_pos+f-3+pos+lwlen
						bestframe=f
		return beststart,bestend,bestframe,LMaxP

# This function validates whether a sequence is pertaining to the type specified. i.e. NT or AA
# Note: This function cannot distinguish amino acids sequence containing only ATCG because it will return True in either NT or AA case
def validateSeq(seq,chartype='NT'):
	if chartype=='NT':
		charset=set('ATCGUN')
	elif chartype=='AA':
		charset=set('ACDEFGHIKLMNPQRSTVWYXZ*')
	else:
		print "I have no clue what type you are talking about. Either NT or AA is allowed."
		return
	diff=set(seq.upper())-charset
	return not diff

# Branch point for requesting inputs
# If DBboolean is True, get MotifTable from MongoDB, else get from local file location
def getMotifInput(DBboolean,species,locus,inputfilename=None):
	if DBboolean=='True':
		filename=motifDB2JSON(convert2TaxID(species),locus)
	else:
		if inputfilename==None:
			sys.exit("No custom Motif Table file was specified")
		filename=motifTXT2JSON(inputfilename,species,locus)
	return filename,species,locus


# Getting Motif Table from Database and output MotifTables file
def motifDB2JSON(species,locus):
	lpos,rpos,lmotif,rmotif,lcommotif,rcommotif,maxlpos,maxrpos=getMotifTable(species,locus)
	pattern=lcommotif+'(.*)'+rcommotif
	lwlen=int(lpos[-1])
	rwlen=int(rpos[-1])
	timestr = time.strftime("%Y%m%d-%H%M%S")
	fout=open('scratch/'+timestr+'.motifT',"w")
	print >>fout,json.dumps(lpos)
	print >>fout,json.dumps(rpos)
	print >>fout,json.dumps(lmotif)
	print >>fout,json.dumps(rmotif)
	print >>fout,json.dumps(lcommotif)
	print >>fout,json.dumps(rcommotif)
	print >>fout,json.dumps(maxlpos)
	print >>fout,json.dumps(maxrpos)
	fout.close()
	return 'scratch/'+timestr+'.motifT'
	
	
# Converting custom Motif Table from user provided location and output MotifTables file
def motifTXT2JSON(inputfilename,species,locus):
	# Processing the user input motif table
	Llist=[]
	Rlist=[]
	L=0
	R=0
	with open(inputfilename) as f:
		for line in f:
			line=line.replace("\r",'')
			line=line.replace("\n",'')
			if line[0]=="\t":
				continue
			if 'LPosition' in line:
				L=1
				R=0
				tmp=line.split("\t")
				pos=tmp[1:]
				pos=[x for x in pos if x != '']
				r=str()
			if 'RPosition' in line:
				L=0
				R=1
				tmp=line.split("\t")
				pos=tmp[1:]
				pos=[x for x in pos if x != '']
				r=str()
			if L==1 and 'LPosition' not in line:
				tmp=line.split("\t")
				aa=tmp[0]
				prob=tmp[1:]
				for i in range(len(pos)):
					r=aa+'_'+pos[i]+'_'+prob[i]
					Llist.append(r)
			if R==1 and 'RPosition' not in line:
				tmp=line.split("\t")
				aa=tmp[0]
				prob=tmp[1:]
				for i in range(len(pos)):
					r=aa+'_'+pos[i]+'_'+prob[i]
					Rlist.append(r)
	f.close()
	# Generating necessary outputs to follow the standards outputs needed for MotifTables file
	Lpos=[]
	Rpos=[]
	Lprob=defaultdict(dict)
	Rprob=defaultdict(dict)
	for i in Llist:
		tmp=i.split('_')
		if tmp[1] not in Lpos:
			Lpos.append(tmp[1])
		Lprob[tmp[1]][tmp[0]]=float(tmp[2])
	for i in Rlist:
		tmp=i.split('_')
		if tmp[1] not in Rpos:
			Rpos.append(tmp[1])
		Rprob[tmp[1]][tmp[0]]=float(tmp[2])
	LcomMotif=str()
	for i in range(int(Lpos[-1])):
		itera=str(i+1)
		maxprob=0.3
		cnt=0
		if (len(Lprob[itera])==0):
			cnt=100
		for k in Lprob[itera]:
			if Lprob[itera][k]>=maxprob:
				cnt+=1
				maxprob=Lprob[itera][k]
				maxaa=k
				if (maxaa=='C' or maxaa=='Y') and maxprob>0.81:
					maxlpos=int(itera)
		if cnt==1:
			LcomMotif=LcomMotif+maxaa
		else:
			LcomMotif=LcomMotif+'.{1}'
	RcomMotif=str()
	for i in range(int(Rpos[-1])):
		itera=str(i+1)
		maxprob=0.3
		cnt=0
		if (len(Rprob[itera])==0):
			cnt=100
		for k in Rprob[itera]:
			if Rprob[itera][k]>=maxprob:
				cnt+=1
				maxprob=Rprob[itera][k]
				maxaa=k
				if (maxaa=='W' or maxaa=='F') and maxprob>0.81:
					maxrpos=int(itera)
		if cnt==1:
			RcomMotif=RcomMotif+maxaa
		else:
			RcomMotif=RcomMotif+'.{1}'
	# Output standards to Motif Tables file
	timestr = time.strftime("%Y%m%d-%H%M%S")
	fout=open('scratch/'+timestr+'.motifT',"w")
	print >>fout,json.dumps(Lpos)
	print >>fout,json.dumps(Rpos)
	print >>fout,json.dumps(Lprob)
	print >>fout,json.dumps(Rprob)
	print >>fout,json.dumps(LcomMotif)
	print >>fout,json.dumps(RcomMotif)
	print >>fout,json.dumps(maxlpos)
	print >>fout,json.dumps(maxrpos)
	fout.close()
	return 'scratch/'+timestr+'.motifT'
	

# This function converts MotifTable .motifT file into variables
def convertMotifTable2Var(inputjson):
	with open(inputjson) as f:
		lines=[f.next() for i in range(8)]
		lpos=[int(i) for i in json.loads(lines[0])]
		rpos=[int(i) for i in json.loads(lines[1])]
		lmotif=json.loads(lines[2])
		rmotif=json.loads(lines[3])
		lcommotif=json.loads(lines[4])
		rcommotif=json.loads(lines[5])
		pattern=lcommotif+'.*'+rcommotif
		lwlen=int(lpos[-1])
		rwlen=int(rpos[-1])
		maxlpos=json.loads(lines[6])
		maxrpos=json.loads(lines[7])
	return lpos,rpos,lmotif,rmotif,lcommotif,rcommotif,pattern,lwlen,rwlen,maxlpos,maxrpos


# Get Motif Table from MongoDB database
# This function requires Species (Unique numbers) and Locus (IGH,IGL,IGK)
def getMotifTable(spee,loc):
	db=connectToIgDatabase()
	motif=db['motifs']
	prob=motif.find_one({u'TaxID':spee,u'Locus':loc})
	Lpos=[]
	Rpos=[]
	Lprob=defaultdict(dict)
	Rprob=defaultdict(dict)
	# Storing the position and the probability tables into dictionaries
	try:
		prob['Lmotif']
	except:
		sys.exit("!!!! Warning: Data does not exist on MongoDB !!!!")

	for i in prob['Lmotif']:
		tmp=i.split('_')
		if tmp[1] not in Lpos:
			Lpos.append(tmp[1])
		Lprob[tmp[1]][tmp[0]]=float(tmp[2])
	for i in prob['Rmotif']:
		tmp=i.split('_')
		if tmp[1] not in Rpos:
			Rpos.append(tmp[1])
		Rprob[tmp[1]][tmp[0]]=float(tmp[2])
	# Finding the most common Left flanking and Right flanking motifs
	LcomMotif=str()
	tempLprob=Lprob.copy()
	for i in range(int(Lpos[-1])):
		itera=str(i+1)
		maxprob=0.3
		cnt=0
		if (len(tempLprob[itera])==0):
			cnt=100
		for k in tempLprob[itera]:
			if tempLprob[itera][k]>=maxprob:
				cnt+=1
				maxprob=tempLprob[itera][k]
				maxaa=k
				if (maxaa=='C' or maxaa=='Y') and maxprob>0.81:
					maxlpos=int(itera)
		if cnt==1:
			LcomMotif=LcomMotif+maxaa
		else:
			LcomMotif=LcomMotif+'.{1}'
	RcomMotif=str()
	tempRprob=Rprob.copy()
	for i in range(int(Rpos[-1])):
		itera=str(i+1)
		maxprob=0.3
		cnt=0
		if (len(tempRprob[itera])==0):
			cnt=100
		for k in tempRprob[itera]:
			if tempRprob[itera][k]>=maxprob:
				cnt+=1
				maxprob=tempRprob[itera][k]
				maxaa=k
				if (maxaa=='W' or maxaa=='F') and maxprob>0.81:
					maxrpos=int(itera)
		if cnt==1:
			RcomMotif=RcomMotif+maxaa
		else:
			RcomMotif=RcomMotif+'.{1}'
	
	return Lpos,Rpos,Lprob,Rprob,LcomMotif,RcomMotif,maxlpos,maxrpos


# This function will fetch all unique species to a list
def fetchSpeciesList():
	db=CDR3MotifDB()
	return db.GetUniqSpecies()
