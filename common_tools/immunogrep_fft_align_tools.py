import numpy as np 
import pyfftw
import math
import matplotlib.pyplot as plt
import re
import immunogrep_useful_functions as useful 

#this lists a series of tools for performing gapless alignment using FFT

#function for converting sequence into complex number. default value is 1.
def seq_convert(sequence,method=1,n=0):
	
	seqLen = len(sequence)
	if n==0:
		vecLen = seqLen
	else:
		vecLen=n
		
	if method == 1:		
		intSeq = pyfftw.n_byte_align_empty(vecLen, 1, 'complex')
		#intSeq = numpy.empty(seqLen,complex);
		for i,let in enumerate(sequence):
			if i<vecLen:
				if let == 'A' or let=='a':
					intSeq[i] = 1j
				elif let=='T' or let == 't':
					intSeq[i]=-1j
				elif let=='C' or let == 'c':
					intSeq[i] = 1				
				elif let =='G' or let == 'g':
					intSeq[i] = -1
				else:
					intSeq[i] = 0
		for i in range(seqLen,vecLen):
			intSeq[i] = 0
	return intSeq
	
#function for converting complex series back int a sequence
def complex_convert(complexN,method=1):
	seqLen = len(complexN)
	
	if method==1:
		seq=""
		for i in range(seqLen):
			if complexN[i]==1j:
				seq=seq+'a'
			elif complexN[i]==-1j:
				seq =seq+'t'
			elif complexN[i]==1:
				seq=seq+'c'
			elif complexN[i]==-1:
				seq=seq+'g'
			else:
				seq=seq+'n'				
	return seq
	

####PRIVATE FUNCTIONS: MOST PEOPLE WILL NOT HAVE TO WORRY ABOUT CALLING OR USING THESE FUNCTIONS##########
#create fft plan
def createFFTplan(algnLen):
	##the following steps will ensure that the fft is fast and efficient
	#forward transform plan THIS IS REQUIRED
	f_signal = pyfftw.n_byte_align_empty(algnLen,1, 'complex')
	fwdPlan = createplan(f_signal,algnLen)		
	return f_signal,fwdPlan

#create reverse (ifft) plan
def createInvPlan(algnLen):
	#reverse transform plan THIS IS REQUIRED
	r_signal = pyfftw.n_byte_align_empty(algnLen,1, 'complex')
	revPlan = createplan(r_signal,algnLen,2)
	## end of optimization	
	return r_signal,revPlan

	
#USE THIS FUNCTION TO CREATE THE PLAN NECESSARY FOR STORING THE PROPER PLAN/PATH FOR COMPUTING FOURIER TRANSFORMS
def createplan(s,n,dirF=1):	
	b = pyfftw.n_byte_align_empty(n,1, 'complex')
	if dirF==1:
		fftwO = pyfftw.FFTW(s,b,flags=('FFTW_PATIENT',))
	else:
		fftwO = pyfftw.FFTW(s,b,flags=('FFTW_PATIENT',),direction="FFTW_BACKWARD")	
	return fftwO
	
#function for performing cross corerelation in fourier space
	#input 1=> fourier series representing query DNA sequence
	#input 2=> fourier series representing germline DNA sequence
	#input 3=> required fft plan for finding transform
	#input 4=> required varaible name needed for computing inverse (this is what we will invert)
def return_cross_corr(complex_query,complex_target,fftPlan,inverseVariable):

	cross_corr_real = pyfftw.n_byte_align_empty(len(complex_query),1, 'complex')
	#find complex conjugate of germline		
	target_conj = np.conjugate(complex_target)			
	#perform corss correlation
	inverseVariable[:] = complex_query*target_conj				
	#take inverse	
	cross_corr_real[:] = fftPlan()		

	cross_corr_real = np.round(np.real(cross_corr_real))
	max_val = np.amax(cross_corr_real) # find max hit
	pos_val = np.argwhere(cross_corr_real==max_val) # find where hits occur	
	return max_val,pos_val,cross_corr_real # return the max valu,e the the max positions, the the vector of alignments

#s complex series representing DNA
def convert_complex_four(s,fftPlan,inverseVariable):
	#if n is smaller than length of array s, array is cropped
	#if n is larger than lengnth of array s, then array is padded with 0	
	inverseVariable[:] = s
	four_series = fftPlan()		
	return four_series
####################################################


#####PUBLIC FUNCTIONS: USE THESE FUNCTIONS IN OTHER SCRIPTS THAT WANT TO TAKE ADVANTAGE OF THIS GROUP OF FUNCTIONS########
#general method to convert a DNA sequence into Fourier series
def seq_to_four(seq1,fftPlan,inverseVariable,method=1,n=0):
	
	if n == 0:
		n = len(seq1)

	sF1 = pyfftw.n_byte_align_empty(n,1,'complex') #create an empty vector
	s1 = seq_convert(seq1,method,n);			
	sF1[:] = convert_complex_four(s1,fftPlan,inverseVariable);
	return sF1

#returns a list of sequences into a list of fourier/complex series
#this is usuefull for database of sequences
def seqDatab_to_four(seqList,fftPlan,inverseVariable,method=1,n=0):		
	four_seq_datab = [];		
	for i in range(len(seqList)):
		if n==0:
			four_seq_datab.append(pyfftw.n_byte_align_empty(len(seqList[i]),1,'complex'));
		else:
			four_seq_datab.append(pyfftw.n_byte_align_empty(n,1,'complex'));					
		four_seq_datab[i][:] = seq_to_four(seqList[i],fftPlan,inverseVariable,method,n)
	return four_seq_datab
	
#use this function to just align two sequences together
#seq1 = input DNA sequence one
#seq2 = input DNA sequence two
#method for converting dna to integer
#this will not take advantage of speed because it makes a planner each time 
def pair_align_fft(seq1,seq2,method=1):
	myLen = len(seq1)+len(seq2)
	
	#lets make a planner
	[g,fwdPlan] = createFFTplan(myLen)
	[gR,revPlan] = createInvPlan(myLen)	
	
	#convert both seqs to fourier	
	sF1 = seq_to_four(seq1,fwdPlan,g,1,myLen);# convert_complex_four(s1,fwdPlan,g);
	sF2 = seq_to_four(seq2,fwdPlan,g,1,myLen);#convert_complex_four(s2,fwdPlan,g);
	
	#align sequences and get result
	[max_val,max_pos,algn] = return_cross_corr(sF1,sF2,revPlan,gR);#,cross_corr_result);
	return max_val,max_pos,algn

#model1 => this defines how the nucleotide sequence that is not shifting is
#currently aligned.  it lists all the positions where a dna sequence
#exists.  For example if this sequence is put first as [ACTG.....] WHERE
#"." represents no sequence or "0 padding" then model1 would be
#[1,2,3,4,0,0,0,0,0] where 4 represents position 4 or nucleotide 4 along
#sequence

#model 2=> this defines how the nucleotide sequence that is shifting along
#the first one is aligned/designed.  For example if model 2 is [.....ACTG]
#where '." represents regions of no sequence or "0" padding, then we right
#it as [0,0,0,0,0,1,2,3,4] 		

#the output is as follows [shift position in fft alignment (index of array), alignment start along query, alignment end along query, alignment start along target, alignment query along query, alignment length)
def modelOverlap(model1,model2):
	model2shift = [0]*len(model2)
	algnShift = []
		
	for shift in range(len(model1)):
		for i in range(len(model2)):
			shiftPos = (i-shift)%len(model2)			
			model2shift[i] = model2[shiftPos]
	
		
		
		algn = [0]*len(model2)
		
		for k in range(len(model1)):
			algn[k] = model1[k]*model2shift[k]

	
		
		currentPos = 0
		
		start = -1	
		batch = 0
		s = 1
		maxLen = 0
		maxP = -1
		endP=-1
		
		posS = [];#[0]*4
		tempPos = [0]*3
		for s in range(len(algn)):			
			if algn[s]!=0:
				if start == -1:
					start =s 
			else:
				if start !=-1:
					endP=s-1	
					tempPos=[start,endP,endP-start+1]					
					posS.append(tempPos)
					
					if posS[batch][2]>maxLen:
						maxP = batch
						maxLen = posS[batch][2]					
					start = -1
					tempPos = [0]*4
					endP = -1
					batch = batch+1
					
					
		if algn[s] != 0:
			if start == -1:
				start = s
				endP = s
			else:
				endP=s
			tempPos[:2] =[start,endP,endP-start+1]
			posS.append(tempPos) 
						
			if posS[batch][2]>maxLen:
				maxP = batch;
				maxLen = posS[batch][2]
			batch =batch+1
			
		if maxP==-1:
			algnInd=[]
		else:
			algnInd = [posS[maxP][0],posS[maxP][1]]
	
		
		if algnInd!=[]:
			algnShift.append([shift,model1[algnInd[0]]-1,model1[algnInd[1]]-1,model2shift[algnInd[0]]-1,model2shift[algnInd[1]]-1,posS[maxP][2]])
		else:
			algnShift.append([shift,-1,-1,-1,-1,0])
			
	return algnShift
		

#model1 => this defines how the nucleotide sequence that is not shifting is
#currently aligned.  it lists all the positions where a dna sequence
#exists.  For example if this sequence is put first as [ACTG.....] WHERE
#"." represents no sequence or "0 padding" then model1 would be
#[1,2,3,4,0,0,0,0,0] where 4 represents position 4 or nucleotide 4 along
#sequence

#model 2=> this defines how the nucleotide sequence that is shifting along
#the first one is aligned/designed.  For example if model 2 is [.....ACTG]
#where '." represents regions of no sequence or "0" padding, then we right
#it as [0,0,0,0,0,1,2,3,4] 		

#the output is as follows [shift position in fft alignment (index of array), alignment start along query, alignment end along query, alignment start along target, alignment query along query, alignment length)
def AlignmentOverlap(model1,model2,shift):
	model2shift = [0]*len(model2)
	algnShift = []
		
	
	for i in range(len(model2)):
		shiftPos = (i-shift)%len(model2)			
		model2shift[i] = model2[shiftPos]

	
	
	algn = [0]*len(model2)
	
	for k in range(len(model1)):
		algn[k] = model1[k]*model2shift[k]


	
	currentPos = 0
	
	start = -1	
	batch = 0
	s = 1
	maxLen = 0
	maxP = -1
	endP=-1
	
	posS = [];#[0]*4
	tempPos = [0]*3
	for s in range(len(algn)):			
		if algn[s]!=0:
			if start == -1:
				start =s 
		else:
			if start !=-1:
				endP=s-1	
				tempPos=[start,endP,endP-start+1]					
				posS.append(tempPos)
				
				if posS[batch][2]>maxLen:
					maxP = batch
					maxLen = posS[batch][2]					
				start = -1
				tempPos = [0]*4
				endP = -1
				batch = batch+1
				
				
	if algn[s] != 0:
		if start == -1:
			start = s
			endP = s
		else:
			endP=s
		tempPos[:,2] =[start,endP,endP-start+1]
		posS.append(tempPos) 
					
		if posS[batch][2]>maxLen:
			maxP = batch;
			maxLen = posS[batch][2]
		batch =batch+1
		
	if maxP==-1:
		algnInd=[]
	else:
		algnInd = [posS[maxP][0],posS[maxP][1]]


	if algnInd!=[]:
		
		algnShift= [model1[algnInd[0]]-1,model1[algnInd[1]]-1,model2shift[algnInd[0]]-1,model2shift[algnInd[1]]-1,posS[maxP][2]]
	else:
		algnShift=[-1,-1,-1,-1,0]
			
	return algnShift
def maxBarcodeLen(barcodeSeqList):
	maxLen = 0
	for bs in barcodeSeqList:				
		if len(bs[0])>maxLen:
			maxLen=len(bs[0])
	
	return maxLen
#determine if a number n is divisble by prime numbers in list 
def IsDivisible(n, primes = [] ):

	for p in primes:
		while n%p==0:
			n = n/p
	
	if n==1:
		return True 
	else:
		return False
		
		

def GetIdealAlgnLen(algnLen):
	#ONLY allow algnLen that is multiple of 2,3,5,7
	possible_lens = [32,35, 36, 40, 42, 45, 48, 49, 50, 54, 56, 60, 63, 64, 70, 72,
					75, 80, 81, 84, 90, 96, 98,64,100, 105, 108, 112, 120, 125, 126, 128, 135, 140, 144, 147, 
					150, 160, 162, 168, 175, 180, 189, 192, 196, 200, 210, 216, 224, 
					225, 240, 243, 245, 250, 252, 256, 270, 280, 288, 294, 300, 
					315, 320, 324, 336, 343, 350, 360, 375, 378, 384, 392, 400,
					405, 420, 432, 441, 448, 450, 480, 486, 490, 500, 504, 512,
					525, 540, 560, 567, 576, 588, 600, 625, 630, 640, 648, 672, 
					675, 686, 700, 720, 729, 735, 750, 756, 768, 784, 800, 810, 840, 
					864, 875, 882, 896, 900, 945, 960, 972, 980,1000, 1008, 1024, 1029, 
					1050, 1080, 1120, 1125, 1134, 1152, 1176, 1200, 1215, 1225, 1250, 1260, 
					1280, 1296, 1323, 1344, 1350, 1372, 1400, 1440, 1458, 1470, 1500, 1512, 
					1536, 1568, 1575, 1600, 1620, 1680, 1701, 1715, 1728, 1750, 1764, 1792, 
					1800, 1875, 1890, 1920, 1944, 1960]  
	
	if algnLen>possible_lens[-1]:
		print("Warning: ALIGNMENT LENGTH IS LARGER THAN HARDCODED POSSIBLE_LENS in function getIdealAlgnLen")
		return algnLen #worst case scenario 
	else:
		#choose the minimum length larger than algnLen 
		for l in possible_lens:
			if l>=algnLen:
				return l

#THESE NEXT LINES ARE ONLY IF INTERESTED IN USING FOURIER GAPLESS ALIGNMENT
#make a plan if you want to use FFT...makes it fast by planning the route
#OVERHEAD FUNCTION/PLANNING FOR THE FFT STEP IF NECESSARY -> functions in file nt_fft_align_tools								
#barcode_list = > list of lists. first list => 0 and 1. element 0 is forward, element 1 is reverse. 
#for each list inside, its the sequence of the barcode 
def MakeFFTStructure(n,barcode_list,maxLenSeq,maxLenBarcode,direction):
		
	[s,p] = createFFTplan(n)	#plan for fft	
	[r_s,r_p] = createInvPlan(n)	#plan for inverse fft	
	
	if direction == 0:
		fourBarCode = seqDatab_to_four(barcode_list[0],p,s,1,n) # fwdPlan,f_signal,1,algnLen)#convert seq database to fourier 		
		fourBarCodeRev = [0]*len(barcode_list[0]) #create an empty set of values for RC 
	elif direction == 1:		
		#generate FFT database using RC 		
		fourBarCodeRev = seqDatab_to_four(barcode_list[1],p,s,1,n)
		fourBarCode = [0]*len(barcode_list[1]) #create an empty set of values for RC 
	elif direction == 2:
		#search both directions of provided barcode		
		fourBarCodeRev = seqDatab_to_four(barcode_list[1],p,s,1,n)
		fourBarCode = seqDatab_to_four(barcode_list[0],p,s,1,n) # fwdPlan,f_signal,1,algnLen)#convert seq database to fourier 		
	else:
		fourBarCode = [0]*len(barcode_list) #create an empty set of values for RC 
		fourBarCodeRev = [0]*len(barcode_list) #create an empty set of values for RC 
	
	modelSeq = [0]*n #we use this to make an alignment model so we knwo what region of each sequences are aligned to one another
	modelSeq[0:maxLenSeq] = [r for r in range(1,maxLenSeq+1)]
	modelBarcode = [0]*n #we use this to make an alignment model so we knwo what region of each sequences are aligned to one another
	modelBarcode[0:maxLenBarcode] = [r for r in range(1,maxLenBarcode+1)]
	aM = modelOverlap(modelSeq,modelBarcode)
	
	all_barcodes_combined = [ fourBarCode,fourBarCodeRev ] 
	
	
	return [s,p,r_s,r_p,all_barcodes_combined,aM]


		
def compare_sequences(seq1,seq2):#barcode_uint8):
	#returns True if sequences are same, False if they are different	
	diff=0 #counter for nucleotide differences	
	
	# Make character array
	# Convert char to uint8 numbers
	a=np.array(seq1.upper(),'c').view(np.uint8)						
	b=np.array(seq2.upper(),'c').view(np.uint8)						
	# find differences simultaneously; if identical c=0
	c=a-b#barcode_uint8
	diff=len(np.nonzero(c)[0])	
	
	return diff 

		
#This is a class for alinging to barcodes using INEXACT GAPPLES MATCH 
#direction = 0,1, or 2 , 2=> search both directions
class BarcodeAligner(object):
	
	def __init__(self,barcodeSeqList,penalize_truncations,direction=0,allowed_mismatches_in_alignment=None,minimum_alignment_length=None,nmax=0,nmin=60):						
		self.barcodeseq_list = barcodeSeqList		
		self.penalize_truncations = penalize_truncations		
		
		self.nmax = nmax 
		self.nmin = nmin	
		self.nbarcodeLen = maxBarcodeLen(barcodeSeqList) #find length of max barcode  
		self.allowed_mismatches_in_alignment = self.nbarcodeLen if not allowed_mismatches_in_alignment else allowed_mismatches_in_alignment #if allowed_mismatches_in_alignment = None, then there is no upper limit on the number of allowed mismatches, it will just return the best alignment
		[self.barcodeseq,self.barcodeheader]=zip(*barcodeSeqList)	
		
		absolute_min_length = 3
		
		#get a minimum alignment length 
		if minimum_alignment_length:
			self.minimum_alignment_length = minimum_alignment_length
		else:
			self.minimum_alignment_length = min([len(s) for s in self.barcodeseq])-self.allowed_mismatches_in_alignment
		#if the final minimum length is too small, adjust it 
		if self.minimum_alignment_length<=absolute_min_length:
			self.minimum_alignment_length = absolute_min_length+1
			
		self.minimum_alignment_length = 0 if not minimum_alignment_length else minimum_alignment_length
		
		self.barcode_min_scores = []
		#in the alignment  base matches get a score of +1, complementary base mismatches get a score of -1 . Therefore the minimum allowed score in the alignment bst be len -2*allowed_mismatches_in_alignment for worst case scenario
		if penalize_truncations == True: #this means we will treat trucanted alignments as mismatches. Therefore the minimum_alignment_length does not matter because it is basically the number of allowed mismatches
			for i in self.barcodeseq:
				#maximum possible score of barcode = > length of barcode (100% alignment)
				#minimum possible score of barcode = > length of barcode - (2*allowed_mismatches_in_alignment)
				#any alignment below this score we will ignore 
				self.barcode_min_scores.append([len(i)-2*self.allowed_mismatches_in_alignment,len(i)])	
		else: #this means truncated alignments to dont get negative mismatch scores 
			for i in self.barcodeseq:
				#maximum possible score of barcode = > length of barcode (100% alignment)
				#minimum possible score of barcode = > MINIMUM ALLOWED LENGTH - (2*allowed_mismatches_in_alignment) => we do this because we allow truncatiosn in alignment so the score of 100% will not be length of barcode			
				self.barcode_min_scores.append([self.minimum_alignment_length-2*self.allowed_mismatches_in_alignment,len(i)])					
							
		if nmax==0:
			self.nmax = 400
		
		maxAlgnLen = GetIdealAlgnLen(self.nbarcodeLen+self.nmax) #this will be the fft plan for the longest sequence
		minAlgnLen = GetIdealAlgnLen(self.nbarcodeLen+self.nmin) #this will be the fft plan for the shorest sequence 
		midLen = (maxAlgnLen-minAlgnLen)/2+minAlgnLen
		
		if midLen>0:	
			midLen_ideal = GetIdealAlgnLen(self.nbarcodeLen+midLen)
			#a nx2 array. elemen 0 => ideal algn len for fft, eleme 1=> actual length of sequence 
			#sort results by actual length of fft plans
			lengths = sorted([[minAlgnLen,self.nmin],[midLen_ideal,midLen],[maxAlgnLen,self.nmax]], key = lambda x:x[0])
		else:
			lengths = sorted([[minAlgnLen,self.nmin],[maxAlgnLen,self.nmax]], key = lambda x:x[0])
		
		
		
		#WE do this because we want to account for both forward and reverse complements of barcode 
		if direction == 0:
			rc_list = ['']*len(self.barcodeseq) #its empty 				
			fd_list = self.barcodeseq
		elif direction == 1:		
			#convert barcodes to rc. ERASE the original 
			rc_list = [useful.Reverse_Complement(bc) for bc in self.barcodeseq]
			fd_list = ['']*len(self.barcodeseq)
		elif direction == 2:
			#keep both directions 
			fd_list = self.barcodeseq
			rc_list = [useful.Reverse_Complement(bc) for bc in self.barcodeseq]
		else:
			fd_list = ['']*len(self.barcodeseq)
			rc_list = ['']*len(self.barcodeseq) #its empty 				
		#update self.barcodeseq to 2D array 
		self.barcodeseq = [fd_list,rc_list]
			
		
		#for each of the lengths provided, make fft plans
		self.fftplans = []
		for each_len in lengths:								
			[a,b,c,d,e,f] = MakeFFTStructure(each_len[0],self.barcodeseq,each_len[1],self.nbarcodeLen,direction)
			self.fftplans.append({
				'f_signal':a,
				'f_plan':b,
				'r_signal':c,
				'r_plan':d,
				'barcode_fft':e,
				'alignmentModel':f,
				'n':each_len[0]
			})
				
		#now make a regular expression for each barcode (use this for exact match)
		self.compiledRegExp = [[],[]]	
		for j in range(2):
			for i in range(len(self.barcodeseq[0])):
				if self.barcodeseq[j][i]:
					self.compiledRegExp[j].append(re.compile(self.barcodeseq[j][i],re.IGNORECASE))
				else:
					self.compiledRegExp[j].append(None)
		
			
				
	
	
	def AlignToSeq(self,seqs):
		
					
		#all sequqences will correpsond to a sequenc eof interest we want to align to barcodes 
		#seqs will be array of either one or two elements 
			#the first element is the forward sequence 
			#the second element is the reverse complement if applicable 
				
		#first check for an exact match 				
		results = {}
		each_seq = seqs
		
		found_results = self.ExactMatch(each_seq)			
		
		if found_results:
			isotypes = []
			scores = []
			mis = []
			algnpos = []
			algnlens = []
			percent=[]
			direction = []
			
			for x in found_results:
				isotypes.append(x['Isotype'])
				mis.append(0)
				algnpos.append(x['AlgnPos'])
				algnlens.append(x['AlgnLen'])
				percent.append(100)
				scores.append(x['MaxScore'])
				direction.append(x['Direction'])
			results = {
				'Direction':direction,
				'MaxScore':scores,
				'Mismatches':mis,
				'Isotype':isotypes,
				'Alignment start position':algnpos,					
				'Percent similarity':percent,
				'Alignment length':algnlens
			}
			
		
		if results:
			return results
		
		#could not find an exact match, 
		#checked to see if we need to do inexact match 
		if self.allowed_mismatches_in_alignment==0 and self.penalize_truncations==True:
			return {}
		
		#revert to FFT match 			
		seq_len = len(seqs)		
		
		#find the FFT plan with the proper alignment length 
		good_plan = None
		for plans in self.fftplans:
			if plans['n']>=(seq_len+self.nbarcodeLen): #find the first plan whose length is >= seq len
				good_plan = plans
				break
		if not good_plan: #have not found the right plan yet 
			#thereforew, we need to make a new fft plan using the current length 				
			n = GetIdealAlgnLen(seq_len+self.nbarcodeLen)
			[a,b,c,d,e,f] = MakeFFTStructure(n,self.barcodeseq,seq_len,self.nbarcodeLen)
			self.fftplans.append({
				'f_signal':a,
				'f_plan':b,
				'r_signal':c,
				'r_plan':d,
				'barcode_fft':e,
				'alignmentModel':f,
				'n':n
			})
			
			good_plan = self.fftplans[-1]
						
		#now actually perform fft alignment using the current plan 
		
		#never consider small/truncated sequqences
		if len(each_seq)<self.minimum_alignment_length:
			return results #ignore sequences with low lengths 
		fftresults = self.FFTMatch(each_seq,good_plan)	
		if fftresults:
			isotypes = []
			mis = []
			scores = []				
			algnpos = []
			algnlens = []
			percent = []
			direction = []
			
			for x in fftresults:
				isotypes.append(x['Isotype'])
				mis.append(x['Mismatches'])
				scores.append(x['MaxScore'])
				algnpos.append(x['AlgnPos'])
				algnlens.append(x['AlgnLen'])
				percent.append(x['PercentSimilarity'])
				direction.append(x['Direction'])
			results = {
				'Direction':direction,
				'Mismatches':mis,
				'MaxScore':scores,
				'Isotype':isotypes,
				'Alignment start position':algnpos,					
				'Percent similarity':percent,
				'Alignment length':algnlens
			}
			
																									
		return results


	#GAPLESS SEARCH
	#this function will find inexact matches without accounting for gaps 
	#min_algn_score => this is the minimum possible alignment score based on user request. The equation is as follows => (min alignment length) -(2*allowed_mismatches)
	#based on minimum number of             
	#barcode_min_scores => list of lists, element 1 => minimum allowed score, elment 2 => barcode length 
	def FFTMatch(self,seqFwd,plan):				
		
		fourBarCode_list  = plan['barcode_fft']
		algnLen = plan['n']
		fwdPlan = plan['f_plan']
		revPlan = plan['r_plan']
		f_signal = plan['f_signal']
		r_signal = plan['r_signal']
		alignmentModel = plan['alignmentModel']
		
		#convert sequences into FFT (note the barcodes are already converted in var fourBarCode)
		fourSeq_fwd = seq_to_four(seqFwd,fwdPlan,f_signal,1,algnLen); #-> function in nt_fft_align_tools
			
		results = []
		barcode_ind = 0
		bestScore = 0
		
		seqFwdLen = len(seqFwd)
		
		#go through each barcode 
		for j in range(len(fourBarCode_list[0])):
			scores = [-100,-100]
			pos_results = [None]*2
			for pos in range(2):				
				barcode_ind = j 
				my_barcode = self.barcodeseq[pos][barcode_ind]
				current_barcode = fourBarCode_list[pos][j] 				
				#loop through each barcode and find best alignment	
				if not my_barcode:
					continue #we aren't considering this as an option for alignment 
				#algnDataF is structured as follows: [1] max alignment, [2] position of max alignment, [3] vector of alignments at each position
				algnDataF= return_cross_corr(fourSeq_fwd,current_barcode,revPlan,r_signal) #align fwd sequence to barcode seq (returns max alignment/all positions with thtat max alingment, a llist of alignemnts)								
						
				if algnDataF[0]<self.barcode_min_scores[barcode_ind][0]:
					continue #this barcode did not align to sequence well (or it did not pass the minimum alignmetn length)
				#the next sequence could potentially match the barcode because it is above the alignmetn filter 									
				
				#First get the length of the alignment		
				algn_pos = algnDataF[1][0][0] 
				#what is the alignmetn between sequence and barcode at algn_pos
				#alignmentResults[0] => should be teh same as algn_pos , just an index
				#alignmentResults[1] => Start of alignment with respect to sequence
				#alignmentResults[2] => End of alignmetn with respenct to sequence 
				#alignmentResults[3] => Start of alingment with respsect to barcode
				#alignmentResults[4] => End of alignment with respect to barcode 
				alignmentResults = alignmentModel[algn_pos] 			
				seqAlgnStart = alignmentResults[1]
				seqAlgnEnd = alignmentResults[2]
				barAlgnEnd = alignmentResults[4]
				barAlgnStart = alignmentResults[3]
			
				if seqAlgnEnd>=seqFwdLen: #=> this means that only substring of the sequence (at the end) aligned to the front of the barcode, and we need to modify its actualy start positon 
					diff = seqAlgnEnd-seqFwdLen+1
					#modify the alignmetn end to match the end of the sequence 
					seqAlgnEnd -= diff 
					barAlgnEnd -= diff 
				
				if barAlgnEnd>=self.barcode_min_scores[barcode_ind][1]:
					diff = barAlgnEnd-self.barcode_min_scores[barcode_ind][1]+1
					seqAlgnEnd-=diff 
					barAlgnEnd-=diff
				
				algnLen = seqAlgnEnd-seqAlgnStart+1
									
			
				#confirm the number of mismatches in the alignment (this is because c-t and a-g are not penalized so there still may be mismatches)
				seqSubString = seqFwd[seqAlgnStart:seqAlgnEnd+1]
			
				#barcodeSubStringNPArray = self.barcode_char_arrays[barcode_ind][barAlgnStart:barAlgnEnd+1]				
				barcode_substring = my_barcode[barAlgnStart:barAlgnStart+algnLen]
									
				num_mismatches = compare_sequences(seqSubString,barcode_substring)
			
				#num_mismatches = 0
						
				if self.penalize_truncations == True:					
					penalties=num_mismatches+(self.barcode_min_scores[barcode_ind][1]-algnLen) #account for the number of mimsatches in the truncation region									
					max_score = self.barcode_min_scores[barcode_ind][1]-penalties
					num_mismatches = penalties				
				else:
					penalties = num_mismatches
					max_score = algnLen-penalties
				
				if penalties<=self.allowed_mismatches_in_alignment:				
					#ok alignment passed filter settings so, store results 
					pos_results[pos] = {
						'Mismatches':num_mismatches,
						'MaxScore':max_score,
						'Isotype':self.barcodeheader[barcode_ind],
						'AlgnPos':seqAlgnStart,
						'PercentSimilarity':round(100*(max_score/float(self.barcode_min_scores[barcode_ind][1])),1),
						'AlgnLen':algnLen
					}			
			#BOTH strands reported an answer
			if pos_results[0] and pos_results[1]:
				#choose the maximum value 
				if pos_results[0]['MaxScore']>=pos_results[1]['MaxScore']:
					pos_results[0]['Direction'] = 0
					results.append(pos_results[0])
				else:
					pos_results[1]['Direction'] = 1
					results.append(pos_results[1])
			#only one strand reported
			elif pos_results[0]:
				pos_results[0]['Direction'] = 0 
				results.append(pos_results[0])
			#second strand repoted 
			elif pos_results[1]:
				pos_results[1]['Direction'] = 1
				results.append(pos_results[1])			
			results.append	
		#sort results by the max score (not the percent similarity)
		results = sorted(results, key = lambda x:x['MaxScore'],reverse=True)
		return results
	
	#search for an exact match of the barcode/isotype in the sequence
	def ExactMatch(self,seq):				
		
		#perform regular expression to search for sequence 
		found_seqs = []
		#check both the forward and reverse complement sequences
		for pos in range(2):
			for barcodeInd,current_seq in enumerate(self.barcodeseq[pos]):	
				if not self.compiledRegExp[pos][barcodeInd]:
					continue
				result = self.compiledRegExp[pos][barcodeInd].search(seq)# current_seq.sub(seq)# re.search(current_seq,seq)	
				if result:		
					found_seqs.append({
						'MaxScore':len(self.barcodeseq[pos][barcodeInd]),
						'Isotype':self.barcodeheader[barcodeInd],
						'AlgnPos':result.start(),		
						'AlgnLen':len(self.barcodeseq[pos][barcodeInd]),
						'Direction':pos
					})
					
		#sort results by maxscore
		found_seqs = sorted(found_seqs, key = lambda x:x['MaxScore'],reverse=True)		
		
		return found_seqs
						
			
		
		
