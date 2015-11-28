import math
import numpy as np
import scipy
from scipy import stats
import re
import json
from pprint import pprint
import timeit
import copy

import json
from collections import defaultdict
# from immunogrep_taxonomy_tools import convert2TaxID
# from immunogrep_query_germline_functions import connectToIgDatabase
# from immunogrep_query_germline_functions import CDR3MotifDB
import time
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from numpy import prod
import timeit

class RegularExpression():
	'''
		Regular Expression Class.
		Creates PWM tables, Find's the entropy and then creates a
		regular expression.

		Contains several instance variables

		==========	   ======================
		Variable		 Definition
		----------	   ----------------------
		lPWM, rPWM	   PWM for both motifs. Numpy Array
		lEnt, rEnt	   Entropies for both motifs. Numpy Array
		lRegEx, rRegEx   Regular Expression for motifs
		RegEx			Combined left and right
		ltrim			distance from START of LEFT motif to cdr3
		rtrim			distance from END of CDR3 to START of RIGHT motif
		==========	   ================

		.. note::
			All of these can be called from simply doing class.Variable where
			class refers to whatever you decided to call the class.

			e.g.::
				>example = RegularExpression(List)
				>print(example.ltrim)
				0

		.. important::
		   The left trim (ltrim) is NOT the distance between the motif and CDR3.
		   The right trim (rtrim) is the distance between the motif and CDR3

		.. image:: CDR3Trimdisambiguation.png
			:align: center
			:alt: Reference for the right and left trims
	'''

	def __init__(self, data=[], filename=None):
		'''
			data is a list. each list must be in the
			following format

			[NAME, lmotif, rmotif, ltrim, rtrim]

			========	   ================
			Variable	   Definition
			--------	   ----------------
			Name		   'Locus'
			lmotif		 Left motif in json format
			rmotif		 Right motif in json format
			ltrim		  Distance from lmotif from CDR3, default is 0
			rtrim		  Distance from rmotif from CDR3, default is 0
			========	   ================

			example::
				RegularExpression(['TRB', lmotif, rmotif])
				or
				RegularExpression(['TRC', lmotif, rmotif, 1, 1])

			.. warning::
				If a the tuple is not 3 or 5 in length, or a string is not the
				first element, then an error is raised

		'''
		#####
		# The only reason this is here is for testing. Technically allows the
		# passing of a single file containing both motifs
		if filename is not None:
			with open(filename) as data_file:
				temp = json.load(data_file)
				data = temp[0]

		# Try catch
		try:

			isinstance(data[0], str)
			if len(data) == 3 or len(data) == 5:
				j = 1
			else:
				j = 0
			h = len(data)/j
		except ValueError:
			print("First value is not a string.")
		except ZeroDivisionError:
			print("input must contain 3 or 5 values")
		else:
			print("input is a-okay, maybe")

		# Class instance variables
		# AA list
		self.aalist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
					   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		self.Table = {'A': 'GC.', 'B': 'AA[CT]|GA[CT]', 'C': 'TG[CT]',
					  'D': 'GA[CT]', 'E': 'GA[AG]', 'F': 'TT[CT]', 'G': 'GG.',
					  'H': 'CA[CT]', 'I': 'AT[ATC]', 'K': 'AA[AG]', 'L': 'CT.',
					  'M': 'ATG', 'N': 'AA[CT]', 'P': 'CC.', 'Q': 'CA[AG]',
					  'R': 'AG[AG]|CG.', 'S': 'AG[CT]|TC.', 'T': 'AC.',
					  'V': 'GT.', 'W': 'TGG', 'X': '...', 'Y': 'TA[CT]',
					  'Z': 'CA[AG]|GA[AG]', '.': '...'}

		# for easy access from methods
		self.data = data
		self.Name = data[0]
		self.lmotif = data[1]
		self.rmotif = data[2]
		# trims are given or not
		if h == 5:
			self.ltrim = data[3]
			self.rtrim = data[4]
		else:
			self.ltrim = 0
			self.rtrim = 0

		# create your PWM
		self.lPWM = self.PWM(self.lmotif)
		# print('pwm')
		# print(self.lPWM)
		self.rPWM = self.PWM(self.rmotif)

		# create Entropies
		self.lEnt = self.entropy(self.lPWM)
		# print('entropy')
		# print(self.lPWM)
		self.rEnt = self.entropy(self.rPWM)

		# create regexpressions
		# self.lRegEx = self.RegExCreator(self.lPWM, self.lEnt)
		# self.rRegEx = self.RegExCreator(self.rPWM, self.rEnt)
		self.RegEx = re.compile(self.RegExCreator(self.lPWM, self.lEnt)+'\S+?'+self.RegExCreator(self.rPWM, self.rEnt))
		# print('regex')
		# print(self.lPWM)

	def PWM(self, motif):
		'''
			Builds the PWM from the json format.
		'''
		# print(max([int(k) for k in motif.keys()])
		PWM = [[] for x in range(0, max([int(k) for k in motif.keys()]))]
		# PWM = np.zeros(max([int(k) for k in motif.keys()]), len(self.aalist))

		for x in motif:
			for y in self.aalist:
				PWM[int(x)-1].append(motif[x][y])

		# print(PWM)
		npPWM = np.asarray(PWM)
		# print('-------------------------------------')
		# print(PWM)
		# print(len(PWM))
		return npPWM

	def entropy(self, PWM):
		'''
			Taking the PWM, finds the entropy at all the different spots
			We find entropy as follows

			2^(-sum(p * log2(p)))

		'''
		Entropies = [[] for x in range(0, len(PWM))]
		l = len(PWM)
		for i in range(0, len(PWM)):
			if PWM[i] == []:
				Entropies[i] = -1
			else:
				Entropies[i] = math.pow(2, scipy.stats.entropy(PWM[i], base=2))

		npEntropies = np.asarray(Entropies)
		# print(npEntropies)
		return npEntropies

	def RegExCreator(self, PWM, Entropies):
		'''
			Creates an Regular Expression string from the PWM
			and the entropies.

			.. note::
				This returns the entire expression including the r. So to use
				this expression, simply run re.search(regex, 'string')
				where 'string' refers to the string you're looking in.
				(not limited to search)
		'''
		TempEntropy = Entropies.tolist()
		TempPWM = copy.deepcopy(PWM.tolist())
		for x in range(0, len(TempEntropy)):
			TempEntropy[x] = int(round(TempEntropy[x], 0))
			if TempEntropy[x] == 0:
				TempEntropy[x] = 1

		# print(Entropies)

		RegExp = [[] for i in range(0, len(Entropies))]

		# finds the maximum value for a spot, finds the corresponding aa and
		# then places that AB into the RegExp list
		# -1 corresponds to no '.'
		for n in range(0, len(RegExp)):
			if TempEntropy[n] == -1:
				RegExp[n] = '.'
			else:
				for m in range(0, TempEntropy[n]):
					if max(TempPWM[n]) > .01:
						RegExp[n] = RegExp[n]+[self.aalist[TempPWM[n].index(max(TempPWM[n]))]]
						# print(max(TempPWM[n]))
						x = max(TempPWM[n])
						TempPWM[n][TempPWM[n].index(max(TempPWM[n]))] = 0
						# print(max(TempPWM[n]))
						# print('------')
						# print(x)
						if (len(RegExp) - n) == 1:
							while max(TempPWM[n]) == x and x > .01:
								# print('in while loop')
								# print(str(max(TempPWM[n])) + ' ' + str(x))
								# if str(RegExp[n]).find(self.aalist[TempPWM[n].index(max(TempPWM[n]))]) == 1:
								# print(RegExp[n])
								RegExp[n] = RegExp[n]+[self.aalist[TempPWM[n].index(max(TempPWM[n]))]]
								# print(RegExp[n])
								TempPWM[n][TempPWM[n].index(max(TempPWM[n]))] = 0
								# if x == 0 or x >= .02:
								#	print('break out of while')

		# print(RegExp)
		RegExpression = ''
		# because each index may contain at least 1 amino acid, this takes the
		# joined list and adds them to the Regular expression
		# The brackets are here for regex notation.
		for x in RegExp:
			RegExpression = RegExpression+'['+''.join(x)+']'

		# since this inserts [.], which looks for a . we need to replace with
		# '.'
		RegExpression = RegExpression.replace('[.]', '.')

		# pprint(RegExpression)

		regex = r"{}".format(RegExpression)
		# print('PWM')
		# print(PWM)
		# print('temp')
		# print(TempPWM)
		# return regex
		print(regex)
		return self.converttoNT(regex)

	def findP(self, substring, side):
		if side.upper() == 'L':
			PWM = self.lPWM
		elif side.upper() == 'R':
			PWM = self.rPWM
		else:
			print('Please provide L or R for left or right side')
			return
		P = 1
		# print(substring)
		# print(len(substring))
		# print(PWM)
		# print(len(PWM))
		for x, l in enumerate(substring):
			# print(PWM[x])
			# print('x = ' + str(x))
			# print('P = ' + str(PWM[x][self.aalist.index(l)]))
			if PWM[x] != []:
				P = P * PWM[x][self.aalist.index(l)]

		return P
	def converttoNT(self, reg):
		regnt = ''
		for c in reg:
			# print(c)
			if c == '[':
				regnt = regnt + '('
			elif c == ']':
				regnt = regnt + ')'
			elif c == '.':
				regnt = regnt + str(self.Table[c])
			else:
				regnt = regnt + str(self.Table[c]) + '|'

		regnt = str.replace(regnt, '|)', ')')
		# regnt.replace('|[', )
		# print(regnt)
		return regnt



class findCDR3():
	def __init__(self, Motif, motif_type='AA', frame=[0, 1, 2, 3, 4, 5]):
		self.Motif = Motif
		self.MotifList = [RegularExpression(data=motif) for motif in self.Motif]
		# for m in self.MotifList:
		#	  print(m.RegEx)
		self.motif_type = motif_type
		self.frame = frame

	def FindCDR3(self, seq, suggest_chain=None,
				 start_pos=0, end_pos = None, strand=None,):
		'''
		look for regex
		for each frame,
			#search all regeular expressions at once
		not found,
		find FindCDR3PWM()
		'''
		[found, results] = self.regexpsearch(seq, start_pos, end_pos, suggest_chain, strand )

		# found = False
		if found is False:
				# print('Running PWM')
				results = self.FindCDR3PWM(seq, start_pos, end_pos, suggest_chain, strand )

		return results

	def FindCDR3PWM(self, seq, start_pos=0, end_pos=None, suggest_chain=None, strand=None,):
		'''
			FindCDR3() function will use a given DNA sequence and
			list of possible LMotif and RMotifs to delimit the CDR3 regions

			.. Note:: the format of the motif is changed because the JQuery is
			reporting the motif selection in its own format

			.. Note:: the spec is that the input seq will be DNA !!!!
		'''
		# start = timeit.default_timer()
		######################################################
		# initialization
		# These are the base values.
		# MaxP is used to find the best P. Rest are pretty self-explanatory
		allscores = []
		MaxP = 1e-99*1e-99
		bestmotifset = 0
		bestscoreset = ()
		bestchain = 'NA'
		cdr3nt = ''
		cdr3aa = ''
		cdr3start = 0
		cdr3end = 0
		motifset = -1

		# this is more or less to insure that the suggest chain is a list
		# if it is not None
		if suggest_chain is not None and not type(suggest_chain) is list:
			suggest_chain = [suggest_chain]

		######################################################
		# Maps motif to variables
		# Chain = 'TRB'
		# Runs this for every motif
		for m in self.MotifList:
			# map a motif set to variables
			motifset += 1
			chain = m.Name
			LMotif = m.lmotif
			RMotif = m.rmotif
			RMotifLength = int(sorted([int(k) for k in RMotif])[-1])
			Ltrim = m.ltrim
			Rtrim = m.rtrim

			if strand is None:
				frame = [0, 1, 2, 3, 4, 5]
			elif strand == '+':
				frame = [0, 1, 2]
			elif strand == '-':
				frame = [3, 4, 5]

			# check for appropriate motif set, if suggest_chain provided then
			# ignore those motif sets that aren't matching
			if (suggest_chain is not None and suggest_chain != []
			   and chain not in suggest_chain):
				scoreset = ('N/A', '-bypassed', 'by', 'user', '!')
				allscores.append(scoreset)
				continue

			[bestLstart, bestLend, Lframe,
			 LMaxP] = PWM(seq, LMotif, frame, start_pos, self.motif_type)
			# scanning for best right flank
			[bestRstart, bestRend, Rframe,
			 RMaxP] = PWM(seq, RMotif, frame, start_pos, self.motif_type)

			# storing scoreset
			scoreset = (LMaxP*RMaxP, bestLstart, bestRend, Lframe, Rframe,
						chain)
			allscores.append(scoreset)

			# checking for best possible CDR3 regions for reporting
			if LMaxP*RMaxP > 1e-16 and LMaxP*RMaxP > MaxP and Lframe == Rframe:
				MaxP = LMaxP*RMaxP
				bestmotifset = motifset
				bestscoreset = scoreset
				bestchain = chain
				if Lframe >= 3:
					seq = str(Seq(seq, generic_dna).reverse_complement())
				# Perform trimming here to report the CDR3
				cdr3ntseq = seq[scoreset[1]:scoreset[2]]
				cdr3start = scoreset[1]+(Ltrim)*3
				cdr3nt = cdr3ntseq[(Ltrim)*3:-1*(RMotifLength-Rtrim+1)*3]
				cdr3end = cdr3start+len(cdr3nt)-1
				cdr3aa = translate_seq(str(Seq(cdr3nt, generic_dna)))

		# stop = timeit.default_timer()
		# print(stop - start)
		# return [MaxP, bestchain, cdr3aa]
		return (bestmotifset, MaxP, cdr3start, cdr3end, cdr3nt, cdr3aa,
				bestchain, bestscoreset, allscores)

	def regexpsearch(self, seq, start_pos=0, end_pos=None,
					 suggest_chain=None, strand=None):
		'''
			Returns:: bestLstart, bestRend, frame, found

			finds the reg exp for each frame. If more than one frame returns a
			regular expression, then this returns (0, 0, 0, False)

			Only returns 1 frame as it looks in the same frame for the L and R
			Motifs.
		'''

		'''
			1) Version 1 :
				Loop through self.motiflist,
				for each motif perform regular expression on each of the frames listed in frame function
			2) The minute you have more than one hit, return False and exit function
			3) After performing all motif searches, if you have a single hit, then return
				True and a set of results
		'''

		if self.motif_type != 'AA':
			print('The regex creator currently only works for AA motif types')
			return False, 0
		else:
			#####################################################
			# initialization
			# These are the base values.
			# MaxP is used to find the best P. Rest are pretty self-explanatory
			allscores = []
			MaxP = 1e-99*1e-99
			bestmotifset = 0
			bestscoreset = ()
			bestchain = 'NA'
			cdr3nt = ''
			cdr3aa = ''
			cdr3start = 0
			cdr3end = 0
			motifset = -1
			bestRegEx = None
			# this is more or less to insure that the suggest chain is a list
			# if it is not None
			if suggest_chain is not None and not type(suggest_chain) is list:
				suggest_chain = [suggest_chain]
			
			######################################################
			# Maps motif to variables
			# Chain = 'TRB'
			# Runs this for every motif
			bestLstart = -1
			bestRend = -1
			bestLtrim = -1
			bestRtrim = -1
			
			for m in self.MotifList:
				# map a motif set to variables
				motifset += 1
				chain = m.Name
				# LMotif = m.lmotif
				
				# LMotifLength = int(sorted([int(k) for k in LMotif])[-1])
				RMotif = m.rmotif
				RMotifLength = int(sorted([int(k) for k in RMotif])[-1])
				Ltrim = m.ltrim
				Rtrim = m.rtrim
				reg = m.RegEx

				if strand is None:
					frame = [0, 1, 2, 3, 4, 5]
				elif strand == '+':
					frame = [0, 1, 2]
				elif strand == '-':
					frame = [3, 4, 5]

				# check for appropriate motif set, if suggest_chain provided then
				# ignore those motif sets that aren't matching
				if (suggest_chain is not None and suggest_chain != []
				   and chain not in suggest_chain):
					scoreset = ('N/A', '-bypassed', 'by', 'user', '!')
					allscores.append(scoreset)
					continue
				if end_pos is None:
					end_pos = len(seq)
				else:
					end_pos = int(end_pos)
				seq = seq[start_pos:end_pos+1]

				# Since we're comparing against the NT sequence now, I've made
				# this run based on strand. If there is no provided strand,
				# run both, otherwise only run the specified part.
				bestframe = 0  # initial best frame
				regexpTEMP = None
				if strand is None or strand == '+':
					seq1 = str(Seq(seq, generic_dna))
					Regexp = list(reg.finditer(seq1))
					# print(Regexp)
					if Regexp != []:
						if regexpTEMP is not None or len(Regexp) > 1:
							# print('More than one regexp match. Use PWM instead')
							return (False, 0)
						elif bestLstart == -1 and bestRend == -1:
							bestLstart = Regexp[0].start()
							bestRend = Regexp[0].end()
							bestchain = chain
							bestframe = bestLstart % 3
							bestRegEx = Regexp[0]
							bestLtrim = Ltrim
							bestRtrim = Rtrim
						elif bestLstart != -1 and bestRend != -1:
							if bestRegEx is Regexp[0]:
								bestchain = [bestchain, chain]
							else:
								# print('Multiple motif matches. Use PWM instead')
								return False, 0
						else:
							# print('Multiple motif matches. Use PWM instead')
							return False, 0
				if strand is None or strand == '-':
					# reverse complement
					rseq = str(Seq(seq, generic_dna).reverse_complement())
					# Regexp = list(re.finditer(lreg+'\S+?'+rreg, aaseq[f]))
					Regexp = list(reg.finditer(rseq))
					if Regexp != []:
						if regexpTEMP is not None or len(Regexp) > 1:
							# print('More than one regexp match. Use PWM instead')
							return (False, 0)
						elif bestLstart == -1 and bestRend == -1:
							bestLstart = Regexp[0].start()
							bestRend = Regexp[0].end()
							bestchain = chain
							bestframe = 3 + bestLstart % 3
							bestRegEx = Regexp[0]
							bestLtrim = Ltrim
							bestRtrim = Rtrim
						elif bestLstart != -1 and bestRend != -1:
							if bestRegEx is Regexp[0]:
								bestchain = [bestchain, chain]
							else:
								# print('Multiple motif matches. Use PWM instead')
								return False, 0
						else:
							# print('Multiple motif matches. Use PWM instead')
							return False, 0
						# may not be the best way to handle this.
						# What this does is if there is a regexp then we give it to
						# these variables. If another regexp is found and these values
						# are not none, then we must use the PWM and find the MaxP
						if Regexp != []:
							regexpTEMP = Regexp[0]
							# print(bestchain)
							# print(Regexp[0])


			# Regex was found, so lets just make all the max = 1
			###################################################
			####### IMPORTANT!!!!
			# Kam's stuff uses everything in terms of NT. But because the
			# PWM tables were given in AA tables, the regexp ONLY applies
			# to AA this means that the following code is a bit repetitive,
			# but it follows the format Kam has already written.

			# now that this code is specific to NT, part of this is necessary.
			if bestframe >= 3:
					seq = str(Seq(seq, generic_dna).reverse_complement())
			MaxP = 999999
			scoreset = (MaxP, bestLstart, bestRend, bestframe, bestframe, chain)
			bestmotifset = motifset
			bestscoreset = scoreset
			#  bestchain = chain
			cdr3ntseq = seq[(scoreset[1]):(scoreset[2])]
			cdr3start = (scoreset[1] + bestLtrim*3)
			cdr3nt = cdr3ntseq[(bestLtrim)*3:-1*(RMotifLength - bestRtrim+1)*3]
			cdr3end = cdr3start + len(cdr3nt) - 1
			cdr3aa = translate_seq(str(Seq(cdr3nt, generic_dna)))

			if bestLstart == -1 and bestRend == -1:
				# print('No regex matches')
				# print('seq')
				# print(seq)
				# print('rseq')
				# print(str(Seq(seq, generic_dna).reverse_complement()))
				# print(m.RegEx)
				return (False, 0)
			else:
				# print('\n>>>>>>>>>>>>>>>' + str(bestLstart) + '\n')
				#return True, [MaxP, bestchain, bestframe, cdr3aa]
				return True, [bestmotifset, MaxP, cdr3start, cdr3end, cdr3nt, cdr3aa, bestchain, bestscoreset, allscores]

# PWM() function will check a sequence for the starting position and the ending
# position of the best motif
# Required inputs: PWM(seq<>text,LMotif<>dict,frame<>[0,1,2,3,4,5],
# start_pos<>intDefault0,motif_type<>AA/NTDefaultAA)
# Note: the start_pos is referencing on the original sequence input
def PWM(seq, LMotif, frame=[0, 1, 2, 3, 4, 5], start_pos=0, motif_type='AA', end_pos=None):
	# if end_pos is not provided, set it as the length of the sequence
	# print('running PWM func')
	if end_pos is None:
		end_pos = len(seq)
	else:
		end_pos = int(end_pos)
	# isolate only the section that needs to be checked
	seq = seq[start_pos:end_pos+1]
	# initialization
	LMaxP = 1e-99  # initial highest probability
	beststart = 0  # initial start of motif position
	bestend = 0  # initial end of motif position
	bestframe = 0  # initial best frame
	lwlen = int(sorted([int(k) for k in LMotif])[-1])
	lpos = [int(k) for k in LMotif]
	aaseq = defaultdict(str)
	# Checking sequence and motif compatibility and proceed with the pwm
	if validateSeq(seq, chartype='AA') and motif_type == 'NT':
		# print "!!!! Critical Error: Cannot compare Amino Acids sequence to a Nucleotide motif !!!!"
		return
	elif validateSeq(seq, chartype='NT') and motif_type == 'AA':
		# Going through all supplied frames stored in the frame list
		for f in frame:
			# Original strand
			if not f > 2:
				aaseq[f] = translate_seq(str(Seq(seq[f:], generic_dna)))
				for pos in range(len(aaseq[f])):
					lw = aaseq[f][pos:pos + lwlen]
					if (len(lw) < lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob = [LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP = prod(lwprob)
					if LP > LMaxP:
						LMaxP = LP
						beststart = start_pos+f+(pos)*3
						bestend = start_pos+f+pos*3+lwlen*3
						bestframe = f
			# Reverse complement strand
			else:
				rseq = str(Seq(seq, generic_dna).reverse_complement())
				aaseq[f] = translate_seq(str(Seq(rseq[(f-3):], generic_dna)))
				for pos in range(len(aaseq[f])):
					lw = aaseq[f][pos:pos+lwlen]
					if (len(lw) < lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob = [LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP = prod(lwprob)
					if LP > LMaxP:
						LMaxP = LP
						beststart = start_pos+f-3+(pos)*3
						bestend = start_pos+f-3+pos*3+lwlen*3
						bestframe = f

		# ---- This section is now obsolete ---- #
		# tempbeststart=defaultdict(int)
		# Note: the AA position is trying its best to give the original frame
		# AA position
		# tempbeststart['AA']=beststart
		# tempbeststart['NT']=int(beststart*3)
		# tempbestend=defaultdict(int)
		# tempbestend['AA']=bestend
		# tempbestend['NT']=int(bestend*3)
		# beststart=tempbeststart
		# bestend=tempbestend
		# ---- The above section is now obsolete ---- #

		# print(str(beststart) + ' ' + str(bestframe) + ' ' + str(LMaxP))
		return beststart, bestend, bestframe, LMaxP
	else:
		# Going through all supplied frames stored in the frame list
		for f in frame:
			if not f > 2:
				for pos in range(len(seq[f:])):
					lw = seq[f:][pos:pos+lwlen]
					if (len(lw) < lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob = [LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP = prod(lwprob)
					if LP > LMaxP:
						LMaxP = LP
						beststart = start_pos+f+pos
						bestend = start_pos+f+pos+lwlen
						bestframe = f
			else:
				rseq = str(Seq(seq, generic_dna).reverse_complement())
				for pos in range(len(rseq[(f-3):])):
					lw = rseq[(f-3):][pos:pos+lwlen]
					if (len(lw) < lwlen) or 'X' in lw or '*' in lw:
						continue
					lwprob = [LMotif[str(c+1)][lw[c]] for c in range(lwlen) if c+1 in lpos]
					LP = prod(lwprob)
					if LP > LMaxP:
						LMaxP = LP
						beststart = start_pos+f-3+pos
						bestend = start_pos+f-3+pos+lwlen
						bestframe = f
		# print(str(beststart) + ' ' + str(bestframe) + ' ' + str(LMaxP))
		return beststart, bestend, bestframe, LMaxP


def translate_seq(nt):
	'''
		Ensure suquence is multiple of 3
	'''
	rem = len(nt) % 3
	addedN = 3 - rem if rem != 0 else 0
	nt += 'N' * addedN
	return str(Seq(nt).translate())


# This function validates whether a sequence is pertaining to the type specified. i.e. NT or AA
# Note: This function cannot distinguish amino acids sequence containing only ATCG because it will return True in either NT or AA case
def validateSeq(seq, chartype='NT'):
	if chartype == 'NT':
		charset = set('ATCGUN')
	elif chartype == 'AA':
		charset = set('ACDEFGHIKLMNPQRSTVWYXZ*')
	else:
		# print "I have no clue what type you are talking about. Either NT or AA is allowed."
		return
	diff = set(seq.upper())-charset
	return not diff


# Branch point for requesting inputs
# If DBboolean is True, get MotifTable from MongoDB, else get from local file location
def getMotifInput(DBboolean, species, locus, inputfilename=None):
	if DBboolean == 'True':
		filename = motifDB2JSON(convert2TaxID(species), locus)
	else:
		if inputfilename is None:
			sys.exit("No custom Motif Table file was specified")
		filename = motifTXT2JSON(inputfilename, species, locus)
	return filename, species, locus


# Getting Motif Table from Database and output MotifTables file
def motifDB2JSON(species, locus):
	lpos, rpos, lmotif, rmotif, lcommotif, rcommotif, maxlpos, maxrpos = getMotifTable(species, locus)
	pattern = lcommotif+'(.*)'+rcommotif
	lwlen = int(lpos[-1])
	rwlen = int(rpos[-1])
	timestr = time.strftime("%Y%m%d-%H%M%S")
	fout = open('scratch/'+timestr+'.motifT', "w")
	print >>fout, json.dumps(lpos)
	print >>fout, json.dumps(rpos)
	print >>fout, json.dumps(lmotif)
	print >>fout, json.dumps(rmotif)
	print >>fout, json.dumps(lcommotif)
	print >>fout, json.dumps(rcommotif)
	print >>fout, json.dumps(maxlpos)
	print >>fout, json.dumps(maxrpos)
	fout.close()
	return 'scratch/'+timestr+'.motifT'


# Converting custom Motif Table from user provided location and output MotifTables file
def motifTXT2JSON(inputfilename,species,locus):
	# Processing the user input motif table
	Llist = []
	Rlist = []
	L = 0
	R = 0
	with open(inputfilename) as f:
		for line in f:
			line = line.replace("\r", '')
			line = line.replace("\n", '')
			if line[0] == "\t":
				continue
			if 'LPosition' in line:
				L = 1
				R = 0
				tmp = line.split("\t")
				pos = tmp[1:]
				pos = [x for x in pos if x != '']
				r = str()
			if 'RPosition' in line:
				L = 0
				R = 1
				tmp = line.split("\t")
				pos = tmp[1:]
				pos = [x for x in pos if x != '']
				r = str()
			if L == 1 and 'LPosition' not in line:
				tmp = line.split("\t")
				aa = tmp[0]
				prob = tmp[1:]
				for i in range(len(pos)):
					r = aa + '_' + pos[i] + '_' + prob[i]
					Llist.append(r)
			if R == 1 and 'RPosition' not in line:
				tmp = line.split("\t")
				aa = tmp[0]
				prob = tmp[1:]
				for i in range(len(pos)):
					r = aa+'_'+pos[i]+'_'+prob[i]
					Rlist.append(r)
	f.close()
	# Generating necessary outputs to follow the standards outputs needed for MotifTables file
	Lpos = []
	Rpos = []
	Lprob = defaultdict(dict)
	Rprob = defaultdict(dict)
	for i in Llist:
		tmp = i.split('_')
		if tmp[1] not in Lpos:
			Lpos.append(tmp[1])
		Lprob[tmp[1]][tmp[0]] = float(tmp[2])
	for i in Rlist:
		tmp = i.split('_')
		if tmp[1] not in Rpos:
			Rpos.append(tmp[1])
		Rprob[tmp[1]][tmp[0]] = float(tmp[2])
	LcomMotif = str()
	for i in range(int(Lpos[-1])):
		itera = str(i+1)
		maxprob = 0.3
		cnt = 0
		if (len(Lprob[itera]) == 0):
			cnt = 100
		for k in Lprob[itera]:
			if Lprob[itera][k] >= maxprob:
				cnt += 1
				maxprob = Lprob[itera][k]
				maxaa = k
				if (maxaa == 'C' or maxaa == 'Y') and maxprob > 0.81:
					maxlpos = int(itera)
		if cnt == 1:
			LcomMotif = LcomMotif + maxaa
		else:
			LcomMotif = LcomMotif + '.{1}'
	RcomMotif = str()
	for i in range(int(Rpos[-1])):
		itera = str(i+1)
		maxprob = 0.3
		cnt = 0
		if (len(Rprob[itera]) == 0):
			cnt = 100
		for k in Rprob[itera]:
			if Rprob[itera][k] >= maxprob:
				cnt += 1
				maxprob = Rprob[itera][k]
				maxaa = k
				if (maxaa == 'W' or maxaa == 'F') and maxprob > 0.81:
					maxrpos = int(itera)
		if cnt == 1:
			RcomMotif = RcomMotif+maxaa
		else:
			RcomMotif = RcomMotif+'.{1}'
	# Output standards to Motif Tables file
	timestr = time.strftime("%Y%m%d-%H%M%S")
	fout = open('scratch/'+timestr+'.motifT', "w")
	print >>fout, json.dumps(Lpos)
	print >>fout, json.dumps(Rpos)
	print >>fout, json.dumps(Lprob)
	print >>fout, json.dumps(Rprob)
	print >>fout, json.dumps(LcomMotif)
	print >>fout, json.dumps(RcomMotif)
	print >>fout, json.dumps(maxlpos)
	print >>fout, json.dumps(maxrpos)
	fout.close()
	return 'scratch/'+timestr+'.motifT'


# This function converts MotifTable .motifT file into variables
def convertMotifTable2Var(inputjson):
	with open(inputjson) as f:
		lines = [f.next() for i in range(8)]
		lpos = [int(i) for i in json.loads(lines[0])]
		rpos = [int(i) for i in json.loads(lines[1])]
		lmotif = json.loads(lines[2])
		rmotif = json.loads(lines[3])
		lcommotif = json.loads(lines[4])
		rcommotif = json.loads(lines[5])
		pattern = lcommotif+'.*'+rcommotif
		lwlen = int(lpos[-1])
		rwlen = int(rpos[-1])
		maxlpos = json.loads(lines[6])
		maxrpos = json.loads(lines[7])
	return (lpos, rpos, lmotif, rmotif, lcommotif, rcommotif,
			pattern, lwlen, rwlen, maxlpos, maxrpos)


# Get Motif Table from MongoDB database
# This function requires Species (Unique numbers) and Locus (IGH,IGL,IGK)
def getMotifTable(spee, loc):
	db = connectToIgDatabase()
	motif = db['motifs']
	prob = motif.find_one({u'TaxID':spee, u'Locus':loc})
	Lpos = []
	Rpos = []
	Lprob = defaultdict(dict)
	Rprob = defaultdict(dict)
	# Storing the position and the probability tables into dictionaries
	try:
		prob['Lmotif']
	except:
		sys.exit("!!!! Warning: Data does not exist on MongoDB !!!!")

	for i in prob['Lmotif']:
		tmp = i.split('_')
		if tmp[1] not in Lpos:
			Lpos.append(tmp[1])
		Lprob[tmp[1]][tmp[0]] = float(tmp[2])
	for i in prob['Rmotif']:
		tmp = i.split('_')
		if tmp[1] not in Rpos:
			Rpos.append(tmp[1])
		Rprob[tmp[1]][tmp[0]] = float(tmp[2])
	# Finding the most common Left flanking and Right flanking motifs
	LcomMotif = str()
	tempLprob = Lprob.copy()
	for i in range(int(Lpos[-1])):
		itera = str(i+1)
		maxprob = 0.3
		cnt = 0
		if (len(tempLprob[itera]) == 0):
			cnt = 100
		for k in tempLprob[itera]:
			if tempLprob[itera][k] >= maxprob:
				cnt += 1
				maxprob = tempLprob[itera][k]
				maxaa = k
				if (maxaa == 'C' or maxaa == 'Y') and maxprob > 0.81:
					maxlpos = int(itera)
		if cnt == 1:
			LcomMotif = LcomMotif + maxaa
		else:
			LcomMotif = LcomMotif+'.{1}'
	RcomMotif = str()
	tempRprob = Rprob.copy()
	for i in range(int(Rpos[-1])):
		itera = str(i+1)
		maxprob = 0.3
		cnt = 0
		if (len(tempRprob[itera]) == 0):
			cnt = 100
		for k in tempRprob[itera]:
			if tempRprob[itera][k] >= maxprob:
				cnt += 1
				maxprob = tempRprob[itera][k]
				maxaa = k
				if (maxaa == 'W' or maxaa == 'F') and maxprob > 0.81:
					maxrpos = int(itera)
		if cnt == 1:
			RcomMotif = RcomMotif + maxaa
		else:
			RcomMotif = RcomMotif+'.{1}'

	return Lpos, Rpos, Lprob, Rprob, LcomMotif, RcomMotif, maxlpos, maxrpos


# This function will fetch all unique species to a list
def fetchSpeciesList():
	db = CDR3MotifDB()
	return db.GetUniqSpecies()





