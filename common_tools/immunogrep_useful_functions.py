import sys
import os
from datetime import datetime
from immunogrep_global_variables import listofextension
from immunogrep_global_variables import descriptor_symbol
import immunogrep_read_file as readwrite
import traceback
import subprocess
from collections import MutableMapping
import glob

dna_codes = {
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A',
	'N': 'N',
	'X': 'X',
	'-': '-',
	'U': 'A',
	'W': 'W',
	'M': 'K',
	'K': 'M',
	'S': 'S',
	'R': 'Y',
	'Y': 'R',
	'a': 't',
	'c': 'g',
	'g': 'c',
	't': 'a',
	'u': 'a',
	'w': 'w',
	'n': 'n',
	'm': 'k',
	'k': 'm',
	's': 's',
	'r': 'y',
	'x': 'x',
	'y': 'r'
}


def purge(filelist):
	"""
		Delete all files with same substring at start of file path
	"""
	for f in filelist:
		f = f.replace('*', '\*')
		for filename in glob.glob(f + "*"):
			os.remove(filename)	

def get_stdout(bash_command):
	"""
		Python wrapper for running python subprocess.call function
		Returns the output from command
	"""
	output = subprocess.check_output(bash_command, shell=True)
	return output


def file_line_count(filepath):
	"""
		Takes a file path as input. Returns the number of lines in the file using wc -l

		.. warning::
			If path is not found, raises an error
	"""
	if os.path.isfile(filepath):
		filepath = os.path.abspath(filepath)
		value = get_stdout("wc -l '{0}'".format(filepath)).split()[0]
		if value:
			return int(value)
		else:
			return 0
	else:
		raise Exception('File does not exist: ' + filepath)


def Reverse_Complement(x):
	"""
		Takes in sequence x and returns its complement (?)
	"""
	rc_list = [dna_codes[c] if c in dna_codes else 'N' if ord(c) < 91 else 'n' for c in reversed(x)]
	return ''.join(rc_list)


def get_parent_dir(path):
	"""
		Uses os.path to return the parent directory of a file
	"""
	return os.path.dirname(os.path.abspath(path))


def split_files_by_seq(filepath, num_files_to_make, number_lines_per_seq, contains_header_row):
	"""
		This function is used for splitting a single file into multiple little files. This will be useful for multithreading purposes.

		.. note:: First make sure files do not start with any documentation fields

		Parameters
		----------
		filepath : string
			Location of the input file
		num_files_to_make : int
			Defines the number of smaller files you want to create
		number_lines_per_seq : int
			Defines the number of lines that correspond to a single sequence
			i.e. fastq files = 4 & fasta files = 2
		contains_header_row : boolean
			Defines whether or not the input file contains a header row (i.e. CSV or TAB files)

		.. important::
			We cannot just use the split command, instead we have to use awk command. This is because we need to copy header lines and IGREP documentation lines to each split file.
	"""
	num_doc_fields = 0
	header_lines = ''

	with open(filepath, 'r') as i:
		while True:
			line = i.readline()
			if line.startswith(descriptor_symbol):
				num_doc_fields += 1
				header_lines += line
			else:
				if contains_header_row:
					num_doc_fields += 1
					# NOW read the top header line
					header_lines += line
				break

	parent_path = get_parent_dir(filepath)
	if num_doc_fields > 0:
		my_temp_file = os.path.join(parent_path, 'header_rows_' + str(datetime.now()).replace(' ', ''))
		with open(my_temp_file, 'w') as e:
			e.write(header_lines)

	# First get the number of lines in the file
	num_lines = file_line_count(filepath)
	# Make sure the number of lines per seq matches what would be expected from num_lines
	if(num_lines % number_lines_per_seq != 0):
		raise Exception('Number of lines in file is not divisible by the number of lines for each sequence: {0}/{1}!=0'.format(str(num_lines), str(number_lines_per_seq)))

	num_seqs = num_lines / number_lines_per_seq
	# Determine how many lines to split the file by
	# Round down, division THEN add 1
	# We add a 1 to ensure that the num_files_to_make is the max number of files made
	num_lines_in_split_files = (int(num_seqs / num_files_to_make) + 1) * number_lines_per_seq

	# Make a temp folder
	subfolder = os.path.join(parent_path, str(datetime.today()).replace(' ', '_').replace(':', '').replace('.', ''))
	os.makedirs(subfolder)

	# Ok THERE are header lines at the top fo the file that we dont watn to split so we need to ignore these lines
	if num_doc_fields > 0:
		# The final '-' is IMPORTANT when combining with tail
		system_command = "tail -n +{5} '{1}'| split -l {0} -a {4} - '{2}/{3}'".format(str(num_lines_in_split_files),filepath,subfolder,os.path.basename(filepath+'.'),num_files_to_make/10+1,num_doc_fields+1)
	else:
		system_command = "split -l {0} -a {4} '{1}' '{2}/{3}'".format(str(num_lines_in_split_files),filepath,subfolder,os.path.basename(filepath+'.'),num_files_to_make/10+1)


	#run the bash split command , split files by lines defined above, generate suffixes whose length is equal to the number of files ot make, output results to temp subfolder, prefix files with inputfile+'.'
	#os.system(system_command)
	subprocess.call(system_command,shell=True)

	#return the contents of files made
	files_created = os.listdir(subfolder)



	if num_doc_fields>0:
		#move files back up one folder, but while moving files also concatenate teh header file
		hack_bash_command = '''
			for file in "{0}/"*
			do
				s=$(basename "$file")
				cat "{1}" "$file" > "{2}/$s"
			done
			rm "{1}"
			rm -r "{0}"
		'''.format(subfolder,my_temp_file,parent_path)
		#os.system(hack_bash_command)
		subprocess.call(hack_bash_command,shell=True)
	else:
		#move files back to starting folder , delete temp folder
		#os.system("mv {0}/*.* {1}/.;rm -r {0}".format(subfolder.replace(' ','\ '),parent_path.replace(' ','\ ')))
		subprocess.call("mv {0}/*.* {1}/.;rm -r {0}".format(subfolder.replace(' ','\ '),parent_path.replace(' ','\ ')),shell=True)

	return sorted([parent_path+'/'+f for f in files_created])



#runs an awk script for merging results from multiple files
def merge_multiple_files(file_list, num_header_lines=1,outfile=None):
	"""
        Runs an awk script for merging results from multiple files

        	=====================   =====================
	        **Input**               **Description**
	        ---------------------   ---------------------
	        file_list               A list of strings corresponding to the path of all files to merge
	        num_header_lines	    The number of lines (x) in each file that correspond to a header row. When, FNR!=NR the first x lines are ignored. In other words only the first x lines from the first file are copied to the new file.
	        outfile			        The location of the output file (optional). When empty it defaults to the first input file +'.merged'
	        =====================   =====================
    """
	if not outfile: #no outfile, then make default path
		outfile = file_list[0]+'.merged'

	file_list = [f.replace(' ','\ ') for f in file_list]
	file_list_strings = ' '.join(file_list)#["'"+f+"'" for f in file_list])

	awk ='''
		awk 'FNR!=NR&&FNR<={0}{{next}};
		{{print $0> "{1}" }}' {2}'''.format(str(num_header_lines),outfile,file_list_strings)

	#os.system(awk)
	output = subprocess.check_output(awk,shell=True)# call(awk)
	return outfile

try:
	from bson.objectid import ObjectId
	oid_type = type(ObjectId())
except:
	print('Cannot convert objectid to str. install pymongo')
	oid_type = str

def flatten_dictionary(d,val={},p='',start=True):
	"""
		This function recursively flattens dictionaries.
		All values in the dictionary that are found to be nested subdictionaries will be converted to keys containing '.'
		This is very useful for moving all keys in a dictionary to the top level rather than having to access each subkey

			=====================   =====================
			**Input**				**Description**
			---------------------	---------------------
			d						A nested dictionary of key:values
			=====================   =====================

		See example below::

			Imagine a dictionary is passed as such
				d = {
					'a':{'hello':5,
						'world':10}
					'ok':10
				}
			Running flatten_dictionary on this function will produce the following
			d = {
				'a.hello':5,
				'a.world':10,
				'ok':10,
			}
			So rather than accesssing d['a']['hello'], We access d['a.hello']

		.. important:: Please refer to the function,DotAccessible, to perform the opposite of this function: Convert a flattened dictionary into a nested dictionary

    """
	if start:
		val = {}
	for k,v in d.iteritems():
		if isinstance(v, dict):
			flatten_dictionary(v,val,p + k + '.', False)
		elif isinstance(v,oid_type):
			val[p+k]=str(v)
		else:
			val[p+k] = v
	return val



class DotAccessible(MutableMapping):
	"""
		A dictionary that allows dot notation for accessing nested subdictionaries. Code [inspired-by/copied-from]:
		http://www.velvetcache.org/2012/03/13/addressing-nested-dictionaries-in-python
		http://stackoverflow.com/questions/3387691/python-how-to-perfectly-override-a-dict
		http://docs.python.org/2/library/collections.html#collections-abstract-base-classes

	"""

	def __init__(self, *args, **kwargs):
		self.store = dict()
		self.update(dict(*args, **kwargs))  # use the free update to set keys

	def __getitem__(self, dotkey):
		value = self.store
		for key in dotkey.split('.'):
			value = value[key]
		return value

	def __setitem__(self, dotkey, value):
		nested_node = self.store
		for key in dotkey.split('.')[:-1]:
			if key not in nested_node:
				nested_node[key] = dict()
			nested_node = nested_node[key]
		nested_node[dotkey.split('.')[-1]] = value

	def __delitem__(self, dotkey):
		nested_node = self.store
		for key in dotkey.split('.')[:-1]:
			nested_node = nested_node[key]
		del nested_node[dotkey.split('.')[-1]]

	def __iter__(self):
		return iter(self.store)

	def __len__(self):
		return len(self.store)

	def __str__(self):
		return str(dict(self))


def RemoveObjId(document):
	"""
        This function will recursively search all the values of a document/dictionary. Any value in the dictionary whose type is equal to a MongoDB ObjectId will be converted to a string.

    		=====================   =====================
	        **Input**               **Description**
	        ---------------------   ---------------------
	        document		A nested or flattened dictionary; or similarly a document returned by MongoDB
	        =====================   =====================
    """
	from bson.objectid import ObjectId
	oid_type=type(ObjectId())
	def remove_oid(document):
		for f,v in document.iteritems():
			if isinstance(v,dict):
				remove_oid(v)
			#its an object id
			elif isinstance(v,oid_type):
				document[f]=str(v)
			elif isinstance(v,list) and isinstance(v[0],dict):
				for sub_vals in v:
					if isinstance(sub_vals,dict):
						remove_oid(sub_vals)
	remove_oid(document)

#function description - this will extract a specific field from a file and write it to the output file as a  single column (without a header)
#if count_field = None, then we assume there are no counts associated with the field of interest.  if count_field =is not None then we assume that column refers to the number of counts a sequence has occurred
def Write_Single_Field(filename=None,outfile_location=None,field=None,count_field = None, file_format=None,contains_header=True):
	"""
			This will extract a specific field from a file and write it to the output file as a single column (without a header.)
			If *count_field* = None, then we assume there are no counts associated with the field of interest.
			If *count_field* =! None, then we assume that column refers to the number of counts a sequence has occurred
	"""

	total_data = 0
	total_field = 0

	if outfile_location == None:
		outfile_location = filename+'_singlefield.txt'

	if (filename):
		isfile = os.path.isfile(filename)
	if (filename==None) or (isfile==False):
		raise Exception("The pathname of the file is invalid")
	if (field==None):
		IF_file = readwrite.immunogrepFile(filelocation=filename,filetype='TAB',contains_header=False,mode='r')
		print("Warning no field name was provided.  This file will be treated as a tab file and the first column will be selected")
		field = 'Column 1'
	else:
		IF_file = readwrite.immunogrepFile(filelocation=filename,filetype=file_format,contains_header=contains_header,mode='r')

		if (file_format==None):
			guessedFiletype = IF_file.getFiletype()
			print("Warning, no file type for this file was provided.  The file was predicted to be a "+guessedFiletype+" file.")
	try:
		outfile = open(outfile_location,'w')
		while not(IF_file.IFclass.eof):
			data = IF_file.IFclass.read() #read in each line as a dictionary
			if data:
				total_data+=1
				if field in data and data[field]:
					value = data[field]
					if count_field!=None and count_field in data and data[count_field]:
						count = data[count_field] #this defines the number of times we will write it to a file
					else:
						count = '1'
					outfile.write(value+'\t'+count+'\n')#write sequence to new file
					total_field+=1
	except Exception as e:
		os.remove(outfile_location)#("rm '{0}'".format(outfile_location))
		print_error(e)

	return [total_data,total_field]


def gunzip_python(path): 
	new_path = path.replace('.gz','')
	subprocess.call('gzip -c {0} > {1}'.format(path, new_path),shell=True)
	return new_path

#Simple python command for counting the occurrences of a sequence ina file.
#assumption 1: sequences are sorted alphabetically
#assumption 2: tab delimited file
#assumption 3: column 1 = sequence of interest, column 2 = counts for that sequence
def count_sorted_seqs(input_file_name, output_file_name):
	"""
		Counts the occurances of a sequence in a file. Takes a few assumptions into account.

		:Assumptions:
			* Sequences are sorted alphabetically
			* tab delimited file
			* column 1 = sequence of interest
			* column 2 = counts for that sequence
	"""
	f_in = open(input_file_name,'r')
	f_out = open(output_file_name,'w')

	line_one = f_in.readline().strip()
	line_one = line_one.split('\t')
	current_seq = line_one[0]
	if len(line_one)>1:
		current_count = int(line_one[1])
	else:
		current_count = 1

	for line in f_in:
		data = line.strip()
		data = data.split('\t')
		temp_seq = data[0]
		if len(data)>1:
			temp_count = int(data[1])
		else:
			temp_count = 1

		if (temp_seq == current_seq): #same sequence as before, so keep adding to teh counts
			current_count+=temp_count
		else:
			f_out.write(current_seq+'\t'+str(current_count)+'\n') # new sequence encountered so output results for old sequence and replace new sequence
			current_seq = data[0]
			current_count = temp_count
	f_out.write(current_seq+'\t'+str(current_count)+'\n')
	f_in.close()
	f_out.close()


#filelocation-> location of file
#field -> name of the column/field in the file that you want to read
#file_format -> JSON/TAB/CSV/ETC...
#contains_header -> whether or not the file contains a header row
def count_unique_values(filelocation=None,output_filelocation=None,field=None,count_field=None,file_format=None,contains_header=True,delete_intermediate_file=True,mem_safe = False):
	"""
	        ===============   ================
	        **Input**         **Description**
	        ---------------   ----------------
	        filelocation      location of file
	        field             name of the column/field in the file that you want to read
	        file_format       JSON/TAB/CSV/etc.
	        contains_header   whether or not the file contains a header row
	        ===============   ================
	"""
	if output_filelocation == None:
		output_filename = removeFileExtension(filelocation)
		output_filename = output_filename+'.unique_vals.txt'
		unique_counts = output_filename+'.counts'
	else:
		output_filename = output_filelocation+'.single_field_list'
		output_filename2 = output_filelocation+'_temp2'
		unique_counts = output_filelocation

	[total_found,total_field] = Write_Single_Field(filelocation,output_filename,field,count_field, file_format,contains_header) #take the field we are interested in and just write it to a temp file

	if mem_safe:
		#bash_command = '''sort '{0}' | uniq -c | awk 'BEGIN{{OFS="\t"}}{{ print $2,$1,$1 }}' > '{1}' '''.format(output_filename,unique_counts)	 --> oldest method, no longer used

		#USE THESE LINES IF WE FIND THAT THERE ARE MEMORY LIMITATIONS/ERRORS WITH THE FUNCTION##
		#print datetime.now()
		bash_command = '''sort '{0}' > '{1}' '''.format(output_filename,output_filename2)
		#os.system(bash_command)
		subprocess.call(bash_command,shell=True)
		count_sorted_seqs(output_filename2,unique_counts)
		os.remove(output_filename2)
		#os.system("rm '{0}'".format(output_filename2))
		#print datetime.now()
		###END OF MEMORY SMART BUT SLOER FUNCTION ################
	else:
		#FAST method, but could have memory limitations for extremely large files)
		bash_command = '''awk '{{OFS = "\t"}} {{arr[$1]+=$2}} END {{for (i in arr) {{print i,arr[i]}}}}' '{0}' > '{1}' '''.format(output_filename,unique_counts)
		#os.system(bash_command)
		subprocess.call(bash_command,shell=True)

	if delete_intermediate_file:
		os.remove(output_filename)

	return [total_found, total_field]


# print an error message if system/program fails
def print_error_string(e=None):
	"""
		Prints an error message if system/program fails.
		Uses traceback to report the last line and module where error occurred

		:exception code:

			if not(e):
				e = "Exception error not passed"
			else:
				print there was an error: str(e)
    """
	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]

	if not(e):
		error = "Exception error not passed\n"
	else:
		error = "There was an error: \n"+str(e)

	error+=	"\nLine Number: {0} , type: {1} , fname: {2}".format(str(exc_tb.tb_lineno),str(exc_type),str(fname))

	tb_error = traceback.format_exc()
	error+="Traceback Error:\n"+str(tb_error)
	return error







# this is a cleanup function to remove any unused fields
# this will iterate the sequence document and remove any subdocuments or fields that have "None" values
def removeEmptyVals(myDict):
	"""
		this is a cleanup function to remove any unused fields
		this will iterate the sequence document/dictionary and remove any subdocuments or fields that have "None" or empty string or empty array/dict values

			================	================
			**Input**			**Description**
			----------------	----------------
			myDict				A python dictionary of keys/values
			================	================

	"""
	if myDict:
		copyDict = myDict.copy();

		for myKeys in myDict:
			if type(myDict[myKeys]) is dict:
				#if it has a subdocument then recursively remove any None elements in subdocument
				subDict = removeEmptyVals(myDict[myKeys])
				if subDict == {}:
					copyDict.pop(myKeys,None)
				else:
					copyDict[myKeys] = subDict;
			else:
				if myDict[myKeys]!=0 and myDict[myKeys]!=False and not myDict[myKeys]:
					copyDict.pop(myKeys,None)
	else:
		copyDict = myDict

	return copyDict

#ASSUMES THAT THE DICTIONARY IS ONLY IN DOT NOTATION!!! I.E. NOT A NESTED DICTIONARY
# this is a cleanup function to remove any unused fields
# this will iterate the sequence document and remove any subdocuments or fields that have "None" values or empty strings or empty lists
#any fields that are 'empty' are retuend as a second variable
def divideEmptyAndNonEmptyVals(myDict):
	"""
		This is a cleanup function that removes any unused fields.
		This will iterate the sequence document and remove any subdocuments or fields that have "None" values, empty strings or empty lists.
		Fields that are empty are returned as a second variable.

	        ===============   ================
	        **Input**         **Description**
	        ---------------   ----------------
	        myDict			  An already flattened dictionary
	        ===============   ================

	        .. note::
			This method assumes that the dictionary is in dot notation. i.e. not a nested dictionary
	"""
	empty_fields={}
	non_empty = {}
	for field,value in myDict.iteritems():
		if value!=False and not(value):
			empty_fields[field]= ""
		else:
			non_empty[field] = value

	return [non_empty,empty_fields]

# this is a cleanup function to remove any unused fields
# this will iterate the sequence document and remove any subdocuments or fields that have "None" values
def removeNoneVals(myDict):
	"""
		this is a cleanup function to remove any unused fields
		this will iterate the sequence document/dictionary and remove any subdocuments or fields that have "None" only

			===============   ================
			**Input**         **Description**
			---------------   ----------------
			myDict            A python dictionary of keys/values
			===============   ================

	"""
	if myDict:
		copyDict = myDict.copy();

		for myKeys in myDict:
			if type(myDict[myKeys]) is dict:
				#if it has a subdocument then recursively remove any None elements in subdocument
				subDict = removeNoneVals(myDict[myKeys])
				if subDict == {}:
					copyDict.pop(myKeys,None)
				else:
					copyDict[myKeys] = subDict;
			else:
				if myDict[myKeys]==None:
					copyDict.pop(myKeys,None)
	else:
		copyDict = myDict

	return copyDict



#spint out a counter of the current status of a process (i.e. percent done)
def LoopStatus(counter,totalSeq,perIndicate,startPer,div='',addedInfo=None):
	"""
        	counter of the current status of a process. (i.e. percent done)
    """
	percentDone = int(counter/float(totalSeq)*100)
	if percentDone%perIndicate==0 and percentDone>startPer:
		stringvar ='{0}% percent done. Time: {1}'.format(str(percentDone),str(datetime.now()))
		if addedInfo:
			stringvar+='\n{0}\n\n'.format(addedInfo)
		print(stringvar)
		startPer = percentDone
	return startPer

#spint out a counter of the current status of a process (i.e. percent done)
#THIS FUNCTION USES GENERATOR INSTEAD
def LoopStatusGen(totalSeq,perIndicate,addedInfo=None):
	"""
		counter of the current status of a process.

		.. note::
			this function uses generator instead

		.. seealso::
			:py:func:`.LoopStatus`
	"""
	counter=0
	startPer = 0
	while True:
		percentDone = int(counter/float(totalSeq)*100)
		if percentDone%perIndicate==0 and percentDone>startPer:
			stringvar ='{0}% percent done. Time: {1}'.format(str(percentDone),str(datetime.now()))
			if addedInfo:
				stringvar+='\n{0}\n\n'.format(addedInfo)

			print(stringvar) #print out the current perecent
			startPer = percentDone
		yield startPer #use a generator to pause
		counter+=1

def removeFileExtension(stringFileName):
	"""
		Removes File Extension
    """
	filename = stringFileName.split('.')
	lastExtension = filename[-1]

	foundExtension = False

	for i in listofextension:
		if i == lastExtension:
			foundExtension = True

	if foundExtension:
		editedFileName = filename[0]
		for i in range(1,len(filename)-1):
			editedFileName+="."+filename[i]
	else:
		editedFileName = stringFileName

	return editedFileName


def fieldsForAnnotatingAb():

	abFields = {
		"FULL_SEQ":None, #dna sequence
		"SEQ_HEADER":None, #header for dna sequence
		"STRAND_DIR":None, #after the alignment is the v/d/j gene aligned to the fwd or reverse complement?
		"VGENE_START":None,
		"JGENE_START":None,
		"FR1_START_NT":None, #nucleotide position along sequence where FR1 starts
		"FR1_END_NT":None,
		"CDR1_START_NT":None,
		"CDR1_END_NT":None,
		"FR2_START_NT":None,
		"FR2_END_NT":None,
		"CDR2_START_NT":None,
		"CDR2_END_NT":None,
		"CDR3_START_NT":None,
		"CDR3_END_NT":None,
		"FR4_START_NT":None,
		"FR4_END_NT":None,
		"SEQ_ALGN":None, #the actual sequence that is aligned to the germline (includes gap characters (-)), to be used for later scripts
		"GERMLINE_ALGN":None, #the actual germline sequence that best aligns to sequence (includes gap characters (-)), to be used for later scripts
		"START_FEATURE":None #what region does the antibody start/along alignment
	};
	return abFields

