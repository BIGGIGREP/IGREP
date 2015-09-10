# tools for repsoma
import os
import time
import errno
import shutil # needed for a probably unnecessary class below
#import appsoma_api
from datetime import datetime
import subprocess
#---A GENERATOR FUNCTION FOR PREVIEWING ZIPPED FILES. THIS FUNCTION IS WITH THE 'SHOWRESULTS' CLASS IN MAKEAPPFUNCTIONS. 
#THIS IS THE PREVIEW FUNCTION USED TO GENERATE A PREVIEW OF THE ZIPPED FILE IN THE APP 
#IT IS CURRENTLY CALLED IN THE FOLLOWING SCRIPTS 1) DOWNLOADSEQUENCESAPP, 2)FASTQ_TO_FASTA
import zipfile
import tarfile
def PreviewCompressedFilesGenerator(filepath):
	if not os.path.isfile(filepath):
		yield None
	else:
		if filepath.endswith('.zip'):
			zip = zipfile.ZipFile(filepath)
			yield {'Files':zip.namelist()}
		elif filepath.endswith('.tar.gz') or filepath.endswith('.tar.bz2') or filepath.endswith('.tar.gz') or filepath.endswith('tgz') or filepath.endswith('tbz'):
			tar = tarfile.open(filepath)
			yield {'Files':tar.getnames()}
		else:
			yield None

#------PARAMETERS liable to change------#
_projName_ = 'igrep' # When we think of something better than 'repsoma'
expDirParams = {'start':'Exp', # Experiment directory format.
				'delim':'_',
				'ndigits':5} # number of digits in incrementing experiment counter
# expDirParams controls how the experiment dirs are named.
#---------------------------------------#
expDirTemplate = "{start}{delim}{{exp_num:0{nd}d}}{delim}{{exp_name}}".format(nd=expDirParams['ndigits'], **expDirParams)
""" This, and expDirParams, define how experiment folders are named.
	Later functions assume expDirTemplate structure doesn't change. Let's 
	discuss convention if anyone cares.
	Usage example:
	>>	expDirTemplate.format(exp_num = 1123, exp_name = Platyhelminthes)
	Exp_01123_Platyhelminthes
	(or something else if expDirParams has changed)
"""
#def url_to_path(path):
#	"""convert a workspace URL to a home directory path.
#	Has no effect if given path is not a workspace URL."""
#	ws_url = appsoma_api.environment_get_workspace_url()
#	home_dir = appsoma_api.environment_get_home()
#	prefix = os.path.commonprefix([ws_url,path])
#	return os.path.join(home_dir,path[len(prefix):])


def InitializeHomeScratch(home):
	"""

		Make sure a 'scratch' folder is found directly at the users' home folder.
		If not, then make the directory

	"""
	scratch = os.path.join(home,'scratch')
	if not os.path.isdir(scratch):
		os.makedirs(scratch)
	return scratch

#get the current users' home directory 
home = os.path.expanduser('~')
#within the home directory, we want to define a default scratch folder to start from
scratch = InitializeHomeScratch(home)

def mkdir_p(path):
	"""same as shell 'mkdir -p'. Will not throw error if dir already exists."""
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else: raise

def dir_files(path):
	"""Get files within top level of given directory.
	Return dict with keys:
		'_all_' : all files
		'<ext1>': files with <ext1>
		'<ext2>': files with <ext2>
	etc. """
	if path[-1]!=os.sep:
		path+=os.sep
	try:		
		files = os.walk(path,followlinks=True).next()[2]		
	except:
		raise Exception("{} is not accessible".format(path))
	contents = {'_all_' : files}
	[contents[ext[1:]].append(base+ext) if contents.has_key(ext[1:]) 
			else contents.update({ext[1:]:[base+ext]})
			for base,ext in [os.path.splitext(f) for f in files]]
	return contents
	
def list_all_dirs_in_path(path):
	found_dirs = [i[0] for i in os.walk(path,followlinks=True)]
	return found_dirs

#recursively searches all files in all folders under parent folder "path". Not just top level
#trim = True => remove everything up until the last directory or path
def get_all_files_recursively(path='scratch',trim=False):
		
	if path[-1] == os.sep:
		path = path[:-1]
	
	basename_path = os.path.basename(path)
	
	files = []
	#try:

	for i in os.walk(path,followlinks=True):
		if i[-1]:		
			folder_name = i[0]+'/'
			if trim:
				path_to_find = '/'+basename_path+'/'
				found_pos =  folder_name.find(path_to_find)
				
				if found_pos>=0:																
					folder_name = folder_name[found_pos+1:]
			
																				
			for j in i[-1]:
				files.append(folder_name+j)

	contents = {'_all_' : files}
	[contents[ext[1:]].append(base+ext) if contents.has_key(ext[1:]) 
			else contents.update({ext[1:]:[base+ext]})
			for base,ext in [os.path.splitext(f) for f in files]]
	return contents

#recursively searches all files in all folders under parent folder "path". Not just top level
#trim = True => remove everything up until the last directory or path
def get_all_folders_recursively(path='scratch',trim=False):
		
	if path[-1] == '/':
		path = path[:-1]
	
	basename_path = os.path.basename(path)
	
	folders = [path]
	#try:

	for i in os.walk(path,followlinks=True):
		if i[1]:
						
			folder_name = i[0]+'/'
			if trim:
				path_to_find = '/'+basename_path+'/'
				found_pos =  folder_name.find(path_to_find)
				
				if found_pos>=0:																
					folder_name = folder_name[found_pos+1:]
			
			for j in i[1]:
				folders.append(folder_name+j)		
									

	#except:
	#	raise Exception("{0} is not accessible".format(path))
	
	#contents = {'_all_' : files}
	#[contents[ext[1:]].append(base+ext) if contents.has_key(ext[1:]) 
	#		else contents.update({ext[1:]:[base+ext]})
	#		for base,ext in [os.path.splitext(f) for f in files]]
	return folders

#if ext_only then will only extract files that have end in '*.*' => useful for IMGT zipped files that may have an individual files folder 
#move_all_files_to_folder_path => if True, then any files created in any subdirectories will be moved to the 'folder_path', all subdirectories will be subsequently deleted
#WIL ALWAYS ZIP INTO THE FOLDER OF THE FIRST FILE !!
def zip_files(file_list,zip_name ='',method='ZIP'):
	
	zip_name = os.path.basename(zip_name)
	
	if not isinstance(file_list,list):
		file_list =[file_list]
	
	method = method.upper().strip()
	allowed_methods = ['ZIP','TAR.GZ','TAR.BZ2']
	if method not in allowed_methods:
		raise Exception("We only allow the following values for the parameter 'method': {0}".format(','.join(allowed_methods)))
	
	for each_file in file_list:
		if not os.path.isfile(each_file):
			raise Exception("The following file, {0}, is not a valid file".format(each_file))
		
	output_path = '/'.join(file_list[0].split('/')[:-1])
		
	
	if not zip_name:
		zip_name = 'compressed_files_{0}'.format(str(datetime.now()).replace(' ','').replace(':','').replace('-','').replace('_','').replace('/',''))
	
	if not(zip_name.endswith('.'+method.lower())):
		zip_name+='.'+method.lower()
	
	if method == 'ZIP':
		#use -j to ignore the folder structure of files 
		command = 'zip -j {2}/{0} {1}'.format(zip_name,' '.join(file_list),output_path)
		subprocess.call(command,shell=True)
	elif method =='TAR.GZ':
		#step1: copy all selected files to base folder ... we do this to ignore folder structure. that is zip files only, perhaps shoudl just be 'ADDING' to targball rather than doing it in one go
		for file in file_list:			
			subprocess.call("cp '{0}' .".format(file),shell=True)
		#step 2: zip all copied files 				
		command = 'tar -czvf {2}/{0} {1} --remove-files'.format(zip_name,' '.join([os.path.basename(file) for file in file_list]),output_path)
		subprocess.call(command,shell=True)		
	elif method == 'TAR.BZ2':
			#step1: copy all selected files to base folder ... we do this to ignore folder structure. that is zip files only, perhaps shoudl just be 'ADDING' to targball rather than doing it in one go
		for file in file_list:			
			subprocess.call("cp '{0}' .".format(file),shell=True)
		#step 2: zip all copied files 				
		command = 'tar -cjvf {2}/{0} {1} --remove-files'.format(zip_name,' '.join([os.path.basename(file) for file in file_list]),output_path)
		subprocess.call(command,shell=True)		
	else:
		raise Exception("HAVE NOT MADE AN ZIP FUNCTION FOR THIS METHOD YET: "+method)	
		
	return output_path+'/'+zip_name


#if ext_only then will only extract files that have end in '*.*' => useful for IMGT zipped files that may have an individual files folder 
#move_all_files_to_folder_path => if True, then any files created in any subdirectories will be moved to the 'folder_path', all subdirectories will be subsequently deleted
def unzip_files_in_folder(folder_path,ext_only=True,move_all_files_to_folder_path=True,keep_zip_filename=False):	
	
	#subfucntion for just using tar
	#for gzip files => decompress = "z", for bzip files, decompress = "j"
	#keep_extension => a hack for unzipping imgt .txz files. it does not store the sample name anymore, so we need to manually add it to the file 
	def use_untar(parent_folder,tar_file,decompress=''):
		#parameters => --wildcards expression => only unzip files matching expression
				#=> -C dir=> output result to the following directory
				#standard gunzip and extract => xzvf => extract, gzip, verbose, filedirectore		
		if parent_folder[-1]!='/':
			parent_folder+='/'		
		tar_folder_path = '/'.join(tar_file.split('/')[:-1])+'/'		
		
		#add a temp folder to extract files 
		tempfolder = tar_folder_path+'tempextract_'+str(datetime.now()).replace(' ','').replace(':','').replace('/','')			
		if tempfolder[-1]!='/':
			tempfolder+='/'
		untar_command = '''tar x{3}vf '{0}' -C '{2}' {1} '''.format(tar_file,"--wildcards '*.*'" if ext_only else '',tempfolder,decompress)		
		#make temp folder 		
		os.mkdir(tempfolder)#system("mkdir '"+tempfolder+"'")
		
		#untar files to temp folder
		subprocess.call(untar_command,shell=True)		
		
		file_substring = '.'.join(os.path.basename(tar_file).split('.')[:-1])
		
		#RENAME FILES USING THE NAEM FROM THE ZIPPED FILE 
		#WE ONLY DO THIS FOR .TXT FILES!!!..BASICALLY THIS IS OUR IMGT HACK
		if keep_zip_filename:	
			zipped_file_name ='.'.join(os.path.basename(tar_file).split('.')[:-1])
			for x in os.listdir(tempfolder):
				if x.endswith('.txt'):
					splits = x.split('.')
					prefix='.'.join(splits[:-1])+'_'
					os.rename(tempfolder+x,tempfolder+prefix+zipped_file_name+'.'+splits[-1])				
								
		if move_all_files_to_folder_path:
			#recursively find all files and move them 			
			subprocess.call("find '{0}' -type f -exec mv {{}} '{1}' \;".format(tempfolder,parent_folder),shell=True)
		else:
			#just move the results from the tar
			subprocess.call("mv '{0}'/* '{1}' ".format(tempfolder,parent_folder),shell=True)			
		#delete temp folder 
		subprocess.call("rm -r '"+tempfolder+"'",shell=True)
		
	
	e = ExperimentDirs()		 	
	#make sure folder provided exists 
	folder_path = e.get_valid_folder_loc(folder_path)
	#currently support extracting files from following zipped file extensions: 
	#currently support .zip, .gz, and .tar.gz/.tgz files, .tar files, .tar.bz2,.tbz,.tb2,.bz2 files 
	
	if folder_path[-1]!='/':
		folder_path+='/'
	
	
	#search all files in folder, identify compressed files 
	#next extract seperate .gz files from .tar.gz files 		
	#1 Get dict of all files in folder	
		#key of dict was the filetype.values in dict are list of files
	files_in_folder = dir_files(folder_path)
	tar_gz = []
	gz_files = []
	tar_files = []
	tar_bz = []
	tar_xz = []
	bz_files = []
	zip_filenames = []
	for file_types,files in files_in_folder.iteritems():
		if file_types == 'zip':
			zip_filenames = files
		if file_types == 'gz':
			tar_gz.extend([f for f in files if f.endswith('.tar.gz')])						
			gz_files = list(set(files)-set(tar_gz))
		elif file_types == 'bz2':
			tar_bz.extend([f for f in files if f.endswith('.tar.bz2')])						
			bz_files = list(set(files)-set(tar_bz))
		elif file_types == 'txz':
			tar_xz.extend(files)
		elif file_types == 'tgz':
			tar_gz.extend(files)
		elif file_types == 'tbz':
			tar_bz.extend(files)
		elif file_types == 'tb2':
			tar_bz.extend(files)
		elif file_types == 'tar':
			tar_files = files
		
	
	#first UNZIP *.ZIP files 
	#parameters => (n) never overwrite OR use (o) to always overwrite files
		#=> j -> extract all files to same folder 	
	if zip_filenames:		
		if keep_zip_filename:
			#append the zipped filename to the end of the file
			for each_file in zip_filenames:
				zipped_file_name = '.'.join(os.path.basename(each_file).split('.')[:-1])
				#make temp folder
				os.makedirs(folder_path+each_file+'_folder')			
				#unzip each individual file to that specific folder
				unzip_command = '''unzip -{0}{1} '{2}' {3} -d '{4}' '''.format('n','j' if move_all_files_to_folder_path else '',folder_path+each_file,'\*.*' if ext_only else '',folder_path+each_file+'_folder')							
				subprocess.call(unzip_command,shell=True)
				#read all of the created files in that folder
				for unzipped_files in os.listdir(folder_path+each_file+'_folder'):
					#only rename .txt files 
					#basically this is our IMGT hack for their current settings 
					if unzipped_files.endswith('.txt'):
						sp = unzipped_files.split('.')					
						new_name = folder_path+'.'.join(sp[:-1])+'_'+zipped_file_name+'.'+sp[-1]								
					else:
						new_name = folder_path+'/'+unzipped_files
					os.rename(folder_path+each_file+'_folder/'+unzipped_files,new_name)						
				shutil.rmtree(folder_path+each_file+'_folder')
		else:
			unzip_command = '''unzip -{0}{1} '{2}'\*.zip {3} -d '{2}' '''.format('n','j' if move_all_files_to_folder_path else '',folder_path,'\*.*' if ext_only else '')							
			subprocess.call(unzip_command,shell=True)
	
	#UNTAR DATA gzip data
	if tar_gz:		
		for each_targz in tar_gz:
			each_targz = folder_path+each_targz			
			use_untar(folder_path,each_targz,'z')
	#UNTAR tar data
	if tar_files:		
		for each_tar in tar_files:		
			each_tar = folder_path+each_tar			
			use_untar(folder_path,each_tar)
	
	if tar_xz:		
		for each_tarxz in tar_xz:
			each_tarxz = folder_path+each_tarxz			
			use_untar(folder_path,each_tarxz,'J')	
	
	if tar_bz:		
		for each_tarbz in tar_bz:
			each_tarbz = folder_path+each_tarbz			
			use_untar(folder_path,each_tarbz,'j')	
	#GUNZIP
	if gz_files:		
		subprocess.call('gunzip {0}'.format(' '.join(["'"+folder_path+f+"'" for f in gz_files])),shell=True)
	
	#BUZNIP
	if bz_files:		
		subprocess.call('bzip2 -d {0}'.format(' '.join(["'"+folder_path+f+"'" for f in bz_files])),shell=True)

#we really shoudl change this code so that it uses the zipfile and tarfile modules provided by python 
#file_name = > a single file or list of files to unzip 
def unzip_files(file_name,output_dir = None,copy_zip_file_as_prefix=False):						
	#if ext_only then will only extract files that have end in '*.*' => useful for IMGT zipped files that may have an individual files folder 
	ext_only = True
	#move_all_files_to_folder_path => if True, then any files created in any subdirectories will be moved to the 'folder_path', all subdirectories will be subsequently deleted
	move_all_files_to_folder_path=True
	if not isinstance(file_name,list):
		file_name = [file_name]
		
	current_time = str(datetime.now()).replace(' ','').replace(':','').replace('-','').replace('_','').replace('/','')
	for f in file_name:
		if not(os.path.isfile(f)):
			raise Exception('The provided file does not exist: '+f)
	
	if not output_dir:
		folder_path = '/'.join(file_name[0].split('/')[:-1])+'/'
	else:
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		folder_path = output_dir
		if folder_path[-1]!='/':
			folder_path+='/'
	
	sub_folder = folder_path+current_time+'/'
	os.makedirs(sub_folder)	
	
	for each_file in file_name:
		shutil.copyfile(each_file,sub_folder+os.path.basename(each_file))
	
	unzip_files_in_folder(sub_folder,ext_only,move_all_files_to_folder_path,copy_zip_file_as_prefix)
	for each_file in file_name:		
		if os.path.isfile(sub_folder+os.path.basename(each_file)):
			os.remove(sub_folder+os.path.basename(each_file))
	
	files_created = []
	
	#move all files to parent folder 
	for files in os.listdir(sub_folder):
		new_file_path = folder_path+files 
		current_file_path = sub_folder+files 
		os.rename(current_file_path,new_file_path)
		files_created.append(files)
	
	shutil.rmtree(sub_folder)
	return [folder_path+f for f in files_created]



class ExperimentDirs(object):
	"""Methods for user experiment directory management.
	_.get_dirs() returns experiment dirs in user's folder.
	_.make_dir(exp_name) creates a new experiment dir and returns its path.
	"""
	def __init__(self,proj=_projName_):
		self.basedir = self.get_basedir(proj)
		mkdir_p(self.basedir) # in case it doesn't yet exist
		# make sure it's accessible
		if not os.access(self.basedir, os.F_OK):
			raise Exception("{} is not accessible".format(self.basedir))

	@staticmethod
	def get_basedir(proj = _projName_):
		return os.path.join(scratch, proj)
	
	
	#given a path, attempts to determine its valid path 
	def get_valid_folder_loc(self,path):	
		if not path:#path provided is empty 
			raise Exception("A file path must be provided. Empty string passed in ")
	
		basedir = os.path.join(scratch, _projName_)
		#check if path is valid directly or file 
		if os.path.isdir(path):
			#since the folder is valid, then we just need to return this path 
			if path[-1]=='/': #remove last '/' 
				path = ''.join(list(path)[:-1])			
		elif os.path.isfile(path):
			#just return the pathlocation of this file (everything before the last '/')
			path = '/'.join(path.split('/')[:-1])							
		#file is not valid 			
		else:
			#assume that the file came from this project
			path = basedir+'/'+path
			
			#repeat steps above 
			#check if path is valid directly or file 
			if os.path.isdir(path):
				#since the folder is valid, then we just need to return this path 
				if path[-1]=='/': #remove last '/' 
					path = ''.join(list(path)[:-1])
					
			elif os.path.isfile(path):
				#just return the pathlocation of this file (everything before the last '/')
				path = '/'.join(path.split('/')[:-1])				
			else:
				raise Exception("Count not resolve the location of the path provided: "+path)
		return path#url_to_path(path)
		
			
			

	def get_dirs(self,exp_format_only=False):
		"""Return the user's experiment directories."""
		self.expdirs = []

		# define how to check experiment folder names. Not a ton of flexibility...
		expname_start = "{start}{delim}".format(**expDirParams)
		digitInds = [len(expname_start)]*2
		digitInds[1] += expDirParams['ndigits']

		#FOR NOW, JUST FIND ANY SUBFOLDER UNDER IGREP
		# find subdirectories matching the format
		for subdir in os.walk(self.basedir,followlinks=True).next()[1]:
			#when exp_format_only = True, then only find folders which match format 
			if exp_format_only:
				if (	len(subdir) > digitInds[1] and 
						subdir.startswith(expname_start) and 
						subdir[digitInds[0]:digitInds[1]].isdigit()):
					self.expdirs.append(subdir)
			else:
				self.expdirs.append(subdir)
		return self.expdirs

	def next_exp_num(self):
		"""Return next higher experiment number."""
		exp_dirs = sorted(self.get_dirs(exp_format_only=True))
		if len(exp_dirs):
		#for exps in exp_dirs[::-1]:
		#	try:
				#use try/except because not all directories will be in 'igrep' format (expDirParams)
			last_exp_dir = exp_dirs[-1]
			return int(last_exp_dir.split(expDirParams['delim'])[1])+1
			#except:
			#	continue
		#return 1
		else:
			return 1

	def make_dir(self, exp_name = time.strftime('%m%d%y_%H%M%S')):
		"""Create a new experiment directory,
		as defined by expDirParams and expDirTemplate"""
		dirname = expDirTemplate.format(exp_num = self.next_exp_num(), exp_name = exp_name)
		pathname = (os.path.join(self.basedir, dirname))
		mkdir_p(pathname)
		return pathname

'''
### --- Below here is an older class that may or may not be useful
job_dir = url_to_path(appsoma_api.environment_get_job_url())
class JobFiles(object):
	"""
	Accepts list of files, either urls or local paths and imports them into workspace.
	If url, pulls file into user's private directory using appsoma_api.resource_pull
		and creates a sym link in 'symlink_dir'.
	If path, just creates a link in 'symlink_dir'.
	self.paths is dict of <file name>:<full path of symlink>
	'symlink_dir' defaults to the current job's directory.
	Files can be added to an already instantiated JobFiles object with .add_files()
	Can set min and max # of files object is allowed to represent.

	Creates object with public attributes:
		paths: <dict> of {file name: file path}
		add_files: method for linking files and populating self.paths
	"""
	def __init__(self, import_these, symlink_dir = job_dir, primary_dir = scratch,
			min_nfiles = 1, max_nfiles = float('Inf')):

		self._realpath = primary_dir
		assert os.path.isdir(self._realpath)
		self._sympath = symlink_dir
		assert os.path.isdir(self._sympath)
		
		_, home, localpath = self._realpath.partition(appsoma_api.environment_get_home())
		assert home == appsoma_api.environment_get_home()
		self._localpath = localpath

		self._base_paths =  dict()
		self.paths =  dict()
		self._min_nfiles = min_nfiles
		self._max_nfiles = max_nfiles

		self.add_files(import_these)

	def add_files(self, files):
		#valid_files = self._validate_files(files) # tried to add validation, difficult with appsoma resource urls
		file_list = [files] if not isinstance(files, list) else files
		for f in file_list:
			t_file_name = os.path.basename(f)
			t_symln_path = os.path.join(self._sympath,t_file_name)

			if os.name == 'nt':  #@LOCAL
				shutil.copyfile(f, t_symln_path)
				self._base_paths[t_file_name] = f
			else:
				if os.access(f, os.R_OK):
					os.symlink(f, t_symln_path)
					self._base_paths[t_file_name] = f
				else:
					t_pullhere_full = os.path.join(self._realpath, t_file_name)
					t_pullhere = t_pullhere_full.partition(appsoma_api.environment_get_home())[2]
					appsoma_api.resource_pull(f,t_pullhere)
					assert os.path.exists(t_pullhere_full)
					os.symlink(t_pullhere_full, t_symln_path)

					self._base_paths[t_file_name] = t_pullhere_full

			self.paths[t_file_name] = t_symln_path
		assert self._min_nfiles <= len(self.paths) <= self._max_nfiles
'''

#This module is not meant to be run
if __name__ == "__main__":	
	print("The answer is 6.")
