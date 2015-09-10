import json
import sys
import os 
import json
import subprocess


parent_of_proxy = os.path.dirname(os.path.dirname(os.path.abspath("__file__")))
sys.path.insert(0,os.path.join(parent_of_proxy,'common_tools'))# '../common_tools/') #add scripts in this file to the proxy list

import immunogrep_db_query_api as query
import immunogrep_useful_functions as useful

try:	
	from simplepam import authenticate
	use_posix = True
except:	
	use_posix = False
	print('IM USING A WINDOWS MACHINE AS MY SERVER WITHOUT POSIX')

def authenticate_user(u,p):
	if use_posix:
		res =  subprocess.check_output('sudo python authenticate_script.py '+u+' '+p,shell=True)
		if int(res)==0:
			return False
		else:
			return True
	else:
		print('im not doing authentication!!!!')
		return True	
	
def load_configuration_file():
	'''
		
		The .cfg file inside the http_service.py folder defines some default parameters we use.
		This will load the file, but if it doesnt exist, then it will set some default parameters.
		
	'''
	if not os.path.isfile('igrep_services.cfg'):
		print('There is no configuration file defining, will set default parameters')
		igrep_params = {
			'igrep_mongoproxy_path':'localhost:6200'
		}
	else:
		with open('igrep_services.cfg','r') as ff:			
			igrep_params = json.load(ff)
	return igrep_params

#def walklevel(some_dir, level=1,followlinks=False):
#	some_dir = some_dir.rstrip(os.path.sep)
#	assert os.path.isdir(some_dir)
#	num_sep = some_dir.count(os.path.sep)
#	for root, dirs, files in os.walk(some_dir,followlinks=followlinks):
#		yield root, dirs, files
#		num_sep_this = root.count(os.path.sep)
#		if num_sep + level <= num_sep_this:
#			del dirs[:]
def splitall(path):
	"""
		take a path split it into lists: 
		/A/B/C/d.txt => [A,B,C,d.txt]
	"""
	allparts = []
	while 1:
		parts = os.path.split(path)
		if parts[0] == path:  # sentinel for absolute paths
			allparts.insert(0, parts[0])
			break
		elif parts[1] == path: # sentinel for relative paths
			allparts.insert(0, parts[1])
			break
		else:
			path = parts[0]
			allparts.insert(0, parts[1])
	return allparts

igrep_params = load_configuration_file()

#THIS MODULE SHOULD SERVE AS A HELPER SCRIPT FOR MAKING OUR IGREP APPS
#JAVASCRIPT FUNCTIONS CAN MAKE AJAX CALLS TO THE PROXY WHICH CAN SUBSEQUEENTLY CALL ANY OF THE
#FUNCTIONS DEFINED IN HTTPIGREPHANDLER

class HtppIgrepHandler():
	def __init__(self,fxn_name,args=[],kargs={}):
		self.fxn = fxn_name
		self.args=args
		self.kargs=kargs
		self.parent_of_proxy = os.path.dirname(os.path.dirname(os.path.abspath("__file__")))
	def run(self):
		#run function passed in by user			
		fxn_to_eval = 'self.{0}(*{1},**{2})'.format(self.fxn,json.dumps(self.args),json.dumps(self.kargs))
		fnx_ans = eval(fxn_to_eval)
		return fnx_ans
	def example(self):
		return 'heere we go'   
	def sum_vals(self,a,b=3):
		print a
		return a+b		
	def authenticate(self,username,password):
		"""
			
			Function for authenticating the current user accessing the computer
			It will check whether the user logging in is a current user in the computer
			
		"""
		return authenticate_user(username,password)

	def get_default_mongoproxy_address(self):
		return igrep_params['igrep_mongoproxy_path']

	def get_apps_list(self):
		"""
			
			Searches the 'apps' folder found in the parent. We want to return a list of possible apps 
			that a user can run. This function will only search three levels of hierarchy of 'projects/app topics'. 
			For each level, we will search for folders containing index.html files. Any folder with an index.html folder
			within the apps folder will be assumed to be an APP for analysis. After the maximum levels of hierarchy, 
			all other subfolders containing html pages will be grouped into the last hierarchy group. 
			
			Example:
				/apps/gsaf download/index.html => there is an app (gsaf download) at LEVEL 1
				/apps/database queries/queries/simple database query/index.html  => there is an app (simple database) at LEVEL 3
				/apps/automated pipelines/ngs processing/index.html => there is an app (ngs processing) at LEVEL 2
				/apps/annotation apps/vgene analysis/pairing/single chain pairing/index.html => in this hypothetical situation the app 'single chain pairing' is found at level 4, so this app will instead be grouped under vgene analysis and the 'pairing' hierarchy is ignored
		
			For each app, a file description.txt and keywords.txt can be used as a brief description for each app and keywords that can be used to search for apps

		"""
		max_hierarchy_level = 3
		app_home = os.path.join(self.parent_of_proxy,'igrep-webpage','apps')
		apps_found = []
		f = []
		d = []
		
		#search for all files named 'index.html' in directories within the apps folder
		all_files = sorted([os.path.relpath(dirpath,start=app_home) for (dirpath,dirnames,filenames) in os.walk(app_home,followlinks=True) for f in filenames if f=='index.html' and dirpath!=app_home])
		app_list = []
		#stores the name of all directories up until the max hierarchy
		directory_dictionary = {}			
		for f in all_files:
			#split the filename into a list
			splitf = splitall(f)			
			
			#search all of the found apps for description files and keywords files
			if os.path.isfile(os.path.join(app_home,f,'description.txt')):
				with open(os.path.join(app_home,f,'description.txt')) as r:
					desc = r.read()					
			else:
				desc = ''
									
			classes = ['-'.join(splitf[:i+1]).replace(' ','-') for i in range(max_hierarchy_level) if i<len(splitf)]
						
			directory_dictionary['.'.join(splitf[:max_hierarchy_level])] = 1
										
			kwd = [splitf[-1]]
			if os.path.isfile(os.path.join(app_home,f,'keywords.txt')):				
				with open(os.path.join(app_home,f,'keywords.txt')) as r:
					for l in r.readlines():
						l = l.strip('\r\n\t')
						kwd.extend(l.split(','))				
			#we hardcode this part because its for the html link 
			html_path = '/'+os.path.basename(app_home)+'/'+'/'.join(splitf)
			app_list.append({'path':html_path,'keywords':kwd,'description':desc,'classes':classes,'name':splitf[-1],'main-group':splitf[0]})
			
								
		directory_dictionary = dict(useful.DotAccessible(directory_dictionary))
		
		return {'directory_dictionary':directory_dictionary,'apps':app_list}

	def get_metadata_on_proxy(self,proxy_path=igrep_params['igrep_mongoproxy_path']):
		"""
			
			Queries the database and returns information necessary for running javascript functions.
			1) Returns metadata for experiments stored in the database.
			2) Returns information regarding the current user running APP
			
	        ===============						================						================
	        **Input**							**Description**							**Default value**
	        ---------------						----------------						----------------
	        proxy_path							String defining the address of			biotseq.icmb.utexas.edu 
												the IGREP mongo db database												
			---------------						----------------						----------------			
	        
			===============						================						================

			:General flow of functions:
			
				1) connect to the mongo database via our mongoproxy
				2) Run a query for getting current user's information (write_access, administrator, etc)
				3) Run a query searching for experiment metadata that the user has access to 
				
			Fields Returned by function

	        ===============						================						================
	        **Output**							**Description**							**Format**
	        ---------------						----------------						----------------
	        users								Information about current user			 Dictionary containing the following keys:
																						 user,write_access,administrator,name
			---------------						----------------						----------------
			metadata_from_experiment_database   All metadata associated with experiments List of dictionaries
	        ===============						================						================
		
		"""
				
		# CREATE an instance of the class made to query data via proxy from database		
		exp_collection_data = query.RunQuery(proxy_path=proxy_path) 		
		# REQUEST information of the current user running program
		# returns user information from database including name, lab, email, administrator access 
		users = exp_collection_data.get_user_access_info()._return(to_file=False).next()
		if not users:
			return {}
		users = users[0]
		#if not(users['user']): #the user was not found in the database
		#	
		#.communicate_javascript_run_function('NoAccess',[])		
		#	return {}
			#os._exit(1) #exit program
								
		# REQUEST all records from experiments collection;
		# first function sets up the parameters/list of functions that will be queried once run; 
		# the second function calls proxy, runs requested functions, and returns results as a list of dictionaries.
		# each document is an index in list, and all metadata for that document is stored as dictionary		
		#metadata_from_experiment_database =  exp_collection_data.get_exp_docs(write_access_only=True)._return()
		metadata_from_experiment_database =  []
		for list_of_docs in exp_collection_data.get_exp_docs()._return(to_file = False,chunk_size = 1000):
			metadata_from_experiment_database.extend(list_of_docs)						
		with open('metadatafile.txt','w') as outme:
			outme.write(json.dumps(metadata_from_experiment_database,indent=4)+'\n')
			outme.write(json.dumps(users))
			
		return {'metadata':metadata_from_experiment_database ,'user_data':users}		
	

