import json
import sys
import sys
sys.path.insert(0,'../common_tools/') #add scripts in this file to the proxy list
import immunogrep_db_query_api as query

try:
	from simplepam import authenticate as authenticate_user
	#pass
except:
	from fakepam import authenticate as authenticate_user
	print('IM USING A WINDOWS MACHINE WITHOUT SIMPLEPAM')


#THIS MODULE SHOULD SERVE AS A HELPER SCRIPT FOR MAKING OUR IGREP APPS
#JAVASCRIPT FUNCTIONS CAN MAKE AJAX CALLS TO THE PROXY WHICH CAN SUBSEQUEENTLY CALL ANY OF THE
#FUNCTIONS DEFINED IN HTTPIGREPHANDLER

class HtppIgrepHandler():
	def __init__(self,fxn_name,args,kargs):
		self.fxn = fxn_name
		self.args=args
		self.kargs=kargs
	def example(self):
		return 'heere we go'   
	def sum_vals(self,a,b=3):
		print a
		return a+b		
	def authenticate(username,password):
		"""

			Function for authenticating the current user accessing the computer
			It will check whether the user logging in is a current user in the computer

		"""
		return authenticate_user(username,password)

	def get_metadata_on_proxy(self,proxy_path=''):
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
		if not(users['user']): #the user was not found in the database
		#	
		.communicate_javascript_run_function('NoAccess',[])		
			return {}
			#os._exit(1) #exit program
								
		# REQUEST all records from experiments collection;
		# first function sets up the parameters/list of functions that will be queried once run; 
		# the second function calls proxy, runs requested functions, and returns results as a list of dictionaries.
		# each document is an index in list, and all metadata for that document is stored as dictionary		
		#metadata_from_experiment_database =  exp_collection_data.get_exp_docs(write_access_only=True)._return()
		metadata_from_experiment_database =  []
		for list_of_docs in exp_collection_data.get_exp_docs()._return(to_file = False,chunk_size = 1000):
			metadata_from_experiment_database.extend(list_of_docs)						
		return {'metadata':metadata_from_experiment_database ,'user_data':users}		
	
	def run(self):
		#run function passed in by user			
		fxn_to_eval = 'self.{0}(*{1},**{2})'.format(self.fxn,json.dumps(self.args),json.dumps(self.kargs))
		fnx_ans = eval(fxn_to_eval)
		return fnx_ans