#THIS VERSION OF MONGOPROXY HAS BEEN MODIFIED TO ALLOW THE NEW QUERY FUNCTIONS
from SocketServer import ThreadingMixIn
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
import threading
import urllib
import urllib2
import json
import sys
import time
import signal
import errno
from bson.json_util import dumps as dumps_mongo
from bson.json_util import loads as loads_mongo
from SimpleXMLRPCServer import SimpleXMLRPCServer
from SimpleXMLRPCServer import SimpleXMLRPCRequestHandler
import os
import traceback
import shutil

#if not(os.path.isdir('appsoma_scripts_here')):
#	os.system("mkdir appsoma_scripts_here")
	
sys.path.insert(0,'../common_tools/') #add scripts in this file to the proxy list
# from contextlib import contextmanager
import immunogrep_handle_proxy as immunogrep_proxy

password_dict = {}

import __builtin__

global modules_imported_by_proxy
modules_imported_by_proxy = sys.modules.keys()

query_objects = {}
__builtin__.query_objects = query_objects

with open( "mongo_password" ) as f:
	password = f.read().strip()
	password_dict = json.loads(password)
	



# The purpose of this context manager is to retrieve and load a module from its Appsoma URL,
# then reload the module in case updates have occured. The context manager handles the automatic reload. 
class RemoteModule():
	'''Context manager that loads a module from a remote location and reloads the module on exit in case it's changed.
	'''
	
	def __init__(self, module_location,new_name=None,import_module=False):

		self.basename = "appsoma_scripts_here/"
		if new_name:
			self.basename += new_name
		else:
			self.basename += os.path.basename(module_location)
		
		urllib.urlretrieve(module_location, self.basename) #downloaded file
		with open(self.basename) as downloaded_file:#lets make sure the file we downloaded is correct/downloaded correctly			
			try:
				first_line = json.loads(downloaded_file.readline().strip())
			except: #if that doesnt work, then thats OK becuase its probalby not a dictionry string
				first_line = {}

			if "error" in first_line:
				raise Exception("\n\n"+json.dumps(first_line)+"\n\n")

		if import_module:
			print self.basename
			self.module = __import__(os.path.splitext(os.path.basename(self.basename))[0])
			reload(self.module)
		else:
			self.module=None

	def __enter__(self):
		return self.module

	def __exit__(self, exception_type, exception_value, traceback):
		pass
		#reload(self.module)

#def set_db_globals():
#	########################################
#	# This is Benni connecting setting up Mongo connections here in the proxy.
#	# Please yell at me if it doesn't work, I'm fucking around, y'all.
#	########################################	
#	import immunogrep_igdbtools as igdbconnect
#	try:
#		_ig_.close()
#		_ig2.close()
#	except:
#		print "cannot close connection"
#		pass	
#	[immunogrep_proxy.db_writer,_ig_]=igdbconnect.connectToIgDatabase("writer",password_dict["writer"])
#	#immunogrep_proxy.db_writer = db_writer # This should set the global (to immunogrep_proxy) variable "db_writer" to the writer connection
#	[immunogrep_proxy.db_reader,_ig2_]=igdbconnect.connectToIgDatabase("reader",password_dict["reader"])
#	#immunogrep_proxy.db_reader = db_reader # This should set the global (to immunogrep_proxy) variable "db_reader" to the reader connection

def return_error_string(e):
	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]					
	error = "There was an error: \n"+str(e)
	error+=	"\nLine Number: {0} , type: {1} , fname: {2}\n".format(str(exc_tb.tb_lineno),str(exc_type),str(fname))
	tb_error = traceback.format_exc()
	error+=str(tb_error)
	return error 

	
def UpdateFiles():
	#first clear any modules currently imported
	current_modules = sys.modules.keys()
	for m in current_modules:
		if m not in modules_imported_by_proxy:
			del(sys.modules[m]) #this will delete any modules imported by downstream files we use in our appsoma_scripts
	
					
	with RemoteModule(handle_proxy_file_name_http,list_of_downloads[handle_proxy_file_name],True) as new_module:
		print "updateing"
		if new_module:
			print "new"
		global immunogrep_proxy
		immunogrep_proxy = new_module

# Threaded http server
#############################################################################################################




class HttpHandler( BaseHTTPRequestHandler ):
	def do_HEAD(self):
		pass

	def do_GET(self):
#		load_code()

		reply_str = "Bad command"
		self.path = urllib.unquote( self.path )


		#reply_str = json.dumps( immunogrep_proxy.handle_GET( password, self.path ) )

		self.send_response( 200 )
		self.send_header( "Content-length", str(len(reply_str)) )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()
		#self.wfile.write( reply_str )


	def do_DELETE(self):
		pass

	def do_POST(self):
		self.path = urllib.unquote( self.path )
		length = int( self.headers["Content-Length"] )
		data = self.rfile.read( length )
		#data_dict = json.loads(data)
		data_dict = loads_mongo(data)
		content_path = 'temp_file_query.txt'	
		try:			
			#run a database function			
			with open(content_path,'w') as out:
				immunogrep_proxy.handle_POST( password_dict, self.path, data_dict,out)			
		except Exception as e:	
			print('An error occurred duing query processing')		
			error_message = return_error_string(e)
			print(error_message)
			#send error back to proxy
			self.send_error(405,error_message)
			return	
														
		self.send_response( 200 )
		self.send_header( "Content-type", "application/json" )		
		self.send_header( "Access-Control-Allow-Origin", "*" )		
		self.end_headers()		
		#transfer the file over proxy
		if os.path.isfile(content_path):
			with open(content_path, 'rb') as content:
				shutil.copyfileobj(content, self.wfile)
		else:
			self.wfile.write('')

	def do_PUT(self):
		pass

	def do_OPTIONS(self):
		pass

	def log_message(self, format, *args):
		return

class ThreadedHTTPServer( ThreadingMixIn, HTTPServer ):
	pass

class HttpServerThread(threading.Thread):
	def __init__(self,port):
		threading.Thread.__init__(self)
		self.port = port
		self.stop_flag = False
		
	def stop(self):
		self.stop_flag = True

	def run(self):
		try:
			server = ThreadedHTTPServer( ('',self.port), HttpHandler )
			server.socket.settimeout(1)	# Important: Without a timeout, handle_request will never return!
			while not self.stop_flag:
				try:
					server.handle_request()
				except Exception as e:
					print str(e), "Server handle request problem"
			print "Stopping http server"
			server.server_close()
			print "Http server stopped"
				# do not use .shutdown() since that is meant to be called from a 
				# separate thread, and this thread handles all it's business...
			sys.exit(0)
		except IOError as e:
			if e.errno != errno.EINTR and e.errno != errno.EPIPE and e.errno != errno.EPERM:
				print str(e), "HTTP server IO error in onRun"
		except Exception as e:
			print str(e), "HTTP server thread onRun"


port = 6200


########### This section is added to test stuff##################################
#import igdbconnect
#[db_reader,_ig_]=igdbconnect.connectToIgDatabase("reader",password_dict["reader"])
#__builtin__.db_reader=db_reader
##################################################################################

print "STARTING NEW MONGOPROXY on ", port

try:
	thread_http_server = HttpServerThread(port)
	thread_http_server.start()
except KeyboardInterrupt:
	thread_http_server.stop()
	print "DIE"
	sys.exit(0)

def signal_handler(signal, frame):
	global thread_http_server
	thread_http_server.stop()
        print('You pressed Ctrl+C!')

signal.signal(signal.SIGINT, signal_handler)



print "SERVING"
#set_db_globals()

