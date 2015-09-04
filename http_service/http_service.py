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
import os
import traceback
import handle_http_service_calls
import re
import httplib
import urlparse
import shutil


def return_error_string(e):
	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]					
	error = "There was an error: \n"+str(e)
	error+=	"\nLine Number: {0} , type: {1} , fname: {2}\n".format(str(exc_tb.tb_lineno),str(exc_type),str(fname))
	tb_error = traceback.format_exc()
	error+=str(tb_error)
	return error 


class HttpHandler( BaseHTTPRequestHandler ):
	
	def do_HEAD(self):
		pass
	
	def get_reply_str(self):
		self.reply_string = None
		self.reply_path = None
		"""
			
			Process self.path in http request to figure out what html files or figures to server to webpage
		
		"""
		# CURRENT ROUTES
		# /apps/:appName
		# HARDCODED ROUTES
			#IGREP/igrep-webpage/assets/images/
			#IGREP/igrep-webpage/assets/scripts/		
			#IGREP/igrep-webpage/apps/

			#to get to harcoded routes from server (http_service.py), use 
				#../igrep-webpage/assets/		
		#with open('test.txt','a') as f:
		#	f.write('original path '+self.path+'\n')
		#print('original path',self.path)
		orig_path  =self.path
		self.homepage_path = 'igrep-webpage'
		self.image_path = 'images'
		self.app_path = 'apps'
		self.scripts_path = 'scripts'
		self.styles_path = 'styles'
		self.parent_of_proxy = os.path.dirname(os.path.dirname(os.path.abspath("__file__")))
		self.reroute = os.path.join(self.parent_of_proxy,self.homepage_path,'assets')
		self.reroute_apps = os.path.join(self.parent_of_proxy,self.homepage_path)
		
		print(self.parent_of_proxy)
		#remove the '/' or '\' from beginning of path
		self.path = self.path.lstrip('/\\')
		if self.path == '':
			self.path = os.path.join(self.parent_of_proxy,self.homepage_path)
		elif self.path == self.homepage_path:
			self.path = os.path.join(self.parent_of_proxy,self.homepage_path)
		elif self.path == self.homepage_path+'/' or self.path == self.homepage_path+'\\':
			self.path = os.path.join(self.parent_of_proxy,self.homepage_path)			
		elif re.match(self.homepage_path+'[/\\\\].',self.path):			
			self.path = os.path.join(self.parent_of_proxy,self.path)
		elif re.match(self.image_path+'[/\\\\].',self.path):			
			self.path = os.path.join(self.reroute,self.path)
		elif re.match(self.scripts_path+'[/\\\\].',self.path):						
			self.path = os.path.join(self.reroute,self.path)
		elif re.match(self.styles_path+'[/\\\\].',self.path):
			self.path=self.path.lstrip('/\\\\')
			self.path = os.path.join(self.reroute,self.path)
		elif re.match(self.app_path+'[/\\\\].',self.path):
			self.path=self.path.lstrip('/\\\\')
			self.path = os.path.join(self.reroute_apps,self.path)
		elif self.path.strip('/\\')=='igrep':
			self.reply_string = "it2sasstart"	
			return
			#return reply_str		
		# Someday you might have other routes				
		#print('newpath',self.path)
		#with open('test.txt','a') as f:
		#	f.write('original path: '+self.path+'\n')
		
		#self.path has been re-routed, now search f
		if os.path.isfile(self.path):
			#with open(self.path,'rb') as w:
			#	reply_str = w.read()		
			self.reply_path = self.path
		elif os.path.isdir(self.path):						
			if os.path.isfile(os.path.join(self.path,'index.html')):
				self.reply_path = os.path.join(self.path,'index.html')
				#with open(self.path+'index.html') as w:
					#reply_str = w.read()
			else:
				#self.reply_string=''
				raise Exception('The following folder provided does not have an index.html file: '+self.path)				
		else:
			#self.reply_string=''
			raise Exception('The following path does not exist: '+self.path)			
		
		#return reply_str

	def do_GET(self):		
		try:
			#reply_str = self.get_reply_str()
			self.get_reply_str()
			
			self.send_response( 200 )	
			#self.send_header( "Content-length", str(len(reply_str)) )
			self.send_header( "Access-Control-Allow-Origin", "*" )
			self.end_headers()
			#self.wfile.write( reply_str )
			if self.reply_string!=None:
				self.wfile.write(self.reply_string)
			elif self.reply_path!=None:				
				self.reply_path =  os.path.abspath(self.reply_path)
				shutil.copyfileobj(open(self.reply_path,'rb'), self.wfile)
		except Exception as e:			
			self.send_error(401,str(e))
			print('Error in html request: '+str(e))
			print return_error_string(e)
	
	def do_DELETE(self):
		pass

	def do_POST(self):
		self.path = urllib.unquote( self.path )
		length = int( self.headers["Content-Length"] )
		#user passes in data/a function they want to run
		data = self.rfile.read( length )
		if not(data):
			return		
		#function requested by user
		try:
			#we load the function as a dict
			data_as_dict = json.loads(data)		
			fxn_to_run = data_as_dict['module']
			#keywordarguments 
			if 'kwargs' in data_as_dict:
				kargs = data_as_dict['kwargs']
			else:
				kargs = {}
			#arguments
			if 'args' in data_as_dict:
				args = data_as_dict['args']
			else:
				args = []
			init = handle_http_service_calls.HtppIgrepHandler(fxn_to_run,args,kargs)
			reply_str = json.dumps(init.run())			
		except Exception as e:
			error_message = return_error_string(e)
			print(error_message)
			#send error back to proxy
			self.send_error(405,error_message)
			return	
					
		#response = 200 => no errors found
		self.send_response( 200 )
		self.send_header( "Content-type", "application/json" )
		self.send_header( "Content-length", str(len(reply_str)) )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()
		#send answer back to webpage/ajax call or proxy call/any http call
		self.wfile.write(reply_str)		
		
	def do_PUT(self):
		pass

	def do_OPTIONS(self):
		self.send_response( 200 )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()

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


port = 6100
print "STARTING on ", port

try:
	thread_http_server = HttpServerThread(port)
	thread_http_server.start()
except KeyboardInterrupt:
	thread_http_server.stop()
	print("DIE")
	sys.exit(0)

def signal_handler(signal, frame):
	global thread_http_server
	thread_http_server.stop()
	print('You pressed Ctrl+C!')

signal.signal(signal.SIGINT, signal_handler)

print("SERVING")

