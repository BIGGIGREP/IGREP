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
	
	def do_GET(self):		
		self.homepage_path = 'igrep-apps'		
		if self.path.startswith('/'):
			self.path=self.path[1:]		
		if self.path=='igrep':
			reply_str = "it2sasstart"	
		elif self.path==self.homepage_path:			
			#html homepage		
			with open('../apps/homepage/index.html') as w:
				reply_str = w.read()					
		elif self.path.startswith(self.homepage_path+'/'):
			#load other app webpages
			folder_path = self.path[len(self.homepage_path):]
			if folder_path[-1]!='/':
				folder_path+='/'				
			with open('../apps'+folder_path+'index.html') as w:
				reply_str = w.read()
		else:
			reply_str ='notsure'
				
		self.send_response( 200 )
		self.send_header( "Content-length", str(len(reply_str)) )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()
		self.wfile.write( reply_str )
	
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
	print "DIE"
	sys.exit(0)

def signal_handler(signal, frame):
	global thread_http_server
	thread_http_server.stop()
	print('You pressed Ctrl+C!')

signal.signal(signal.SIGINT, signal_handler)

print "SERVING"

