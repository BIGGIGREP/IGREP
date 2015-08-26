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
import os
import traceback



class HttpHandler( BaseHTTPRequestHandler ):
	def do_HEAD(self):
		pass

	def do_GET(self):
    reply_str = ""
		self.send_response( 200 )
		self.send_header( "Content-length", str(len(reply_str)) )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()
		self.wfile.write( reply_str )
    #for monogproxy:
      #zip all files created
      #then send wfie.write the zipped file
      #then unzip on client side
    #f = open( filename )
    #l = get the ength of the file
    #while sent < l:
    #  to_send = f.read( 100 )
		#self.wfile.write( to_send )

	def do_DELETE(self):
		pass

	def do_POST(self):
		self.path = urllib.unquote( self.path )
		length = int( self.headers["Content-Length"] )
		data = self.rfile.read( length )

    reply_str = "{}"
		self.send_response( 200 )
		self.send_header( "Content-type", "application/json" )
		self.send_header( "Content-length", str(len(reply_str)) )
		self.send_header( "Access-Control-Allow-Origin", "*" )
		self.end_headers()
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


port = 5998
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

