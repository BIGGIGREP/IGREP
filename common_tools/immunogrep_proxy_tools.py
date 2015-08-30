## Tools for working with the proxy server between Appsoma and the Mongo database.

import inspect
import os

import json
import urlparse
import httplib
import ssl

from urllib2 import HTTPError
from httplib import BadStatusLine

import urllib2
import urllib

def on_proxy():
	"""Check whether the code calling this function is running on the proxy or not.
	Warning: this assumes that some function that is always running while the proxy is active is named 'mongoproxy.py'.
	"""
	stack_list = inspect.stack()
	stack_function_names = [os.path.basename(frame[1]) for frame in stack_list]
	if 'mongoproxy.py' in stack_function_names:
		return True
	else:
		return False
		
		
def http(url, params={}, data="", action="GET", headers={}, progressCallback=None, toFilename=None, retrys=0, retryDelay=5, returnHeaders=False,onlySendCommand=False):
	attempts = 0
	while attempts <= retrys:
		try:
			parsed = urlparse.urlparse(url)
			if url.startswith("https://"):
				# WARNING: This does not do any checking of the SSL certificate
				h = httplib.HTTPSConnection(parsed.netloc)
			else:
				h = httplib.HTTPConnection(parsed.netloc)

			# The URL may have a query string which is combined with the params
			if parsed.query:
				inlineQueryAsDict = urlparse.parse_qs(parsed.query)

				# FLATTEN because parse_qs annoyingly returns a list for each value
				for k in inlineQueryAsDict:
					inlineQueryAsDict[k] = inlineQueryAsDict[k][0]

				# UNION the sent in params with the query string
				params = dict(params.items() + inlineQueryAsDict.items())

			headers["Content-Length"] = str(len(data))
			headers["User-Agent"] = "fake browser"
			if data and data != "" and ("Content-Type" not in headers) and not(action == "GET" or action == "DELETE"):
				raise HTTPError(url, 0, "Trying to send data without Content-Type to url:" + url, None, None)

			paramString = ""
			if params:
				paramString = "?"+urllib.urlencode(params).encode('utf-8')

			try:
				escapePath = urllib.quote(parsed.path).encode('utf-8')
				h.request(action, escapePath + paramString, data, headers)
				resp = h.getresponse()
			except ssl.SSLError as e:
				raise HTTPError(url, 406, "SSL failed: " + url + ", " + str(e), None, None)
			except BadStatusLine as e:
				raise HTTPError(url, 404, "URL returned bad status:" + url + ", " + str(e), None, None)
			except Exception as e:
				raise HTTPError(url, 0, "Can not reach url: " + url + ", " + str(e), None, None)

#			if resp.status != 200:
#				raise HTTPError(url, resp.status, "Return code not 200", None, None)
			
			if onlySendCommand:
				return None			
			readData = ""
			f = None			
			
			try:				
				if toFilename and resp.status == 200:
					try:
						os.makedirs('/'.join(toFilename.split('/')[:-1]))
					except:
						pass
					f = open(toFilename, "wb")				

				totalSize = int( resp.getheader("content-length", 0))
				blockSize = 1024 * 1024
				count = 0
				lengthRead = 0
				while True:
					chunk = resp.read(blockSize)
					if chunk and progressCallback:
						lengthRead += len(chunk)
						progressCallback(url, lengthRead, totalSize)
					if not chunk:
						break
					if f:
						f.write(chunk)
					else:
						readData += chunk

					count += 1
			except Exception as e:
				raise HTTPError(url, resp.status, "Bad read or write: " + url + " " + str(e), None, None)
			finally:
				if f is not None:
					f.close()
			print resp.status
			if resp.status != 200:
				print('error is here!!!!!')
				raise HTTPError(url, resp.status, readData, None, None)

			if returnHeaders:
				replyHeaders = {}
				print resp.getheaders()
				for k, v in resp.getheaders():
					replyHeaders[k.lower()] = v
				return readData, replyHeaders
			else:
				return readData

		except Exception as e:
			if attempts >= retrys:
				raise e
			attempts += 1
			if retryDelay:
				time.sleep(retryDelay)

	# We should never get here, but in case we throw an exception
	raise HTTPError(url, 0, "http fatal: " + url, None, None)



