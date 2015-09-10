
import json
import urllib2



def getFullName(user_auth_key=None):	
	print('YOU NEED TO FIX AUTHENTICATION!!!')
	username='cchrysostomou'
	#auth = json.loads( urllib2.urlopen( "https://appsoma.com/users/key_authenticate?key="+user_auth_key ).read() )		   
	#username = auth['user_name'] #appsoma_api.environment_get_username()
	#print username
	#if username in Administrators:
	#	admin = True
	#	user = True
	#	name = Administrators[username]["name"] 	
	#	lab = Administrators[username]["lab"]
	#elif username in Members:
	#	admin = False
	#	user = True
	#	name = Members[username]["name"]
	#	lab = Members[username]["lab"]
	#else:
	#	admin = False
	#	user = False
	#	name = "nouser_zz#"
	#	lab = ""
	
	#return {"user":user,"admin":admin,"name":name,"lab":lab}
	return username