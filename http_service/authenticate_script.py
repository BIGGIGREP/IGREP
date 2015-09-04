from simplepam import authenticate

user = sys.argv[1]
pass = sys.argv[2]

try:
	print simplepam.authenticate(user,pass)
except:
	print 0
