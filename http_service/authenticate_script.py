import sys
from simplepam import authenticate

user = sys.argv[1]
passwd = sys.argv[2]

if authenticate(user,passwd):
	#0 => means success in exit code
	sys.exit(0)	
else:
	#=>anything else is an error of somekind
	sys.exit(1)
	
