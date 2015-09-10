import sys
from simplepam import authenticate

user = sys.argv[1]
passwd = sys.argv[2]

if authenticate(user,passwd):
	print 1
else:
	print 0 
