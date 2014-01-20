#!/usr/bin/env python

import sys
import string
import commands
true  = 1
false = 0

def do (command):
	print command
	output = commands.getoutput(command)
	if (len(output) > 0): print output

binDir     = "~/py/bin"
currentDir = commands.getoutput("pwd").strip()

args = sys.argv[1:]
remove = False
if (args[0] == "--remove"):
	remove = True
	args.pop(0)

binDir = args.pop(0)

for f in args:

	if f.endswith (".py"):
		fPy   = f
		fNoPy = f[:-3]
	else:
		fPy   = "%s.py" % f
		fNoPy = f

	fileHere  = "%s/%s" % (currentDir, fPy)
	fileInBin = "%s/%s" % (binDir,     fNoPy)

	if (remove):
		do ("rm %s" % (fileInBin))
	else:
		do ("ln -s %s %s" % (fileHere,fileInBin))
		do ("chmod +x %s" % (fileInBin))
