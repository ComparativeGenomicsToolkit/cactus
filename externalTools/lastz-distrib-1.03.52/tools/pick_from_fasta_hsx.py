#!/usr/bin/env python
"""
Select a subset of sequences from a fasta file indexed by an hsx file
---------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys
from hsx_file import HsxFile


def usage(s=None):
	message = """
pick_from_fasta_hsx hsx_file [--names=<file>] [name1 name2 ...]
  --names=<file>  read sequence names from a file
  --nowarn        don't warn about sequences that aren't found
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	##########
	# parse the command line
	##########

	hsxFileName   = None
	seqNames      = []
	warnOnMissing = True
	showProgress  = False
	debug         = []

	args = sys.argv[1:]
	while (len(args) > 0):
		arg = args.pop(0)
		val = None
		fields = arg.split("=",1)
		if (len(fields) == 2):
			arg = fields[0]
			val = fields[1]
			if (val == ""):
				usage("missing a value in %s=" % arg)

		if (arg == "--names") and (val != None):
			f = file(val)
			seqNames += [line.strip() for line in f]
			f.close()
		elif (arg == "--nowarn") and (val == None):
			warnOnMissing = False
		elif (arg == "--progress") and (val == None):
			showProgress = True
		elif (arg == "--debug") and (val == None):
			debug += ["debug"]
		elif (arg == "--debug") and (val != None):
			debug += [val]
		elif (arg.startswith("--")):
			usage("unknown argument: %s" % arg)
		elif (hsxFileName == None) and (val == None):
			hsxFileName = arg
		elif (val == None):
			seqNames += [arg]
		else:
			usage("unknown argument: %s" % arg)

	if (hsxFileName == None): usage("you must give me an hsx file!")
	if (seqNames    == []):   usage("you must give me some sequence names!")

	##########
	# fetch the sequences
	##########

	hsx = HsxFile(hsxFileName,debug=debug)
	for name in seqNames:
		seq = hsx.get_sequence(name)
		if (seq != None):
			print seq
			if (showProgress):
				print >>sys.stderr, name
		elif (warnOnMissing):
			print >>sys.stderr, "WARNING: %s not found" % name
	hsx.close()


if __name__ == "__main__": main()

