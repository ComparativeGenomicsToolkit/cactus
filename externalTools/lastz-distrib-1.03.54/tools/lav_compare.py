#!/usr/bin/env python
"""
Compare two lav files, reporting differences but ignoring some trivial ones
---------------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys

def usage(s=None):
	message = """
lav_diff lav_file1 lav_file2
"""

	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	if (len(sys.argv) < 3):
		usage("you must specify two lav files")
	elif (len(sys.argv) > 3):
		usage("wrong number of arguments")

	lav1Filename = sys.argv[1]
	lav2Filename = sys.argv[2]

	# compare the files

	lav1 = file(lav1Filename,"rt")
	lav2 = file(lav2Filename,"rt")

	different = True
	stanza    = None
	lineNum   = 0

	while (True):
		lineNum += 1
		line1 = lav1.readline()
		line2 = lav2.readline()
		if (line1 == "") and (line2 == ""):
			different = False
			break
		line1 = line1.rstrip()
		line2 = line2.rstrip()

		if (stanza != None):
 			if (line1 == "}") != (line2 == "}"): break
 			if (line1 == "}") and (line2 == "}"):
	 			stanza = None
				continue
			stanzaIx += 1

		if (stanza == "d") and (stanzaIx == 1):
			continue	# ignore command line differences

		elif (stanza == "s") and (stanzaIx <= 2):
			line1 = line1.strip()
			line2 = line2.strip()

		elif (stanza == "h") and (stanzaIx <= 2):
			line1 = header_strip(line1)
			line2 = header_strip(line2)

		if (line1 != line2):
			# print >>sys.stderr,"%s\n%s" % (line1,line2)
			break

		if (stanza != None) and (line1 == "}"):
			stanza = None
			continue

		if (line1.endswith("{")):
			stanza  = line1[:-1].strip()
			stanzaIx = 0

	if (different):
		print >>sys.stderr,"FAILURE: %s and %s are different (line %d)" \
		                 % (lav1Filename,lav2Filename,lineNum)
		sys.exit(1)

	print >>sys.stderr,"SUCCESS: %s and %s are equivalent" \
					 % (lav1Filename,lav2Filename)


def header_strip(s):
	s = s.strip()
	if (s.startswith("\"")) and (s.endswith("\"")):
		s = s[1:-1].strip()
	if (s.startswith(">")):
		s = s[1:].strip()
	return s


if __name__ == "__main__": main()
