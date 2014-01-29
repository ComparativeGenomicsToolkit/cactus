#!/usr/bin/env python
"""
Add scoring-related parameters to a lastz scores file
-----------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)

Typical input scores file:

	# (a LASTZ scoring set, created by "LASTZ --infer")

	bad_score          = X:-1910 # used for sub[X][*] and sub[*][X]
	fill_score         = -191    # used when sub[*][*] not otherwise defined
	gap_open_penalty   = 400
	gap_extend_penalty = 30

	      A     C     G     T
	A    85  -164   -70  -191
	C  -164   100  -151   -70
	G   -70  -151   100  -164
	T  -191   -70  -164    85
"""

import sys

def usage(s=None):
	message = """
expand_scores_file [options]< scores_file > scores_file
  --overridegaps  ignore gap scores already set 
"""

	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	overrideGaps = False

	for arg in sys.argv[1:]:
		if (arg == "--overridegaps"):
				overrideGaps = True
				continue
		raise "unrecognized argument: %s" % arg

	# read the scores file

	lines = []
	numValueLines = None
	valuesFinished = False
	nameToVal = {}
	subs = subRows = subColumns = None

	lineNumber = 0
	for line in sys.stdin:
		lineNumber += 1
		line = line.rstrip()
		lines += [line]
		if (line == ""): continue
		if (line.startswith("#")): continue
		if ("#" in line): line = line.split("#",1)[0].strip()

		if ("=" in line):
			if (valuesFinished):
				raise "in scores file, unexpected assignment (line %d): %s" \
					% (lineNumber,line)
			fields = line.split("=",1)
			name = fields[0].strip()
			val  = fields[1].strip()
			if   (name == "gap_open_penalty"):   name = "O"
			elif (name == "gap_extend_penalty"): name = "E"
			if (name in nameToVal):
				raise "in scores file, %s is assigned twice (line %d): %s" \
					% (name,lineNumber,line)
			if (overrideGaps):
				if (name in ["O","E"]):
					lines.pop()
					continue
			try:
				nameToVal[name] = int_or_float(fields[1])
			except:
				if (name in ["O","E"]):
					raise "in scores file, bad assignment value (line %d): %s" \
						% (lineNumber,line)
		elif (not valuesFinished):
			numValueLines = len(lines) - 1
			valuesFinished = True
			subColumns = line.split()
			subRows    = []
			subs       = {}
		else:
			fields  =  line.split()
			rowCh   =  fields.pop(0)
			subRows += [rowCh]
			if (len(fields) != len(subColumns)):
				raise "in scores file, inconsistent matrix (line %d): %s" \
					% (lineNumber,line)
			for ix in range(len(fields)):
				colCh = subColumns[ix]
				subs[rowCh+colCh] = int_or_float(fields[ix])

	if (subs == None):
		raise "scores file is missing a matrix"

	if ("AA" not in subs):
		raise "scores file lacks A-to-A score"

	# compute a few values from the scores matrix

	bestSub  = float(max([subs[digram] for digram in subs]))
	worstSub = float(min([subs[digram] for digram in subs]))
	aaSub    = float(subs["AA"])

	# add expanded values

	knownVals = [name for name in nameToVal]

	if ("O" not in nameToVal):
		nameToVal["O"] = -int(3.25 * worstSub)
		
	if ("E" not in nameToVal):
		nameToVal["E"] = -int(0.25 * worstSub)
		
	if ("X" not in nameToVal):
		nameToVal["X"] = int(10 * aaSub)
		
	if ("Y" not in nameToVal):
		nameToVal["Y"] = int(nameToVal["O"] + 100*nameToVal["E"])
		
	if ("K" not in nameToVal):
		nameToVal["K"] = int(30 * bestSub)
		
	if ("L" not in nameToVal):
		nameToVal["L"] = int(30 * bestSub)

	if ("T" not in nameToVal) and (worstSub/bestSub < -1.5):
		nameToVal["T"] = "2"

	if ("Z" not in nameToVal) and (worstSub/bestSub < -3.0):
		nameToVal["Z"] = "3"

	# figure out what values we've added, and in what order to print them

	addedNames =  [name for name in ["T","Z","O","E","X","Y","K","L"] \
	                 if    (name in nameToVal) \
	                   and (name not in knownVals)]
	addedNames += [name for name in nameToVal \
	                 if    (name not in addedNames) \
	                   and (name not in knownVals)]

	# print the new scores file

	blankLine = False

	for ix in range(numValueLines):
		print lines[ix]
		blankLine = (lines[ix] == "")

	if (addedNames != []):
		if (not blankLine): print ""
		print "# (score parameters added by expand_scores_file)"
		print ""

		for name in addedNames:
			print "%s=%s" % (name,nameToVal[name])

		blankLine = (lines[numValueLines] == "")
		if (not blankLine): print ""

	for ix in range(numValueLines,len(lines)):
		print lines[ix]


def int_or_float(s):
	try:    return int(s)
	except: return float(s)


if __name__ == "__main__": main()
