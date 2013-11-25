#!/usr/bin/env python
"""
Convert probabilities to a LASTZ scores file (including quantum scores)
-----------------------------------------------------------------------

Given background probabilities, probabilities of each DNA substitution event,
and an optional list of quantum symbols, we create a log-odds scoring matrix
suitable for LASTZ.

Typical command line:

	probabilities_to_scores --scaleto=100 \
	   A:.26585  C:.23415  G:.23415  T:.26585 \   <--- background probabilties
	  AA:.18204 AC:.01903 AG:.04510 AT:.01967 \
	  CA:.01903 CC:.15508 CG:.01495 CT:.04510 \   <--- substitution probabilties
	  GA:.04510 GC:.01495 GG:.15508 GT:.01903 \
	  TA:.01967 TC:.04510 TG:.01903 TT:.18204 \
	  R=G:.5,A:.5 Y=T:.5,C:.5                     <--- quantum symbols

An equivalent command line that takes advantage of the usual symmetry:

	probabilities_to_scores --scaleto=100 \
	  --symmetric                             \
	   A:.26585  C:.23415                     \   <--- background probabilties
	  AA:.18204 AC:.01903 AG:.04510 AT:.01967 \   <--- substitution probabilties
	            CC:.15508 CG:.01495           \
	  R=G:.5,A:.5 Y=T:.5,C:.5                     <--- quantum symbols

The resulting scores file would look like this:

	     A     C     G     T     R     Y
	A   91  -114   -31  -123    52  -119
	C -114   100  -125   -31  -119    52
	G  -31  -125   100  -114    52  -119
	T -123   -31  -114    91  -119    52
	R   52  -119    52  -119    52  -119
	Y -119    52  -119    52  -119    52

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys
from math import log

def usage(s=None):
	message = """
probabilities_to_scores [options] > lastz_score_file
  --scaleto=<max>    scale scores to give desired max
  --symmetric        map probabilities symmetrically
  --nodna            don't include A,G,C,T in the alphabets
  --dnarows          (target) row alphabet is A,C,G,T
  --dnacol[umn]s     (query) column alphabet is A,C,G,T
  --hoxd70           use HOXD70 (lastz default scores) for probabilities
  --iupac            alphabets are IUPAC 15-letter code
  --writecode=<file> write quantum code to a file
  --creator=<string> set name of creator to write as a comment in output
  --nocreator        inhibit creator comment in output
  <base>=<prob>      set background probability of a nucleotide
  <basepair>=<prob>  set basepair substitution probability
  <symbol>=<profile> define the profile for a quantum symbol
                     .. e.g. Y=T:.5,C:.5 or 07=A:0.311,C:0.228,G:0.422,T:0.039
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


bases      = "ACGT"
basePairs  = ["AA","AC","AG","AT",
              "CA","CC","CG","CT",
              "GA","GC","GG","GT",
              "TA","TC","TG","TT"]

symmetries = [["A","T"],["C","G"],
              ["AA","TT"],["CC","GG"],["AT","TA"],["CG","GC"],
              ["AC","CA","GT","TG"],["AG","GA","CT","TC"]]

hoxd70     = [("A", .26585),("C", .23415),
              ("AA",.18204),("AC",.01903),("AG",.04510),("AT",.01967),
                            ("CC",.15508),("CG",.01495)]

iupac      = [("R","G,A"),
              ("Y","T,C"),
              ("K","G,T"),
              ("M","A,C"),
              ("S","G,C"),
              ("W","A,T"),
              ("B","G,T,C"),
              ("D","G,A,T"),
              ("H","A,C,T"),
              ("V","G,C,A"),
              ("N","A,C,G,T")]


def main():

	##########
	# parse the command line
	##########

	prob       = {}
	scaleTo    = None
	symmetric  = False
	dnaQuery   = True
	symbols    = []
	symProb    = {}
	symGroup   = {}
	settings   = []
	rowsAreDNA = False
	colsAreDNA = False
	creator    = "probabilities_to_scores"
	codeName   = None
	debug      = []

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

		if (arg == "--scaleto") and (val != None):
			try:               scaleTo = int(val)
			except ValueError: scaleTo = float(val)
		elif (arg == "--symmetric") and (val == None):
			symmetric = True
		elif (arg == "--nodna") and (val == None):
			dnaQuery = False
		elif (arg == "--dnarows") and (val == None):
			rowsAreDNA = True
		elif (arg in ["--dnacols","--dnacolumns"]) and (val == None):
			colsAreDNA = True
		elif (arg in ["--hoxd70","--HOXD70"]) and (val == None):
 			symmetric = True
			for (s,p) in hoxd70:
				assert (s not in prob), "duplicate DNA event: %s" % s
				prob[s] = p
		elif (arg in ["--iupac","--IUPAC"]) and (val == None):
			for (sym,val) in iupac:
				assert (sym not in symProb), "duplicate quantum symbol: %s" % sym
				symbols += [sym]
				symProb[sym]  = {}
				symGroup[sym] = ""
				vals = val.split(",")
				for s in vals:
					symProb[sym][s] = 1.0/len(vals)
					symGroup[sym] += s
		elif (arg == "--writecode") and (val != None):
			codeName = val
		elif (arg == "--nocreator") and (val == None):
			creator = None
		elif (arg == "--creator") and (val != None):
			creator = val
		elif (arg == "--debug") and (val != None):
			debug.append(val)
		elif (arg == "--debug") and (val == None):
			debug.append("debug")
		elif (arg.startswith("--")) and (val != None):
			settings += [(arg[2:],val)]
		elif (arg.startswith("--")):
			usage("unknown argument: %s" % arg)
		elif (val == None) and (":" in arg):
			(s,p) = dna_event(arg)
			assert (s not in prob), "duplicate DNA event: %s" % s
			prob[s] = p
		elif (valid_quantum_symbol(arg)) and (val != None):
			sym = arg
			assert (sym not in symProb), "duplicate quantum symbol: %s" % sym
			symbols += [sym]
			symProb[sym]  = {}
			symGroup[sym] = ""
			vals = val.split(",")
			haveProbs = False
			for val in vals:
				if (":" in val):
					haveProbs = True
					break
			if (haveProbs):
				for val in vals:
					(s,p) = dna_event(val)
					assert (len(s) == 1), \
						   "invalid DNA event for %s: %s" % (sym,s)
					assert (s not in symProb[sym]), \
						   "duplicate DNA event for %s: %s" % (sym,s)
					symProb[sym][s] = p
					symGroup[sym] += s
			else:
				for s in vals:
					assert (len(s) == 1) and (s in bases), \
						   "invalid DNA event for %s: %s" % (sym,s)
					assert (s not in symProb[sym]), \
						   "duplicate DNA event for %s: %s" % (sym,s)
					symProb[sym][s] = 1.0/len(vals)
					symGroup[sym] += s
		else:
			usage("unknown argument: %s" % arg)

	##########
	# sanity check
	##########

	if (symmetric):
		for group in symmetries:
			present = len([x for x in group if (x in prob)])
			assert (present == 1), \
			       "need a probability for exactly one of %s" \
			     % (",".join(group))
			val = None
			for x in group:
				if (x in prob):
					val = prob[x]
					break
			for x in group:
				if (x not in prob): prob[x] = val

	for nuc in bases:
		assert (nuc in prob), \
		       "need a probability for %s" % nuc

	for xy in basePairs:
		assert (xy in prob), \
		       "need a probability for %s" % (xy)

	p = sum([prob[nuc] for nuc in bases])
	assert (abs(p-1) < .000001), \
	       "base probabilities sum to %f" % p

	p = sum([prob[xy] for xy in basePairs])
	assert (abs(p-1) < .000001), \
	       "base pair probabilities sum to %f" % p

	for sym in symProb:
		p = sum([symProb[sym][nuc] for nuc in symProb[sym]])
		assert (abs(p-1) < .000001), \
		       "probabilities for %s sum to %f" % (sym,p)
		for nuc in bases:
			if (nuc not in symProb[sym]):
				symProb[sym][nuc] = 0

	if (dnaQuery):
		for sym in bases:
			if (sym in symProb): continue
			symbols += [sym]
			symProb[sym]  = {}
			symGroup[sym] = sym
			for nuc in bases:
				if (nuc == sym): symProb[sym][nuc] = 1
				else:            symProb[sym][nuc] = 0
		symbols = [sym for sym in bases] \
		        + [sym for sym in symbols if (sym not in bases)]

	if (rowsAreDNA): rowSymbols = bases
	else:            rowSymbols = symbols

	if (colsAreDNA): colSymbols = bases
	else:            colSymbols = symbols

	##########
	# print what we got
	##########

	if ("debug" in debug):
		print "  ".join([" %s:%.5f" % (nuc,prob[nuc]) for nuc in bases])

		for x in bases:
			print "  ".join(["%s:%.5f" % (x+y,prob[x+y]) for y in bases])

		print
		for sym in symbols:
			p = symProb[sym]
			print "%s -> %s" \
			    % (sym,"  ".join([" %s:%.5f" % (nuc,p[nuc]) for nuc in bases]))

	##########
	# write quantum code file
	##########
	
	if (codeName != None):
		codeF = file(codeName,"wt")
		for sym in symbols:
			p = symProb[sym]
			print >>codeF, "%s\t%s" \
			             % (sym,"\t".join(["%.6f" % p[nuc] for nuc in bases]))
		codeF.close()

	##########
	# assign scores
	##########

	sub = {}
	maxSub = None

	for row in rowSymbols:
		u = symProb[row]
		sub[row] = {}
		for col in colSymbols:
			v = symProb[col]
			numer = sum([u[y]*v[x]*prob[y+x]       for (y,x) in basePairs])
			denom = sum([u[y]*v[x]*prob[y]*prob[x] for (y,x) in basePairs])
			sub[row][col] = log (float(numer) / float(denom))
			if (maxSub == None) or (sub[row][col] > maxSub):
				maxSub = sub[row][col]

	if (scaleTo != None):
		scale = scaleTo / maxSub
		for row in rowSymbols:
			for col in colSymbols:
				sub[row][col] *= scale
				if (type(scaleTo) == int):
					sub[row][col] = round(sub[row][col])

	##########
	# print the settings, if there are any
	##########

	if (creator != None):
		print "# created by %s" % creator
		print

	if (settings != []):
		sLen = max([len(s) for (s,val) in settings])
		for (s,val) in settings:
			print "%-*s = %s" % (sLen,s,val)
		print

	##########
	# print the substitution matrix
	##########

	if (scaleTo != None) and (type(scaleTo) == int):
		wSub = 4
		for row in rowSymbols:
			for col in colSymbols:
				wSub = max(wSub,len("%d" % sub[row][col]))

		print "%s %s" \
		    % ("#","  ".join(["%*s" % (wSub,non_single(symGroup[col])) for col in colSymbols]))

		print "%s %s" \
		    % (" ","  ".join(["%*s" % (wSub,col) for col in colSymbols]))

		for row in rowSymbols:
			print "%s %s%s" \
			    % (row,
			       "  ".join(["%*d" % (wSub,sub[row][col]) for col in colSymbols]),
			       non_single_comment(symGroup[row]))

	else:
		wSub = 4
		for row in rowSymbols:
			for col in colSymbols:
				wSub = max(wSub,len("%.6f" % sub[row][col]))

		print "%s %s" \
		    % ("#","  ".join(["%*s" % (wSub,non_single(symGroup[col])) for col in colSymbols]))

		print "%s %s" \
		    % (" ","  ".join(["%*s" % (wSub,col) for col in colSymbols]))

		for row in rowSymbols:
			print "%s %s%s" \
			    % (row,
			       "  ".join(["%*.6f" % (wSub,sub[row][col]) for col in colSymbols]),
			       non_single_comment(symGroup[row]))


def dna_event(s):
	(s,p) = s.split(":",1)
	assert (valid_dna_event(s)), "invalid DNA event: %s" % s
	try:
		p = float(p)
		if (not (0 <= p <= 1)): raise ValueError
	except ValueError:
		assert (False), "invalid probability for %s: %s" % (s,p)
	return (s,p)


def valid_dna_event(s):
	if (len(s) == 0):
		return False
	if (len(s) == 1):
		return (s in bases)
	if (len(s) == 2):
		return (s[0] in bases) and (s[1] in bases)
	return False


def valid_quantum_symbol(s):
	if (len(s) == 0):
		return False
	if (len(s) == 1):
		return (s in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
	if (len(s) == 2):
		if (s == "00"): return False
		return (s[0] in "0123456789ABCDEF") and (s[1] in "0123456789ABCDEF")
	return False


def non_single_comment(s):
	if (len(s) == 1): return ""
	else:             return " # " + s


def non_single(s):
	if (len(s) == 1): return ""
	else:             return s


if __name__ == "__main__": main()
