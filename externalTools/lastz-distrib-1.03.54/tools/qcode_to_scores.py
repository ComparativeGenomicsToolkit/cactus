#!/usr/bin/env python
"""
Convert quantum-code files to a LASTZ scores file
-------------------------------------------------

Given background probabilities, probabilities of each DNA substitution event,
and one (or two) quantum code files, we create a log-odds scoring matrix
suitable for LASTZ.

Typical command line:

	qcode_to_scores --scaleto=100 \
	   A:.26585  C:.23415  G:.23415  T:.26585 \   <--- background probabilties
	  AA:.18204 AC:.01903 AG:.04510 AT:.01967 \
	  CA:.01903 CC:.15508 CG:.01495 CT:.04510 \   <--- substitution probabilties
	  GA:.04510 GC:.01495 GG:.15508 GT:.01903 \
	  TA:.01967 TC:.04510 TG:.01903 TT:.18204 \
	  --code.target=<codefile> --code.query=<codefile>

An equivalent command line that takes advantage of the usual symmetry:

	qcode_to_scores --scaleto=100 \
	  --symmetric                             \
	   A:.26585  C:.23415                     \   <--- background probabilties
	  AA:.18204 AC:.01903 AG:.04510 AT:.01967 \   <--- substitution probabilties
	            CC:.15508 CG:.01495           \
	  --code.target=<codefile> --code.query=<codefile>

Quantum code files look something like the one below.  Each row represents a
quantum symbol.  The first value is the code value, either a single ascii
character or a two character hex value.  The remaining four values are the
probability of that symbol being A, C, G, or T.  Lines beginning with a # are
comments, and anything other than five columns is an error.

	#   p(A)      p(C)      p(G)      p(T)
	01  0.125041  0.080147  0.100723  0.694088
	02  0.111162  0.053299  0.025790  0.809749
	03  0.065313  0.007030  0.004978  0.922679
	 ...

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys
from math import log

def usage(s=None):
	message = """
qcode_to_scores [options] > lastz_score_file
  --scaleto=<max>         scale scores to give desired max
  --symmetric             map probabilities symmetrically
  --hoxd70                use HOXD70 (lastz default scores) for probabilities
  --code.target=<codefile> specify the quantum code for rows (LASTZ target)
  --code.query=<codefile>  specify the quantum code for columns (LASTZ query)
  --code=<codefile>       specify the quantum code for both rows *and* columns
  --creator=<string>      set name of creator to write as a comment in output
  --nocreator             inhibit creator comment in output
  <base>.target:<prob>    set target background probability of a nucleotide
  <base>.query:<prob>     set query background probability of a nucleotide
  <base>:<prob>           set background probability of a nucleotide for *both*
                          target and query 
  <basepair>:<prob>       set basepair substitution probability;  first base is
                          for target, second for query
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


bases          = ["A","C","G","T"]
basePairs      = ["AA","AC","AG","AT",
                  "CA","CC","CG","CT",
                  "GA","GC","GG","GT",
                  "TA","TC","TG","TT"]

baseSymmetries = [["A","T"],["C","G"]]
pairSymmetries = [["AA","TT"],["CC","GG"],["AT","TA"],["CG","GC"],
                  ["AC","CA","GT","TG"],["AG","GA","CT","TC"]]

hoxd70         = [("A", .26585),("C", .23415),
                  ("AA",.18204),("AC",.01903),("AG",.04510),("AT",.01967),
                                ("CC",.15508),("CG",.01495)]


def main():

	##########
	# parse the command line
	##########

	rProb       = {}
	cProb       = {}
	rcProb      = {}
	scaleTo     = None
	symmetric   = False
	dnaQuery    = True
	symbols     = []
	settings    = []
	rowCodeName = None
	colCodeName = None
	creator     = "qcode_to_scores"
	debug       = []

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
				assert (s not in rProb) and (s not in cProb), \
				       "duplicate DNA event: %s" % s
				rProb[s] = cProb[s] = p
		elif (arg in ["--code.row","--code.target"]) and (val != None):
			assert (rowCodeName == None), \
			       "can't have more than one row/target code"
			rowCodeName = val
		elif (arg in ["--code.column","--code.col","--code.query"]) and (val != None):
			assert (colCodeName == None), \
			       "can't have more than one column/target code"
			colCodeName = val
		elif (arg == "--code") and (val != None):
			assert (rowCodeName == None), \
			       "can't have more than one row/target code"
			assert (colCodeName == None), \
			       "can't have more than one column/target code"
			rowCodeName = colCodeName = val
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
			(s,which,p) = dna_event(arg)
			if   (which == "target"): w = "row"
			elif (which == "query"):  w = "col"
			elif (which == "column"): w = "col"
			else:                     w = which
			assert (w in ["row","col",None]), \
			       "can't decipher \"%s\" (in %s)" % (which,arg)
			if (w == "row"):
				assert (s in bases), \
				       "can't specify %s for %s (in %s)" % (which,s,arg)
				assert (s not in rProb), \
					   "duplicate DNA event: %s.target" % s
				rProb[s] = p
			elif (w == "col"):
				assert (s in bases), \
				       "can't specify %s for %s (in %s)" % (which,s,arg)
				assert (s not in cProb), \
					   "duplicate DNA event: %s.query" % s
				cProb[s] = p
			elif (s in bases):
				assert (s not in rProb) and (s not in cProb), \
				       "duplicate DNA event: %s" % s
				rProb[s] = cProb[s] = p
			else:
				assert (s not in rcProb), \
				       "duplicate DNA pair event: %s" % s
				rcProb[s] = p
		else:
			usage("unknown argument: %s" % arg)

	##########
	# sanity check
	##########

	if (symmetric):
		conProb = {}
		for nuc in bases:
			if (nuc in rProb) and (nuc not in cProb):
				conProb[nuc] = rProb[nuc]
			elif (nuc in cProb) and (nuc not in rProb):
				conProb[nuc] = cProb[nuc]
			elif (nuc in cProb) and (nuc in rProb):
				assert (rProb[nuc] == cProb[nuc]), \
				       "can't use --symmetric with %s.target != %s.query" \
				     % (nuc,nuc)
				conProb[nuc] = rProb[nuc]

		for group in baseSymmetries:
			present = len([x for x in group if (x in conProb)])
			assert (present == 1), \
			       "need a probability for exactly one of %s" \
			     % (",".join(group))
			val = None
			for x in group:
				if (x in conProb):
					val = conProb[x]
					break
			for x in group:
				if (x not in conProb): conProb[x] = val
		rProb = cProb = conProb

		for group in pairSymmetries:
			present = len([x for x in group if (x in rcProb)])
			assert (present == 1), \
			       "need a probability for exactly one of %s" \
			     % (",".join(group))
			val = None
			for x in group:
				if (x in rcProb):
					val = rcProb[x]
					break
			for x in group:
				if (x not in rcProb): rcProb[x] = val

	for nuc in bases:
		assert (nuc in rProb), \
		       "need a target probability for %s" % nuc
		assert (nuc in cProb), \
		       "need a query probability for %s" % nuc

	for xy in basePairs:
		assert (xy in rcProb), \
		       "need a probability for %s" % (xy)

	p = sum([rProb[nuc] for nuc in bases])
	assert (abs(p-1) < .00001), \
	       "target base probabilities sum to %f" % p

	p = sum([cProb[nuc] for nuc in bases])
	assert (abs(p-1) < .00001), \
	       "query base probabilities sum to %f" % p

	p = sum([rcProb[yx] for yx in basePairs])
	assert (abs(p-1) < .00001), \
	       "base pair probabilities sum to %f" % p

	##########
	# read code files
	##########

	# read row code

	if (rowCodeName == None):
		rowCode = simple_dna_quantum_code()
	else:
		rowCode = read_quantum_code(rowCodeName)

	if (".order" in rowCode):
		rowSymbols = rowCode[".order"]
	else:
		rowSymbols = [sym for sym in rowCode]
		rowSymbols.sort()

	# read column code

	if (colCodeName == None):
		colCode = simple_dna_quantum_code()
	elif (colCodeName == rowCodeName):
		colCode = rowCode
	else:
		colCode = read_quantum_code(colCodeName)

	if (".order" in colCode):
		colSymbols = colCode[".order"]
	else:
		colSymbols = [sym for sym in colCode]
		colSymbols.sort()

	##########
	# print what we got
	##########

	if ("debug" in debug):
		print "target" \
		    + "  ".join([" %s:%.5f" % (nuc,rProb[nuc]) for nuc in bases])
		print "query" \
		    + "  ".join([" %s:%.5f" % (nuc,cProb[nuc]) for nuc in bases])
		for y in bases:
			print "  ".join(["%s:%.5f" % (y+x,rcProb[y+x]) for x in bases])

	##########
	# assign scores
	##########

	sub = {}
	maxSub = None

	for row in rowSymbols:
		u = rowCode[row]
		sub[row] = {}
		for col in colSymbols:
			v = colCode[col]
			numer = sum([u[y]*v[x]*rcProb[y+x]       for (y,x) in basePairs])
			denom = sum([u[y]*v[x]*rProb[y]*cProb[x] for (y,x) in basePairs])
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

	wRow = max([len(row) for row in rowSymbols])

	if (scaleTo != None) and (type(scaleTo) == int):
		wCol = 4
		for row in rowSymbols:
			for col in colSymbols:
				wCol = max(wCol,len("%d" % sub[row][col]))

		print "%-*s %s" \
		    % (wRow," ","  ".join(["%*s" % (wCol,col) for col in colSymbols]))

		for row in rowSymbols:
			print "%-*s %s" \
			    % (wRow,row,
			       "  ".join(["%*d" % (wCol,sub[row][col]) for col in colSymbols]))

	else:
		wCol = 4
		for row in rowSymbols:
			for col in colSymbols:
				wCol = max(wCol,len("%.6f" % sub[row][col]))

		print "%-*s %s" \
		    % (wRow," ","  ".join(["%*s" % (wCol,col) for col in colSymbols]))

		for row in rowSymbols:
			print "%-*s %s" \
			    % (wRow,row,
			       "  ".join(["%*.6f" % (wCol,sub[row][col]) for col in colSymbols]))


def simple_dna_quantum_code():
	symToProfile = {}
	for nuc1 in bases:
		symToProfile[nuc1] = {}
		for nuc2 in bases:
			if (nuc2 == nuc1): symToProfile[nuc1][nuc2] = 1
			else:              symToProfile[nuc1][nuc2] = 0
	return symToProfile


def read_quantum_code(codeName):
	codeF = file (codeName, "rt")

	symToProfile = {}
	codeNumUsed  = {}
	symOrder     = []

	lineNum = 0
	for line in codeF:
		lineNum += 1
		line = line.strip()
		if ("#" in line):
			line = line.split("#",1)[0].strip()
		if (line == ""):
			continue

		fields = line.split()

		assert (len(fields) >= 5), \
		       "fewer than four probabilities (%s line %d)" \
		     % (codeName,lineNum)
		assert (len(fields) <= 5), \
		       "more than four probabilities (%s line %d)" \
		     % (codeName,lineNum)

		try:
			sym     = fields[0]
			codeNum = quantum_code_num(sym)
		except ValueError:
			assert (False), \
			       "%s is not a valid quantum symbol (%s line %d)" \
			     % (sym,codeName,lineNum)

		if (codeNum in codeNumUsed):
			assert (False), \
			       "%s (or equivalent) appears more than once (%s line %d)" \
			     % (sym,codeName,lineNum)

		try:
			profile = {}
			for ix in range(4):
				p = float_or_fraction(fields[ix+1])
				if (not (0 <= p <= 1)): raise ValueError
				profile[bases[ix]] = p
		except:
			assert (False), \
			       "%s is a bad probability value (%s line %d)" \
			     % (fields[ix+1],codeName,lineNum)

		symToProfile[sym] = profile
		codeNumUsed[codeNum] = True
		symOrder += [sym]

	codeF.close ()

	# sanity check

	assert (len(symToProfile) >= 1), \
	       "%s contains no code vectors!" % codeName

	for sym in symToProfile:
		p = sum([symToProfile[sym][nuc] for nuc in bases])
		assert (abs(p-1) < .00001), \
		       "probabilities for %s sum to %f (in %s)" % (sym,p,codeName)

	symToProfile[".order"] = symOrder

	return symToProfile


def dna_event(s):
	(s,p) = s.split(":",1)
	if ("." in s): (s,which) = s.split(".",1)
	else:          which     = None
	assert (valid_dna_event(s)), "invalid DNA event: %s" % s
	try:
		p = float_or_fraction(p)
		if (not (0 <= p <= 1)): raise ValueError
	except ValueError:
		assert (False), "invalid probability for %s: %s" % (s,p)
	return (s,which,p)


def valid_dna_event(s):
	if (len(s) == 0):
		return False
	if (len(s) == 1):
		return (s in bases)
	if (len(s) == 2):
		return (s[0] in bases) and (s[1] in bases)
	return False


def float_or_fraction(s):
	if ("/" in s):
		(n,d) = s.split("/",1)
		return float(n)/float(d)
	else:
		return float(s)


def quantum_code_num(s):
	if (len(s) == 0):
		raise ValueError
	if (len(s) == 1):
		if (0x21 <= ord(s) <= 0x7E): return ord(s)
		else:                        raise ValueError
	if (len(s) == 2):
		if (s == "00"): raise ValueError
		try:    return int(s,16)
		except: raise ValueError
	raise ValueError


if __name__ == "__main__": main()
