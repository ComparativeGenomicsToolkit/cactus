#!/usr/bin/env python
"""
Break a fasta file into fragments.

$$$ todo: spread out the fragment starts so that the last fragment ends at the
$$$       .. end of a sequence, if possible

$$$ todo: find runs of N and reset the fragment start position to skip past
$$$       .. such runs
"""

from sys import argv,stdin,exit
import math

def usage(s=None):
	message = """fasta_fragments [options] < fasta_file > fasta_file
  Split a fasta file into overlapping fragments.

  options:
    --fragment=<length>       length of each fragment
                              (default is 100)
    --step=<length>           distance between the start of each fragment
                              (default is 50)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	fragmentLength = 100
	stepLength     = 50

	for arg in argv[1:]:
		if (arg.startswith("--fragment=")):
			fragmentLength = int(arg.split("=",1)[1])
		elif (arg.startswith("--step=")):
			stepLength = int(arg.split("=",1)[1])
		elif (arg.startswith("--")):
			usage("can't understand %s" % arg)
		else:
			usage("can't understand %s" % arg)

	# process the sequences

	for (name,seq) in chunk_fasta_sequences(stdin, stepLength, fragmentLength):
		print ">%s\n%s\n" % (name,seq)

# fasta_sequences--
#	Read the fasta sequences from a file	

def chunk_fasta_sequences(f, stepLength, fragmentLength):
	"""Chunks up a fasta stream into overlapping fragments.
	
	>>> for (seqName, seq) in chunk_fasta_sequences([">1", "ACTG", "AGGG", "TGCTGC", ">2", "AT", ">3", "CCCGCCT", ">4" ], 3, 6):
	...     print seqName, seq                                                                                                       
	1_1 ACTGAG
	1_4 GAGGGT
	1_7 GGTGCT
	1_10 GCTGC
	2_1 AT
	3_1 CCCGCC
	3_4 GCCT
	>>> 
	"""
	allN = "N" * fragmentLength
	seqName = None
	seqNucs = ""
	offset = 0
	
	def chunk(eatWholeSequence):
		for ix in xrange(0,len(seqNucs),stepLength):
			frag = seqNucs[ix:ix+fragmentLength]
			if frag == allN: 
				continue
			if len(frag) < fragmentLength:
				if eatWholeSequence and len(frag) > 0:
					yield ("%s_%d" % (seqName,ix+1+offset), frag)
				break
			yield ("%s_%d" % (seqName,ix+1+offset), frag)
			
	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				for seqPair in chunk(True):
					yield seqPair
			seqName = line[1:].strip()
			seqNucs = ""
			offset = 0
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += line
			for seqPair in chunk(False):
				yield seqPair
			#Cut off the trailing sequence
			if len(seqNucs) >= fragmentLength:
				i = int(math.ceil((float(len(seqNucs)) - fragmentLength)/stepLength))*stepLength
				offset += i
				seqNucs = seqNucs[i:]

	if (seqName != None):
		for seqPair in chunk(True):
			yield seqPair

if __name__ == "__main__": 
	import doctest
	doctest.testmod()
	main()
