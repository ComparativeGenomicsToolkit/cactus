#!/usr/bin/env python
"""
Break a fasta file into fragments.

$$$ todo: spread out the fragment starts so that the last fragment ends at the
$$$       .. end of a sequence, if possible

$$$ todo: find runs of N and reset the fragment start position to skip past
$$$       .. such runs
"""

from sys import argv,stdin,exit


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

	allN = "N" * fragmentLength

	# process the sequences

	for (name,seq) in fasta_sequences(stdin):
		seq = seq.upper()
		for ix in xrange(0,len(seq)-fragmentLength,stepLength):
			frag = seq[ix:ix+fragmentLength]
			if (frag == allN): continue
			print ">%s_%d" % (name,ix+1)
			print frag


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = line[1:].strip()
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


if __name__ == "__main__": main()

