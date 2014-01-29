#!/usr/bin/env python
"""
Sort alignment blocks in a maf file, according to the user's choice of key
--------------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys,re

validKeys = ["score","pos1","pos2","beg1","beg2","end1","end2","diag","name1","name2"]

def usage(s=None):
	message = """
maf_sort --key=[-]<score|beg1|beg2|end1|end2|diag|name1|name2> < maf_file > maf_file
"""

	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	if (len(sys.argv) < 2):
		usage("you must specify a key")
	elif (len(sys.argv) > 2):
		usage("wrong number of arguments")

	arg = sys.argv[1]
	if (not arg.startswith("--key=")):
		usage("unrecognized argument: \"%s\"" % arg)

	keyName    = arg[arg.find("=")+1:]
	keyReverse = False
	if (keyName.startswith("-")):
		keyName    = keyName[1:]
		keyReverse = True
	if (keyName.startswith("+")):
		keyName    = keyName[1:]
		keyReverse = False
	if (keyName not in validKeys):
		usage("unrecognized key: \"%s\"" % keyName)

	# process the blocks

	blocks = []
	for (block,comments) in read_blocks(sys.stdin):
		key = get_key_value(keyName,block)
		blocks += [(key,block,comments)]

	if (len(blocks) > 0):
		blocks.sort()
		if (keyReverse): blocks.reverse()
		for (key,block,comments) in blocks:
			if (comments != []):
				print "\n".join([line for line in comments])
			print "\n".join([line for line in block])
			print


# read_blocks--
#	Collect the lines that belong to the next alignment block.  A block has the
#	form shown below.
#
#		a score=19951
#		s apple  23871 367 + 70000 CCCCCGC...
#		s orange    13 390 -   408 CTCCTGC...

def read_blocks(f):
	comments = []
	block    = []
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.rstrip()
		if (line.startswith("#")):
			comments += [line]
			continue
		if (line == ""):
			if (len(block) == 3):
				yield (block,comments)
				comments = []
				block    = []
				continue
			elif (len(block) == 0):
				continue
			else:
				assert (False), "premature end of block at line %d" % lineNumber
		if (len(block) == 3): "long block at line %d" % lineNumber
		block += [line]

	if (len(block) == 3):
		yield (block,comments)
	elif (len(block) != 0):
		assert (False), "premature end of file"


# get_key_value--
#	Extract the specied key value from a maf block
#
#		a score=19951
#		s apple  23871 367 + 70000 CCCCCGC...
#		s orange    13 390 -   408 CTCCTGC...

scoreRe = re.compile("^a score=(?P<score>.+)$")
textRe  = re.compile("^s"
                   + " +(?P<name>[^ ]+)"
                   + " +(?P<pos>[0-9]+)"
                   + " +(?P<len>[0-9]+)"
                   + " +(?P<strand>[-+])"
                   + " +[0-9]+"
                   + " +[-ACGTacgtNn]+$")


def get_key_value(keyName,block):
	try:
		line = block[0]
		m = scoreRe.match(line)
		if (m == None): raise ValueError
		score = float(m.group("score"))
	except ValueError:
		assert (False), "bad score line: %s" % line

	try:
		line = block[1]
		m = textRe.match(line)
		if (m == None): raise ValueError
		name1   = m.group("name")
		pos1    = int(m.group("pos"))
		len1    = int(m.group("len"))
		strand1 = m.group("strand")
	except ValueError:
		assert (False), "bad line: %s" % line

	try:
		line = block[2]
		m = textRe.match(line)
		if (m == None): raise ValueError
		name2   = m.group("name")
		pos2    = int(m.group("pos"))
		len2    = int(m.group("len"))
		strand2 = m.group("strand")
	except ValueError:
		assert (False), "bad line: %s" % line

	if (keyName == "score"):
		return (score,pos1,strand1,pos2,strand2,len1,len2,name1,name2)

	if (keyName in ["pos1","beg1"]):
		return (pos1,strand1,pos2,strand2,len1,len2,score,name1,name2)

	if (keyName in ["pos2","beg2"]):
		return (pos2,strand2,pos1,strand1,len2,len1,score,name1,name2)

	if (keyName in ["end1"]):
		return (pos1+len1,strand1,pos2+len2,strand2,len1,len2,score,name1,name2)

	if (keyName in ["end2"]):
		return (pos2+len2,strand2,pos1+len1,strand1,len2,len1,score,name1,name2)

	if (keyName in ["diag"]):
		return (strand1,strand2,pos1-pos2,pos1,len1,len2,score,name1,name2)

	if (keyName in ["name1"]):
		return (name1,score,len1,strand1,pos1,name2,len2,strand2,pos2)

	if (keyName in ["name2"]):
		return (name2,score,len2,strand2,pos2,name1,len1,strand1,pos1)

	assert False


if __name__ == "__main__": main()
