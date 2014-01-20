#!/usr/bin/env python
"""
Sort the a-stanzas in a lav file, according to the user's choice of key
-----------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys

validKeys = ["score","pos1","pos2","beg1","beg2","end1","end2"]

def usage(s=None):
	message = """
lav_sort --key=[-]<score|beg1|beg2|end1|end2> < lav_file > lav_file
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

	# process the stanzas

	blocks = []
	for (kind,stanza) in read_stanzas(sys.stdin):
		if (kind == "a"):
			key = get_key_value(keyName,stanza)
			blocks += [(key,stanza)]
			continue
		if (len(blocks) > 0):
			blocks.sort()
			if (keyReverse): blocks.reverse()
			for (key,s) in blocks:
				print "\n".join(s)
			blocks = []
		print "\n".join(stanza)

	if (len(blocks) > 0):
		blocks.sort()
		if (keyReverse): blocks.reverse()
		for (key,s) in blocks:
			print "\n".join(s)

# read_stanzas--
#	Collect the lines that belong to the next stanza.  A stanza has the form
#	shown below.  It consists of several lines bracketed by a pair of curlies,
#	and has a type indicated by a single letter.
#
#		x {
#		  ...
#		}
#
#	In this routine we generalize the stanza concept to include lines not
#	strictly with a pair of curlies.  First, lines beginning with a "#:" are
#	considered to be single line stanzas with no type (e.g. the "#:lav" and
#	"#:eof" lines).  Second, any other blank lines are appended to whatever
#	stanza preceeded them.  This allows for lav+text and other debugging output
#	from lastz to be carried around with the appropriate stanza.

def read_stanzas(f):
	kind    = None
	stanza  = []
	inCurly = False
	for line in f:
		line = line.rstrip()
		if (not inCurly):
			isWaffle = line.startswith("#:")
			inCurly  = (len(line) == 3) and (line.endswith(" {"))
			if (isWaffle) or (inCurly):
				if (len(stanza) > 0):
					yield (kind,stanza)
					stanza = []
				if (isWaffle):
					yield (line[2:],[line])
					kind = None
					continue
				kind = line[0]
			stanza += [line]
		else: # (inCurly)
			stanza += [line]
			if (line == "}"): inCurly = False

	assert (len(stanza) == 0), "premature end of file"

# get_key_value--
#	Extract the specied key value from an a-stanza.  A typical a-stanza looks
#	like this one:
#
#		a {
#		  s 14400
#		  b 425 4438
#		  e 697 4714
#		  l 425 4438 448 4461 96
#		  l 449 4464 579 4594 83
#		  l 581 4595 604 4618 96
#		  l 605 4627 609 4631 100
#		  l 617 4632 648 4663 91
#		  l 649 4666 697 4714 90
#		}

def get_key_value(keyName,aStanza):
	if (keyName == "score"):
		assert (len(aStanza) >= 2) and (aStanza[1].startswith("  s"))
		score = aStanza[1].split()[1]
		try:
			return int(score)
		except:
			try:
				return float(score)
			except:
				pass
		return score

	if (keyName in ["pos1","beg1"]):
		assert (len(aStanza) >= 3) and (aStanza[2].startswith("  b"))
		beg1 = aStanza[2].split()[1]
		return int(beg1)

	if (keyName in ["pos2","beg2"]):
		assert (len(aStanza) >= 3) and (aStanza[2].startswith("  b"))
		beg2 = aStanza[2].split()[2]
		return int(beg2)

	if (keyName in ["end1"]):
		assert (len(aStanza) >= 4) and (aStanza[3].startswith("  e"))
		end1 = aStanza[3].split()[1]
		return int(end1)

	if (keyName in ["end2"]):
		assert (len(aStanza) >= 4) and (aStanza[3].startswith("  e"))
		end2 = aStanza[3].split()[2]
		return int(end2)

	assert False


if __name__ == "__main__": main()
