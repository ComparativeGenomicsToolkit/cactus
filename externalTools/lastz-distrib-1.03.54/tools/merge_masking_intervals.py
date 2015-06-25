#!/usr/bin/env python
"""
Given a file of masking intervals, combine overlapping intervals.

Masking intervals are as would be used as lastz's softmask, xmask, or nmask
sequence specifier actions.  They can be produced by lastz using the --masking
and --outputmasking options.

	input:            output:
	555110 555310     555110 555310
	555941 556479     555941 556663
	555966 556402
	555976 556402
	555977 556402
	556125 556479
	556153 556663
	557674 558206     557674 558278
	557802 558278
	559509 559769     559509 559769
	798462 798922     798462 799008
	798462 798963
	798614 799008
	799495 799603     799495 799603

Intervals are origin 1, closed.  They needn't be sorted.

We consider adjoining intervals to be overlapping.
"""

__author__ = "Bob Harris (rsharris@bx.psu.edu)"


from sys import argv,stdin


def main():
	global origin,adjoining

	assert (len(argv) == 1), "give me no arguments"

	# collect the intervals
	# nota bene: internally we work with them as origin-zero, half-open

	intervals = []

	lineNumber = 0
	for line in stdin:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		try:
			s = int(fields[0]) - 1
			e = int(fields[1])
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		intervals += [(s,e)]

	# merge 'em

	intervals.sort()

	start = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s > end):
			print "%d\t%d" % (start+1,end)
			(start,end) = (s,e)
			continue
		elif (e > end):
			end = e

	if (start != None):
		print "%d\t%d" % (start+1,end)


if __name__ == "__main__": main()
