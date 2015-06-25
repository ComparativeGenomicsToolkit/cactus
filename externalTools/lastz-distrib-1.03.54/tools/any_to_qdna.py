#!/usr/bin/env python
"""
Convert any file to a LASTZ quantum dna file, just by appending qdna headers

Qdna file format is shown below (omitting "named properties", which we don't
use).  We simply create all the headers and copy the file as the "data
sequence".
	
	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
	offset 0x04: 00 00 02 00   version 2.0 (fourth byte is sub version)
	offset 0x08: 00 00 00 14   header length (in bytes, including this field)
	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
	offset 0x18: 00 00 00 00   (offset to named properties, not used)
	offset    N: ...           name (zero-terminated string)
	offset    S: ...           data sequence

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

from sys import argv,stdin,stdout,exit


def usage(s=None):
	message = """any_to_qdna [options] < any_file > qdna_file
  Convert any file to a LASTZ quantum dna file.

  options:
    --name=<string>    the name of the sequence
                       (by default, the sequence is unnamed)
    --striplinebreaks  strip line breaks from the file
                       (default is to include line breaks in the qdna file)
    --simple           create an "old-style" qdna file
                       (default is to create a version 2 qda file)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	qdnaOldMagic = 0xF656659EL	# big endian magic number for older qdna files
	qdnaMagic    = 0xC4B47197L	# big endian magic number for qdna files
	qdnaVersion  = 0x00000200L 

	# parse args

	name   = None
	strip  = False
	simple = False

	for arg in argv[1:]:
		if (arg.startswith("--name=")):
			name = arg.split("=",1)[1]
		elif (arg == "--striplinebreaks") or (arg == "--strip"):
			strip = True
		elif (arg == "--simple") or (arg == "--old"):
			simple = True
		elif (arg.startswith("--")):
			usage("can't understand %s" % arg)
		else:
			usage("can't understand %s" % arg)

	if (simple) and (name != None):
		uaseg("simple qdna file cannot carry a sequence name")

	# === read the input file ===

	seq = []
	for line in stdin:
		if (strip): line = line.rstrip()
		seq += [line]
	seq = "".join(seq)

	# === write the qdna file ===

	if (not simple):
		headerLen = 20
		if (name == None):
			nameOffset = 0
			seqOffset  = headerLen + 8;
		else:
			nameOffset = headerLen + 8;
			seqOffset  = nameOffset + len(name) + 1

	# prepend magic number

	if (simple): write_4(stdout,qdnaOldMagic)
	else:        write_4(stdout,qdnaMagic)

	# write the rest of the header

	if (not simple):
		write_4(stdout,qdnaVersion)
		write_4(stdout,headerLen)
		write_4(stdout,seqOffset)
		write_4(stdout,nameOffset)
		write_4(stdout,len(seq))
		write_4(stdout,0)

		if (name != None):
			stdout.write(name)
			stdout.write(chr(0))

	# write the sequence

	stdout.write(seq)


def write_4(f,val):
	f.write (chr((val >> 24) & 0xFF))
	f.write (chr((val >> 16) & 0xFF))
	f.write (chr((val >>  8) & 0xFF))
	f.write (chr( val        & 0xFF))


if __name__ == "__main__": main()
