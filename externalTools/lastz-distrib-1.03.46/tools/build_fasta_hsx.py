#!/usr/bin/env python
"""
Build a "hashed sequence index" (hsx) file for a fasta file
-----------------------------------------------------------

(see the header of hsx_file.py for file format details)

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys
from math import ceil
from hsx_file import HsxFile


def usage(s=None):
	message = """
build_fasta_hsx [options] [fasta_file ...] > hsx_file
  (if no fasta_files are present, we read fasta from stdin)
  --bucketsize=<N>  set the average hash bucket size
  --numbuckets=<N>  set the number of hash buckets (overrides avg size)
  --anonymous       don't copy fasta_file name into the index
  --secondary       use secondary hash in file instead of sequence names
  --skipheader      point to sequence data rather than header
  --windows         the fasta file has two-byte line feeds (this is what
                    microsoft windows uses)
  --bigendian       write fields as big endian (default is little endian)
  --oddbuckets      force the number of hash buckets to be odd
  --keepempties     don't discard empty sequences
  --progress[=<N>]  print progress reports on stderr
"""
	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():
	global write4,write5,write6

	##########
	# parse the command line
	##########

	avgBucket   = 10
	numBuckets  = None
	anonymous   = False
	doSecondary = False
	skipHeader  = False
	isWindows   = False
	fileNames   = []
	bigEndian   = False
	oddBuckets  = False
	keepEmpties = False
	screwup     = [] # (so that we can verify that validation works!)
	debug       = []
	progress    = None

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

		if (arg == "--bucketsize") and (val != None):
			try:
				avgBucket = int(val)
				if (avgBucket < 1): raise ValueError
			except ValueError:
				assert (False), "invalid bucket size: %s" % val
		elif (arg == "--numbuckets") and (val != None):
			try:
				numBuckets = int(val)
				if (numBuckets < 1): raise ValueError
			except ValueError:
				assert (False), "invalid number of buckets: %s" % val
		elif (arg == "--secondary") and (val == None):
			doSecondary = True
			assert (False), "secondary hash is not implemented yet (sorry)"
		elif (arg == "--anonymous") and (val == None):
			anonymous = True
		elif (arg == "--skipheader") and (val == None):
			skipHeader = True
		elif (arg == "--windows") and (val == None):
			isWindows = True
		elif (arg == "--bigendian") and (val == None):
			bigEndian = True
		elif (arg == "--oddbuckets") and (val == None):
			oddBuckets = True
		elif (arg == "--keepempties") and (val == None):
			keepEmpties = True
		elif (arg == "--screwup") and (val != None):
			screwup += [val]
		elif (arg == "--debug") and (val == None):
			debug += ["debug"]
		elif (arg == "--debug") and (val != None):
			debug += [val]
		elif (arg == "--progress") and (val == None):
			debug += ["progress"]
			progress = None
		elif (arg == "--progress") and (val != None):
			debug += ["progress"]
			progress = int(val)
		elif (arg.startswith("--")):
			usage("unknown argument: %s" % arg)
		elif (val == None):
			fileNames += [arg]
		else:
			usage("unknown argument: %s" % arg)

	# sanity check on file name

	if (fileNames != []):
		for fileName in fileNames:
			try:
				slash = fileName.rfind("/")
				dot   = fileName.rfind(".")
				if (dot < 0):     raise ValueError
				if (dot < slash): raise ValueError
				if (fileName[dot:] not in [".fa",".fasta"]): raise ValueError
			except ValueError:
				assert (False), \
				       "bad fasta file name (it has to end with .fa or .fasta)" \
				     % fileName

	if (anonymous) and (len(fileNames) > 1):
		assert (False), "can't use anonymous when you have multiple fasta files"

	assert (len(fileNames) <= 255), "too many input files (max is 255)"

	# set up big- or little-endian

	if (bigEndian):
		write4 = write4_big_endian
		write5 = write5_big_endian
		write6 = write6_big_endian
	else:
		write4 = write4_little_endian
		write5 = write5_little_endian
		write6 = write6_little_endian

	##########
	# read the fasta file(s)
	##########

	fileNameToNum = {}

	# read the fasta file(s), collecting names, etc.

	if (fileNames == []):
		fileNames += [""]

	sequences = []
	nameSeen  = {}

	for (fileNum,fileName) in enumerate(fileNames):
		assert (fileName not in fileNameToNum), \
		       "can't use the same file twice (%s)" % fileName
		fileNameToNum[fileName] = fileNum

		if (fileName == ""):
			f = sys.stdin
		else:
			try:
				f = file(fileName,"rt")
			except IOError:
				assert (False), "unable to open %s" % fileName

		seqNum = 0
		for seqInfo in fasta_sequences(f,twoByteLFs=isWindows):
			(name,length,lineNum,headerOffset,seqOffset) = seqInfo
			seqNum += 1

			assert (name not in nameSeen), \
				"%s is used for two sequences (at %s and %s)" \
			  % (name,
			     line_reference(nameSeen[name]),
			     line_reference((fileName,lineNum)))
			nameSeen[name] = (fileName,lineNum)

			if (length == 0):
				if (keepEmpties):
					print >>sys.stderr, "WARNING: keeping empty sequence %s (%s)" \
					                  % (name,line_reference((fileName,lineNum)))
				else:
					print >>sys.stderr, "WARNING: discarding empty sequence %s (%s)" \
					                  % (name,line_reference((fileName,lineNum)))
					continue

			if (skipHeader): sequences += [(name,length,fileNum,seqOffset)]
			else:            sequences += [(name,length,fileNum,headerOffset)]

			if ("progress" in debug) and (progress != None) and (seqNum % progress == 0):
				print >>sys.stderr, "read sequence %d (%s)" % (seqNum,name)

		if (fileName != ""): f.close()

		if ("progress" in debug):
			if (fileName != ""):
				print >>sys.stderr, "finished reading %s" % fileName
			else:
				print >>sys.stderr, "finished reading input file"

	# scan collected sequence info and assign hash values

	numSequences = len(sequences)
	assert (numSequences > 0), "input file contains no sequences!"
	if (numBuckets == None):
		numBuckets = int(ceil(numSequences / avgBucket))
	if (oddBuckets) and (numBuckets % 1 == 0):
		numBuckets += 1

	sequences = [(HsxFile.hash(name) % numBuckets,name,length,fileNum,offset) \
	                             for (name,length,fileNum,offset) in sequences] 
	sequences.sort()

	if ("progress" in debug):
		print >>sys.stderr, "finished computing hashes"

	if ("info" in debug):
		for (hash,name,length,fileNum,offset) in sequences:
			print >>sys.stderr, "%10d==%08X %2d:%08X %s %d" \
			                  % (HsxFile.hash(name),hash,fileNum,offset,name,length)

	##########
	# write the index
	##########

	# decide how we will write the file names

	fileNumToOffset    = {}
	fileNumToFastaName = {}
	fileNumToFastaExt  = {}
	fileInfoLength = 0

	for fileName in fileNames:
		fileNum = fileNameToNum[fileName]

		fastaName = ""
		fastaExt  = "fa"
		if (fileName != ""):
			dot = fileName.rfind(".")
			fastaExt = fileName[dot+1:]
			if (not anonymous):
				fastaName = fileName[:dot]

		fileNumToOffset   [fileNum] = fileInfoLength
		fileNumToFastaName[fileNum] = fastaName
		fileNumToFastaExt [fileNum] = fastaExt
		fileInfoLength += len(fastaExt)+1 + len(fastaName)+1 

	# determine header and table sizes

	headerLength    = 0x1C
	headerPad       = pad_for_16(8+headerLength)
	headerSize      = headerLength + headerPad

	numFiles        = len(fileNames)
	fileTableOffset = 0x08 + headerSize
	fileTableLength = numFiles * 4
	fileTablePad    = pad_for_16(fileTableLength)
	fileTableSize   = fileTableLength + fileTablePad

	fileInfoOffset  = fileTableOffset + fileTableSize
	fileInfoPad     = pad_for_16(fileInfoLength)
	fileInfoSize    = fileInfoLength + fileInfoPad

	hashTableOffset = fileInfoOffset + fileInfoSize
	hashTableLength = (numBuckets+1) * 5
	hashTablePad    = pad_for_16(hashTableLength)
	if ("hashpad" in screwup): hashTablePad = -1
	hashTableSize   = hashTableLength + hashTablePad

	seqTableOffset  = hashTableOffset + hashTableSize

	if ("file" in debug):
		print >>sys.stderr, "fileTableOffset = %08X (%08X)" % (fileTableOffset,fileTableSize)
		print >>sys.stderr, "fileInfoOffset  = %08X (%08X)" % (fileInfoOffset,fileInfoSize)
		print >>sys.stderr, "hashTableOffset = %08X (%08X)" % (hashTableOffset,hashTableSize)
		print >>sys.stderr, "seqTableOffset  = %08X"        % seqTableOffset

	# determine offsets into the sequence table

	nameToOffset = {}

	prevHash = None
	for (hash,name,length,fileNum,offset) in sequences:
		if (hash == prevHash): continue
		nameToOffset[name] = True

	seqOffset = seqTableOffset
	for (hash,name,length,fileNum,offset) in sequences:
		if (name in nameToOffset):
			nameToOffset[name] = seqOffset
		seqOffset += 12 + len(name) + 1
	nameToOffset[""] = seqOffset

	# write header

	write4(HsxFile.magicBig)
	write4(HsxFile.version)

	write4(headerLength)
	write4(numFiles)
	write4(fileTableOffset)
	write4(numBuckets)
	write4(hashTableOffset)
	write4(numSequences)
	write4(seqTableOffset)
	writeZeros(headerPad)

	if ("progress" in debug):
		print >>sys.stderr, "finished writing header"

	# write file table and file info

	for fileName in fileNames:
		fileNum = fileNameToNum[fileName]
		write4(fileInfoOffset + fileNumToOffset[fileNum])
	writeZeros(fileTablePad)

	for fileName in fileNames:
		fileNum = fileNameToNum[fileName]
		writeString(fileNumToFastaExt [fileNum])
		writeString(fileNumToFastaName[fileNum])
	writeZeros(fileInfoPad)

	if ("progress" in debug):
		print >>sys.stderr, "finished writing file table"

	# write hash table

	msBit5 = 0x80 << (4*8)

	prevHash = None

	for (hash,name,length,fileNum,offset) in sequences:
		if (hash == prevHash):
			bucketSize += 1
			continue

		if (prevHash != None):
			# output previous bucket
			write5(seqOffset)
			if ("progress" in debug) and (progress != None) and ((hash+1) % progress == 0):
				print >>sys.stderr, "wrote hash bucket %d" % (hash+1)
			# output intervening empty buckets
			prevHash += 1
			while (prevHash < hash):
				write5(msBit5 + nameToOffset[name])
				prevHash += 1
				if ("progress" in debug) and (progress != None) and (prevHash % progress == 0):
					print >>sys.stderr, "wrote hash bucket %d" % (prevHash)

		bucketSize = 1
		seqOffset  = nameToOffset[name]
		prevHash   = hash

	# output previous bucket
	write5(seqOffset)
	seqOffset = nameToOffset[""] # offset past end of sequence index table
	# output intervening empty buckets
	prevHash += 1
	while (prevHash < numBuckets):
		write5(msBit5 + seqOffset)
		prevHash += 1
	# output extra bucket
	write5(msBit5 + seqOffset)

	writeZeros(hashTablePad)

	if ("progress" in debug):
		print >>sys.stderr, "finished writing hash table"

	# write sequence table

	for (seqNum,(hash,name,length,fileNum,offset)) in enumerate(sequences):
		write5(length)			# length of the sequence
		write1(fileNum)			# file number (index into file table)
		write6(offset)			# offset to the sequence data
		writeString(name)		# name of sequence
		if ("progress" in debug) and (progress != None) and ((seqNum+1) % progress == 0):
			print >>sys.stderr, "wrote sequence entry %d" % (seqNum+1)

	if ("progress" in debug):
		print >>sys.stderr, "finished writing index"


def pad_for_16(n):
	return (16 - (n % 16)) % 16

def write1(val):
	sys.stdout.write(chr(val & 0xFF))

def writeString(s):
	assert (len(s) <= 255)
	sys.stdout.write(chr(len(s)))
	sys.stdout.write(s)

def writeZeros(n):
	for i in range(n): sys.stdout.write(chr(0))

def write4_little_endian(val):
	sys.stdout.write(chr( val        & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >> 24) & 0xFF))

def write5_little_endian(val):
	sys.stdout.write(chr( val        & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >> 24) & 0xFF))
	sys.stdout.write(chr((val >> 32) & 0xFF))

def write6_little_endian(val):
	sys.stdout.write(chr( val        & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >> 24) & 0xFF))
	sys.stdout.write(chr((val >> 32) & 0xFF))
	sys.stdout.write(chr((val >> 40) & 0xFF))

def write4_big_endian(val):
	sys.stdout.write(chr((val >> 24) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr( val        & 0xFF))

def write5_big_endian(val):
	sys.stdout.write(chr((val >> 32) & 0xFF))
	sys.stdout.write(chr((val >> 24) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr( val        & 0xFF))

def write6_big_endian(val):
	sys.stdout.write(chr((val >> 40) & 0xFF))
	sys.stdout.write(chr((val >> 32) & 0xFF))
	sys.stdout.write(chr((val >> 24) & 0xFF))
	sys.stdout.write(chr((val >> 16) & 0xFF))
	sys.stdout.write(chr((val >>  8) & 0xFF))
	sys.stdout.write(chr( val        & 0xFF))

# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f,nameParse=None,twoByteLFs=False):

	lineNum    = 0
	fileOffset = 0
	seqName    = None
	seqLength  = 0

	for line in f:
		lineNum += 1
		lineOffset = fileOffset
		fileOffset += len(line)
		if (twoByteLFs): fileOffset += 1
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				if (seqOffset == None): seqOffset = lineOffset
				yield (seqName,seqLength,seqLine,headerOffset,seqOffset)
			seqLine      = lineNum
			headerOffset = lineOffset
			seqName      = sequence_name(line)
			seqLength    = 0
			seqOffset    = None
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			if (seqOffset == None): seqOffset = lineOffset
			seqLength += len(line)

	if (seqName != None):
		if (seqOffset == None): seqOffset = fileOffset
		yield (seqName,seqLength,seqLine,headerOffset,seqOffset)


# sequence_name--
#	Extract the sequence name from a fasta header.
#	$$$ this needs to use nameParse

def sequence_name(s,nameParse=None):
	s = s[1:].strip()
	if (s == ""): return ""
	else:         return s.split()[0]


def line_reference((fileName,lineNum)):
	if (fileName == ""): return "line %d" % lineNum
	else:                return "line %s:%d" % (fileName,lineNum)


if __name__ == "__main__": main()

