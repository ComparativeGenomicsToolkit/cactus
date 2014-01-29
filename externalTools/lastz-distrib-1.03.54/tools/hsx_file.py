#!/usr/bin/env python
"""
"Hashed sequence index" (hsx) file reader (for a fasta file)
-------------------------------------------------------------------

offset 0x00: D2 52 70 95        big endian magic number
								.. (95 70 52 D2 => little endian)
offset 0x04: 00 00 01 xx        version 1.0 (see note 1)
offset 0x08: 00 00 00 1C        header length (in bytes, including this
								.. field)
offset 0x0C: xx xx xx xx        FN, number of files (see note 2)
offset 0x10: xx xx xx xx        FO, offset to file table
offset 0x14: xx xx xx xx        HN, number of hash buckets (see notes 3 and 4)
offset 0x18: xx xx xx xx        HO, offset to hash table
offset 0x1C: xx xx xx xx        SN, number of sequences
offset 0x20: xx xx xx xx        SO, offset to sequence index table (see
								.. note 5)

offset FO:   xx xx xx xx        FIO0, offset to file info for file 0
			  ...               (FN-1 more entries, at 4 bytes per)

offset FIOn: LL xx ..           type of file (ascii "fa", "2bit", etc., see
                                note 6)
			 LL xx ..           name of file (see note 7)
			  ...               (FN-1 more entries, variable length)

offset HO:   xx xx xx xx xx     SIOn, offset into sequence index table (see
								.. notes 8, 9 and 10)
			  ...               (HN-1 more entries, at 5 bytes per)
			 xx xx xx xx xx     offset past end of sequence index table

offset SO:   xx xx xx xx xx     length of the sequence (see note 11)
			 xx                 file number (index into file table)
			 xx xx xx xx xx xx  offset to the sequence data (see note 12)
			 LL xx ..           name of sequence (see note 13)
			  ...               (SN-1 more entries, variable length)
 
Notes:

	(1)  The least significant byte of the version is the "sub version".
	     For version 1, this is 00 (secondary hashes are not in use) or 01
	     (secondary hashes are in use).
	(2)  The number of files is limited to 255.
	(3)  It is assumed that the number of buckets is set so that the average
	     number of sequences per bucket (SN/HN) is reasonably small (e.g. 10).
	(4)  The hash table actually includes HN+1 buckets.  The extra bucket has
	     size zero and gives the offset to just past the end of the sequence
	     index table.
	(5)  Entries in the sequence index table are necessarily stored in hash
	     order.  Entries with the same hash are stored in alphabetical order;
	     actually, in lexicographic order over the bytes of their names.
	(6)  Strings are stored as a length byte followed by ascii text.
	(7)  If a file info record contains an empty name, the name of the file is
	     the same as the index file itself, with the file type used as the
	     extension (e.g. "reads.hsx" becomes "reads.fa").  This allows files to
	     be renamed without rebuilding the index.
	(8)  SIOn is the file offset for the nth entry in the sequence index table.
	     When this is in a hash table entry, it is the index for the first
	     sequence in that hash's bucket.
	(9)  The most significant bit in a bucket's SIOn value is used to indicate
	     whether the bucket is empty or not.  If a bucket is empty, this bit is
	     set (1), otherwise it is clear.
	(10) The end of a bucket can be determined from the SIOn entry for the
	     start of the next bucket.
	(11) A sequence may be empty, so zero is a legitimate value for the
	     sequence length.
	(12) The offset to the sequence data is an offset into the sequence file.
	     For fasta it can point to the ">" at the start of the sequence's
	     header, or directly to the sequence data.
	(13) When secondary hashes are in use, the sequence name (including the
	     terminating zero) is replaced by the four-byte secondary hash.

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys,struct
import hassock_hash


class HsxFile(object):

	def __init__(self,fileName,debug=None):
		self.fileName = fileName
		self.file     = None
		self.numFiles = 0
		if (debug == None): self.debug = []
		else:               self.debug = debug
		self.open()

	magicBig    = 0xD2527095L
	magicLittle = 0x957052D2L
	version     = 0x00000100L
	msBit5      = 0x80 << (4*8)

	def open(self):
		self.file = file(self.fileName,"rb")

		self.magic = magic = struct.unpack(">L",self.file.read(4))[0]
		if   (magic == HsxFile.magicBig):    self.byteOrder = ">" # (big endian)
		elif (magic == HsxFile.magicLittle): self.byteOrder = "<" # (little endian)
		else:
			assert (False), \
			       "%s is not an hsx file (magic = %08X)" \
			     % (self.fileName,magic)
		self.struct4 = "%sL" % self.byteOrder

		self.version = self.read4()
		assert (self.version == HsxFile.version), \
		       "%s is hsx version %08X, which is not supported" \
		     % (self.fileName,self.version)

		self.read_header()
		self.load_file_table()

	def close(self):
		self.file.close()
		for fileIx in range(self.numFiles):
			(name,file) = self.fileTable[fileIx]
			if (file != None): file.close()

	def read_header(self):
		self.headerLength = self.read4()
		assert (self.headerLength >= 0x1C), \
		       "%s has unsupported header length (%08X)" \
		     % (self.fileName,self.headerSize)
		(self.numFiles,
		 self.fileTableOffset,
		 self.numBuckets,
		 self.hashTableOffset,
		 self.numSequences,
		 self.seqTableOffset) = struct.unpack("%sLLLLLL" % self.byteOrder,self.file.read(24))
		assert (self.numBuckets != 0), \
		       "%s has corrupt header (numBuckets = 0)" % (self.fileName)

	def load_file_table(self):
		self.file.seek(self.fileTableOffset)
		offsetTable = self.file.read(4*self.numFiles)
		offsetTable = struct.unpack("%s%s" % (self.byteOrder,"L"*self.numFiles),offsetTable)
		self.fileTable = [None] * self.numFiles

		basePath = baseName = None
		for fileIx in range(self.numFiles):
			self.file.seek(offsetTable[fileIx])
			extension = self.readString()
			name      = self.readString()
			if (name == ""):
				if (baseName == None):
					baseName = self.base_file_name()
				name = baseName + "." + extension
			else:
				if (basePath == None):
					basePath = self.base_file_path()
				name = basePath + name + "." + extension
			self.fileTable[fileIx] = (name,None) # (second field holds file when opened)
			#.. print "fileTable[%d] = %s" % (fileIx,name)

	def base_file_name(self):
		slash = self.fileName.rfind("/")
		dot   = self.fileName.rfind(".")
		if (dot < 0):     return self.fileName
		if (dot < slash): return self.fileName
		return self.fileName[:dot]

	def base_file_path(self):
		slash = self.fileName.rfind("/")
		if (slash < 0): return ""
		return self.fileName[:slash+1]

	def get_sequence(self,name):
		if ("fetch" in self.debug):
			print >>sys.stderr, "[fetching %s]" % name
		# read hash bucket for this name
		bucket = HsxFile.hash(name) % self.numBuckets
		if ("fetch" in self.debug):
			print >>sys.stderr, "[  bucket = %d (file offset %08X)]" \
			                  % (bucket,self.hashTableOffset+5*bucket)
		self.file.seek(self.hashTableOffset + 5*bucket)
		bucketOffset = self.read5()
		if (bucketOffset & HsxFile.msBit5 != 0):
			if ("fetch" in self.debug):
				print >>sys.stderr, "[  bucket is empty]"
			return None
		bucketEnd = self.read5() & ~HsxFile.msBit5
		if ("fetch" in self.debug):
			print >>sys.stderr, "[  bucket offset = %010X..%010X ]" \
			                  % (bucketOffset,bucketEnd)
		# scan the bucket until we find this sequence
		self.file.seek(bucketOffset)
		seqIx   = 1
		seqName = None
		while (bucketOffset < bucketEnd):
			seqLength = self.read5()
			fileIx    = self.read1()
			seqOffset = self.read6()
			seqName   = self.readString()
			if ("fetch" in self.debug):
				print >>sys.stderr, "[  (%010X) name %d = %s]" \
				                    % (bucketOffset,seqIx,seqName)
			if (seqName == name): break
			if (seqName >  name): return None
			bucketOffset += 1 + 6 + 5 + len(seqName) + 1
			seqIx += 1
		if (seqName != name):
			if ("fetch" in self.debug):
				print >>sys.stderr, "[  %s not in bucket]" % name
			return None
		# open the sequence file (if it isn't already open)
		assert (fileIx < len(self.fileTable)), \
		       "file index for %s is out of bounds (%d > %d)" \
		     % (name,fileIx,len(self.fileTable))
		(seqFileName,seqFile) = self.fileTable[fileIx]
		if (seqFile == None):
			if ("fetch" in self.debug):
				print >>sys.stderr, "[  opening %s]" % seqFileName
			seqFile = file(seqFileName,"rt")
			self.fileTable[fileIx] = (seqFileName,seqFile)
		if ("fetch" in self.debug):
			print >>sys.stderr, "[  reading from %s:%012X]" \
			                  % (seqFileName,seqOffset)
		# read the sequence
		seqFile.seek(seqOffset)
		seqLines = []
		seqRead = 0
		while (True):
			line = seqFile.readline()
			if (line == ""): break
			line = line.strip()
			if ("fetch" in self.debug):
				print >>sys.stderr, "[  read %s]" % line
			if (line.startswith(">")):
				if (len(seqLines) != 0): break
				seqLines += [line]
				continue
			seqRead += len(line)
			if (seqRead > seqLength):
				line = line[:-seqLength-seqRead]
				seqRead = seqLength
			seqLines += [line]
			if (seqRead == seqLength):
				break
		assert (seqRead == seqLength), \
		       "sequence for %s is short (%d < %d)" \
		     % (name,seqRead,seqLength)
		return "\n".join(seqLines)

	def read1(self):
		return ord(self.file.read(1))

	def read4(self):
		return struct.unpack(self.struct4,self.file.read(4))[0]

	def read5(self):
		return self.read_and_unpack(5)

	def read6(self):
		return self.read_and_unpack(6)

	def readString(self):
		ch = self.file.read(1)
		s = self.file.read(ord(ch))
		return "".join(s)

	def read_and_unpack(self,bytes):
		data = self.file.read(bytes)
		if (self.byteOrder == "<"): # (make data big endian)
			data = [ch for ch in data]
			data.reverse()
		val = 0
		for ch in data: val = (val << 8) + ord(ch)
		return val

	# hash

	def hash(name):
		return hassock_hash.hassock_hash(name)
	hash = staticmethod(hash)


if __name__ == "__main__": main()

