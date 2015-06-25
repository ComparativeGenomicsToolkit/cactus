//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: capsule.c
//
//----------
//
// capsule--
//	Support for "capsule" files and sharing of target data structures between
//	multiple processes.
//
// Sharing is achieved through the unix mmap function.  We mmap to a single
// file containing (a) the target sequence and (b) the seed word position
// table.
//
//----------

//----------
//
// other files
//
//----------

#include <unistd.h>				// standard UNIX stuff (non ANSI)
#include <stdio.h>				// standard C i/o stuff
#include <sys/stat.h>			// some mysterious UNIX voodoo
#include <sys/mman.h>			// UNIX memory manager stuff (non ANSI)
#include <fcntl.h>				// UNIX file control stuff (non ANSI)
#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <string.h>				// standard C string stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed matching stuff
#include "pos_table.h"			// position table stuff

#define  capsule_owner			// (make this the owner of its globals)
#include "capsule.h"			// interface to this module

// debugging defines

//#define snoopBytesWritten		// if this is defined, extra code is added to
								// .. track the number of bytes written to the
								// .. capsule file

//----------
//
// write_capsule_file--
//	Write a Target Sequence Capsule File corresponding to the current target
//	sequence and seeding state.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to write.  The caller should have already
//							.. opened this, with "wb" access.
//	char*		filename:	The name of the file being written to.  This is
//							.. only used for error reporting, and may be NULL.
//	seq*		seq:		The target sequence.
//	u8*			revNucs:	The reverse of the target sequence (NOT reverse
//							.. complement);  this may be NULL, in which case it
//							.. is left out of the file.
//	postable*	pt:			A table of positions of words in target.
//	seed*		seed:		The seed used to build the position table.
//
// Returns:
//	The number of bytes written to the file.
//
//----------

//=== stuff for snoopBytesWritten ===

#ifndef snoopBytesWritten
#define debugSnoopBytesWritten_1 ;
#define debugSnoopBytesWritten_2 ;
#define debugSnoopBytesWritten_3(length,bytes) ;
#define debugSnoopBytesWritten_4 ;
#define debugSnoopBytesWritten_5 ;
#define debugSnoopBytesWritten_6 ;
#endif // not snoopBytesWritten

#ifdef snoopBytesWritten

#define debugSnoopBytesWritten_1                                               \
	fprintf (stderr,"write_field(0x%016" PRIX64 ",%s)\n",                      \
	                bytesToWrite, reason);

#define debugSnoopBytesWritten_2                                               \
	fprintf (stderr,"write_sized_field(0x%016" PRIX64 ",%s)\n",                \
	                bytesToWrite, reason);

#define debugSnoopBytesWritten_3(length,bytes)                                 \
	fprintf (stderr,"write_padding(0x%016" PRIX64 ",0x%016" PRIX64 ",%s)\n",   \
	                (u64) length, (u64) bytes, reason);

#define debugSnoopBytesWritten_4                                               \
	fprintf (stderr,"  bytesWritten      = 0x%016" PRIX64 "\n", bytesWritten); \
	fprintf (stderr,"  totalBytesWritten = 0x%016" PRIX64 "\n", totalBytesWritten);

#define debugSnoopBytesWritten_5                                               \
	fprintf (stderr,"  sizeof(size_t) = %d\n", (int) sizeof(size_t));          \
	fprintf (stderr,"  sizeof(last[0]) = %d\n", (int) sizeof(pt->last[0]));    \
	fprintf (stderr,"  allocLast  = 0x%016" PRIX64 "\n", (u64) pt->allocLast); \
	fprintf (stderr,"  lastLength = 0x%016" PRIX64 "\n", lastLength);          \
	fprintf (stderr,"  lastBytes  = 0x%016" PRIX64 "\n", lastBytes);           \
	fprintf (stderr,"  sizeof(prev[0]) = %d\n", (int) sizeof(pt->prev[0]));    \
	fprintf (stderr,"  allocPrev  = 0x%016" PRIX64 "\n", (u64) pt->allocPrev); \
	fprintf (stderr,"  prevLength = 0x%016" PRIX64 "\n", prevLength);          \
	fprintf (stderr,"  prevBytes  = 0x%016" PRIX64 "\n", prevBytes);

#define debugSnoopBytesWritten_6                                               \
	fprintf (stderr,"  header length: 0x%016" PRIX64 "\n", (u64) headerLength);\
	fprintf (stderr,"   header bytes: 0x%016" PRIX64 "\n", (u64) headerBytes); \
	fprintf (stderr,"    name offset: 0x%016" PRIX64 "\n", (u64) nameOffset);  \
	fprintf (stderr,"    name length: 0x%016" PRIX64 "\n", (u64) nameLength);  \
	fprintf (stderr,"     name bytes: 0x%016" PRIX64 "\n", (u64) nameBytes);   \
	fprintf (stderr,"    nucs offset: 0x%016" PRIX64 "\n", (u64) nucsOffset);  \
	fprintf (stderr,"    nucs length: 0x%016" PRIX64 "\n", (u64) nucsLength);  \
	fprintf (stderr,"     nucs bytes: 0x%016" PRIX64 "\n", (u64) nucsBytes);   \
	fprintf (stderr,"    rvrs offset: 0x%016" PRIX64 "\n", (u64) rvrsOffset);  \
	fprintf (stderr,"    rvrs length: 0x%016" PRIX64 "\n", (u64) rvrsLength);  \
	fprintf (stderr,"     rvrs bytes: 0x%016" PRIX64 "\n", (u64) rvrsBytes);   \
	fprintf (stderr,"    bits offset: 0x%016" PRIX64 "\n", (u64) bitsOffset);  \
	fprintf (stderr,"    bits length: 0x%016" PRIX64 "\n", (u64) bitsLength);  \
	fprintf (stderr,"     bits bytes: 0x%016" PRIX64 "\n", (u64) bitsBytes);   \
	fprintf (stderr,"    last offset: 0x%016" PRIX64 "\n", (u64) lastOffset);  \
	fprintf (stderr,"    last length: 0x%016" PRIX64 "\n", (u64) lastLength);  \
	fprintf (stderr,"     last bytes: 0x%016" PRIX64 "\n", (u64) lastBytes);   \
	fprintf (stderr,"    prev offset: 0x%016" PRIX64 "\n", (u64) prevOffset);  \
	fprintf (stderr,"    prev length: 0x%016" PRIX64 "\n", (u64) prevLength);  \
	fprintf (stderr,"     prev bytes: 0x%016" PRIX64 "\n", (u64) prevBytes);   \
	fprintf (stderr,"    info offset: 0x%016" PRIX64 "\n", (u64) infoOffset);  \
	fprintf (stderr,"    info length: 0x%016" PRIX64 "\n", (u64) infoLength);  \
	fprintf (stderr,"     info bytes: 0x%016" PRIX64 "\n", (u64) infoBytes);   \
	fprintf (stderr,"    part offset: 0x%016" PRIX64 "\n", (u64) partOffset);  \
	fprintf (stderr,"    part length: 0x%016" PRIX64 "\n", (u64) partLength);  \
	fprintf (stderr,"     part bytes: 0x%016" PRIX64 "\n", (u64) partBytes);   \
	fprintf (stderr,"    pool offset: 0x%016" PRIX64 "\n", (u64) poolOffset);  \
	fprintf (stderr,"    pool length: 0x%016" PRIX64 "\n", (u64) poolLength);  \
	fprintf (stderr,"     pool bytes: 0x%016" PRIX64 "\n", (u64) poolBytes);   \
	fprintf (stderr,"    seed offset: 0x%016" PRIX64 "\n", (u64) seedOffset);  \
	fprintf (stderr,"    seed length: 0x%016" PRIX64 "\n", (u64) seedLength);  \
	fprintf (stderr,"     seed bytes: 0x%016" PRIX64 "\n", (u64) seedBytes);   \
	fprintf (stderr,"     end offset: 0x%016" PRIX64 "\n", (u64) endOffset);

#endif // snoopBytesWritten


//=== macros to write fields ===

#define write_field(fieldName)                                               \
	bytesToWrite = sizeof(fieldName);                                        \
	debugSnoopBytesWritten_1;                                                \
	bytesWritten = fwrite (&fieldName, 1, bytesToWrite, f);                  \
	if (bytesWritten != bytesToWrite) goto write_failure;                    \
	totalBytesWritten += bytesWritten;                                       \
    debugSnoopBytesWritten_4;

#define write_sized_field(fieldName,bytes)                                   \
	bytesToWrite = bytes;                                                    \
	debugSnoopBytesWritten_2;                                                \
	bytesWritten = fwrite (fieldName, 1, bytesToWrite, f);                   \
	if (bytesWritten != bytesToWrite) goto write_failure;                    \
	totalBytesWritten += bytesWritten;                                       \
    debugSnoopBytesWritten_4;

#define write_padding(length,bytes)                                          \
	debugSnoopBytesWritten_3(length,bytes);                                  \
	if (bytes > length)                                                      \
		{                                                                    \
		bytesToWrite = bytes - length;                                       \
		bytesWritten = fwrite (zeroes, 1, bytesToWrite, f);                  \
		if (bytesWritten != bytesToWrite) goto write_failure;                \
		totalBytesWritten += bytesWritten;                                   \
		debugSnoopBytesWritten_4;                                            \
		}


//=== write_capsule_file ===

u64 write_capsule_file
   (FILE*			f,
	char*			filename,
	seq*			seq,
	u8*				revNucs,
	postable*		pt,
	seed*			seed)
	{
	seqpartition*	sp = &seq->partition;
	u8				zeroes[32];
	u64				totalBytesWritten = 0;
	u64				bytesToWrite, bytesWritten;
	u32				headerLength, headerBytes; // (yes, u32 not u64)
	u64				nameLength, nameBytes, nameOffset;
	u64				nucsLength, nucsBytes, nucsOffset;
	u64				rvrsLength, rvrsBytes, rvrsOffset;
	u64				bitsLength, bitsBytes, bitsOffset;
	u64				lastLength, lastBytes, lastOffset;
	u64				prevLength, prevBytes, prevOffset;
	u64				infoLength, infoBytes, infoOffset;
	u64				partLength, partBytes, partOffset;
	u64				poolLength, poolBytes, poolOffset;
	u64				seedLength, seedBytes, seedOffset;
	u64				endOffset, badOffset;
	u64				magic;
	u32				version, headerEntries;
	u32				dataTypeCode, dataItem, extraInfo;
	char*			seqName, *reason;
	u32*			flipScan;
	int				numFlips, flipIx, partIx;

	if (filename == NULL) filename = "unnamed capsule file";

	memset (zeroes, 0, sizeof(zeroes));

	if (sizeof(unspos) != sizeof(u32))
		suicide ("internal error, capsule expects positions to be 32 bits");

	if (sizeof(partition) != sizeof(cappartition))
		suicidef ("internal error, capsule expects partition records to be %d bytes",
		          sizeof(cappartition));

	//////////
	// write magic number
	//////////

	reason = "magic";
	magic = (((u64) refcapMagicABig) << 32) + refcapMagicBBig;
	write_field (magic);

	//////////
	// figure out block sizes and offsets
	//////////

	// figure out the header size (not including the 12 bytes for magic number
	// and version)

	headerEntries = 6;
	if (revNucs    != NULL) headerEntries++;
	if (pt->asBits != NULL) headerEntries++;
	if (sp->p      != NULL) headerEntries += 2;

	headerLength  = sizeof(u32)									// length field
	              + (headerEntries * capsuleHeaderEntrySize)	// entries
	              + sizeof(u32);								// terminator

	headerBytes = round_up_32 (headerLength + capsulePreHeaderSize)
	            - capsulePreHeaderSize;

	// figure out the name block's size

	seqName = (seq->useFullNames)? seq->header : seq->shortHeader;
	if ((seqName == NULL) || (seqName[0] == 0)) seqName = "(unnamed)";
	nameLength = strlen(seqName) + 1;
	nameBytes  = round_up_32 (nameLength);

	// figure out the nucleotide blocks' size

	nucsLength = seq->len + 1;
	nucsBytes  = round_up_32 (nucsLength);

	rvrsLength = (revNucs == NULL)? 0 : nucsLength;
	rvrsBytes  = round_up_32 (rvrsLength);

	bitsLength = (pt->asBits == NULL)? 0 : round_up_16((nucsLength+3) / 4);
	bitsBytes  = round_up_32 (bitsLength);

	infoLength = sizeof(capseqinfo);
	infoBytes  = round_up_32 (infoLength);

	// figure out the position table blocks' size

	lastLength = ((size_t) pt->allocLast) * sizeof(pt->last[0]);
	lastBytes  = round_up_32 (lastLength);

	prevLength = ((size_t) pt->allocPrev) * sizeof(pt->prev[0]);
	prevBytes  = round_up_32 (prevLength);

	debugSnoopBytesWritten_5;

	// figure out the partition blocks' sizes

	partBytes = poolBytes = 0;
	if (sp->p != NULL)
		{
		partLength = (sp->len + 1) * sizeof(partition);
		partBytes  = round_up_32 (partLength);

		poolLength = sp->poolLen;
		poolBytes  = round_up_32 (poolLength);
		}

	// figure out the seed block's size

	numFlips = 0;
	for (flipScan=seed->transFlips ; *flipScan!=0 ; flipScan++)
		numFlips++;

	seedLength = (sizeof(capseed) - sizeof(u32))		// standard fields
	           + (seed->numParts * sizeof(u32))			// shift[] array
	           + (seed->numParts * sizeof(u32))			// mask[] array
	           + ((numFlips+1) * sizeof(u32));			// transFlips[] array
	seedBytes  = round_up_32 (seedLength);

	// figure out all the offsets

	nameOffset = capsulePreHeaderSize + headerBytes;
	nucsOffset = nameOffset + nameBytes;
	rvrsOffset = nucsOffset + nucsBytes;
	bitsOffset = rvrsOffset + rvrsBytes;
	lastOffset = bitsOffset + bitsBytes;
	prevOffset = lastOffset + lastBytes;
	infoOffset = prevOffset + prevBytes;
	partOffset = infoOffset + infoBytes;
	poolOffset = partOffset + partBytes;
	seedOffset = poolOffset + poolBytes;
	endOffset  = seedOffset + seedBytes;

	capsule_set_stat  (display, true);
	capsule_copy_stat (headerLength);
	capsule_copy_stat (headerBytes);
	capsule_copy_stat (nameOffset);
	capsule_copy_stat (nameLength);
	capsule_copy_stat (nameBytes);
	capsule_copy_stat (nucsOffset);
	capsule_copy_stat (nucsLength);
	capsule_copy_stat (nucsBytes);
	capsule_copy_stat (rvrsOffset);
	capsule_copy_stat (rvrsLength);
	capsule_copy_stat (rvrsBytes);
	capsule_copy_stat (bitsOffset);
	capsule_copy_stat (bitsLength);
	capsule_copy_stat (bitsBytes);
	capsule_copy_stat (lastOffset);
	capsule_copy_stat (lastLength);
	capsule_copy_stat (lastBytes);
	capsule_copy_stat (prevOffset);
	capsule_copy_stat (prevLength);
	capsule_copy_stat (prevBytes);
	capsule_copy_stat (infoOffset);
	capsule_copy_stat (infoLength);
	capsule_copy_stat (infoBytes);
	capsule_copy_stat (partOffset);
	capsule_copy_stat (partLength);
	capsule_copy_stat (partBytes);
	capsule_copy_stat (poolOffset);
	capsule_copy_stat (poolLength);
	capsule_copy_stat (poolBytes);
	capsule_copy_stat (seedOffset);
	capsule_copy_stat (seedLength);
	capsule_copy_stat (seedBytes);
	capsule_copy_stat (endOffset);

	debugSnoopBytesWritten_6;

	//////////
	// finish writing the pre-header
	//////////

	// write the file size

	reason = "file size";
	write_field (endOffset);

	// write the version

	reason = "version";
	version = refcapVersion;
	write_field (version);

	// write the header length

	reason = "header";
	write_field (headerLength);

	//////////
	// write the header
	//////////

	// write the sequence name

	reason       = "name entry";
	dataTypeCode = cap_seqName;       write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (nameOffset);
	                                  write_field (nameLength);

	// write the nucleotides

	reason       = "nucs entry";
	dataTypeCode = cap_seqForward;    write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (nucsOffset);
	                                  write_field (nucsLength);

	// write the reversed nucleotides

	if (rvrsBytes > 0)
		{
		reason       = "rvrs entry";
		dataTypeCode = cap_seqReverse;  write_field (dataTypeCode);
		extraInfo    = 0;               write_field (extraInfo);
		                                write_field (rvrsOffset);
		                                write_field (rvrsLength);
		}

	// write the nucleotide bits

	if (bitsBytes > 0)
		{
		reason       = "nuc bits entry";
		dataTypeCode = cap_seqBits;     write_field (dataTypeCode);
		extraInfo    = 0;               write_field (extraInfo);
		                                write_field (bitsOffset);
		                                write_field (bitsLength);
		}

	// write the last[] array

	reason       = "last entry";
	dataTypeCode = cap_lastPosTable;  write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (lastOffset);
	                                  write_field (lastLength);

	// write the prev[] array

	reason       = "prev entry";
	dataTypeCode = cap_prevPosTable;  write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (prevOffset);
	                                  write_field (prevLength);

	// write the sequence info

	reason       = "info entry";
	dataTypeCode = cap_seqInfo;       write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (infoOffset);
	                                  write_field (infoLength);

	// write the partition[] array

	if (partBytes > 0)
		{
		reason       = "parititon entry";
		dataTypeCode = cap_partitions; write_field (dataTypeCode);
		extraInfo    = 0;             write_field (extraInfo);
		                              write_field (partOffset);
		                              write_field (partLength);
		}

	// write the partition names[] array

	if (poolBytes > 0)
		{
		reason       = "parititon entry";
		dataTypeCode = cap_partitionNames; write_field (dataTypeCode);
		extraInfo    = 0;             write_field (extraInfo);
		                              write_field (poolOffset);
		                              write_field (poolLength);
		}

	// write the seed

	reason       = "seed entry";
	dataTypeCode = cap_seed;          write_field (dataTypeCode);
	extraInfo    = 0;                 write_field (extraInfo);
	                                  write_field (seedOffset);
	                                  write_field (seedLength);

	// write the terminator and padding

	reason       = "terminator";
	dataTypeCode = cap_terminator;    write_field (dataTypeCode);

	reason       = "header padding";
	write_padding (headerLength, headerBytes);

	//////////
	// write the data blocks
	//////////

	badOffset = 0; // (placate complier)

	// write the sequence name

	reason = "name";
	if (totalBytesWritten != nameOffset)
		{ badOffset = nameOffset;  goto wrong_offset; }

	write_sized_field (seqName, nameLength);
	write_padding     (nameLength, nameBytes);

	// write the nucleotides

	reason = "nucs";
	if (totalBytesWritten != nucsOffset)
		{ badOffset = nucsOffset;  goto wrong_offset; }

	write_sized_field (seq->v, nucsLength);
	write_padding     (nucsLength, nucsBytes);

	// write the reversed nucleotides

	if (rvrsBytes > 0)
		{
		reason = "rvrs";
		if (totalBytesWritten != rvrsOffset)
			{ badOffset = rvrsOffset;  goto wrong_offset; }

		write_sized_field (revNucs, rvrsLength);
		write_padding     (rvrsLength, rvrsBytes);
		}

	// write the nucleotide bits

	if (bitsBytes > 0)
		{
		reason = "bits";
		if (totalBytesWritten != bitsOffset)
			{ badOffset = bitsOffset;  goto wrong_offset; }

		write_sized_field (pt->asBits, bitsLength);
		write_padding     (bitsLength, bitsBytes);
		}

	// write the last[] array

	reason = "last";
	if (totalBytesWritten != lastOffset)
		{ badOffset = lastOffset;  goto wrong_offset; }

	write_sized_field (pt->last, lastLength);
	write_padding     (lastLength, lastBytes);

	// write the prev[] array

	reason = "prev";
	if (totalBytesWritten != prevOffset)
		{ badOffset = prevOffset;  goto wrong_offset; }

	write_sized_field (pt->prev, prevLength);
	write_padding     (prevLength, prevBytes);

	// write the sequence info

	reason = "info";
	if (totalBytesWritten != infoOffset)
		{ badOffset = infoOffset;  goto wrong_offset; }

	dataItem = seq->startLoc;        write_field (dataItem);
	dataItem = seq->trueLen;         write_field (dataItem);
	dataItem = seq->revCompFlags;    write_field (dataItem);
	dataItem = seq->contig;          write_field (dataItem);

	dataItem = (seq->partition.p == NULL)? 0 : seq->partition.len;
	write_field (dataItem);

	write_padding (infoLength, infoBytes);

	// write the partitions and names

	if (partBytes > 0)
		{
		reason = "part";
		if (totalBytesWritten != partOffset)
			{ badOffset = partOffset;  goto wrong_offset; }

		write_sized_field (sp->p, partLength);
		write_padding     (partLength, partBytes);
		}

	if (poolBytes > 0)
		{
		reason = "pool";
		if (totalBytesWritten != poolOffset)
			{ badOffset = poolOffset;  goto wrong_offset; }

		write_sized_field (sp->pool, poolLength);
		write_padding     (poolLength, poolBytes);
		}

	// write the seed

	reason = "seed";
	if (totalBytesWritten != seedOffset)
		{ badOffset = seedOffset;  goto wrong_offset; }

	dataItem = pt->step;             write_field (dataItem);
	dataItem = seed->type;           write_field (dataItem);
	dataItem = seed->length;         write_field (dataItem);
	dataItem = seed->weight;         write_field (dataItem);
	dataItem = seed->resolvingMask;  write_field (dataItem);
	dataItem = seed->revComp;        write_field (dataItem);
	dataItem = seed->isHalfweight;   write_field (dataItem);
	dataItem = seed->numParts;       write_field (dataItem);

	for (partIx=0 ; partIx<seed->numParts ; partIx++)
		{ dataItem = seed->shift[partIx];  write_field (dataItem); }

	for (partIx=0 ; partIx<seed->numParts ; partIx++)
		{ dataItem = seed->mask[partIx];  write_field (dataItem); }

	for (flipIx=0 ; flipIx<numFlips ; flipIx++)
		{ dataItem = seed->transFlips[flipIx];  write_field (dataItem); }

	dataItem = 0;  write_field (dataItem);

	write_padding (seedLength, seedBytes);

	// sanity check on file length

	if (totalBytesWritten != endOffset)
		goto wrong_file_length;

	// success!

	return endOffset;

	//////////
	// failure exits
	//////////

wrong_offset:
	suicidef ("internal error writing to %s (offset for %s = 0x%s, actual is 0x%s)",
	          filename, reason, hex_64_string(badOffset), hex_64_string(totalBytesWritten));

wrong_file_length:
	suicidef ("internal error writing to %s (file length = 0x%s, actual is 0x%s)",
	          filename, hex_64_string(endOffset), hex_64_string(totalBytesWritten));

write_failure:
	suicidef_with_perror ("unable to write to %s (attempted %d bytes, wrote %d, for %s)",
	                      filename, bytesToWrite, bytesWritten, reason);
	return 0; // (never gets here)
	}

//----------
//
// open_capsule_file--
//	Open a Target Sequence Capsule File and map it for sharing.
//
//----------
//
// Arguments:
//	char*	filename:		The name of the capsule file.
//
// Returns:
//	A pointer to a capsule info record, allocated from the heap.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The caller must eventually call close_capsule_file() to unmap and
//		.. reliquish the mapped memory, as well as disposing of the capsule
//		.. info record.
//
// (2)	Our use of open/mmap follows an example from chapter four of "Linux
//		System Programming: Talking Directly to the Kernel and C Library", by
//		Robert Love, as presented at
//			www.devshed.com/c/a/BrainDump/Using-mmap-for-Advanced-File-IO
//
//----------

capinfo* open_capsule_file
   (char*		filename)
	{
	int			fdes;
	struct stat	sb;
	void*		mappedData;
	size_t		dataSize;
	u64			magic, fileSize;
	u32			magicA, magicB;
	int			swap64halves, littleEndian;
	capinfo*	cap;

	// open the file

	fdes = open (filename, O_RDONLY);
	if (fdes < 0)
		goto open_failed;

	if (fstat (fdes, &sb) == -1)
		goto fstat_failed;

	if (!S_ISREG (sb.st_mode))
		goto non_regular_file;

	// map the pre-header and read the file size;  note that we have to check
	// the magic number to figure out how to descramble the file size

	mappedData = mmap (0, capsulePreHeaderSize, PROT_READ, MAP_SHARED, fdes, 0);
	if ((mappedData == NULL) || (mappedData == MAP_FAILED))
		{ fileSize = capsulePreHeaderSize;  goto mmap_failed; }

	magic    = ((u64*) mappedData)[0];
	fileSize = ((u64*) mappedData)[1];
	munmap (mappedData, capsulePreHeaderSize);

	swap64halves = littleEndian = false;
	magicA = (u32) (magic >> 32);
	magicB = (u32)  magic;

	if (((magicA == refcapMagicABig)    && (magicB == refcapMagicBBig))
	 || ((magicA == refcapMagicALittle) && (magicB == refcapMagicBLittle)))
		; // ok, and no half swapping needed
	else if (((magicA == refcapMagicBBig)    && (magicB == refcapMagicABig))
	      || ((magicA == refcapMagicBLittle) && (magicB == refcapMagicALittle)))
		{ // ok, but half swapping needed
		magic    = swap_64_halves (magic);
		fileSize = swap_64_halves (fileSize);
		magicA = (u32) (magic >> 32);
		magicB = (u32)  magic;
		swap64halves = true;
		}
	else
		goto bad_magic;

	if (magicA == refcapMagicALittle)
		{
		fileSize = swap_two32_endian (fileSize);
		littleEndian = true;
		}

	dataSize = fileSize;
	if (dataSize != fileSize)
		goto file_size_overlow;

	// now that we know the size, map the whole thing;  note that there is
	// no point in memory-mapping a file on an architecture with different
	// endianness that the one the file was created on

	if ((littleEndian) || (swap64halves))
		goto architecture_mismatch;

	mappedData = mmap (0, dataSize, PROT_READ, MAP_SHARED, fdes, 0);
	if ((mappedData == NULL) || (mappedData == MAP_FAILED))
		goto mmap_failed;

	capsule_set_stat (sharedAddress, mappedData);

	// close the file;  suprisingly, we can close it even though it is mapped

	if (close (fdes) == -1)
		goto close_failed;

	// success!

	cap = (capinfo*) zalloc_or_die ("open_capsule_file", sizeof(capinfo));

	cap->dataSize     = dataSize;
	cap->mappedData   = mappedData;
	cap->swap64halves = swap64halves;
	cap->littleEndian = littleEndian;

	return cap;

	//////////
	// failure exits
	//////////

open_failed:
	suicidef_with_perror ("open(%s) failed (returned file descriptor = %d)",
	                      filename, fdes);
	return NULL; // (never gets here)

close_failed:
	suicidef_with_perror ("close() for %s failed", filename);
	return NULL; // (never gets here)

fstat_failed:
	suicidef_with_perror ("fstat() for %s failed", filename);
	return NULL; // (never gets here)

non_regular_file:
	suicidef ("%s is not a regular file, so is not mappable as a capsule file",
	          filename);
	return NULL; // (never gets here)

mmap_failed:
	suicidef_with_perror ("mmap() for %s failed (attempted 0x%s bytes, return value = %p)",
	                      filename, hex_64_string(fileSize), mappedData);
	return NULL; // (never gets here)

bad_magic:
	suicidef ("%s is not a capsule file (magic = 0x%s)",
	          filename, hex_64_string(magic));
	return NULL; // (never gets here)

file_size_overlow:
	suicidef ("file size overflow in %s (fileSize = 0x%s, dataSize = 0x%s)",
	          filename, hex_64_string(fileSize), hex_64_string(dataSize));
	return NULL; // (never gets here)

#define suggestions " rebuild it using --writecapsule"

architecture_mismatch:
	if ((littleEndian) && (!swap64halves))
		suicidef ("architecture mismatch for %s (8-byte words have halves swapped);"
		          suggestions,
		          filename);
	else if ((!littleEndian) && (swap64halves))
		suicidef ("architecture mismatch for %s (4-byte words are wrong endian);"
		          suggestions,
		          filename);
	else // if ((littleEndian) && (swap64halves))
		suicidef ("architecture mismatch for %s (8-byte words are wrong endian);"
		          suggestions,
		          filename);
	return NULL; // (never gets here)
	}

//----------
//
// close_capsule_file--
//	Close (== unmap) a Target Sequence Capsule File.
//
//----------
//
// Arguments:
//	capinfo*	cap:	The capsule info record, as returned by
//						.. open_capsule_file.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void close_capsule_file
   (capinfo*	cap)
	{
	if (cap == NULL) return;
	if (cap->mappedData != NULL) munmap (cap->mappedData, cap->dataSize);
	free_if_valid ("close_capsule_file", cap);
	}

//----------
//
// locate_capsule_data--
//	Loacte a data block in a mapped Target Sequence Capsule.
//
//----------
//
// Arguments:
//	capinfo*	cap:		The capsule info record.
//	u32			dataType:	The desired block's data type code (one of cap_xxx).
//	u32*		blockInfo:	Place to return the block info word.  This can be
//							.. NULL.
//	u32*		blockSize:	Place to return the block's size.  This can be NULL.
//
// Returns:
//	A pointer to the mapped block (NULL if the block is not found).
//
//----------

void* locate_capsule_data
   (capinfo*	cap,
	u32			blockType,
	u32*		blockInfo,
	u64*		blockSize)
	{
	char*		scan;
	u32			headerLength, numEntries, ix;
	u32			dataTypeCode;
	u64			blockOffset;

	scan = ((char*) cap->mappedData) + capsulePreHeaderSize;

	headerLength = *((u32*) scan);  scan += sizeof(u32);
	if ((headerLength % capsuleHeaderEntrySize) != 8)
		suicidef ("bad capsule header (length = %08X)", headerLength);
	numEntries = (headerLength - 8) / capsuleHeaderEntrySize;

	for (ix=0 ; ix<numEntries ; ix++)
		{
		dataTypeCode = *((u32*) scan);  scan += sizeof(u32);
		if (dataTypeCode == cap_terminator)
			suicide ("bad capsule header (premature terminator)");
		if (dataTypeCode != blockType)
			{ scan += capsuleHeaderEntrySize - sizeof(u32);  continue; }

		// block found

		if (blockInfo != NULL) *blockInfo = *((u32*) scan);
		scan += sizeof(u32);

		blockOffset = *((u64*) scan);  scan += sizeof(u64);

		if (blockSize != NULL) *blockSize = *((u64*) scan);

		return (void*) (((char*) cap->mappedData) + blockOffset);
		}

	// block not found

	return NULL;
	}

//----------
//
// capsule_zero_stats--
//	Clear the statistics for this module.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

void capsule_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&refSharingStats, 0, sizeof(refSharingStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// capsule_show_stats:
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to print the stats to.
//
// Returns:
//	(nothing)
//
//----------

void capsule_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats

	if (f == NULL) return;
	if (!refSharingStats.display) return;

	fprintf (f, "     header length: %s\n", hex_64_string(refSharingStats.headerLength));
	fprintf (f, "      header bytes: %s\n", hex_64_string(refSharingStats.headerBytes));
	fprintf (f, "       name offset: %s\n", hex_64_string(refSharingStats.nameOffset));
	fprintf (f, "       name length: %s\n", hex_64_string(refSharingStats.nameLength));
	fprintf (f, "        name bytes: %s\n", hex_64_string(refSharingStats.nameBytes));
	fprintf (f, "       nucs offset: %s\n", hex_64_string(refSharingStats.nucsOffset));
	fprintf (f, "       nucs length: %s\n", hex_64_string(refSharingStats.nucsLength));
	fprintf (f, "        nucs bytes: %s\n", hex_64_string(refSharingStats.nucsBytes));
	fprintf (f, "       rvrs offset: %s\n", hex_64_string(refSharingStats.rvrsOffset));
	fprintf (f, "       rvrs length: %s\n", hex_64_string(refSharingStats.rvrsLength));
	fprintf (f, "        rvrs bytes: %s\n", hex_64_string(refSharingStats.rvrsBytes));
	fprintf (f, "       bits offset: %s\n", hex_64_string(refSharingStats.bitsOffset));
	fprintf (f, "       bits length: %s\n", hex_64_string(refSharingStats.bitsLength));
	fprintf (f, "        bits bytes: %s\n", hex_64_string(refSharingStats.bitsBytes));
	fprintf (f, "       last offset: %s\n", hex_64_string(refSharingStats.lastOffset));
	fprintf (f, "       last length: %s\n", hex_64_string(refSharingStats.lastLength));
	fprintf (f, "        last bytes: %s\n", hex_64_string(refSharingStats.lastBytes));
	fprintf (f, "       prev offset: %s\n", hex_64_string(refSharingStats.prevOffset));
	fprintf (f, "       prev length: %s\n", hex_64_string(refSharingStats.prevLength));
	fprintf (f, "        prev bytes: %s\n", hex_64_string(refSharingStats.prevBytes));
	fprintf (f, "       info offset: %s\n", hex_64_string(refSharingStats.infoOffset));
	fprintf (f, "       info length: %s\n", hex_64_string(refSharingStats.infoLength));
	fprintf (f, "        info bytes: %s\n", hex_64_string(refSharingStats.infoBytes));
	fprintf (f, "       part offset: %s\n", hex_64_string(refSharingStats.partOffset));
	fprintf (f, "       part length: %s\n", hex_64_string(refSharingStats.partLength));
	fprintf (f, "        part bytes: %s\n", hex_64_string(refSharingStats.partBytes));
	fprintf (f, "       pool offset: %s\n", hex_64_string(refSharingStats.poolOffset));
	fprintf (f, "       pool length: %s\n", hex_64_string(refSharingStats.poolLength));
	fprintf (f, "        pool bytes: %s\n", hex_64_string(refSharingStats.poolBytes));
	fprintf (f, "       seed offset: %s\n", hex_64_string(refSharingStats.seedOffset));
	fprintf (f, "       seed length: %s\n", hex_64_string(refSharingStats.seedLength));
	fprintf (f, "        seed bytes: %s\n", hex_64_string(refSharingStats.seedBytes));
	fprintf (f, "        end offset: %s\n", hex_64_string(refSharingStats.endOffset));
	fprintf (f, "    shared address: %p\n", refSharingStats.sharedAddress);
	fprintf (f, "-------------------\n");

#endif // collect_stats
	}
