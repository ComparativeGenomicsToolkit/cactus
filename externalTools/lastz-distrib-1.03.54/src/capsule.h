//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: capsule.h
//
//----------

#ifndef capsule_H				// (prevent multiple inclusion)
#define capsule_H

// other files

#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff

// establish ownership of global variables

#ifdef capsule_owner
#define global
#else
#define global extern
#endif

//----------
//
// Reference Sequence Capsule File:
//
// A capsule file encapsulates the information needed for a reference sequence,
// including its seed word index table.
//
//----------
//
// Reference capsule file format:
//
//	Fields can be in big- or little-endian format;  they must match the
//	endianess of the magic number.  Similarly, 64-bit fields must match the
//	order indicated by the 64-bit magic number.  If the data blocks containing
//	multi-byte numeric fields are directly mapped by the program (e.g. with
//	mmap), the program will reject a file fi the magic number doesn't match
//	what is expected on the running-on platfrom.
//
//	Version 1
//
//	offset 0x00:   DA C8 9D 8E   big endian magic number (8E 9D C8 DA => little endian)
//	offset 0x04:   60 11 EF 1B   big endian magic number (1B EF 11 60 => little endian)
//	offset 0x08:   xx xx xx xx   file size (one half) (in bytes)
//	offset 0x0C:   xx xx xx xx   file size (other half)
//	offset 0x10:   00 00 01 00   version 1.0 (fourth byte is sub version)
//	offset 0x14:   xx xx xx xx   header length H (in bytes, including this field)
//	offset 0x18:    ...          header
//	 ...
//	offset H+0x14: 68 45 6E 64   ('hEnd') header terminator
//	offset ...:     ...          other data, pointed to by header entries
//
//	Each header entry (at offset E) consists of the following fields:
//
//	offset E+0x00: xx xx xx xx   data type code (one of cap_xxx below)
//	offset E+0x04: xx xx xx xx   extra info, per type code
//	offset E+0x08: xx xx xx xx   data offset (one half)
//	offset E+0x0C: xx xx xx xx   data offset (other half)
//	offset E+0x10: xx xx xx xx   data length (one half)
//	offset E+0x14: xx xx xx xx   data length (other half)
//
//	Type codes indicate the type of data block.  Valid types are shown here.
//	and are usually four ascii characters.
//
//	'name'	The name of the reference sequence.
//
//		offset 0x00:  ...			zero-terminated ascii string
//
//	'nucs'	The reference sequence nucleotides, one byte per nucleotide.  This
//			.. corresponds to seq.v.
//
//		offset 0x00:  ...           zero-terminated ascii string of length L+1,
//		                            .. where L is the length of the sequence
//
//	'rvrs'	The reference sequence nucleotides in reverse order (*not*
//			complemented).  This corresponds to targetRev in lastz's main().
//
//		offset 0x00:  ...           zero-terminated ascii string of length L+1.
//
//	'bits'	The reference sequence nucleotides in forward order, encoded as
//			bits.  This corresponds to postable.asBits
//
//		offset 0x00:  ...           bit data of length L/4 (rounded up to a
//									.. multiple of 16 bytes).
//
//	'last'	The last position, in the reference, of each seed word.  This
//			.. corresponds to postable.last.  It has has 4^W entries with each
//			.. entry being an index into the reference nucleotides.
//
//		offset 0x00: xx xx xx xx    postable.last[0]
//		offset 0x04: xx xx xx xx    postable.last[1]
//		offset 0x08:  ...
//
//	'prev'	The previous position, in the reference, of each seed word.  This
//			.. corresponds to postable.prev.  It has has ceil(L/Z) entries with
//			.. entry being an 'adjusted' index into the reference nucleotides.
//
//		offset 0x00: xx xx xx xx    postable.prev[0]
//		offset 0x04: xx xx xx xx    postable.prev[1]
//		offset 0x08:  ...
//
//	'info'	additional information about the reference sequence.
//
//		offset 0x00:    xx xx xx xx   seq.start
//		offset 0x04:    xx xx xx xx   seq.trueLen
//		offset 0x08:    xx xx xx xx   seq.revCompFlags
//		offset 0x0C:    xx xx xx xx   seq.contig
//		offset 0x10:    xx xx xx xx   seqpartition.len
//
//	'part'	sequence partitions
//
//		offset 0x00: xx xx xx xx    seqpartition.p[0].sepBefore
//		offset 0x04: xx xx xx xx    seqpartition.p[0].sepAfter
//		offset 0x08: xx xx xx xx    seqpartition.p[0].contig
//		offset 0x0C: xx xx xx xx    seqpartition.p[0].startLoc
//		offset 0x10: xx xx xx xx    seqpartition.p[0].trueLen
//		offset 0x14: xx xx xx xx    seqpartition.p[0].header
//		offset 0x18: xx xx xx xx    seqpartition.p[1].sepBefore
//		offset 0x1C: xx xx xx xx    seqpartition.p[1].sepAfter
//		offset 0x20: xx xx xx xx    seqpartition.p[1].contig
//		offset 0x24: xx xx xx xx    seqpartition.p[1].startLoc
//		offset 0x28: xx xx xx xx    seqpartition.p[1].trueLen
//		offset 0x2C: xx xx xx xx    seqpartition.p[1].header
//		offset 0x30: xx xx xx xx    seqpartition.p[2].sepBefore
//		offset 0x34:  ...
//
//	'pNam'	names for sequence partitions
//
//		offset 0x00: xx xx xx xx    seqpartition.pool[0] .. seqpartition.pool[3]
//		offset 0x04: xx xx xx xx    seqpartition.pool[4] .. seqpartition.pool[7]
//		offset 0x08:  ...
//
//	'seed'	The seed pattern used to build the seed word position table.  This
//			.. (mostly) corresponds to struct seed, but is not directly mapped.
//			.. Instead an instance of a struct seed is built from these values.
//			.. (seed capseed struct)
//
//		offset 0x00:    xx xx xx xx   postable.step
//		offset 0x04:    xx xx xx xx   seed.type
//		offset 0x08:    xx xx xx xx   seed.length
//		offset 0x0C:    xx xx xx xx   seed.weight
//		offset 0x10:    xx xx xx xx   seed.resolvingMask
//		offset 0x14:    xx xx xx xx   seed.revComp
//		offset 0x18:    xx xx xx xx   seed.isHalfweight
//		offset 0x1C:    xx xx xx xx   (P) seed.numParts
//		offset 0x20:    xx xx xx xx   seed.shift[0]
//		offset 0x24:     ..
//		offset 0x1C+4P: xx xx xx xx   seed.shift[P-1]
//		offset 0x20+4P: xx xx xx xx   seed.mask[0]
//		offset 0x24+4P:  ..
//		offset 0x1C+8P: xx xx xx xx   seed.mask[P-1]
//		offset 0x20+8P: xx xx xx xx   seed.transFlips[0]
//		offset 0x24+8P:  ...
//		 ...
//		offset ...:     00 00 00 00
//
//----------

// capsule file magic number(s)

static const u32 refcapMagicABig    = 0xDAC89D8E;	// in big endian format
static const u32 refcapMagicALittle = 0x8E9DC8DA;	// in little endian format

static const u32 refcapMagicBBig    = 0x6011EF1B;	// in big endian format
static const u32 refcapMagicBLittle = 0x1BEF1160;	// in little endian format

static const u32 refcapVersion      = 0x00000100;	// version 1.0

// miscellaneous capsule file block sizes

#define capsulePreHeaderSize   0x14
#define capsuleHeaderEntrySize 0x18

// capsule file data type codes
// $$$ add a creator field

#define cap_seqName        0x6E616D65	// 'name'
#define cap_seqForward     0x6E756373	// 'nucs'
#define cap_seqReverse     0x72767273	// 'rvrs'
#define cap_seqBits        0x62697473	// 'bits'
#define cap_lastPosTable   0x6C617374	// 'last'
#define cap_prevPosTable   0x70726576	// 'prev'
#define cap_seqInfo        0x696E666F	// 'info'
#define cap_seed           0x73656564	// 'seed'
#define cap_partitions     0x70617274	// 'part'
#define cap_partitionNames 0x704E616D	// 'pNam'
#define cap_terminator     0x68456E64	// 'hEnd'

typedef struct capinfo
	{
	void*	mappedData;			// Pointer to the file's mapped data.
	size_t	dataSize;			// Size of the mapped data.
	int		swap64halves;		// true => perform swap_64_halves on 64-bit
								//         .. values
	int		littleEndian;		// true => perform swap_32_endian on 32-bit
								//         .. values
	} capinfo;

typedef struct capseqinfo		// (corresponds to cap_info block)
	{
	u32		startLoc;			// seq.startLoc
	u32		trueLen;			// seq.trueLen
	u32		revCompFlags;		// seq.revCompFlags
	u32		contig;				// seq.contig
	u32		numPartitions;		// seqpartition.len
	} capseqinfo;

typedef struct cappartition		// (corresponds to cap_partitions block)
	{							//  .. layout must match struct partition (sequences.h)
	u32		sepBefore;			// partition.sepBefore
	u32		sepAfter;			// partition.sepAfter
	u32		contig;				// partition.contig
	u32		startLoc;			// partition.startLoc
	u32		trueLen;			// partition.trueLen
	u32		header;				// partition.header
	} cappartition;

typedef struct cappartitionold	// (corresponds to old cap_partitions block,
	{							//  .. lacking sepAfter and startLoc fields)
	u32		sepBefore;			// partition.sepBefore
	u32		contig;				// partition.contig
	u32		trueLen;			// partition.trueLen
	u32		header;				// partition.header
	} cappartitionold;

typedef struct capseed			// (corresponds to cap_seed block)
	{
	u32		step;				// postable.step
	u32		type;				// seed.type
	u32		length;				// seed.length
	u32		weight;				// seed.weight
	u32		resolvingMask;		// seed.resolvingMask
	u32		revComp;			// seed.revComp
	u32		isHalfweight;		// seed.isHalfweight
	u32		numParts;			// seed.numParts
	u32		shift0;				// first entry of seed.shift[]
	} capseed;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int   display;
	u64   headerLength;
	u64   headerBytes;
	u64   nameOffset;
	u64   nameLength;
	u64   nameBytes;
	u64   nucsOffset;
	u64   nucsLength;
	u64   nucsBytes;
	u64   rvrsOffset;
	u64   rvrsLength;
	u64   rvrsBytes;
	u64   bitsOffset;
	u64   bitsLength;
	u64   bitsBytes;
	u64   lastOffset;
	u64   lastLength;
	u64   lastBytes;
	u64   prevOffset;
	u64   prevLength;
	u64   prevBytes;
	u64   infoOffset;
	u64   infoLength;
	u64   infoBytes;
	u64   partOffset;
	u64   partLength;
	u64   partBytes;
	u64   poolOffset;
	u64   poolLength;
	u64   poolBytes;
	u64   seedOffset;
	u64   seedLength;
	u64   seedBytes;
	u64   endOffset;
	void* sharedAddress;
	} refSharingStats;

// stats macros

#define capsule_count_stat(field)   ++refSharingStats.field
#define capsule_uncount_stat(field) --refSharingStats.field
#define capsule_set_stat(field,val) (refSharingStats.field = val)
#define capsule_copy_stat(field)    (refSharingStats.field = field)
#define capsule_add_stat(field,val) (refSharingStats.field += val)
#else
#define capsule_count_stat(field)
#define capsule_uncount_stat(field)
#define capsule_set_stat(field,val)
#define capsule_copy_stat(field)
#define capsule_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void capsule_zero_stats (void);
void capsule_show_stats (FILE* f);

//----------
//
// prototypes for routines in capsule.c
//
//----------

u64      write_capsule_file  (FILE* f, char* filename, seq* seq, u8* revNucs,
                              postable* pt, seed* seed);
capinfo* open_capsule_file   (char* filename);
void     close_capsule_file  (capinfo* cap);
void*    locate_capsule_data (capinfo* cap, u32 blockType,
                              u32* blockInfo, u64* blockSize);

#undef global
#endif // capsule_H
