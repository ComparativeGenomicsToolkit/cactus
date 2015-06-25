//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sequences.c
//
//----------
//
// sequences--
//	Support for DNA sequences.
//
// FASTA format stores DNA sequences as plain text.  A single file can store
// multiple sequences, which we call contigs.  Each contig begins with a header
// line, which begins with a ">".
//
// FASTQ format stores DNA sequences and associated base-call quality scores as
// plain text.  A single file can store multiple sequences (contigs).  Each
// contig consists of four lines-- a header line (which begins with a "@"), a
// sequence line (nucleotides), a separator line (which begins with a "+"), and
// a qualities line.  As of Apr/2011, a spec for FASTQ files can be found at
//	http://maq.sourceforge.net/fastq.shtml
// Note that there has been much confusion with FASTQ quality values, since
// different variants of the format encode them in different ways.  As of
// Apr/2011, the only thing we don't interpret quality values, so we are immune
// to that confusion.
//
// NIB format stores a single DNA sequence, containing {A,C,G,T,a,c,g,t,N,N} in
// two bases per byte.  As of Jan/2008, a spec for NIB files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format8
//
// 2BIT format stores multiple DNA sequences encoded as four bases per byte
// with some additional information describing runs of masked bases or Ns.  As
// of Jan/2008, a spec for 2BIT files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format7
//
// HSX format is a hashed sequence index that consists of references to
// sequences in other files.  This allows random access.  See the file format
// spec in the lastz readme file for more information.
//
// Quantum-dna format stores each base as a byte, with the meaning of the
// byte values essentially defined by the scoring matrix.  See the file format
// spec in the lastz readme file for more information.
//
//----------
//
// WARNING:	As of this writing (Apr/2008), the code to read sequences is a
//			mess.  The additions of rewindability and contigs-of-interest did
//			not fit well with the original design, and have been patched in
//			with bandaids.  The author hopes to rectify this in the future.
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <ctype.h>				// standard C upper/lower stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff

#define  sequences_owner		// (make this the owner of its globals)
#include "sequences.h"			// interface to this module


#define pathSlash '/'
#ifdef compileForWindows
#undef pathSlash
#define pathSlash '\\'
#endif

//----------
//
// debug set ups
//
//----------

#define sequence_filename(seq)  (((seq)->filename != NULL)? (seq)->filename : "(unnamed sequence file)")

//--- debug set up for debugging the contigs-of-interest names file ---
//    .. and also chores file

//#define debugNamesFile

#ifndef debugNamesFile
#define debugNamesFile_1 ;
#define debugNamesFile_2 ;
#define debugNamesFile_3 ;
#define debugNamesFile_4 ;
#define debugNamesFile_5 ;
#define debugNamesFile_6 ;
#define debugNamesFile_7 ;
#define debugNamesFile_8 ;
#define debugNamesFile_9 ;
#define debugNamesFile_10 ;
#define debugNamesFile_11 ;
#define debugNamesFile_12 ;
#define debugNamesFile_13 ;
#define debugNamesFile_14 ;
#define debugNamesFile_15 ;
#define debugNamesFile_16 ;
#define debugNamesFile_17 ;
#define debugNamesFile_18 ;
#define debugNamesFile_19 ;
#define debugNamesFile_20 ;
#define debugNamesFile_21 ;
#define debugNamesFile_22 ;
#define debugNamesFile_23 ;
#endif // not debugNamesFile

#ifdef debugNamesFile

#define debugNamesFile_1                                                     \
	fprintf (stderr, "load_sequence (%s):\n", sequence_filename(_seq));      \
	if ((_seq != NULL) && (_seq->preLoaded))                                 \
		fprintf (stderr, "  preloaded, header: %s\n", _seq->header);

#define debugNamesFile_2                                                     \
	fprintf (stderr, "  (back in load_sequence)\n");

#define debugNamesFile_3                                                     \
	fprintf (stderr, "load_sequence_core (%s):\n", sequence_filename(_seq)); \
	if (!keeper)                                                             \
		fprintf (stderr, "  (not a keeper)\n");

#define debugNamesFile_4                                                     \
	fprintf (stderr, "  header: %s\n", _seq->header);

#define debugNamesFile_5                                                    \
	fprintf (stderr, "locate_hsx_first_sequence:"                            \
	                 " single contig-of-interest"                            \
	                 " hsx.contigFilePos <- %010lX\n",                       \
	                 (long unsigned int) _seq->hsx.contigFilePos);

#define debugNamesFile_6                                                     \
	fprintf (stderr, "locate_hsx_first_sequence:"                            \
	                 " no contig names"                                      \
	                 " hsx.contigFilePos <- %010lX\n",                       \
	                 (long unsigned int) _seq->hsx.contigFilePos);

#define debugNamesFile_7                                                     \
	fprintf (stderr, "load_hsx_sequence:"                                    \
	                 " hsx.contigFilePos == %010lX\n",                       \
	                 (long unsigned int) _seq->hsx.contigFilePos);

#define debugNamesFile_8                                                     \
	fprintf (stderr, "load_hsx_sequence:"                                    \
	                 " hsx.contigFilePos <- %010lX\n",                       \
	                 (long unsigned int) _seq->hsx.contigFilePos);

#define debugNamesFile_9                                                     \
	fprintf (stderr, "lookup_hsx_sequence(%s):\n", name);                    \
	fprintf (stderr, "  bucket      == %08X\n  hash entry  == %010lX\n",     \
	                 (unsigned int) bucket, (long unsigned int) fileOffset);

#define debugNamesFile_10                                                    \
	fprintf (stderr, "  bucketStart == %010lX\n",                            \
	                 (long unsigned int) bucketStart);

#define debugNamesFile_11                                                    \
	fprintf (stderr, "another_sequence (%s):\n", sequence_filename(_seq));

#define debugNamesFile_12                                                    \
	long int	fpos;                                                        \
	fprintf (stderr, "find_next_general_fasta_coi (%s):\n", sequence_filename(_seq));

#define debugNamesFile_13                                                    \
	long int	fpos;                                                        \
	fprintf (stderr, "find_next_fastq_coi (%s):\n", sequence_filename(_seq));

#define debugNamesFile_14                                                    \
	fpos = ftell (_seq->f) - _seq->pendingLen - 1;

#define debugNamesFile_15                                                    \
	{                                                                        \
	char headerSaveCh = header[headerLen];                                   \
	header[headerLen] = 0;                                                   \
	fprintf (stderr, "  test:  [%08lX] %s\n", fpos, header);                 \
	header[headerLen] = headerSaveCh;                                        \
	}

#define debugNamesFile_16                                                    \
	fprintf (stderr, "  found: [%08lX] %s\n", fpos, _seq->nextContigName);

#define debugNamesFile_17                                                    \
	fprintf (stderr, "find_next_2bit_coi (%s): [%s]\n",                      \
	                 sequence_filename(_seq), _seq->nextContigName);

#define debugNamesFile_18                                                    \
	fprintf (stderr, "find_next_hsx_coi (%s): [%s]\n",                       \
	                 sequence_filename(_seq), _seq->nextContigName);

#define debugNamesFile_19                                                    \
	fprintf (stderr, "load_hsx_sequence:"                                    \
	                 " hsx.contigFilePos <- %010lX\n",                       \
	                 (long unsigned int) _seq->hsx.contigFilePos);

#define debugNamesFile_20                                                    \
	fprintf (stderr, "read_contig_name: %s\n", line);

#define debugNamesFile_21                                                    \
	static int choreNumber = 0;

#define debugNamesFile_22                                                    \
	choreNumber++;                                                           \
	fprintf (stderr, "read_chore #%d: %s\n", choreNumber, line);

#define debugNamesFile_23                                                    \
	fprintf (stderr, "read_chore -->");                                      \
                                                                             \
	if (_seq->chore.tSubrange)                                               \
		fprintf (stderr, " %s " unsposFmt " " unsposFmt,                     \
		                 _seq->chore.tName,                                  \
		                 _seq->chore.tStart, _seq->chore.tEnd);              \
	else                                                                     \
		fprintf (stderr, " %s * *",                                          \
		                 _seq->chore.tName);                                 \
                                                                             \
	if (_seq->chore.qSubrange)                                               \
		fprintf (stderr, " %s " unsposFmt " " unsposFmt,                     \
		                 _seq->nextContigName,                               \
		                 _seq->chore.qStart, _seq->chore.qEnd);              \
	else                                                                     \
		fprintf (stderr, " %s * *",                                          \
		                 _seq->nextContigName);                              \
                                                                             \
	if      (_seq->chore.qStrand == 0) fprintf (stderr, " +");               \
	else if (_seq->chore.qStrand <  0) fprintf (stderr, " -");               \
                                                                             \
	if      (_seq->chore.idTag[0] != 0)                                      \
		fprintf (stderr, " id=%s", _seq->chore.idTag);                       \
                                                                             \
	if ((header != NULL) && (strcmp (header, _seq->nextContigName) == 0))    \
		fprintf (stderr, " (same as existing header)");                      \
                                                                             \
	fprintf (stderr, "\n");

#endif // debugNamesFile


//--- debug set up for debugging partitioned sequences ---

//#define debugPartitions

#ifndef debugPartitions
#define debugPartitions_1a ;
#define debugPartitions_1b ;
#define debugPartitions_2 ;
#define debugPartitions_3 ;
#define debugPartitions_4 ;
#endif // not debugPartitions

#ifdef debugPartitions

#define debugPartitions_1a                                                  \
	if (sp->p != NULL)                                                      \
		{                                                                   \
		print_partition_table (stderr, _seq);                               \
		print_sequence (stderr, _seq, "", 100);                             \
		}

#define debugPartitions_1b                                                  \
	if (sp->p != NULL)                                                      \
		{                                                                   \
		print_partition_table (stderr, _seq);                               \
		}

#define debugPartitions_2                                                   \
	if (_seq->filename == NULL)                                             \
		fprintf (stderr, "lookup_partition(.," unsposFmt ")\n",             \
		                 pos);                                              \
	else                                                                    \
		fprintf (stderr, "lookup_partition(%s," unsposFmt ")\n",            \
		                 _seq->filename, pos);                              \
	fprintf (stderr, "   p[%d].sepBefore=" unsposFmt                        \
	                   " p[%d].sepBefore=" unsposFmt "\n",                  \
	                 lo, p[lo].sepBefore, hi, p[hi].sepBefore);

#define debugPartitions_3                                                   \
	fprintf (stderr, "   p[%d].sepBefore=" unsposFmt                        \
	                   " p[%d].sepBefore=" unsposFmt                        \
	                   " p[%d].sepBefore=" unsposFmt "\n",                  \
	                 lo, p[lo].sepBefore, ix, p[ix].sepBefore, hi, p[hi].sepBefore);

#define debugPartitions_4                                                   \
	fprintf (stderr, "   -> %u." unsposFmt "\n", ix, pos);

#endif // debugPartitions


//#define debugSeparation

#ifndef debugSeparation
#define debugSeparation_1 ;
#define debugSeparation_2 ;
#define debugSeparation_3 ;
#define debugSeparation_4 ;
#define debugSeparation_5 ;
#define debugSeparation_6 ;
#define debugSeparation_7 ;
#define debugSeparation_8 ;
#define debugSeparation_9 ;
#define debugSeparation_10 ;
#endif // not debugSeparation

#ifdef debugSeparation

#define debugSeparation_1                                                   \
	if (sp->p != NULL)                                                      \
		{                                                                   \
		print_partition_table (stderr, _seq);                               \
		print_sequence (stderr, _seq, "", 100);                             \
		}

#define debugSeparation_2                                                   \
	fprintf (stderr, "   incoming separator at " unsposFmt "\n",            \
	                 (unspos) ((scan-1)-_seq->v));

#define debugSeparation_3                                                   \
	fprintf (stderr, "   extraPieces=%d\n",                                 \
	                 extraPieces);

#define debugSeparation_4                                                   \
	fprintf (stderr, "   copying sentinel from partition %d"                \
	                   " to partition %d\n",                                \
	                 fromIx, toIx);

#define debugSeparation_5                                                   \
	fprintf (stderr, "   scanning partition %d,"                            \
	                   " from " unsposFmt                                   \
	                   " thru " unsposFmt                                   \
	                   " (to partition = %d)\n",                            \
	                 fromIx, sepPrefix+1, sepSuffix-1, toIx);

#define debugSeparation_6                                                   \
	fprintf (stderr, "   seq->v[" unsposFmt "]=%c\n",                       \
	                 scan-_seq->v, ch);

#define debugSeparation_7                                                   \
	fprintf (stderr, "   sepAfter=" unsposFmt                               \
	                   " (ch=%c)\n",                                        \
	                 sepAfter, ch);

#define debugSeparation_8                                                   \
	fprintf (stderr, "   copying from partition %d"                         \
	                   " to partition %d"                                   \
	                   " sepBefore=" unsposFmt                              \
	                   " sepAfter=" unsposFmt                               \
	                   " startLoc=" unsposFmt "\n",                         \
	                 fromIx, toIx, sepBefore, sepAfter, startLoc);

#define debugSeparation_9                                                   \
	fprintf (stderr, "   done scanning partitions,"                         \
	                   " (to partition = %d)\n",                            \
	                 toIx);

#define debugSeparation_10                                                  \
	print_partition_table (stderr, _seq);                                   \
	print_sequence (stderr, _seq, "", 100);


#endif // debugSeparation


//--- debug set up for debugging problems with reading binary files ---

//#define debugBinaryFile

#ifdef debugBinaryFile

static FILE* dbg_fopen_or_die (const char* name, const char* mode);
static FILE* dbg_fopen_or_die (const char* name, const char* mode)
	{
	FILE* f = fopen_or_die (name, mode);
	fprintf (stderr, "fopen (\"%s\", \"%s\") = %08X\n",
	                 name, mode, (u32)f);
	return f;
	}

static void dbg_rewind (FILE *f);
static void dbg_rewind (FILE *f)
	{
	fprintf (stderr, "rewind (%08X)\n", (u32)f);
	rewind (f);
	}

static int dbg_fseek (FILE *f, long int offset, int whence);
static int dbg_fseek (FILE *f, long int offset, int whence)
	{
	fprintf (stderr, "fseek (%08X, %08X, %d)\n",
	                 (u32)f, (u32)offset, whence);
	return fseek (f, offset, whence);
	}

static size_t dbg_fread (void *ptr, size_t size, size_t nmemb, FILE* f);
static size_t dbg_fread (void *ptr, size_t size, size_t nmemb, FILE* f)
	{
	fprintf (stderr, "fread (%08X, %08X, %08X, %08X)\n",
	                 (u32)ptr, (u32)size, (u32)nmemb, (u32)f);
	return fread (ptr, size, nmemb, f);
	}

#define fopen_or_die dbg_fopen_or_die
#define rewind       dbg_rewind
#define fseek        dbg_fseek
#define fread        dbg_fread

#endif // debugBinaryFile


//--- debug set up for debugging problems with reading text files ---

//#define debugTextFile

#ifndef debugTextFile
#define debugTextFile_1 ;
#define debugTextFile_2 ;
#define debugTextFile_3 ;
#define debugTextFile_4 ;
#define debugTextFile_5 ;
#define debugTextFile_6 ;
#define debugTextFile_7 ;
#endif // not debugTextFile

#ifdef debugTextFile

#define debugTextFile_1                                                    \
	fprintf (stderr, "(for %s) rewinding\n",                               \
	                 _seq->filename);

#define debugTextFile_2                                                    \
	fprintf (stderr, "(for %s) add_partition\n",                           \
	                 _seq->filename);

#define debugTextFile_3                                                    \
	fprintf (stderr, "(for %s) from file: %d --> %s\n",                    \
	                 _seq->filename, ch, char_to_description(ch));

#define debugTextFile_4                                                    \
	fprintf (stderr, "(for %s) from buff: %d --> %s\n",                    \
	                 _seq->filename, ch, char_to_description(ch));

#define debugTextFile_5                                                    \
	fprintf (stderr, "(for %s) to   buff: %d --> %s\n",                    \
	                 _seq->filename, ch, char_to_description(ch));

#define debugTextFile_6                                                    \
	fprintf (stderr, "(for %s) saving file state\n",                       \
	                 _seq->filename);

#define debugTextFile_7                                                    \
	fprintf (stderr, "(for %s) restoring file state\n",                    \
	                 _seq->filename);

#endif // debugTextFile


//--- debug set up for debugging problems with reading text files ---

//#define debugFastaFile

#ifndef debugFastaFile
#define debugFastaFile_1 ;
#define debugFastaFile_2 ;
#endif // not debugFastaFile

#ifdef debugFastaFile

#define debugFastaFile_1                                                    \
	if (keeper)                                                             \
		fprintf (stderr, "load_fasta_sequence(%s,keeper)\n",                \
		                 _seq->filename);                                   \
	else                                                                    \
		fprintf (stderr, "load_fasta_sequence(%s,NOT keeper)\n",            \
		                 _seq->filename);

#define debugFastaFile_2                                                    \
	fprintf (stderr, "  header = \"%s\"\n", _seq->header);

#endif // debugFastaFile


//--- debug set up for debugging sequence cloning ---

//#define debugCloning

#ifndef debugCloning
#define debugCloning_1 ;
#define debugCloning_2 ;
#endif // not debugCloning

#ifdef debugCloning

#define debugCloning_1                                                      \
	fprintf (stderr, "cloning this seq:\n-----------------\n");             \
	dump_sequence_state (stderr, _seq);                                     \
	fprintf (stderr, "\n");

#define debugCloning_2                                                      \
	fprintf (stderr, "newSeq:\n-------\n");                                 \
	dump_sequence_state (stderr, newSeq);                                   \
	fprintf (stderr, "\n");

#endif // debugCloning


//--- debug set up for debugging "fencing" used to facilitate chores ---

//#define debugFencing

#ifndef debugFencing
#define debugFencing_1 ;
#define debugFencing_2 ;
#define debugFencing_3 ;
#define debugFencing_4 ;
#define debugFencing_5 ;
#define debugFencing_6 ;
#define debugFencing_7 ;
#define debugFencing_8 ;
#endif // not debugFencing

#ifdef debugFencing

#define debugFencing_1                                                      \
	fprintf (stderr, "fence:  ");

#define debugFencing_2                                                      \
	fprintf (stderr, " left: [" unsposFmt "]=%02X",                         \
	                 _seq->leftFencePos, _seq->leftFenceCh);

#define debugFencing_3                                                      \
	fprintf (stderr, " right: [" unsposFmt "]=%02X",                        \
	                 _seq->rightFencePos, _seq->rightFenceCh);

#define debugFencing_4                                                      \
	fprintf (stderr, "\n");

#define debugFencing_5                                                      \
	fprintf (stderr, "unfence:");

#define debugFencing_6                                                      \
	fprintf (stderr, " left: [" unsposFmt "]=%02X",                         \
	                 _seq->leftFencePos, _seq->leftFenceCh);

#define debugFencing_7                                                      \
	fprintf (stderr, " right: [" unsposFmt "]=%02X",                        \
	                 _seq->rightFencePos, _seq->rightFenceCh);

#define debugFencing_8                                                      \
	fprintf (stderr, "\n");

#endif // debugFencing

//----------
//
// private global data relating to fasta and csfasta format
//
//----------

// tables to map 8-bit ascii character to fasta or csfasta character type
//   "nucleotide" characters are the A, C, G, T and N
//   "ambiguous"  characters are the remaining IUPAC 15-letter alphabet

enum
	{ _bad = 0, _whitespace, _newline, _nucleotide, _ambiguous, _color };

#define __ _bad
#define _w _whitespace
#define _l _newline
#define _n _nucleotide
#define _a _ambiguous
#define _c _color

static const u8 char_to_fasta_type[256] =
	{
	__,__,__,__,__,__,__,__,__,_w,_l,__,_w,_l,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	_w,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	_w,_w,_w,_w,_w,_w,_w,_w,_w,_w,__,__,__,__,__,__, // 3x (numbers)
	__,_n,_a,_n,_a,__,__,_n,_a,__,__,_a,__,_a,_n,__, // 4x (upper case)
	__,__,_a,_a,_n,__,_a,_a,_n,_a,__,__,__,__,__,__, // 5x (upper case)
	__,_n,_a,_n,_a,__,__,_n,_a,__,__,_a,__,_a,_n,__, // 6x (lower case)
	__,__,_a,_a,_n,__,_a,_a,_n,_a,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__  // Fx
	};

static const u8 char_to_csfasta_type[256] =
	{
	__,__,__,__,__,__,__,__,__,_w,_l,__,_w,_l,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	_w,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	_c,_c,_c,_c,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,_n,__,_n,__,__,__,_n,__,__,__,__,__,__,__,__, // 4x (upper case)
	__,__,__,__,_n,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 6x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__  // Fx
	};

#undef __
#undef _w
#undef _l
#undef _n
#undef _a
#undef _c

//----------
//
// private global data relating to nib format
//
//----------

// nib file magic number(s)

static const u32 nibMagicBig     = 0x6BE93D3A;	// in big endian format
static const u32 nibMagicLittle  = 0x3A3DE96B;	// in little endian format

// tables to map 4-bit nybbles from nib file to a character
//
// nybbles map as follows:
//		nybble:    0 1 2 3 4 5 6 7 8 9 A B C D E F
//		character: T C A G N X X X t c a g n x x x
// For (alleged) efficiency's sake, we use separate lookup tables for the left
// and right nybble mapping.

static const unsigned char nibTo1stChar[256] = 
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"tttttttttttttttt"
	"cccccccccccccccc"
	"aaaaaaaaaaaaaaaa"
	"gggggggggggggggg"
	"nnnnnnnnnnnnnnnn"
	"xxxxxxxxxxxxxxxx"
	"xxxxxxxxxxxxxxxx"
	"xxxxxxxxxxxxxxxx";

static const unsigned char nibTo2ndChar[256] = 
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx"
	"TCAGNXXXtcagnxxx";

static const unsigned char nibTo1stCharUnmasked[256] = 
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"TTTTTTTTTTTTTTTT"
	"CCCCCCCCCCCCCCCC"
	"AAAAAAAAAAAAAAAA"
	"GGGGGGGGGGGGGGGG"
	"NNNNNNNNNNNNNNNN"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX"
	"XXXXXXXXXXXXXXXX";

static const unsigned char nibTo2ndCharUnmasked[256] = 
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX"
	"TCAGNXXXTCAGNXXX";

//----------
//
// private global data relating to 2bit format
//
//----------

// 2bit file magic number(s)

static const u32 twobitMagicBig    = 0x1A412743;	// in big endian format
static const u32 twobitMagicLittle = 0x4327411A;	// in little endian format

static const char* twobitToChars[256] =
	{
	"TTTT","TTTC","TTTA","TTTG","TTCT","TTCC","TTCA","TTCG",
	"TTAT","TTAC","TTAA","TTAG","TTGT","TTGC","TTGA","TTGG",
	"TCTT","TCTC","TCTA","TCTG","TCCT","TCCC","TCCA","TCCG",
	"TCAT","TCAC","TCAA","TCAG","TCGT","TCGC","TCGA","TCGG",
	"TATT","TATC","TATA","TATG","TACT","TACC","TACA","TACG",
	"TAAT","TAAC","TAAA","TAAG","TAGT","TAGC","TAGA","TAGG",
	"TGTT","TGTC","TGTA","TGTG","TGCT","TGCC","TGCA","TGCG",
	"TGAT","TGAC","TGAA","TGAG","TGGT","TGGC","TGGA","TGGG",

	"CTTT","CTTC","CTTA","CTTG","CTCT","CTCC","CTCA","CTCG",
	"CTAT","CTAC","CTAA","CTAG","CTGT","CTGC","CTGA","CTGG",
	"CCTT","CCTC","CCTA","CCTG","CCCT","CCCC","CCCA","CCCG",
	"CCAT","CCAC","CCAA","CCAG","CCGT","CCGC","CCGA","CCGG",
	"CATT","CATC","CATA","CATG","CACT","CACC","CACA","CACG",
	"CAAT","CAAC","CAAA","CAAG","CAGT","CAGC","CAGA","CAGG",
	"CGTT","CGTC","CGTA","CGTG","CGCT","CGCC","CGCA","CGCG",
	"CGAT","CGAC","CGAA","CGAG","CGGT","CGGC","CGGA","CGGG",

	"ATTT","ATTC","ATTA","ATTG","ATCT","ATCC","ATCA","ATCG",
	"ATAT","ATAC","ATAA","ATAG","ATGT","ATGC","ATGA","ATGG",
	"ACTT","ACTC","ACTA","ACTG","ACCT","ACCC","ACCA","ACCG",
	"ACAT","ACAC","ACAA","ACAG","ACGT","ACGC","ACGA","ACGG",
	"AATT","AATC","AATA","AATG","AACT","AACC","AACA","AACG",
	"AAAT","AAAC","AAAA","AAAG","AAGT","AAGC","AAGA","AAGG",
	"AGTT","AGTC","AGTA","AGTG","AGCT","AGCC","AGCA","AGCG",
	"AGAT","AGAC","AGAA","AGAG","AGGT","AGGC","AGGA","AGGG",

	"GTTT","GTTC","GTTA","GTTG","GTCT","GTCC","GTCA","GTCG",
	"GTAT","GTAC","GTAA","GTAG","GTGT","GTGC","GTGA","GTGG",
	"GCTT","GCTC","GCTA","GCTG","GCCT","GCCC","GCCA","GCCG",
	"GCAT","GCAC","GCAA","GCAG","GCGT","GCGC","GCGA","GCGG",
	"GATT","GATC","GATA","GATG","GACT","GACC","GACA","GACG",
	"GAAT","GAAC","GAAA","GAAG","GAGT","GAGC","GAGA","GAGG",
	"GGTT","GGTC","GGTA","GGTG","GGCT","GGCC","GGCA","GGCG",
	"GGAT","GGAC","GGAA","GGAG","GGGT","GGGC","GGGA","GGGG"
	};

//----------
//
// private global data relating to hsx format
//
//----------

// hsx file magic number(s)

static const u32 hsxMagicBig    = 0xD2527095;	// in big endian format
static const u32 hsxMagicLittle = 0x957052D2;	// in little endian format

static const u64 hsxMsBit5      = ((u64) 0x80) << (4*8);
static const u64 hsxMaxFilePos  = (u64) ((long int) -1);

//----------
//
// private global data relating to quantum-dna format
//
//----------

// quantum-dna file magic number(s)

static const u32 qdnaMagicBig       = 0xC4B47197;	// in big endian format
static const u32 qdnaMagicLittle    = 0x9771B4C4;	// in little endian format
static const u32 oldQdnaMagicBig    = 0xF656659E;	// in big endian format
static const u32 oldQdnaMagicLittle = 0x9E6556F6;	// in little endian format

//----------
//
// prototypes for private functions
//
//----------

static seq*   alloc_sequence_record (char* id);
static void   skip_sequences        (seq* _seq, int skipCount);
static void   load_sequence_core    (seq* _seq, int keeper);
static void   load_fasta_sequence   (seq* _seq, int keeper);
static void   parse_fasta_header    (seq* _seq);
static unspos parse_fasta           (seq* _seq, int storeEm);
static void   load_fastq_sequence   (seq* _seq, int keeper);
static void   parse_fastq_header    (seq* _seq);
static unspos parse_fastq           (seq* _seq);
static int    fastq_skip_content    (seq* _seq);
static void   load_csfasta_sequence (seq* _seq, int keeper);
static void   parse_csfasta_header  (seq* _seq);
static unspos parse_csfasta         (seq* _seq, int storeEm);
static void   load_nib_sequence     (seq* _seq, int keeper);
static void   read_2bit_header      (seq* _seq);
static void   load_2bit_sequence    (seq* _seq, int keeper);
static void   read_hsx_header       (seq* _seq);
static void   locate_hsx_first_sequence (seq* _seq);
static void   load_hsx_sequence     (seq* _seq, int keeper);
static void   load_qdna_sequence    (seq* _seq, int keeper);
static int    another_sequence_core (seq* _seq);
static void   create_short_header   (seq* _seq);
static int    find_next_fasta_coi   (seq* _seq);
static int    find_next_fastq_coi   (seq* _seq);
static int    find_next_csfasta_coi (seq* _seq);
static int    find_next_general_fasta_coi (seq* _seq, int allowComments);
static int    find_next_2bit_coi    (seq* _seq);
static int    find_next_hsx_coi     (seq* _seq);
static int    read_contig_name      (seq* _seq);
static int    read_chore            (seq* _seq);
static void   shorten_header        (char* src, int nameParseType, int skipPath,
                                     char** dst, u32* dstSize);
static void   whitespace_to_under   (char* s, int sLen);
static void   expand_nickname       (char* src, u32 contigNumber,
                                     char** dst, u32* dstSize);
static void   separate_sequence     (seq* _seq, char sepCh);
static void   add_partition         (seq* _seq, unspos sepPos,
                                     unspos startLoc, unspos trueLenOffset);
static void   copy_partitions       (seq* seqTo, seq* seqFrom);
static void   enough_partitions     (seq* _seq, u32 numPartitions, u32 poolSize,
                                     int anticipate, int roundUp);
static void   parse_sequence_name   (const char* name,
                                     char** filename, char** nickname,
                                     char** contigOfInterest,
                                     char** namesFilename,
                                     char** choresFilename,
                                     int* subsampleK, int* subsampleN,
                                     char** softMaskFilename, int* softMaskComplement,
                                     char** xMaskFilename, int* xMaskComplement,
                                     char** nMaskFilename, int* nMaskComplement,
                                     int* nameParseType,
                                     char** nameTrigger,
                                     int* doRevCompFlags,
                                     int* doUnmask,
                                     int* doPartitioning, int* doJoin,
                                     char* separatorCh,
                                     int* useFullNames,
                                     int* fileType,
                                     int* isQuantum, char** qCodingFilename,
                                     unspos* _start, unspos* _end,
                                     int* endIsSoft);
static int    detect_file_type      (seq* _seq);
static u32    read_4                (seq* _seq, int asBigEndian);
static u32    read_4_big            (seq* _seq);
static u32    read_4_little         (seq* _seq);
static u64    read_5                (seq* _seq, int asBigEndian);
static u64    read_5_big            (seq* _seq);
static u64    read_5_little         (seq* _seq);
static u64    read_6                (seq* _seq, int asBigEndian);
static u64    read_6_big            (seq* _seq);
static u64    read_6_little         (seq* _seq);
static int    skip_seq_whitespace   (seq* _seq);
static int    seq_getc              (seq* _seq);
static void   seq_ungetc            (char ch, seq* _seq);
static int    skip_chars            (seq* _seq, u32 toSkip);
static int    test_rewindability    (seq* _seq);
static void   save_fstate           (seq* _seq);
static void   restore_fstate        (seq* _seq);

//----------
//
// open_sequence_file--
//	Open a sequence file for read operations.
// open_rewindable_sequence_file--
//	Open a sequence file for read operations, which may be rewound later.  Note
//	that if the actual file is not rewindable, but only contains one sequence,
//	we can still satisfy the caller's desire for rewindability.
//
//----------
//
// Arguments:
//	char*	name:			The name of the file that sequence data is to be
//							.. read from.  This can be NULL, which indicates
//							.. that data is to be read from stdin.  Note that
//							.. the name may have actions attached as a suffix.
// 	int		fileType:		The type of file, e.g. seq_type_nib.  This can be
//							.. seq_type_unknown if the caller would like the
//							.. type to be determined from the file' contents.
//	int		choresAllowed:	true => an "alignment chores" action is allowed.
//	char*	choresFilename:	The name of the file to read "alignment chores"
//							.. from.  This can be NULL, in which case there
//							.. are no *specific* chores to perform.  Note that,
//							.. even if this is NULL, a chores file may be
//							..  specified by an action attached to the file
//							.. name.
//	unspos	allocLen:		If non-zero, pre-allocate for a sequence of this
//							.. length.  Zero indicates that the caller doesn't
//							.. have any pre-allocation preference.
//	int		needTrueLen:	true  => set seq->trueLen correctly, even if this
//							         .. means reading additional characters
//							         .. outside the desired (sub)interval
//							false => the value written to trueLen is unimportant
//	int		allowAmbiDNA:	true  => permit ambiguous DNA characters
//							          .. B,D,H,K,M,R,S,V,W,Y
//							false => only A,C,G,T,N,X permitted
//	u8*		qToComplement:	(similar to nuc_to_complement) array to map a
//							.. quantum base to its complement.  This is only
//							.. used if the sequence is quantum DNA.  This may
//							.. be NULL (in which case the sequence can not be
//							.. reverse-complemented).
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------
//
// Implementation notes:
//
// To satisfy the caller's request that the file be rewindable, we perform the
// initial sequence load here (and set a flag to let load_sequence know this
// has happened).  Then we check whether the file contains any additional data.
// If it doesn't, then we treat the file as rewindable regardless of whether
// the underlying file actualy is.  Only if the file contains additional data
// do we then perform a test for rewindability, by attempting to set the file's
// position back to the end of the first sequence.
//
//----------

static seq* private_open_sequence_file (char* name, int fileType,
               int choresAllowed, char* choresFilename,
               unspos allocLen,
               int beRewindable, int needTrueLen, int allowAmbiDNA,
               u8* qToComplement);

seq* open_sequence_file (char* name, int fileType, int choresAllowed, char* choresFilename, unspos allocLen, int needTrueLen, int allowAmbiDNA, u8* qToComplement)
	{ return private_open_sequence_file (name, fileType, choresAllowed, choresFilename, allocLen, false, needTrueLen, allowAmbiDNA, qToComplement); }

seq* open_rewindable_sequence_file (char* name, int fileType, int choresAllowed, char* choresFilename, unspos allocLen, int needTrueLen, int allowAmbiDNA, u8* qToComplement)
	{ return private_open_sequence_file (name, fileType, choresAllowed, choresFilename, allocLen, true, needTrueLen, allowAmbiDNA, qToComplement); }

static seq* private_open_sequence_file
   (char*	name,
	int		fileType,
	int		choresAllowed,
	char*	choresFilename,
	unspos	allocLen,
	int		beRewindable,
	int		needTrueLen,
	int		allowAmbiDNA,
	u8*		qToComplement)
	{
	seq*	_seq;
	int		isQuantum = false;
	char*	qCodingFilename = NULL;
	int		forcedfileType = seq_type_unknown;
	char*	header;
	int		searchForContig;
	int		err;

	// allocate the sequence tracking structure

	_seq = alloc_sequence_record ("open_sequence");
	_seq->vOwner  = true;  				// (even though _seq->v is NULL, we
	_seq->vcOwner = true;				//  .. still will be considered as the
	_seq->vqOwner = true;				//  .. 'owner' so we can resize it)
	_seq->partition.poolOwner = true;	// (similarly, we are owner of partition
										//  .. names pool so we can resize it)
	_seq->headerOwner         = true;  	// (and we're owner of header[],
	_seq->shortHeaderOwner    = true;	//  .. shortHeader[] and trueHeader[]
	_seq->trueHeaderOwner     = true;	//  .. so we can resize them)

	_seq->pendingChars = zalloc_or_die ("open_sequence (pendingChars)",
                                        seqBufferSize);
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen   = 0;

	_seq->lockedHeader = false;
	_seq->needTrueLen  = needTrueLen;
	_seq->allowAmbiDNA = allowAmbiDNA;

	// if we have no name, use stdin

	if (name == NULL)
		{
		_seq->filename = copy_string ("(stdin)");
		_seq->f        = stdin;
		}

	// otherwise, open the file
	// nota bene:  we'd like to open the file as "rt" instead of "rb" if it is
	//             a fasta file;  unfortunately we don't know what the file
	//             type is until we open it

	else
		{
		if (choresFilename != NULL) _seq->choresFilename = copy_string (choresFilename);
		                       else _seq->choresFilename = NULL;

		parse_sequence_name (name,
		                     &_seq->filename, &_seq->header,
		                     &_seq->contigOfInterest,
		                     &_seq->namesFilename,
		                     &_seq->choresFilename,
                             &_seq->subsampleK, &_seq->subsampleN,
		                     &_seq->softMaskFilename, &_seq->softMaskComplement,
		                     &_seq->xMaskFilename, &_seq->xMaskComplement,
		                     &_seq->nMaskFilename, &_seq->nMaskComplement,
		                     &_seq->nameParseType,
		                     &_seq->nameTrigger,
		                     &_seq->doRevCompFlags,
		                     &_seq->doUnmask,
		                     &_seq->doPartitioning, &_seq->doJoin,
		                     &_seq->separatorCh,
		                     &_seq->useFullNames,
		                     &forcedfileType,
		                     &isQuantum, &qCodingFilename,
		                     &_seq->startLimit, &_seq->endLimit,
		                     &_seq->endIsSoft);
		_seq->f = fopen_or_die (_seq->filename, "rb");
		if (_seq->header != NULL)
			{
			_seq->headerSize   = strlen (_seq->header) + 1;
			_seq->lockedHeader = true;
			_seq->hasNickname  = true;
			}
		}

	if ((!choresAllowed) && (_seq->choresFilename != NULL))
		suicidef ("can't use [chores] for the target file (%s)\n"
		          "move [chores] to the query file, or use the --chores option",
		          sequence_filename(_seq));

	if ((_seq->doJoin) && (_seq->choresFilename != NULL))
		{
		if (choresAllowed) suicidef ("can't use --chores with [multiple]");
		              else suicidef ("can't use [chores] with [multiple]");
		}

	// init any non-zero fields

	_seq->hasSavedState = false;
	_seq->rewindable    = -1;	// (rewindability unknown at this point)
	_seq->contig        = 0;

	if (isQuantum)
		{
		if ((fileType != seq_type_qdna) && (fileType != seq_type_unknown))
			suicidef ("clashing file type for %s ([quantum] used for %s file)",
			          sequence_filename(_seq), seqTypeNames[fileType]);
		_seq->fileType = seq_type_qdna;
		}
	else if (forcedfileType != seq_type_unknown)
		{
		if ((fileType != forcedfileType) && (fileType != seq_type_unknown))
			suicidef ("clashing file type for %s (%s used for %s file)",
			          sequence_filename(_seq), seqTypeNames[forcedfileType], seqTypeNames[fileType]);
		_seq->fileType = forcedfileType;
		}
	else
		_seq->fileType = fileType;

	// initialize subsampling

	if (_seq->subsampleN == 1)
		_seq->subsampleN = 0; // (subsampling 1 of 1 is meaningless)

	if (_seq->subsampleN == 0)
		_seq->subsampleSkip = 0;
	else
		_seq->subsampleSkip = _seq->subsampleK - 1;

	// if the sequences in this file are to be joined into a parititioned
	// sequence, initialize that

	if (_seq->doPartitioning)
		{
		enough_partitions (_seq, /*numPartitions*/ 100, /*poolSize*/ 0,
		                   /*anticipate*/ false, /*round up*/ true);
		_seq->partition.state = seqpart_empty;
		}

	// if the file type is not yet known, figure out what it is

	if (_seq->fileType == seq_type_unknown)
		_seq->fileType = detect_file_type (_seq);

	if ((!sequences_dbgAllowColors)
	 && (_seq->fileType == seq_type_csfasta))
		suicidef ("sorry, color space is not fully implemented yet");

	if ((_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx)
	 && (_seq->contigOfInterest != NULL))
		suicidef ("specific contig-of-interest only valid for 2bit or hsx files (%s)",
		          _seq->contigOfInterest);

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx)
	 && (_seq->namesFilename != NULL))
		suicidef ("sequence-subset file only valid for fasta, fastq, csfasta, 2bit or hsx files\n(%s)",
		          _seq->namesFilename);

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx)
	 && (_seq->choresFilename != NULL))
		suicidef ("chores file only valid for fasta, fastq, csfasta, 2bit or hsx files\n(%s)",
		          _seq->choresFilename);

	if ((_seq->fileType == seq_type_hsx)
	 && (parse_type(_seq->nameParseType) != name_parse_type_core))
		suicidef ("\"nameparse=\" is not valid for hsx files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (parse_type(_seq->nameParseType) == name_parse_type_alnum))
		suicidef ("\"nameparse=alphanum\" only valid for fasta, fastq, csfasta or 2bit files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->fileType != seq_type_2bit)
	 && (parse_type(_seq->nameParseType) == name_parse_type_darkspace))
		suicidef ("\"nameparse=darkspace\" only valid for fasta, fastq, csfasta or 2bit files");

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_csfasta)
	 && (_seq->nameTrigger != NULL))
		suicidef ("\"nameparse=tag:%s\" only valid for fasta, fastq or csfasta files", _seq->nameTrigger);

	if ((_seq->fileType == seq_type_csfasta)
	 && (_seq->separatorCh != 0))
		suicidef ("[separator=%c] is not allowed for csfasta files", _seq->separatorCh);

	// for fastq files, bind allocation for vq (qualities) to v (nucleotides)

	if (_seq->fileType == seq_type_fastq)
		_seq->needsVq = true;

	// pre-allocate;  note that we can't do this earlier, since we wouldn't
	// have known whether to allocate quality values

	if (allocLen != 0)
		sequence_long_enough (_seq, allocLen, false);

	// for quantum DNA files, attach the complement mapping

	if (_seq->fileType == seq_type_qdna)
		_seq->qToComplement = qToComplement;

	// make sure the file is rewindable if the caller requires it to be

	if (beRewindable)
		{
		if (load_sequence (_seq))
			{
			_seq->preLoaded = true;
			if (another_sequence (_seq))
				{
				err = test_rewindability (_seq);
				if (err != 0) goto not_rewindable;
				_seq->rewindable = true;
				}
			}
		}

	// for 2bit or hsx files, we need to read the header
	// $$$ should this be moved to before the pre-load above?

	if (_seq->fileType == seq_type_2bit)
		read_2bit_header (_seq);
	else if (_seq->fileType == seq_type_hsx)
		read_hsx_header (_seq);

	// if we have a contigs-of-interest file, open it, read the first line, and
	// decide whether we'll need to advance to that contig

	_seq->contigPending = false;
	searchForContig = false;

	if (_seq->namesFilename != NULL)
		{
		_seq->namesFile = fopen_or_die (_seq->namesFilename, "rt");
		if (!read_contig_name (_seq))
			suicidef ("contigs-of-interest file is empty: %s", _seq->namesFilename);

		searchForContig = true;
		if (_seq->preLoaded)
			{
			header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
			if (strcmp (header, _seq->nextContigName) == 0)
				searchForContig = false;
			}
		}

	// if we have a chores file, open it, read the first chore, and decide
	// whether we'll need to advance to that contig

	else if (_seq->choresFilename != NULL)
		{
		_seq->choresFile = fopen_or_die (_seq->choresFilename, "rt");
		_seq->choresLineNum = 0;
		if (!read_chore (_seq))
			suicidef ("chores file is empty: %s", _seq->choresFilename);

		header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;

		searchForContig = true;
		if (_seq->preLoaded)
			{
			if (strcmp (header, _seq->nextContigName) == 0)
				searchForContig = false;
			}
		else // if (!_seq->preLoaded)
			{
			if ((header != NULL) && (strcmp (header, _seq->nextContigName) == 0))
				{
				searchForContig = false;
				_seq->preLoaded = true;
				validate_rev_comp (_seq);
				}
			}
		}

	// if necessary, advance to the first specified contig

	if (searchForContig)
		{
		if (_seq->fileType == seq_type_fasta)
			find_next_fasta_coi (_seq);
		else if (_seq->fileType == seq_type_fastq)
			find_next_fastq_coi (_seq);
		else if (_seq->fileType == seq_type_csfasta)
			find_next_csfasta_coi (_seq);
		else if (_seq->fileType == seq_type_2bit)
			find_next_2bit_coi (_seq);
		else // if (_seq->fileType == seq_type_hsx)
			find_next_hsx_coi (_seq);
		}

	// if we have a quantum coding file, read it

	if (qCodingFilename != NULL)
		{
		_seq->qCoding = read_quantum_code_by_name (qCodingFilename);
		free_if_valid   ("open_sequence_file (qCodingFilename)", qCodingFilename);
		}

	return _seq;

// failure exits

not_rewindable:
	if (name == NULL) name = "(stdin)";
	suicidef_with_perror ("sequence file \"%s\" is not rewindable"
						  " (fseek returned %d): %s",
						  name, err, sequence_filename(_seq));
	return NULL; // (never gets here)
	}

//----------
//
// rewind_sequence_file--
//	Rewind a sequence file, so that the sequence(s) within it can be read again.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to rewind.
//
// Returns:
//	(nothing;  failures cause program termination)
//
//----------

void rewind_sequence_file
   (seq*	_seq)
	{
	int		err;

	if (_seq->f == NULL)  return;	// (no file to rewind)

	// decide whether we can take a short cut and ignore the rewind request

	if ((_seq->fileType == seq_type_2bit)		// (we have only one
	 && (_seq->contigOfInterest != NULL))		//  .. sequence in the file)
		{
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if ((_seq->fileType == seq_type_hsx)		// (we have only one
	 && (_seq->contigOfInterest != NULL))		//  .. sequence in the file)
		{
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if (_seq->contig < 2)						// (we haven't read more than
		{										//  .. one sequence from file)
		if (_seq->contig != 0) _seq->preLoaded = true;
		return;
		}

	if (_seq->fileType == seq_type_2bit)
		{
		_seq->twoBit.contigFilePos = _seq->twoBit.indexFilePos;
		_seq->twoBit.contigLoaded  = false;
		goto reset_file_data;
		}

	if (_seq->fileType == seq_type_hsx)
		{
		_seq->hsx.contigLoaded = false;
		goto reset_file_data;
		}

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if (_seq->rewindable == false)
		suicidef ("sequence file is not rewindable: %s", sequence_filename(_seq));

	// rewind the file and reset the data

	debugTextFile_1;
	rewind (_seq->f);

reset_file_data:
	_seq->len           = 0;
	_seq->contig        = 0;
	_seq->preLoaded     = false;
	_seq->pendingStack  = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen    = 0;
	_seq->hasSavedState = false;

	// re-initialize subsampling

	if (_seq->subsampleN == 0)
		_seq->subsampleSkip = 0;
	else
		_seq->subsampleSkip = _seq->subsampleK - 1;

	// rewind the contigs-of-interest file

	if (_seq->namesFile != NULL)
		{
		err = fseek (_seq->namesFile, 0, SEEK_SET);
		if (err != 0)
			suicidef ("failed to seek to position in file\n"
			          "in rewind_sequence_file for %s, index fseek(%08lX) returned %d",
			          _seq->namesFilename, 0, err);
		_seq->contigPending = false;
		}

	// rewind the chores file

	if (_seq->choresFile != NULL)
		{
		err = fseek (_seq->choresFile, 0, SEEK_SET);
		if (err != 0)
			suicidef ("failed to seek to position in file\n"
			          "in rewind_sequence_file for %s, index fseek(%08lX) returned %d",
			          _seq->choresFilename, 0, err);
		_seq->contigPending = false;
		}

	if (_seq->fileType == seq_type_hsx)
		locate_hsx_first_sequence (_seq);
	}

//----------
//
// clone_sequence--
//	Make a copy of a sequence, so that the copy is in the same state as if
//	open_rewindable_sequence_file() had been called.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to make a clone of.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

seq* clone_sequence
   (seq*	_seq)
	{
	seq*	newSeq;

	debugCloning_1;
	newSeq = copy_sequence (_seq);
	newSeq->preLoaded = true;
	debugCloning_2;

	return newSeq;
	}

//----------
//
// copy_sequence--
//	Make a copy of a sequence.  Note that the copy is not as functional as the
//	original.  For example, file operations cannot be performed on the copy.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to make a copy of.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

seq* copy_sequence
   (seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	seq*			newSeq;
	unspos			ix;

	// allocate the sequence tracking structure

	newSeq = alloc_sequence_record ("copy_sequence");

	newSeq->pendingChars = zalloc_or_die ("copy_sequence (pendingChars)",
	                                      seqBufferSize);
	newSeq->pendingStack = newSeq->pendingChars + seqBufferSize;
	newSeq->pendingLen   = 0;

	// allocate the sequence content

	if ((_seq->v == NULL) || (_seq->len < 1))
		{
		newSeq->v       = NULL;
		newSeq->vOwner  = true;
		newSeq->size    = 0;
		newSeq->len     = 0;
		newSeq->trueLen = 0;
		}
	else
		{
		if (_seq->len > maxSequenceLen)
			suicidef ("in copy_sequence, "
			          "sequence length " unsposFmt " exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		newSeq->v       = malloc_or_die ("copy_sequence (v)", _seq->len+1);
		newSeq->vOwner  = true;
		newSeq->size    = _seq->len+1;
		newSeq->len     = _seq->len;
		newSeq->trueLen = _seq->trueLen;
		}

	if (_seq->vc == NULL) newSeq->vc = NULL;
	                 else newSeq->vc = malloc_or_die ("copy_sequence (vc)", _seq->len+1);
	newSeq->vcOwner = true;

	if (_seq->vq == NULL) newSeq->vq = NULL;
	                 else newSeq->vq = malloc_or_die ("copy_sequence (vq)", _seq->len+1);
	newSeq->vqOwner = true;

	// set up file info

	newSeq->fileType = _seq->fileType;

	if (_seq->filename != NULL)
		newSeq->filename = copy_string (_seq->filename);

	// copy the sequence content (including a terminating zero)

	if (newSeq->v != NULL)
		{
		for (ix=0 ; ix<=newSeq->len ; ix++)
			newSeq->v[ix] = _seq->v[ix];
		}

	newSeq->preLoaded = _seq->preLoaded;
	newSeq->contig    = _seq->contig;

	// copy the partition info

	if (sp->p != NULL)
		{
		if (sp->state != seqpart_ready)
			suicidef ("internal error, attempt to copy sequence (%s) in partition state %d",
			          _seq->filename, sp->state);
		newSeq->partition.state = seqpart_ready;
		copy_partitions (/*to*/ newSeq, /*from*/ _seq);
		}

	// copy other fields

	newSeq->startLoc     = _seq->startLoc;
	newSeq->revCompFlags = _seq->revCompFlags;
	newSeq->contig       = _seq->contig;
	newSeq->lockedHeader = _seq->lockedHeader;
	newSeq->allowAmbiDNA = _seq->allowAmbiDNA;

	if (_seq->header == NULL)
		{
		newSeq->header     = NULL;
		newSeq->headerSize = 0;
		}
	else
		{
		newSeq->header     = copy_string (_seq->header);
		newSeq->headerSize = strlen (newSeq->header) + 1;
		}
	newSeq->headerOwner = true;

	if (_seq->shortHeader == NULL)
		{
		newSeq->shortHeader     = NULL;
		newSeq->shortHeaderSize = 0;
		}
	else
		{
		newSeq->shortHeader      = copy_string (_seq->shortHeader);
		newSeq->shortHeaderSize  = strlen (newSeq->shortHeader) + 1;
		}
	newSeq->shortHeaderOwner = true;

	if (_seq->trueHeader == NULL)
		{
		newSeq->trueHeader     = NULL;
		newSeq->trueHeaderSize = 0;
		}
	else
		{
		newSeq->trueHeader      = copy_string (_seq->trueHeader);
		newSeq->trueHeaderSize  = strlen (newSeq->trueHeader) + 1;
		}
	newSeq->trueHeaderOwner = true;

	newSeq->startLimit   = _seq->startLimit;
	newSeq->endLimit     = _seq->endLimit;
	newSeq->endIsSoft    = _seq->endIsSoft;
	newSeq->useFullNames = _seq->useFullNames;

	if (_seq->contigOfInterest == NULL)
		newSeq->contigOfInterest = NULL;
	else
		newSeq->contigOfInterest = copy_string (_seq->contigOfInterest);

	return newSeq;
	}

//----------
//
// new_sequence--
//	Create a new, empty, sequence.
//
//----------
//
// Arguments:
//	unspos	allocLen:	The sequence length to allocate for.  The special
//						.. value seqposInfinity indicates that no memory should
//						.. be allocated for the sequence vector.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

seq* new_sequence
   (unspos	allocLen)
	{
	seq*	_seq;

	// allocate the sequence tracking structure

	_seq = alloc_sequence_record ("new_sequence");

	_seq->pendingChars = zalloc_or_die ("new_sequence (pendingChars)",
	                                   seqBufferSize);
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	_seq->pendingLen   = 0;

	// allocate space for sequence data (but leave it empty)

	if (allocLen == seqposInfinity)
		{
		_seq->v       = NULL;
		_seq->vc      = NULL;
		_seq->vq      = NULL;
		_seq->vOwner  = false;
		_seq->vcOwner = false;
		_seq->vqOwner = false;
		_seq->size    = 0;
		_seq->len     = 0;
		}
	else
		{
		_seq->v       = malloc_or_die ("new_sequence (v)", allocLen+1);
		_seq->vc      = NULL;
		_seq->vq      = NULL;
		_seq->vOwner  = true;
		_seq->vcOwner = true;
		_seq->vqOwner = true;
		_seq->size    = allocLen+1;
		_seq->len     = 0;
		_seq->v[0]    = 0;
		}

	// initialize the other fields

	_seq->fileType = seq_type_nofile;
	_seq->contig   = 1;
	_seq->startLoc = 1;
	_seq->trueLen  = 0;

	return _seq;
	}

//----------
//
// alloc_sequence_record--
//	Allocate a new sequence tracking structure, and make sure all pointer
//	fields are NULL.
//
//----------
//
// Arguments:
//	char*	id:	an identifying string to be used when trackMemoryUsage is
//				.. turned on;  this can be NULL.
//
// Returns:
//	A pointer to the sequence;  failures result in fatality.  The caller must
//	eventually de-allocate this by calling free_sequence().
//
//----------

static seq* alloc_sequence_record
   (arg_dont_complain(char* id))
	{
	seq*	_seq;

	_seq = zalloc_or_die (id, sizeof(*_seq));

	_seq->v                   = NULL;
	_seq->vc                  = NULL;
	_seq->vq                  = NULL;
	_seq->needsVq             = false;
	_seq->pendingChars        = NULL;
	_seq->filename            = NULL;
	_seq->header              = NULL;
	_seq->shortHeader         = NULL;
	_seq->trueHeader          = NULL;
	_seq->f                   = NULL;
	_seq->namesFile           = NULL;
	_seq->namesFilename       = NULL;
	_seq->choresFile          = NULL;
	_seq->choresFilename      = NULL;
	_seq->subsampleK          = 0;
	_seq->subsampleN          = 0;
	_seq->softMaskFilename    = NULL;
	_seq->xMaskFilename       = NULL;
	_seq->nMaskFilename       = NULL;
	_seq->nameTrigger         = NULL;
	_seq->contigOfInterest    = NULL;
	_seq->twoBit.nBlockStarts = NULL;
	_seq->twoBit.nBlockSizes  = NULL;
	_seq->twoBit.mBlockstarts = NULL;
	_seq->twoBit.mBlocksizes  = NULL;
	_seq->partition.p         = NULL;
	_seq->partition.pool      = NULL;
	_seq->allowAmbiDNA        = false;
	_seq->qToComplement       = NULL;

	return _seq;
	}

//----------
//
// sequence_long_enough--
//	Make sure a sequence has enough room (including an extra byte for a
//	terminating zero).
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to check.
//	unspos	allocLen:	The sequence length to allocate for (not including the
//						.. terminator).
//	int		anticipate:	true  => allocate extra, anticipating the need for more
//						false => don't
//
// Returns:
//	nothing;  the sequence's v[] (and possible vq, see note 1) may be modified;
//	failures result in fatality.
//
//----------
//
//	Notes:
//
//	(1)	Normally we only allocate v[].  But if _seq->needsVq is true, we also
//		allocate vq[], and keep it the same length as v[].
//
//----------

void sequence_long_enough
   (seq*	_seq,
	unspos	allocLen,
	int		anticipate)
	{
	char*	name;

	if (_seq->size >= allocLen+1)
		return;

	allocLen += 2;						// (add space for a terminating zero,
	if (anticipate)						//  .. etc.)
		allocLen += 30 + allocLen / 8;	// anticipatory, grow by about 13%
	allocLen = round_up_16K (allocLen);	// we expect that allocation in
										// .. multiples of 16K is better for
										// .. the heap manager

	if (!_seq->vOwner)
		{
		name = (_seq->filename != NULL)? _seq->filename
		                               : _seq->header;
		suicidef ("internal error, attempt to resize external sequence (%s)",
		          name);
		}

	if ((_seq->needsVq) && (!_seq->vqOwner))
		{
		name = (_seq->filename != NULL)? _seq->filename
		                               : _seq->header;
		suicidef ("internal error, attempt to resize external qualities (%s)",
		          name);
		}

	// fprintf (stderr, "re-allocating " unsposFmt " bytes\n", allocLen);
	_seq->v    = realloc_or_die ("sequence_long_enough", _seq->v, allocLen);
	_seq->size = allocLen;

	if (_seq->needsVq)
		_seq->vq = realloc_or_die ("sequence_long_enough (qualities)", _seq->vq, allocLen);
	}

//----------
//
// free_sequence--
//	Deallocate a sequence, along with any associated memory or files.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_sequence
   (seq*	_seq)
	{
	if (_seq == NULL) return;

	if (_seq->vOwner)           free_if_valid ("free_sequence (v)",                _seq->v);
	if (_seq->vcOwner)          free_if_valid ("free_sequence (vc)",               _seq->vc);
	if (_seq->vqOwner)          free_if_valid ("free_sequence (vq)",               _seq->vq);
	                            free_if_valid ("free_sequence (pendingChars)",     _seq->pendingChars);
	                            free_if_valid ("free_sequence (filename)",         _seq->filename);
	if (_seq->headerOwner)      free_if_valid ("free_sequence (header)",           _seq->header);
	if (_seq->shortHeaderOwner) free_if_valid ("free_sequence (shortHeader)",      _seq->shortHeader);
	if (_seq->trueHeaderOwner)  free_if_valid ("free_sequence (trueHeader)",       _seq->trueHeader);
	                            free_if_valid ("free_sequence (qCoding)",          _seq->qCoding);

	                            free_if_valid ("free_sequence (contigOfInterest)", _seq->contigOfInterest);

	if (_seq->fileType != seq_type_nofile)
		{
		fclose_if_valid (_seq->f);
		fclose_if_valid (_seq->namesFile);
		free_if_valid   ("free_sequence (namesFilename)",     _seq->namesFilename);
		fclose_if_valid (_seq->choresFile);
		free_if_valid   ("free_sequence (choresFilename)",    _seq->choresFilename);
		free_if_valid   ("free_sequence (softMaskFilename)",  _seq->softMaskFilename);
		free_if_valid   ("free_sequence (xMaskFilename)",     _seq->xMaskFilename);
		free_if_valid   ("free_sequence (nMaskFilename)",     _seq->nMaskFilename);
		free_if_valid   ("free_sequence (nameTrigger)",       _seq->nameTrigger);
		}

	free_if_valid ("free_sequence (twoBit.nBlockStarts)",     _seq->twoBit.nBlockStarts);
	free_if_valid ("free_sequence (twoBit.nBlockSizes)",      _seq->twoBit.nBlockSizes);
	free_if_valid ("free_sequence (twoBit.mBlockstarts)",     _seq->twoBit.mBlockstarts);
	free_if_valid ("free_sequence (twoBit.mBlocksizes)",      _seq->twoBit.mBlocksizes);

	if (_seq->hsx.fileInfo != NULL)
		{
		u32 fileNum;
		for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
			fclose_if_valid (_seq->hsx.fileInfo[fileNum].f);
		free_if_valid ("free_sequence (hsx.fileInfo)",        _seq->hsx.fileInfo);
		}

	                               free_if_valid ("free_sequence (partition.p)",    _seq->partition.p);
	if (_seq->partition.poolOwner) free_if_valid ("free_sequence (partition.pool)", _seq->partition.pool);

	free_if_valid ("free_sequence (_seq)", _seq);
	}

//----------
//
// load_sequence--
//	Load the next sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//
// Returns:
//	true if there was another sequence to load, false if not;  failures other
//	than eof result in fatality.
//
//----------

int load_sequence
   (seq*			_seq)
	{
	seqpartition*	sp;
	partition*		p;
	unspos			sepPos;
	unspos			oldTrueLen;

	debugNamesFile_1;

	if (_seq == NULL)                      suicide ("load_sequence(NULL)");
	if (_seq->preLoaded)                   { _seq->preLoaded = false;  return true; }
	if (_seq->fileType == seq_type_nofile) return false;
	if (!another_sequence (_seq))          return false;
	// (another_sequence() may set _seq->preLoaded, so we must check a second time)
	if (_seq->preLoaded)                   { _seq->preLoaded = false;  return true; }

	debugNamesFile_2;

	// get rid of sequence data from previous load

	_seq->len     = 0;
	_seq->trueLen = 0;

	//////////
	// read the sequence data, either
	//	o sequence is not partitioned => read just one sequence from the file
	//	o sequence is partitioned     => read all sequences from the file
	//////////

	sp = &_seq->partition;
	if (!_seq->doPartitioning)
	 	{
		// not partitioned, just read the next sequence

		if (_seq->subsampleN > 0)
			{ // we're subsampling, skip sequences as appropriate
			if (_seq->subsampleSkip > 0)
				skip_sequences (_seq, _seq->subsampleSkip);
			_seq->subsampleSkip = _seq->subsampleN-1;
			}

		load_sequence_core (_seq, /*keeper*/ true);
		}

	else if (_seq->doJoin)
		{
		// partitioned, read all the sequences

		sp->state = seqpart_loading;

		// write first separator

		sequence_long_enough (_seq, 1, false);
		_seq->v[0] = 0;
		_seq->len  = 0;
		sp->len    = 0;

		// read all the sequences

		while (another_sequence (_seq))
			{
			// load the next sequence

			sepPos = _seq->len++;		// (advance length past the separator)
			oldTrueLen = _seq->trueLen;

			if (_seq->subsampleN > 0)
				{ // we're subsampling, skip sequences as appropriate
				if (_seq->subsampleSkip > 0)
					skip_sequences (_seq, _seq->subsampleSkip);
				_seq->subsampleSkip = _seq->subsampleN-1;
				}
			load_sequence_core (_seq, /*keeper*/ true);
            add_partition (_seq, sepPos, _seq->startLoc, _seq->trueLen-oldTrueLen);
			}

		// add final separator to the table

		p = &sp->p[sp->len-1];
		p->sepAfter = _seq->len;

		p = &sp->p[sp->len];
		p->sepBefore = _seq->len; 

		sp->state = seqpart_ready;
		}

	else if (!_seq->doJoin)
		{
		// partitioned, but only so as to allow separators, just read the next
		// sequence

		sp->state = seqpart_loading;

		// write first separator

		sequence_long_enough (_seq, 1, false);
		_seq->v[0] = 0;
		_seq->len  = 0;
		sp->len    = 0;

		// load the sequence

		sepPos = _seq->len++;		// (advance length past the separator)
		oldTrueLen = _seq->trueLen;

		if (_seq->subsampleN > 0)
			{ // we're subsampling, skip sequences as appropriate
			if (_seq->subsampleSkip > 0)
				skip_sequences (_seq, _seq->subsampleSkip);
			_seq->subsampleSkip = _seq->subsampleN-1;
			}
		load_sequence_core (_seq, /*keeper*/ true);
		add_partition (_seq, sepPos, _seq->startLoc, _seq->trueLen-oldTrueLen);

		// add final separator to the table

		p = &sp->p[sp->len-1];
		p->sepAfter = _seq->len;

		p = &sp->p[sp->len];
		p->sepBefore = _seq->len; 

		sp->state = seqpart_ready;
		}

	// apply any required operators to it;  note that nib and 2bit sequences
	// are unmasked (if desired) during the earlier call to load_nib_sequence
	// or load_2bit_sequence, so there is no need to unmask them here;  further,
	// csfasta files do not have the concept of masking

	if ((_seq->doUnmask)
	 && ((_seq->fileType == seq_type_fasta)
	  || (_seq->fileType == seq_type_fastq)
	  || (_seq->fileType == seq_type_hsx)))
		upper_sequence (_seq);

	if (_seq->fileType == seq_type_qdna)
		{
		if ((_seq->softMaskFilename != NULL)
		 || (_seq->xMaskFilename    != NULL)
		 || (_seq->nMaskFilename    != NULL))
			suicidef ("masking not allowed for %s", sequence_filename(_seq));
		if (_seq->doUnmask)
			suicidef ("unmasking not allowed for %s", sequence_filename(_seq));
		if (((_seq->doRevCompFlags & rcf_comp) != 0)
		  && (_seq->qToComplement == NULL))
			suicidef ("reverse complement not allowed for %s\n",
			          "(the score file lacks complements)",
			           sequence_filename(_seq));
		}

	if (_seq->softMaskFilename != NULL)
		{
		if (_seq->softMaskComplement)
			mask_sequence_keep (_seq, _seq->softMaskFilename, -1);
		else
			mask_sequence      (_seq, _seq->softMaskFilename, -1);
		}
	if (_seq->xMaskFilename != NULL)
		{
		if (_seq->xMaskComplement)
			mask_sequence_keep (_seq, _seq->xMaskFilename, 'X');
		else
			mask_sequence      (_seq, _seq->xMaskFilename, 'X');
		}
	if (_seq->nMaskFilename != NULL)
		{
		if (_seq->nMaskComplement)
			mask_sequence_keep (_seq, _seq->nMaskFilename, 'N');
		else
			mask_sequence      (_seq, _seq->nMaskFilename, 'N');
		}

	if (_seq->separatorCh != 0)
		separate_sequence (_seq, _seq->separatorCh);

	if (_seq->doRevCompFlags == rcf_revcomp)
		rev_comp_sequence (_seq, _seq->qToComplement);
	else if (_seq->doRevCompFlags == rcf_rev)
		backward_sequence (_seq);
	else if (_seq->doRevCompFlags == rcf_comp)
		{
		backward_sequence (_seq);
		rev_comp_sequence (_seq, _seq->qToComplement);
		}

	//debugPartitions_1a;
	debugPartitions_1b;

	if (sequences_dbgDumpSequence)
		dump_sequence (stderr, _seq);

	return true;
	}


//-- skip_sequences--

static void skip_sequences
   (seq*	_seq,
	int		skipCount)
	{
	while ((skipCount-- > 0) && another_sequence_core (_seq))
		load_sequence_core (_seq, /*keeper*/ false);
	}


//-- load_sequence_core --

static void load_sequence_core
   (seq*	_seq,
	int		keeper)
	{
	debugNamesFile_3;

	// get rid of header data from previous load

	if (!_seq->lockedHeader)
		{
		if ((_seq->header != NULL) && (_seq->headerSize != 0))
			_seq->header[0] = 0;
		if ((_seq->shortHeader != NULL) && (_seq->shortHeaderSize != 0))
			_seq->shortHeader[0] = 0;
		if ((_seq->trueHeader != NULL) && (_seq->trueHeaderSize != 0))
			_seq->trueHeader[0] = 0;
		}

	// read the next sequence for this type

	_seq->revCompFlags = rcf_forward;
	_seq->contig++;

	switch (_seq->fileType)
		{
		case seq_type_fasta:
			load_fasta_sequence (_seq, keeper);
			break;

		case seq_type_fastq:
			load_fastq_sequence (_seq, keeper);
			break;

		case seq_type_csfasta:
			load_csfasta_sequence (_seq, keeper);
			break;

		case seq_type_nib:
			load_nib_sequence (_seq, keeper);
			break;

		case seq_type_2bit:
			if (_seq->contigOfInterest != NULL)
				_seq->contig--; // (cancel earlier increment)
			load_2bit_sequence (_seq, keeper);
			break;

		case seq_type_hsx:
			if (_seq->contigOfInterest != NULL)
				_seq->contig--; // (cancel earlier increment)
			load_hsx_sequence (_seq, keeper);
			break;

		case seq_type_qdna:
			load_qdna_sequence (_seq, keeper);
			break;

		default:
			suicidef ("unknown sequence type: %X", _seq->fileType);
		}

	debugNamesFile_4;

	_seq->contigPending = false;

	if ((_seq->header != NULL)
	 && (_seq->headerSize != 0)
	 && ((_seq->shortHeader == NULL)
	  || (_seq->shortHeaderSize == 0)
	  || (_seq->shortHeader[0] == 0)
	  || (_seq->hasNickname)))
		create_short_header (_seq);

	if ((_seq->header != NULL)
	 && (_seq->headerSize != 0)
	 && ((_seq->nameParseType & name_parse_fill_white) != 0))
		whitespace_to_under (_seq->header, strlen(_seq->header));
	}

//----------
//
// load_fasta_sequence--
//	Load the next fasta sequence from the associated file.
//
// A typical file looks like this:
//
//	> some header for the first sequence
//	GCGGTATCGCGCACAAGATTTAGGGATAGATCGTTTTGATGACCTCTCGCCACCTGGCAA
//	  ...
//	AAAAAAGGTAGGCCCATTAGCCCCCC
//
// The header line is optional.  However, if several sequences are included in
// the same file, the header lines are necessary as separators.  Nucleotides
// can be upper or lower case.  X can be used to indicate a masked position.
// whitespace can be added as desired (and so line breaks can be anywhere).
// digits are ignored so it is easy to use numerically annotated sequences.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_fasta_sequence
   (seq*	_seq,
	int		keeper)
	{
	int		ch;
	unspos	length;

	if (_seq == NULL) suicide ("load_fasta_sequence(NULL)");

	debugFastaFile_1;

	//////////
	// read the header
	//////////

	ch = skip_seq_whitespace (_seq);

	if (ch != '>') seq_ungetc (ch, _seq);
	          else parse_fasta_header (_seq);

	debugFastaFile_2;

	//////////
	// read ahead to determine the length of the sequence and to pre-allocate
	// the vector
	//////////

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if ((_seq->rewindable == true) && (keeper))
		{
		// read ahead, counting chars needed

		save_fstate (_seq);
		length = parse_fasta (_seq, /*storeEm*/ false);
		restore_fstate (_seq);

		// allocate the vector

		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_fasta_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));

		sequence_long_enough (_seq, _seq->len+length, false);
		}

	//////////
	// read the sequence
	//////////

	parse_fasta (_seq, /*storeEm*/ keeper);
	}

//----------
//
// parse_fasta_header--
//	Parse a fasta header from the associated file.  This assumes that the
//	sequence's file is positioned at the first character in the header, *after*
//	the '>' character.
//
//	Upon return, the file is positioned at the start of the first line
//	following the header.
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence being parsed.
//
// Returns:
//  (nothing; the header is written to _seq->header)
//
//----------

static void parse_fasta_header
   (seq*	_seq)
	{
	u32		headerLen;
	int		ch;

	if (_seq->lockedHeader)
		{
		ch = seq_getc (_seq);
		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			ch = seq_getc (_seq);
		}
	else
		{
		headerLen = 0;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		if (sequences_keepFastaArrow)
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, '>');
			ch = seq_getc (_seq);
			while ((ch == ' ') || (ch == '\t'))
				{
				append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
				ch = seq_getc (_seq);
				}
			}
		else
			ch = skip_seq_whitespace (_seq);

		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
			ch = seq_getc (_seq);
			}
		append_char (&_seq->header, &_seq->headerSize, &headerLen, 0);

		if (_seq->nameTrigger != NULL)
			{
			char* triggerFound, *src, *dst;
			triggerFound = strstr (_seq->header, _seq->nameTrigger);
			if (triggerFound != NULL)
				{
				triggerFound += strlen (_seq->nameTrigger);
				for (src=triggerFound,dst=_seq->header ; *src!=0 ; )
					{
					ch = *(src++);
					if ((!isalnum(ch)) && (ch != '_'))
						break;
					*(dst++) = ch;
					}
				*dst = 0;
				}
			}
		}

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	}

//----------
//
// parse_fasta--
//	Parse a fasta sequence from the associated file.  This assumes that the
//	sequence's file is positioned at the first character in the sequence,
//	*after* the sequence header line.
//
// (see load_fasta_sequence() for info about the file format)
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to parse.
//	int		storeEm:	true  => store the results (in the sequence)
//						false => just count
//
// Returns:
//  The number of characters read into the sequence;  failure causes program
//	fatality.
//
//----------

static unspos parse_fasta
   (seq*	_seq,
	int		storeEm)
	{
	unspos	index, startLimit, endLimit, count;
	int		prevCh, ch;
	char*	description;

	index      = 0;
	count      = 0;
	startLimit = _seq->startLimit;
	endLimit   = _seq->endLimit;

	// scan the file, keeping characters that are (a) nucleotides and (b) are
	// within our index limits

	prevCh = '\n';
	ch     = skip_seq_whitespace (_seq);
	while (ch != EOF)
		{
		if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
			{ seq_ungetc (ch, _seq);  break; }

		if ((_seq->separatorCh == 0) || (ch != _seq->separatorCh))
			{
			switch (char_to_fasta_type[(u8)ch])
				{
				case _nucleotide:
					break;
				case _ambiguous:
					if (!_seq->allowAmbiDNA)
						goto bad_char;
					break;
				case _newline:
					ch = '\n'; // (allow for unix, mac, or pc line ends)
					goto next_char;
				case _bad:
					goto bad_char;
				}
			}

		// this is a nucleotide (or separator), do we want it?

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// we want it;  are we just counting?

		if ((!storeEm) && (count+1 < count))
			suicidef ("in parse_fasta, "
			          "sequence length " unsposFmt "+1 overflows internal data type",
			          count);

		count++;
		if (!storeEm) goto next_char;

		// ok, let's store it

		if (_seq->len > maxSequenceLen - 1)
			suicidef ("in parse_fasta, "
			          "sequence length " unsposFmt "+1 exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		sequence_long_enough (_seq, _seq->len+1, true);
		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		ch = skip_seq_whitespace (_seq);
		}

	if (storeEm)
		{
		_seq->v[_seq->len] = 0;				// (set the terminating zero)
		_seq->trueLen += index;				// (account for the characters
											//  .. we've read so far)
		}

	// make sure we got somethin' useful

	if ((startLimit != 0) && (startLimit > index))
		goto beyond_start;

	if ((endLimit != 0) && (endLimit > index))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			goto beyond_end;
		}

	if ((count == 0) && (storeEm))
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		}

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

// (no longer needed;  the loop above exits only at end-of-file or end-of-seq)
//	// skip to the next sequence
//
//	if ((storeEm) && (_seq->needTrueLen))
//		{
//		prevCh = '\n';
//		ch     = skip_seq_whitespace (_seq);
//		while (ch != EOF)
//			{
//			if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
//				{ seq_ungetc (ch, _seq);  break; }
//
//			switch (char_to_fasta_type[(u8)ch])
//				{
//				case _nucleotide:
//				case _ambiguous:
//					_seq->trueLen++;
//					break;
//				case _newline:
//					ch = '\n'; // (allow for unix, mac, or pc line ends)
//					break;
//				case _bad:
//					goto bad_char;
//				}
//
//			// go try the next character
//
//			prevCh = ch;
//			ch = skip_seq_whitespace (_seq);
//			}
//		}

	return count;

	// failure exits
	// $$$ report line number here

bad_char:
	description = char_to_description (ch);
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("bad fasta character in %s (%s)\n"
		          "remove or replace non-ACGTN characters or consider using --ambiguous=iupac",
		          sequence_filename(_seq), description);
	else
		suicidef ("bad fasta character in %s, %s (%s)\n"
		          "remove or replace non-ACGTN characters or consider using --ambiguous=iupac",
				  sequence_filename(_seq), _seq->header, description);

beyond_start:
	if (_seq->header == NULL)
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), startLimit, index);
	else
		suicidef ("beyond end in %s, %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->header, startLimit, index);

beyond_end:
	if (_seq->header == NULL)
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), endLimit, index);
	else
		suicidef ("beyond end in %s, %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->header, endLimit, index);

	return 0; // (never gets here)
	}

//----------
//
// load_fastq_sequence--
//	Load the next fastq sequence from the associated file.
//
// A typical file looks like this:
//
//	@HWI-ST407_110227_0090_A80FT9ABXX:1:1:1190:2064#0/1
//	AGCTAAGGAATGACACAATTTGTCCTAATGGCAAATGCAGGGATTGTGATAAATATATCCNATATCTTA
//	+
//	cccWXccbcc[aZ_bRaaa`edadeefffff\cdaIXPWZdd`adcXaXN_BBBBBBBBBBBBBBBBBB
//	@HWI-ST407_110227_0090_A80FT9ABXX:1:1:1120:2089#0/1
//	TAGAACATGAGGGAAAGGAACAACCCTGCTGACTGACATGAGGCTGCCTGCCGCGGGGGGATGGGCAGG
//	+
//	c_aaeaaWcdddcaNZa`baaXa^VYX\V_[[[]TVJJOYbbcBBBBBBBBBBBBBBBBBBBBBBBBBB
//	 ...
//
// In each four-line block, the first line is @name.  The third line may also
// contain the name, in which case it much match the first line except that it
// starts with a plus sign.  The second line is the nucleotides.  The fourth
// line is qualities.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_fastq_sequence
   (seq*	_seq,
	int		keeper)
	{
	int		ch;
	char*	description;

	if (_seq == NULL) suicide ("load_fastq_sequence(NULL)");

	//////////
	// read the header
	//////////

	ch = seq_getc (_seq);
	if (ch == EOF) goto end_of_file;
	if (ch != '@') goto bad_fastq_header;

	parse_fastq_header (_seq);

	//////////
	// read the sequence
	//////////

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	sequence_long_enough (_seq, _seq->len+maxFastqSequenceLen, false);

	if (!keeper) fastq_skip_content (_seq);
	        else parse_fastq (_seq);

	return;

	// failure exits
	// $$$ report line number here

end_of_file:
	suicidef ("premature end of fastq file %s\n",
	          sequence_filename(_seq));

bad_fastq_header:
	description = char_to_description (ch);
	suicidef ("bad fastq header character in %s (expected \"@\" but read \"%s\")\n",
	          sequence_filename(_seq), description);

	return; // (never gets here)
	}

//----------
//
// parse_fastq_header--
//	Parse a fastq header from the associated file.  This assumes that
//	the sequence's file is positioned at the first character in the header,
//	*after* the '@' character.
//
//	Upon return, the file is positioned at the start of the first line
//	following the header.
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence being parsed.
//
// Returns:
//  (nothing; the header is written to _seq->header)
//
//----------

static void parse_fastq_header
   (seq*	_seq)
	{
	u32		headerLen;
	int		ch;
	char*	s;

	// parse the complete header and save it in trueHeader

	_seq->trueHeaderOwner = true;

	headerLen = 0;
	ch = seq_getc (_seq);
	while ((ch != '\n') && (ch != '\r') && (ch != EOF))
		{
		append_char (&_seq->trueHeader, &_seq->trueHeaderSize, &headerLen, ch);
		ch = seq_getc (_seq);
		}
	append_char (&_seq->trueHeader, &_seq->trueHeaderSize, &headerLen, 0);

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	// copy from trueHeader into the header (unless the header is locked), then
	// if a name trigger is active and matched, handle it

	if (_seq->lockedHeader)
		{
		ch = seq_getc (_seq);
		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			ch = seq_getc (_seq);
		}
	else
		{
		_seq->headerOwner = _seq->shortHeaderOwner = true;

		headerLen = 0;
		for (s=_seq->trueHeader ; (*s)!=0 ; s++)
			append_char (&_seq->header, &_seq->headerSize, &headerLen, *s);
		append_char (&_seq->header, &_seq->headerSize, &headerLen, 0);

		if (_seq->nameTrigger != NULL)
			{
			char* triggerFound, *src, *dst;
			triggerFound = strstr (_seq->header, _seq->nameTrigger);
			if (triggerFound != NULL)
				{
				triggerFound += strlen (_seq->nameTrigger);
				for (src=triggerFound,dst=_seq->header ; *src!=0 ; )
					{
					ch = *(src++);
					if ((!isalnum(ch)) && (ch != '_'))
						break;
					*(dst++) = ch;
					}
				*dst = 0;
				}
			}
		}

	}

//----------
//
// parse_fastq--
//	Parse a fastq sequence from the associated file.  This assumes that the
//	sequence's file is positioned at the first character in the sequence,
//	*after* the sequence header line.
//
// (see load_fastq_sequence() for info about the file format)
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to parse.
//
// Returns:
//  The number of characters read into the sequence;  failure causes program
//	fatality.
//
//----------

static unspos parse_fastq
   (seq*	_seq)
	{
	unspos	index, startLimit, endLimit, nucCount;
	unspos	qualCount, qualLen;
	int		ch;
	int		checkVsHeader;
	u32		headerIx;
	char*	description;

	headerIx = qualCount = 0; // (placate compiler)

	qualLen    = _seq->len;
	startLimit = _seq->startLimit;
	endLimit   = _seq->endLimit;

	//////////
	// read nucleotides and keep any that are within our index limits
	//////////

	index    = 0;
	nucCount = 0;

	ch = seq_getc (_seq);
	while (ch != EOF)
		{
		if ((_seq->separatorCh == 0) || (ch != _seq->separatorCh))
			{
			switch (char_to_fasta_type[(u8)ch])
				{
				case _nucleotide:
					break;
				case _ambiguous:
					if (!_seq->allowAmbiDNA) goto bad_nucleotide;
					break;
				case _newline:
					goto end_of_nucleotides;
				case _bad:
					goto bad_nucleotide;
				}
			}

		// this is a nucleotide (or separator), do we want it?

		index++;
		if ((startLimit != 0) && (index < startLimit)) goto next_nucleotide;
		if ((endLimit   != 0) && (index > endLimit))   goto next_nucleotide;

		// we want it;  let's store it

		nucCount++;

		if (_seq->len > maxSequenceLen - 1)
			suicidef ("in parse_fastq, "
			          "sequence length " unsposFmt "+1 exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		sequence_long_enough (_seq, _seq->len+1, true);
		_seq->v[_seq->len++] = ch;

		// go try the next nucleotide character

	next_nucleotide:
		ch = seq_getc (_seq);
		}

	// we've read all the nucleotides, finish up the sequence info

end_of_nucleotides:

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	_seq->v[_seq->len] = 0;				// (set the terminating zero)
	_seq->trueLen += index;				// (account for the characters
										//  .. we've read so far)

	// make sure we got somethin' useful

	if ((startLimit != 0) && (startLimit > index))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), startLimit, index);

	if ((endLimit != 0) && (endLimit > index))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), endLimit, index);
		}

	if (nucCount == 0)
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		}

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

	//////////
	// read and validate the 3rd line
	//////////

	// make sure it starts with a plus sign

	ch = seq_getc (_seq);
	if (ch == EOF) goto end_of_file;
	if (ch != '+') goto bad_fastq_third_line;

	// if there's nothing else on the line, we're done with it

	ch = seq_getc (_seq);
	if ((ch == '\n') || (ch == '\r')) goto end_of_third_line;

	// otherwise, it has to match the name that was given in the header

	checkVsHeader = (_seq->trueHeader != NULL) && (_seq->trueHeader[0] != 0);

	headerIx = 0;
	while ((ch != '\n') && (ch != '\r'))
		{
		if (ch == EOF) goto end_of_file;
		if (checkVsHeader)
			{ if (ch != _seq->trueHeader[headerIx]) goto third_line_mismatch; }
		headerIx++;
		ch = seq_getc (_seq);
		}

	if (checkVsHeader)
		{ if (headerIx != strlen(_seq->trueHeader)) goto third_line_short; }

end_of_third_line:

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	//////////
	// read qualities and keep those that are within our index limits
	//////////

	index     = 0;
	qualCount = 0;

	ch = seq_getc (_seq);
	while (ch != EOF)
		{
		if ((ch == '\n') || (ch == '\r')) goto end_of_qualities;
		if ((ch < minFastqCh) || (ch > maxFastqCh)) goto bad_quality;

		// this is a quality character, do we want it?

		index++;
		if ((startLimit != 0) && (index < startLimit)) goto next_quality;
		if ((endLimit   != 0) && (index > endLimit))   goto next_quality;

		// we want it;  let's store it;  note that we don't need to check
		// whether the array is long enough, since vq[] has been allocated in
		// lock step with v[]

		qualCount++;
		if (qualCount > nucCount) goto too_many_qualities;
		_seq->vq[qualLen++] = ch;

		// go try the next quality character

	next_quality:
		ch = seq_getc (_seq);
		}

	// we've read all the qualities, finish up the sequence info

end_of_qualities:

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	if (qualCount < nucCount) goto not_enough_qualities;
	_seq->vq[qualLen] = 0;				// (set the terminating zero)

	return nucCount;

	// failure exits
	// $$$ report line number here

bad_nucleotide:
	description = char_to_description (ch);
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("bad fastq nucleotide character in %s (%s)\n"
		          "remove or replace non-ACGTN characters or consider using --ambiguous=iupac",
		          sequence_filename(_seq), description);
	else
		suicidef ("bad fastq nucleotide character in %s, %s (%s)\n"
		          "remove or replace non-ACGTN characters or consider using --ambiguous=iupac",
				  sequence_filename(_seq), _seq->header, description);

bad_quality:
	description = char_to_description (ch);
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("bad fastq quality character in %s (%s)\n",
		          sequence_filename(_seq), description);
	else
		suicidef ("bad fastq quality character in %s, %s (%s)\n",
				  sequence_filename(_seq), _seq->header, description);

not_enough_qualities:
	description = char_to_description (ch);
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("not enough fastq quality characters in %s\n"
		          unsposFmt " nucleotides and only " unsposFmt " quality characters\n"
		          "(this may be a line-wrapped FASTQ file, which is not supported)",
		          sequence_filename(_seq), nucCount, qualCount);
	else
		suicidef ("not enough fastq quality characters in %s, %s\n"
		          unsposFmt " nucleotides and only " unsposFmt " quality characters\n"
		          "(this may be a line-wrapped FASTQ file, which is not supported)",
				  sequence_filename(_seq), _seq->header, nucCount, qualCount);

too_many_qualities:
	description = char_to_description (ch);
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("too many fastq quality characters in %s\n"
		          unsposFmt " nucleotides and at least " unsposFmt " quality characters\n",
		          sequence_filename(_seq),nucCount,qualCount);
	else
		suicidef ("too many fastq quality characters in %s, %s\n"
		          unsposFmt " nucleotides and at least " unsposFmt " quality characters\n",
				  sequence_filename(_seq), _seq->header, nucCount, qualCount);

bad_fastq_third_line:
	description = char_to_description (ch);
	suicidef ("bad fastq third line character in %s (expected \"+\" but read \"%s\")\n"
	          "(this may be a line-wrapped FASTQ file, which is not supported)",
	          sequence_filename(_seq), description);

third_line_mismatch:
	if (headerIx >= strlen(_seq->trueHeader)) goto third_line_long;
	description = char_to_description (ch);
	suicidef ("fastq third line mismatch in %s (character %d is \"%s\")\n(expected \"+%s\")\n",
	          sequence_filename(_seq), headerIx+2, description, _seq->trueHeader);

third_line_long:
	suicidef ("fastq third line mismatch in %s (line has more than %d characters)\n(expected \"+%s\")\n",
	          sequence_filename(_seq), strlen(_seq->trueHeader)+1, _seq->trueHeader);

third_line_short:
	suicidef ("fastq third line mismatch in %s (line has only %d characters)\n(expected \"+%s\")\n",
	          sequence_filename(_seq), headerIx+1, _seq->trueHeader);

end_of_file:
	if ((_seq->header == NULL) || (_seq->header[0] == 0))
		suicidef ("premature end of fastq file %s\n",
		          sequence_filename(_seq));
	else
		suicidef ("premature end of fastq file %s, %s\n",
				  sequence_filename(_seq), _seq->header);

	return 0; // (never gets here)
	}

//----------
//
// fastq_skip_content--
//	Skip over the content of one sequence in a fastq file.  This assumes that
//	the sequence's file is positioned at the start of the line following the
//	header.
//
//	Upon return, the file is positioned at the start of the header for the next
//	contig (or at the end of file).
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence being parsed.
//
// Returns:
//	true if we were successful;  false if we hit the end-of-file before
//	reading all the content.
//
//----------

static int fastq_skip_content
   (seq*	_seq)
	{
	int		linesToSkip;
	int		prevCh, ch;

	linesToSkip = 3;
	prevCh = 0;
	while (linesToSkip > 0)
		{
		ch = seq_getc (_seq);
		if (ch == EOF) return false;
		if ((ch == '\n') && (prevCh == '\r'))
			{ prevCh = 0;  continue; }
		if ((ch == '\n') || (ch == '\r')) linesToSkip--;
		prevCh = ch;
		}

	return true;
	}

//----------
//
// load_csfasta_sequence--
//	Load the next fasta color sequence from the associated file.
//
// A typical file looks like this:
//
//	# Wed Apr 22 15:07:58 2009 ...
//	>538_743_229_F7
//	T013131021212033022020113200231003030002
//	>538_4021_559_F7
//	T002120310210323111000110101233231231210
//	>534_6488_139_F7
//	T112211320333111020130303120302210313113
//
// Line beginning with '#' are comments and are ignored, but they can only
// occur immediately in front of a header line. Lines beginning with ">" are
// header lines.  If the file contains only one sequence, the header line is
// optional.  Sequences must begin with a nucleotide and thereafter consist
// only of the digits '0', '1', '2' and '3'.  Sequences may occupy multiple
// lines.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_csfasta_sequence
   (seq*	_seq,
	int		keeper)
	{
	int		ch;
	unspos	length;

	if (_seq == NULL) suicide ("load_csfasta_sequence(NULL)");

	//////////
	// read the header
	//////////

	while (true) // skip comment lines
		{
		ch = skip_seq_whitespace (_seq);
		if (ch != '#') break;

		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			ch = seq_getc (_seq);
		}

	if (ch != '>') seq_ungetc (ch, _seq);
	          else parse_csfasta_header (_seq);

	//////////
	// read ahead to determine the length of the sequence and to pre-allocate
	// the vector
	//////////

	if (_seq->rewindable == -1)
		_seq->rewindable = (test_rewindability (_seq) == 0);

	if ((_seq->rewindable == true) && (keeper))
		{
		// read ahead, counting chars needed

		save_fstate (_seq);
		length = parse_csfasta (_seq, /*storeEm*/ false);
		restore_fstate (_seq);

		// allocate the vector

		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_csfasta_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));

		sequence_long_enough (_seq, _seq->len+length, false);
		}

	//////////
	// read the sequence
	//////////

	parse_csfasta (_seq, /*storeEm*/ keeper);
	}

//----------
//
// parse_csfasta_header--
//	Parse a csfasta header from the associated file.  This assumes that the
//	sequence's file is positioned at the first character in the header, *after*
//	the '>' character.
//
//	Upon return, the file is positioned at the start of the first line
//	following the header.
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence being parsed.
//
// Returns:
//  (nothing; the header is written to _seq->header)
//
//----------

static void parse_csfasta_header
   (seq*	_seq)
	{
	u32		headerLen;
	int		ch;

	if (_seq->lockedHeader)
		{
		ch = seq_getc (_seq);
		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			ch = seq_getc (_seq);
		}
	else
		{
		headerLen = 0;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		ch = seq_getc (_seq);

		while ((ch != '\n') && (ch != '\r') && (ch != EOF))
			{
			append_char (&_seq->header, &_seq->headerSize, &headerLen, ch);
			ch = seq_getc (_seq);
			}
		append_char (&_seq->header, &_seq->headerSize, &headerLen, 0);

		if (_seq->nameTrigger != NULL)
			{
			char* triggerFound, *src, *dst;
			triggerFound = strstr (_seq->header, _seq->nameTrigger);
			if (triggerFound != NULL)
				{
				triggerFound += strlen (_seq->nameTrigger);
				for (src=triggerFound,dst=_seq->header ; *src!=0 ; )
					{
					ch = *(src++);
					if ((!isalnum(ch)) && (ch != '_'))
						break;
					*(dst++) = ch;
					}
				*dst = 0;
				}
			}
		}

	if (ch == '\r') // handle possible DOS CR-LF line ending
		{
		ch = seq_getc (_seq);
		if (ch != '\n') seq_ungetc (ch, _seq);
		}

	}

//----------
//
// parse_csfasta--
//	Parse a csfasta sequence from the associated file.  This assumes that the
//	sequence's file is positioned at the first character in the sequence,
//	*after* the sequence header line.
//
// (see load_csfasta_sequence() for info about the file format)
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to parse.
//	int		storeEm:	true  => store the results (in the sequence)
//						false => just count
//
// Returns:
//  The number of characters read into the sequence;  failure causes program
//	fatality.
//
//----------
//
//	Notes:
//
//	(1)	Unlike fasta and fastq, we do not allow csfasta to include separator
//		characters.  The reasoning is that following a separator we would not
//		have the primer nucleotide that we have at the start of the sequence.
//
//----------

static unspos parse_csfasta
   (seq*	_seq,
	int		storeEm)
	{
	unspos	index, startLimit, endLimit, count;
	int		prevCh, ch;
	u8		chType;

	index      = 0;
	count      = 0;
	startLimit = _seq->startLimit;
	endLimit   = _seq->endLimit;

	// scan the file, keeping characters that are (a) colors (or an initial
	// nucleotide) and (b) are within our index limits

	prevCh = '\n';
	ch     = skip_seq_whitespace (_seq);
	while (ch != EOF)
		{
		if ((prevCh == '\n')
		 && ((ch == '#') || (ch == '>'))) 			// (start of next sequence)
			{ seq_ungetc (ch, _seq);  break; }

		chType = char_to_csfasta_type[(u8)ch];
		switch (chType)
			{
			case _nucleotide:
			case _color:
				break;
			case _newline:
				ch = '\n'; // (allow for unix, mac, or pc line ends)
				goto next_char;
			case _bad:
				goto bad_char;
			}

		// this is a color or nucleotide, do we want it?

		if ((index == 0) != (chType == _nucleotide))
			{
			if (index == 0) goto bad_nucleotide;
			           else goto bad_color;
			}

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// we want it;  are we just counting?

		if ((!storeEm) && (count+1 < count))
			suicidef ("in parse_csfasta, "
			          "sequence length " unsposFmt "+1 overflows internal data type",
			          count);

		count++;
		if (!storeEm) goto next_char;

		// ok, let's store it

		if (_seq->len > maxSequenceLen - 1)
			suicidef ("in parse_csfasta, "
			          "sequence length " unsposFmt "+1 exceeds maximum (" unsposFmt ")",
			          _seq->len, maxSequenceLen);

		sequence_long_enough (_seq, _seq->len+1, true);
		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		ch = skip_seq_whitespace (_seq);
		}

	if (storeEm)
		{
		_seq->v[_seq->len] = 0;				// (set the terminating zero)
		_seq->trueLen += index;				// (account for the characters
											//  .. we've read so far)
		}

	// make sure we got somethin' useful

	if ((startLimit != 0) && (startLimit > index))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), startLimit, index);

	if ((endLimit != 0) && (endLimit > index))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), endLimit, index);
		}

	if ((count == 0) && (storeEm))
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		}

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

// (no longer needed;  the loop above exits only at end-of-file or end-of-seq)
//	// skip to the next sequence
//
//	if ((storeEm) && (_seq->needTrueLen))
//		{
//		prevCh = '\n';
//		ch     = skip_seq_whitespace (_seq);
//		while (ch != EOF)
//			{
//			if ((prevCh == '\n')
//			 && ((ch == '#') || (ch == '>')))		// (start of next sequence)
//				{ seq_ungetc (ch, _seq);  break; }
//
//			switch (char_to_csfasta_type[(u8)ch])
//				{
//				case _nucleotide:
//					goto bad_color;
//				case _color:
//					_seq->trueLen++;
//					break;
//				case _newline:
//					ch = '\n'; // (allow for unix, mac, or pc line ends)
//					break;
//				case _bad:
//					goto bad_char;
//				}
//
//			// go try the next character
//
//			prevCh = ch;
//			ch = skip_seq_whitespace (_seq);
//			}
//		}

	return count;

	// failure exits
	// $$$ report line number and sequence name here

bad_char:
	if (dna_isprint(ch))
		suicidef ("bad csfasta character in %s: %c",
		          sequence_filename(_seq), (int) ch);
	else
		suicidef ("bad csfasta character in %s (ascii %02X)",
		          sequence_filename(_seq), (u8) ch);
	return 0; // (never gets here)

bad_nucleotide:
	if (dna_isprint(ch))
		suicidef ("bad csfasta nucleotide in %s: %c",
		          sequence_filename(_seq), (u8) ch);
	else
		suicidef ("bad csfasta nucleotide in %s (ascii %02X)",
		          sequence_filename(_seq), (int) ch);
	return 0; // (never gets here)

bad_color:
	if (dna_isprint(ch))
		suicidef ("bad csfasta color in %s: %c",
		          sequence_filename(_seq), (int) ch);
	else
		suicidef ("bad csfasta color in %s (ascii %02X)",
		          sequence_filename(_seq), (u8) ch);
	return 0; // (never gets here)
	}

//----------
//
// load_nib_sequence--
//	Load a nib sequence from the associated file.
//
// A nib file stores each nucleotide in four bits (one nybble).  The file
// consists of a 4 byte magic number, followed by a 4 byte length, followed by
// the nucleotides.  The magic number is in the file as
//	(first byte) 3A 3D E9 6B (third byte)
// The length field is in little-endian order, so
//	(first byte) C0 E1 E4 00 (third byte)
// means 0x00E4E1C0 bytes (15 million).  The length is the number of
// nucleotides.  The first nucleotide is in the most significant nybble of the
// 9th byte, the second one is in the least significant nybble, the third in
// the 10th byte (msnybble), and so on.  Nybble bits are mapped to characters
// as per the tables nibTo1stChar[] and nibTo2ndChar[].
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void load_nib_sequence
   (seq*		_seq,
	int			keeper)
	{
	u32			magic, length;
	unspos		newSeqLen, ix;
	u32			startLimit, endLimit, startIndex;
	u32			bytesLeft, bytesToRead, bytesRead;
	u8			ch;
	const u8*	to1stChar, *to2ndChar;
	u8*			dst;

	if (!keeper) return;
	if (_seq == NULL) suicide ("load_nib_sequence(NULL)");

	//////////
	// get the sequence length
	//////////

	// check the magic number and decide if it's little or big endian

	magic = read_4_big (_seq);

	// read the length

	if (magic == nibMagicLittle)
		length = read_4_little (_seq);
	else if (magic == nibMagicBig)
		length = read_4_big (_seq);
	else
		{
		length = 0; // (placate compiler)
		suicidef ("bad nib magic number in %s (%08lX)",
		          sequence_filename(_seq), magic);
		}

	if ((length == 0) || (((s32) length) == -1))
		suicidef ("bad nib length in %s (%08lX)", sequence_filename(_seq), length);

	// validate sequence limits

	if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) length))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->startLimit, length);

	if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) length))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->endLimit, length);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	_seq->trueLen += length;

	// skip ahead to the first desired base, and determine how many bases
	// we'll read

	if (startLimit == 0) startLimit = 1;
	startIndex = startLimit - 1;

	bytesLeft = length;
	if (startIndex > 0)
		{ skip_chars (_seq, startIndex/2);  bytesLeft -= 2*(startIndex/2); }

	length = bytesLeft;
	if ((startIndex&1) != 0)	// start offset is odd
		length--;
	if (endLimit != 0)
		{
		if (length > endLimit - startIndex)
			length = endLimit - startIndex;
		}

	//////////
	// allocate the vector, including an extra byte since we may overshoot by
	// 1 when unpacking
	//////////

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
		suicidef ("in load_nib_sequence for %s, "
		          "sequence length %s+%s exceeds maximum (%s)",
		          sequence_filename(_seq),
		          commatize(_seq->len), commatize(length),
		          commatize(maxSequenceLen));
#endif

	newSeqLen = _seq->len + length;
	sequence_long_enough (_seq, newSeqLen+1, false);

	//////////
	// read the sequence
	//////////

	// decide which lookup tables we'll use

	if (_seq->doUnmask)
		{
		to1stChar = nibTo1stCharUnmasked;
		to2ndChar = nibTo2ndCharUnmasked;
		}
	else
		{
		to1stChar = nibTo1stChar;
		to2ndChar = nibTo2ndChar;
		}

	// read the first, partial, byte

	ix = _seq->len;
	if ((startIndex&1) != 0)	// start offset is odd
		{
		ch = seq_getc (_seq);  bytesLeft -= 2;
		_seq->v[ix++] = to2ndChar[ch];
		}

	// process any bytes in the pending buffer, one at a time

	while ((ix < newSeqLen) && (_seq->pendingLen > 0))
		{
		ch = seq_getc (_seq);  bytesLeft -= 2;
		_seq->v[ix++] = to1stChar[ch];
		_seq->v[ix++] = to2ndChar[ch];
		}

	// read the remaining bytes to the tail end of the buffer

	bytesToRead = ((newSeqLen-ix) + 1) / 2;
	dst = _seq->v + _seq->size - bytesToRead;
	if (bytesToRead > 0)
		{
		bytesRead = fread (dst, 1, bytesToRead, _seq->f);
		if (bytesRead != bytesToRead)
			suicidef ("in load_nib_sequence(%s), block read\n"
			          "wanted %d bytes, only got %d",
			          sequence_filename(_seq), bytesToRead, bytesRead);
		}

	// unpack those bytes;  note that although we are writing into the same
	// buffer that we are reading from, and writing two bytes for each one read,
	// the write pointer will not overtake the read pointer;  further note that
	// we may unpack an extra nybble (this will be overwritten when we set the
	// terminating zero)

	while (ix < newSeqLen)
		{
		ch = *(dst++);  bytesLeft--;
		_seq->v[ix++] = to1stChar[ch];
		_seq->v[ix++] = to2ndChar[ch];
		}

	_seq->v[newSeqLen] = 0;					// (set the terminating zero)
	_seq->len = newSeqLen;

	skip_chars (_seq, bytesLeft);			// skip the rest of the sequence
	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

	//////////
	// create a header
	//////////

	if (!_seq->lockedHeader)
		{
		length = snprintf (_seq->header, 0, "%s:" unsposDashFmt,
		                   sequence_filename(_seq),
		                   _seq->startLoc,
		                   _seq->startLoc + _seq->len-1);

		if (_seq->headerSize < length+1)
			{
			_seq->header      = realloc_or_die ("load_nib_sequence (header)",
			                                   _seq->header, length+1);
			_seq->headerSize  = length+1;
			}
		_seq->headerOwner = _seq->shortHeaderOwner = true;

		snprintf (_seq->header, length+1, "%s:" unsposDashFmt,
		          sequence_filename(_seq),
		          _seq->startLoc,
		          _seq->startLoc + _seq->len-1);
		}
	}

//----------
//
// read_2bit_header, load_2bit_sequence--
//	Load a 2bit sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	(load_2bit_sequence only)
//					true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static int find_2bit_sequence    (seq* _seq, char* name);
static u32 read_2bit_index_entry (seq* _seq, char seqName[256], u32 seqNum);

//--- read_2bit_header ---

static void read_2bit_header
   (seq*	_seq)
	{
	u32		magic, version, reserved;

	// read and validate the header

	magic = read_4_big (_seq);

	if (magic == twobitMagicLittle)
		_seq->twoBit.bigEndian = false;
	else if (magic == twobitMagicBig)
		_seq->twoBit.bigEndian = true;
	else
		suicidef ("bad 2bit magic number in %s (%08lX)",
				  sequence_filename(_seq), magic);

	version                 = read_4 (_seq, _seq->twoBit.bigEndian);
	_seq->twoBit.numContigs = read_4 (_seq, _seq->twoBit.bigEndian);
	reserved                = read_4 (_seq, _seq->twoBit.bigEndian);

	if (version != 0)
		suicidef ("bad 2bit version in %s (%08lX)",
		          sequence_filename(_seq), version);
	if (reserved != 0)
		suicidef ("bad 2bit header word 4 in %s (%08lX)",
		          sequence_filename(_seq), reserved);

	if (_seq->twoBit.numContigs == 0)
		suicidef ("empty 2bit file %s", sequence_filename(_seq));

	// save index's file position

	_seq->twoBit.indexFilePos = _seq->twoBit.contigFilePos = ftell (_seq->f);

	// if we have a single contig-of-interest, locate it

	if (_seq->contigOfInterest != NULL)
		{
		if (!find_2bit_sequence (_seq, _seq->contigOfInterest))
			suicidef ("2bit file %s doesn't contain %s",
					  sequence_filename(_seq), _seq->contigOfInterest);
		}

	}

//--- load_2bit_sequence ---

static void load_2bit_sequence
   (seq*	_seq,
	int		keeper)
	{
	char	seqName[maxSequenceName+1];
	int		numChars;
	u32		dnaSize, reserved;
	u32		nBlockCount, maskBlockCount;
	u32		seqDataPos, seekPos;
	unspos	oldSeqLen, ix;
	u32		startLimit, endLimit, length, startIndex, endIndex;
	u32		basesToGo, bytesToSkip, bytesToRead, bytesRead;
	u32		blockIx, s, e, scanIx;
	u8		ch;
	u8*		data, *dst;
	char*	seekType;
	int		err;

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	//////////
	// read the sequence's index table entry
	//////////

	err = fseek (_seq->f, _seq->twoBit.contigFilePos, SEEK_SET);
	if (err != 0)
		{ seekType = "index";  seekPos = _seq->twoBit.contigFilePos;  goto fseek_failed; }

	seqDataPos = read_2bit_index_entry (_seq, seqName, _seq->contig);
	_seq->twoBit.contigFilePos = ftell (_seq->f);

	if (!keeper) return;

	// copy the sequence name as our header (unless the header is locked)

	if (!_seq->lockedHeader)
		{
		numChars = strlen (seqName);
		if (_seq->headerSize < (unsigned) (numChars+1))
			{
			_seq->header      = realloc_or_die ("load_2bit_sequence (header)",
											  _seq->header, numChars+1);
			_seq->headerSize  = numChars+1;
			_seq->headerOwner = _seq->shortHeaderOwner = true;
			}

		strcpy (/*to*/ _seq->header, /*from*/ seqName);
		}

	//////////
	// make sure we have enough room for the sequence's data
	//////////

	err = fseek (_seq->f, seqDataPos, SEEK_SET);
	if (err != 0)
		{ seekType = "header data";  seekPos = seqDataPos;  goto fseek_failed; }

	dnaSize = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if ((dnaSize == 0) || (((s32) dnaSize) == -1))
		suicidef ("bad 2bit length in %s (%08lX)", sequence_filename(_seq), dnaSize);

	if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) dnaSize))
		suicidef ("beyond end in %s (%ld > %ld)",
		          sequence_filename(_seq), _seq->startLimit, dnaSize);

	if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) dnaSize))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->endLimit, dnaSize);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	if (endLimit   == 0) endLimit   = dnaSize;

	_seq->trueLen += dnaSize;

	// allocate the vector; we ask for an additional three characters because
	// during unpacking we write 4 bytes at a time, and thus may overshoot the
	// end by 3
	// $$$ note that this may be way more than needed, if start and end limits
	// $$$ .. reduce the number of characters we actually want;  this could be
	// $$$ .. improved if it causes problems

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((dnaSize+3 > maxSequenceLen) || (_seq->len > maxSequenceLen - (dnaSize+3)))
		goto sequence_too_big;
#endif

	sequence_long_enough (_seq, _seq->len + dnaSize+3, false);

	//////////
	// read and save the intervening block-marking fields
	//////////

	// make sure we have enough room for the n-blocks

	nBlockCount = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if (nBlockCount > _seq->twoBit.nBlocksSize)
		{
		_seq->twoBit.nBlockStarts = (u32*) realloc_or_die ("nBlockStarts", _seq->twoBit.nBlockStarts, nBlockCount * sizeof(u32));
		_seq->twoBit.nBlockSizes  = (u32*) realloc_or_die ("nBlockSizes",  _seq->twoBit.nBlockSizes,  nBlockCount * sizeof(u32));
		_seq->twoBit.nBlocksSize  = nBlockCount;
		}

	// read the n-blocks

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		_seq->twoBit.nBlockStarts[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		_seq->twoBit.nBlockSizes[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	seqDataPos += 4 * (2 * nBlockCount);

	// make sure we have enough room for the mask-blocks;  note that if we are
	// unmasking the sequence, then we skip over the mask-blocks

	maskBlockCount = read_4 (_seq, _seq->twoBit.bigEndian);
	seqDataPos += 4;

	if (_seq->doUnmask)
		{
		seqDataPos += 4 * (2 * maskBlockCount);
		err = fseek (_seq->f, seqDataPos, SEEK_SET);
		if (err != 0)
			{ seekType = "mask data";  seekPos = seqDataPos;  goto fseek_failed; }
		maskBlockCount = 0;
		}
	else
		{
		if (maskBlockCount > _seq->twoBit.mBlocksSize)
			{
			_seq->twoBit.mBlockstarts = (u32*) realloc_or_die ("mBlockstarts", _seq->twoBit.mBlockstarts, maskBlockCount * sizeof(u32));
			_seq->twoBit.mBlocksizes  = (u32*) realloc_or_die ("mBlocksizes",  _seq->twoBit.mBlocksizes,  maskBlockCount * sizeof(u32));
			_seq->twoBit.mBlocksSize  = maskBlockCount;
			}
		}

	// read the mask-blocks

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		_seq->twoBit.mBlockstarts[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		_seq->twoBit.mBlocksizes[blockIx] = read_4 (_seq, _seq->twoBit.bigEndian);

	seqDataPos += 4 * (2 * maskBlockCount);

	// skip the reserved data prefix

	reserved = read_4 (_seq, _seq->twoBit.bigEndian);
	if (reserved != 0)
		suicidef ("bad 2bit reserved data prefix in %s\n"
		          "         (data at %08lX is %08lX)",
		          sequence_filename(_seq), seqDataPos, reserved);

	//////////
	// read the sequence's data
	//////////

	// skip to the first byte containing data of interest

	startIndex = startLimit-1;
	length     = basesToGo = endLimit+1 - startLimit;

	bytesToSkip = startIndex / 4;
	if (bytesToSkip != 0)
		{
		err = fseek (_seq->f, bytesToSkip, SEEK_CUR);
		if (err != 0)
			{ seekType = "data";  seekPos = _seq->twoBit.contigFilePos;  goto fseek_failed; }
		startIndex -= 4*bytesToSkip;
		}

	// read the leading partial byte (if any)

	ix = oldSeqLen = _seq->len;
	if (startIndex > 0)
		{
		ch = (u8) seq_getc (_seq);
		data = (u8*) (twobitToChars[ch] + startIndex);
		while (*data != 0)
			{
			if (basesToGo-- <= 0) break;
			_seq->v[ix++] = *(data++);
			}
		}

	// process any bytes in the pending buffer;  note that we may end up writing
	// as many as 3 bytes beyond the end of the sequence, but will correct this
	// when we write the terminating zero;  also note that we have to separate
	// the last iteration of this loop, since basesToGo is unsigned

	for ( ; (basesToGo>=4)&&(_seq->pendingLen>0) ; basesToGo-=4)
		{
		ch   = (u8)  seq_getc (_seq);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		}

	if ((basesToGo > 0) && (_seq->pendingLen > 0))
		{
		ch   = (u8)  seq_getc (_seq);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		basesToGo = 0;
		}

	// read the remaining bytes to the tail end of the buffer

	bytesToRead = (basesToGo + 3) / 4;
	dst = _seq->v + _seq->size - bytesToRead;
	if (bytesToRead > 0)
		{
		bytesRead = fread (dst, 1, bytesToRead, _seq->f);
		if (bytesRead != bytesToRead) goto read_failed;
		}

	// unpack those bytes;  note that although we are writing into the same
	// buffer that we are reading from, and writing four bytes for each one
	// read, the write pointer will not overtake the read pointer;  as above,
	// we may end up writing as many as 3 bytes beyond the end of the sequence

	for ( ; basesToGo>=4 ; basesToGo-=4)
		{
		ch   = *(dst++);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		}

	if (basesToGo > 0)
		{
		ch   = *(dst++);
		data = (u8*) twobitToChars[ch];
		_seq->v[ix++] = data[0];
		_seq->v[ix++] = data[1];
		_seq->v[ix++] = data[2];
		_seq->v[ix++] = data[3];
		basesToGo = 0;
		}

	_seq->len += length;
	_seq->v[_seq->len] = 0;					// (set the terminating zero)

	//////////
	// mark the Ns and masked bases
	//////////

	startIndex = startLimit-1;
	endIndex   = endLimit;

	for (blockIx=0 ; blockIx<nBlockCount ; blockIx++)
		{
		s = _seq->twoBit.nBlockStarts[blockIx];
		e = s + _seq->twoBit.nBlockSizes[blockIx];
		if (e <= startIndex) continue;
		if (s >= endIndex)   continue;
		if (s <  startIndex) s = startIndex;
		if (e >  endIndex)   e = endIndex;
		s -= startIndex;
		e -= startIndex;
		for (scanIx=s ; scanIx<e ; scanIx++)
			_seq->v[oldSeqLen+scanIx] = 'N';
		}

	for (blockIx=0 ; blockIx<maskBlockCount ; blockIx++)
		{
		s = _seq->twoBit.mBlockstarts[blockIx];
		e = s + _seq->twoBit.mBlocksizes[blockIx];
		if (e <= startIndex) continue;
		if (s >= endIndex)   continue;
		if (s <  startIndex) s = startIndex;
		if (e >  endIndex)   e = endIndex;
		s -= startIndex;
		e -= startIndex;
		for (scanIx=s ; scanIx<e ; scanIx++)
			_seq->v[oldSeqLen+scanIx] = dna_tolower (_seq->v[oldSeqLen+scanIx]);
		}

	_seq->twoBit.contigLoaded = true;

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;
	return;

// failure exits

fseek_failed:
	suicidef ("failed to seek to position in \"%s\"\n"
	          "in load_2bit_sequence, %s fseek(%08lX) returned %d",
	          sequence_filename(_seq), seekType, seekPos, err);
	return; // (never gets here)

sequence_too_big:
	suicidef ("in load_2bit_sequence for %s, "
			  "sequence length %s+%s exceeds maximum (%s)",
			  sequence_filename(_seq),
			  commatize(_seq->len),commatize(dnaSize+3),
			  commatize(maxSequenceLen));
	return; // (never gets here)

read_failed:
	suicidef ("in load_2bit_sequence for %s,"
	          " block read for sequence %u\n"
	          "wanted %d bytes, only got %d",
	          sequence_filename(_seq), _seq->contig, bytesToRead, bytesRead);
	return; // (never gets here)
	}

//--- find_2bit_sequence ---

static int find_2bit_sequence
   (seq*	_seq,
	char*	name)
	{
	char	seqName[maxSequenceName+1];
	u32		ix;

	for (ix=0 ; ix<_seq->twoBit.numContigs ; ix++)
		{
		_seq->twoBit.contigFilePos = ftell (_seq->f);
		read_2bit_index_entry (_seq, seqName, ix+1);
		if (strcmp (seqName, name) == 0)
			{ _seq->contig = ix + 1;  return true; }
		}

	return false;
	}

//--- read_2bit_index_entry ---

static u32 read_2bit_index_entry
   (seq*			_seq,
	char			seqName[maxSequenceName+1],
	u32				seqNum)
	{
	unsigned int	nameSize;
	size_t			bytesRead;

	// read the name

	nameSize = getc_or_die (_seq->f, _seq->filename);
	if (nameSize > 0)
		{
		bytesRead = fread (seqName, 1, nameSize, _seq->f);
		if (bytesRead != nameSize)
			suicidef ("in load_2bit_sequence for %s, short read for sequence %u\n"
			          "wanted %d bytes, only got %d",
			          sequence_filename(_seq), seqNum, nameSize, bytesRead);
		}
	seqName[nameSize] = 0;

	// read the data offset

	return read_4 (_seq, _seq->twoBit.bigEndian);
	}

//----------
//
// read_hsx_header, load_hsx_sequence--
//	Load a sequence from the associated hsx file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	(load_hsx_sequence only)
//					true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static u64   lookup_hsx_sequence  (seq* _seq, char* name);
static u64   find_hsx_sequence    (seq* _seq, char* name,
                                   u64 bucketStart, u64 bucketEnd);
static char* read_hsx_index_entry (seq* _seq);
static char* read_hsx_string      (seq* _seq, FILE* f);

//--- read_hsx_header ---

static void read_hsx_header
   (seq*		_seq)
	{
	u32			fileInfoOffset[255];
	u32			magic, headerLength;
	char*		s;
	char		extension[10];
	u32			infoBytes, nameBytes;
	char*		slash, *dot, *nameScan;
	int  		baseLen, pathLen;
	u32			fileNum;
	int			err;

	// read and validate the header

	magic = read_4_big (_seq);

	if (magic == hsxMagicLittle)
		_seq->hsx.bigEndian = false;
	else if (magic == hsxMagicBig)
		_seq->hsx.bigEndian = true;
	else
		suicidef ("bad hsx magic number in %s (%08lX)",
				  sequence_filename(_seq), magic);

	_seq->hsx.version = read_4 (_seq, _seq->hsx.bigEndian);
	if (_seq->hsx.version != 0x00000100L)
		suicidef ("bad hsx version in %s (%08lX)",
		          sequence_filename(_seq), _seq->hsx.version);

	headerLength = read_4 (_seq, _seq->hsx.bigEndian);
	if (headerLength != 0x1C)
		suicidef ("bad hsx header length in %s (%08lX)",
		          sequence_filename(_seq), headerLength);

	_seq->hsx.numFiles        =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.fileTableOffset = (u64) read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.numBuckets      =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.hashTableOffset = (u64) read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.numContigs      =       read_4 (_seq, _seq->hsx.bigEndian);
	_seq->hsx.seqTableOffset  = (u64) read_4 (_seq, _seq->hsx.bigEndian);

	if (_seq->hsx.numFiles == 0)
		suicidef ("empty file table in hsx file %s", sequence_filename(_seq));

	if (_seq->hsx.numFiles > 255)
		suicidef ("corrupt header in hsx file %s (numFiles > 255;  %d)",
		          sequence_filename(_seq), _seq->hsx.numFiles);

	if (_seq->hsx.numBuckets == 0)
		suicidef ("corrupt header in hsx file %s (numBuckets = 0)",
		          sequence_filename(_seq));

	// read and validate the file table

	err = fseek (_seq->f, (long int) _seq->hsx.fileTableOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
		          sequence_filename(_seq), _seq->hsx.fileTableOffset, err);

	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		fileInfoOffset[fileNum] = (u64) read_4 (_seq, _seq->hsx.bigEndian);

	slash = strrchr (_seq->filename, pathSlash);
	dot   = strrchr (_seq->filename, '.');
	if ((dot == NULL) || ((slash != NULL) && (dot < slash)))
		baseLen = strlen(_seq->filename);
	else
		baseLen = dot - _seq->filename;
	if (slash == NULL)
		pathLen = 0;
	else
		pathLen = slash+1 - _seq->filename;

	infoBytes = sizeof(hsxfileinfo) * _seq->hsx.numFiles;
	nameBytes = 0;
	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		{
		err = fseek (_seq->f, (long int) fileInfoOffset[fileNum], SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
			          sequence_filename(_seq), fileInfoOffset[fileNum], err);

		s = read_hsx_string (_seq, _seq->f);
		if ((strcmp (s, "fa")    != 0)
		 && (strcmp (s, "fasta") != 0))
			suicidef ("in read_hsx_header for %s, unsupported file type: %s",
					  sequence_filename(_seq), s);
		strncpy (/*to*/ extension, /*from*/ s, sizeof(extension));

		s = read_hsx_string (_seq, _seq->f);
		if (s[0] != 0)
			nameBytes += pathLen + strlen(s) + 1 + strlen(extension) + 1;
		else
			nameBytes += baseLen + 1 + strlen(extension) + 1;
		}

	_seq->hsx.fileInfo = (hsxfileinfo*) zalloc_or_die ("read_hsx_header", infoBytes + nameBytes);

	nameScan = ((char*) _seq->hsx.fileInfo) + infoBytes;
	for (fileNum=0 ; fileNum<_seq->hsx.numFiles ; fileNum++)
		{
		_seq->hsx.fileInfo[fileNum].name = nameScan;
		_seq->hsx.fileInfo[fileNum].f    = NULL;

		err = fseek (_seq->f, (long int) fileInfoOffset[fileNum], SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s, file table fseek(%08lX) returned %d",
			          sequence_filename(_seq), fileInfoOffset[fileNum], err);

		s = read_hsx_string (_seq, _seq->f);
		strncpy (/*to*/ extension, /*from*/ s, sizeof(extension));

		s = read_hsx_string (_seq, _seq->f);
		if (s[0] != 0)
			{
			strncpy (/*to*/    nameScan,
					 /*from*/  _seq->filename,
					 /*limit*/ pathLen);
			strcpy  (/*to*/   nameScan + pathLen,
					 /*from*/ s);
			nameScan[pathLen+strlen(s)] = '.';
			strcpy  (/*to*/   nameScan + pathLen+strlen(s) + 1,
					 /*from*/ extension);
			nameScan += pathLen + strlen(s) + 1 + strlen(extension) + 1;
			}
		else
			{
			strncpy (/*to*/    nameScan,
					 /*from*/  _seq->filename,
					 /*limit*/ baseLen);
			nameScan[baseLen] = '.';
			strcpy  (/*to*/   nameScan + baseLen + 1,
					 /*from*/ extension);
			nameScan += baseLen + 1 + strlen(extension) + 1;
			}
		}

	// locate the first contig

	locate_hsx_first_sequence (_seq);
	}

//--- locate_hsx_first_sequence ---

static void locate_hsx_first_sequence
   (seq*		_seq)
	{
	u64			fileOffset;
	long int	bucketOffset;
	int			err;

	// if we have a single contig-of-interest, locate it

	if (_seq->contigOfInterest != NULL)
		{
		fileOffset = lookup_hsx_sequence (_seq, _seq->contigOfInterest);
		if ((fileOffset & hsxMsBit5) != 0)
			suicidef ("hsx file %s doesn't contain %s",
					  sequence_filename(_seq), _seq->contigOfInterest);
		if (fileOffset > hsxMaxFilePos)
			suicidef ("in read_hsx_header for %s,"
			          " file pos for %s (%010lX) exceeds max (%010lX)",
					  sequence_filename(_seq), _seq->contigOfInterest,
					  fileOffset, hsxMaxFilePos);
		_seq->hsx.contigFilePos = fileOffset;
		debugNamesFile_5;
		}

	// otherwise, if we have no contig names or chores, locate the first
	// sequence in the index

	else if ((_seq->namesFilename == NULL) && (_seq->choresFilename == NULL))
		{
		bucketOffset = (long int) _seq->hsx.hashTableOffset;
		err = fseek (_seq->f, bucketOffset, SEEK_SET);
		if (err != 0)
			suicidef ("in read_hsx_header for %s,"
			          " file table fseek(%010lX) returned %d",
					  sequence_filename(_seq), 0, err);

		fileOffset = read_5 (_seq, _seq->hsx.bigEndian) & ~hsxMsBit5;
		if (fileOffset > hsxMaxFilePos)
			suicidef ("in read_hsx_header for %s,"
			          " file pos for index 0 (%010lX) exceeds max (%010lX)",
					  sequence_filename(_seq), fileOffset, hsxMaxFilePos);
		_seq->hsx.contigFilePos = fileOffset;
		debugNamesFile_6;
		}

	}

//--- load_hsx_sequence ---

static void load_hsx_sequence
   (seq*	_seq,
	int		keeper)
	{
	int		err;
	char*	seqName;
	int		numChars;
	unspos	dnaSize;
	unspos	index, startLimit, endLimit;
	int		prevCh, ch;
	char*	seqFName;
	FILE*	seqF;
	char*	description;

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	debugNamesFile_7;

	//////////
	// read the sequence's index table entry
	//////////

	if (_seq->hsx.contigFilePos > hsxMaxFilePos)
		suicidef ("in load_hsx_sequence for %s,"
		          " file pos for contig %u (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), _seq->contig,
				  _seq->hsx.contigFilePos, hsxMaxFilePos);
	err = fseek (_seq->f, (long int) _seq->hsx.contigFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in load_hsx_sequence for %s, index fseek(%010lX) returned %d",
		          sequence_filename(_seq), _seq->hsx.contigFilePos, err);

	seqName = read_hsx_index_entry (_seq);
	_seq->hsx.contigFilePos = (u64) ftell (_seq->f);
	debugNamesFile_8;

	if (!keeper) return;

	if (_seq->hsx.seqLength > (u64) maxSequenceLen)
		suicidef ("in load_hsx_sequence for %s, "
		          "sequence length " unsposFmt " for %s "
		          "exceeds maximum (" unsposFmt ")",
		          sequence_filename(_seq), _seq->hsx.seqLength, seqName,
		          maxSequenceLen);

	// copy the sequence name as our header

	numChars = strlen (seqName);
	if (_seq->headerSize < (unsigned) (numChars+1))
		{
		_seq->header      = realloc_or_die ("load_hsx_sequence (header)",
						 				  _seq->header, numChars+1);
		_seq->headerSize  = numChars+1;
		_seq->headerOwner = _seq->shortHeaderOwner = true;
		}

	strcpy (/*to*/ _seq->header, /*from*/ seqName);

	//////////
	// make sure we have enough room for the sequence's data
	// 
	// $$$ note that the allocated vector may be way more than needed, if start
	// $$$ .. and end limits reduce the number of characters we actually want;
	// $$$ .. this could be improved if it causes problems
	//////////

	dnaSize = (unspos) _seq->hsx.seqLength;

	if ((_seq->startLimit != 0) && (_seq->startLimit > dnaSize))
		suicidef ("beyond end in %s/%s (%ld > " unsposFmt ")",
		          sequence_filename(_seq), seqName, _seq->startLimit, dnaSize);

	if ((_seq->endLimit != 0) && (_seq->endLimit > dnaSize))
		{
		if (_seq->endIsSoft)
			{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
		else
			suicidef ("beyond end in %s/%s (%ld > " unsposFmt ")",
			          sequence_filename(_seq), seqName, _seq->endLimit, dnaSize);
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	if (endLimit   == 0) endLimit   = dnaSize;

	_seq->trueLen += dnaSize;

#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
	if ((dnaSize > maxSequenceLen) || (_seq->len > maxSequenceLen - dnaSize))
		suicidef ("in load_hsx_sequence for %s/%s, "
		          "sequence length %s+%s exceeds maximum (%s)",
		          sequence_filename(_seq), seqName,
		          commatize(_seq->len),commatize(dnaSize),
		          commatize(maxSequenceLen));
#endif

	sequence_long_enough (_seq, _seq->len + dnaSize, false);

	// if the sequence is empty, warn the user but return the empty sequence to
	// our caller

	if (_seq->hsx.seqLength == 0)
		{
		if (_seq->header == NULL)
			fprintf (stderr, "WARNING. %s contains an empty sequence\n",
			                 sequence_filename(_seq));
		else
			fprintf (stderr, "WARNING. %s contains an empty sequence:\n%s\n",
			                 sequence_filename(_seq), _seq->header);
		return;
		}

	//////////
	// read the sequence's data
	//////////

	seqFName = _seq->hsx.fileInfo[_seq->hsx.seqFileIx].name;
	seqF     = _seq->hsx.fileInfo[_seq->hsx.seqFileIx].f;
	if (seqF == NULL)
		{
		// $$$ we should probably keep track of the number of open files and
		// $$$ close some (by LRU) if too many are open
		seqF = fopen_or_die (seqFName, "rb");
		_seq->hsx.fileInfo[_seq->hsx.seqFileIx].f = seqF;
		}

	if (_seq->hsx.seqFilePos > hsxMaxFilePos)
		suicidef ("in load_hsx_sequence for %s/%s,"
		          " file pos for sequence %s (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), seqName, _seq->header,
				  _seq->hsx.seqFilePos, hsxMaxFilePos);
	err = fseek (seqF, _seq->hsx.seqFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in load_hsx_sequence for %s/s,"
		          " data fseek(%s,%08lX) returned %d",
		          sequence_filename(_seq), seqName,
		          seqFName, _seq->hsx.seqFilePos, err);

	// if the first character is a '>' (and the length is non-zero), we have to
	// skip this sequence header
	// $$$ it might be a good idea to validate that the header we are reading
	// $$$ .. matches the name of the sequence we think we're going to be reading

	prevCh = '\n';

	ch = getc_or_die (seqF, seqFName);
	if ((ch == '>') && (dnaSize != 0))
		{
		while (ch != '\n') // (skip line)
			ch = getc_or_die (seqF, seqFName);
		// get first character of next line
		ch = getc_or_die (seqF, seqFName);
		}

	while ((ch == ' ') || (ch == '\t')) // (skip whitespace)
		ch = getc_or_die (seqF, seqFName);

	// scan the file, keeping characters that are (a) nucleotides and (b) are
	// within our index limits

	index = 0;
	while (ch != EOF)
		{
		if ((prevCh == '\n') && (ch == '>')) // (start of next sequence)
			break;

		if ((_seq->separatorCh == 0) || (ch != _seq->separatorCh))
			{
			switch (char_to_fasta_type[(u8)ch])
				{
				case _nucleotide:
					break;
				case _ambiguous:
					if (!_seq->allowAmbiDNA) goto bad_char;
					break;
				case _newline:
					ch = '\n'; // (allow for unix, mac, or pc line ends)
					goto next_char;
				case _bad:
					goto bad_char;
				}
			}

		// this is a nucleotide (or separator), do we want it?

		index++;

		if ((startLimit != 0) && (index < startLimit)) goto next_char;
		if ((endLimit   != 0) && (index > endLimit))   goto next_char;

		// ok, let's store it

		_seq->v[_seq->len++] = ch;

		// go try the next character

	next_char:
		prevCh = ch;
		do // (skip whitespace)
			{
			ch = getc_or_die (seqF, seqFName);
			} while ((ch == ' ') || (ch == '\t'));
		}

	_seq->v[_seq->len] = 0;				// (set the terminating zero)

	// $$$ we should make sure the sequence was as long as it said it was

	_seq->hsx.contigLoaded = true;

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

	return;

	// failure exits
	// $$$ report line number and sequence name here

bad_char:
	description = char_to_description (ch);
	suicidef ("bad fasta character in %s, %s (%s)\n"
	          "remove or replace non-ACGTN characters or consider using --ambiguous=iupac",
	          sequence_filename(_seq), seqName, description);
	} 

//--- lookup_hsx_sequence ---

static u64 lookup_hsx_sequence
   (seq*	_seq,
	char*	name)
	{
	u32		bucket;
	u64		fileOffset;
	u64		bucketStart, bucketEnd;
	int		err;

	bucket = hassock_hash (name, strlen(name)) % _seq->hsx.numBuckets;
	fileOffset = _seq->hsx.hashTableOffset + (5 * (u64) bucket);
	debugNamesFile_9;
	if (fileOffset > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s hash bucket %d (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucket, fileOffset, hsxMaxFilePos);
	err = fseek (_seq->f, (long int) fileOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in lookup_hsx_sequence for %s, file table fseek(%010lX) returned %d",
		          sequence_filename(_seq), fileOffset, err);

	bucketStart = read_5 (_seq, _seq->hsx.bigEndian);
	debugNamesFile_10;
	if ((bucketStart & hsxMsBit5) != 0) // (bucket is empty)
		return hsxMsBit5; // (not found)
	if (bucketStart > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s bucket start (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucketStart, hsxMaxFilePos);

	bucketEnd = read_5 (_seq, _seq->hsx.bigEndian) & ~hsxMsBit5;
	if (bucketEnd > hsxMaxFilePos)
		suicidef ("in lookup_hsx_sequence for %s,"
				  " file pos for %s bucket end (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), bucketEnd, hsxMaxFilePos);

	return find_hsx_sequence (_seq, name, bucketStart, bucketEnd);
	}

//--- find_hsx_sequence ---

static u64 find_hsx_sequence
   (seq*	_seq,
	char*	name,
	u64		bucketStart,
	u64		bucketEnd)
	{
	u64		bucketOffset = bucketStart;
	char*	seqName;
	int		diff, err;

	err = fseek (_seq->f, (unsigned long) bucketOffset, SEEK_SET);
	if (err != 0)
		suicidef ("in find_hsx_sequence for %s,"
		          " file table fseek(%010lX) returned %d",
		          sequence_filename(_seq), bucketOffset, err);

	while (bucketOffset < bucketEnd)
		{
		seqName = read_hsx_index_entry (_seq);
		diff = strcmp (seqName, name);
		if (diff == 0) return bucketOffset; // (sequence name found)
		if (diff >  0) break;               // (sequence name not found)
		bucketOffset += 1 + 6 + 5 + strlen(seqName) + 1;
		}

	return hsxMsBit5; // (not found)
	}

//--- read_hsx_index_entry ---

static char* read_hsx_index_entry
   (seq* _seq)
	{
	_seq->hsx.seqLength  = read_5   (_seq, _seq->hsx.bigEndian);
	_seq->hsx.seqFileIx  = seq_getc (_seq);
	_seq->hsx.seqFilePos = read_6   (_seq, _seq->hsx.bigEndian);
	return read_hsx_string (_seq, _seq->f);
	}

//--- read_hsx_string ---

static char* read_hsx_string
   (seq*			_seq,
	FILE*			f)
	{
	static char		s[256];
	unsigned int	stringSize;
	size_t			bytesRead;

	// read the name

	stringSize = getc_or_die (_seq->f, _seq->filename);
	if (stringSize == 0)
		{ s[0] = 0;  return s; }

	bytesRead = fread (s, 1, stringSize, f);
	if (bytesRead != stringSize)
		suicidef ("in read_hsx_string for %s, short read\n"
				  "wanted %d bytes, only got %d",
				  sequence_filename(_seq), stringSize, bytesRead);

	s[stringSize] = 0;
	return s;
	}

//----------
//
// load_qdna_sequence--
//	Load a quantum-dna sequence from the associated file.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to load.
//	int		keeper:	true  => actually load the sequence
//					false => just skip the sequence
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------
//
// Qdna file format:
//
//	Fields can be in big- or little-endian format;  they must match the
//	endianess of the magic number.
//
//	Version 2 (name ignored; files with named properites not supported):
//
//	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
//	offset 0x04: 00 00 02 00   version 2.0 (fourth byte is sub version)
//	offset 0x08: 00 00 00 14   header length (in bytes, including this field)
//	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
//	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
//	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
//	offset 0x18: xx xx xx xx   (for version >= 2.0) P, offset to named
//							   .. properties, 0 indicates no properties
//	offset    N: ...           name (zero-terminated string)
//	offset    S: ...           data sequence
//	offset    P: ...           named properties (see below)
//
//	The named properties section is not allowed in this implementation.
//
//	Version 1 (name ignored):
//
//	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
//	offset 0x04: 00 00 01 00   version (fourth byte will be sub version)
//	offset 0x08: 00 00 00 10   header length (in bytes, including this field)
//	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
//	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
//	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
//	offset    N:  ...          name (zero-terminated string)
//	offset    S:  ...          data sequence
//
//	Version 0:
//
//	offset 0x00: 9E 65 56 F6   magic number
//	offset 0x04:  ...          data sequence
//
//	Additionally, we will accept any binary file and interpret it as the data
//	sequence.  Note that if the data sequence happens to begin with one of the
//	magic numbers above, we will fail to read the file properly.  Further, if
//	the file contains newlines that are not part of the sequence, we will fail
//	to read the file properly.
//
//----------

// $$$ why don't we use the name from the file?????

static void load_qdna_sequence
   (seq*	_seq,
	int		keeper)
   	{
	u32		magic, version, seqOffset, propOffset;
	u32		length, startLimit, startIndex, endLimit;
	unspos	newSeqLen;
	int		oldFormat, bigEndian, lengthKnown;
	int		ch;
	int		err, numChars;

	if (!keeper) return;
	if (_seq == NULL) suicide ("load_qdna_sequence(NULL)");

	//////////
	// process the header
	//////////

	// validate the magic number

	oldFormat = bigEndian = false;

	magic = read_4_big (_seq);
	if      (magic == qdnaMagicLittle)    { ; }
	else if (magic == qdnaMagicBig)       { bigEndian = true; }
	else if (magic == oldQdnaMagicLittle) { oldFormat = true; }
	else if (magic == oldQdnaMagicBig)    { oldFormat = bigEndian = true; }
	else
		{
		seq_ungetc ((magic >> 24) & 0xFF, _seq);
		seq_ungetc ((magic >> 16) & 0xFF, _seq);
		seq_ungetc ((magic >> 8)  & 0xFF, _seq);
		seq_ungetc ( magic        & 0xFF, _seq);
		oldFormat = true;
		}

	// skip the header (unless it's the old format)

	if (oldFormat)
		{
		lengthKnown = false;
		length      = 0;
		}
	else
		{
		version = read_4 (_seq, bigEndian);
		if (((version >> 8) != 1) && ((version >> 8) != 2))
			suicidef ("unsupported qdna version in %s (%08lX)",
			          sequence_filename(_seq), version);

		/*headerLen=*/  read_4 (_seq, bigEndian);
		seqOffset     = read_4 (_seq, bigEndian);
		/*nameOffset=*/ read_4 (_seq, bigEndian);
		length        = read_4 (_seq, bigEndian);  lengthKnown = true;

		if ((version >> 8) == 1)
			skip_chars (_seq, seqOffset - 0x18);
		if ((version >> 8) == 2)
			{
			propOffset = read_4 (_seq, bigEndian);
			if (propOffset != 0)
				suicidef ("qdna named properties are not supported in %s",
				          sequence_filename(_seq));
			skip_chars (_seq, seqOffset - 0x1C);
			}

		_seq->trueLen += length;
		}

	//////////
	// skip ahead to the first desired base, and try to determine how many
	// bases we'll read
	//////////

	if (lengthKnown)
		{
		if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) length))
			suicidef ("beyond end in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->startLimit, length);

		if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) length))
			{
			if (_seq->endIsSoft)
				{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
			else
				suicidef ("beyond end in %s (%ld > %ld)",
				          sequence_filename(_seq), _seq->endLimit, length);
			}
		}
	else
		{
		if ((_seq->startLimit != 0) && (_seq->startLimit > (unspos) 0xFFFFFFFF))
			suicidef ("invalid start limit in %s (%ld > %ld)",
			          sequence_filename(_seq), _seq->startLimit, 0xFFFFFFFF);

		if ((_seq->endLimit != 0) && (_seq->endLimit > (unspos) 0xFFFFFFFF))
			{
			if (_seq->endIsSoft)
				{ _seq->endLimit = 0;  _seq->endIsSoft = false; }
			else
				suicidef ("invalid end limit in %s (%ld > %ld)",
				          sequence_filename(_seq), _seq->endLimit, 0xFFFFFFFF);
			}
		}

	startLimit = (u32) _seq->startLimit;
	endLimit   = (u32) _seq->endLimit;

	if (startLimit == 0) startLimit = 1;
	startIndex = startLimit - 1;

	if (startIndex > 0)
		{
		if (!skip_chars (_seq, startIndex))
			suicidef ("bad start index for %s: %d",
			          sequence_filename(_seq), startIndex);
		}

	if (endLimit != 0)
		{ length = endLimit - startIndex;  lengthKnown = true; }

	//////////
	// allocate the vector (if we know the length)
	//////////

	newSeqLen = 0;
	if (lengthKnown)
		{
#if (maxSequenceIndex <= 32)	// otherwise compiler complains that this test is
								// .. always false
		if ((length > maxSequenceLen) || (_seq->len > maxSequenceLen - length))
			suicidef ("in load_qdna_sequence for %s, "
			          "sequence length %s+%s exceeds maximum (%s)",
			          sequence_filename(_seq),
			          commatize(_seq->len), commatize(length),
			          commatize(maxSequenceLen));
#endif

		newSeqLen = _seq->len + length;
		sequence_long_enough (_seq, newSeqLen, false);
		}

	//////////
	// read the sequence
	//////////

	while (true)
		{
		if ((newSeqLen != 0) && (_seq->len >= newSeqLen))
			break;

		// read the next character from the sequence

		ch = seq_getc (_seq);
		if (ch == EOF) break;

		if (ch == 0)
			suicidef ("in load_qdna_sequence(), file contains a zero");

		// allocate more room in the vector if we need it, and deposit the
		// character in the sequence

		if (!lengthKnown)
			sequence_long_enough (_seq, _seq->len+1, true);

		_seq->v[_seq->len++] = (u8) ch;
		}

	_seq->v[_seq->len] = 0;

	if (oldFormat)
		_seq->trueLen += _seq->len + startIndex;
	else
		{ ; } // (for new format, we already added the file length to _seq->trueLen

	if ((newSeqLen != 0) && (_seq->len < newSeqLen))
		suicidef ("beyond end in %s (%ld > end of file)",
		          sequence_filename(_seq), endLimit);

	_seq->pendingLen   = 0;					// (discard any pending chars)
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;

	if (startLimit == 0) _seq->startLoc = 1;
	                else _seq->startLoc = startLimit;

	// skip the rest of the sequence

	if ((oldFormat) && (_seq->needTrueLen))
		{
		while ((ch = seq_getc (_seq)) != EOF)
			_seq->trueLen++;
		}
	else
		{
		err = fseek (_seq->f, 0, SEEK_END);	// skip the rest of the sequence
		if (err != 0)
			suicidef_with_perror ("in load_qdna_sequence(), fseek returned %d",
			                      err);
		}

	//////////
	// create a sequence header
	//////////

	if (!_seq->lockedHeader)
		{
		numChars = snprintf (_seq->header, 0, "%s:" unsposDashFmt,
		                     sequence_filename(_seq),
		                     _seq->startLoc,
		                     _seq->startLoc + _seq->len-1);

		if (_seq->headerSize < (unsigned) numChars+1)
			{
			_seq->header      = realloc_or_die ("load_qdna_sequence (header)",
			                                   _seq->header, numChars+1);
			_seq->headerSize  = numChars+1;
			_seq->headerOwner = _seq->shortHeaderOwner = true;
			}

		snprintf (_seq->header, numChars+1, "%s:" unsposDashFmt,
		          sequence_filename(_seq),
		          _seq->startLoc,
		          _seq->startLoc + _seq->len-1);
		}
	}

//----------
//
// another_sequence--
//	Determine if the associated file has another sequence (see note 1 for
//	clarification).
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to dispose of.
//
// Returns:
//	true if there is another sequence to load, false if not.
//
//----------
//
//	Notes:
//
//	(1)	If a chores file is active, "another sequence" means "another chore",
//		even if any remaining chores are for the current sequence.
//
//----------

int another_sequence
   (seq*			_seq)
	{
	debugNamesFile_11;

	if (_seq == NULL)                      suicide ("another_sequence(NULL)");
	if (_seq->fileType == seq_type_nofile) return false;

	// if we're subsampling the file's sequences, skip past sequences as needed

	if (_seq->subsampleN > 0)
		{
		if (_seq->subsampleSkip > 0)
			skip_sequences (_seq, _seq->subsampleSkip);
		_seq->subsampleSkip = 0;
		}

	return another_sequence_core (_seq);
	}

static int another_sequence_core
   (seq*			_seq)
	{
	seqpartition*	sp;
	int				ch;
	char*			header;
	int				haveNamesFile, inhibitSearch;

	// if this is a partitioned sequence and we've finished loading it, there's
	// never another sequence

	sp = &_seq->partition;
	if ((sp->p != NULL) && (sp->state >= seqpart_ready))
		{
		if (_seq->doJoin) return false;
		sp->state = seqpart_reusable;
		}

	// if we've previously positioned the file to a contig, but have yet to
	// load that contig, then we have another sequence

	if (_seq->contigPending) return true;

	// if we pre-loaded a sequence then that sequence counts as "another"
	// sequence since the caller doesn't know we did so

	if (_seq->preLoaded)
		{
		if ((sp->p != NULL) && (sp->state == seqpart_reusable))
			sp->state = seqpart_ready;
		return true;
		}

	// if we have a contigs-of-interest file, get the next contig name ...

	inhibitSearch = false;

	if (_seq->namesFile != NULL)
		{
		if (!read_contig_name (_seq))
			return false;
		}

	// ... or, if we have a chores file, get the next chore

	else if (_seq->choresFile != NULL)
		{
		if (!read_chore (_seq))
			return false;

		header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
		if ((header != NULL) && (strcmp (header, _seq->nextContigName) == 0))
			{
			_seq->preLoaded = true;
			if ((sp->p != NULL) && (sp->state == seqpart_reusable))
				sp->state = seqpart_ready;
			inhibitSearch = true;
			validate_rev_comp (_seq);
			}
		}

	// for 2bit or hsx files it's a matter of whether we're at the end of the
	// index list

	haveNamesFile = (_seq->namesFile != NULL) || (_seq->choresFile != NULL);

	if (_seq->fileType == seq_type_2bit)
		{
		if ((!haveNamesFile)
		 && (!_seq->twoBit.contigLoaded))            return (_seq->twoBit.numContigs > 0);
		if (_seq->contigOfInterest != NULL)          return false;
		if ((haveNamesFile)
		 && (!inhibitSearch))                        return find_next_2bit_coi (_seq);
		if ((_seq->contig >= _seq->twoBit.numContigs)
		 && (!inhibitSearch))                        return false;
		return true;
		}
	else if (_seq->fileType == seq_type_hsx)
		{
		if ((!haveNamesFile)
		 && (!_seq->hsx.contigLoaded))               return (_seq->hsx.numContigs > 0);
		if (_seq->contigOfInterest != NULL)          return false;
		if ((haveNamesFile)
		 && (!inhibitSearch))                        return find_next_hsx_coi (_seq);
		if ((_seq->contig >= _seq->hsx.numContigs)
		 && (!inhibitSearch))                        return false;
		return true;
		}

	// otherwise it's a matter of having data left in the file

	if (_seq->f == NULL)  return false;		// we've have no file to read from
	if (feof   (_seq->f)) return false;		// we've previously hit end of file
	if (ferror (_seq->f)) return false;		// we've previously had a problem

	if ((_seq->fileType == seq_type_fasta)	// we've got the next contig-of-
	 && (haveNamesFile))					// .. interest
		return find_next_fasta_coi (_seq);
	if ((_seq->fileType == seq_type_fastq)   && (haveNamesFile))
		return find_next_fastq_coi (_seq);
	if ((_seq->fileType == seq_type_csfasta) && (haveNamesFile))
		return find_next_csfasta_coi (_seq);

	if (_seq->pendingLen > 0) return true;	// we have characters to process

	ch = getc_or_die (_seq->f,				// take a peek and see what's left
	                  _seq->filename);
	if (ch == EOF) return false;			// we're at end of file now

	seq_ungetc (ch, _seq);					// save what we peeked at
	return true;							// we have characters to process
	}


// find_next_fasta_coi, find_next_csfasta_coi--
//	advance to the next contig-of-interest in fasta or csfasta file
//  (always returns true)

static int find_next_fasta_coi (seq* _seq)
	{ return find_next_general_fasta_coi (_seq, false); }

static int find_next_csfasta_coi (seq* _seq)
	{ return find_next_general_fasta_coi (_seq, true); }

static int find_next_general_fasta_coi
   (seq*	_seq,
	int		allowComments)
	{
	char	buffer[maxSequenceHeader+1];
	char*	header;
	int		headerLen;
	int		mustBeHeader;
	int		leadingWhite;
	char	ch, *s;
	int		ix;

	debugNamesFile_12;

	mustBeHeader = true;

	while (true)
		{
		// find the next header

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;

		if ((allowComments) && (ch == '#'))
			{ // comment, skip to end-of-line and go back and try again
			while (ch != '\n')
				{
				ch = seq_getc (_seq);
				if (ch == EOF) goto failure;
				}
			continue;
			}

		if (ch != '>')
			{
			if (mustBeHeader)
				suicidef ("internal error in find_next_fasta_coi\n"
				          "processing %s, looking for \"%s\"\n",
				          sequence_filename(_seq), _seq->nextContigName);
			continue;
			}

		if (!mustBeHeader) _seq->contig++;
		mustBeHeader = false;

		// skip leading white space

		debugNamesFile_14;

		leadingWhite = 0;

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;
		while ((ch != '\n') && (isspace (ch)))
			{
			leadingWhite++;
			ch = seq_getc (_seq);
			if (ch == EOF) goto failure;
			}

		if (ch == '\n')
			continue;	// (unnamed sequence)

		// read the header

		s = buffer;
		while (ch != '\n')
			{
			if (s - buffer >= maxSequenceHeader)	// (overflow;
				break;								//  .. truncate the header)
			*(s++) = ch;
			ch = seq_getc (_seq);
			if (ch == EOF) goto failure;
			}
		*s = 0;

		// if we have a name trigger, locate the sequence's name

		if (_seq->nameTrigger != NULL)
			{
			header = strstr (buffer, _seq->nameTrigger);
			if (header == NULL) continue;	// (effectively unnamed sequence)
			header += strlen (_seq->nameTrigger);

			s = header;
			while ((*s != 0) && ((isalnum(*s)) || (*s == '_')))
				s++;
			headerLen = s-header;
			}
		else if (!_seq->useFullNames)
			{
			shorten_header (/* from */ buffer, _seq->nameParseType, false,
			                /* to   */ NULL, NULL);
			header    = buffer;
			headerLen = strlen(buffer);
			}
		else
			{
			header    = buffer;
			headerLen = strlen(buffer);
			}

		if ((_seq->nameParseType & name_parse_fill_white) != 0)
			whitespace_to_under (header, headerLen);

		// compare header to the contig-of-interest

		debugNamesFile_15;

		if (strncmp (header, _seq->nextContigName, headerLen) != 0)
			continue;
		if ((int) strlen (_seq->nextContigName) != headerLen)
			continue;

		break; // found a match!
		}

	debugNamesFile_16;

	// unget the header

	seq_ungetc (ch, _seq);	// (ch terminated the header)

	for (ix=strlen(buffer) ; ix>0 ; )
		seq_ungetc (buffer[--ix], _seq);

	while (leadingWhite-- > 0) seq_ungetc (' ', _seq);
	seq_ungetc ('>', _seq);

	_seq->contigPending = true;
	return true;

	// failure, the contig name was not found

failure:
	suicidef ("%s does not contain (or contains out of order)\n"
	          "         the sequence \"%s\"",
			  sequence_filename(_seq), _seq->nextContigName);
	return false; // (will never reach here)
	}

// find_next_fastq_coi--
//	advance to the next contig-of-interest in fastq file;  note that we don't
//	completely validate the file format here-- if we can locate a suitable
//	sequence header we use it;  validation is left to the sequence parser
//  (always returns true)

static int find_next_fastq_coi
   (seq*	_seq)
	{
	char	buffer[maxSequenceHeader+1];
	char*	header;
	int		headerLen;
	char	ch, *s;
	int		ix;
	int		ok;

	debugNamesFile_13;

	while (true)
		{
		// find the next header

		debugNamesFile_14;

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;

		if (ch != '@')
			suicidef ("internal error in find_next_fastq_coi\n"
			          "processing %s, looking for \"%s\"\n",
			          sequence_filename(_seq), _seq->nextContigName);

		// read the header

		ch = seq_getc (_seq);
		if (ch == EOF) goto failure;

		s = buffer;
		while ((ch != '\n') && (ch != '\r'))
			{
			if (s - buffer >= maxSequenceHeader)	// (overflow;
				break;								//  .. truncate the header)
			*(s++) = ch;
			ch = seq_getc (_seq);
			if (ch == EOF) goto failure;
			}
		*s = 0;

		if (ch == '\r') // handle possible DOS CR-LF line ending
			{
			ch = seq_getc (_seq);
			if (ch != '\n') seq_ungetc (ch, _seq);
			}

		// if we have a name trigger, locate the sequence's name

		if (_seq->nameTrigger != NULL)
			{
			header = strstr (buffer, _seq->nameTrigger);
			if (header == NULL) goto skip_content;	// (effectively unnamed sequence)
			header += strlen (_seq->nameTrigger);

			s = header;
			while ((*s != 0) && ((isalnum(*s)) || (*s == '_')))
				s++;
			headerLen = s-header;
			}
		else if (!_seq->useFullNames)
			{
			shorten_header (/* from */ buffer, _seq->nameParseType, false,
			                /* to   */ NULL, NULL);
			header    = buffer;
			headerLen = strlen(buffer);
			}
		else
			{
			header    = buffer;
			headerLen = strlen(buffer);
			}

		if ((_seq->nameParseType & name_parse_fill_white) != 0)
			whitespace_to_under (header, headerLen);

		// compare header to the contig-of-interest;  if this is not the
		// contig-of-interest, skip the sequence content

		debugNamesFile_15;

		if ((strncmp (header, _seq->nextContigName, headerLen) != 0)
		 || ((int) strlen (_seq->nextContigName) != headerLen))
			{
		skip_content:
			ok = fastq_skip_content (_seq);
			if (!ok) goto failure;
			continue;
			}

		break; // found a match!
		}

	debugNamesFile_16;

	// unget the header

	seq_ungetc (ch, _seq);	// (ch terminated the header)

	for (ix=strlen(buffer) ; ix>0 ; )
		seq_ungetc (buffer[--ix], _seq);

	seq_ungetc ('@', _seq);

	_seq->contigPending = true;
	return true;

	// failure, the contig name was not found

failure:
	suicidef ("%s does not contain (or contains out of order)\n"
	          "         the sequence \"%s\"",
			  sequence_filename(_seq), _seq->nextContigName);
	return false; // (will never reach here)
	}

// find_next_2bit_coi--
//	advance to the next contig-of-interest in 2bit header
//  (always returns true)

static int find_next_2bit_coi
   (seq*		_seq)
	{
	char		seqName[maxSequenceName+1];
	long int	savedContigFilePos;
	int			err;

	debugNamesFile_17;

	// position to the sequence's next index table entry

	err = fseek (_seq->f, _seq->twoBit.contigFilePos, SEEK_SET);
	if (err != 0)
		suicidef ("in find_next_2bit_coi(%s), index fseek(%08lX) returned %d",
		          sequence_filename(_seq), _seq->twoBit.contigFilePos, err);

	// read index table entries until we find the one we're looking for

	while (true)
		{
		if (_seq->contig >= _seq->twoBit.numContigs)
			suicidef ("%s does not contain (or contains out of order)\n"
			          "         the sequence \"%s\"",
			          sequence_filename(_seq), _seq->nextContigName);

		// read the sequence's next index table entry

		savedContigFilePos = ftell (_seq->f);
		/*seqDataPos=*/ read_2bit_index_entry (_seq, seqName, _seq->contig);
		if (strcmp (seqName, _seq->nextContigName) == 0) break;

		_seq->contig++;
		}

	_seq->twoBit.contigFilePos = savedContigFilePos;
	_seq->contigPending        = true;
	return true;
	}


// find_next_hsx_coi--
//	advance to the next contig-of-interest in hsx header
//  (always returns true)

static int find_next_hsx_coi
   (seq*	_seq)
	{
	u64		fileOffset;

	debugNamesFile_18;

	fileOffset = lookup_hsx_sequence (_seq, _seq->nextContigName);

	if ((fileOffset & hsxMsBit5) != 0)
		suicidef ("hsx file %s doesn't contain %s",
				  sequence_filename(_seq), _seq->nextContigName);
	if (fileOffset > hsxMaxFilePos)
		suicidef ("in find_next_hsx_coi for %s,"
				  " file pos for %s (%010lX) exceeds max (%010lX)",
				  sequence_filename(_seq), _seq->nextContigName, fileOffset);

	_seq->hsx.contigFilePos = fileOffset;
	_seq->contigPending     = true;
	debugNamesFile_19;
	return true;
	}

//----------
//
// read_contig_name--
//	Read the next name from a contigs-of-interest file.
//
// The file will contain one name per line.  Any leading whitespace is ignored,
// any comment lines are ignored (# is the comment character), and the name is
// only up to the first whitespace character.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	true if we were successful;  false if there are no more names in the file.
//
//----------

static int read_contig_name
   (seq*	_seq)
	{
	char*	line     = _seq->nextContigName;
	int		lineSize = sizeof(_seq->nextContigName);
	char	discard[maxSequenceName+1];
	int		len;
	int		missingEol;
	char*	waffle, *s;

	while (fgets (line, lineSize, _seq->namesFile) != NULL)
		{
		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line (and discard the extra)

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), _seq->namesFile) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines
		// nota bene: since illumina read names contain our comment character
		// (#), and to maintain backward compatibility for lines that contain
		// a contig name *and* a comment, the comment character now requires a
		// space or tab just before it (unless it is at the start of the line)

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		while (waffle != NULL)
			{
			if (waffle == line)
				{ *waffle = 0;  break; }
			else if ((waffle[-1] == ' ') || (waffle[-1] == '\t'))
				{ *waffle = 0;  break; }
			waffle = strchr (waffle+1, '#');
			}

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it

		s = skip_darkspace(line);
		*s = 0;

		debugNamesFile_20;

		return true;
		}

	return false;
	}

//----------
//
// read_chore--
//	Read the next chore from a chores file.
//
// The file contains one chore per line, and any comment lines are ignored. "#"
// is the comment character, but only if it appears at the beginning of a line
// or immediately after whitespace.
//
// A chore contains the following fields (some of which are optional):
//
//	<tName> <tStart> <tEnd> <qName> [<qStart> <qEnd>] [<qStrand>] [id=<tag>]
//
// In cases where the target name is irrelevant (i.e. there is only one name in
// the target sequence file), "*" can replace <tName>.  Similarly, if we don't
// have a target (or query) subrange, "* *" can be used.  Note that the query
// subrange and strand are optional, as is the tag.
//
// The tag can be any short string (but without whitespace) the user wants to
// associate with the chore (maximum length is defined by maxChoreTagLen).
// This tag can be reported along with alignments for the chore, in the general
// tab-delimited format.
//
// Note that intervals in the file are origin-one half open, and they are
// *not* altered as written into the chore struct.
//
// Typical lines:
//
//	chr11 5931512 5931843  APPLE_READ_00009
//	chr11 5931512 5931843  APPLE_READ_00009 +
//	chr11 5931512 5931843  APPLE_READ_00009 -
//	*     5931512 5931843  APPLE_READ_00036
//	chr22 *       *        APPLE_READ_00087
//	chr11 2878300 1933292  chr11 1486276 1741268 +
//	chr11 2878300 1933292  chr2  6865671 7149925 -
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	true if we were successful;  false if there are no more chores in the file.
//
//----------

static int read_chore
   (seq*	_seq)
	{
	char	line[511+1], discard[511+1];
	int		len;
	int		missingEol;
	char*	scan, *waffle, *field;
	int		numItems, charsUsed;
	char*	tName, *tStart, *tEnd;
	char*	qName, *qStart, *qEnd, *qStrand;
	char*	idTag;
	char*	header;

	debugNamesFile_21;

	while (fgets (line, sizeof(line), _seq->choresFile) != NULL)
		{
		_seq->choresLineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line (and discard the extra)

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), _seq->choresFile) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines
		// nota bene: since illumina read names contain our comment character
		// (#), and to allow lines that contain a chore *and* a comment, the
		// comment character now requires a space or tab just before it (unless
		// it is at the start of the line)

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		while (waffle != NULL)
			{
			if (waffle == line)
				{ *waffle = 0;  break; }
			else if ((waffle[-1] == ' ') || (waffle[-1] == '\t'))
				{ *waffle = 0;  break; }
			waffle = strchr (waffle+1, '#');
			}

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it;  parse to find the chore fields

		debugNamesFile_22;

		tName = scan = line;

		scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
		if (*scan == 0) goto not_enough_fields;
		tStart = scan;

		scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
		if (*scan == 0) goto not_enough_fields;
		tEnd = scan;

		scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
		if (*scan == 0) goto not_enough_fields;
		qName = scan;

		qStart = qEnd = qStrand = idTag = NULL;
		scan = skip_darkspace (scan);  scan = skip_whitespace (scan);

		if ((*scan != 0)
		 && ((scan[0] != '+') || ((scan[1] != 0) && (!isspace (scan[1]))))
		 && ((scan[0] != '-') || ((scan[1] != 0) && (!isspace (scan[1]))))
		 && (strncmp (scan, "id=", strlen("id+")) != 0))
			{
			qStart = scan;

			scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
			if (*scan == 0) goto missing_query_end;
			qEnd = scan;

			scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
			}

		if (((scan[0] == '+') && ((scan[1] == 0) || (isspace (scan[1]))))
		 || ((scan[0] == '-') && ((scan[1] == 0) || (isspace (scan[1])))))
			{
			qStrand = scan;
			scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
			}

		if (strncmp (scan, "id=", strlen("id+")) == 0)
			{
			idTag = scan + strlen("id+");
			scan = skip_darkspace (scan);  scan = skip_whitespace (scan);
			}

		if (*scan != 0) goto extra_fields;

		// interpret the chore fields

		_seq->chore.tSubrange = false;
		_seq->chore.qSubrange = false;

		field = tStart;
		scan = skip_darkspace (field);  *scan = 0;
		if (strcmp (field, "*") != 0)
			{
			charsUsed = -1;
			numItems = sscanf (field, unsposFmtScanf "%n", &_seq->chore.tStart, &charsUsed);
			if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
			if (_seq->chore.tStart == 0) goto bad_target_start1;
			_seq->chore.tSubrange = true;
			}

		field = tEnd;
		scan = skip_darkspace (field);  *scan = 0;
		if (strcmp (field, "*") == 0)
			{ if (_seq->chore.tSubrange) goto bad_target_end; }
		else
			{
			charsUsed = -1;
			numItems = sscanf (field, unsposFmtScanf "%n", &_seq->chore.tEnd, &charsUsed);
			if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
			if (!_seq->chore.tSubrange) goto bad_target_start2;
			if (_seq->chore.tEnd <= _seq->chore.tStart) goto bad_target_interval;
			}

		if (qStart != NULL)
			{
			field = qStart;
			scan = skip_darkspace (field);  *scan = 0;
			if (strcmp (field, "*") != 0)
				{
				charsUsed = -1;
				numItems = sscanf (field, unsposFmtScanf "%n", &_seq->chore.qStart, &charsUsed);
				if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
				if (_seq->chore.qStart == 0) goto bad_query_start1;
				_seq->chore.qSubrange = true;
				}

			field = qEnd;
			scan = skip_darkspace (field);  *scan = 0;
			if (strcmp (field, "*") == 0)
				{ if (_seq->chore.qSubrange) goto bad_query_end; }
			else
				{
				charsUsed = -1;
				numItems = sscanf (field, unsposFmtScanf "%n", &_seq->chore.qEnd, &charsUsed);
				if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
				if (!_seq->chore.qSubrange) goto bad_query_start2;
				if (qEnd <= qStart) goto bad_query_interval;
				}
			}

		if (qStrand == NULL)
			_seq->chore.qStrand = 1;										// (both strands)
		else
			{
			scan = skip_darkspace (qStrand);  *scan = 0;
			if (strcmp (qStrand, "+") == 0) _seq->chore.qStrand = 0;	// (forward strand only)
			                           else _seq->chore.qStrand = -1;	// (reverse strand only)
			}

		scan = skip_darkspace (tName);  *scan = 0;
		if (strcmp (tName, "*") == 0)
			_seq->chore.tName[0] = 0;
		else if (strlen(tName) >= sizeof(_seq->chore.tName))
			goto target_name_too_long;
		else
			strcpy (/*to*/ _seq->chore.tName, /*from*/ tName);

		scan = skip_darkspace (qName);  *scan = 0;
		if (strlen(qName) >= sizeof(_seq->nextContigName))
			goto query_name_too_long;
		else
			strcpy (/*to*/ _seq->nextContigName, /*from*/ qName);

		header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
		if ((header == NULL) || (strcmp (header, _seq->nextContigName) != 0))
			_seq->chore.num = 1;
		else
			_seq->chore.num++;

		if (idTag == NULL)
			_seq->chore.idTag[0] = 0;
		else
			{
			scan = skip_darkspace (idTag);  *scan = 0;
			if (strlen(idTag) >= sizeof(_seq->chore.idTag))
				goto id_tag_too_long;
			strcpy (/*to*/ _seq->chore.idTag, /*from*/ idTag);
			}

		debugNamesFile_23;

		return true;
		}

	return false;

// failure exits

not_enough_fields:
	suicidef ("bad chore (in %s, line %d): \"%s\"\n"
			  "not enough fields in line",
			  _seq->choresFilename, _seq->choresLineNum, line);
	return false; // (never gets here)

extra_fields:
	suicidef ("bad chore (in %s, line %d): \"%s\"\n"
			  "extra fields in line: \"%s\"",
			  _seq->choresFilename, _seq->choresLineNum, line, scan);
	return false; // (never gets here)

missing_query_end:
	suicidef ("bad chore (in %s, line %d): \"%s\"\n"
			  "has start of query subrange but not end",
			  _seq->choresFilename, _seq->choresLineNum, line);
	return false; // (never gets here)

bad_field:
	suicidef ("bad chore field (in %s, line %d): \"%s\"",
			  _seq->choresFilename, _seq->choresLineNum, field);
	return false; // (never gets here)

bad_target_start1:
	suicidef ("bad chore target interval (in %s, line %d)\n"
	          "start cannot be zero",
			  _seq->choresFilename, _seq->choresLineNum);
	return false; // (never gets here)

bad_target_start2:
	suicidef ("bad chore target interval (in %s, line %d): * " unsposFmt "\n"
	          "can't wildcard start and not end",
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.tEnd);
	return false; // (never gets here)

bad_target_end:
	suicidef ("bad chore target interval (in %s, line %d): " unsposFmt " *\n"
	          "can't wildcard end and not start",
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.tStart);
	return false; // (never gets here)

bad_target_interval:
	suicidef ("bad chore target interval (in %s, line %d): " unsposFmt ">=" unsposFmt,
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.tStart, _seq->chore.tEnd);
	return false; // (never gets here)

bad_query_start1:
	suicidef ("bad chore query interval (in %s, line %d)\n"
	          "start cannot be zero",
			  _seq->choresFilename, _seq->choresLineNum);
	return false; // (never gets here)

bad_query_start2:
	suicidef ("bad chore query interval (in %s, line %d): * " unsposFmt "\n"
	          "can't wildcard start and not end",
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.qEnd);
	return false; // (never gets here)

bad_query_end:
	suicidef ("bad chore query interval (in %s, line %d): " unsposFmt " *\n"
	          "can't wildcard end and not start",
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.qStart);
	return false; // (never gets here)

bad_query_interval:
	suicidef ("bad chore query interval (in %s, line %d): " unsposFmt ">=" unsposFmt,
			  _seq->choresFilename, _seq->choresLineNum, _seq->chore.qStart, _seq->chore.qEnd);
	return false; // (never gets here)

target_name_too_long:
	suicidef ("chore target name too long (in %s, line %d): \"%s\"",
			  _seq->choresFilename, _seq->choresLineNum, tName);
	return false; // (never gets here)

query_name_too_long:
	suicidef ("chore query name too long (in %s, line %d): \"%s\"",
			  _seq->choresFilename, _seq->choresLineNum, qName);
	return false; // (never gets here)

id_tag_too_long:
	suicidef ("chore id tag too long, allowed length is %d (in %s, line %d): \"%s\"",
			   sizeof(_seq->chore.idTag)-1, _seq->choresFilename, _seq->choresLineNum, idTag);
	return false; // (never gets here)
	}

//----------
//
// create_short_header--
//	Convert a sequence's header into a shorter version.  The shorter version is
//	intended to be useful as a sequence's name in maf or axt files.
//
// Examples:
//
//	>~username/human/hg18/_seq/chr16.nib:120000-190000           chr16
//  owl_monkey 122000-180000 of owl_monkey.ENm008.fa             owl_monkey
//  > armadillo|ENm001|JAN-2006|9361|NISC|...|1|1|.              armadillo
//	>reverse complement of ~username/human/hg18/_seq/chr14.nib   chr14
//	>positions 180000-250000 of armadillo|ENm008|...             armadillo
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//
// Returns:
//	(nothing;  _seq->shortHeader and _seq->shortHeaderSize are modified)
//
//----------
//
// Note:	It is possible for the resulting short header name to be an empty
//			string.
//
//----------

static void create_short_header
   (seq*	_seq)
	{
	int		skipPath;

	if (_seq->header == NULL)
		{
		if ((_seq->shortHeader != NULL) && (_seq->shortHeaderSize != 0))
			_seq->shortHeader[0] = 0;
		return;
		}

	if ((!_seq->headerOwner) || (!_seq->shortHeaderOwner))
		{
		char* name = (_seq->filename != NULL)? _seq->filename
		                                     : _seq->header;
		suicidef ("internal error, attempt to shorten external sequence header (%s)",
		          name);
		}

	if (strstr(_seq->header,"{number}") != NULL)
		expand_nickname( /* from */ _seq->header, _seq->contig,
	                     /* to   */ &_seq->shortHeader, &_seq->shortHeaderSize);
	else
		{
		skipPath = (_seq->fileType == seq_type_nib);
		shorten_header (/* from */ _seq->header, _seq->nameParseType, skipPath,
		                /* to   */ &_seq->shortHeader, &_seq->shortHeaderSize);
		}
	}


static void shorten_header
   (char*	src,
	int		nameParseType,
	int		skipPath,
	char**	_dst, // (NULL => write it in place)
	u32*	_dstSize)
	{
	char*	dst;
	u32		dstSize;
	char*	h, *hh, *s;
	u32		len, sLen;

	// skip fasta '>', leading whitespace, and/or a path

	h = src;
	if (h[0] == '>') h++;
	h = skip_whitespace (h);

	// skip "reverse complement of" and/or "positions A-B of"

	s = "reverse complement of ";
	if (strcmp_prefix (h, s) == 0)
		h = skip_whitespace (h + strlen (s));

	s = "positions ";
	if (strcmp_prefix (h, s) == 0)
		{
		hh = skip_whitespace (h + strlen (s));
		hh = skip_darkspace (hh);
		hh = skip_whitespace (hh);
		s = "of ";
		if (strcmp_prefix (hh, s) == 0)				// (we only change h if
			h = skip_whitespace (hh + strlen (s));	//  .. "of" is present)
		}

	// skip a path

	if (skipPath)
		{
		while (true)
			{
			hh = strchr (h, pathSlash);
			if (hh == NULL) break;
			h = hh + 1;
			}
		}

	h = skip_whitespace (h);

	// figure out the length to copy;  we'll truncate at the first space or
	// "funny" character

	if (parse_type(nameParseType) == name_parse_type_alnum)
		{
		len = strspn (h, "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		                 "abcdefghijklmnopqrstuvwxyz"
		                 "0123456789" "_");
		goto skip_suffix_trim;
		}
	else if (parse_type(nameParseType) == name_parse_type_darkspace)
		len = strcspn (h, " \t");
	else // if (parse_type(nameParseType) == name_parse_type_core)
		len = strcspn (h, " \t|:");

	// if the suffix is ".nib", ".fasta", etc., remove it

	s = ".nib";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".2bit";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".hsx";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".fasta";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

	s = ".fa";
	sLen = strlen(s);
	if ((len > sLen) && (strncmp (h+len-sLen,s,sLen) == 0))
		len -= sLen;

skip_suffix_trim:

	// create the header

	if (_dst == NULL)
		{
		strncpy (src, h, len);
		src[len] = 0;

		if ((nameParseType & name_parse_fill_white) != 0)
			whitespace_to_under (src, strlen(src));
		}
	else
		{
		dst     = *_dst;
		dstSize = *_dstSize;

		if (len+1 > dstSize)
			{
			dst = realloc_or_die ("shorten_header", dst, len+1);
			dstSize = len+1;
			}

		strncpy (dst, h, len);
		dst[len] = 0;
		*_dst     = dst;
		*_dstSize = dstSize;

		if ((nameParseType & name_parse_fill_white) != 0)
			whitespace_to_under (dst, strlen(dst));
		}

	}


static void whitespace_to_under (char* s, int sLen)
	{ for ( ; sLen-->0 ; s++) { if (isspace (*s)) *s = '_'; } }


static void expand_nickname
   (char*	src,
	u32		contigNumber,
	char**	_dst,
	u32*	_dstSize)
	{
	char*	dst;
	u32		dstSize;
	char*	s, *d, *expand;
	u32		len;

	// determine the size of the resulting header

	len = strlen (src)
	    - strlen ("{number}")
	    + snprintf (NULL, 0, unsposFmt, contigNumber);

	// allocate the header (if necessary)

	dst     = *_dst;
	dstSize = *_dstSize;

	if (len+1 > dstSize)
		{
		dst = realloc_or_die ("expand_nickname", dst, len+1);
		dstSize = len+1;
		}

	// create the header

	s = src;
	d = dst;

	expand = strstr (src, "{number}");
	if (expand > src)
		{
		strncpy (d, s, expand-src);
		d += expand-src;
		s =  expand + strlen("{number}");
		}

	sprintf (d, unsposFmt, contigNumber);
	d += strlen(d);

	strcpy (d, s);

	*_dst     = dst;
	*_dstSize = dstSize;
	}

//----------
//
// separate_sequence--
//	Separate each of a sequence's subsequences, spliting them into pieces
//	wherever a specified character occurs.
//
// We require that the sequence is already partitioned (though it may have only
// one subsequence), and we introduce additional partition blocks whenever a
// subsequence is split into pieces.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence to separate.
//	char	sepCh:	The character at which to separate the sequence.
//
// Returns:
//	nothing
//
//----------
//
//	Notes:
//
//	(1)	The term 'separator' has two similar meanings in this routine.  The
//		NUL characters between partitions are called separators, as are the
//		user-defined separator characters in the incoming sequence.  One of the
//		results of this routine is to turn those latter characters into the
//		former.
//
//----------

static void separate_sequence
   (seq*			_seq,
	char			_sepCh)
	{
	seqpartition*	sp = &_seq->partition;
	u8				sepCh = (u8) _sepCh;
	u32				extraPieces;
	int				inSepRun;
	u8				ch, chBefore, chAfter;
	u32				newSpLen;
	partition*		p, *pFrom, *pTo;
	u32				fromIx, toIx;
	unspos			sepPrefix, sepSuffix, sepBefore, sepAfter;
	unspos			startLoc;
	u8*				scan;

	if (sp->p == NULL)
		suicidef ("internal error in separate_sequence\n"
		          "sequence has no partition table");
	if (sp->state == seqpart_empty)
		suicidef ("internal error in separate_sequence\n"
		          "partition table is in empty state");
	if (sp->state == seqpart_reusable)
		suicidef ("internal error in separate_sequence\n"
		          "partition table is in re-usable state");
	if (sp->state == seqpart_loading)
		suicidef ("internal error in separate_sequence\n"
		          "partition table is in loading state");

	debugSeparation_1;

	// scan the sequence, looking for runs of separators, to see how many
	// additional partitions we'll need;  in rough terms, we will need an
	// extra partition whenever we see a transition from nuc to sep or sep to
	// nuc;  however, if the sep has a NUL on the other end (marking the start
	// or end of the current partition), then we will not need an extra
	// partition
	//
	// examples (0 is NUL, X is a separator, and n is a nucleotide):
	//    0XXXXXnnnnnnnnnnXXXXX0                   no extras needed
	//    0XXXXXnnnnnnnnnnXXXXXXnnnnnnnnnnXXXXX0   one extra needed

	extraPieces = 0;
	inSepRun = false;
	chBefore = 0;
	for (scan=_seq->v ; scan<_seq->v+_seq->len ; scan++)
		{
		ch = *scan;
		if (ch == 0)
			inSepRun = false;
		else if (ch == sepCh)
			{ if ((chBefore != 0) && (chBefore != sepCh)) inSepRun = true; }
		else
			{
			if (inSepRun)
				{
				extraPieces++;
				debugSeparation_2;
				}
			inSepRun = false;
			}
		chBefore = ch;
		}

	debugSeparation_3;

	// allocate the extra partitions;  note that we won't need any additional
	// space in the header pool, since the new partitions will share existing
	// names

	newSpLen = sp->len + extraPieces;
	if ((newSpLen < sp->len) || (newSpLen == u32max))
		suicidef ("in separate_sequence, "
		          "number of partitions overflows internal data type");

	if (extraPieces != 0)
		enough_partitions (_seq, newSpLen, sp->poolLen,
		                   /*anticipate*/ false, /*round up*/ true);

	// scan the current partitions from last to first, dowing the following:
	//  (1) expanding any partition that contains multiple pieces
	//  (2) adjusting partition bounds whenever a separator exists at an end
	//  (3) replacing all separators with NUL characters

	p = sp->p;
	fromIx = sp->len;  pFrom = &p[fromIx];
	toIx   = newSpLen;   pTo   = &p[toIx];
	pTo->sepBefore = pFrom->sepBefore;		// copy sentinel
	debugSeparation_4;

	sepPrefix = pFrom->sepBefore;
	while (fromIx-- > 0)
		{
		pFrom--;
		sepSuffix = sepPrefix;			// (separator at end of this partition)
		sepPrefix = pFrom->sepBefore;	// (separator at start of this partition)

		debugSeparation_5;

		sepAfter = 0;
		chAfter  = 0;
		for (scan=_seq->v+sepSuffix-1 ; scan>_seq->v+sepPrefix ; scan--)
			{
			ch = *scan;
			//debugSeparation_6;
			if (ch == 0)
				{
				suicidef ("internal error in separate_sequence\n"
				          "seq->v[" unsposFmt "]=0x00", (unspos) (scan-_seq->v));
				}
			else if (ch == sepCh)
				{
				*scan = 0; // replace separator with NUL
				if ((chAfter != 0) && (chAfter != sepCh) && (sepAfter != 0))
					{
					toIx--;  pTo--;

					sepBefore = scan - _seq->v;
					startLoc  = pFrom->startLoc + (sepBefore - pFrom->sepBefore);

					pTo->sepBefore = sepBefore;
					pTo->sepAfter  = sepAfter;
					pTo->contig    = pFrom->contig;
					pTo->startLoc  = startLoc;
					pTo->trueLen   = pFrom->trueLen;
					pTo->header    = pFrom->header;
					debugSeparation_8;

					sepAfter = 0;
					}
				}
			else if (sepAfter == 0)
				{
				sepAfter = (scan - _seq->v) + 1;
				debugSeparation_7;
				}
			chAfter = ch;
			}

		if (sepAfter != 0)
			{
			toIx--;  pTo--;

			sepBefore = scan - _seq->v;
			startLoc  = pFrom->startLoc + (sepBefore - pFrom->sepBefore);

			pTo->sepBefore = sepBefore;
			pTo->sepAfter  = sepAfter;
			pTo->contig    = pFrom->contig;
			pTo->startLoc  = startLoc;
			pTo->trueLen   = pFrom->trueLen;
			pTo->header    = pFrom->header;
			debugSeparation_8;
			}
		}
	debugSeparation_9;

	if (toIx != 0)
		suicidef ("internal error in separate_sequence\n"
		          "toIx=%d", toIx);

	sp->len = newSpLen;
	debugSeparation_10;
	}

//----------
//
// add_partition--
//	Add a partition block for a sequence.
//
//----------
//
// Arguments:
//	seq*	_seq:		The sequence to add the partition to.
//	unspos	sepPos:		The position of the partition's prefix separator;  see
//						.. note (1) below.
//	unspos	startLoc:	The partition's startLoc field.
//	unspos	trueLen:	The partition's trueLen field.
//
// Returns:
//	nothing
//
//----------
//
//	Notes:
//
//	(1)	We use sepPos as this partition's prefix separator (sepBefore) *and* as
//		the previous partition's suffix separator (sepAfter). 
//
//	(2)	We initially write a zero to this partition's sepAfter.  It is replaced
//		with the correct on the next call (the call to add the next partition).
//		However, this leaves the final partition's sepAfter field unset.  So the
//		caller *must* set it.
//
//----------

static void add_partition
   (seq*			_seq,
	unspos			sepPos,
	unspos			startLoc,
	unspos			trueLen)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	char*			header;
	int				headerLen;

	debugTextFile_2;

	header = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
	headerLen = strlen(header);

	enough_partitions (_seq, sp->len+1, sp->poolLen+headerLen+1,
					   /*anticipate*/ true, /*round up*/ true);

	if (sp->len > 0)
		{
		p = &sp->p[sp->len-1];
		p->sepAfter = sepPos;
		}

	p = &sp->p[sp->len];
	p->sepBefore = sepPos;
	p->sepAfter  = 0;  // (see note 2)
	p->contig    = _seq->contig; 
	p->startLoc  = startLoc;
	p->trueLen   = trueLen;
	p->header    = sp->poolLen;
	strcpy (/*to*/ &sp->pool[p->header], /*from*/ header);
	sp->poolLen += headerLen+1;
	sp->len++;
	}

//----------
//
// copy_partitions--
//	Make a copy of a sequence's parititions.
//
//----------
//
// Arguments:
//	seq*	seqTo:		The sequence to copy partition info to.
//	seq*	seqFrom:	The sequence to copy partition info from.
//
// Returns:
//	nothing;  failures result in fatality.
//
//----------

static void copy_partitions
   (seq*			seqTo,
	seq*			seqFrom)
	{
	seqpartition*	spFrom = &seqFrom->partition;
	seqpartition*	spTo   = &seqTo->partition;
	u32				len, poolLen;
	u32				ix;

	// make sure we have space for the partitions

	poolLen = spFrom->poolLen;
	spTo->poolOwner = true;
	enough_partitions (seqTo, spFrom->len, poolLen,
	                   /*anticipate*/ false, /*round up*/ false);

	// copy the pool directly

	spTo->poolLen = poolLen;
	memcpy (spTo->pool, spFrom->pool, poolLen);

	// copy the partition records

	len = spTo->len = spFrom->len;
	for (ix=0 ; ix<=len ; ix++)
		spTo->p[ix] = spFrom->p[ix];
	}

//----------
//
// enough_partitions--
//	Make sure a sequence has enough room for a specified number of partitions.
//
//----------
//
// Arguments:
//	seq*	_seq:			The sequence to check.
//	u32		numPartitions:	The number of partitions to allocate for (not
//							.. including the extra partition used as a
//							.. sentinel).
//	u32		poolSize:		The number of bytes to allocate for a pool of
//							.. headers.  If this is zero, we will estimate the
//							.. pool size from the number of partitions.
//	int		anticipate:		true  => allocate extra, anticipating the need for
//							         .. more
//							false => don't
//	int		roundUp:		true  => round up the allocation size to some
//							         .. convenient size
//							false => don't
//
// Returns:
//	nothing;  _seq->partition.p and _seq->partition.pool may be modified;
//	failures result in fatality.
//
//----------

#define averageHeaderSize 20

static void enough_partitions
   (seq*			_seq,
	u32				numPartitions,
	u32				poolSize,
	int				anticipate,
	int				roundUp)
	{
	seqpartition*	sp = &_seq->partition;
	u32				bytesNeeded;

	// if we have enough already, just return

	if ((sp->p != NULL)
	 && (sp->size > numPartitions)
	 && (sp->pool != NULL)
	 && ((poolSize > 0) && (sp->poolSize >= poolSize)))
		return;

	if (!sp->poolOwner)
		{
		char* name = (_seq->filename != NULL)? _seq->filename
		                                     : _seq->header;
		suicidef ("internal error, attempt to resize external partition names pool (%s)",
		          name);
		}

	if (sp->p    == NULL) sp->len     = 0;
	if (sp->pool == NULL) sp->poolLen = 0;

	if (poolSize == 0) poolSize = numPartitions * (averageHeaderSize + 1);

	// allocate partition array;  note that we bump up the number of records
	// allocated to as many as can fit in a multiple of 16K

	if (sp->size <= numPartitions)
		{
		numPartitions++;					// (extra one for a sentinel)
		if (anticipate)						// anticipatory, grow by about 13%
			numPartitions += 30 + numPartitions / 8;
		bytesNeeded = numPartitions * sizeof(partition);

		if (roundUp)
			{
			bytesNeeded   = round_up_16K (bytesNeeded);
			numPartitions = bytesNeeded / sizeof(partition);
			bytesNeeded   = numPartitions * sizeof(partition);
			}

		sp->p    = realloc_or_die ("enough_partitions (p)", sp->p, bytesNeeded);
		sp->size = numPartitions;
		}

	// allocate pool for partition headers

	if (sp->poolSize < poolSize)
		{
		if (anticipate)						// anticipatory, grow by about 13%
			poolSize += 30*(averageHeaderSize+1) + poolSize / 8;
		bytesNeeded   = round_up_16K (poolSize);
		sp->pool      = realloc_or_die ("enough_partitions (pool)",
		                                sp->pool, bytesNeeded);
		sp->poolSize  = bytesNeeded;
		sp->poolOwner = true;
		}

	}

//----------
//
// lookup_partition, lookup_partition_no_die--
//	Map a position in a partitioned sequence to its partition record.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//	unspos	pos:	The position to look up.  This is an index into _seq->v,
//					.. (origin zero).
//
// Returns:
//	A pointer to the partition record (see note (1)); lookup_partition_no_die
//	returns NULL on failure;  for lookup_partition, failures result in fatality.
//
//----------
//
//	Notes:
//
//	(1)	The pointer p returned points to the partition record which contains the
//		contig, header, sepBefore and sepAfter for the partition.  sepBefore is
//		the bounding NUL at the left end (lower index) of the sequence, and
//		sepAfter is the bounding NUL at the right end (higher index) of the
//		sequence.
//
//	(2)	It has been suggested that there are cases for which caching the latest
//		lookup is beneficial.  Code to do this can be enabled by #defining
//		cache_partition_lookups.  However, in the author's tests the cached
//		lookup was rarely hit.
//
//----------

static partition* lookup_partition_core (seq* _seq, unspos pos, int dieOnFailure);

partition* lookup_partition_no_die (seq* _seq, unspos pos)
	{ return lookup_partition_core (_seq, pos, /* dieOnFailure */ false); }

partition* lookup_partition (seq* _seq, unspos pos)
	{ return lookup_partition_core (_seq, pos, /* dieOnFailure */ true); }


partition* lookup_partition_core
   (seq*				_seq,
	unspos				pos,
	int					dieOnFailure)
	{
#ifdef cache_partition_lookups // see note (2)
	static seq*			cachedSeq    = NULL;
	static unspos		cachedPos    = ((unspos) -1);
	static partition*	cachedResult = NULL;
#endif // cache_partition_lookups
	seqpartition*		sp = &_seq->partition;
	partition*			p;
	u32					hi, lo, ix;
	char*				reason;

	if (sp->p   == NULL) goto no_partitions;
	if (sp->len == 0)    goto no_partitions;

	sequence_count_stat (partitionLookups);

#ifdef cache_partition_lookups
	if ((_seq == cachedSeq) && (pos == cachedPos))
		{
		sequence_count_stat (partitionHits);
		return cachedResult;
		}
#endif // cache_partition_lookups

	p = sp->p;

	lo = 0;
	hi = sp->len;

	debugPartitions_2;

	if (pos <= p[lo].sepBefore) goto before_first;
	if (pos >= p[hi].sepBefore) goto after_last;

	// perform binary search for the position
	//
	// loop invariants:					loop termination:
	//		0 <= lo < hi <= len				ix = lo = hi-1
	//		pos > p[lo].sepBefore			pos > p[ix].sepBefore
	//		pos < p[hi].sepBefore			pos < p[ix+1].sepBefore

	while (true)
		{
		sequence_count_stat (lookupIterations);
		ix = (lo + hi) / 2;			// when hi==lo+1, ix==lo
		debugPartitions_3;
		if (hi == lo+1) break;
		if      (pos < p[ix].sepBefore) hi = ix;
		else if (pos > p[ix].sepBefore) lo = ix;
		else goto on_separator;		// pos == p[ix].sepBefore, which is illegal
		}

	// make sure the position was within the actual partition
	// nota bene:  we use ">" rather than ">=" to allow the caller to position
	//             on the open end of the sequence

	if (pos > p[ix].sepAfter) goto not_in_partition;

	// success

	debugPartitions_4;

#ifdef cache_partition_lookups
	cachedSeq    = _seq;
	cachedPos    = pos;
	cachedResult = &sp->p[ix];
#endif // cache_partition_lookups

	return &sp->p[ix];

	// failure

no_partitions:
	reason = "there are no partitions";
	goto failure;

before_first:
	reason = "before first partition";
	goto failure;

after_last:
	reason = "after last partition";
	goto failure;

on_separator:
	reason = "on partition separator prefix";
	goto failure;

not_in_partition:
	reason = "not within partition";
	goto failure;

failure:
	if (!dieOnFailure) return NULL;
	if (_seq->filename == NULL)
		suicidef ("lookup_partition could not locate position " unsposFmt "\n%s",
		          pos, reason);
	else
		suicidef ("lookup_partition could not locate position " unsposFmt " in %s\n%s",
		          pos, _seq->filename, reason);


	return NULL; // (never gets here)
	}

//----------
//
// lookup_named_partition--
//	Map a name to the corresponding partition record in a partitioned sequence.
//
//----------
//
// Arguments:
//	seq*	_seq:	The sequence.
//	name	name:	The name to look up.
//
// Returns:
//	A pointer to the partition record (see note 1);  NULL if not found.
//
//----------
//
//	Notes:
//
//	(1)	The parition record returned is the first with the given name.  If
//		there are other partitions with the same name, lookup_partition_seq_pos
//		should then be used to determine the partition containing that position.
//		Or, last_partition_with_name can then be used to locate the last
//		partition with that same given name.
//
//----------

partition* lookup_named_partition
   (seq*			_seq,
	char*			name)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		part;
	u32				ix;
	int				found;

	if (sp->p   == NULL) return NULL;
	if (sp->len == 0)    return NULL;

	// perform linear search for the name

	found = false;
	for (ix=0 ; ix<sp->len ; ix++)
		{
		part = &sp->p[ix];
		if (strcmp (&sp->pool[part->header], name) == 0)
			{ found = true;  break; }
		}
	if (!found) return NULL;

	return part;
	}

//----------
//
// last_partition_with_name--
//	Find the last partition record, in a partitioned sequence, that has the
//	same name as a given partition record.
//
//----------
//
// Arguments:
//	seq*		_seq:	The sequence.
//	partition*	part:	The first partition with the desired name.
//
// Returns:
//	A pointer to the partition record.
//
//----------

partition* last_partition_with_name
   (seq*			_seq,
	partition*		firstPart)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		scanPart, *prevPart;
	u32				ix;
	char*			name;

	if (sp->p   == NULL) return NULL;
	if (sp->len == 0)    return NULL;

	// determine the index into the partition list
	// $$$ perhaps to be safe we should use integer arithmetic to make sure
	//     .. firstPart - sp->p is a multiple of sizeof(partition)

	ix = firstPart - sp->p;
	if (ix >= sp->len)
		suicidef ("internal error in lookup_named_partition\n"
		          "invalid partition pointer");

	// get the name of this partition

	name = &sp->pool[firstPart->header];

	// scan forward along the partition list until we find a different name

	prevPart = firstPart;
	for (ix=ix+1 ; ix<sp->len ; ix++)
		{
		scanPart = &sp->p[ix];
		if (strcmp (&sp->pool[scanPart->header], name) != 0) break;
		prevPart = scanPart;
		}

	return prevPart;
	}

//----------
//
// lookup_partition_seq_pos--
//	Given the first partition for a given name, map a sequence position to the
//	corresponding partition record.
//
//----------
//
// Arguments:
//	seq*		_seq:	The sequence.
//	partition*	part:	The first partition for a given sequence (by name).
//	unspos		pos:	The position (within the sequence) to look for.  This
//						.. is origin-one.
//
// Returns:
//	A pointer to the partition record;  NULL if not found.
//
//----------

partition* lookup_partition_seq_pos
   (seq*			_seq,
	partition*		_part,
	unspos			pos)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		part = _part;
	char*			name;
	unspos			endLoc;
	u32				ix;
	int				found;

	if (sp->p   == NULL) return NULL;
	if (sp->len == 0)    return NULL;
	if (_part   == NULL) return NULL;

	name = &sp->pool[part->header];
	ix   = part - sp->p;

	found = false;
	while (true)
		{
		if (pos >= part->startLoc)
			{
			endLoc = part->startLoc + part->sepAfter - part->sepBefore;
			if (pos < endLoc)
				{ found = (pos < endLoc-1);  break; }
			}
		if (++ix >= sp->len) break;  // (not found)
		part = &sp->p[ix];
		if (strcmp (&sp->pool[part->header], name) != 0) break;  // (not found)
		}
	if (!found) return NULL;
	
	return part;
	}

//----------
//
// print_sequence--
//	Print a sequence to a file.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to write to.
//	seq*	seq:		The sequence to print.
//	char*	header:		The header to give the sequence.  If this is NULL,
//						.. _seq->header is used.  If this is empty, we assume
//						.. the caller has already taken care of printing the
//						.. header.
//	int		perLine:	The number of nucleotides to print per line.  Zero
//						.. indicates that they should all be printed on one
//						.. line.
//
// Returns:
//	(nothing)
//
//----------

void print_sequence
   (FILE*	f,
	seq*	_seq,
	char*	header,
	int		perLine)
	{
	unspos	ix;

	if (_seq == NULL)
		{ fprintf (f, "(null sequence)\n");  return; }

	if ((header == NULL) || (header[0] != 0))
		{
		if (header == NULL)
			{
			header = _seq->header;
			if (header == NULL)
				header = "";
			}

		if (header[0] == '>')
			header = skip_whitespace(header+1);

		if (header[0] == 0) fprintf (f, ">\n");
		               else fprintf (f, "> %s\n", header);
		}

	if (_seq->fileType == seq_type_qdna)
		{
		for (ix=0 ; ix<_seq->len ; ix++)
			{
			if ((ix != 0) && ((ix % perLine) == 0)) fprintf (f, "\n");
			fprintf (f, " %02X", _seq->v[ix]);
			}
		fprintf (f, "\n");
		}
	else
		{
		for (ix=0 ; ix<_seq->len ; ix++)
			{
			if ((ix != 0) && ((ix % perLine) == 0)) fprintf (f, "\n");
			if (_seq->v[ix] == 0) fprintf (f, "*");
			                 else fprintf (f, "%c", _seq->v[ix]);
			}
		fprintf (f, "\n");
		}
	}

//----------
//
// print_partition_table--
//	Print a sequence's partition table.  (for debugging only)
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write to.
//	seq*	seq:	The sequence to print.
//
// Returns:
//	(nothing)
//
//----------

void print_partition_table
   (FILE*			f,
	seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;

	if (sp->p == NULL)
		{ fprintf (f, "sequence has no partition\n");  return; }
	if (sp->state == seqpart_empty)
		{ fprintf (f, "partition table is in empty state\n");  return; }
	if (sp->state == seqpart_reusable)
		{ fprintf (f, "partition table is in re-usable state\n");  return; }
	if (sp->state == seqpart_loading)
		{ fprintf (f, "partition table is in loading state\n");  return; }

	fprintf (f, "     %8s %8s %8s %8s %6s %s\n",
	            "sepBefore", "sepAfter", "startLoc", "trueLen", "contig", "header");

	p = sp->p;
	for (ix=0 ; ix<=sp->len ; ix++)
		{
		fprintf (f, "[%2d] %8u %8u " unsposStarFmt " " unsposStarFmt,
		            ix, p[ix].sepBefore, p[ix].sepAfter,
		            8, p[ix].startLoc,
		            8, p[ix].trueLen);
		if (ix < sp->len)
			fprintf (f, " %6d %s", p[ix].contig, &sp->pool[p[ix].header]);
		fprintf (f, "\n");
		}

	}

//----------
//
// mask_sequence, mask_sequence_keep--
//	Mask a sequence, in place, as prescribed by some file.
//
// mask_sequence() replaces any base in the prescribed intervals.
// mask_sequence_keep() replaces any base NOT in the prescribed intervals.
//
// A typical masking file looks like this:
//
//	1527933 3184039
//	4165389 6877343
//	7374477 7902860
//
// Each line describes a region to be masked.  Indexes are one-based, and
// inclusive on both ends.  Numbers are free format.  Comment lines (beginning
// with #) are ignored, as are blank lines.  Additional information after the
// first two columns is also ignored.
//
// Note that if the sequence has been reverse complemented, the masking
// intervals are relative to the reverse strand.
//
//----------
//
// Arguments:
//	seq*	seq:			The sequence to mask.
//	char*	maskFilename:	The name of the file containing masking
//							.. information.
//	int		maskChar:		The character to ask as a mask.  Normally this is a
//							.. character (e.g. 'X' or 'N').  However, the value
//							.. -1 means that we should mask by changing to
//							.. lowercase.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

void mask_sequence
   (seq*			_seq,
	char*			maskFilename,
	int				_maskChar)
	{
	char			line[511+1], discard[511+1];
	char			maskChar = 0;
	int				toLower  = false;
	FILE*			maskF;
	int				len;
	int				lineNum, missingEol;
	char*			waffle;
	char			extra;
	int				numItems;
	seqpartition*	sp = &_seq->partition;
	partition*		part;
	unspos			b, e, pB, pE, pOffset, pLen;
	u32				ix;

	if (_maskChar == -1) toLower  = true;
	else                 maskChar = (u8) _maskChar;

	// read the masking intervals and deposit the masking character at every
	// contained base

	maskF = fopen_or_die (maskFilename, "rt");

	lineNum = 0;
	while (fgets (line, sizeof(line), maskF) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line (and discard the extra)

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), maskF) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it;  parse it as an interval

		numItems = sscanf (line, unsposFmtScanf " " unsposFmtScanf "%c",
		                   &b, &e, &extra);

		if ((numItems == 3) && ((extra == ' ') || (extra == '\t')))
			numItems = 2;

		if (numItems != 2)
			{
			char* scan = skip_whitespace (line);
			numItems = 0;
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }

			if (numItems == 3)
				suicidef ("bad interval (in %s, line %d): \"%s\"\n"
				          "three-column masking intervals are not supported for this operation",
				          maskFilename, lineNum, line);
			else
				suicidef ("bad interval (in %s, line %d): \"%s\"",
				          maskFilename, lineNum, line);
			}

		// trim the left end of the interval to our subsequence 

		if (e < _seq->startLoc) continue;
		if (b < _seq->startLoc) b = _seq->startLoc;
		b -= _seq->startLoc;    // (nota bene, b is zero-based afterwards)
		e -= _seq->startLoc-1;  // (nota bene, e is open interval afterwards)

		if (sp->p == NULL)				//=== sequence is not partitioned ===
			{
			// trim the right end of the interval to our subsequence

			if (b >= _seq->len) continue;
			if (e >= _seq->len) e = _seq->len;
			if (e <= b)         continue;

			// mask the interval

			if (toLower)
				{
				for ( ; b<e ; b++)
					{
					if ((_seq->v[b] >= 'A') && (_seq->v[b] <= 'Z'))
						_seq->v[b] = _seq->v[b] + 'a' - 'A';
					}
				}
			else
				memset (_seq->v+b, maskChar, (size_t) (e-b));
			}
		else							//=== sequence is partitioned ===
			{
			for (ix=0 ; ix<sp->len ; ix++)
				{
				part    = &sp->p[ix];
				pOffset = part->sepBefore + 1;
				pLen    = (part+1)->sepBefore - pOffset;

				// trim the right end of the interval to our subsequence

				pB = b;
				pE = e;

				if (pB >= pLen) continue;
				if (pE >= pLen) pE = pLen;
				if (pE <= pB)   continue;

				// mask the interval

				pB += pOffset;
				pE += pOffset;

				if (toLower)
					{
					for ( ; pB<pE ; pB++)
						{
						if ((_seq->v[pB] >= 'A') && (_seq->v[pB] <= 'Z'))
							_seq->v[pB] = _seq->v[pB] + 'a' - 'A';
						}
					}
				else
					memset (_seq->v+pB, maskChar, (size_t) (pE-pB));
				}
			}
		}

	fclose(maskF);
	}

void mask_sequence_keep
   (seq*			_seq,
	char*			maskFilename,
	int				_maskChar)
	{
	char			line[511+1], discard[511+1];
	char			maskChar = 0;
	int				toLower  = false;
	FILE*			maskF;
	int				len;
	int				lineNum, missingEol;
	char*			waffle;
	char			extra;
	int				numItems;
	seqpartition*	sp = &_seq->partition;
	partition*		part;
	unspos			b, e, pB, pE, pOffset, pLen;
	u32				ix;

	if ((_seq->fileType != seq_type_fasta)
	 && (_seq->fileType != seq_type_fastq)
	 && (_seq->fileType != seq_type_nib)
	 && (_seq->fileType != seq_type_2bit)
	 && (_seq->fileType != seq_type_hsx))
		suicidef ("masking of interval complements only valid for DNA sequences\n"
		          " (%s is %s file)",
		          sequence_filename(_seq), seqTypeNames[_seq->fileType]);

	if (_maskChar == -1) toLower  = true;
	else                 maskChar = (u8) _maskChar;

	// read the masking intervals and mark the most-significant bit of every
	// contained base

	maskF = fopen_or_die (maskFilename, "rt");

	lineNum = 0;
	while (fgets (line, sizeof(line), maskF) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way);
		// if the line was split we simply read ahead until we find the end of
		// the line (and discard the extra)

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		if (missingEol)
			{
			while (fgets (discard, sizeof(discard), maskF) != NULL)
				{
				len = strlen(discard);
				if (len == 0) break;
				if (discard[len-1] == '\n') break;
				}
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// ok, the line has something in it;  parse it as an interval

		numItems = sscanf (line, unsposFmtScanf " " unsposFmtScanf "%c",
		                   &b, &e, &extra);

		if ((numItems == 3) && ((extra == ' ') || (extra == '\t')))
			numItems = 2;

		if (numItems != 2)
			{
			char* scan = skip_whitespace (line);
			numItems = 0;
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }
			if (*scan != 0) { scan = skip_darkspace (scan);  scan = skip_whitespace (scan);  numItems++; }

			if (numItems == 3)
				suicidef ("bad interval (in %s, line %d): \"%s\"\n"
				          "three-column masking intervals are not supported for this operation",
				          maskFilename, lineNum, line);
			else
				suicidef ("bad interval (in %s, line %d): \"%s\"",
				          maskFilename, lineNum, line);
			}

		// trim the left end of the interval to our subsequence 

		if (e < _seq->startLoc) continue;
		if (b < _seq->startLoc) b = _seq->startLoc;
		b -= _seq->startLoc;    // (nota bene, b is zero-based afterwards)
		e -= _seq->startLoc-1;  // (nota bene, e is open interval afterwards)

		if (sp->p == NULL)				//=== sequence is not partitioned ===
			{
			// trim the right end of the interval to our subsequence

			if (b >= _seq->len) continue;
			if (e >= _seq->len) e = _seq->len;

			// mark the interval

			while (b<e) _seq->v[b++] |= 0x80;
			}
		else							//=== sequence is partitioned ===
			{
			for (ix=0 ; ix<sp->len ; ix++)
				{
				part    = &sp->p[ix];
				pOffset = part->sepBefore + 1;
				pLen    = (part+1)->sepBefore - pOffset;

				// trim the right end of the interval to our subsequence

				pB = b;
				pE = e;

				if (pB >= pLen) continue;
				if (pE >= pLen) pE = pLen;

				// mark the interval

				pB += pOffset;
				pE += pOffset;

				while (pB<pE) _seq->v[pB++] |= 0x80;
				}
			}

		}

	fclose(maskF);

	// scan the sequence, replacing unmarked bases with the masking character,
	// and removing the marks

	for (b=0 ; b<_seq->len ; b++)
		{
		if (_seq->v[b] == 0) continue;
		if ((_seq->v[b] & 0x80) != 0)	// marked => erase mark
			_seq->v[b] &= ~0x80;
		else if (!toLower)				// unmarked => mask it
			_seq->v[b] = maskChar;
		else							// unmarked => change to lower case
			{
			if ((_seq->v[b] >= 'A') && (_seq->v[b] <= 'Z'))
				_seq->v[b] = _seq->v[b] + 'a' - 'A';
			}
		}

	}

//----------
//
// colorize_sequence--
//	Create a sequence's color-space equivalent (see description below).
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to modify.
//
// Returns:
//	(nothing)
//
//----------
//
// We colorize by mapping each pair of nucleotides to a single "color" value
// according to the table below.  This extends the normal color space
// definition to account for lower case, N, and other characters, in a way that
// makes sense in the context of lastz's alignment operations.  The DNA sequence
// may contain upper and lower case nucleotides, N's, and other stuff.  Further,
// a partitioned sequence will contain NUL character separators (shown as '*' in
// the table below).  We also treat the DNA sequence as if it had an 'X'
// prepended to it.
//
//		                  second in pair
//		      | A | C | G | T | a | c | g | t | N/n | other | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		   A  | 0 | 1 | 2 | 3 | 0 | 1 | 2 | 3 |  N  |   X   | * |
//		   C  | 1 | 0 | 3 | 2 | 1 | 0 | 3 | 2 |  N  |   X   | * |
//		   G  | 2 | 3 | 0 | 1 | 2 | 3 | 0 | 1 |  N  |   X   | * |
//	first  T  | 3 | 1 | 2 | 0 | 3 | 1 | 2 | 0 |  N  |   X   | * |
//	 in	   a  | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |  N  |   X   | * |
//	  pair c  | 1 | 0 | 3 | 2 | 5 | 4 | 7 | 6 |  N  |   X   | * |
//		   g  | 2 | 3 | 0 | 1 | 6 | 7 | 4 | 5 |  N  |   X   | * |
//		   t  | 3 | 1 | 2 | 0 | 7 | 6 | 5 | 4 |  N  |   X   | * |
//		  N/n | N | N | N | N | N | N | N | N |  N  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		other | X | X | X | X | X | X | X | X |  X  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//		   *  | X | X | X | X | X | X | X | X |  X  |   X   | * |
//		------+---+---+---+---+---+---+---+---+-----+-------+---+
//
// Example:
//
//	 dna:   G T C G A A C C C G * C A A C C G T A T T * T A A T A A G T T T
//	 color: X 1 2 3 2 0 1 0 0 3 * X 1 0 1 0 3 1 3 3 0 * X 3 0 3 3 0 2 1 0 0 
//
//----------

static void do_colorize (u8* seq, u8* colorSeq, unspos seqLen);

void colorize_sequence
   (seq*	_seq)
	{
	char*	name;

	if (_seq == NULL) suicide ("colorize_sequence(NULL)");
	if (_seq->len < 1) return;

	if (!_seq->vcOwner)
		{
		name = (_seq->filename != NULL)? _seq->filename : _seq->header;
		suicidef ("internal error, attempt to colorize external sequence (%s)",
		          name);
		}

	if (_seq->vc != NULL)
		{
		name = (_seq->filename != NULL)? _seq->filename : _seq->header;
		suicidef ("internal error, attempt to re-colorize sequence (%s)",
		          name);
		}

	_seq->vc = malloc_or_die ("colorize_sequence (vc)", _seq->len+1);
	do_colorize (_seq->v, _seq->vc, _seq->len);
	}

#define A_ 0
#define C_ 1
#define G_ 2
#define T_ 3
#define a_ 4
#define c_ 5
#define g_ 6
#define t_ 7
#define N_ 8
#define X_ 9
#define __ X_
#define Z_ 10

const s8 nuc_to_color_bits[256] =
	{
	Z_,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,A_,__,C_,__,__,__,G_,__,__,__,__,__,__,N_,__, // 4x (upper case)
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,a_,__,c_,__,__,__,g_,__,__,__,__,__,__,N_,__, // 6x (lower case)
	__,__,__,__,t_,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__	 // Fx
	};

static void do_colorize
   (u8*		dnaSeq,
	u8*		colorSeq,
	unspos	seqLen)
	{
	u8		bits1, bits2;

	bits1 = nuc_to_color_bits[X_];

	while (seqLen-- > 0)
		{
		bits2 = nuc_to_color_bits[*(dnaSeq++)];

		if (bits1 == Z_)
			*(colorSeq++) = 0;
		else if ((bits1 == X_) || (bits2 == X_))
			*(colorSeq++) = 'X';
		else if ((bits1 == N_) || (bits2 == N_))
			*(colorSeq++) = 'N';
		else if ((bits1 <= T_) || (bits2 <= T_))
			*(colorSeq++) = '0' + ((bits1 ^ bits2) & 3);
		else
			*(colorSeq++) = '4' + ((bits1 ^ bits2) & 3);

		bits1 = bits2;
		}

	*colorSeq = 0;
	}

#undef A_
#undef C_
#undef G_
#undef T_
#undef a_
#undef c_
#undef g_
#undef t_
#undef N_
#undef X_
#undef __
#undef Z_

//----------
//
// validate_rev_comp--
//	Validate that a sequence is in the same direction as it would normally be
//	after operations were performed at initial load.  If it is not, it is
//	reversed and/or complemented as needed to bring it back to that state.
//
//		v IS v      should be ->> | F  | C  | R  | RC |
//		--------------------------+----+----+----+----+
//		forward              (F)  | -  | C  | R  | RC | <<- correction
//		complemented         (C)  | C  | -  | RC | R  |
//		reversed             (R)  | R  | RC | -  | C  |
//		reverse-complemented (RC) | RC | R  | C  | -  |
//		--------------------------+----+----+----+----+
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to operate upon.
//
// Returns:
//	(nothing)
//
//----------

void validate_rev_comp
   (seq*	_seq)
	{
	int		operation = _seq->revCompFlags ^ _seq->doRevCompFlags;

	if (operation == rcf_revcomp)
		rev_comp_sequence (_seq, _seq->qToComplement);
	else if (operation == rcf_rev)
		backward_sequence (_seq);
	else if (operation == rcf_comp)
		{
		backward_sequence (_seq);
		rev_comp_sequence (_seq, _seq->qToComplement);
		}
	}

//----------
//
// rev_comp_sequence--
//	Convert a sequence to its reverse-complement, in place.  Note that in
//	partitioned sequences, each partition is reverse-complemented separately
//	(independent of the others), so that the partition table can remain
//	unchanged.
//
//----------
//
// Arguments:
//	seq*		seq:				The sequence to reverse-complement.
//	const u8*	nucToComplement:	Array to map a nucleotide to its complement.
//									.. If this is NULL, nuc_to_complement is
//									.. used for DNA sequences.  For quantum
//									.. sequences this cannot be NULL.
//
// Returns:
//	(nothing)
//
//----------

static void revcomp_in_place (u8* seq, unspos seqLen, const u8* nucToComplement);
static void reverse_in_place (u8* seq, unspos seqLen);


void rev_comp_sequence
   (seq*			_seq,
	const u8*		_nucToComplement)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;
	const u8*		nucToComplement;

	if (_seq == NULL) suicide ("rev_comp_sequence(NULL)");
	if (_seq->len < 1) return;

	if ((_seq->fileType == seq_type_qdna)
	 && (_nucToComplement == NULL))
		{
		suicidef ("quantum DNA cannot be complemented (%s)\n"
		          "(the score file lacks complements)",
		          sequence_filename(_seq));
		return; // (we never reach here)
		}

	if (_nucToComplement != NULL) nucToComplement = _nucToComplement;
	                         else nucToComplement = nuc_to_complement;

	if (sp->p == NULL)
		{
		revcomp_in_place (_seq->v, _seq->len, nucToComplement);
		if (_seq->vq != NULL)
			reverse_in_place (_seq->vq, _seq->len);
		}
	else
		{
		p = sp->p;
		for (ix=0 ; ix<sp->len ; ix++)
			{
			revcomp_in_place (/*start*/  _seq->v + p[ix].sepBefore+1,
			                  /*length*/ p[ix].sepAfter - (p[ix].sepBefore+1),
			                  /*how*/    nucToComplement);
			if (_seq->vq != NULL)
				reverse_in_place (/*start*/  _seq->vq + p[ix].sepBefore+1,
				                  /*length*/ p[ix].sepAfter - (p[ix].sepBefore+1));
			}
		}

	if (_seq->vc != NULL)
		do_colorize (_seq->v, _seq->vc, _seq->len);

	_seq->revCompFlags = _seq->revCompFlags ^ rcf_revcomp;
	}


static void revcomp_in_place
   (u8*			_seq,
	unspos		seqLen,
	const u8*	nucToComplement)
	{
	u8*		left, *right;
	u8		nuc;

	left  = _seq;
	right = left + seqLen-1;
	for ( ; left<=right ; left++,right--)
		{
		nuc    = nucToComplement[*left ];
		*left  = nucToComplement[*right];
		*right = nuc;
		}
	}


static void reverse_in_place
   (u8*		_seq,
	unspos	seqLen)
	{
	u8*		left, *right;
	u8		nuc;

	left  = _seq;
	right = left + seqLen-1;
	for ( ; left<=right ; left++,right--)
		{ nuc = *left;  *left = *right;  *right = nuc; }
	}

//----------
//
// backward_sequence--
//	Convert a sequence to its reverse (without complement), in place.  Note
//	that in partitioned sequences, each partition is reversed separately.
//
// [[ see also copy_reverse_of_string and strncpy_reverse ]]
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to reverse.
//
// Returns:
//	(nothing)
//
//----------

static void backward_in_place (u8* seq, unspos seqLen);

void backward_sequence
   (seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	partition*		p;
	u32				ix;

	if (_seq->fileType == seq_type_csfasta)
		suicidef ("color space cannot be reversed (%s)", sequence_filename(_seq));

	if (_seq == NULL) suicide ("backward_sequence(NULL)");
	if (_seq->len < 1) return;

	if (sp->p == NULL)
		backward_in_place (_seq->v, _seq->len);
	else
		{
		p = sp->p;
		for (ix=0 ; ix<sp->len ; ix++)
			backward_in_place (/*start*/  _seq->v + p[ix].sepBefore+1,
			                   /*length*/ p[ix].sepAfter - (p[ix].sepBefore+1));
		}

	_seq->revCompFlags = _seq->revCompFlags ^ rcf_rev;
	}


static void backward_in_place
   (u8*		_seq,
	unspos	seqLen)
	{
	u8*		left, *right;
	u8		nuc;

	left  = _seq;
	right = left + seqLen-1;
	for ( ; left<=right ; left++,right--)
		{ nuc = *left;  *left = *right;  *right = nuc; }
	}

//----------
//
// upper_sequence--
//	Convert a sequence to its upper-case equivalent, in place.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to modify.
//
// Returns:
//	(nothing)
//
//----------

static void upper_in_place (u8* seq, unspos seqLen);

void upper_sequence
   (seq*	_seq)
	{
	if (_seq == NULL) suicide ("upper_sequence(NULL)");
	if (_seq->len < 1) return;

	if (_seq->fileType == seq_type_csfasta)
		suicidef ("color space cannot be upper-cased (%s)", sequence_filename(_seq));
	else if (_seq->fileType == seq_type_qdna)
		suicidef ("quantum DNA cannot be upper-cased (%s)", sequence_filename(_seq));

	upper_in_place (_seq->v, _seq->len);
	}

static void upper_in_place
   (u8*		_seq,
	unspos	seqLen)
	{
	u8*		left, *right;

	left  = _seq;
	right = left + seqLen;
	for ( ; left<right ; left++)
		{
		if (('a' <= *left) && (*left <= 'z'))
			*left += 'A'-'a';
		}
	}

//----------
//
// copy_reverse_of_string--
//	Make a copy of a string, in reversed order.
//
// [[ see also backward_sequence and strncpy_reverse ]]
//
//----------
//
// Arguments:
//	char*	s:		The string to copy.
//	unspos	len:	Its length (not including the terminating zero).
//
// Returns:
//	The copy, newly allocated from the heap (including a terminating zero).
//	Caller must eventually dispose of this with a call to free().
//
//----------

char* copy_reverse_of_string
   (char*	s,
	unspos	len)
	{
	char* r = malloc_or_die ("copy_reverse_of_string", (size_t)(len+1));
	char* t;

	for (t=r+len-1 ; t>=r ; t--)
		*t = *s++;

	r[len] = 0;
	return r;
	}

//----------
//
// strncpy_reverse--
//	Copy of a string, in reversed order.
//
// [[ see also backward_sequence and copy_reverse_of_string ]]
//
//----------
//
// Arguments:
//	char*	d:		The place to copy the string to.
//	char*	s:		The string to copy.
//	unspos	len:	Its length (not including the terminating zero).
//
// Returns:
//	(nothing)
//	The copy, newly allocated from the heap (including a terminating zero).
//	Caller must eventually dispose of this with a call to free().
//
//----------

void strncpy_reverse
   (char*	d,
	char*	s,
	unspos	len)
	{
	char* t;

	for (t=d+len-1 ; t>=d ; t--)
		*t = *s++;

	d[len] = 0;
	}

//----------
//
// fence_sequence_interval--
//	Place two markers in a sequence, bracketing the ends of an interval.
//
// Note: the markers can later be removed by calling unfence_sequence_interval.
//
//----------
//
// Arguments:
//	seq*		seq:		The sequence to modify.
//	interval	interval:	The interval to mark.  We expect this to be origin-
//							.. zero, half-open.  We mark the positions s-1
//							.. and e.  It is legal for the interval to indicate
//							.. points beyond the sequence.
//	u8			ch:			The character to mark with.
//
// Returns:
//	(nothing)
//
//----------

void fence_sequence_interval
   (seq*		_seq,
	interval	_interval,
	u8			ch)
	{
	unspos		s, e;
	u8			sCh, eCh;

	if (_seq == NULL) suicide ("fence_sequence_interval(NULL)");

	if ((_seq->hasLeftFence) || (_seq->hasRightFence))
		suicide ("INTERNAL ERROR-- sequence already has fences");

	debugFencing_1

	s = _interval.s;
	if (s >= 1)
		{
		s--;
		sCh = _seq->v[s];
		_seq->v[s] = ch;
		_seq->hasLeftFence = true;
		_seq->leftFencePos = s;
		_seq->leftFenceCh  = sCh;
		debugFencing_2
		}

	e = _interval.e;
	if (e <= _seq->len)
		{
		eCh = _seq->v[e];
		_seq->v[e] = ch;
		_seq->hasRightFence = true;
		_seq->rightFencePos = e;
		_seq->rightFenceCh  = eCh;
		debugFencing_3
		}

	debugFencing_4
	}

//----------
//
// unfence_sequence_interval--
//	Remove the markers placed by fence_sequence_interval.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to modify.
//
// Returns:
//	(nothing)
//
//----------

void unfence_sequence_interval
   (seq*		_seq)
	{
	if ((!_seq->hasLeftFence) && (!_seq->hasRightFence))
		suicide ("INTERNAL ERROR-- sequence has no fences to tear down");

	debugFencing_5

	if (_seq->hasLeftFence)
		{
		_seq->v[_seq->leftFencePos] = _seq->leftFenceCh;
		_seq->hasLeftFence = false;
		debugFencing_6
		}

	if (_seq->hasRightFence)
		{
		_seq->v[_seq->rightFencePos] = _seq->rightFenceCh;
		_seq->hasRightFence = false;
		debugFencing_7
		}

	debugFencing_8
	}

//----------
//
// parse_sequence_name--
//	Parse a sequence name
//
// The seqence name is the name of a file, plus some control options.  The
// basic form is
//
//	nickname::sequence_file/contig_name{mask_file}[actions]-
//
//----------
//
// Arguments:
//	const char*	name:				The name string to parse
//	char**		filename:			Place to return the file name.
//	char**		nickname:			Place to return the nickname, if any.
//	char**		contigOfInterest:	Place to return the name of the particular
//									contig-of-interest, if any.
//	char**		namesFilename:		Place to return the contigs-of-interest
//									.. file name, if any.
//	char**		choresFilename:		Place to return the chores file name, if
//									.. any.  Note that, unlike the other
//									.. filename arguments, this *must* be set
//									.. to either NULL or a filename upon entry.
//	int*		subsampleK:			Place to return K-of-N subsampling
//	int*		subsampleN:			.. specification.
//	char**		softMaskFilename:	Place to return the soft-mask file name, if
//									.. any.
//	int*		softMaskComplement:	Place to return true/false for soft-masking.
//	char**		xMaskFilename:		Place to return the x-mask file name, if any.
//	int*		xMaskComplement:	Place to return true/false for x-masking.
//	char**		nMaskFilename:		Place to return the n-mask file name, if any.
//	int*		xMaskComplement:	Place to return true/false for n-masking.
//	int*		nameParseType:		Place to return the name parse type, if any.
//	char**		nameTrigger:		Place to return the mask name trigger, if
//									.. any.
//	int*		doRevCompFlags:		Place to return rcf_forward/rcf_revcomp.
//									.. Actually can return any combination of
//									.. rcf_xxx flags.
//	int*		doUnmask:			(see field in the sequence structure).
//	int*		doPartitioning:		(see field in the sequence structure).
//	int*		doJoin:				(see field in the sequence structure).
//	char*		separatorCh:		Place to return a separator character if
//									.. any is specified.
//	int*		useFullNames:		(see field in the sequence structure).
// 	int*		fileType:			Place to return file type if any is
//									.. specified (one of seq_type_xxx).
//	int*		isQuantum:			Place to return true/false for whether the
//									.. sequence is to be quantum DNA.
//	char**		qCodingFilename:	Place to return the quantum coding file
//									.. name, if any.
//	unspos*		start:				Place to return the starting index (origin-1).
//									.. If no start is specified, 0 is placed
//									.. here. (see note 1 below)
//	unspos*		end:				Place to return the ending index (inclusive).
//									.. If no end is specified, 0 is placed here.
//									.. (see note 1 below)
//	int*		endIsSoft:			Place to return a flag telling the caller
//									.. that, while the end has been specified,
//									.. it is a soft specification and the caller
//									.. can trim it to the actual end of the
//									.. sequence.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------
//
// Notes:
//
// (1)	Start and end are origin-1, inclusive on both ends.  So, for example,
//		start=10 and end=15 defines a 6 letter sequence, the 10th thru 15th
//		letters of the file.  The first 9 characters from the file should be
//		skipped.
//
//----------

void print_file_actions (FILE* f)
	{
	int   typeIx, indent, nameLen, needComma, lineWidth;
	char* action, *description, *name;

	fprintf (f, "Supported actions:\n");
	fprintf (f, "  <subrange>          only process a subrange of the file (see below)\n");
	fprintf (f, "  revcomp             reverse complement\n");
	fprintf (f, "  multiple            file's sequences are internally treated as a single\n");
	fprintf (f, "                      sequence\n");
	fprintf (f, "  separator=<ch>      file's sequences are internally separated by the given\n");
	fprintf (f, "                      character;  no alignments will cross a separator\n");
	fprintf (f, "                      (this forces multiple)\n");
	fprintf (f, "  subset=<namesfile>  process only the sequences listed in namesfile\n");
	fprintf (f, "                      (only valid for fasta, fastq, 2bit and hsx)\n");
	fprintf (f, "  chores=<choresfile> process \"alignment chores\" listed in choresfile\n");
	fprintf (f, "                      (only valid for fasta, fastq, 2bit and hsx)\n");
	fprintf (f, "  subsample=<k>/<n>   process only the kth sequence of every group of n\n");
	fprintf (f, "                      sequences.  k ranges from 1 to n\n");
	fprintf (f, "                      (only valid for fasta, 2bit and hsx)\n");
	fprintf (f, "  unmask              convert any lowercase bases to uppercase\n");
	fprintf (f, "  softmask=<file>     mask segments specified in <file>, replacing them with\n");
	fprintf (f, "                      lowercase equivalents\n");
	fprintf (f, "  softmask=keep:<file> mask bases NOT in segments specified in <file>, with Xs\n");
	fprintf (f, "  xmask=<file>        mask segments specified in <file>, replacing them with Xs\n");
	fprintf (f, "  xmask=keep:<file>   mask bases NOT in segments specified in <file>, with Xs\n");
	fprintf (f, "  nmask=<file>        mask segments specified in <file>, replacing them with Ns\n");
	fprintf (f, "  nmask=keep:<file>   mask bases NOT in segments specified in <file>, with Ns\n");
	fprintf (f, "  nickname=<name>     name to use for this sequence in any output files\n");
	fprintf (f, "  nameparse=full      report full names in alignments instead of short names\n");
	fprintf (f, "  nameparse=alphanum  pull short name from sequence header, alphanumeric only\n");
	fprintf (f, "  nameparse=darkspace pull short name from sequence header, non-whitespace only\n");
	fprintf (f, "  nameparse=tag:<marker> pull a short name from sequence header, starting from\n");
	fprintf (f, "                      marker (only valid for fasta)\n");
	fprintf (f, "  quantum             the sequence contains quantum DNA\n");
	fprintf (f, "  quantum=<codefile>  the sequence contains quantum DNA, and <codefile>\n");
	fprintf (f, "                      describes the mapping from symbols to probabilities (only\n");
	fprintf (f, "                      meaningful for --format=text)\n");

	action      = "  format=<type>       ";
	description = "override auto-format detect;  <type> is one of ";
	fprintf (f, "%s%s", action, description);
	indent    = strlen (action);
	lineWidth = indent + strlen (description);
	for (typeIx=seq_type_unknown+1 ; typeIx<seq_type_max ; typeIx++)
		{
		name    = seqTypeNames[typeIx];
		nameLen = strlen(name);

		needComma = true;
		if (typeIx == seq_type_unknown+1)
			needComma = false;
		else if (lineWidth + 2 + nameLen >= 79)
			{
			fprintf (f, ",\n%*s", indent, " ");
			lineWidth = indent;
			needComma = false;
			}

		if (needComma) fprintf (f, ", ");
		fprintf (f, "%s", name);
		if (needComma) lineWidth += 2;
		lineWidth += nameLen;
		}
	fprintf (f, "\n\n");

	fprintf (f, "Subranges:\n");
	fprintf (f, "  start,end           same as start..end (for BLASTZ compatibility)\n");
	fprintf (f, "  start..end          process from start thru end, inclusive\n");
	fprintf (f, "  start..             process from given start thru the end of the sequence\n");
	fprintf (f, "  ..end               process from the start of the sequence thru given end\n");
	fprintf (f, "  start#length        same as start..start+length-1\n");
	fprintf (f, "  center^length       same as center-length/2..center+length/2-1\n");
	fprintf (f, "  start..end+zoom%%    process from start thru end, zoomed out by zoom%%\n");
	fprintf (f, "  (subrange indices begin with 1 and are inclusive)\n");
	}


//--- parse_sequence_name--

static void parse_sequence_name
   (const char*	name,
	char**		filename,
	char**		nickname,
	char**		contigOfInterest,
	char**		namesFilename,
	char**		choresFilename,
	int*		subsampleK,
	int*		subsampleN,
	char**		softMaskFilename,
	int*		softMaskComplement,
	char**		xMaskFilename,
	int*		xMaskComplement,
	char**		nMaskFilename,
	int*		nMaskComplement,
	int*		nameParseType,
	char**		nameTrigger,
	int*		doRevCompFlags,
	int*		doUnmask,
	int*		doPartitioning,
	int*        doJoin,
	char*		separatorCh,
	int*		useFullNames,
	int*		fileType,
	int*		isQuantum,
	char**		qCodingFilename,
	unspos*		_start,
	unspos*		_end,
	int*		endIsSoft)
	{
	int			len;
	char*		fname, *bracket, *mask, *actions, *action, *actionName;
	char*		parse, *slashParse, *extParse;
	int			numItems, charsUsed;
	unspos		start, end, pendingStart, temp, length, mid;
	int			tempInt;
	float		size, zoom, fLength;
	int			parsed;

	*namesFilename    = NULL;
	*contigOfInterest = NULL;
	*qCodingFilename  = NULL;

	*_start = *_end = 0;
	*endIsSoft = false;

	*doRevCompFlags = rcf_forward;
	*doUnmask       = false;
	*doPartitioning = false;
	*doJoin         = false;
	*separatorCh    = 0;
	*useFullNames   = false;
	*fileType       = seq_type_unknown;
	*isQuantum      = false;

	//////////
	// copy the name, splitting out the nickname if present;  we will shorten
	// this copy if other components are present, so we are potentially
	// allocating more memory for the copy than is really needed
	//////////

	if (name == NULL) suicide ("parse_sequence_name(NULL)");

	parse   = strstr (name, "::");
	actions = strchr (name, '[');
	if ((parse == NULL)								// no "::"
	 || ((actions != NULL) && (parse > actions)))	// "::" is after "["
		{
		*nickname = NULL;
		*filename = fname = copy_string (name);
		}
	else
		{
		if (parse-name == 0) goto empty_species_name;
		*nickname = copy_prefix (name, parse-name);
		*filename = fname = copy_string (parse+2);
		}

	len = strlen (fname);
	if (len < 1) goto empty_file_name;

	//////////
	// see if we are to reverse the sequence
	//////////

	switch (fname[len-1])
		{
		case '-': *doRevCompFlags = rcf_revcomp;  fname[--len] = 0;  break;
		case '+':                                 fname[--len] = 0;  break;
		}

	if (len < 1) goto empty_file_name;

	//////////
	// split the file name string into its components
	//////////

	mask    = strchr (fname, '{');
	bracket = strchr (fname, '[');
	if ((bracket != NULL) && (mask != NULL))
		{ if (mask > bracket) mask = NULL; }

	if (mask == fname) goto empty_file_name;

	if (mask != NULL)
		*(mask++) = 0;

	parse = (mask == NULL)? fname : mask;

	actions = strchr (parse, '[');
	if (actions == parse) goto empty_file_name;

	if (actions != NULL)
		*(actions++) = 0;

	//////////
	// parse the mask file name
	//////////

	*softMaskFilename = NULL;
	*xMaskFilename    = NULL;
	*nMaskFilename    = NULL;

	if (mask != NULL)
		{
		len = strlen (mask);
		if (mask[--len] != '}') goto bad_mask;
		if (len == 0)           goto empty_mask_file_name;
		mask[len] = 0;
		*xMaskFilename = copy_string (mask);
		}

	//////////
	// split out the contig-of-interest if present
	//////////

	slashParse = NULL;

	extParse = strstr (fname, ".2bit/");
	if (extParse != NULL)
		slashParse = extParse+5;
	else
		{
		extParse = strstr (fname, ".hsx/");
		if (extParse != NULL)
			slashParse = extParse+4;
		}

	if (slashParse != NULL)
		{
		if      (strchr (slashParse+1, pathSlash) != NULL) extParse = NULL;
		else if (strstr (slashParse+1, ".2bit")   != NULL) extParse = NULL;
		}

	if (extParse != NULL)
		{
		*contigOfInterest = copy_string (slashParse+1);
		*slashParse = 0;
		}

	//////////
	// parse the actions list
	//////////

	if (actions != NULL)
		{
		len = strlen (actions);
		if (len == 0)
			goto bad_action_list;
		else if (actions[len-1] != ']')
			{
			if (strchr (actions, ']') == NULL) goto bad_action_list;
			                              else goto actions_not_at_end;
			}
		actions[--len] = 0;

		start = end = pendingStart = 0;
		while (actions != NULL)
			{
			//fprintf(stderr,"actions=\"%s\"\n", actions);
			action  = actions;
			actions = strchr (actions, ',');
			if (actions != NULL)
				*(actions++) = 0;
			else
				{
				actions = strstr (action, "][");
				if (actions != NULL)
					{ *(actions++) = 0;  actions++; }
				}

			//fprintf(stderr,"  action=\"%s\"\n", action);
			len = strlen (action);
			if (len == 0) goto blank_action;

			// parse simple actions

			if (strcmp (action, "unmask") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doUnmask = true;
				continue;
				}

			if (strcmp (action, "revcomp") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doRevCompFlags ^= rcf_revcomp;
				continue;
				}

			if ((strcmp (action, "backward")  == 0)  // (unadvertised)
			 || (strcmp (action, "backwards") == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doRevCompFlags ^= rcf_rev;
				continue;
				}

			if ((strcmp (action, "multi")    == 0)
			 || (strcmp (action, "multiple") == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*doPartitioning = *doJoin = true;
				continue;
				}

			if (strcmp_prefix (action, "sep=") == 0)
				{
				actionName = action + strlen("sep=");
				goto action_separator;
				}

			if (strcmp_prefix (action, "separator=") == 0)
				{
				actionName = action + strlen("separator=");
			action_separator:
				if (pendingStart != 0) goto unsatisfied_start;
				if (actionName[0] == 0) goto bad_separator;
				if (actionName[1] != 0) goto bad_separator;
				if (*separatorCh != 0) goto many_separators;
				*separatorCh = actionName[0];
				*doPartitioning = true;
				continue;
				}

			if ((strcmp (action, "nameparse=full") == 0)
			 || (strcmp (action, "fullname")       == 0)
			 || (strcmp (action, "fullnames")      == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				*useFullNames = true;
				continue;
				}

			if (action[0] == '@')
				{
				actionName = action + strlen("@");
				goto action_subset;
				}

			if (strcmp_prefix (action, "subset=") == 0)
				{
				actionName = action + strlen("subset=");
			action_subset:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*namesFilename != NULL) goto many_name_files;
				if (strlen(actionName) == 0) goto bad_name_file;
				*namesFilename = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "chores=") == 0)
				{
				actionName = action + strlen("chores=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*choresFilename != NULL) goto many_chore_files;
				if (strlen(actionName) == 0) goto bad_chore_file;
				*choresFilename = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "subsample=") == 0)
				{
				actionName = action + strlen("subsample=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*subsampleN != 0) goto many_subsamples;
				if (strlen(actionName) == 0) goto bad_subsample;
				slashParse = strchr (actionName, '/');
				if (slashParse == NULL) goto bad_subsample;

				len = slashParse - actionName;
				*slashParse = ']';	// (write a sentinel for parsing K)
				charsUsed = -1;
				numItems = sscanf (actionName, "%d]%n", &tempInt, &charsUsed);
				if ((numItems != 1) || (charsUsed != len+1) || (tempInt < 1))
					{ *slashParse = '/';  goto bad_subsample; }
				*subsampleK = tempInt;

				*(slashParse++) = '/';
				len = strlen(slashParse);
				slashParse[len] = ']';	// (write a sentinel for parsing N)
				charsUsed = -1;
				numItems = sscanf (slashParse, "%d]%n", &tempInt, &charsUsed);
				if ((numItems != 1) || (charsUsed != len+1) || (tempInt < *subsampleK))
					{ slashParse[len] = 0;  goto bad_subsample; }
				*subsampleN = tempInt;
				continue;
				}

			if ((strcmp (action, "nameparse=alnum")    == 0)
			 || (strcmp (action, "nameparse=alphanum") == 0)
			 || (strcmp (action, "name:alnum")         == 0)
			 || (strcmp (action, "name:alphanum")      == 0))
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_parse_types;
				*nameParseType = name_parse_type_alnum
				               | (*nameParseType & name_parse_fill_white);
				continue;
				}

			if (strcmp (action, "nameparse=darkspace") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_parse_types;
				*nameParseType = name_parse_type_darkspace
				               | (*nameParseType & name_parse_fill_white);
				continue;
				}

			if (strcmp_prefix (action, "nickname=") == 0)
				{
				actionName = action + strlen("nickname=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nickname != NULL)  goto many_nicknames;
				if (strlen(actionName) == 0) goto bad_nickname;
				*nickname = copy_string (actionName);
				continue;
				}

			if (strcmp_prefix (action, "name=") == 0)
				{
				actionName = action + strlen("name=");
				goto action_tag;
				}

			if (strcmp_prefix (action, "nameparse=tag:") == 0)
				{
				actionName = action + strlen("nameparse=tag:");
			action_tag:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nameTrigger != NULL) goto many_name_triggers;
				if (strlen(actionName) == 0) goto bad_name_trigger;
				if (parse_type(*nameParseType) != name_parse_type_core)
					goto many_name_parse_types;
				*nameParseType = name_parse_type_trigger 
				               | (*nameParseType & name_parse_fill_white);
				*nameTrigger   = copy_string (actionName);
				continue;
				}

			if (strcmp (action, "namejoin") == 0)
				{
				*nameParseType |= name_parse_fill_white;
				continue;
				}

			if (strcmp_prefix (action, "soft=keep:") == 0)
				{
				actionName = action + strlen("soft=keep:");
				goto action_softmaskkeep;
				}

			if (strcmp_prefix (action, "softmask=keep:") == 0)
				{
				actionName = action + strlen("softmask=keep:");
			action_softmaskkeep:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*softMaskFilename != NULL) goto many_soft_mask_files;
				if (strlen(actionName) == 0) goto bad_soft_mask_file;
				*softMaskFilename   = copy_string (actionName);
				*softMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "soft=") == 0)
				{
				actionName = action + strlen("soft=");
				goto action_softmask;
				}

			if (strcmp_prefix (action, "softmask=") == 0)
				{
				actionName = action + strlen("softmask=");
			action_softmask:
				if (pendingStart != 0) goto unsatisfied_start;
				if (*softMaskFilename != NULL) goto many_soft_mask_files;
				if (strlen(actionName) == 0) goto bad_soft_mask_file;
				*softMaskFilename   = copy_string (actionName);
				*softMaskComplement = false;
				continue;
				}

			if (strcmp_prefix (action, "xmask=keep:") == 0)
				{
				actionName = action + strlen("xmask=keep:");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*xMaskFilename != NULL) goto many_x_mask_files;
				if (strlen(actionName) == 0) goto bad_x_mask_file;
				*xMaskFilename   = copy_string (actionName);
				*xMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "xmask=") == 0)
				{
				actionName = action + strlen("xmask=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*xMaskFilename != NULL) goto many_x_mask_files;
				if (strlen(actionName) == 0) goto bad_x_mask_file;
				*xMaskFilename   = copy_string (actionName);
				*xMaskComplement = false;
				continue;
				}

			if (strcmp_prefix (action, "nmask=keep:") == 0)
				{
				actionName = action + strlen("nmask=keep:");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nMaskFilename != NULL) goto many_n_mask_files;
				if (strlen(actionName) == 0) goto bad_n_mask_file;
				*nMaskFilename   = copy_string (actionName);
				*nMaskComplement = true;
				continue;
				}

			if (strcmp_prefix (action, "nmask=") == 0)
				{
				actionName = action + strlen("nmask=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*nMaskFilename != NULL) goto many_n_mask_files;
				if (strlen(actionName) == 0) goto bad_n_mask_file;
				*nMaskFilename   = copy_string (actionName);
				*nMaskComplement = false;
				continue;
				}

			if (strcmp (action, "quantum") == 0)
				{
				if (pendingStart != 0) goto unsatisfied_start;
				if (*isQuantum) goto many_quantums;
				if (*fileType != seq_type_unknown) goto many_file_types;
				*isQuantum = true;
				*fileType = seq_type_qdna;
				continue;
				}

			if (strcmp_prefix (action, "quantum=") == 0)
				{
				actionName = action + strlen("quantum=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*isQuantum) goto many_quantums;
				if (strlen(actionName) == 0) goto bad_code_file;
				*qCodingFilename = copy_string (actionName);
				*isQuantum       = true;
				continue;
				}

			if (strcmp_prefix (action, "format=") == 0)
				{
				int fType, typeIx;
				actionName = action + strlen("format=");
				if (pendingStart != 0) goto unsatisfied_start;
				if (*fileType != seq_type_unknown) goto many_file_types;
				fType = seq_type_unknown;
				for (typeIx=seq_type_unknown+1 ; typeIx<seq_type_max ; typeIx++)
					{
					if (strcmp (actionName, seqTypeNames[typeIx]) == 0)
						{ fType = typeIx;  break; }
					}
				if (fType == seq_type_unknown) goto bad_file_format;
				*fileType = typeIx;
				continue;
				}

			// any other option would be a sequence limit;  so if we already
			// have sequence limits, this is an error

			if ((start != 0) || (end != 0))
				goto bad_action;

			// try to parse the field as a single number

			action[len] = ']';		// (write a sentinel;  this location is
									//  guaranteed to be part of the string)

			charsUsed = -1;
			numItems = sscanf (action, unsposFmtScanf "]%n", &temp, &charsUsed);
			if ((numItems == 1) && (charsUsed == len+1))
				{
				//fprintf(stderr,"  parsed as one (" unsposFmt ") %d\n", temp, len);
				if (temp == 0)
					{ action[len] = 0;  goto bad_sequence_position; }
				if (pendingStart == 0)
					pendingStart = temp;
				else
					{ start = pendingStart;  end = temp;  pendingStart = 0; }
				continue;
				}
			else if (pendingStart != 0)
				{ action[len] = 0;  goto unsatisfied_start; }

			// try to parse the field as sequence limits

			parsed = false;

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposDotsFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, dots\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposDotsFmtScanf "+%f%%]%n", &start, &end, &zoom, &charsUsed);
				if ((numItems == 3) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, dots\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					if (end < start)
						{ action[len] = 0;  goto bad_limits; }

					length  = end+1 - start;
					mid     = start + length/2;	// (avoids overflow, instead of
					zoom    = 1 + (zoom/100.0);	//  .. (start+end)/2 )
					fLength = length * zoom;
					length  = fLength;			// (ok if it overflows)
					if (mid - fLength/2 <= 1) start = 1;
					                     else start = 1 + (mid - length/2);
					if (mid + fLength/2 >= maxSequenceLen) end = maxSequenceLen;
					                                  else end = mid + length/2;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposCommaFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, comma\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, waffle\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "K]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end *= 1000;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#%fK]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end =  (size * 1000) + 1;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#" unsposFmtScanf "M]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end *= 1000 * 1000;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "#%fM]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, waffle with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end =  (size * 1000 * 1000) + 1;
					end += start-1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as two, caret\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "K]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   *= 1000;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^%fK]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with K\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   =  (size * 1000) + 1;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^" unsposFmtScanf "M]%n", &start, &end, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   *= 1000 * 1000;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "^%fM]%n", &start, &size, &charsUsed);
				if ((numItems == 2) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, caret with M\n");
					if ((start == 0) || (end == 0))
						{ action[len] = 0;  goto bad_limits; }
					end   =  (size * 1000 * 1000) + 1;
					start -= end / 2;
					end   += start-1;
					if (start < 1) start = 1;
					*endIsSoft = true;
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, unsposFmtScanf "..]%n", &start, &charsUsed);
				if ((numItems == 1) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, dots at end\n");
					if (start == 0)
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				charsUsed = -1;
				numItems = sscanf (action, ".." unsposFmtScanf "]%n", &end, &charsUsed);
				if ((numItems == 1) && (charsUsed == len+1))
					{
					//fprintf(stderr,"  parsed as one, dots at start\n");
					if (end == 0)
						{ action[len] = 0;  goto bad_limits; }
					parsed = true;
					}
				}

			if (!parsed)
				{
				action[len] = 0;		// (clear the sentinel)
				goto bad_action;
				}
			}

		if (pendingStart != 0) goto unsatisfied_start_2;

		if ((start != 0) && (end != 0) && (start > end))
			{
			temp = start;  start = end;  end = temp;
			*doRevCompFlags ^= rcf_revcomp;
			}

		if ((start != 0) || (end != 0))
			{
			*_start = start;
			*_end   = end;
			}
		}

	//////////
	// sanity checks
	//////////

	if ((*contigOfInterest != NULL) && (*namesFilename != NULL))
		suicidef ("(for %s) can't use these together:\n  %s\n  %s",
		          fname, *contigOfInterest, *namesFilename);

	if ((*contigOfInterest != NULL) && (*choresFilename != NULL))
		suicidef ("(for %s) can't use these together:\n  %s\n  %s",
		          fname, *contigOfInterest, *choresFilename);

	if ((*namesFilename != NULL) && (*choresFilename != NULL))
		suicidef ("(for %s) can't use these together:\n  %s\n  %s",
		          fname, *namesFilename, *choresFilename);

// (these are no longer prohibited combinations)
//
//	if ((*doPartitioning) && (*softMaskFilename != NULL))
//		{
//		if (*softMaskComplement)
//			suicidef ("(for %s) can't use [multi] with [softmask=keep:%s]",
//			          fname, *softMaskFilename);
//		else
//			suicidef ("(for %s) can't use [multi] with [softmask=%s]",
//			          fname, *softMaskFilename);
//		}
//
//	if ((*doPartitioning) && (*xMaskFilename != NULL))
//		{
//		if (*xMaskComplement)
//			suicidef ("(for %s) can't use [multi] with [xmask=keep:%s]",
//			          fname, *xMaskFilename);
//		else
//			suicidef ("(for %s) can't use [multi] with [xmask=%s]",
//			          fname, *xMaskFilename);
//		}
//
//	if ((*doPartitioning) && (*nMaskFilename != NULL))
//		{
//		if (*nMaskComplement)
//			suicidef ("(for %s) can't use [multi] with [nmask=keep:%s]",
//			          fname, *nMaskFilename);
//		else
//			suicidef ("(for %s) can't use [multi] with [nmask=%s]",
//			          fname, *nMaskFilename);
//		}

	return;

	//////////
	// failure exits
	//////////

empty_file_name:
	suicidef ("sequence file name is absent from \"%s\"", name);
	return; // (can't reach here)

empty_species_name:
	suicidef ("(for %s) empty nickname", parse+2);
	return; // (can't reach here)

bad_mask:
	suicidef ("(for %s) mask name needs closing } (%s)", fname, mask);
	return; // (can't reach here)

empty_mask_file_name:
	suicidef ("(for %s) use [unmask] instead of {}", fname);
	return; // (can't reach here)

actions_not_at_end:
	suicidef ("(for %s[%s)\n"
	          "The action list is not at the end of the sequence specifier.  See the README\n"
	          "section on sequence specifiers.  Perhaps you forgot a space after the closing\n"
	          "square bracket?",
	          fname,actions);
	return; // (can't reach here)

bad_action_list:
	suicidef ("(for %s) bad action list", fname);
	return; // (can't reach here)

blank_action:
	suicidef ("(for %s) blank action", fname);
	return; // (can't reach here)

bad_action:
	suicidef ("(for %s) bad action \"%s\"", fname, action);
	return; // (can't reach here)

bad_sequence_position:
	suicidef ("(for %s) bad limit \"%s\"", fname, action);
	return; // (can't reach here)

bad_limits:
	suicidef ("(for %s) bad limits \"%s\"", fname, action);
	return; // (can't reach here)

unsatisfied_start:
	suicidef ("(for %s) incomplete limits (%d,%s)", fname, pendingStart, action);
	return; // (can't reach here)

unsatisfied_start_2:
	suicidef ("(for %s) incomplete limits (%d)", fname, pendingStart);
	return; // (can't reach here)

many_separators:
	suicidef ("(for %s) only one separator allowed:\n"
	          "         separator=%c\n"
	          "         separator=%c",
	          fname, *separatorCh, actionName[0]);
	return; // (can't reach here)

bad_separator:
	if (actionName[0] == 0)
		suicidef ("(for %s) separator=<ch> requires a character", fname, actionName);
	else
		suicidef ("(for %s) bad separator, separator=%s, only one character allowed", fname);
	return; // (can't reach here)

many_name_files:
	suicidef ("(for %s) only one names file allowed:\n"
	          "         subset=%s\n"
	          "         subset=%s",
	          fname, *namesFilename, actionName);
	return; // (can't reach here)

bad_name_file:
	suicidef ("(for %s) subset= requires a names file", fname);
	return; // (can't reach here)

many_chore_files:
	suicidef ("(for %s) only one chores file allowed:\n"
	          "         chores=%s\n"
	          "         chores=%s",
	          fname, *choresFilename, actionName);
	return; // (can't reach here)

bad_chore_file:
	suicidef ("(for %s) chores= requires a chores file", fname);
	return; // (can't reach here)

many_subsamples:
	suicidef ("(for %s) only one subsampling allowed:\n"
	          "         subsample=%d/%d\n"
	          "         subsample=%s",
	          fname, *subsampleK, *subsampleN, actionName);
	return; // (can't reach here)

bad_subsample:
	suicidef ("(for %s) bad subsample \"%s\"", fname, actionName);
	return; // (can't reach here)

many_nicknames:
	suicidef ("(for %s) only one nickname allowed:\n"
	          "         nickname=%s\n"
	          "         nickname=%s",
	          fname, *nickname, actionName);
	return; // (can't reach here)

bad_nickname:
	suicidef ("(for %s) nickname= requires a non-empty string", fname);
	return; // (can't reach here)

many_name_parse_types:
	suicidef ("(for %s) only one name parsing allowed\n", fname);
	return; // (can't reach here)

many_name_triggers:
	suicidef ("(for %s) only one name trigger allowed:\n"
	          "         nameparse=tag:%s\n"
	          "         nameparse=tag:%s",
	          fname, *nameTrigger, actionName);
	return; // (can't reach here)

bad_name_trigger:
	suicidef ("(for %s) nameparse=tag: requires a non-empty string", fname);
	return; // (can't reach here)

many_soft_mask_files:
	suicidef ("(for %s) only one softmask allowed:\n"
	          "         softmask=%s\n"
	          "         softmask=%s",
	          fname, *softMaskFilename, actionName);
	return; // (can't reach here)

bad_soft_mask_file:
	suicidef ("(for %s) softMask= or softMask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_x_mask_files:
	suicidef ("(for %s) only one xmask allowed:\n"
	          "         xmask=%s\n"
	          "         xmask=%s",
	          fname, *xMaskFilename, actionName);
	return; // (can't reach here)

bad_x_mask_file:
	suicidef ("(for %s) xmask= or xmask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_n_mask_files:
	suicidef ("(for %s) only one nmask allowed:\n"
	          "         nmask=%s\n"
	          "         nmask=%s",
	          fname, *nMaskFilename, actionName);
	return; // (can't reach here)

bad_n_mask_file:
	suicidef ("(for %s) nmask= or nmask=keep: require a non-empty string", fname);
	return; // (can't reach here)

many_file_types:
	suicidef ("(for %s) more than one file type is defined", fname);
	return; // (can't reach here)

bad_file_format:
	suicidef ("(for %s) unknown file format: %s", fname, actionName);
	return; // (can't reach here)

many_quantums:
	suicidef ("(for %s) only one instance of quantum allowed", fname);
	return; // (can't reach here)

bad_code_file:
	suicidef ("(for %s) quantum= requires a non-empty string", fname);
	return; // (can't reach here)
	}

//----------
//
// detect_file_type--
//	Attempt to determine what type of file we are dealing with (e.g. fasta,
//	nib, quantum, etc.).
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence.
//
// Returns:
//	The type of file being read (one of seq_type_xxx);  failure causes program
//	fatality.
//
//----------

static int detect_file_type
   (seq*	_seq)
	{
	int		type = seq_type_unknown;
	char	buffer[maxSequenceHeader+3];
	u32		bufferLen;
	u32		magic;
	u8		ch;
	int		intCh;

	//////////
	// determine if it's a nib file (from the magic number)
	//////////

	// read the first four bytes;  if it's a recognizable magic number then it
	// must be the corresponding type of file;  otherwise we assume it must be
	// a fasta file (if not, it will die later)

	magic = read_4_little (_seq);

	if ((magic == nibMagicLittle)          || (magic == nibMagicBig))
		type = seq_type_nib;
	else if ((magic == twobitMagicLittle)  || (magic == twobitMagicBig))
		type = seq_type_2bit;
	else if ((magic == hsxMagicLittle)     || (magic == hsxMagicBig))
		type = seq_type_hsx;
	else if ((magic == qdnaMagicLittle)    || (magic == qdnaMagicBig))
		type = seq_type_qdna;
	else if ((magic == oldQdnaMagicLittle) || (magic == oldQdnaMagicBig))
		type = seq_type_qdna;

	// put those four bytes back in the file (in reverse of the read order)

	seq_ungetc (magic >> 24, _seq);
	seq_ungetc (magic >> 16, _seq);
	seq_ungetc (magic >>  8, _seq);
	seq_ungetc (magic      , _seq);

	if (type != seq_type_unknown)
		return type;

	//////////
	// determine if it's a fastq file
	//////////

	ch = seq_getc (_seq);
	seq_ungetc (ch, _seq);
	if (ch == '@')
		return seq_type_fastq;

	//////////
	// determine if it's a fasta or csfasta file;  these are very similar
	// formats;  if the first character is a '#' we know it is csfasta (because
	// regular fasta doesn't allow comments);  if the first character is a '>'
	// then it can still be either, in which case we need to skip the header
	// line and read the first two sequence characters
	//////////

	ch = seq_getc (_seq);
	if (ch == '#')
		{
		seq_ungetc (ch, _seq);
		return seq_type_csfasta;
		}

	if (ch != '>')
		seq_ungetc (ch, _seq);
	else
		{
		// read header

		bufferLen = 0;
		buffer[bufferLen++] = ch;

		do
			{
			intCh = seq_getc (_seq);
			if (intCh == EOF) goto unknown;
			if (bufferLen >= sizeof(buffer)-3) goto unknown;
			buffer[bufferLen++] = (char) intCh;
			} while (intCh != '\n');

		// first character must be a nucleotide or it is not fasta
		// $$$ we are ignoring the slim chance of an empty first sequence, or
		// $$$ .. that the first line contains just a single nucleotide

		intCh = seq_getc (_seq);
		if (intCh == EOF) goto unknown;
		buffer[bufferLen++] = intCh;
		if (ustrchr ("ACGTacgtNn", intCh) != NULL)
			{
			intCh = seq_getc (_seq);
			if (intCh == EOF) goto unknown;
			buffer[bufferLen++] = intCh;
			if (ustrchr ("ACGTacgtNn", intCh) != NULL)
				type = seq_type_fasta;
			else if (ustrchr ("0123", intCh) != NULL)
				type = seq_type_csfasta;
			}

	unknown:
		while (bufferLen > 0)
			seq_ungetc (buffer[--bufferLen], _seq);

		if (type != seq_type_unknown)
			return type;
		}

	//////////
	// if all else fails, assume it's a fasta file (for compatibility with
	// blastz)
	//////////

	if (type == seq_type_unknown)
		type = seq_type_fasta;

	return type;
	}

//----------
//
// read_4, read_4_big, read_4_little--
//	Read four bytes from a file, in big or little endian order.
// read_5, read_5_big, read_5_little--
//	Read five bytes from a file, in big or little endian order.
// read_6, read_6_big, read_6_little--
//	Read six bytes from a file, in big or little endian order.
//
//----------
//
// Arguments:
//	seq*	seq:			The sequence.
//	int		asBigEndian:	true  => read 'em as big endian
//							false => read 'em as little endian
//
// Returns:
//	The magic number read.
//
//----------

static u32 read_4
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_4_big    (_seq);
	            else return read_4_little (_seq);
	}

static u32 read_4_big
   (seq*	_seq)
	{
	u32		val;

	val =  seq_getc (_seq) << 24;
	val |= seq_getc (_seq) << 16;
	val |= seq_getc (_seq) << 8;
	val |= seq_getc (_seq);

	return val;
	}

static u32 read_4_little
   (seq*	_seq)
	{
	u32		val;

	val =  seq_getc (_seq);
	val |= seq_getc (_seq) << 8;
	val |= seq_getc (_seq) << 16;
	val |= seq_getc (_seq) << 24;

	return val;
	}


static u64 read_5
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_5_big    (_seq);
	            else return read_5_little (_seq);
	}

static u64 read_5_big
   (seq*	_seq)
	{
	u64		val;

	val =  ((u64) seq_getc (_seq)) << 32;
	val |=        seq_getc (_seq)  << 24;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq);

	return val;
	}

static u64 read_5_little
   (seq*	_seq)
	{
	u64		val;

	val =         seq_getc (_seq);
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 24;
	val |= ((u64) seq_getc (_seq)) << 32;

	return val;
	}


static u64 read_6
   (seq*	_seq,
	int		asBigEndian)
	{
	if (asBigEndian) return read_6_big    (_seq);
	            else return read_6_little (_seq);
	}

static u64 read_6_big
   (seq*	_seq)
	{
	u64		val;

	val =  ((u64) seq_getc (_seq)) << 40;
	val |= ((u64) seq_getc (_seq)) << 32;
	val |=        seq_getc (_seq)  << 24;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq);

	return val;
	}

static u64 read_6_little
   (seq*	_seq)
	{
	u64		val;

	val =         seq_getc (_seq);
	val |=        seq_getc (_seq)  << 8;
	val |=        seq_getc (_seq)  << 16;
	val |=        seq_getc (_seq)  << 24;
	val |= ((u64) seq_getc (_seq)) << 32;
	val |= ((u64) seq_getc (_seq)) << 40;

	return val;
	}

//----------
//
// skip_seq_whitespace--
//	Read characters from the associated sequence until we get something that
//	ain't whitespace.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to read.
//
// Returns:
//	(same as for getc())
//
//----------

static int skip_seq_whitespace
   (seq*	_seq)
	{
	int		ch;

	do
		{
		ch = seq_getc (_seq);
		} while ((ch == ' ') || (ch == '\t'));

	return ch;
	}

//----------
//
// seq_getc--
//	Read the next character from the associated file.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to read.
//
// Returns:
//	(same as for getc();  the character read, or EOF)
//
//----------

static int seq_getc
   (seq*	_seq)
	{
	int		ch;

	// if there are characters pending, get one from the pending buffer;
	// otherwise, feed one straight from the file

	if (_seq->pendingLen == 0)
		{
		ch = getc_or_die (_seq->f, _seq->filename);
		debugTextFile_3;
		return ch;
		}

	_seq->pendingLen--;
	ch = (int) (u8) *(_seq->pendingStack++);

	debugTextFile_4;
	return ch;
	}

//----------
//
// seq_ungetc--
//	Give back a character to the associated file.
//
// WARNING: This routine is not an exact drop in replacement for the standard
//	c routine ungetc().
//
//----------
//
// Arguments:
//	char	ch:		The character to return.  Characters should be returned in
//					.. the opposite order of that read (ie newest char returned
//					.. first).
//	seq*	seq:	The sequence being read.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void seq_ungetc
   (char	ch,
	seq*	_seq)
	{
	debugTextFile_5;

	if (_seq->pendingLen >= seqBufferSize)
		suicide ("seq_ungetc() buffer is already full");

	_seq->pendingLen++;
	*(--_seq->pendingStack) = ch;
	}

//----------
//
// skip_chars--
//	Skip the next so many characters from the associated file.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence being read.
//	u32		toSkip:	The number of characters to skip.
//
// Returns:
//	true if successful, false if there's a problem (such as reaching premature
//	end-of-file).
//
//----------

static int skip_chars
   (seq*	_seq,
	u32		toSkip)
	{
	int		ch;

	// see if we have any to skip in the pending buffer

	if (_seq->pendingLen >= toSkip)
		{
		// we have all we need in the pending buffer

		_seq->pendingLen   -= toSkip;
		_seq->pendingStack += toSkip;
		return true;
		}

	if (_seq->pendingLen > 0)
		{
		// everything in the pending buffer will be skipped

		toSkip            -= _seq->pendingLen;
		_seq->pendingLen   =  0;
		_seq->pendingStack = _seq->pendingChars + seqBufferSize;
		}

	if (toSkip == 0) return true; // (none left to skip)

	// skip the rest by seeking past the characters

	if (fseek (_seq->f, toSkip, SEEK_CUR) != 0)
		{
		// seek failed, so let's try reading instead

		while (toSkip-- > 0)
			{
			ch = getc (_seq->f);
			if ((ch == EOF) || (ferror (_seq->f)))
				return false;
			}
		}

	return true;
	}

//----------
//
// test_rewindability--
//	Test whether a sequence's underflying file is rewinable.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence being read.
//
// Returns:
//	An fseek error code;  zero indicates the underlying file is rewindable;
//	any other value indicates that it is not.
//
//----------

static int test_rewindability
   (seq*		_seq)
	{
	long int	savedFilePos;

	savedFilePos = ftell (_seq->f);
	return fseek (_seq->f, savedFilePos, SEEK_SET);
	}

//----------
//
// save_fstate, restore_fstate--
//	Save and restore the state of the associated file.  This lets the caller
//	read ahead in a file, then return to the original point.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence being read.
//
// Returns:
//  (nothing;  failure causes program fatality)
//
//----------

static void save_fstate
   (seq*	_seq)
	{
	if (_seq == NULL) suicide ("save_fstate(NULL)");

	// save read head

	debugTextFile_6;

	_seq->savedFilePos  = ftell (_seq->f) - _seq->pendingLen;
	_seq->hasSavedState = true;
	}

static void restore_fstate
   (seq*	_seq)
	{
	int		err;

	if (_seq == NULL)         suicide ("restore_fstate(NULL)");
	if (!_seq->hasSavedState) suicide ("restore_fstate(), no state saved");

	// restore read head

	debugTextFile_7;

	err = fseek (_seq->f, _seq->savedFilePos, SEEK_SET);
	if (err != 0)
		suicidef_with_perror ("restore_fstate(), fseek returned %d", err);

	// restore pending character state

	_seq->pendingLen   = 0;
	_seq->pendingStack = _seq->pendingChars + seqBufferSize;
	}

//----------
//
// match_composition--
//	Count the number of matched DNA letter pairs in a gap free alignment.
//
//----------
//
// Arguments:
//	seq*	seq1:		The first sequence.
//	unspos 	pos1:		The subsequence start position in seq1 (origin-0).
//	seq*	seq2:		The second sequence.
//	unspos 	pos2:		The subsequence start position in seq2 (origin-0).
//	unspos 	length:		The length of the subsequence.
//	unspos	count[4][4]:Place to return the counts of each matched DNA letter
//						.. pair.  Indexing is as per nuc_to_bits.
//
// Returns:
//	nothing;  composition is returned in the count[][] array
//
//----------
//
// Note:  Masked (lowercase) bp do not contribute to the results.
//
//----------

void match_composition
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	unspos	count[4][4])
	{
	u8*		s1 = seq1->v + pos1;
	u8*		s2 = seq2->v + pos2;
	unspos	ix;
	int		r, c;

	for (r=0 ; r<4 ; r++)
		for (c=0 ; c<4 ; c++)
			count[r][c] = 0;

	for (ix=0 ; ix<length ; ix++)
		{
		r = upper_nuc_to_bits[*(s1++)];
		c = upper_nuc_to_bits[*(s2++)];
		if ((r >= 0) && (c >= 0)) 
			count[r][c]++;
		}
	}

//----------
//
// percent_identical--
//	Determine the percentage of bases that match in two subsequences.
//
//----------
//
// Arguments:
//	seq*	seq1:	The first sequence.
//	unspos 	pos1:	The subsequence start position in seq1 (origin-0).
//	seq*	seq2:	The second sequence.
//	unspos 	pos2:	The subsequence start position in seq2 (origin-0).
//	unspos 	length:	The length of the subsequence.
//
// Returns:
//	The percentage (an integer in the range 0..100).
//
//----------
//
// Note:  Masked (lowercase) bp *do* contribute to the results, but illegal
//        values like N do not (and they are not counted in the denominator
//        either).
//
//----------

int percent_identical
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	u8*		s1 = seq1->v + pos1;
	u8*		s2 = seq2->v + pos2;
	s8		c1, c2;
	unspos	numMatches = 0;
	unspos	denom = 0;
	unspos	ix;

	if (length == 0)
		return 0;

	if ((seq1->fileType == seq_type_qdna)
	 || (seq2->fileType == seq_type_qdna))
		return 0;

	for (ix=0 ; ix<length ; ix++)
		{
		c1 = nuc_to_bits[*(s1++)];
		c2 = nuc_to_bits[*(s2++)];
		if ((c1 >= 0) && (c2 >= 0))
			{
			if (c1 == c2) numMatches++;
			denom++;
			}
		}

	if (denom == 0)
		return 0;
	else
		return (200*numMatches + denom) / (2*denom); // 100*numMatches/denom, rounded
	}

//----------
//
// score_match--
//	Determine the substitution score of aligned bases in two subsequences.
//
//----------
//
// Arguments:
//	scoreset*	scoring:	The scoring scheme to use.
//	seq*		seq1:		The first sequence.
//	unspos		pos1:		The subsequence start position in seq1 (origin-0).
//	seq*		seq2:		The second sequence.
//	unspos		pos2:		The subsequence start position in seq2 (origin-0).
//	unspos		length:		The length of the subsequence.  Note that this may
//							.. be zero.
//
// Returns:
//	The substitution score of the two subsequences.
//
//----------

score score_match
   (scoreset*	scoring,
	seq*		seq1,
	unspos		pos1,
	seq*		seq2,
	unspos		pos2,
	unspos		length)
	{
	u8*			s1 = seq1->v + pos1;
	u8*			s2 = seq2->v + pos2;
	u8*			stop = s1 + length;
	score		similarity = 0;

	if (length == 0)
		return (score) 0;

	while (s1 < stop)
		similarity += scoring->sub[*(s1++)][*(s2++)];

	return similarity;
	}

//----------
//
// dump_aligned_nucleotides--
//	Dump the nucleotides (from each sequence) for a gap-free alignment.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The first aligned position in sequence 1.
//	seq*	seq2:	The other sequence.
//	unspos	pos2:	The first aligned position in sequence 2.
//	unspos	length:	The length of the alignment.
//
// Returns:
//	(nothing)
//
//----------

void dump_aligned_nucleotides
   (FILE*	f,
	seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	int		isRev1 = ((seq1->revCompFlags & rcf_rev) != 0);
	int		isRev2 = ((seq2->revCompFlags & rcf_rev) != 0);
	char*	start1 = (char*) seq1->v + pos1;
	char*	start2 = (char*) seq2->v + pos2;
	int		digits = 10;

	fprintf      (f, unsposStarFmt "%c:", digits, pos1+1, (isRev1)?'-':'+');
	print_prefix (f, (char*) seq1->v + pos1, (int) length);
	fprintf      (f, "\n");

	fprintf      (f, "%*s  ", digits, "");
	print_dna_similarities
	             (f, start1, start2, (int) length);
	fprintf      (f, "\n");

	fprintf      (f, unsposStarFmt "%c:", digits, pos2+1, (isRev2)?'-':'+');
	print_prefix (f, (char*) seq2->v + pos2, (int) length);
	fprintf      (f, "\n");
	}

//----------
//
// dump_sequence--
//	Write a sequence to a file (for debugging).
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq:	The sequence to print.
//
// Returns:
//	(nothing)
//
//----------

void dump_sequence
   (FILE*	f,
	seq*	_seq)
	{
	char	buffer[101];
	unspos	ix, start = 0;
	int		width;
	char	ch;
	int		bx;
	int		needSeparator;

	start = _seq->len;
	width = 1;
	while (start > 9) { start/=10;  width++; }

	needSeparator = false;

	bx = 0;
	for (ix=0 ; ix<_seq->len ; ix++)
		{
		ch = (char) _seq->v[ix];

		if (ch == 0)
			{
			if (bx > 0)
				{
				if (needSeparator)
					{ fprintf (f, "%*s  =====\n", width, "");  needSeparator = false; }
				buffer[bx] = 0;
				fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
				}
			bx = 0;
			needSeparator = true;
			continue;
			}

		if (bx == sizeof(buffer)-1)
			{
			if (needSeparator)
				{ fprintf (f, "%*s  =====\n", width, "");  needSeparator = false; }
			buffer[bx] = 0;
			fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
			bx = 0;
			}

		if (bx == 0) start = ix;
		buffer[bx++] = ch;
		}

	if (bx > 0)
		{
		if (needSeparator)
			{ fprintf (f, "%*s  =====\n", width, "");  needSeparator = false; }
		buffer[bx] = 0;
		fprintf (f, unsposStarFmt ": %s\n", width, start, buffer);
		}
	}

//----------
//
// dump_sequence_state--
//	Write a sequence's state information to a file (for debugging).
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq:	The sequence to print.
//
// Returns:
//	(nothing)
//
//----------

void dump_sequence_state
   (FILE*			f,
	seq*			_seq)
	{
	seqpartition*	sp = &_seq->partition;
	u32				ix;

	fprintf (f, "size:               %s\n", commatize(_seq->size));
	fprintf (f, "len:                %s\n", commatize(_seq->len));
	fprintf (f, "needsVq:            %d\n", _seq->needsVq);

	fprintf (f, "v:                  %s%p\n", (_seq->vOwner)?"[owner] ":"", _seq->v);
	if (_seq->vc != NULL) fprintf (f, "vc:                 %s%p\n", (_seq->vcOwner)?"[owner] ":"", _seq->vc);
	if (_seq->vq != NULL) fprintf (f, "vq:                 %s%p\n", (_seq->vqOwner)?"[owner] ":"", _seq->vq);

	fprintf (f, "startLoc:           %s\n", commatize(_seq->startLoc));
	fprintf (f, "trueLen:            %s%s\n", (_seq->needTrueLen)?"[need] ":"", commatize(_seq->trueLen));
	fprintf (f, "revCompFlags:       %d\n", _seq->revCompFlags);

	if (_seq->contigOfInterest != NULL) fprintf (f, "vc:                 \"%s\"\n", _seq->contigOfInterest);
	fprintf (f, "contig:             %u\n", _seq->contig);
	fprintf (f, "preLoaded:          %d\n", _seq->preLoaded);

	fprintf (f, "lockedHeader:       %d\n", _seq->lockedHeader);
	fprintf (f, "headerSize:         %s\n", commatize(_seq->headerSize));
	if (_seq->header != NULL) fprintf (f, "header:             %s\"%s\"\n", (_seq->headerOwner)?"[owner] ":"", _seq->header);
	                     else fprintf (f, "header:             %s(null)\n", (_seq->headerOwner)?"[owner] ":"");
	fprintf (f, "shortHeaderSize:    %s\n", commatize(_seq->shortHeaderSize));
	if (_seq->shortHeader != NULL) fprintf (f, "shortHeader:        %s\"%s\"\n", (_seq->shortHeaderOwner)?"[owner] ":"", _seq->shortHeader);
	                          else fprintf (f, "shortHeader:        %s(null)\n", (_seq->shortHeaderOwner)?"[owner] ":"");
	fprintf (f, "hasNickname:        %d\n", _seq->hasNickname);
	fprintf (f, "trueHeaderSize:     %s\n", commatize(_seq->trueHeaderSize));
	if (_seq->trueHeader != NULL) fprintf (f, "trueHeader:         %s\"%s\"\n", (_seq->trueHeaderOwner)?"[owner] ":"", _seq->trueHeader);
	                         else fprintf (f, "trueHeader:         %s(null)\n", (_seq->trueHeaderOwner)?"[owner] ":"");

	if (_seq->hasLeftFence)  fprintf (f, "leftFence:          " unsposFmt " %02X\n", _seq->leftFencePos,  _seq->leftFenceCh);
	if (_seq->hasRightFence) fprintf (f, "rightFence:         " unsposFmt " %02X\n", _seq->rightFencePos, _seq->rightFenceCh);

	if (_seq->filename != NULL) fprintf (f, "filename:           \"%s\"\n", _seq->filename);
	                       else fprintf (f, "filename:           (null)\n");
	if (_seq->f        != NULL) fprintf (f, "f:                  %p\n", _seq->f);
	                       else fprintf (f, "f:                  (null)\n");
	fprintf (f, "fileType:           %s (%d)\n", seqTypeNames[_seq->fileType], _seq->fileType);
	fprintf (f, "rewindable:         %d\n", _seq->rewindable);

	fprintf (f, "pending:            ");
	if (_seq->pendingLen == 0)
		fprintf (f, "(empty)\n");
	else
		{
		fprintf (f, "\"");
		for (ix=0 ; ix<_seq->pendingLen ; ix++)
			fprintf (f, "%c", _seq->pendingStack[ix]);
		fprintf (f, "\"\n");
		}

	if (_seq->hasSavedState) fprintf (f, "savedFilePos:       %016lX\n", _seq->savedFilePos);

	if (_seq->namesFile != NULL)
		{
		if (_seq->namesFilename != NULL) fprintf (f, "namesFilename:      \"%s\"\n", _seq->namesFilename);
		                            else fprintf (f, "namesFilename:      (null)\n");
		fprintf (f, "namesFile:          %p\n", _seq->namesFile);
		fprintf (f, "nextContigName:     \"%s\"\n", _seq->nextContigName);
		fprintf (f, "contigPending:      %d\n",     _seq->contigPending);
		}

	if (_seq->choresFile != NULL)
		{
		if (_seq->choresFilename != NULL) fprintf (f, "choresFilename:     \"%s\"\n", _seq->choresFilename);
		                             else fprintf (f, "choresFilename:     (null)\n");
		fprintf (f, "choresFile:         %p\n", _seq->choresFile);
		fprintf (f, "choresLineNum:      %d\n", _seq->choresLineNum);
		fprintf (f, "chore.num:          %d\n",     _seq->chore.num);
		fprintf (f, "chore.tName:        \"%s\"\n", _seq->chore.tName);
		if (_seq->chore.tSubrange) fprintf (f, "chore.tSubrange:    " unsposDotsFmt "\n", _seq->chore.tStart, _seq->chore.tEnd);
		                      else fprintf (f, "chore.tSubrange:    (whole sequence)\n");
		if (_seq->chore.qSubrange) fprintf (f, "chore.qSubrange:    " unsposDotsFmt "\n", _seq->chore.qStart, _seq->chore.qEnd);
		                      else fprintf (f, "chore.qSubrange:    (whole sequence)\n");
		fprintf (f, "chore.qStrand:      %s\n", (_seq->chore.qStrand<0)?"- strand":((_seq->chore.qStrand>0)?"+ strand":"both strands"));
		fprintf (f, "chore.idTag:        \"%s\"\n", _seq->chore.idTag);

		fprintf (f, "chore.targetInt:    " unsposDotsFmt "\n", _seq->chore.targetInterval.s, _seq->chore.targetInterval.e);
		fprintf (f, "chore.queryInt:     " unsposDotsFmt "\n", _seq->chore.queryInterval.s,  _seq->chore.queryInterval.e);
		}

	if (_seq->subsampleN > 0) fprintf (f, "subsample:          %d/%d [skip %d]\n", _seq->subsampleK, _seq->subsampleN, _seq->subsampleN);

	if (_seq->softMaskFilename != NULL) fprintf (f, "softMaskFilename:   %s\"%s\"\n", (_seq->softMaskComplement)?"[keep] ":"", _seq->softMaskFilename);
	if (_seq->xMaskFilename    != NULL) fprintf (f, "xMaskFilename:      %s\"%s\"\n", (_seq->xMaskComplement)?"[keep] ":"",    _seq->xMaskFilename);
	if (_seq->nMaskFilename    != NULL) fprintf (f, "nMaskFilename:      %s\"%s\"\n", (_seq->nMaskComplement)?"[keep] ":"",    _seq->nMaskFilename);

	if ((_seq->startLimit != 0) || (_seq->endLimit != 0))
		fprintf (f, "limits:             " unsposDotsFmt "%s\n", _seq->startLimit, _seq->endLimit, (_seq->endIsSoft)?"[soft] ":"");

	fprintf (f, "doRevCompFlags:     %d\n", _seq->doRevCompFlags);
	fprintf (f, "doUnmask:           %d\n", _seq->doUnmask);
	fprintf (f, "doPartitioning:     %d\n", _seq->doPartitioning);
	fprintf (f, "doJoin:             %d\n", _seq->doJoin);
	if (_seq->separatorCh != 0)
		{
		if ((' ' < _seq->separatorCh) && (_seq->separatorCh <= '~')) fprintf (f, "separatorCh:        '%c' %02X\n", _seq->separatorCh, _seq->separatorCh);
		                                                        else fprintf (f, "separatorCh:        %02X\n",      _seq->separatorCh);
		}
	fprintf (f, "useFullNames:       %d\n", _seq->useFullNames);
	fprintf (f, "nameParseType:      %d\n", _seq->nameParseType);
	if (_seq->nameTrigger != NULL) fprintf (f, "nameTrigger:        \"%s\"\n", _seq->nameTrigger);
	fprintf (f, "allowAmbiDNA:       %d\n", _seq->allowAmbiDNA);
	if (_seq->qToComplement != NULL) fprintf (f, "qToComplement:      %p\n", _seq->qToComplement);
	if (_seq->qCoding       != NULL) fprintf (f, "qCoding:            %p\n", _seq->qCoding);

	if (_seq->fileType == seq_type_2bit)
		{
		; // $$$ dump _seq->twoBit
		}

	else if (_seq->fileType == seq_type_hsx)
		{
		; // $$$ dump _seq->hsx
		}

	if (sp->p != NULL)
		{
		fprintf (f, "partition.state:    %d\n", sp->state);
		fprintf (f, "partition.size:     %s\n", commatize(sp->size));
		fprintf (f, "partition.len:      %s\n", commatize(sp->len));
		fprintf (f, "partition.p:        %p\n", sp->p);

		fprintf (f, "partition.poolSize: %s\n", commatize(sp->poolSize));
		fprintf (f, "partition.poolLen:  %s\n", commatize(sp->poolLen));
		fprintf (f, "partition.pool:     %s%p\n", (sp->poolOwner)?"[owner] ":"", sp->pool);

		print_partition_table (f, _seq);
		}

	}

//----------
//
// sequence_zero_stats--
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

void sequence_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&sequenceStats, 0, sizeof(sequenceStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// sequence_show_stats,
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

void sequence_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;

	fprintf (f, " partition lookups: %s\n", commatize (sequenceStats.partitionLookups));
	if (sequenceStats.partitionHits != 0)
		fprintf (f, "    partition hits: %s\n", commatize (sequenceStats.partitionHits));
	if (sequenceStats.partitionLookups != 0)
		fprintf (f, " lookup iterations: %s (%.1f per)\n",
		            commatize (sequenceStats.lookupIterations),
		            sequenceStats.lookupIterations / ((double) sequenceStats.partitionLookups));
	fprintf (f, "-------------------\n");

#endif // collect_stats
	}
