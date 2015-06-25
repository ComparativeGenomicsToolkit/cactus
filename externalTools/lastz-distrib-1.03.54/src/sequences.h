//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sequences.h
//
//----------

#ifndef sequences_H				// (prevent multiple inclusion)
#define sequences_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff

// establish ownership of global variables

#ifdef sequences_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef sequences_owner
int sequences_keepFastaArrow = false;	// true => keep ">" on header in
										// .. load_fasta_sequence()
int sequences_dbgDumpSequence = false;	// true => dump sequences to stderr
										// .. after input, masking, etc.
int sequences_dbgAllowColors = false;	// true => allow color space
#else
global int sequences_keepFastaArrow;
global int sequences_dbgDumpSequence;
global int sequences_dbgAllowColors;
#endif

//----------
//
// sequence sizing types
//	(see also "malloc sizing range" in utilties.h)
//	These types control the range of sequence lengths we can handle before some
//	values will overflow.  For the sake of efficiency, we generally do not check
//	for overflow when indexing sequences.
//
//	Sequence lengths are normally assumed to be small enough to fit into a
//	31-bit integer.  This gives a maximum length of about 2.1 billion bp, which
//	is half the length of a (hypothetical) monoploid human genome.  The
//	programmer can override this at compile time by defining max_sequence_index
//	as 32 or 63.
//
//	maxSequenceLen (L) is the longest sequence that we intend to support.
//
//	unspos is an unsigned type that can represent 0..L.
//
//	sgnpos is a signed type that can represent -L..L.
//
//	possum is an unsigned type that can represent a sum of values in the 0..L
//	range.
//
//----------

#define maxSequenceOverrun 10

#if defined(max_sequence_index)
#define maxSequenceIndex max_sequence_index
#else
#define maxSequenceIndex 31
#endif

#if (maxSequenceIndex == 31)
#define unspos_sz 4
typedef u32 unspos;
typedef s32 sgnpos;
typedef u64 possum;
#elif (maxSequenceIndex == 32)
#define unspos_sz 4
typedef u32 unspos;
typedef s64 sgnpos;
typedef u64 possum;
#elif (maxSequenceIndex == 63)
#define unspos_sz 8
typedef u64 unspos;
typedef s64 sgnpos;
typedef u64 possum;
#else
#error ***** undecipherable max sequence length definition *****
#endif

#define maxSequenceLen ((unspos)(((1ULL)<<maxSequenceIndex)-(maxSequenceOverrun+1)))
#define seqposInfinity ((unspos) -1LL)

#ifdef override_inttypes

#if (maxSequenceIndex == 31)
#define unsposFmt     "%u"
#define unsposStarFmt "%*u"
#define sgnposFmt     "%d"
#define possumFmt     "%ju"
#elif (maxSequenceIndex == 32)
#define unsposFmt     "%u"
#define unsposStarFmt "%*u"
#define sgnposFmt     "%jd"
#define possumFmt     "%ju"
#elif (maxSequenceIndex == 63)
#define unsposFmt     "%ju"
#define unsposStarFmt "%*ju"
#define sgnposFmt     "%jd"
#define possumFmt     "%ju"
#endif

#define unsposFmtScanf     unsposFmt
#define unsposStarFmtScanf unsposStarFmt
#define sgnposFmtScanf     sgnposFmt
#define possumFmtScanf     possumFmt

#else

#include <inttypes.h>
#if (maxSequenceIndex == 31)
#define unsposFmt      "%"  PRIu32
#define unsposStarFmt  "%*" PRIu32
#define sgnposFmt      "%"  PRId32
#define possumFmt      "%"  PRIu64
#define unsposFmtScanf "%"  SCNu32
#define sgnposFmtScanf "%"  SCNd32
#define possumFmtScanf "%"  SCNu64
#elif (maxSequenceIndex == 32)
#define unsposFmt      "%"  PRIu32
#define unsposStarFmt  "%*" PRIu32
#define sgnposFmt      "%"  PRId64
#define possumFmt      "%"  PRIu64
#define unsposFmtScanf "%"  SCNu32
#define sgnposFmtScanf "%"  SCNd64
#define possumFmtScanf "%"  SCNu64
#elif (maxSequenceIndex == 63)
#define unsposFmt      "%"  PRIu64
#define unsposStarFmt  "%*" PRIu64
#define sgnposFmt      "%"  PRId64
#define possumFmt      "%"  PRIu64
#define unsposFmtScanf "%"  SCNu64
#define sgnposFmtScanf "%"  SCNd64
#define possumFmtScanf "%"  SCNu64
#endif

#endif // override_inttypes

#define unsposSlashFmt       unsposFmt      "/"   unsposFmt
#define unsposSlashSFmt      unsposFmt      "%s/" unsposFmt      "%s"
#define unsposSlashCFmt      unsposFmt      "%c/" unsposFmt      "%c"
#define unsposCommaFmt       unsposFmt      ","   unsposFmt
#define unsposDashFmt        unsposFmt      "-"   unsposFmt
#define unsposDotsFmt        unsposFmt      ".."  unsposFmt

#define unsposSlashFmtScanf  unsposFmtScanf "/"   unsposFmtScanf
#define unsposSlashSFmtScanf unsposFmtScanf "%s/" unsposFmtScanf "%s"
#define unsposSlashCFmtScanf unsposFmtScanf "%c/" unsposFmtScanf "%c"
#define unsposCommaFmtScanf  unsposFmtScanf ","   unsposFmtScanf
#define unsposDashFmtScanf   unsposFmtScanf "-"   unsposFmtScanf
#define unsposDotsFmtScanf   unsposFmtScanf ".."  unsposFmtScanf

// interval (start-end position pair)
// note that we do not define here whether the the interval is origin-zero/one
// .. or open/closed

typedef struct interval { unspos s;  unspos e; } interval;

//----------
//
// sequence data structures and types
//
//----------
//
// sequences--
//	A sequence is a vector of characters (unsigned 8-bit values).
//
//----------
//
// Notes:
//	(1)	In the struct definitions, an (H) indicates pointers to items allocated
//		separately, and which therefore must be deallocated when the sequence is
//		deallocated.  (O) indicates pointers to items for which another field
//		tells us whether or not should be deallocated here.  (I) indicates
//		pointers to items within the same structure.  (X) indicates pointers to
//		other things but which do not need to be deallocated.  (F) indicates
//		files which must be closed when the sequence is deallocated.
//
//	(2)	For partitioned sequences, p[0].sepBefore is always 0 and
//		p[len].sepBefore is always seq->len.  seq->v[p[ix].sepBefore] is a NUL,
//		for all ix (even for 0 and len).  All partitions are bounded on both
//		sides by a NUL, even the first partition.
//
//----------

#define maxSequenceName     100
#define maxSequenceHeader   992
#define seqBufferSize       (maxSequenceHeader+32)
#define maxFastqSequenceLen 10000
#define maxChoreTagLen      15

typedef struct chore
	{
	int		num;				// the number of this chore among chores on the
								// .. same query;  the first chore is 1

	char	tName			 	// name of the target sequence;  an empty
			  [maxSequenceName+1];// .. string indicates a wildcard name;
			  					// .. note that the query sequence name is
			  					// .. stored in seq.nextContigName

	int		tSubrange;			// true  => tStart,tEnd are meaningful
								// false => the entire target is to be used 
	unspos	tStart, tEnd;		// origin-one half-open

	int		qSubrange;			// true  => qStart,qEnd are meaningful
								// false => the entire query is to be used 
	unspos	qStart, qEnd;		// origin-one half-open
	int		qStrand;			// 0   => search + strand only
								// < 0 => search - strand only
								// > 0 => search both strands

	char	idTag				// user-specified tag to identify this chore
			  [maxChoreTagLen+1];//.. an empty string indicates "no tag"

	interval targetInterval;	// corresponding index range in target->v[] or
	interval queryInterval;		// and query->v[];  the intervals are origin-
								// .. zero closed
	} chore;


typedef struct partition
	{							//  .. layout must match struct cappartition (capsule.h)
	unspos	sepBefore;			//     the position (in seq->v) of the
								//     .. separating NUL preceding this
								//     .. partition;  valid only for entries
								//     .. 0..len in seqpartition->p;  for
								//     .. entry seqpartition.len, it is the
								//     .. position of the final NUL in seq->v;
								//     .. see note (2) above
	unspos	sepAfter;			//     the position (in seq->v) of the
								//     .. separating NUL following this
								//     .. partition;  valid only for entries
								//     .. 0..len-1 in seqpartition->p
	u32		contig;				//     the contig number of this partition
								//     .. in the actual sequence file;  valid
								//     .. only for entries 0..len-1 in
								//     .. seqpartition->p
	unspos	startLoc;			//     1-based starting location;  this is only
								//     .. valid after a sequence is loaded
	unspos	trueLen;			//     the number of characters in the actual
								//     .. sequence for this partition (similar
								//     .. to seq->trueLen);  valid only for
								//     .. entries 0..len-1 in seqpartition->p
	u32		header;				//     the header of this partition;  valid
								//     .. only for entries 0..len-1 in
								//     .. seqpartition->p;  points into
								//     .. seqpartition->pool
	} partition;


typedef struct seqpartition
	{
	// information for a partitioned sequence-- a single sequence which holds
	// several actual sequences, each separated by a NUL character (a zero)
	
	int		state;				//     the current initialization state;  one
								//     .. of seqpart_xxx below

	u32		size;				//     the number of entries allocated for
								//     .. p[];  note that this must always be
								//     .. at least 1 more than len
	u32		len;				//     the number of partitions of the sequence

	partition*	p;				// (H) an array, indexed by 0..len, of the
								//     .. partitions

	u32		poolSize;			//     the number of bytes allocated (size) and
	u32		poolLen;			//     .. actually used (len) for pool[]
	char*	pool;				// (O) an array containing the characters for
								//     .. header[]
	int		poolOwner;			//     true  => we must deallocate pool[]
								//     false => we must *not* deallocate it
	} seqpartition;

enum							// states for a partitioned sequence (only valid
	{							// .. if seqpartition.p != NULL
	seqpart_empty = 0,			// created but no partitions loaded yet
	seqpart_reusable,			// all partitions loaded, ready for use if same
								// .. sequence is to be used again
	seqpart_loading,			// some partitions loaded, but not all
	seqpart_ready				// all partitions loaded
	};


typedef struct twobit
	{
	// information for 2bit files only

	int		bigEndian;			//     true  => file was written as big endian
								//     false => written as little endian
	u32		numContigs;			//     the number of contigs in the file (this
								//     .. is the sequenceCount field)
	long int indexFilePos;		//     position (as per ftell) of the beginning
								//     .. of the file index
	long int contigFilePos;		//     position (as per ftell) of the next entry
								//     .. in the file index

	int		contigLoaded;		//     true => the sequence corresponding to
								//             .. seq.contig has been loaded
								//             .. into v[]

	u32		nBlocksSize;		//     "temporary" arrays to hold n-blocks;
	u32*	nBlockStarts;		// (H) .. nBlocksSize is the number of entries
	u32*	nBlockSizes;		// (H) .. *allocated* for nBlockStarts[] and
								//     .. nBlockSizes[]

	u32		mBlocksSize;		//     "temporary" arrays to hold mask-blocks;
	u32*	mBlockstarts;		// (H) .. mBlocksSize is the number of entries
	u32*	mBlocksizes;		// (H) .. *allocated* for mBlockstarts[] and
								//     .. mBlocksizes[]
	} twobit;


typedef struct hsxfileinfo		// (see hsx.fileInfo)
	{
	char*	name;				//     file's name
	FILE*	f;					//     (if non-NULL), pointer to open file
	} hsxfileinfo;

typedef struct hsx
	{
	// information for hsx files only

	int		bigEndian;			//     true  => file was written as big endian
								//     false => written as little endian

	u32		version;
	u32		numContigs;			//     the number of contigs in the file (this
								//     .. is the numSequences field)
	u32		numFiles;			//     the number of files in the file index
	u32		numBuckets;			//     the number of hash buckets (not including
								//     .. the sentinel bucket)

	u64		fileTableOffset;	//     position of the beginning of the file
								//     .. table
	u64		hashTableOffset;	//     position of the beginning of the hash
								//     .. table
	u64		seqTableOffset;		//     position of the beginning of the sequence
								//     .. index table

	int		contigLoaded;		//     true => the sequence corresponding to
								//             .. seq.contig has been loaded
								//             .. into v[]

	u64		contigFilePos;		//     position of the next entry in the
								//     .. sequence table

	hsxfileinfo* fileInfo;		// (H) internal copy of the file info table;
								//     .. indexed by 0..numFiles-1;
								//     .. *seqFileInfo[i].name is the
								//     .. reconstructed name of file i;
								//     .. *seqFileInfo[i].f is the file pointer
								//     .. to file i, if non-NULL;
	int     seqFileIx;			//     index (into fileInfo) of the current
								//     .. sequence
	u64		seqFilePos;			//     position of the current sequence, in
								//     .. seqFile
	u64		seqLength;			//     length (in bytes) of the current sequence
	} hsx;


typedef struct seq
	{
	// sequence content

	unspos	size;				//     the number of bytes allocated for v[]
	unspos	len;				//     the number of characters in the sequence,
								//     .. not including a terminating zero
	int		needsVq;			//     true  => we allocate for vq with v
								//     fasle => we don't
	u8*		v;					// (O) the sequence content;  if v is not NULL,
								//     .. v[len] is a terminating zero
	u8*		vc;					// (O) the sequence content in color space;
								//     .. if vc is not NULL, vc[len] is a
								//     .. terminating zero
	u8*		vq;					// (O) the sequence's base-call qualities (ascii)
								//     .. if vq is not NULL, vq[len] is a
								//     .. terminating zero

	int		vOwner;				//     true  => we must deallocate v[]
								//     false => we must *not* deallocate it
	int		vcOwner;			//     true  => we must deallocate vc[]
								//     false => we must *not* deallocate it
	int		vqOwner;			//     true  => we must deallocate vq[]
								//     false => we must *not* deallocate it

	unspos	startLoc;			//     1-based starting location;  this is only
								//     .. valid after a sequence is loaded
	unspos	trueLen;			//     the number of characters in the actual
								//     .. sequence, including those not stored
								//     .. in v[]
	int		needTrueLen;		//     true  => trueLen must be set correctly,
								//              .. even if this means reading
								//              .. additional characters outside
								//              .. the desired (sub)interval
	int		revCompFlags;		//     two bits describing how this sequence
								//     .. relates to what was read from the
								//     .. file; the four values are the rcf_xxx
								//     .. values defined below
	char*	contigOfInterest;	// (H) the name of the only sequence of
								//     .. interest;  NULL means we're interested
								//     .. in every sequence
	u32		contig;				//     the number of the subsequence of the
								//     .. actual file (provided for fasta, 2bit
								//     .. and hsx files);  the first contig is 1
	int		preLoaded;			//     true => the data currently in v[] has
								//             .. been pre-loaded from the file
								//             .. behind the caller's back

	int		lockedHeader;		//     true => don't change header or
								//             .. shortHeader
	u32		headerSize;			//     number of bytes allocated for header[]
	char*	header;				// (O) the sequence's header (e.g. for each
								//     .. sequence in a fasta file);  this may
								//     .. be NULL
	int		headerOwner;		//     true  => we must deallocate header[]
								//     false => we must *not* deallocate it
	u32		shortHeaderSize;	//     number of bytes allocated for
								//     .. shortHeader[]
	char*	shortHeader;		// (O) short version of the header;  this may be
								//     .. NULL
	int		shortHeaderOwner;	//     true  => we must deallocate shortHeader[]
								//     false => we must *not* deallocate it
	int		hasNickname;		//     true => header is a nickname, and should
								//             .. be copied into shortHeader
								//             .. without stripping paths
	u32		trueHeaderSize;		//     number of bytes allocated for trueHeader[]
	char*	trueHeader;			// (O) the sequence's true header (primarily
								//     .. used to validate fastq files);  this
								//     .. may be NULL
	int		trueHeaderOwner;	//     true  => we must deallocate fauxHeader[]
								//     false => we must *not* deallocate it

	int		hasLeftFence;		//     a 'fence' at the left end of some
	unspos	leftFencePos;		//     .. interval;  if hasLeftFence is true,
	u8		leftFenceCh;		//     .. a marker has been written to
								//     .. v[leftFencePos];  leftFenceCh is the
								//     .. character that was there before
	int		hasRightFence;		//     a 'fence' at the right end of some
	unspos	rightFencePos;		//     .. interval (similar ot the left fence)
	u8		rightFenceCh;		//

	// file containing the sequence

	char*	filename;			// (H) the name of the file associated with this
								//     .. sequence;  this may be NULL
	FILE*	f;					// (F) the file associated with this sequence;
								//     .. this can be NULL;  with NULL, a
								//     .. fileType of seq_type_nofile indicates
								//     .. this sequence was manufactured in
								//     .. memory;  any other fileType indicates
								//     .. a cloned seqeunce
	int		fileType;			//     the type of file being read (one of
								//     .. seq_type_xxx)
	int		rewindable;			//     true  => we can rewind the file
								//     false => we cannot rewind the file
								//     -1    => we do not know if we can rewind
								//              .. the file or not

	u32		pendingLen;			//     characters that have been read from the
	char*	pendingChars;		// (H) .. file but which have not been consumed;
	char*	pendingStack;		// (X) .. pendingChars is an array of size
			    				//     .. seqBufferSize;  pendingLen gives the
			    				//     .. number of non-consumed characters
								//     .. remaining;  pendingStack points to the
								//     .. top of the stack, which builds down
								//     .. from pendingChars+seqBufferSize

	twobit	twoBit;				//     additional info for 2bit files (only
								//     .. valid if fileType == seq_type_2bit)
	hsx		hsx;				//     additional info for hsx files (only
								//     .. valid if fileType == seq_type_hsx)
	seqpartition partition;		//     additional info for partitioned
								//     .. sequences (only valid if partition.p
								//     .. is not NULL)
	qcode*	qCoding;			// (H) table to map quantum symbols to
								//     .. probabilities  (this may be NULL)

	// saved file state

	int		hasSavedState;		//     true => savedFilePos is meaningful
	long int savedFilePos;		//     (as per ftell)

	// sequence read control options

	char*	namesFilename;		// (H) the name of a file containing a list of
								//     .. contig names;  this may be NULL
	FILE*	namesFile;			// (F) file corresponding to namesFilename
	char	nextContigName	 	//     the name of the next contig-of-interest;
			    [maxSequenceName+1];// .. only valid if namesFile is not NULL
	int		contigPending;		//     true => the sequence file has already
								//             .. been scanned to locate the
								//             .. sequence for nextContigName,
								//             .. is positioned appropriately,
								//             .. and that sequence hasn't been
								//             .. loaded yet;  only valid if
								//             .. namesFile is not NULL

	char*	choresFilename;		// (H) the name of a file containing a list of
								//     .. "alignment chores";  this may be NULL
	FILE*	choresFile;			// (F) file corresponding to choresFilename
	int		choresLineNum;		//     line number of the current chore
	chore	chore;				//     the current alignment chore (valid only
								//     .. if choresFile is non-NULL)

	int		subsampleK;			//     specifier for K-of-N subsampling;  only
	int		subsampleN;			//     the Kth sequence of every group of N are
								//     .. processed, with 1<=K<=N;  only valid
								//     .. if subsampleN is greater than zero
	int		subsampleSkip;		//     the current subsampling state;  this is
								//     .. the number of sequences that we will
								//     .. skip before the next one that we keep

	char*	softMaskFilename;	// (H) the name of a file containing
								//     .. soft-masking info to apply;  this
								//     .. may be NULL
	int		softMaskComplement;	//     false => soft-masking replaces any bases
								//              .. in the intervals
								//     true  => soft-masking replaces any bases
								//              .. NOT in the intervals
	char*	xMaskFilename;		// (H) the name of a file containing X-masking
								//     .. info to apply;  this may be NULL
	int		xMaskComplement;	//     false => x-masking replaces any bases in
								//              .. the intervals
								//     true  => x-masking replaces any bases
								//              .. NOT in the intervals
	char*	nMaskFilename;		// (H) the name of a file containing N-masking
								//     .. info to apply;  this may be NULL
	int		nMaskComplement;	//     (similar to xMaskComplement)
	unspos	startLimit,			//     1-based starting and ending locations
			endLimit;			//     .. (inclusive) limiting the part of the
								//     .. sequence we are to read;  these are
								//     .. zero if there are no limits
	int		endIsSoft;			//     true => while the endLimit has been
								//             .. specified, it is a soft limit
								//             .. and can be trimmed to the
								//             .. actual end of the sequence

	int		doRevCompFlags;		//     two bits describing whether any sequence
								//     .. we read is to be reversed or
								//     .. complemented; the four values are the
								//     .. rcf_xxx values defined below
	int		doUnmask;			//     true  => any sequence we read is to be
								//              .. unmasked (converted to upper
								//              .. case)
	int		doPartitioning;		//     true  => the sequence will be partitioned
	int		doJoin;				//     true  => combine the file's sequences
								//              .. into a partitioned sequence
	char	separatorCh;		//     (if not NUL) the file's sequence should
								//     .. be separated into partitions at any
								//     .. run of this character
	int		useFullNames;		//     true  => report full names in alignments
								//     false => report short names instead
	int		nameParseType;		//     how to parse sequence headers, if at all;
								//     .. one of name_parse_type_xxx
	char*	nameTrigger;		//     a string to trigger the fetch of the
								//     .. short sequence name from a fasta
								//     .. header line;  e.g. "name:" means to
								//     .. fetch the name stating after the ":";
								//     .. NULL indicates an absent trigger
	int		allowAmbiDNA;		//     (this applies only to reading fasta
								//      .. files)
								//     true  => permit ambiguous DNA characters
								//              .. B,D,H,K,M,R,S,V,W,Y
								//     false => only A,C,G,T,N,X permitted
	u8*		qToComplement;		// (X) (similar to nuc_to_complement) array to
								//     .. map a quantum base to its complement;
								//     .. this may be NULL
	} seq;


enum
	{
	seq_type_nofile=0,
	seq_type_unknown,
	seq_type_fasta,				// dna, permitting ambiguity
	seq_type_fastq,				// dna with base-call quality scores
	seq_type_csfasta,			// color-space dna
	seq_type_nib,				// nybble-coded dna
	seq_type_2bit,				// 2bit-coded dna
	seq_type_hsx,				// hashed sequence index
	seq_type_qdna,				// quantum-dna (byte-coded prob'listic ambiguity)
	seq_type_max				// (sentinel for enum of seq_type_xxx)
	};

#ifdef sequences_owner
char* seqTypeNames[seq_type_max] =
    {"(no file)", "unknown",
     "fasta", "fastq", "csfasta", "nib", "2bit", "hsx", "qdna" };
#else
extern char* seqTypeNames[];
#endif

enum
	{
	rcf_forward = 0,			// sequence is in normal order, uncomplemented
	rcf_comp    = 1,			// sequence is complemented but not reversed
	rcf_rev     = 2,			// sequence is reversed but not complemented
	rcf_revcomp = 3				// sequence is reversed and complemented
	};

enum
	{
	name_parse_fill_white=1,	// (modifier for types below) convert any
								// whitespace in the resulting name to
								// underline character
	name_parse_type_core=0,		// trim path on left, junk and ext on right
	name_parse_type_alnum=2,	// trim path on left, use only alphanumeric and
								// .. underscore
	name_parse_type_darkspace=4,// trim path on left, use anything but
								// .. whitespace
	name_parse_type_trigger=6	// use nameTrigger
	};

#define parse_type(nameParseType) ((nameParseType) & (~name_parse_fill_white))

#define minFastqCh '!'
#define maxFastqCh '~'

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	u64   partitionLookups;
	u64   partitionHits;
	u64   lookupIterations;
	} sequenceStats;

// stats macros

#define sequence_count_stat(field)   ++sequenceStats.field
#define sequence_uncount_stat(field) --sequenceStats.field
#define sequence_set_stat(field,val) (sequenceStats.field = val)
#define sequence_add_stat(field,val) (sequenceStats.field += val)
#else
#define sequence_count_stat(field)
#define sequence_uncount_stat(field)
#define sequence_set_stat(field,val)
#define sequence_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void sequence_zero_stats (void);
void sequence_show_stats (FILE* f);

//----------
//
// prototypes for routines in sequences.c
//
//----------

seq*  open_sequence_file            (char* name, int fileType,
                                     int choresAllowed, char* choresFilename,
                                     unspos allocLen,
                                     int needTrueLen, int prohibitAmbiDNA,
                                     u8* qToComplement);
seq*  open_rewindable_sequence_file (char* name, int fileType,
                                     int choresAllowed, char* choresFilename,
                                     unspos allocLen,
                                     int needTrueLen, int prohibitAmbiDNA,
                                     u8* qToComplement);
void  rewind_sequence_file          (seq* seq);
seq*  clone_sequence                (seq* seq);
seq*  copy_sequence                 (seq* seq);
seq*  new_sequence                  (unspos allocLen);
void  sequence_long_enough          (seq* seq, unspos allocLen, int anticipate);
void  free_sequence                 (seq* seq);
int   load_sequence                 (seq* seq);
int   another_sequence              (seq* seq);
partition* lookup_partition_no_die  (seq* seq, unspos pos);
partition* lookup_partition         (seq* seq, unspos pos);
partition* lookup_named_partition   (seq* seq, char* name);
partition* last_partition_with_name (seq* seq, partition* firstPart);
partition* lookup_partition_seq_pos (seq* seq, partition* part, unspos pos);
void  print_sequence                (FILE* f, seq* seq, char* header, int perLine);
void  print_partition_table         (FILE* f, seq* _seq);
void  mask_sequence                 (seq* seq, char* maskFilename, int maskChar);
void  mask_sequence_keep            (seq* seq, char* maskFilename, int maskChar);
void  colorize_sequence             (seq* seq);
void  validate_rev_comp             (seq* seq);
void  rev_comp_sequence             (seq* seq, const u8* nucToComplement);
void  backward_sequence             (seq* seq);
void  upper_sequence                (seq* seq);
char* copy_reverse_of_string        (char* s, unspos len);
void  strncpy_reverse               (char* d, char* s, unspos len);
void  fence_sequence_interval       (seq* seq, interval interval, u8 ch);
void  unfence_sequence_interval     (seq* seq);
void  print_file_actions            (FILE* f);
void  match_composition             (seq* seq1, unspos pos1,
                                     seq* seq2, unspos pos2, unspos length,
                                     unspos count[4][4]);
int   percent_identical             (seq* seq1, unspos pos1,
                                     seq* seq2, unspos pos2, unspos length);
score score_match                   (scoreset* scoring,
                                     seq* seq1, unspos pos1,
                                     seq* seq2, unspos pos2, unspos length);
void  dump_aligned_nucleotides      (FILE* f,
                                     seq* seq1, unspos pos1,
                                     seq* seq2, unspos pos2, unspos length);
void  dump_sequence                 (FILE* f, seq* _seq);
void  dump_sequence_state           (FILE* f, seq* _seq);

#undef global
#endif // sequences_H
