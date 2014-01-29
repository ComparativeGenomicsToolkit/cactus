//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: gapped_extend.c
//
//----------
//
// gapped_extend--
//	Support for extending anchors to alignments (gaps allowed).
//
// The Y-drop variant of dynamic programming is applied to create a gapped
// alignment by extending in both directions from each 'anchor point".  We use
// the term Y-drop to distinguish this from the similar X-drop technique for
// ungapped alignments.
//
// The underlying DP algorithm here is the one shown in figure 4 of reference
// [1].  However, where the algorithm in [1] computes the entire DP matrix, we
// only compute horizontal slices of the DP matrix, with the bounds of each row
// determined by (1) any neighboring alignment segments, and (2) cells scoring
// less than Y from the max score.
//
// Throughout this module, we consider sequence 1 ("target") on the vertical
// edge of the dynamic programming matrix, and sequence 2 ("query") horizontal.
// Deletions are vertical, insertions horizontal.
//
// More detailed information can be found in the headers of the functions.
// Specifically, the function ydrop_one_sided_align() is the heart of the DP
// algorithm, and contains more details.
//
// References:
//	[1] Approximate Matching of Regular Expressions.  Myers and Miller, Bull.
//	    Math. Biol. 51 (1989), pp. 5-37.
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
#include <string.h>				// standard C string stuff
#include <limits.h>				// standard C value limit stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "identity_dist.h"		// identity distribution stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "continuity_dist.h"	// query continuity distribution stuff
#include "output.h"				// alignment outout format stuff

#define  gapped_extend_owner	// (make this the owner of its globals)
#include "gapped_extend.h"		// interface to this module

// debugging defines

//#define snoopAnchors			// if this is defined, extra code is added to
								// .. track anchors through the alignment
								// .. process
//#define snoopAnchorToGapped	// if this is defined, extra code is added to
								// .. track anchors through the alignment
								// .. process (different than snoopAnchors)
//#define snoopBlocks			// if this is defined, extra code is added to
								// .. track alignment blocks through the process
//#define snoopSubprobs			// if this is defined, extra code is added to
								// .. report alignment sub-problems
//#define snoopAlgorithm		// if this is defined, extra code is added to
								// .. track the dynamic programming algorithm
//#define snoopAlgorithmTrap	// if this is defined, extra code is added to
								// .. allow trapping the dynamic programming
								// .. algorithm in a debugger, at a specific
								// .. cell
//#define debugPosA1 525189		// if defined (and snoopAlgorithm is also),
//#define debugPosA2 501468		// .. only information about extension of this
//#define debugPosB1 55189		// .. particular seed is output;  note that
//#define debugPosB2 91468		// .. this position is origin-zero

//#define snoopTraceback		// if this is defined, extra code is added to
								// .. track the dynamic programming algorithm's
								// .. traceback process

//#define debugHspImmediate		// if this is defined, extra code is added to
								// .. aid debugging of gappily_extend_hsps()

//#define snoopBatches			// if this is defined, extra code is added to
								// .. track the processing of HSPs in batches;
								// .. this is not meaningful if 
//#define snoopEditScripts		// if this is defined, extra code is added to
								// .. examine the edit scripts created by
								// .. format_alignment()
//#define snoopAlignioInput		// if this is defined, extra code is added to
								// .. examine the alignio data structure, which
								// .. is used to pass data into ydrop_align()
//#define snoopAlignioOutput	// if this is defined, extra code is added to
								// .. examine the results in the alignio data
								// .. structure, which is used to pass data out
								// .. of ydrop_align()

//#define snoopSpecialHsp		// if this is defined, extra code is added to
//								// .. so that only one particular HSP is
//								// .. processed
//#define specialPosA 439		// if defined (and snoopSpecialHsp is also),
//#define specialPosB 16		// .. only this HSP is processed;  note that
//								// .. this position is origin-zero

//#define snoopBounds			// if this is defined, extra code is added to
								// .. examine bounding alignments

//----------
//
// stats to augment crude profiling
//
//----------

//--- this module's contribution to profiling for the whole program --

#ifndef dbgTiming
#define dbg_timing_set_stat(field,val)     ;
#define dbg_timing_count_stat(field)       ;
#define dbg_timing_report_stat(field,name) ;
#endif // not dbgTiming

#ifdef dbgTiming
struct
	{
	int numExtensions;
	} gappedExtendTimingStats;

#define dbg_timing_set_stat(field,val)     (gappedExtendTimingStats.field = val)
#define dbg_timing_count_stat(field)       ++gappedExtendTimingStats.field
#define dbg_timing_report_stat(field,name) fprintf(stderr,"%-26s %d\n",name":",gappedExtendTimingStats.field)
#endif // dbgTiming

//--- profiling for this module only --

#ifndef dbgTimingGappedExtend
#define dbg_timing_gapped_extend_sub(v) ;
#define dbg_timing_gapped_extend_add(v) ;
#define dbg_timing_gapped_extend_copy(dst,src) ;
#define dbg_timing_gapped_extend_report(f,v,s) ;
#endif // not dbgTimingGappedExtend

#ifdef dbgTimingGappedExtend
#define read_clock() microsec_clock()
#define clocksPerSec 1000000
u64 microsec_clock (void);  // (from lastz.c)

s64 debugClockAboveBelow         = 0,
    debugClockLeftRight          = 0,
    debugClockYdropAlign         = 0,
    debugClockYdropOneSidedAlign = 0,
    debugClockUpdateLrBounds     = 0,
    debugClockNextSweepSeg       = 0,
    debugClockPrevSweepSeg       = 0,
    debugClockUpdateActiveSegs   = 0,
    debugClockFilterActiveSegs   = 0;

#define dbg_timing_gapped_extend_sub(v)  { v -= (s64) read_clock();  }
#define dbg_timing_gapped_extend_add(v)  { v += (s64) read_clock();  }
#define dbg_timing_gapped_extend_copy(dst,src) { dst = src; }

#define dbg_timing_gapped_extend_report(f,v,s) { fprintf(f,"%-26s %.3f\n",s":",((float)(v))/clocksPerSec); }
#endif // dbgTimingGappedExtend

//----------
//
// private data
//
//----------

// miscellany

#define negInf negInfinity

#define anchorPeakLen 31

#undef  min
#define min(x,y) ((x)<(y)?(x):(y))
#undef  max
#define max(x,y) ((x)>(y)?(x):(y))

#define signed_difference(a,b) (((sgnpos)(a))-((sgnpos)(b)))

// straight alignment segments

typedef struct aliseg
	{
	char		type;			// one of diagSeg, horzSeg, vertSeg
	unspos		b1, b2, e1, e2;
	struct aliseg* nextSeg;
	struct aliseg* prevSeg;
	} aliseg;

#define diagSeg 0
#define horzSeg 1				// parallel to sequence 2
#define vertSeg 2				// parallel to sequence 1

// gapped alignment;  this starts as a single point anchor, and is then extended
// into an alignment

typedef struct galign
	{
	unspos		pos1, pos2;		// anchor location or start of alignment (these
								// .. are origin-zero)
	unspos		end1, end2;		// end of alignment (these are inclusive)

	aliseg*		firstSeg;		// the alignment, in diagonal, vertical and
	aliseg*		lastSeg;		// .. horizontal segments

	alignel*	align;			// second form of the alignment, for external
								// .. consumption only

	struct galign *leftAlign1;	// at the alignment's beginning and ending
	struct galign *rightAlign1;	// .. points, we keep pointers to the alignments
	struct galign *leftAlign2;	// .. immediately left and right of the terminus
	struct galign *rightAlign2;

	aliseg*		leftSeg1;		// these correspond to leftAlign1, etc., and
	aliseg*		rightSeg1;		// .. are the closest segments to the alignment
	aliseg*		leftSeg2;		// .. (e.g. leftSeg1 is a segment in leftAlign1)
	aliseg*		rightSeg2;

	struct galign* next;		// alignments are linked both by increasing
	struct galign* prev;		// .. start-point and by decreasing end-point
	} galign;

// input/output for ydrop_align()

typedef struct alignio
	{
	// the following must be supplied by caller

	u8*		seq1, *seq2;		// target and query sequences
	u8*		rev1, *rev2;		// reverse of target and query sequences (NOT
								// .. reverse complement)
	unspos	len1, len2;			// total length of sequences
	unspos	low1, low2;			// limits of sub-interval in each sequence;  low
	unspos	high1, high2;		// .. is the leftmost position allowed;  high is
								// .. one past the rightmost position  allowed;
								// .. both are origin-0 indexes into seq->v[]
	unspos	anchor1, anchor2;	// position to start the alignment

	scoreset* scoring;			// substitution and gap scores to use
	score	yDrop;				// value of Y-dropoff parameter
	int		trimToPeak;			// whether y-drop should be trimmed (see
								// description in ydrop_one_sided_align)
	tback*	tb;					// block of memory for alignment trace-back

	galign* leftAlign;			// closest alignment to left of anchor pt
	galign* rightAlign;			// closest alignment to right of anchor pt
	aliseg* leftSeg;			// specific segment of leftAlign
	aliseg* rightSeg;			// specific segment of rightAlign

	galign* aboveList;			// alignments starting above the anchor point
	galign* belowList;			// alignments ending below the anchor point

	// the following are returned by ydrop_align()

	score	s;					// alignment's score
	unspos	start1, start2;		// alignment's start on target and query
	unspos	stop1,  stop2;		// alignment's end on target and query
	editscript* script;			// alignment's edit script
	} alignio;

// dynamic programming structure
//
// We (try to) declare the basic cell with a power-of-two size;  doing so MAY
// improve speed on some platforms, presumably because random indexing into the
// array of dpCells would involve a shift rather than a multiply;  however, if
// we can't figure out what size the padding should be (or if we don't need any,
// as would be the case for dpCell_padding_sz==16), we don't add any

#define unpadded_dpCell_sz (score_sz+score_sz+unspos_sz)

#if (unpadded_dpCell_sz == 12)
#define dpCell_padding_sz 4
#elif (unpadded_dpCell_sz == 20)
#define dpCell_padding_sz 12
#elif (unpadded_dpCell_sz == 24)
#define dpCell_padding_sz 8
#endif


typedef struct dpCell
	{
	score	DD, CC;
	unspos	mask;		// mask out grid points that are part of previous
						// .. alignments
#ifdef dpCell_padding_sz
	char	padding[dpCell_padding_sz];
#endif
	} dpCell;

typedef struct dpMatrix
	{
	dpCell*	p;
	u32		len;
	} dpMatrix;

// linked list of active segments

typedef struct activeseg
	{
	aliseg*	seg;
	unspos	x;			// column position where segment intersects sweep row;
						// .. this is relative to the current slice of the DP
						// .. matrix, e.g. ranging from LY to RY in the context
						// .. of the alignment sweep
	unspos	lastRow;	// (also relative to the alignment sweep)
	char	type;		// one of diagSeg, horzSeg, vertSeg
	char	filter;
	struct activeseg* next;
	} activeseg;


#ifdef snoopSubprobs
seq* snoopSubprobsSeq1;
seq* snoopSubprobsSeq2;
#endif // snoopSubprobs


// segment batches-- partitioning of a segment table

typedef struct segbatch
	{
	u32		start;				// index (into a segment table) of the first
								// .. entry in a batch
	u32		end;				// index (into a segment table) of the first
								// .. entry NOT in a batch (i.e. the one after
								// .. the last entry).
	partition*	part;			// sequence partition that "contains" this
								// .. batch;  this can be NULL if we aren't
								// .. dealing with partitions
	} segbatch;

typedef struct sbtable
	{
	u32		size;				// the number of entries allocated for batch[]
	u32		len;				// the number of batches (the number of entries
								// .. actually used)
	segbatch batch[1];			// the batch table (variable-length array)
	} sbtable;

#define sbtable_bytes(size) (sizeof(sbtable) + (((size)-1)*sizeof(segbatch)))

//----------
//
// prototypes for private functions
//
//----------

static unspos   segment_peak       (u8* s1, u8* s2, unspos segLength,
                                    scoreset* scoring);
static sbtable* batched_segments   (segtable* anchors, seqpartition* sp);
static galign** init_from_anchors  (segtable* anchors, u32 numExtraSlots);
static int      identical_sequences (seq* seq1, seq* seq2,
                                    scoreset* scoring, score* s);
static int      identical_partitioned_sequences
                                   (seq* seq1, seq* seq2);
static int      identical_partition_of_sequence
                                   (seq* seq1, seq* seq2);
static score    score_identical_partition
                                   (seq* seq1, seq* seq2,
                                    partition* p1, partition* p2,
                                    scoreset* scoring);
static score    score_identical_partition_of
                                   (seq* seq1, seq* seq2,
                                    partition* p1,
                                    scoreset* scoring);
static void     ydrop_align        (alignio* io);
static score    ydrop_one_sided_align (alignio* io, int reversed,
                                    u8* A, u8* B, unspos M, unspos N,
                                    int trimToPeak,
                                    editscript** script,
                                    unspos* end1, unspos* end2);
static void     dp_ready           (dpMatrix* dynProg, unspos needed);
static int      msp_left_right     (galign* obi, galign* m);
static void     get_above_below    (alignio* io, galign* obi, galign* oed);
static void     align_left_right   (galign* obi, galign* m);
static void     insert_align       (galign* m, galign** obi, galign** oed);
static void     update_LR_bounds   (int reversed,
                                    aliseg** rightSeg,   aliseg** leftSeg,
                                    galign** rightAlign, galign** leftAlign,
                                    unspos row, unspos anchor1, unspos anchor2,
                                    sgnpos* L, sgnpos* R,
                                    unspos* LY, unspos* RY);
static sgnpos   next_sweep_seg     (int lookRight, aliseg** bp, galign** mp,
                                    unspos row, unspos anchor1, unspos anchor2);
static sgnpos   prev_sweep_seg     (int lookRight, aliseg** bp, galign** mp,
                                    unspos row, unspos anchor1, unspos anchor2);
static void     update_active_segs (int reversed, activeseg** active,
                                    galign** alignList, dpCell* dp,
                                    unspos row, unspos anchor1, unspos anchor2,
                                    unspos LY, unspos RY);
static void     build_active_seg   (int reversed, activeseg* act, dpCell* dp,
                                    unspos row, unspos anchor1, unspos anchor2,
                                    unspos LY, unspos RY);
static activeseg* add_new_active   (int reversed, activeseg* active,
                                    galign* alignList, dpCell* dp,
                                    unspos row, unspos anchor1, unspos anchor2,
                                    unspos LY, unspos RY);
static void     filter_active_segs (activeseg** active, int filter);
static alignel* format_alignment   (alignio* io, galign* m);
static void     save_seg           (galign* m,
                                    unspos b1, unspos b2, unspos e1, unspos e2);

#ifdef snoopAlignioInput
static void     dump_alignio_input (FILE* f, alignio* io);
#endif // snoopAlignioInput
#ifdef snoopAlignioOutput
static void     dump_alignio_output(FILE* f, alignio* io);
#endif // snoopAlignioOutput

static u64      count_paired_bases (galign* mp);
static void     warn_for_paired_bases_limit (seq* seq2, u64 maxPairedBases,
                                    int overlyPairedKeep);

//----------
//
// reduce_to_points--
//	Convert each segment in a table to it's "peak".  The definition of peak is
//	the midpoint of the highest scoring subsegment of a given length.
//
//----------
//
// Arguments:
//	seq*		seq1:		The first sequence.
//	seq*		seq2:		The second sequence.
//	scoreset*	scoring:	The scoring scheme to use.
//	segtable*	anchors:	The segment table to modify.
//
// Returns:
//	(nothing)
//
//----------

void reduce_to_points
   (seq*		seq1,
	seq*		seq2,
	scoreset*	scoring,
	segtable*	anchors)
	{
	u32			ix;
	segment*	seg;
	unspos		peak;

	if (gapped_extend_dbgShowAnchors)
		{
		if (anchors->len == 0)
			fprintf (stderr, "reduce_to_points: no anchors\n");
		else
			fprintf (stderr, "reduce_to_points: %u anchors\n", anchors->len);
		}

	for (ix=0,seg=anchors->seg ; ix<anchors->len ; ix++,seg++)
		{
		if (gapped_extend_dbgShowAnchors)
			fprintf (stderr, "reduce_to_points: reducing ("
			                 unsposSlashFmt " " unsposFmt
			                 " diag=" sgnposFmt
			                 " score=" scoreFmtSimple ")",
			                 seg->pos1, seg->pos2, seg->length,
			                 diagNumber(seg->pos1,seg->pos2),
			                 seg->s);

		peak = segment_peak (seq1->v+seg->pos1, seq2->v+seg->pos2, seg->length,
		                     scoring);
		seg->pos1   += peak;
		seg->pos2   += peak;
		seg->length =  0;

		if (gapped_extend_dbgShowAnchors)
			fprintf (stderr, " to (" unsposSlashFmt " " unsposFmt ")\n",
			                 seg->pos1, seg->pos2, seg->length);
		}

	gapped_extend_add_stat (numAnchors, anchors->len);
	}


static unspos segment_peak
   (u8*			s1,
	u8*			s2,
	unspos		segLength,
	scoreset*	scoring)
	{
	u8*			t1 = s1;
	u8*			t2 = s2;
	score		similarity, best;
	unspos		ix, peak;

	if (segLength <= anchorPeakLen)
		peak = segLength / 2;
	else
		{
		similarity = 0;
		for (ix=0 ; ix<anchorPeakLen; ix++)
			{
			similarity += scoring->sub[*t1++][*t2++];
			//fprintf (stderr, "%c %c " scoreFmtSimple " " scoreFmtSimple "\n",
			//                 t1[-1], t2[-1], scoring->sub[t1[-1]][t2[-1]], similarity);
			}
		best = similarity;
		peak = anchorPeakLen / 2;

		for ( ; ix<segLength; ix++)
			{
			similarity -= scoring->sub[*s1++][*s2++];
			similarity += scoring->sub[*t1++][*t2++];
			//fprintf (stderr, "%c %c " scoreFmtSimple " " scoreFmtSimple "\n",
			//                 t1[-1], t2[-1], scoring->sub[t1[-1]][t2[-1]], similarity);
			if (similarity > best)
				{
				best = similarity;
				peak = ix - (anchorPeakLen / 2);
				}
			}

		gapped_extend_count_stat (numPeaks);
		gapped_extend_add_stat   (totalPeakScore, best);
		}

	return peak;
	}

//----------
//
// gapped_extend--
//	performed gapped extension given a set of anchor segments.
//
//----------
//
// Arguments:
//	seq*		seq1:			The sequence being searched.
//	u8*			rev1:			The reverse (NOT reverse complement) of seq1,
//								.. as a zero-terminated string.
//	seq*		seq2:			The sequence being searched for.
//	u8*			rev2:			The reverse of the seq2 (analagous to rev1).
//	int			inhibitTrivial: true => don't output the trivial self-alignment.
//	scoreset*	scoring:		The scoring scheme to use.
//	segtable*	anchors:		The anchor segments.
//	void*		tb:				Memory in which to track gapped alignment
//								.. traceback.
//	int			allBounds:		true  => bound gapped alignments by *all* gapped
//								         .. extensions of higher-scoring HSPs (a
//								         .. la blastz)
//								false => bound gapped alignments only by gapped
//								         .. extensions that meet the score
//								         .. threshold
//	score		yDrop:			Threshold to stop gapped extensions;  if the
//								.. score drops off by more than yDrop, extension
//								.. stops
//	int			trimToPeak:		Whether y-drop should be trimmed (see
//								.. description in ydrop_one_sided_align).
//	sthresh		scoreThresh:	Minimum score required;  gapped alignments are
//								.. discarded if they score less than this.
//	u64			maxPairedBases:	Maximum number of "paired bases" we'll allow.
//								.. A paired base is a match or substitution in
//								.. the DP matrix.  If we exceed this limit, we
//								.. abort processing.  However, any gapped
//								.. alignments we find up to that point are part
//								.. of what we return.  A value of zero inidcates
//								.. there is no limit.
//	int			overlyPairedWarn: true => write a warning message when a query
//								..        .. exceeds the maxPairedBases
//								..        .. threshold.
//	int			overlyPairedKeep:
//								How we should treat alignments for queries that
//								.. exceed the maxPairedBases threshold.
//								..   false => we discard all the alignments
//								..   true  => we output whatever alignments we
//								..            .. happened to find prior to
//								..            .. exceeding the limit
//
// Returns:
//	A linked list of alignments.  The caller is responsible for disposing of
//	these.  The caller must also deallocate other memory;  see note (1) below.
//
//----------
//
// Notes:
//
// (1)	This routine calls other routines that allocate long-term memory.  The
//		Calling program is repsonsible for making sure that this memory is
//		disposed of at the appropriate time (usually at program exit).  This
//		should be done by calling free_segment_batches().
//
//----------

//=== stuff for gapped_extend_verbosity ===

#define debugGappedExtendVerbosity_1                                         \
	if (gapped_extend_verbosity >= 2)                                        \
		{                                                                    \
		pos1 = mp->pos1;                                                     \
		pos2 = mp->pos2;                                                     \
		                                                                     \
		if (sp1->p != NULL)                                                  \
			{                                                                \
			part =  lookup_partition (seq1, pos1);                           \
			pos1 += part->sepBefore + 1;                                     \
			}                                                                \
		if (sp2->p != NULL)                                                  \
			{                                                                \
			part =  lookup_partition (seq2, pos2);                           \
			pos2 += part->sepBefore + 1;                                     \
			}                                                                \
		                                                                     \
		pos1 += seq1->startLoc;                                              \
		pos2 += seq2->startLoc;                                              \
		                                                                     \
		fprintf (stderr, "processing anchor #%u (of %u)"                     \
		                 " (" unsposSlashFmt ")"                             \
		                 " "  unsposSlashFmt "\n",                           \
						 i+1, anchors->len,                                  \
						 mp->pos1, mp->pos2,                                 \
						 pos1,     pos2);                                    \
		}

#define debugGappedExtendVerbosity_2                                         \
	if (gapped_extend_verbosity >= 2)                                        \
		{                                                                    \
		pos1 = mp->pos1;      len1 = mp->end1 - pos1;                        \
		pos2 = mp->pos2;      len2 = mp->end2 - pos2;                        \
		                                                                     \
		if (sp1->p != NULL)                                                  \
			{                                                                \
			part =  lookup_partition (seq1, pos1);                           \
			pos1 += part->sepBefore + 1;                                     \
			}                                                                \
		if (sp2->p != NULL)                                                  \
			{                                                                \
			part =  lookup_partition (seq2, pos2);                           \
			pos2 += part->sepBefore + 1;                                     \
			}                                                                \
		                                                                     \
		pos1 += seq1->startLoc;                                              \
		pos2 += seq2->startLoc;                                              \
		                                                                     \
		fprintf (stderr, "alignment block"                                   \
						 " score=" scoreFmtSimple                            \
						 " at (" unsposSlashFmt ") " unsposSlashFmt          \
						 " length " unsposSlashFmt "\n",                     \
						 mp->align->s,                                       \
						 mp->pos1, mp->pos2, pos1, pos2,                     \
						 len1, len2);                                        \
	}


//=== stuff for gapped_extend_dbgShowIdentity ===
// nota bene: positions reported below are 1-based (not 0-based)

#define debugGappedExtendDbgShowIdentity_1                                   \
	if (gapped_extend_dbgShowIdentity)                                       \
		printf ("discarding " unsposSlashSFmt " score=" scoreFmtSimple "\n", \
				mp->pos1+1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->pos2+1, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->align->s);

#define debugGappedExtendDbgShowIdentity_2                                   \
	if (gapped_extend_dbgShowIdentity)                                       \
		printf ("discarding " unsposSlashSFmt " score=" scoreFmtSimple "\n", \
				mp->pos1+1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->pos2+1, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->align->s);

#define debugGappedExtendDbgShowIdentity_3                                   \
	if (gapped_extend_dbgShowIdentity)                                       \
		printf ("keeping " unsposSlashSFmt " score=" scoreFmtSimple "\n",    \
				mp->pos1+1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->pos2+1, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+", \
				mp->align->s);


//=== stuff for snoopAnchors ===

#ifndef snoopAnchors
#define debugSnoopAnchors_1 ;
#endif // not snoopAnchors

#ifdef snoopAnchors

#define debugSnoopAnchors_1                                                  \
	{                                                                        \
	u32 ix;                                                                  \
	fprintf (stderr,"===== msp =====\n");                                    \
	for (ix=0 ; ix<anchors->len ; ix++)                                      \
		{                                                                    \
		fprintf (stderr,"anchors[%3d] " unsposSlashFmt " " unsposFmt " " scoreFmt, \
		         ix, anchors->seg[ix].pos1+1,                                \
		             anchors->seg[ix].pos2+1,                                \
		             anchors->seg[ix].length,                                \
		             anchors->seg[ix].s);                                    \
		fprintf (stderr,"  msp[%3d] " unsposSlashFmt " " unsposSlashFmt "\n",  \
		         ix, msp[ix]->pos1+1, msp[ix]->pos2+1,                       \
		             msp[ix]->end1,   msp[ix]->end2);                        \
		}                                                                    \
	}

#endif // snoopAnchors


//=== stuff for snoopAnchorToGapped ===

#ifndef snoopAnchorToGapped
#define debugSnoopAnchorToGapped_1 ;
#define debugSnoopAnchorToGapped_2 ;
#define debugSnoopAnchorToGapped_3 ;
#define debugSnoopAnchorToGapped_4 ;
#define debugSnoopAnchorToGapped_5 ;
#endif // not snoopAnchorToGapped

#ifdef snoopAnchorToGapped

#define debugSnoopAnchorToGapped_1                                           \
	fprintf (stderr, "processing anchor #%u (of %u)\n",                      \
					 i+1, anchors->len);

#define debugSnoopAnchorToGapped_2                                           \
	fprintf (stderr, "anchor: " unsposSlashFmt " (diag " sgnposFmt ")\n",    \
	                  mp->pos1, mp->pos2, diagNumber(mp->pos1,mp->pos2));

#define debugSnoopAnchorToGapped_3                                           \
	fprintf (stderr, "  gapped: " unsposSlashFmt " " unsposSlashFmt " " scoreFmt "\n", \
	                 io.start1, io.start2, io.stop1, io.stop2, io.s);

#define debugSnoopAnchorToGapped_4                                           \
	fprintf (stderr, "finished processing %u anchors\n", anchors->len);      \
	fprintf (stderr, "  head=%p\n", head);                                   \
	for (a=head,i=0 ; a!=NULL ; a=a->next,i++)                               \
		fprintf (stderr, "  [%u] %p"                                         \
		                 " %s:" unsposDotsFmt " %s:" unsposDotsFmt "\n",     \
			             i, a,                                               \
			             seq1->header, a->beg1, a->end1,                     \
			             seq2->header, a->beg2, a->end2);

#define debugSnoopAnchorToGapped_5                                           \
	fprintf (stderr, "finished processing %u anchors\n", anchors->len);      \
	fprintf (stderr, "  paired threshold was exceeded\n");

#endif // snoopAnchorToGapped

//=== stuff for snoopBlocks ===

#ifndef snoopBlocks
#define debugSnoopBlocks_1 ;
#define debugSnoopBlocks_2 ;
#define debugSnoopBlocks_3 ;
#define debugSnoopBlocks_3b ;
#define debugSnoopBlocks_4 ;
#define debugSnoopBlocks_5 ;
#endif // not snoopBlocks

#ifdef snoopBlocks

#define debugSnoopBlocks_1                                                   \
	fprintf (stderr, "===== searching for alignment blocks =====\n");

#define debugSnoopBlocks_2                                                   \
	fprintf (stderr, "discarding alignment block [%8p]"                      \
	                 " b " unsposFmt " " unsposFmt                           \
	                 " e " unsposFmt " " unsposFmt                           \
	                 " s " scoreFmtSimple "\n",                              \
	                 mp, mp->pos1, mp->pos2, mp->end1, mp->end2, mp->align->s);

#define debugSnoopBlocks_3                                                   \
	fprintf (stderr, "===== collecting alignment blocks =====\n");

#define debugSnoopBlocks_3b                                                  \
	fprintf (stderr, "===== discarding alignment blocks =====\n");

#define debugSnoopBlocks_4                                                   \
	fprintf (stderr, "keeping alignment block [%8p -> %8p]"                  \
	                 " b " unsposFmt " " unsposFmt                           \
	                 " e " unsposFmt " " unsposFmt                           \
	                 " s " scoreFmtSimple "\n",                              \
	                 mp, mp->align,                                          \
	                 mp->pos1, mp->pos2, mp->end1, mp->end2, mp->align->s);

#define debugSnoopBlocks_5                                                   \
	fprintf (stderr, "discarding alignment block [%8p]"                      \
	                 " b " unsposFmt " " unsposFmt                           \
	                 " e " unsposFmt " " unsposFmt                           \
	                 " s " scoreFmtSimple "\n",                              \
	                 mp, mp->pos1, mp->pos2, mp->end1, mp->end2, mp->align->s);

#endif // snoopBlocks


//=== stuff for snoopSubprobs ===

#ifndef snoopSubprobs
#define debugSnoopSubprobs_1 ;
#endif // not snoopSubprobs

#ifdef snoopSubprobs

#define debugSnoopSubprobs_1                                                 \
	snoopSubprobsSeq1 = seq1;                                                \
	snoopSubprobsSeq2 = seq2;

#endif // snoopSubprobs


//=== stuff for snoopBatches ===

#ifndef snoopBatches
#define debugSnoopBatches_1 ;
#endif // not snoopBatches

#ifdef snoopBatches

#define debugSnoopBatches_1                                                  \
	if (doHspsInBatches)                                                     \
		{                                                                    \
		partition* batPart = segBatches->batch[batIx].part;                  \
		fprintf (stderr, "batch[%u] %u..%u",                                 \
		                 batIx, startSegIx, endSegIx-1);                     \
		if (batPart != NULL)                                                 \
			fprintf (stderr, " " unsposFmt ".." unsposFmt,                   \
			                 batPart->sepBefore+1, batPart->sepAfter);       \
		fprintf (stderr, "\n");                                              \
		}

#endif // snoopBatches


//=== stuff for tryout ===

#ifndef tryout
#define debugTriviality_1 ;
#define debugTriviality_2 ;
#define debugTriviality_3 ;
#define debugTriviality_4 ;
#define debugTriviality_5 ;
#define debugTriviality_6 ;
#endif // not tryout

#ifdef tryout

#define debugTriviality_1                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "trivial?: \"%s\" " unsposFmt " \"%s\" " unsposFmt "\n", \
						 name1, len1, name2, len2);                          \
		}

#define debugTriviality_2                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "  sequence lengths differ\n");                     \
		}

#define debugTriviality_3                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "  alignment lengths differ\n");                    \
		}

#define debugTriviality_4                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "  sequence names differ\n");                       \
		}

#define debugTriviality_5                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "  alignment content differs\n");                   \
		}

#define debugTriviality_6                                                    \
	if (gapped_extend_dbgTriviality)                                         \
		{                                                                    \
		fprintf (stderr, "  it's trivial!\n");                               \
		}

#endif // tryout


//=== stuff for dbgTimingGappedExtend ===

#ifdef dbgTimingGappedExtend
void gapped_extend_timing_report
   (arg_dont_complain(FILE* f))
	{
	dbg_timing_gapped_extend_report (f, debugClockAboveBelow,         "total time in above_below()");
	dbg_timing_gapped_extend_report (f, debugClockLeftRight,          "total time in left_right()");
	dbg_timing_gapped_extend_report (f, debugClockYdropAlign,         "total time in ydrop_align()");
	dbg_timing_gapped_extend_report (f, debugClockYdropOneSidedAlign, "    ydrop_one_sided_align()");
	dbg_timing_gapped_extend_report (f, debugClockUpdateLrBounds,     "    update_lr_bounds()");
	dbg_timing_gapped_extend_report (f, debugClockNextSweepSeg,       "    next_sweep_seg()");
	dbg_timing_gapped_extend_report (f, debugClockPrevSweepSeg,       "    prev_sweep_seg()");
	dbg_timing_gapped_extend_report (f, debugClockUpdateActiveSegs,   "    update_active_segs()");
	dbg_timing_gapped_extend_report (f, debugClockFilterActiveSegs,   "    filter_active_segs()");
	}
#endif // dbgTimingGappedExtend


//=== stuff for snoopEditScripts ===

#ifndef snoopEditScripts
#define debugSnoopEditScripts_1 ;
#define debugSnoopEditScripts_2 ;
#define debugSnoopEditScripts_3 ;
#define debugSnoopEditScripts_4 ;
#define debugSnoopEditScripts_5 ;
#endif // not snoopEditScripts

#ifdef snoopEditScripts

#define debugSnoopEditScripts_1                                              \
	fprintf (stderr, "(adding full length trivial alignment)\n");

#define debugSnoopEditScripts_2                                              \
	fprintf (stderr, "(adding full length trivial alignment vs. partition %u)\n", \
	                 trivialPartIx);

#define debugSnoopEditScripts_3                                              \
	fprintf (stderr, "(adding trivial alignment for partition %u)\n",        \
	                 partIx);

#define debugSnoopEditScripts_4                                              \
	{                                                                        \
	char*    opName[4] = { "?", "I", "D", "S" };                             \
	alignel* aa = mp->align;                                                 \
	u32      opIx, op, rpt;                                                  \
	fprintf (stderr, "  " unsposDotsFmt " vs " unsposDotsFmt                 \
	                 " (diag " sgnposFmt ")"                                 \
	                 " score=" scoreFmt,                                     \
	                 aa->beg1, aa->end1, aa->beg2, aa->end2,                 \
	                 diagNumber(aa->beg1,aa->beg2),                          \
	                 aa->s);                                                 \
	if (aa->isTrivial) fprintf (stderr, " (trivial)");                       \
	fprintf (stderr, "\n   ");                                               \
	for (opIx=0 ; opIx<aa->script->len ; opIx++)                             \
		{                                                                    \
		op  = edit_op_operation (aa->script->op[opIx]);                      \
		rpt = edit_op_repeat    (aa->script->op[opIx]);                      \
		fprintf (stderr, " %dx%s", rpt, opName[op]);                         \
		}                                                                    \
	fprintf (stderr, "\n");                                                  \
	}

#define debugSnoopEditScripts_5                                              \
	fprintf (stderr, "  (alignment is empty)\n");

#endif // snoopEditScripts


//=== stuff for snoopSpecialHsp ===

#ifndef snoopSpecialHsp
#define debugsnoopSpecialHsp_1 ;
#endif // not snoopSpecialHsp

#ifdef snoopSpecialHsp

#define debugsnoopSpecialHsp_1                                               \
	{                                                                        \
	if ((mp->pos1 != specialPosA) || (mp->pos2 != specialPosB))              \
		{                                                                    \
		fprintf (stderr,"  Ignoring msp[%3d] " unsposSlashFmt " " unsposSlashFmt "\n", \
		                i, mp->pos1+1, mp->pos2+1, mp->end1, mp->end2);      \
		continue;                                                            \
		}                                                                    \
	}

#endif // snoopSpecialHsp


//=== finally, the actual function gapped_extend() ===

alignel* gapped_extend
   (seq*			seq1,
	u8*				rev1,
	seq*			seq2,
	u8*				rev2,
	int				inhibitTrivial,
	scoreset*		scoring,
	segtable*		anchors,
	tback*			tb,
	int				allBounds,
	score			yDrop,
	int				trimToPeak,
	sthresh			scoreThresh,
	u64				maxPairedBases,
	int				overlyPairedWarn,
	int				overlyPairedKeep)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		p1, *p2;
	galign**		msp;
	alignel*		a;
	score			s;
	alignio			io;
	galign*			mp, *mpNext;
	galign*			orderBegInc, *orderEndDec;
	alignel*		head, *last;
	aliseg*			bp, *bq;
	u32				i;
	unspos			pos1, pos2, len1, len2;
	partition*		part;
	int				partitionedTriviality, delayedCheckForTrivial;
	int				doHspsInBatches;
	sbtable*		segBatches;
	u32				batIx, startSegIx, endSegIx, partIx;
	int				trivialPartIx = -1;
	u32				freeSlot;
	u64				pairedBases = 0, newPairedBases;

	orderBegInc = orderEndDec = NULL;	// (compiler appeasement)

	if (scoreThresh.t != 'S')
		suicidef ("gapped_extend can't handle score threshold %s",
		          score_thresh_to_string (&scoreThresh));

	// create a gapped alignment table containing one entry for each HSP (plus
	// an additional slot for each possible trivial self-alignment-- see
	// note (1) in init_from_anchors);  note that batched_segments sorts the
	// HSPs in descreasing score order (per batch)

	doHspsInBatches = false;
	if (gapped_extend_dbgAllowBatches)
		doHspsInBatches = (sp1->p != NULL);  // (seq1 is partitioned)

	if (!doHspsInBatches) segBatches = batched_segments (anchors, NULL);
	                 else segBatches = batched_segments (anchors, sp1);

	if ((sp1->p == NULL) || (sp2->p == NULL))
	 	msp = init_from_anchors (anchors, 1);
	else
		{
		u32 numExtraSlots = sp1->len;
		if (sp2->len < numExtraSlots) numExtraSlots = sp2->len;
		msp = init_from_anchors (anchors, numExtraSlots);
		}

	debugSnoopAnchors_1;
	debugSnoopBlocks_1;
	debugSnoopSubprobs_1;

	// set up the "io" block, which is used for communication with lower-level
	// routines.

	io.seq1 = seq1->v;
	io.seq2 = seq2->v;
	io.rev1 = rev1;
	io.rev2 = rev2;
	io.low1 = 0;   io.len1 = io.high1 = seq1->len;
	io.low2 = 0;   io.len2 = io.high2 = seq2->len;

	io.scoring    = scoring;
	io.yDrop      = yDrop;
	io.trimToPeak = trimToPeak;

	if (tb == NULL)
		suicide ("gapped_extend was given a NULL traceback pointer.");
	io.tb = tb;

	// special case check for trivial alignments of identical sequences; if the
	// two sequences are identical, we add the trivial self alignment to the
	// table immediately (so it can prevent some other alignment near the main
	// diagonal from merging with major portions of the trivial alignment);
	// note that if trivial alignments aren't desired we will discard them
	// later, after we've extended all the other anchors
	//
	// see note (1) in init_from_anchors-- the gapped alignment table has
	// extra slots to allow for the addition of these trivial self-alignments;
	// there is one extra slot for min(nTarget,nQuery), where e.g. nTarget is
	// the number of partitions in the target (or 1, if it is not partitioned).
	//
	// the check for trivialty becomes more complicated if we are dealing with
	// partitioned sequences
	//
	// also note that the more complicated test provides *no guarantee* that
	// it will prevent other other alignments from merging into the diagonal;
	// it is possible that two identical sequences will end up with an aligment
	// that only has part of the main diagonal, and no alignment covering the
	// complete main diagonal;  the author does not consider this a significant
	// failure because it can only happen when the user does *not* tell the
	// program she is doing a self-alignment

	partitionedTriviality  = false;
	delayedCheckForTrivial = ((inhibitTrivial)
	                       && ((sp1->p != NULL) || (sp2->p != NULL)));

	if ((sp1->p != NULL)  // (seq1 is partitioned)
	 && (sp2->p == NULL)) // (seq2 is not partitioned)
	 	{
		// $$$ this could be modified to check seq2 vs *every* partition in
		// $$$ .. seq1, but would require some restructuring of the code;  we
		// $$$ .. don't expect there are several identical sequences in seq1
		trivialPartIx = identical_partition_of_sequence (seq1, seq2);
		partitionedTriviality  = (trivialPartIx != -1);
		delayedCheckForTrivial = (inhibitTrivial) && (!partitionedTriviality);
		}
	else if ((sp1->p != NULL)  // (seq1 is partitioned)
	      && (sp2->p != NULL)) // (seq2 is partitioned)
	 	{
		// $$$ this could be modified to check every partition in seq1 vs every
		// $$$ .. partition in seq2, but would require some restructuring of the
		// $$$ .. code;  we don't expect seq2 to be partitioned in normal use;
		// $$$ .. moreover, this would require *many* more slots in the gapped
		// $$$ .. alignment table
		partitionedTriviality  = identical_partitioned_sequences (seq1, seq2);
		delayedCheckForTrivial = (inhibitTrivial) && (!partitionedTriviality);
		}

	if (!doHspsInBatches)					// empty the list of bounding
		orderBegInc = orderEndDec = NULL;	// ..  alignments

	if ((sp1->p == NULL)  // (seq1 is not partitioned)
	 && (sp2->p == NULL)  // (seq2 is not partitioned)
	 && (identical_sequences (seq1, seq2, scoring, &s)))
		{
		if (doHspsInBatches)
			suicidef ("internal error, attempt to add trivial self-alignment with batches in play");

		freeSlot = anchors->len;
		mp = msp[freeSlot];
		mp->pos1 = mp->pos2 = 0;
		mp->end1 = mp->end2 = seq1->len-1;
	    mp->leftAlign1 =
		    mp->leftAlign2 =
		    mp->rightAlign1 =
		    mp->rightAlign2 = NULL;
		mp->leftSeg1 =
		    mp->leftSeg2 =
		    mp->rightSeg1 =
		    mp->rightSeg2 = NULL;
		mp->firstSeg = NULL;
		save_seg (mp, mp->pos1, mp->pos2, mp->end1, mp->end2);
		insert_align(mp, &orderBegInc, &orderEndDec);
		mp->lastSeg = mp->firstSeg;
		mp->firstSeg->prevSeg = mp->lastSeg->nextSeg = NULL;
		a = mp->align = malloc_or_die ("gapped_extend", sizeof(alignel));
		a->script = edit_script_new();
		edit_script_sub (&a->script, seq1->len);
		a->beg1 = a->beg2 = 1;
		a->end1 = a->end2 = seq1->len;
		a->seq1 = seq1->v;
		a->seq2 = seq2->v;
		if (s < scoreThresh.s) a->s = scoreThresh.s; // so it won't be discarded
		                  else a->s = s;
		a->next = NULL;
		a->isTrivial = true;
		debugSnoopEditScripts_1;
		debugSnoopEditScripts_4;
		}

	else if ((partitionedTriviality)
	      && (sp2->p == NULL)) // (seq2 is not partitioned)
		{
		if (doHspsInBatches)
			suicidef ("internal error, attempt to add trivial self-alignment with batches in play");

		freeSlot = anchors->len;
		p1 = &sp1->p[trivialPartIx];
		s = score_identical_partition_of (seq1, seq2, p1, scoring);
		mp = msp[freeSlot];
		mp->pos1 = p1->sepBefore + 1;
		mp->pos2 = 0;
		mp->end1 = p1->sepAfter  - 1;
		mp->end2 = seq2->len - 1;
		mp->leftAlign1 =
			mp->leftAlign2 =
			mp->rightAlign1 =
			mp->rightAlign2 = NULL;
		mp->leftSeg1 =
			mp->leftSeg2 =
			mp->rightSeg1 =
			mp->rightSeg2 = NULL;
		mp->firstSeg = NULL;
		save_seg (mp, mp->pos1, mp->pos2, mp->end1, mp->end2);
		insert_align(mp, &orderBegInc, &orderEndDec);
		mp->lastSeg = mp->firstSeg;
		mp->firstSeg->prevSeg = mp->lastSeg->nextSeg = NULL;
		a = mp->align = malloc_or_die ("gapped_extend", sizeof(alignel));
		a->script = edit_script_new();
		edit_script_sub (&a->script, seq2->len);
		a->beg1 = p1->sepBefore + 2;
		a->beg2 = 1;
		a->end1 = p1->sepAfter;
		a->end2 = seq2->len;
		a->seq1 = seq1->v;
		a->seq2 = seq2->v;
		if (s < scoreThresh.s) a->s = scoreThresh.s; // so it won't be discarded
						  else a->s = s;
		a->next = NULL;
		a->isTrivial = true;
		debugSnoopEditScripts_2;
		debugSnoopEditScripts_4;
		}

	else if ((partitionedTriviality)
	      && (sp2->p != NULL)) // (seq2 is partitioned)
		{
		if (doHspsInBatches)
			suicidef ("internal error, attempt to add trivial self-alignment with batches in play");

		freeSlot = anchors->len;
		for (partIx=0 ; partIx<sp1->len ; partIx++)
			{
			p1 = &sp1->p[partIx];
			p2 = &sp2->p[partIx];
			s = score_identical_partition (seq1, seq2, p1, p2, scoring);
			mp = msp[freeSlot++];
			mp->pos1 = p1->sepBefore + 1;
			mp->pos2 = p2->sepBefore + 1;
			mp->end1 = p1->sepAfter  - 1;
			mp->end2 = p2->sepAfter  - 1;
			mp->leftAlign1 =
			    mp->leftAlign2 =
			    mp->rightAlign1 =
			    mp->rightAlign2 = NULL;
			mp->leftSeg1 =
			    mp->leftSeg2 =
			    mp->rightSeg1 =
			    mp->rightSeg2 = NULL;
			mp->firstSeg = NULL;
			save_seg (mp, mp->pos1, mp->pos2, mp->end1, mp->end2);
			insert_align(mp, &orderBegInc, &orderEndDec);
			mp->lastSeg = mp->firstSeg;
			mp->firstSeg->prevSeg = mp->lastSeg->nextSeg = NULL;
			a = mp->align = malloc_or_die ("gapped_extend", sizeof(alignel));
			a->script = edit_script_new();
			edit_script_sub (&a->script, p1->sepAfter - (p1->sepBefore+1));
			a->beg1 = p1->sepBefore + 2;
			a->beg2 = p2->sepBefore + 2;
			a->end1 = p1->sepAfter;
			a->end2 = p2->sepAfter;
			a->seq1 = seq1->v;
			a->seq2 = seq2->v;
			if (s < scoreThresh.s) a->s = scoreThresh.s; // so it won't be discarded
			                  else a->s = s;
			a->next = NULL;
			a->isTrivial = true;
			debugSnoopEditScripts_3;
			debugSnoopEditScripts_4;
			}
		}

	// process each batch of anchors

	head = last = NULL;	// (initialize the list of qualifying gapped alignments)

	for (batIx=0 ; batIx<segBatches->len ; batIx++)
		{
		startSegIx = segBatches->batch[batIx].start;
		endSegIx   = segBatches->batch[batIx].end;
		debugSnoopBatches_1;

		if (doHspsInBatches)					// empty the list of bounding
			orderBegInc = orderEndDec = NULL;	// .. alignments, for this batch

		// convert each anchor in this batch to a gapped extension, processing
		// the anchors from high score to low (they've previously been sorted
		// into that order)

		for (i=startSegIx ; i<endSegIx ; i++)
			{
			mp = msp[i];

			debugGappedExtendVerbosity_1;
			debugSnoopAnchorToGapped_1;
			debugsnoopSpecialHsp_1;

			// find the horizontal bounding alignments/segments of this anchor

			if (!msp_left_right (orderBegInc, mp))
				continue;		// (an earlier alignment contains this anchor)

			io.leftAlign  = mp->leftAlign1;
			io.rightAlign = mp->rightAlign1;
			io.leftSeg    = mp->leftSeg1;
			io.rightSeg   = mp->rightSeg1;

			// find the closest alignments ending before or starting after this
			// anchor

			debugSnoopAnchorToGapped_2;
			io.anchor1 = mp->pos1;
			io.anchor2 = mp->pos2;
			dbg_timing_gapped_extend_sub (debugClockAboveBelow);
			get_above_below (&io, orderBegInc, orderEndDec);
			dbg_timing_gapped_extend_add (debugClockAboveBelow);

			// if either sequence is partitioned, figure out the limits of the
			// partition containing this anchor

			if (segBatches->batch[batIx].part != NULL)
				{ p1 = segBatches->batch[batIx].part;  goto set_limits1; }
			else if (sp1->p != NULL)
				{
				p1 = lookup_partition (seq1, io.anchor1);
			set_limits1:
				io.low1  = p1->sepBefore + 1;
				io.high1 = p1->sepAfter;
				}

			if (sp2->p != NULL)
				{
				p2 = lookup_partition (seq2, io.anchor2);
				io.low2  = p2->sepBefore + 1;
				io.high2 = p2->sepAfter;
				}

			// if we have a chore, further restrict the limits to the chore

			if (seq2->choresFile != NULL)
				{
				interval tInt = seq2->chore.targetInterval;
				interval qInt = seq2->chore.queryInterval;
				if (tInt.s > io.low1)  io.low1  = tInt.s;
				if (tInt.e < io.high1) io.high1 = tInt.e;
				if (qInt.s > io.low2)  io.low2  = qInt.s;
				if (qInt.e < io.high2) io.high2 = qInt.e;
				}

			// extend this anchor into a gapped alignment, in both directions

			dbg_timing_gapped_extend_sub (debugClockYdropAlign);
			ydrop_align (&io);
			dbg_timing_gapped_extend_add (debugClockYdropAlign);

			debugSnoopAnchorToGapped_3;
			mp->align = format_alignment (&io, mp);
			mp->pos1  = io.start1;
			mp->pos2  = io.start2;
			mp->end1  = io.stop1;
			mp->end2  = io.stop2;

			if (mp->firstSeg == NULL)	// (the gapped alignment is empty,
				{						//  .. so skip it)
				debugSnoopEditScripts_5;
				continue;
				}

			debugSnoopEditScripts_4;

			// record the alignment's tail and detach the circular pointer

			mp->lastSeg = mp->firstSeg->prevSeg;
			mp->firstSeg->prevSeg = mp->lastSeg->nextSeg = NULL;

			// if this alignment doesn't meet the score threshold, discard it
			// now;  otherwise, save it and use it to bound subsequent gapped
			// extensions;  note that in blastz, the alignment was *always*
			// used as a bound (regardless of whether it met score threshold),
			// and we provide that functionality if allBounds is true (in which
			// case the low-scoring alignments are discarded later)

			if ((!allBounds) && (mp->align->s < scoreThresh.s))
				{
				debugGappedExtendDbgShowIdentity_1;
				debugSnoopBlocks_2;
				free_align_list (mp->align);

				for (bp=mp->firstSeg ; bp!=NULL ; bp=bq)
					{ bq = bp->nextSeg;  free_if_valid ("gapped_extend seg", bp); }

				continue;
				}

			// record the horizontal bounding alignments/segments of this
			// anchor's gapped extension, and insert it into the vertically
			// ordered alignment lists

			dbg_timing_gapped_extend_sub (debugClockLeftRight);
			align_left_right (orderBegInc, mp);
			insert_align (mp, &orderBegInc, &orderEndDec);
			dbg_timing_gapped_extend_add (debugClockLeftRight);

			// if we have a limit on the number of paired bases we'll
			// accept, check for that now;  if we've exceeded the limit,
			// we just stop processing HSPs

			if (maxPairedBases > 0)
				{
				newPairedBases = count_paired_bases (mp);
				pairedBases += newPairedBases;
				//fprintf (stderr, "paired bases: %s of %s\n",
				//                  commatize(newPairedBases),
				//                  commatize(pairedBases));
				if (pairedBases > maxPairedBases)
					{
					if (overlyPairedWarn)
						warn_for_paired_bases_limit (seq2, maxPairedBases,
													 overlyPairedKeep);
					if (!overlyPairedKeep) goto discard_alignments;
					break;  // exit the HSP loop
					}
				}

			debugGappedExtendVerbosity_2;
			}

		// link the high scoring alignments together into a list;  discard the
		// trivial self-alignment if it isn't wanted, and discard any alignments
		// that don't meet the score threshold
		// note that unless allBounds is true (blastz compatibility), there will
		// be no low-scoring alignments to discard here
		// also note that if we have more than one batch, the self-alignment
		// won't be here, so we don't have to worry that we'll discard it too
		// early (but since we like to worry, we do check for that case)

		debugSnoopBlocks_3

		for (mp=orderBegInc; mp!=NULL ; mp=mpNext)
			{
			mpNext = mp->next;
			if (mp->align->s < scoreThresh.s)
				{
				debugGappedExtendDbgShowIdentity_2;
				goto free_up_all;
				}
			if ((inhibitTrivial) && (mp->align->isTrivial))
				goto free_up_all;
			else if (delayedCheckForTrivial)
				{
				aliseg* seg = mp->firstSeg;
				char*   name1, *name2;

				if (mp->lastSeg != seg) goto not_trivial;
				if (seg->type != diagSeg) goto not_trivial;

				if (doHspsInBatches)
					suicidef ("internal error, attempt to discard trivial self-alignment with batches in play");

				if (sp1->p == NULL)
					{
					name1 = seq1->header;
					len1  = seq1->trueLen;
					if ((name1 != NULL) && (name1[0] == '>'))
						name1 = skip_whitespace(name1+1);
					}
				else if (segBatches->batch[batIx].part != NULL)
					{ p1 = segBatches->batch[batIx].part;  goto set_name1; }
				else
					{
					p1 = lookup_partition (seq1, mp->pos1);
				set_name1:
					name1 = &sp1->pool[p1->header];
					len1  = p1->trueLen;
					}

				if (sp2->p == NULL)
					{
					name2 = seq2->header;
					len2  = seq2->trueLen;
					if ((name2 != NULL) && (name2[0] == '>'))
						name2 = skip_whitespace(name2+1);
					}
				else
					{
					p2    = lookup_partition (seq2, mp->pos2);
					name2 = &sp2->pool[p2->header];
					len2  = p2->trueLen;
					}

				debugTriviality_1;

				if (len1 != len2)
					{ debugTriviality_2;  goto not_trivial; }
				if (mp->end1+1 - mp->pos1 != len1)
					{ debugTriviality_3;  goto not_trivial; }
				if (strcmp (name1, name2) != 0)
					{ debugTriviality_4;  goto not_trivial; }

				for (pos1=mp->pos1,pos2=mp->pos2 ; pos1<=mp->end1 ; pos1++,pos2++)
					{
					if (seq1->v[pos1] != seq2->v[pos2])
						{ debugTriviality_5;  goto not_trivial; }
					}

				debugTriviality_6;
				goto free_up_all;
				}
		not_trivial:

			if (head == NULL) head = last       = mp->align;
						 else last = last->next = mp->align;
			debugGappedExtendDbgShowIdentity_3;
			debugSnoopBlocks_4;
			goto free_up_segments;

			// discard an alignment block

		free_up_all:

			debugSnoopBlocks_5;
			free_align_list (mp->align);
			// (fall thru to free_up_segments)

			// discard an alignment block's segments

		free_up_segments:
			for (bp=mp->firstSeg ; bp!=NULL ; bp=bq)
				{ bq = bp->nextSeg;  free_if_valid ("gapped_extend seg", bp); }
			}

		}

	debugSnoopAnchorToGapped_4;

	free_if_valid ("gapped_extend msp[]", msp);

	return head;

	//////////
	// failure exit
	//////////

	// discard all alignments we found

discard_alignments:

	debugSnoopBlocks_3b

	for (mp=orderBegInc; mp!=NULL ; mp=mpNext)
		{
		mpNext = mp->next;

		// discard an alignment block and its segments

		debugSnoopBlocks_5;
		free_align_list (mp->align);

		for (bp=mp->firstSeg ; bp!=NULL ; bp=bq)
			{ bq = bp->nextSeg;  free_if_valid ("gapped_extend seg", bp); }
		}

	debugSnoopAnchorToGapped_5;

	free_if_valid ("gapped_extend msp[]", msp);

	return NULL;
	}

//----------
//
// batched_segments--
//	Separate a list of segments in batches relative to a partitioned sequence.
//	Conceptually we produce one batch per partition.  But if a partition
//	contains no segment, it is just left out of the batch table.
//
//----------
//
// Arguments:
//	segtable*		anchors:	The anchors.  Note that these well be reordered
//								.. by this routine, so that every batch is
//								.. sorted by descreasing score.
//	seqpartition*	sp:			The partitioning to separate by.  This should
//								.. relate to sequence 1.  A special value of
//								.. NULL indicates that partitioning is not,
//								.. wanted, so we should treat the whole segment
//								.. table as a single batch.
//
// Returns:
//	A pointer to the batch table.  This is allocated data, which the caller
//	must eventually dispose of by calling free_segment_batches().
//
//----------

static sbtable* _segBatches = NULL;

static sbtable* batched_segments
   (segtable*		anchors,
	seqpartition*	sp)
	{
	u32				entriesNeeded;
	size_t			bytesNeeded;
	partition*		part;
	segment*		anc;
	unspos			pEnd;
	u32				batIx, partIx, ancIx;
	u32				startIx, endIx;

	// allocate the batch table, or resize it if it's not big enough
	// nota bene: (as of this writing), we don't expect realloc to ever occur,
	//            because the sequence that we're partitioning with is the same
	//            one through the entire program;  that may cahnge in the future

	if (sp == NULL) entriesNeeded = 1;
	           else entriesNeeded = sp->len;
	bytesNeeded = sbtable_bytes (entriesNeeded);

	if (_segBatches == NULL)
		{
		_segBatches = (sbtable*) malloc_or_die ("batched_segments", bytesNeeded);
		_segBatches->size = entriesNeeded;
		}
	else if (entriesNeeded > _segBatches->size)
		{
		_segBatches = (sbtable*) realloc_or_die ("batched_segments", _segBatches, bytesNeeded);
		_segBatches->size = entriesNeeded;
		}

	// if we don't have a partitioned sequence, create one batch convering the
	// entire segment list

	if (sp == NULL)
		{
		_segBatches->batch[0].part  = NULL;
		_segBatches->batch[0].start = 0;
		_segBatches->batch[0].end   = anchors->len;
		_segBatches->len = 1;

		sort_segments (anchors, qSegmentsByDecreasingScore);
		}

	// otherwise, create batches, dividing the segments by partition;  note
	// that we assume each segment is contained within a single partition

	else
		{
		sort_segments (anchors, qSegmentsByPos1);

		batIx = 0;
		ancIx = 0;  anc = &anchors->seg[ancIx];
		for (partIx=0 ; partIx<sp->len ; partIx++)
			{
			if (ancIx >= anchors->len) break;

			part = &sp->p[partIx];
			pEnd = part->sepAfter;
			if (pEnd < anc->pos1 + anc->length) continue;

			_segBatches->batch[batIx].part  = part;
			_segBatches->batch[batIx].start = startIx = ancIx++;  anc++;
			while ((ancIx < anchors->len)
			    && (pEnd >= anc->pos1 + anc->length))
				{ ancIx++;  anc++; }

			_segBatches->batch[batIx].end = endIx = ancIx;
			batIx++;

			sort_some_segments (anchors, startIx, endIx,
			                    qSegmentsByDecreasingScore);
			}

		_segBatches->len = batIx;
		}

	return _segBatches;
	}

//----------
//
// free_segment_batches--
//	Dispose of memory alloceted by batched_segments().
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

void free_segment_batches (void)
	{
	free_if_valid ("free_segment_batches", _segBatches);
	_segBatches = NULL;
	}

//----------
//
// init_from_anchors--
//	Initialize a gapped alignments set from a set of anchors.
//
//----------
//
// Arguments:
//	segtable*	anchors:		The anchors, which we will expand into gapped
//								.. alignment records, intitally as single
//								.. points.
//	u32			numExtraSlots:	The number of extra slots to add (see note 1
//								.. below)
//
// Returns:
//	A pointer to newly allocated data (see description below);  failures result
//	in program fatality.  The caller must eventually dispose of the table, with
//	a call to free().
//
// Data Description:
//	The allocated data is an array of pointers to galign blocks.  It is
//	allocated as a single block from the heap (which is why it can be freed wih
//	a simple call to free).  The first N locations form the array, and the data
//	following that contains the galign blocks.
//
//----------
//
// Notes:
//
// (1)	We allocate extra slots in the gapped alignment table, to allow a
//		trivial self-alignments to be added later.  If the target sequence is
//		*not* partitioned, there will be at most one trivial alignment and we
//		add only one extra slot.  The same is true if the query is *not*
//		partitioned.  If both target and query are partitioned, there can be
//		a trivial alignment for each partition in whichever has fewer
//		partitions, so we add that many extra slots.  We push this determination
//		back to the caller, who must provide the value for numExtraSlots.
//
// (2)	We only copy pos1/pos2 from the anchors.  We ignore the length and
//		don't set end1/end2.  The expectation is that the anchors are single
//		points without length-- upstream code should have reduced them to
//		single points, and initial downstream code only uses pos1/pos2.
//
//----------

static galign** init_from_anchors
   (segtable*	anchors,
	u32			numExtraSlots)
	{
	int			numAnchors, numSlots, ix;
	u32			bytesNeeded, bytesPointers, bytesBlocks;
	galign**	snakes, *snake;
	segment*	anchor;

	numAnchors = anchors->len;
	numSlots   = numAnchors + numExtraSlots;

	// allocate;  note that we allocate extra slots as per note (1) above
	// $$$ Why is zalloc used here instead of malloc?  If the table is large,
	// $$$ .. zalloc will spec significant time zeroing it, but why is that
	// $$$ .. necessary?  This sets the anchor's length to zero, but would it
	// $$$ .. be more efficient to set that explicitely, so that we don't waste
	// $$$ .. time clearing out the extra slots which are often not used?

	bytesPointers =  round_up_16 (numSlots * sizeof(galign*));
	bytesBlocks   =  round_up_16 (numSlots * sizeof(galign));
	bytesNeeded   =  bytesPointers;
	bytesNeeded   += bytesBlocks;  if (bytesNeeded < bytesBlocks) goto overflow;

	snakes = (galign**) zalloc_or_die ("init_from_anchors", bytesNeeded);
	gapped_extend_count_stat (zallocCallsA);
	gapped_extend_add_stat   (zallocTotalA, bytesNeeded);

	// hook up pointers;  this links the array of N pointers to the actual
	// blocks

	snakes[0] = (galign*) (((char*) snakes) + bytesPointers);

	for (ix=1 ; ix<numSlots ; ix++)
		snakes[ix] = snakes[ix-1]+1;

	// initialize the pos1/pos2 values of the first N blocks, corresponding
	// to the anchors

	for (ix=0 ; ix<numAnchors ; ix++)
		{
		snake  = snakes[ix];
		anchor = &anchors->seg[ix];
		snake->pos1 = anchor->pos1;
		snake->pos2 = anchor->pos2;
		}

	return snakes;

overflow:
	suicidef ("in init_from_anchors(), structure size exceeds 2^32 (%u+%u)",
	          bytesPointers, bytesBlocks);

	return NULL; // (doesn't get here)
	}

//----------
//
// identical_sequences--
//	Determine if two sequences are identical for practical purposes.
//
//----------
//
// Arguments:
//	seq*		seq1:		One sequence
//	seq*		seq2:		The other sequence.
//	scoreset*	scoring:	The scoring scheme to use.  This may be NULL if the
//							.. caller isn't interested in the score.
//	score*		s:			Place to return the score.  This may be NULL.
//
// Returns:
//  true if the sequences are identical, but allowing differences in upper/lower
//	case;  false otherwise.
//
//----------
//
// Notes:
//
// (1)	(this note also applies to score_identical_partition and
//		score_identical_partition_of)
//		When computing the score below we have to consider the possibility
//		of overflow.  For example, if matches are worth 100, then a sequence
//		of length 21.5M will exceed 2^31 and overflow (assuming scores are
//		32-bit ints).  We prevent this by checking whether adding a positve
//		score will overflow.  However, for negative match scores (e.g. N vs N),
//		we simply subtract.  Thus we are assuming that that we won't have so
//		many negative scores that we underflow.
//
//----------

static int identical_sequences
   (seq*			seq1,
	seq*			seq2,
	scoreset*		scoring,
	score*			_s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	u8*				a,   *b;
	u8				aNuc, bNuc;
	score			s, sub;
	unspos			ix;

	if ((sp1->p != NULL) || (sp2->p != NULL))	// one or the
		return false;							// .. other is partitioned

	if ((seq1->fileType == seq_type_qdna) != (seq2->fileType == seq_type_qdna))
		return false;

	if (seq1->len != seq2->len)
		return false;

	if (seq1->revCompFlags != seq2->revCompFlags)
		return false;

	a = seq1->v;
	b = seq2->v;
	s = 0;

	for (ix=0 ; ix<seq1->len ; ix++)
		{
		aNuc = (u8) dna_toupper(a[ix]);
		bNuc = (u8) dna_toupper(b[ix]);
		if (aNuc != bNuc) return false;
		if (scoring == NULL) continue;

		sub = scoring->sub[aNuc][bNuc];		// (see note above about overflow)
		if (s == bestPossibleScore)
			;
		else if ((sub <= 0) || (s < bestPossibleScore - sub))
			s += sub;
		else
			s =  bestPossibleScore;
		}

	if (_s != NULL) *(_s) = s;
	return true;
	}

//----------
//
// identical_partitioned_sequences--
//	Determine if two partitioned sequences are identical for practical purposes.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence.
//	seq*	seq2:	The other sequence.
//
// Returns:
//  true if the sequences are identical across their partitions, but allowing
//	differences in upper/lower case;  false otherwise.
//
//----------

static int identical_partitioned_sequences
   (seq*			seq1,
	seq*			seq2)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		p1, *p2;
	u32				partIx;
	u8*				a,   *b;
	u8				aNuc, bNuc;
	unspos			len1, len2, ix;

	if ((sp1->p == NULL) || (sp2->p == NULL))	// one or the other
		return false;							// .. is *not* partitioned

	if ((seq1->fileType == seq_type_qdna) != (seq2->fileType == seq_type_qdna))
		return false;

	if (sp1->len != sp2->len)
		return false;

	if (seq1->revCompFlags != seq2->revCompFlags)
		return false;

	for (partIx=0 ; partIx<sp1->len ; partIx++)
		{
		p1 = &sp1->p[partIx];
		p2 = &sp2->p[partIx];
		len1 = p1->sepAfter - (p1->sepBefore+1);
		len2 = p2->sepAfter - (p2->sepBefore+1);
		if (len1 != len2) return false;

		a = seq1->v + p1->sepBefore+1;
		b = seq2->v + p2->sepBefore+1;
		for (ix=0 ; ix<len1 ; ix++)
			{
			aNuc = (u8) dna_toupper(a[ix]);
			bNuc = (u8) dna_toupper(b[ix]);
			if (aNuc != bNuc) return false;
			}
		}

	return true;
	}

//----------
//
// identical_partition_of_sequence--
//	Determine if a sequence is identical (for practical purposes) to any
//	partition of another sequence.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence (partitioned).
//	seq*	seq2:	The other sequence (not partitioned).
//
// Returns:
//	The number of a partition of seq1 that is identical to the entirety of
//	seq2;  -1 if there is no matching partition.
//
//----------
//
//	Notes:
//
//	(1)	In many use cases, this function is useless.  For example, if we are
//		mapping a million short reads (seq2) to a reference genome with tens
//		of chromosomes (seq1), there is no chance that we will get a match.
//		Worse, if the reference is a 100 thousand medium sized coding exons,
//		the time spent search this partition table for a length match, even if
//		no lengths match, can be significant.
//
//		To circumvent this, we #define cache_partition_lengths.  This enables
//		code that computes the range of partition lengths on the first call.
//		Any sequences outside that length range are quickly reported as "no
//		match".
//		
//----------

#define cache_partition_lengths // see note (1)


static int identical_partition_of_sequence
   (seq*			seq1,
	seq*			seq2)
	{
#ifdef cache_partition_lengths // see note (1)
	static seqpartition*	cachedSp1    = NULL;
	static unspos			cachedMinLen = ((unspos) -1);
	static unspos			cachedMaxLen = ((unspos) -1);
#endif // cache_partition_lengths
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		p1;
	u32				partIx;
	u8*				a,   *b;
	u8				aNuc, bNuc;
	unspos			len1, len2, ix;
	int				isMatch;

	if (sp1->p == NULL) return -1;		// seq1 is *not* partitioned
	if (sp2->p != NULL) return -1;		// seq2 is partitioned

	if ((seq1->fileType == seq_type_qdna) != (seq2->fileType == seq_type_qdna))
		return -1;

	if (seq1->revCompFlags != seq2->revCompFlags)
		return -1;

	// if this is a previously unseen partitioning, prescan to find the min and
	// max lengths

#ifdef cache_partition_lengths
	if (sp1 != cachedSp1)
		{
		cachedSp1 = sp1;

		partIx = 0;
		p1 = &sp1->p[partIx];
		len1 = p1->sepAfter - (p1->sepBefore+1);
		cachedMinLen = cachedMaxLen = len1;

		for (partIx=1 ; partIx<sp1->len ; partIx++)
			{
			p1 = &sp1->p[partIx];
			len1 = p1->sepAfter - (p1->sepBefore+1);
			if      (len1 < cachedMinLen) cachedMinLen = len1;
			else if (len1 > cachedMaxLen) cachedMaxLen = len1;
			}

		//fprintf (stderr, "partition lengths: " unsposFmt ".." unsposFmt "\n",
		//                 cachedMinLen, cachedMaxLen);
		}
#endif // cache_partition_lengths

	// if the second sequence is shorter or longer than all partitions, there
	// can be no match

	len2 = seq2->len;

#ifdef cache_partition_lengths
	if (sp1 == cachedSp1)
		{
		if ((len2 < cachedMinLen) || (len2 > cachedMaxLen))
			{
			// fprintf (stderr, "length out of range: " unsposFmt "\n", len2);
			return -1;
			}
		}
#endif // cache_partition_lengths

	// otherwise, scan all partitions for a match

	for (partIx=0 ; partIx<sp1->len ; partIx++)
		{
		p1 = &sp1->p[partIx];
		len1 = p1->sepAfter - (p1->sepBefore+1);
		if (len1 != len2) continue;

		a = seq1->v + p1->sepBefore+1;
		b = seq2->v;
		isMatch = true;
		for (ix=0 ; ix<len1 ; ix++)
			{
			aNuc = (u8) dna_toupper(a[ix]);
			bNuc = (u8) dna_toupper(b[ix]);
			if (aNuc != bNuc) { isMatch = false;  break; }
			}
		if (isMatch) return partIx;
		}

	return -1; // (no matching partition was found)
	}

//----------
//
// score_identical_partition--
//	Compute the score of two partitions, which are assumed to have indentical
//	length.
//
//----------
//
// Arguments:
//	seq*		seq1:		The sequence containing p1.
//	seq*		seq2:		The sequence containing p2.
//	partition*	p1:			One partition.
//	partition*	p2:			The other partition.
//	scoreset*	scoring:	The scoring scheme to use.
//
// Returns:
//	the score.
//
//----------
//
// Notes:
//
// (1)	We assume, WITHOUT CHECKING, that the partitions have the same length.
//
// (2)	(see note 1 in identical_sequences, regarding overflow)
//
//----------

static score score_identical_partition
   (seq*			seq1,
	seq*			seq2,
	partition*		p1,
	partition*		p2,
	scoreset*		scoring)
	{
	u8*				a,   *b;
	u8				aNuc, bNuc;
	score			s, sub;
	unspos			len, ix;

	len = p1->sepAfter - (p1->sepBefore+1);

	a = seq1->v + p1->sepBefore+1;
	b = seq2->v + p2->sepBefore+1;
	s = 0;
	for (ix=0 ; ix<len ; ix++)
		{
		aNuc = (u8) dna_toupper(a[ix]);
		bNuc = (u8) dna_toupper(b[ix]);
		if (aNuc != bNuc) return false;

		sub = scoring->sub[aNuc][bNuc];		// (see note above about overflow)
		if (s == bestPossibleScore)
			;
		else if ((sub <= 0) || (s < bestPossibleScore - sub))
			s += sub;
		else
			s =  bestPossibleScore;
		}

	return s;
	}

//----------
//
// score_identical_partition_of--
//	Compute the score of a partition versus a sequence, which are assumed to
//	have indentical length.
//
//----------
//
// Arguments:
//	seq*		seq1:		The sequence containing p1.
//	seq*		seq2:		The other sequence (not partitioned).
//	partition*	p1:			The partition.
//	scoreset*	scoring:	The scoring scheme to use.
//
// Returns:
//	the score.
//
//----------
//
// Notes:
//
// (1)	We assume, WITHOUT CHECKING, that the partition and sequence have the
//		same length.
//
// (2)	(see note 1 in identical_sequences, regarding overflow)
//
//----------

static score score_identical_partition_of
   (seq*			seq1,
	seq*			seq2,
	partition*		p1,
	scoreset*		scoring)
	{
	u8*				a,   *b;
	u8				aNuc, bNuc;
	score			s, sub;
	unspos			len, ix;

	len = p1->sepAfter - (p1->sepBefore+1);

	a = seq1->v + p1->sepBefore+1;
	b = seq2->v;
	s = 0;
	for (ix=0 ; ix<len ; ix++)
		{
		aNuc = (u8) dna_toupper(a[ix]);
		bNuc = (u8) dna_toupper(b[ix]);
		if (aNuc != bNuc) return false;

		sub = scoring->sub[aNuc][bNuc];		// (see note above about overflow)
		if (s == bestPossibleScore)
			;
		else if ((sub <= 0) || (s < bestPossibleScore - sub))
			s += sub;
		else
			s =  bestPossibleScore;
		}

	return s;
	}

//----------
//
// new_traceback--
//	Allocate some traceback data.
//
//----------
//
// Arguments:
//	u32		size:	The maximum number of bytes to allocate.
//
// Returns:
//	A pointer to a newly allocated traceback data;  failures result in program
//	fatality.  The caller must eventually dispose of the table, with a call to
//	free_traceback().
//
//----------

tback* new_traceback
   (u32		size)
	{
	tback*	tb;
	int		cells;

	// sanity check

	if (size < sizeof(tback))
		suicidef ("in new_traceback(), size can't be %u", size);

	// allocate

	cells = 1 + (size - sizeof(tback)) / sizeof(tb->space[0]);
	tb = (tback*) malloc_or_die ("new_traceback", size);

	// initialize

	tb->size = cells;

	return tb;
	}

//----------
//
// free_traceback--
//	Deallocate traceback data.
//
//----------
//
// Arguments:
//	tback*	tb:	The traceback data to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_traceback
   (tback*	tb)
	{
	free_if_valid ("free_traceback", tb);
	}

//----------
//
// ydrop_align--
//	The Y-drop variant of dynamic programming is applied to create a gapped
//	alignment by extending in both directions from an "anchor point".
//
//----------
//
// Arguments:
//	alignio*	io:		The collection of input arguments, and a place to store
//						.. results.  The input arguments describe the alignment
//						.. problem to be solved.  See the definition of the type
//						.. alignio for more details.
//
// Returns:
//  nothing;  actual result values are in io: s, start1, start2, stop1, stop2
//	and script.
//
//----------

#ifndef snoopSubprobs
#define snoopSubprobsA_V ;
#define snoopSubprobsA_1 ;
#define snoopSubprobsA_2 ;
#define snoopSubprobsA_3 ;
#endif // not snoopSubprobs

#ifdef snoopSubprobs

#define snoopSubprobsA_V                                                     \
	partition*	part;                                                        \
	char*		name1, *name2

#define snoopSubprobsA_1                                                     \
	if (snoopSubprobsSeq1->partition.p == NULL)                              \
		{                                                                    \
		name1 = (snoopSubprobsSeq1->useFullNames)? snoopSubprobsSeq1->header \
		                                         : snoopSubprobsSeq1->shortHeader; \
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";              \
		}                                                                    \
	else                                                                     \
		{                                                                    \
		part  = lookup_partition (snoopSubprobsSeq1, anchor1);               \
		name1 = &snoopSubprobsSeq1->partition.pool[part->header];            \
		}                                                                    \
                                                                             \
	if (snoopSubprobsSeq2->partition.p == NULL)                              \
		{                                                                    \
		name2 = (snoopSubprobsSeq2->useFullNames)? snoopSubprobsSeq2->header \
		                                         : snoopSubprobsSeq2->shortHeader; \
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";              \
		}                                                                    \
	else                                                                     \
		{                                                                    \
		part  = lookup_partition (snoopSubprobsSeq2, anchor2);               \
		name2 = &snoopSubprobsSeq2->partition.pool[part->header];            \
		}                                                                    \
                                                                             \
	fprintf (stderr, "align:\tbck\t%s\t" unsposFmt "\t" unsposFmt            \
	                         "\t%s\t" unsposFmt "\t" unsposFmt,              \
	                 name1, io->len1-anchor1-2, (anchor1+1)-io->low1,        \
	                 name2, io->len2-anchor2-2, (anchor2+1)-io->low2);

#define snoopSubprobsA_2                                                     \
	if (snoopSubprobsSeq1->partition.p == NULL)                              \
		{                                                                    \
		name1 = (snoopSubprobsSeq1->useFullNames)? snoopSubprobsSeq1->header \
		                                         : snoopSubprobsSeq1->shortHeader; \
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";              \
		}                                                                    \
	else                                                                     \
		{                                                                    \
		part  = lookup_partition (snoopSubprobsSeq1, anchor1);               \
		name1 = &snoopSubprobsSeq1->partition.pool[part->header];            \
		}                                                                    \
                                                                             \
	if (snoopSubprobsSeq2->partition.p == NULL)                              \
		{                                                                    \
		name2 = (snoopSubprobsSeq2->useFullNames)? snoopSubprobsSeq2->header \
		                                         : snoopSubprobsSeq2->shortHeader; \
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";              \
		}                                                                    \
	else                                                                     \
		{                                                                    \
		part  = lookup_partition (snoopSubprobsSeq2, anchor2);               \
		name2 = &snoopSubprobsSeq2->partition.pool[part->header];            \
		}                                                                    \
                                                                             \
	fprintf (stderr, "\nalign:\tfwd\t%s\t" unsposFmt "\t" unsposFmt          \
	                           "\t%s\t" unsposFmt "\t" unsposFmt,            \
	                 name1, anchor1, io->high1-(anchor1+1),                  \
	                 name2, anchor2, io->high2-(anchor2+1));

#define snoopSubprobsA_3                                                     \
	fprintf (stderr, "\n");

#endif // snoopSubprobs


static void ydrop_align
   (alignio*	io)
	{
	unspos		anchor1, anchor2;
	unspos		end1, end2;
	score		scoreLeft, scoreRight;
	editscript* script, *scriptRight;
#ifdef snoopAlgorithm
	int			snoop = true;
#endif // snoopAlgorithm
	snoopSubprobsA_V;

	if (io == NULL)
		suicide ("ydrop_align() called with NULL pointer.");

#ifdef snoopAlignioInput
	dump_alignio_input (stderr, io);
#endif // snoopAlignioInput

	gapped_extend_count_stat (numAnchorsExtended);

	anchor1 = io->anchor1;
	anchor2 = io->anchor2;

#ifdef snoopAlgorithm
#if ((defined debugPosA1) && (defined debugPosB1))
	snoop = ((anchor1 == debugPosA1) && (anchor2 == debugPosA2))
	     || ((anchor1 == debugPosB1) && (anchor2 == debugPosB2));
#elif (defined debugPosA1)
	snoop = (anchor1 == debugPosA1) && (anchor2 == debugPosA2);
#elif (defined debugPosB1)
	snoop = (anchor1 == debugPosB1) && (anchor2 == debugPosB2);
#endif // debugPosA1,debugPosB1
#endif // snoopAlgorithm

#ifdef snoopAlgorithm
	if (snoop)
		{
		fprintf (stderr, "gapE  = " scoreFmtSimple "\n", io->scoring->gapExtend);
		fprintf (stderr, "gapOE = " scoreFmtSimple "\n", io->scoring->gapOpen + io->scoring->gapExtend);
		fprintf (stderr, "yDrop = " scoreFmtSimple "\n", io->yDrop);
		dump_score_set (stderr, io->scoring, (u8*)"ACGTacgtNn", (u8*)"ACGTacgtNn");
		}
#endif // snoopAlgorithm

	snoopSubprobsA_1;

	script = edit_script_new();
	dbg_timing_gapped_extend_sub (debugClockYdropOneSidedAlign);
	scoreLeft = ydrop_one_sided_align (io, /*reverse*/ true,
	                                   io->rev1 + io->len1 - anchor1 - 2,
	                                   io->rev2 + io->len2 - anchor2 - 2,
	                                   (anchor1+1) - io->low1,
	                                   (anchor2+1) - io->low2,
	                                   io->trimToPeak,
	                                   &script, &end1, &end2);
	dbg_timing_gapped_extend_add (debugClockYdropOneSidedAlign);
	io->start1 = anchor1 + 1 - end1;
	io->start2 = anchor2 + 1 - end2;

	snoopSubprobsA_2;

	scriptRight = edit_script_new();
	dbg_timing_gapped_extend_sub (debugClockYdropOneSidedAlign);
	scoreRight = ydrop_one_sided_align (io, /*reverse*/ false,
	                                    io->seq1 + anchor1,
	                                    io->seq2 + anchor2,
	                                    io->high1 - (anchor1+1),
	                                    io->high2 - (anchor2+1),
	                                    io->trimToPeak,
	                                    &scriptRight, &end1, &end2);
	dbg_timing_gapped_extend_add (debugClockYdropOneSidedAlign);
	io->stop1 = anchor1 + end1;
	io->stop2 = anchor2 + end2;

	snoopSubprobsA_3;

#ifdef snoopAlgorithm
	if (snoop)
		{
		fprintf          (stderr, "==== left edit script:\n");
		dump_edit_script (stderr, script);
		fprintf          (stderr, "==== right edit script:\n");
		dump_edit_script (stderr, scriptRight);
		}
#endif // snoopAlgorithm

	edit_script_reverse (scriptRight);
	edit_script_append  (&script, scriptRight);
	free_if_valid       ("ydrop_align, edit script", scriptRight);

	io->s      = scoreRight + scoreLeft;
	io->script = script;

#ifdef snoopAlgorithm
	if (snoop)
		{
		fprintf          (stderr, "==== combined edit script:\n");
		dump_edit_script (stderr, script);
		fprintf          (stderr, "\n");
		}
#endif // snoopAlgorithm

#ifdef snoopAlignioOutput
	dump_alignio_output (stderr, io);
#endif // snoopAlignioOutput
	}

//----------
//
// ydrop_one_sided_align--
//	Find an gapped alignment in one direction from a starting point.
//
//----------
//
// Arguments:
//	alignio*	io:			Other input arguments for the alignment problem
//							.. being solved.  Of interest to this routine are
//							.. the scoring parameters (scoring and yDrop) and
//							.. the rightSeg/leftSeg stuff, which is used to
//							.. track sweep row bounds caused by previous
//							.. alignments.
//	int			reversed:	true  => search to the lower-left
//							false => search to the upper-right
//	u8*			A:			Query sequence (conceptually vertical)
//							.. (indexed by 1..M, see note below)
//	u8*			B:			Subject sequence (conceptually horizontal)
//							.. (indexed by 1..N, see note below)
//	unspos		M:			Length of query sequence.
//	unspos		N:			Length of subject sequence.
//	int			trimToPeak:	true  => always trim end of alignment to score peak.
//							false => don't trim if extension runs into end of
//							         .. sequence (see note 12 below).
//	editscript** script:	Place to return the resulting edit script.
//	unspos*		end1:		Place to return the alignment end, in A.
//	unspos*		end2:		Place to return the alignment end, in B.
//
// Returns:
//  The score of the alignment.
//
//----------
//
// Notes:
//	- A and B are indexed starting from 1, so the pointers must point to one
//	  character before the start of the string.  That character is ignored.
//
//	- Whether reversed is true or false, the DP matrix is conceptually scanned
//	  from bottom-up, and from left-to-right (see diagram below).  What changes
//	  is the equivalence between DP row,col and sequence row,col.  
//
//			           ------ right ----->
//			      +----+------------------+
//			row M |A[M]|                  |  ^
//			 ...  |... |                  |  |
//			row 1 |A[1]|    DP Matrix     |  up
//			row 0 |    |                  |  |
//			      |    +------------------+  |
//			      |          B[1] ... B[N]|
//			      +-----------------------+
//			            col  col  ...  col
//			             0    1         N
//
//----------
//
// Implementation Notes:
//
// (1)	We manage the external bounds of the current sweep row by keeping track
//		of the alignment and the specific segments that constrain the feasible
//		region on the left (leftAlign,leftSeg) and right (rightAlign,rightSeg),
//		if they exist.  These are initially set to the constraints of the anchor
//		point, and information is maintained to allow efficient updating as the
//		sweep row advances.
//
//		For a given row, the values L and R bracket the legal DP columns to
//		compute.  Neither column L nor R is a valid column in that row, being
//		part of a previous alignment or off the left or right edge of the
//		matrix.  Similarly, LY and RY bracket the intersection of the y-drop
//		feasible region and the left/right bounds.
//
// (2)  The new alignment should not be allowed to intersect any previously
//		found alignments, so we have to keep track of "active segments"-- gap-
//		free segments of earlier alignments that intersect the sweep row within
//		the feasible region.  These are kept in a list that supports insertion
//		and deletion when the sweep row reaches either end of an earlier
//		alignment.
//
//		As we advance the sweep row, update_active_segs identifies all columns
//		(within the feasible range) that intersect an existing alignment.  It
//		marks the corresponding DP cells as 'masked' (the mask indicator value
//		is the current row number, so that we don't have to erase them for the
//		next sweep row).  When the sweep encounters a masked cell, it refuses
//		to step to it, and sets all three scores to -infinity so that no step
//		will be taken from it.
//
// (3)	Memory for traceback cells is provided by the caller (in io->tb).  That
//		memory is carved into pieces matching the computed rows/slices of the DP
//		matrix.  This routine allocates an array of indexes to the start of
//		each row (with the index offset to the left so that columns index the
//		correct cells).
//
//		The memory for this index array is allocated here, and lives across
//		calls as a static variable.  We allocate it in multiples of 512K in the
//		hope that this improves performance on some platforms (specifically, to
//		avoid expensive mmap/munmap activity in linux with gnu malloc);  it is
//		not clear whether this helps.
//
//		We use an array of indexes, rather than an array of pointers, to reduce
//		the memory requirement on 64-bit machines.  Pointers on such machines
//		would require 8 bytes per entry, whereas indexes require only 4 bytes.
//
//		The BLASTZ version of this function allocated a (potentially) very large
//		pointer array.  It allocated at least M+1 entries, one for each possible
//		row.  For long sequences such as human chromosome 1 (~ 250 MB), this
//		required 1GB of memory if the anchor was close to either end of the
//		sequence.  Typical alignments are very unlikely to use that many rows,
//		especially if B is much shorter than A.  So in this implementation, the
//		pointer array is initially allocated as 1% of the the sequence length,
//		but not more than the length of B, and rounded up to a multiple of 512K
//		bytes.  It is "sort of" expanded at a rate of 6% whenever the DP scan
//		requires more rows (in practice, it just allocates the next chunk of
//		512K bytes).  We don't want to have to reallocate very often, because
//		realloc may have to copy the entire contents to the new buffer (even if
//		we don't need it to, as when we are starting a new alignment;  too bad
//		malloc/realloc doesn't provide a routine that lets us say "I need more,
//		but I don't need you to preserve the contents).
//
// (4)	We only use DP memory corresponding to one row/slice of the matrix.  At
//		any time, we need information about cells in one row and the row above
//		it, but we only need both rows at the scan point.  We use three local
//		variables to represent the extra cell there.
//
//		One pointer (dp) scans the cells for the lower row (the row being read),
//		another (dq) scans the upper row.  These pointers are the same when
//		the two rows have the same left bound.  When the upper row has a
//		different bound (which is necessarily further to the right), dq will
//		lag behind dp.
//
// (5)	The inner loop scans the current row from left to the right, starting
//		with column LY and ending before column RY.  Note that column LY is an
//		invalid column (it is column zero, or it is on a segment, or it is
//		beyond ydrop);  this corresponds to the traditional full DP algorithm
//		computing column zero even though column one is the first character of
//		B.
//
//		The following are true for each pass through the inner loop (see note
//		(4) for information about dp and dq):
//
//			We enter the loop each time with
//				dp points to cell to read  DP[*][col]
//				dq points to cell to write DP[*][col]
//				C[row-1][col] in dp->CC
//				C[row  ][col] in c (proposed, not final)
//				D[row  ][col] in dp->DD
//				I[row  ][col] in i
//				A[row  ][*  ] reflected in sub[*]
//				B[col+1]      in *b
//			And during the loop, we will write
//				C[row  ][col  ] to dq->CC, and check/update bestScore
//				C[row  ][col+1] to c (proposed)
//				I[row  ][col+1] to i
//				D[row+1][col  ] to dq->DD
//				link into C[row][col] to *tbp
//
//		It is important to realize that the stored values in one cell are for
//		the C node in one row and the D node for the row above it.
//
// (6)	At the beginning of each row, we grab a pointer (sub) to the row of the
//		substitution scores matrix for the character A[row].  At the end of each
//		pass through the inner loop, we look up the score for the character
//		B[col+1].  Thus we are computing the value of a substitution into
//		node C[row,col+1].
//
// (7)	The trace-back information, i.e., to trace backwards along an optimal
//		path, is packed into one byte per dynamic-programming grid point.  At
//		each grid point there are three nodes: C which is entered by one
//		diagonal edge,  D which is entered by two vertical edges (from a D node
//		and a C node), and I which is entered by two horizontal edges (from an
//		I node and a C node).  See ref [1] for further details.
//
//		The right-most two bits hold cFromC (0), cFromI (1) or cFromD (2),
//		telling how this cell's C node gets the maximum score from an entering
//		edge.
//
//		iExtend (4) and dExtend (8) relate to edges *exiting* this cell to the I
//		or D node in the next column or row, respectively.  iExtend is set iff
//		it is better to take a horizontal edge from this cell's I node (a gap
//		extend) than a horizontal edge from the C node (a gap open).  dExtend
//		has a similar meaning for vertical edges.
//
// (8)	The loop to compute the first row of the DP matrix is limited by col<=N,
//		but instead this ought to be col<R.  However, that might necessitate
//		some code to allow for any overhang on the next row (similar to the code
//		for note (9)).  Further, the other limit, c>=-yDrop, will keep us from
//		going very far (with typical values, this condition will stop us after
//		around 300 columns).
//
// (9)	When we stop computing one row at its right bound, we may still need
//		to prolong computation in the row to support an overhang on the row
//		above.  We do so by allowing a series of insertions.
//
//		The following are true for each pass through the row prolongation loop:
//			dq points to cell to write DP[*][col]
//			I[row][col] in i
//
// (10) We use the loop body below left even though the one below right seems
//		more mathematically correct.
//
//			dq->CC = c;						dq->CC = c;
//			(dq++)->DD = c - gapOE;			c -= gapE;
//			c -= gapE;						(dq++)->DD = c;
//			*(tbp++) = cFromI;				*(tbp++) = cFromI;
//
//		However, using the one on the right caused problems aligning
//
//			A=GTTTTTTTTTTTTTTTCTT
//			B=TTTTTTTTTTTTTTTTTTCTT
//
//		It preferred the alignment on the left to the one on the right, because
//		it failed to charge a second gap open penalty.
//
//			---GTTTTTTTTTTTTTTTCTT			--GTTTTTTTTTTTTTTTCTT
//			TTT-TTTTTTTTTTTTTTTCTT			TTTTTTTTTTTTTTTTTTCTT
//
// (11) Blastz had a potential problem with the feasibility bounds LY and RY.
//		update_LR_bounds occasionally will create the condition RY < LY.  This
//		means that no columns of this dp matrix row will be computed, because
//		the col<RY clause in the column loop will halt the loop, and the post-
//		loop check for LY>=RY will halt the routine.  However, the calculation
//		of tbNeeded is suspect, since it could conceivably be negative.  I don't
//		believe it ever goes negative in practice (because I don't think RY-LY
//		is ever less than -2), but for safety sake we clip it at zero.
//
// (12) Correction for short read "mismatch shadow".  The original intent of
//		this routine is that the two sequences would be long (e.g. chromosomes).
//		In that context, it made sense to trim any negatively-scoring suffix
//		from the end of an alignment.  However, if one of the sequences is a
//		short read (e.g. 50 bases), we are usually interested in aligning the
//		entire read if we can.  If we don't, we will be creating a bias against
//		reporting mismatches near the end of the reads.
//
//		Suppose scores are +1 for a match and -7 for a mismatch.  Then any
//		mismatch within 7 bases of the end will create a negatively-scoring
//		suffix, even if all the bases between it and the end are matches.  If
//		this is trimmed, that mismatch is never reported as part of an
//		alignment.
//
//		The argument trimToPeak is added to control this behavior.  If
//		trimToPeak is true, we *always* trim the alignment end back to the
//		maximum scoring location.  If trimToPeak is false, we don't trim if we
//		happen to encounter the end of the sequence.
//
// (13) The BLASTZ loop to compute the first row of dp matrix was equivalent to
//		this:
//
//			dq = dynProg->p;
//			dq->CC = 0;
//			c = (dq++)->DD = -gapOE;
//
//			for (col=1 ; (col<=N)&&(c>=-yDrop) ; col++)
//				{
//				dq->CC = c;
//				(dq++)->DD = c - gapOE;
//				c -= gapE;
//				}
//
//		The (c>=-yDrop) part of this test is not correct, because it is testing
//		the value of the D entry in the cell, instead of the C entry.
//
// (14)	Bounding segment calculations are initially done w.r.t. the full DP
//		matrix with forward orientation of the sequences.  When we are doing a
//		reversed alignment, L and R (the left and right bounds) must be
//		swapped.  This swap appeared to be performed wrong in BLASTZ, and has
//		been changed, as shown in this table:
//			
//			left    right  | BLASTZ  | LASTZ     |
//			bound?  bound? | L'   R' | L'    R'  |
//			---------------+---------+-----------+
//			  no      no   |  L   R  |  L    R   |
//			  no      yes  | -R   R  | -R+1  N+1 |
//			  yes     no   |  L  -L  |  0   -L-1 |
//			  yes     yes  | -R  -L  | -R+1 -L-1 |
//			---------------+---------+-----------+
//
//----------

//=== macros for ydrop_one_sided_align ===

#define cFromC  0		// (c bit is no bits) see note (7)
#define cFromI  1
#define cFromD  2
#define iExtend 4
#define dExtend 8
#define cidBits (cFromC | cFromI | cFromD)

#define op_string(op) ((op==cFromC)? "SUB" : (op==cFromI)? "INS" : (op==cFromD)? "DEL" : "???")


#define prune                                                                \
	c = dp->CC + sub[*(b++)];         /* propose C[row][col+1]         */    \
	if (col == LY)                    /* (we're at first col in sweep) */    \
		LY++;                                                                \
	else                                                                     \
		{                             /* 'set' I[row][col+1]           */    \
		i = dq->DD = dq->CC = negInf; /* .. and set D[row+1][col]      */    \
		dq++;                         /* .. and set C[row+1][col]      */    \
		}                                                                    \
	dp++;                                                                    \
	*(tbp++) = 0;


//=== stuff for link_to_string ===

#define include_link_to_string                                               \
static char* link_to_string (int link)                                       \
	{                                                                        \
	static char  str[100];                                                   \
	char* s = str;                                                           \
	*(s++) = '{';                                                            \
	if ((link & cidBits) == 0) { *(s++) = 'c';  *(s++) = ','; }              \
	if ((link & cFromI)  != 0) { *(s++) = 'i';  *(s++) = ','; }              \
	if ((link & cFromD)  != 0) { *(s++) = 'd';  *(s++) = ','; }              \
	if ((link & iExtend) != 0) { *(s++) = 'i';  *(s++) = 'i';  *(s++) = ','; } \
	if ((link & dExtend) != 0) { *(s++) = 'd';  *(s++) = 'd';  *(s++) = ','; } \
	if (s > str+1) s--;                                                      \
	*(s++) = '}';                                                            \
	*s = 0;                                                                  \
	return str;                                                              \
	}


//=== stuff for snoopAlgorithm ===

#ifndef snoopAlgorithm
#define snoopAlgorithm_1    ;
#define snoopAlgorithm_2    ;
#define snoopAlgorithm_3    ;
#define snoopAlgorithm_4A   ;
#define snoopAlgorithm_4B   ;
#define snoopAlgorithm_4C   ;
#define snoopAlgorithm_5    ;
#define snoopAlgorithm_5A   ;
#define snoopAlgorithm_6    ;
#define snoopAlgorithm_7A   ;
#define snoopAlgorithm_7B   ;
#define snoopAlgorithm_7C   ;
#define snoopAlgorithm_7D   ;
#define snoopAlgorithm_7E   ;
#define snoopAlgorithm_8    ;
#endif // not snoopAlgorithm

#ifdef snoopAlgorithm

static char ccSnoop[100];
static char ddSnoop[100];
static char iiSnoop[100];

static char* relative_to_infinity (score s)
	{
	static char  str[100];
	if      (s <  negInf)      sprintf (str, "-inf"  scoreFmtSimple, s-negInf);
	else if (s == negInf)      sprintf (str, "-inf");
	else if (s <= negInf+2000) sprintf (str, "-inf+" scoreFmtSimple, s-negInf);
	else                       sprintf (str,         scoreFmtSimple, s);
	return str;
	}

include_link_to_string

#define snoopAlgorithm_1                                                     \
	if (snoop)                                                               \
		{                                                                    \
		char A25[26], B25[26];                                               \
		                                                                     \
		strncpy (A25, (char*) A+1, sizeof(A25));  A25[sizeof(A25)-1] = 0;    \
		strncpy (B25, (char*) B+1, sizeof(B25));  B25[sizeof(B25)-1] = 0;    \
		                                                                     \
		fprintf (stderr, "\nydrop_one_sided_align(" unsposSlashFmt ") %s"    \
						 "  A=%s (" unsposFmt ")"                            \
						 "  B=%s (" unsposFmt ")\n",                         \
						 io->anchor1, io->anchor2,                           \
						 (reversed)?" reversed":"forward",                   \
						 A25, M, B25, N);                                    \
		}

#define snoopAlgorithm_2                                                     \
	if (snoop)                                                               \
		{                                                                    \
		strcpy (ccSnoop, relative_to_infinity ((dq-1)->CC));                 \
		strcpy (ddSnoop, relative_to_infinity ((dq-1)->DD));                 \
		fprintf (stderr,"init edge, dynProg=%08lX dpLen=%08X\n",             \
					 (long) dynProg, dynProg->len);                          \
		fprintf (stderr, "\n[%3u,%3u]"                                       \
						 "  [%d].cc=%s  [%d].DD=%s"                          \
						 "  dq=%d\n",                                        \
						 0, 0, 0, ccSnoop, 0, ddSnoop, 0);                   \
		}

#define snoopAlgorithm_3                                                     \
	if (snoop)                                                               \
		{                                                                    \
		strcpy (ccSnoop, relative_to_infinity ((dq-1)->CC));                 \
		strcpy (ddSnoop, relative_to_infinity ((dq-1)->DD));                 \
		fprintf (stderr, "[%3u,%3u]"                                         \
						 "  [" unsposFmt "].cc=%s"                           \
						 "  [" unsposFmt "].DD=%s"                           \
						 "  dq=%d\n",                                        \
						 0, col, col, ccSnoop, col, ddSnoop, col);           \
		}

#define snoopAlgorithm_4A                                                    \
	if (snoop)                                                               \
		fprintf (stderr, "\nrow " unsposFmt                                  \
		                 "  LY=" unsposFmt                                   \
		                 " RY=" unsposFmt,                                   \
		                 row, LY, RY);

#define snoopAlgorithm_4B                                                    \
	if (snoop)                                                               \
		fprintf (stderr, " ->  L=" sgnposFmt                                 \
		                 " R=" sgnposFmt                                     \
		                 "  LY=" unsposFmt                                   \
		                 " RY=" unsposFmt,                                   \
		                 L, R, LY, RY);

#define snoopAlgorithm_4C                                                    \
	if (snoop)                                                               \
		fprintf (stderr, "  tbNeeded=%d  dp=%ld  dq=%ld\n",                  \
		                 tbNeeded, (long) (dp-dynProg->p), (long) (dq-dynProg->p));

#define snoopAlgorithm_5                                                     \
	if (snoop)                                                               \
		{                                                                    \
		strcpy (ccSnoop, relative_to_infinity (c));                          \
		strcpy (ddSnoop, relative_to_infinity (dp->DD));                     \
		strcpy (iiSnoop, relative_to_infinity (i));                          \
		fprintf (stderr, "[%3u,%3u]  %c%c  B[%ld]=%c"                        \
						 "  c=%s  d=%s  i=%s"                                \
						 "  dp=%ld  dq=%ld",                                 \
						 row, col, A[row], B[col], (long) (b-B), *b,         \
						 ccSnoop, ddSnoop, iiSnoop,                          \
						 (long) (dp-dynProg->p), (long) (dq-dynProg->p));    \
		}

#define snoopAlgorithm_5A                                                    \
	if (snoop)                                                               \
		fprintf (stderr, "\n");

#define snoopAlgorithm_6                                                     \
	if (snoop)                                                               \
		{                                                                    \
		strcpy (ccSnoop, relative_to_infinity ((dq-1)->CC));                 \
		strcpy (ddSnoop, relative_to_infinity ((dq-1)->DD));                 \
		fprintf (stderr, " [%lu].CC=%s [%lu].DD=%s link=%s\n",               \
						 (long) ((dq-1)-dynProg->p), ccSnoop,                \
						 (long) ((dq-1)-dynProg->p), ddSnoop,                \
						 link_to_string (link));                             \
		}

#define snoopAlgorithm_7A                                                    \
	if (snoop)                                                               \
		{                                                                    \
		if ((rightSeg != NULL) && (R > 0))                                   \
			fprintf (stderr, "NN <- " sgnposFmt " (from R)\n", R-1);         \
		else                                                                 \
			fprintf (stderr, "NN <- " sgnposFmt " (from N)\n", (sgnpos)N);   \
		}

#define snoopAlgorithm_7B                                                    \
	if (snoop)                                                               \
		fprintf (stderr, "RY <- " unsposFmt                                  \
		                 " (hit ydrop prior to RY)\n",                       \
		                 RY);

#define snoopAlgorithm_7C                                                    \
	if (snoop)                                                               \
		prevRY = RY;

#define snoopAlgorithm_7D                                                    \
	if ((snoop) && (RY != prevRY))                                           \
		fprintf (stderr, "RY <- " unsposFmt                                  \
		                 " (feasible overhang)\n",                           \
		                 RY);

#define snoopAlgorithm_7E                                                    \
	if (snoop)                                                               \
		fprintf (stderr, "RY <- " unsposFmt                                  \
		                 " (room at right boundary)\n",                      \
		                 RY);

#define snoopAlgorithm_8                                                     \
	if (snoop)                                                               \
		fprintf (stderr, "alignment ends at [%3u,%3u]\n", *_end1, *_end2);

#endif // snoopAlgorithm

#ifdef snoopAlgorithmTrap
static unspos trapRow = 497;
static unspos trapCol = (unspos) -1;
#endif // snoopAlgorithmTrap


//=== stuff for snoopTraceback ===

#ifndef snoopTraceback
#define snoopTraceback_1    ;
#define snoopTraceback_2    ;
#define snoopTraceback_3    ;
#endif // not snoopTraceback

#ifdef snoopTraceback

#define snoopTraceback_1                                                     \
	fprintf (stderr, "(tbp=%08X) tb[%3d,%3d] <- %-2d %-9s\n",                \
					 (u32) tbp-1,                                            \
					 row, col, link, link_to_string(link));

#define snoopTraceback_2                                                     \
	fprintf (stderr, "(tbp=%08X) tb[%3d,%3d] is %-2d %-9s"                   \
	                 " prevOp=%s --> op=%d",                                 \
					 (u32) &tb->space[tbRow[row] + col],                     \
					 row, col, link, link_to_string(link),                   \
					 op_string(prevOp), op);

#define snoopTraceback_3                                                     \
	fprintf (stderr, " -> %s -> row %d, col %d\n",                           \
					 op_string(op), row, col);

#ifndef snoopAlgorithm
include_link_to_string
#endif //snoopAlgorithm

#endif // snoopTraceback


//=== stuff for snoopSubprobs ===

#ifndef snoopSubprobs
#define snoopSubprobsB_1 ;
#define snoopSubprobsB_2 ;
#define snoopSubprobsB_3 ;
#endif // not snoopSubprobs

#ifdef snoopSubprobs

#define snoopSubprobsB_1                                                     \
	if (leftSeg == NULL)                                                     \
		fprintf (stderr, "\n    leftSeg=(none)");                            \
	else                                                                     \
		fprintf (stderr, "\n    leftSeg=%s"                                  \
		                 " " unsposSlashFmt                                  \
		                 " " unsposSlashFmt,                                 \
		                   (leftSeg->type==diagSeg)? "diag"                  \
		                 : (leftSeg->type==horzSeg)? "horz"                  \
		                 : (leftSeg->type==vertSeg)? "vert"                  \
		                                           : "????",                 \
		                 leftSeg->b1, leftSeg->b2,                           \
		                 leftSeg->e1, leftSeg->e2);                          \
	if (rightSeg == NULL)                                                    \
		fprintf (stderr, "\n    rightSeg=(none)");                           \
	else                                                                     \
		fprintf (stderr, "\n    rightSeg=%s"                                 \
		                 " " unsposSlashFmt                                  \
		                 " " unsposSlashFmt,                                 \
		                   (rightSeg->type==diagSeg)? "diag"                 \
		                 : (rightSeg->type==horzSeg)? "horz"                 \
		                 : (rightSeg->type==vertSeg)? "vert"                 \
		                                           : "????",                 \
		                 rightSeg->b1, rightSeg->b2,                         \
		                 rightSeg->e1, rightSeg->e2);                        \
	if (!reversed)                                                           \
		fprintf (stderr, "\n    L=" sgnposFmt " R=" sgnposFmt,               \
		                 L, R);                                              \
	else                                                                     \
		{                                                                    \
		sgnpos tempL = L;                                                    \
		sgnpos tempR = R;                                                    \
		sgnpos tempT = 0;                                                    \
		if      ((leftSeg == NULL) && (rightSeg != NULL)) {                tempL = -R+1;  tempR = N+1;   } \
		else if ((leftSeg != NULL) && (rightSeg == NULL)) { tempR = -L-1;  tempL = 0;                    } \
		else if ((leftSeg != NULL) && (rightSeg != NULL)) { tempT = -L-1;  tempL = -R+1;  tempR = tempT; } \
		fprintf (stderr, "\n    L=" sgnposFmt " R=" sgnposFmt               \
		                   " -> L=" sgnposFmt " R=" sgnposFmt,               \
		                 L, R, tempL, tempR);                                \
		}                                                                    \
	fprintf (stderr, "\n    ----");

#define snoopSubprobsB_2                                                     \
	fprintf (stderr, "\n    row " unsposFmt                                  \
					 " ->  L=" sgnposFmt                                     \
					 " R=" sgnposFmt                                         \
					 "  LY=" unsposFmt                                       \
					 " RY=" unsposFmt,                                       \
					 row, L, R, LY, RY);

#ifndef collect_stats
#define snoopSubprobsB_3                                                     \
	fprintf (stderr, "\n");                                                  \
	fprintf (stderr, "\trows=" unsposFmt, row);
#endif // not collect_stats

#ifdef collect_stats
#define snoopSubprobsB_3                                                     \
	fprintf (stderr, "\n");                                                  \
	fprintf (stderr, "\trows=" unsposFmt, row);                              \
	fprintf (stderr, "\tmaxDpRows=%s", commatize(gappedExtendStats.maxDpRows));
#endif // collect_stats

#endif // snoopSubprobs


//=== stuff for snoopBounds ===

#ifndef snoopBounds
#define debugSnoopBounds_1 ;
#define debugSnoopBounds_2 ;
#define debugSnoopBounds_3 ;
#endif // not snoopBounds

#ifdef snoopBounds

static char* segmentTypeName[3] = {"diag", "horz", "vert"};

#define debugSnoopBounds_1                                                   \
	{                                                                        \
	fprintf (stderr, "bounds: leftSeg(%s)=(" unsposFmt ")/" unsposFmt        \
	                 " anchor=(" unsposFmt ")/" unsposFmt "\n",              \
	                 segmentTypeName[(u8)leftSeg->type],                     \
	                 leftSeg->b1, leftSeg->b2,                               \
	                 anchor1, anchor2);                                      \
	if (leftSeg->type == diagSeg)                                            \
		fprintf (stderr, "bounds: L=(" unsposFmt "-" unsposFmt ")"           \
		                 "-(" unsposFmt "-" unsposFmt ")"                    \
		                 "=" sgnposFmt "\n",                                 \
		                 leftSeg->b2, anchor2, leftSeg->b1, anchor1, L);     \
	else                                                                     \
		fprintf (stderr, "bounds: L=(" unsposFmt "-" unsposFmt ")"           \
		                 "=" sgnposFmt "\n",                                 \
		                 leftSeg->b2, anchor2, L);                           \
	}

#define debugSnoopBounds_2                                                   \
	{                                                                        \
	fprintf (stderr, "bounds: rightSeg(%s)=(" unsposFmt ")/" unsposFmt       \
	                 " anchor=(" unsposFmt ")/" unsposFmt "\n",              \
	                 segmentTypeName[(u8)rightSeg->type],                    \
	                 rightSeg->b1, rightSeg->b2,                             \
	                 anchor1, anchor2);                                      \
	if (rightSeg->type == diagSeg)                                           \
		fprintf (stderr, "bounds: R=(" unsposFmt "-" unsposFmt ")"           \
		                 "-(" unsposFmt "-" unsposFmt ")"                    \
		                 "=" sgnposFmt "\n",                                 \
		                 rightSeg->b2, anchor2, rightSeg->b1, anchor1, R);   \
	else                                                                     \
		fprintf (stderr, "bounds: R=(" unsposFmt "-" unsposFmt ")"           \
		                 "=" sgnposFmt "\n",                                 \
		                 rightSeg->b2, anchor2, R);                          \
	}

#define debugSnoopBounds_3                                                   \
	fprintf (stderr, "bounds: swapped to L=" sgnposFmt " R=" sgnposFmt "\n", \
	                 L, R);

#endif // snoopBounds


//=== traceback row memory (see note (3)) ===
//
// minTbRowsNeeded is the inverse of the round_up result in tbrow_needed

#define minTbRowsNeeded (((((512*1024)/sizeof(u32))-1)*16)/17)

static u32*	tbRow = NULL;			// memory to track row positions in
static u32	tbRowLen = 0;			// .. traceback array

static void tbrow_needed (u32 rowsNeeded);
static void tbrow_needed (u32 rowsNeeded)
	{
	size_t needed;

	if (rowsNeeded <= tbRowLen) return;

	needed = round_up(sizeof(u32)*(rowsNeeded+1+rowsNeeded/16), 512*1024);
	tbRow = realloc_or_die ("ydrop_one_sided_align tbRow", tbRow, needed);
	tbRowLen = needed / sizeof(u32);
	}

void free_traceback_rows (void)
	{
	free_if_valid ("free_traceback_rows", tbRow);
	tbRow = NULL;
	tbRowLen = 0;
	}


//=== ydrop_one_sided_align ===

static score ydrop_one_sided_align
   (alignio*	io,
	int			reversed,			// true => search to the lower-left
	u8*			A,					// vertical sequence text
	u8*			B,					// horizontal sequence text
	unspos		M,					// vertical sequence length
	unspos		N,					// horizontal sequence length
	int			trimToPeak,
	editscript** script,
	unspos*		_end1,
	unspos*		_end2)
	{
	tback*		tb;					// tb is space provided to record traceback;
	u8*			tbp;				// .. tbp is the current traceback cell
	int			tbLen;				// .. and steps linearly from tbp->space;
	int			tbNeeded;			// .. tbRow[r] is the conceptual start of
									// .. the traceback cells for row r;  it is
									// .. indexed as tbRow[r][c] for c=LY..R',
									// .. where R' is the last cell used in
									// .. row r
	unspos		anchor1, anchor2;	// anchor positions in A and B
	scorerow*	allSub;				// substitution scores matrix
	score*		sub;				// substitution scores vector, current row
	score		gapOE, gapE;		// gap penalties
	score		yDrop;				// score drop threshold
	int			yDropTail;			// length of shortest score fall >= yDrop

									// nota bene:  row, col, leftCol, L, R,
									// .. LY, RY, prevLY, NN, npCol, end1
									// .. and end2 are all relative to the DP
									// .. matrix

	unspos		row, col = 0;		// current DP grid point
	unspos		leftCol;            // (copy of left column, for stats)
	sgnpos		L,  R;				// external column limits for current row
	unspos		LY, RY;				// actual column limits for current row;
									// .. ("Y" is for y-drop, not cartesian Y)
	unspos		prevLY;				// left column limit for previous row
	sgnpos		NN;					// truncated right side bound, current row
	unspos		npCol;				// last non-pruned cell in current row
	unspos		end1, end2;			// end of optimal alignment
	int			endIsBoundary;		// true => report boundaryScore instead of
									//         .. bestScore
	score		bestScore;			// score of best alignment seen so far
	score		boundaryScore;		// score of best alignment seen so far that
									// .. ends at the end of either sequence
	dpMatrix*	dynProg;			// DP cells
	dpCell*		dp;					// scans previous row of dp matrix
	dpCell*		dq;					// scans current row of dp matrix
	u8*			b;					// scans horizontal sequence
	score		c, d, i;			// running scores for DP cells
	score		cOpen, cNext, cTemp;// scratch values for cell scores
	activeseg*	active;				// list of segments that intersect the
									// .. sweep row within the feasible region
	galign*		alignList;
	galign*		rightAlign, *leftAlign;
	aliseg*		rightSeg,   *leftSeg;
	u8			link = 0;			// traceback link
	u8			op, prevOp;			// edit operations
#if ((defined snoopAlgorithm) || (defined snoopAlgorithmTrap))
	int			snoop = true;
#ifdef snoopAlgorithm
	unspos		prevRY = 0;
#endif // snoopAlgorithm
#endif // snoopAlgorithm or snoopAlgorithmTrap

	// sanity check;  if either sequence is empty there's no alignment

	if ((N <= 0) || (M <= 0))
		{ *(_end1) = *(_end2) = 0;  return 0; }

#if ((defined snoopAnchors) || ((defined snoopAlgorithm) && (defined debugPosA1)))
	fprintf (stderr,"ydrop_one_sided_align(" unsposSlashFmt ") %s\n",
	                 io->anchor1, io->anchor2, (reversed)?" reversed":"");
#endif // snoopAnchors

#ifdef snoopAlgorithm
#if ((defined debugPosA1) && (defined debugPosB1))
	snoop = ((io->anchor1 == debugPosA1) && (io->anchor2 == debugPosA2))
	     || ((io->anchor1 == debugPosB1) && (io->anchor2 == debugPosB2));
#elif (defined debugPosA1)
	snoop = (io->anchor1 == debugPosA1) && (io->anchor2 == debugPosA2);
#elif (defined debugPosB1)
	snoop = (io->anchor1 == debugPosB1) && (io->anchor2 == debugPosB2);
#endif // debugPosA1,debugPosB1
#endif // snoopAlgorithm

	gapped_extend_count_stat (numExtensions);
	dbg_timing_count_stat    (numExtensions);

	snoopAlgorithm_1;

	// extract scoring constants

	allSub = io->scoring->sub;
	gapE   = io->scoring->gapExtend;
	gapOE  = io->scoring->gapOpen + gapE;	// (any new gap gets both penalties)
	yDrop  = io->yDrop;

	tb    = io->tb;
	tbLen = tb->size;

	if (gapE != 0)
		yDropTail = (yDrop/gapE) + 6;
	else
		{
		// when gapE is zero, the above results would be infinite;  but we can
		// limit yDropTail to the distance from the length of the sequence;
		// this can increase the amount of memory needed;  the solution here is
		// not completely sufficient;  "truncating alignment" reports are still
		// likely

		int maxYDropTail = 500*1000;
		if (N < (unsigned int) maxYDropTail) yDropTail = N+1;
		                                else yDropTail = maxYDropTail;
		}

	// determine initial left and right constraints

	L = 0;
	R = N+1;								// (in blastz this was R=N)
	anchor1 = io->anchor1;
	anchor2 = io->anchor2;

	leftSeg = io->leftSeg;
	if (leftSeg != NULL)
		{
		L = signed_difference (leftSeg->b2, anchor2);
		if (leftSeg->type == diagSeg)
			L -= signed_difference (leftSeg->b1, anchor1);
		debugSnoopBounds_1;
		}

	rightSeg = io->rightSeg;
	if (rightSeg != NULL)
		{
		R = signed_difference (rightSeg->b2, anchor2);
		if (rightSeg->type == diagSeg)
			R -= signed_difference (rightSeg->b1, anchor1);
		debugSnoopBounds_2;
		}

	snoopSubprobsB_1;

	// if we're doing a reversed alignment we need to swap the L-R bounds (see
	// note (14))

	if (reversed)
		{
		sgnpos temp = 0;	// (placate compiler)
		if      ((leftSeg == NULL) && (rightSeg != NULL)) {               L = -R+1;  R = N+1;  }
		else if ((leftSeg != NULL) && (rightSeg == NULL)) { R    = -L-1;  L = 0;               }
		else if ((leftSeg != NULL) && (rightSeg != NULL)) { temp = -L-1;  L = -R+1;  R = temp; }
		debugSnoopBounds_3;
		}

	active     = NULL;
	rightAlign = io->rightAlign;
	leftAlign  = io->leftAlign;
	alignList  = (!reversed)? io->aboveList
	                        : io->belowList;

	// make sure we have a reasonable number of traceback rows to start with
	// (see note (3))

	tbrow_needed (minTbRowsNeeded);

	tbRow[0] = 0;
	tbp = tb->space;

	//////////
	// compute first row of dp matrix
	//////////

	// make sure we have enough traceback and dp space for the first row

	tbNeeded = yDropTail;
	if (tbNeeded > tbLen)
		suicide ("not enough space in trace_back array");

	dynProg = zalloc_or_die ("ydrop_one_sided_align dynProg", sizeof(dpMatrix));
	dp_ready (dynProg, tbNeeded);
	gapped_extend_count_stat (zallocCallsB);
	gapped_extend_add_stat   (zallocTotalB, sizeof(dpMatrix));

	// compute first row of dp matrix (see notes (8), (10) and (13))

	dq = dynProg->p;
	dq->CC = cTemp = 0;						// set C[0][0]
	c = (dq++)->DD = -gapOE;				// set D[1][0]
	*(tbp++) = 0;

	snoopAlgorithm_2;

	for (col=1 ; (col<=N)&&(cTemp>=-yDrop) ; col++)
		{
		dq->CC = cTemp = c;					// set C[0][col]
		(dq++)->DD = c - gapOE;				// set D[1][col]
		c -= gapE;
		*(tbp++) = cFromI;

		snoopAlgorithm_3;
		}

	gapped_extend_add_stat (dpCellsVisited, col);

	LY = 0;
	RY = col; // (1 column beyond the feasible region)

	//////////
	// compute additional rows of DP matrix
	//////////

	end1 = end2 = 0;
	bestScore = 0;
	boundaryScore = negInf;
	endIsBoundary = false;

	for (row=1; row<=M ; row++)
		{
		snoopAlgorithm_4A;

#ifdef snoopAlgorithmTrap
		if ((snoop) && (row == trapRow) && (trapCol == (unspos) -1))
			cTemp = 0;
#endif

		// update sweep row bounds, active segments, masking

		prevLY = LY;
		update_LR_bounds   (reversed,
		                    &rightSeg, &leftSeg, &rightAlign, &leftAlign,
		                    row, anchor1, anchor2, &L, &R, &LY, &RY);
		update_active_segs (reversed, &active, &alignList, dynProg->p-prevLY,
		                    row, anchor1, anchor2, LY, RY);

		snoopAlgorithm_4B;
		snoopSubprobsB_2;

		// make sure we have enough traceback and dp space for this row (see
		// note (3))

		tbrow_needed (row+1);
		gapped_extend_max_stat (maxDpRows, row+1);

		if (RY < LY) RY = LY;	// (see note 11)
		tbNeeded = RY - LY + yDropTail;
		if ((tbp - tb->space) + tbNeeded >= tbLen)
			{
			if (gapped_extend_inhibitTruncationReport)
				goto dp_finished;

			if (!reversed)
				fprintf (stderr, "truncating alignment ending at (" unsposCommaFmt ");",
				                 end1 + anchor1 + 1, end2 + anchor2 + 1);
			else
				fprintf (stderr, "truncating alignment starting at (" unsposCommaFmt ");",
				                 anchor1 + 2 - end1, anchor2 + 2 - end2);
			fprintf(stderr, "  anchor at (" unsposCommaFmt ")\n", anchor1, anchor2);
			goto dp_finished;
			}
		tbRow[row] = (tbp - tb->space) - LY;

		// set up DP pointers for this sweep row (see note (5))

		dp_ready (dynProg, tbNeeded);	// make sure we have enough DP space
		dq = dynProg->p;				// dq cells start at col == LY
		dp = dq + LY - prevLY;			// dp cells start at col == prevLY

		snoopAlgorithm_4C;

		// compute DP values for all bounded columns in this row (see note (4))

		sub = allSub[A[row]];

		col = leftCol = LY;
		b = B + col + 1;	// (b scans horizontal sequence, one column ahead)
		npCol = col;		// npCol records the last non-pruned position

		i = negInf;			// 'set' I[row][col]
		c = negInf;			// propose C[row][col]

		for ( ; (col<RY)&&((unspos)(b-B)<=N+1) ; col++)
			{
			d = dp->DD;						// get D[row][col]

#ifdef snoopAlgorithmTrap
			if ((snoop) && (row == trapRow) && (col == trapCol))
				cTemp = 0;
#endif

			// at this point d, i and c contain the DP values for the cell at
			// [row][col];  the latter is the value for reaching C from previous
			// grid points, but since C also has edges entering it from the
			// current D and I nodes, we might improve it
			// nota bene: when we *can* improve c, we make an arbitrary choice
			// to prefer deletion to insertion (when i and d are equal)
			// nota bene 2: all paths through this series of ifs assign a value
			// to link

			if ((active != NULL) && (dp->mask == row))
				{ snoopAlgorithm_5;  prune;  snoopAlgorithm_5A;  continue; }

			if ((d > c) || (i > c))			// === we CAN improve C ===
				{
				// nota bene: both iExtend and dExtend are set here because
				// the value of the C and I (or C and D) are equal, so traceback
				// may as well take a gap extend into this cell
				if (d >= i) { c = d;  link = cFromD | iExtend | dExtend; }
				       else { c = i;  link = cFromI | iExtend | dExtend; }
				snoopAlgorithm_5;
				if (c < bestScore - yDrop)
					{ prune;  snoopAlgorithm_5A;  continue; }

#ifndef allowBackToBackGaps
				// not allowing back-to-back gaps, so we don't need to consider
				// opening a gap here
				i -= gapE;					// 'set' I[row][col+1]
				dq->DD = d - gapE;			// set D[row+1][col]
#else
				// back-to-back gaps are allowed, so we must consider gap opens
				cOpen = c - gapOE;
				d -= gapE;					// set D[row+1][col]
				if (cOpen > d) { dq->DD = cOpen;  link &= ~dExtend; }
				          else   dq->DD = d;

				i -= gapE;					// 'set' I[row][col+1]
				if (cOpen > i) { i = cOpen;  link &= ~iExtend; }
#endif // allowBackToBackGaps
				}
			else							// === we CANNOT improve C ===
				{
				snoopAlgorithm_5;
				if (c < bestScore - yDrop)
					{ prune;  snoopAlgorithm_5A;  continue; }

				if (c >= bestScore)
					{ bestScore = c;  end1 = row;  end2 = col;  endIsBoundary = false; }
				if ((!trimToPeak)
				      && (c >= boundaryScore)
				      && ((row == M) || (col == N)))
					{ boundaryScore = c;  end1 = row;  end2 = col;  endIsBoundary = true; }

				cOpen = c - gapOE;
				d -= gapE;					// set D[row+1][col]
				if (cOpen > d) { dq->DD = cOpen;  link = cFromC; }
				          else { dq->DD = d;      link = cFromC | dExtend; }

				i -= gapE;					// 'set' I[row][col+1]
				if (cOpen > i) i = cOpen;
				          else link |= iExtend;
				}

			npCol = col;					// save as last non-pruned position

			// save C for this column, and compute proposed C for the next
			// column (see note (6))

			cNext = (dp++)->CC+sub[*(b++)];	// propose C[row][col+1]
			(dq++)->CC = c;					// set C[row][col]
			c = cNext;
			*(tbp++) = link;				// set link into C[row][col]

			snoopAlgorithm_6;
			//snoopTraceback_1;
			}

		gapped_extend_add_stat (dpCellsVisited, col-leftCol);

		// if the feasible region is empty, we're done

		if (LY >= RY)
			goto dp_finished;

		// finish up this row, by either moving the right bound left or
		// prolonging the row to support an overhang on the row above

		snoopAlgorithm_7A;
		NN = ((rightSeg != NULL) && (R > 0))? (R-1) : ((sgnpos) N);

		if (RY > npCol+1)					// we hit ydrop prior to RY
			{
			RY = npCol+1;
			snoopAlgorithm_7B;
			}
		else
			{
			// the current row reached its right bound, but the row above may
			// still have a feasible overhang so we prolong this row with
			// insertions (see note (9))

			snoopAlgorithm_7C;
			while ((i >= bestScore - yDrop) && (((sgnpos)RY) <= NN))
				{
				if (((u32)(dq - dynProg->p)) >= dynProg->len)
					suicidef("(in ydrop_one_sided_align:%d, dq-dynProg->p==%d, dynProg->len=" unsposFmt ")",
					         __LINE__, dq - dynProg->p, dynProg->len);
				dq->CC = i;					// set C[row][col]
				(dq++)->DD = i - gapOE;		// set D[row+1][col]
				i -= gapE;					// 'set' I[row][col+1]
				*(tbp++) = cFromI;
				RY++;
				}
			snoopAlgorithm_7D;
			}

		// terminate the cell at the right boundary, so that nothing will
		// step from it (termination occurs if the if is false and we thus
		// *fail* to increment RY)

		if (((sgnpos)RY) <= NN)
			{
			if (((u32)(dq - dynProg->p)) >= dynProg->len)
				suicidef("(in ydrop_one_sided_align:%d, dq-dynProg->p==%d, dynProg->len=" unsposFmt ")",
				         __LINE__, dq - dynProg->p, dynProg->len);
			dq->DD = dq->CC = negInf;		// set D[row+1][col]
			RY++;							// .. and set C[row][col]
			snoopAlgorithm_7E;
			}
		}

dp_finished:
	snoopSubprobsB_3;
	*(_end1) = row = end1;
	*(_end2) = col = end2;

	snoopAlgorithm_8;

	//////////
	// traceback the alignment to create the edit script
	//////////

	for (prevOp=0 ; (row>=1) || (col>0) ; prevOp=op)
		{
		link = tb->space[tbRow[row] + col];
		op = link & cidBits;
		if ((prevOp == cFromI) && ((link & iExtend) != 0)) op = cFromI;
		if ((prevOp == cFromD) && ((link & dExtend) != 0)) op = cFromD;
		snoopTraceback_2

		if      (op == cFromI) {         col--;  edit_script_ins (script, 1); }
		else if (op == cFromD) { row--;          edit_script_del (script, 1); }
		else                   { row--;  col--;  edit_script_sub (script, 1); }
		snoopTraceback_3
		}

	filter_active_segs (&active, 2);	// (disposes of everything in the list)

	free_if_valid ("ydrop_one_sided_align dynProg->p", dynProg->p);
	free_if_valid ("ydrop_one_sided_align dynProg",    dynProg);

	if (endIsBoundary) return boundaryScore;
	              else return bestScore;
	}

//----------
//
// dp_ready--
//	Allocate the dynamic-programming "sweep row" array, and ensure that the
//	first n elements exist and have been initialized.
//
//----------
//
// Arguments:
//	dpMatrix*	dynProg:	The current sweep array.
//	unspos		needed:		The number of elements required.
//
// Returns:
//  nothing;  the contents of dynProg are (potentially) modified by adding
//	enough empty cells to satisfy the number needed
//
//----------

static void dp_ready
   (dpMatrix*	dynProg,
	unspos		needed)
	{
	u32			oldLen;
	int			added;

	if (dynProg == NULL)
		suicide ("dp_ready called with NULL pointer");

	if (needed > 0xFFFFFFFF)
		suicidef ("in dp_ready, number of DP cells needed exceeds limit\n"
		          "  " unsposFmt " cells requested", needed);

	// make sure that there's room for n cells

	if (dynProg->p == NULL)
		{
		oldLen = 0;
		dynProg->len = needed + 1000;
		dynProg->p = malloc_or_die ("dp_ready", dynProg->len * sizeof(dpCell));
		}
	else
		{
		oldLen = dynProg->len;
		if (needed > oldLen)
			{
			dynProg->len = needed + dynProg->len/16 + 1000;
			dynProg->p = realloc_or_die ("dp_ready", dynProg->p,
			                             dynProg->len * sizeof(dpCell));
			}
		}

	// 0..oldLen-1 are already in use, but zero the rest to clear their mask
	// fields

	added = dynProg->len - oldLen;
	if (added > 0)
		memset (dynProg->p + oldLen, 0, sizeof(dynProg->p[0])*added);

	gapped_extend_max_stat (maxDpColumns, needed);
	}

//----------
//
// msp_left_right--
//	Given MSP m, we search obi, the list of alignments ordered by beginning
//	point (in vertical/seq1), to find the alignments and their contained
//	segments (gap-free pieces) that are the closest to the left and right of
//	the MSP's anchor-point.
//
//----------
//
// Arguments:
//	galign*	obi:	The list of (previously computed) alignments, ordered by
//					.. increasing beginning point.  This may be NULL.
//	galign*	m:		The MSP to check against the list.
//
// Returns:
//	false if the anchor-point is in an already-computed alignment;  true
//	otherwise;  actual result values (the segment pointers) are deposited in
//	the MSP.
//
//----------

static int msp_left_right
   (galign*	obi,
	galign*	m)
	{
	unspos	pos1 = m->pos1;
	unspos	pos2 = m->pos2;
	unspos	right,  left;
	galign* mRight, *mLeft;
	aliseg* bRight, *bLeft;
	aliseg* bp;
	sgnpos	x;

	right  = left  = seqposInfinity;
	mRight = mLeft = NULL;
	bRight = bLeft = NULL;

	// process all alignments in the obi that overlap (along vertical axis)
	// m's anchor point

	for ( ; (obi!=NULL)&&(obi->pos1<=pos1) ; obi=obi->next)
		{
		if (obi->end1 < pos1)
			continue;

		// invariant: obi->pos1 <= pos1 <= obi->end1

		for (bp=obi->firstSeg ; bp!=NULL ; bp=bp->nextSeg)
			{ if (bp->e1 >= pos1) break; }
		if (bp == NULL)
			continue;

		// invariant: bp is the first segment (in this alignment) that
		// intersects the line y=pos1

		if (bp->type == horzSeg)
			suicide ("msp_left_right: cannot be horizontal");

		if (bp->type == diagSeg)
			x = signed_difference (bp->b2, pos2)		// x is how far to the
			  + signed_difference (pos1,   bp->b1);		// .. right of m the
		else // (bp->type == vertSeg)					// .. segment is, along
			x = signed_difference (bp->b2, pos2);		// .. the line y=pos1

		if (x == 0)
			return false;

		if ((x > 0) && ((unspos) x < right))
			{ // a new closest segment to the right
			right  = (unspos) x;
			mRight = obi;
			bRight = bp;
			}
		else if ((x < 0) && ((unspos) -x < left))
			{ // a new closest segment to the left
			left  = (unspos) -x;
			mLeft = obi;
			bLeft = bp;
			}
		}

	m->rightAlign1 = m->rightAlign2 = mRight;
	m->rightSeg1   = m->rightSeg2   = bRight;
	m->leftAlign1  = m->leftAlign2  = mLeft;
	m->leftSeg1    = m->leftSeg2    = bLeft;

	return true;
	}

//----------
//
// get_above_below--
//	In preparation for extending the anchor-point of an MSP, find the closest
//	alignment ending below the anchor-point and the closest alignment starting
//	above the anchor-point.
//
//----------
//
// Arguments:
//	alignio*	io:		io->anchor1 is the anchor's vertical/seq1 position;
//						this routine fills in io->belowList and io->aboveList
//	galign*		obi:	The list of (previously computed) alignments, ordered by
//						.. increasing beginning point.  This may be NULL.
//	galign*		oed:	The list of (previously computed) alignments, ordered by
//						.. decreasing ending point.  This  may be NULL.
//
// Returns:
//	nothing;  values are written to io->belowList and io->aboveList.
//
//----------

static void get_above_below
   (alignio*	io,
	galign*		obi,
	galign*		oed)
	{
	unspos		pos1 = io->anchor1;
	galign*		mp;

	for (mp=oed ; mp!=NULL ; mp=mp->prev)
		{ if (mp->end1 < pos1) break; }
	io->belowList = mp;

	for (mp=obi ; mp!=NULL ; mp=mp->next)
		{ if (mp->pos1 > pos1) break; }
	io->aboveList = mp;
	}

//----------
//
// align_left_right--
//	Given a new alignment, determine the alignments and their segments that are
//	the closest in either horizontal direction at both ends of the alignment.
//
//----------
//
// Arguments:
//	galign*	obi:	The list of (previously computed) alignments, ordered by
//					.. increasing beginning point.  This may be NULL.
//	galign*	m:		The alignment to check.
//
// Returns:
//	(nothing)
//
//----------

static void align_left_right
   (galign*	obi,
	galign*	m)
	{
	unspos	pos1 = m->pos1, pos2 = m->pos2;
	unspos	end1 = m->end1, end2 = m->end2;
	unspos	rightOfBottom,  rightOfTop,   leftOfBottom,   leftOfTop;
	galign*	mRightOfBottom, *mRightOfTop, *mLeftOfBottom, *mLeftOfTop;
	aliseg*	bRightOfBottom, *bRightOfTop, *bLeftOfBottom, *bLeftOfTop;
	aliseg*	bp;
	sgnpos	x;

	rightOfBottom  = rightOfTop    = leftOfBottom  = leftOfTop  = seqposInfinity;
	mRightOfBottom = mLeftOfBottom = mRightOfTop   = mLeftOfTop = NULL;
	bRightOfBottom = bRightOfTop   = bLeftOfBottom = bLeftOfTop = NULL;

	for ( ; obi!=NULL ; obi=obi->next)
		{
		if ((obi->pos1 > end1) || (obi->end1 < pos1))
			continue;

		// invariant: obi->pos1 <= end1   and   obi->end1 >= pos1
		// meaning:   obi overlaps m along vertical/seq1

		for (bp=obi->firstSeg ; bp!=NULL ; bp=bp->nextSeg)
			{ if ((bp->type != horzSeg) && (bp->e1 >= pos1)) break; }

		// invariant: bp is either NULL or the first non-horizontal segment (in
		// this alignment) that intersects the line y=pos1

		if ((bp != NULL) && (bp->b1 <= pos1))
			{
			if (bp->type == diagSeg)
				x = signed_difference (bp->b2, pos2)	// x is how far to the
				  + signed_difference (pos1,   bp->b1);	// .. right of m the
			else // (bp->type == vertSeg)				// .. segment is, along
				x = signed_difference (bp->b2, pos2);	// .. the line y=pos1

			if ((x > 0) && ((unspos) x < rightOfBottom))
				{
				rightOfBottom  = (unspos) x;
				mRightOfBottom = obi;
				bRightOfBottom = bp;
				}
			else if ((x < 0) && ((unspos) -x < leftOfBottom))
				{
				leftOfBottom  = (unspos) -x;
				mLeftOfBottom = obi;
				bLeftOfBottom = bp;
				}
			}

		for ( ; bp!=NULL ; bp=bp->nextSeg)
			{ if (bp->type != horzSeg && bp->e1 >= end1) break; }

		if ((bp != NULL) && (bp->type != horzSeg) && (bp->e1 >= end1))
			{
			if (bp->type == diagSeg)
				x = signed_difference (bp->b2, end2)	// x is how far to the
				  + signed_difference (end1,   bp->b1);	// .. right of m the
			else // (bp->type == vertSeg)				// .. segment is, along
				x = signed_difference (bp->b2, end2);	// .. the line y=end1

			if ((x > 0) && ((unspos) x < rightOfTop))
				{
				rightOfTop  = (unspos) x;
				mRightOfTop = obi;
				bRightOfTop = bp;
				}
			else if ((x < 0) && ((unspos) -x < leftOfTop))
				{
				leftOfTop  = (unspos) -x;
				mLeftOfTop = obi;
				bLeftOfTop = bp;
				}
			}
		}

	m->rightAlign1 = mRightOfBottom;
	m->rightSeg1   = bRightOfBottom;
	m->rightAlign2 = mRightOfTop;
	m->rightSeg2   = bRightOfTop;
	m->leftAlign1  = mLeftOfBottom;
	m->leftSeg1    = bLeftOfBottom;
	m->leftAlign2  = mLeftOfTop;
	m->leftSeg2    = bLeftOfTop;
	}

//----------
//
// insert_align--
//	Insert a new alignment into the two lists of alignments, one ordered by
//	increasing beginning point (in vertical/seq1) and the other ordered by
//	decreasing end point.
//
//----------
//
// Arguments:
//	galign*		m:		The MSP that was used to anchor the alignment.
//	galign**	obi:	Pointer to the list of (previously computed) alignments,
//						.. ordered by increasing beginning point.  The list may
//						.. be empty (*obi == NULL).  This is updated by this
//						.. function.
//	galign**	oed:	Pointer to the list of (previously computed) alignments,
//						.. ordered by decreasing ending point.  The list may
//						.. be empty (*oed == NULL).  This is updated by this
//						.. function.
//
// Returns:
//	(nothing)
//
//----------

//=== stuff for snoopBlocks ===

#ifndef snoopBlocks
#define debugSnoopBlocks_6 ;
#endif // not snoopBlocks

#ifdef snoopBlocks

#define debugSnoopBlocks_6                                                   \
	fprintf (stderr, "adding alignment block [%8p]"                          \
	                 " b " unsposFmt " " unsposFmt                           \
	                 " e " unsposFmt " " unsposFmt                           \
	                 " s " scoreFmtSimple "\n",                              \
	                 m, m->pos1, m->pos2, m->end1, m->end2,                  \
					 (m->align == NULL)? ((score) 0) : m->align->s);

#endif // snoopBlocks


static void insert_align
   (galign*		m,
   	galign**	_obi,
   	galign**	_oed)
	{
	galign* 	mp;		// (mp scans list, 
	galign*		mq;		//  .. mq is predecessor in scan direction)
	galign* 	obi = *_obi;
	galign* 	oed = *_oed;

	if (m->firstSeg == NULL)
		suicide ("insert_align: null first segment");

	debugSnoopBlocks_6;

	for (mq=NULL,mp=obi ; mp!=NULL ; mq=mp,mp=mp->next)
		{ if (mp->pos1 >= m->pos1) break; }

	if (mq != NULL) { mq->next = m;   m->next = mp; }
	           else { m->next = obi;  obi = m;      }

	for (mq=NULL,mp=oed ; mp!=NULL ; mq=mp,mp=mp->prev)
		{ if (mp->end1 <= m->end1) break; }

	if (mq != NULL) { mq->prev = m;   m->prev = mp; }
	           else { m->prev = oed;  oed = m;      }

	*(_obi) = obi; *(_oed) = oed;
	}

//----------
//
// update_LR_bounds--
//	As we move to a new DP row, update LY and RY, the actual column limits
//	for dynamic programming.
//
//----------
//
// Arguments:
//	int			reversed:	true  => the DP row advances downward
//							false => it advances upward
//	aliseg**	rightSeg:	Current constraining segments;  these may be
//	aliseg**	leftSeg:	.. updated by this function.
//	galign**	rightAlign:	Alignments containing those segments;  these may be
//	galign**	leftAlign:	.. updated by this function.
//	unspos		row:		The sweep row.
//	unspos		anchor1:	The position at which the alignment
//	unspos		anchor2:	.. began.
//	sgnpos*		L,			(pointer to) most recent bounding constraints;
//	sgnpos*		R:			.. these may be updated by this function.
//	unspos*		LY:			(pointer to) most recent Y-drop constraints;  these
//	unspos*		RY:			.. may lie strictly inside L and R;  these may be
//							.. updated by this function.
//
// Returns:
//  (nothing)
//
//----------

#if ((!defined snoopAlgorithm)&&(!defined snoopSubprobs))
#define snoopLR_forward_left_1  ;
#define snoopLR_forward_left_2  ;
#define snoopLR_forward_left_3  ;
#define snoopLR_forward_right_1 ;
#define snoopLR_forward_right_2 ;
#define snoopLR_forward_right_3 ;
#define snoopLR_reverse_left_1  ;
#define snoopLR_reverse_left_2  ;
#define snoopLR_reverse_left_3  ;
#define snoopLR_reverse_right_1 ;
#define snoopLR_reverse_right_2 ;
#define snoopLR_reverse_right_3 ;
#endif // not snoopAlgorithm and not snoopSubprobs


#if (defined snoopAlgorithm)

#define snoopLR_forward_left_1                                               \
	if (snoop)                                                               \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		if ((*leftSeg)->type == diagSeg)                                     \
			fprintf (stderr, "update_LR: L <- " sgnposFmt                    \
			                 " (on same diag seg)\n",                        \
			                 L);                                             \
		else                                                                 \
			fprintf (stderr, "update_LR: L <- " sgnposFmt                    \
			                 " (on same vert seg)\n",                        \
			                 L);                                             \
		}

#define snoopLR_forward_left_2                                               \
	if ((snoop) && (*leftSeg != NULL))                                       \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: L <- " sgnposFmt                        \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")\n",                                              \
		                 L,                                                  \
		                 (*leftSeg)->b1, (*leftSeg)->b2,                     \
		                 (*leftSeg)->e1, (*leftSeg)->e2);                    \
		}

#define snoopLR_forward_left_3                                               \
	if ((snoop) && (L > (sgnpos) LY))                                        \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: LY <- " unsposFmt                       \
		                 "\n",                                               \
		                 (unspos) max ((sgnpos) LY, L));                     \
		}


#define snoopLR_forward_right_1                                              \
	if (snoop)                                                               \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		if ((*rightSeg)->type == diagSeg)                                    \
			fprintf (stderr, "update_LR: L <- " sgnposFmt                    \
			                 " (on same diag seg)\n",                        \
			                 R);                                             \
		else                                                                 \
			fprintf (stderr, "update_LR: R <- " sgnposFmt                    \
			                 " (on same vert seg)\n",                        \
			                 R);                                             \
		}

#define snoopLR_forward_right_2                                              \
	if ((snoop) && (*rightSeg != NULL))                                      \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: R <- " sgnposFmt                        \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")\n",                                              \
		                 R,                                                  \
		                 (*rightSeg)->b1, (*rightSeg)->b2,                   \
		                 (*rightSeg)->e1, (*rightSeg)->e2);                  \
		}
#define snoopLR_forward_right_3                                              \
	if ((snoop) && (R < (sgnpos) RY))                                        \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: RY <- " unsposFmt                       \
		                 "\n",                                               \
		                 (unspos) min ((sgnpos) RY, R));                     \
		}

#define snoopLR_reverse_left_1                                               \
	if (snoop)                                                               \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		if ((*rightSeg)->type == diagSeg)                                    \
			fprintf (stderr, "update_LR: L <- " sgnposFmt                    \
			                 " (on same diag seg)\n",                        \
			                 L);                                             \
		else                                                                 \
			fprintf (stderr, "update_LR: L <- " sgnposFmt                    \
			                 " (on same vert seg)\n",                        \
			                 L);                                             \
		}

#define snoopLR_reverse_left_2                                               \
	if ((snoop) && (*rightSeg != NULL))                                      \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: L <- " sgnposFmt                        \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")\n",                                              \
		                 L,                                                  \
		                 (*rightSeg)->b1, (*rightSeg)->b2,                   \
		                 (*rightSeg)->e1, (*rightSeg)->e2);                  \
		}

#define snoopLR_reverse_left_3                                               \
	if ((snoop) && (L > (sgnpos) LY))                                        \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: LY <- " unsposFmt                       \
		                 "\n",                                               \
		                 (unspos) max ((sgnpos) LY, L));                     \
		}

#define snoopLR_reverse_right_1                                              \
	if (snoop)                                                               \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		if ((*leftSeg)->type == diagSeg)                                     \
			fprintf (stderr, "update_LR: R <- " sgnposFmt                    \
			                 " (on same diag seg)\n",                        \
			                 R);                                             \
		else                                                                 \
			fprintf (stderr, "update_LR: R <- " sgnposFmt                    \
			                 " (on same vert seg)\n",                        \
			                 R);                                             \
		}

#define snoopLR_reverse_right_2                                              \
	if ((snoop) && (*leftSeg != NULL))                                       \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: R <- " sgnposFmt                        \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")\n",                                              \
		                 R,                                                  \
		                 (*leftSeg)->b1, (*leftSeg)->b2,                     \
		                 (*leftSeg)->e1, (*leftSeg)->e2);                    \
		}

#define snoopLR_reverse_right_3                                              \
	if ((snoop) && (R < (sgnpos) RY))                                        \
		{                                                                    \
		if (!havePrinted)                                                    \
			{ fprintf (stderr, "\n");  havePrinted = true; }                 \
		fprintf (stderr, "update_LR: RY <- " unsposFmt                       \
		                 "\n",                                               \
		                 (unspos) min ((sgnpos) RY, R));                     \
		}

#endif // snoopAlgorithm


#if ((!defined snoopAlgorithm)&&(defined snoopSubprobs))

#define snoopLR_forward_left_1                                               \
	if (true)                                                                \
		{                                                                    \
		if ((*leftSeg)->type == diagSeg)                                     \
			fprintf (stderr, "\n    update_LR: L <- " sgnposFmt              \
			                 " (on same diag seg)",                          \
			                 L);                                             \
		else                                                                 \
			fprintf (stderr, "\n    update_LR: L <- " sgnposFmt              \
			                 " (on same vert seg)",                          \
			                 L);                                             \
		}

#define snoopLR_forward_left_2                                               \
	if (*leftSeg != NULL)                                                    \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: L <- " sgnposFmt                  \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")",                                                \
		                 L,                                                  \
		                 (*leftSeg)->b1, (*leftSeg)->b2,                     \
		                 (*leftSeg)->e1, (*leftSeg)->e2);                    \
		}

#define snoopLR_forward_left_3                                               \
	if (L > (sgnpos) LY)                                                     \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: LY <- " unsposFmt,                \
		                 (unspos) max ((sgnpos) LY, L));                     \
		}


#define snoopLR_forward_right_1                                              \
	if (true)                                                                \
		{                                                                    \
		if ((*rightSeg)->type == diagSeg)                                    \
			fprintf (stderr, "\n    update_LR: L <- " sgnposFmt              \
			                 " (on same diag seg)",                          \
			                 R);                                             \
		else                                                                 \
			fprintf (stderr, "\n    update_LR: R <- " sgnposFmt              \
			                 " (on same vert seg)",                          \
			                 R);                                             \
		}

#define snoopLR_forward_right_2                                              \
	if (*rightSeg != NULL)                                                   \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: R <- " sgnposFmt                  \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")",                                                \
		                 R,                                                  \
		                 (*rightSeg)->b1, (*rightSeg)->b2,                   \
		                 (*rightSeg)->e1, (*rightSeg)->e2);                  \
		}

#define snoopLR_forward_right_3                                              \
	if (R < (sgnpos) RY)                                                     \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: RY <- " unsposFmt,                \
		                 (unspos) min ((sgnpos) RY, R));                     \
		}

#define snoopLR_reverse_left_1                                               \
	if (true)                                                                \
		{                                                                    \
		if ((*rightSeg)->type == diagSeg)                                    \
			fprintf (stderr, "\n    update_LR: L <- " sgnposFmt              \
			                 " (on same diag seg)",                          \
			                 L);                                             \
		else                                                                 \
			fprintf (stderr, "\n    update_LR: L <- " sgnposFmt              \
			                 " (on same vert seg)",                          \
			                 L);                                             \
		}

#define snoopLR_reverse_left_2                                               \
	if (*rightSeg != NULL)                                                   \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: L <- " sgnposFmt                  \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")",                                                \
		                 L,                                                  \
		                 (*rightSeg)->b1, (*rightSeg)->b2,                   \
		                 (*rightSeg)->e1, (*rightSeg)->e2);                  \
		}

#define snoopLR_reverse_left_3                                               \
	if (L > (sgnpos) LY)                                                     \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: LY <- " unsposFmt,                \
		                 (unspos) max ((sgnpos) LY, L));                     \
		}

#define snoopLR_reverse_right_1                                              \
	if (true)                                                                \
		{                                                                    \
		if ((*leftSeg)->type == diagSeg)                                     \
			fprintf (stderr, "\n    update_LR: R <- " sgnposFmt              \
			                 " (on same diag seg)",                          \
			                 R);                                             \
		else                                                                 \
			fprintf (stderr, "\n    update_LR: R <- " sgnposFmt              \
			                 " (on same vert seg)",                          \
			                 R);                                             \
		}

#define snoopLR_reverse_right_2                                              \
	if (*leftSeg != NULL)                                                    \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: R <- " sgnposFmt                  \
		                 " (on new seg " unsposSlashFmt " " unsposSlashFmt   \
		                 ")",                                                \
		                 R,                                                  \
		                 (*leftSeg)->b1, (*leftSeg)->b2,                     \
		                 (*leftSeg)->e1, (*leftSeg)->e2);                    \
		}

#define snoopLR_reverse_right_3                                              \
	if (R < (sgnpos) RY)                                                     \
		{                                                                    \
		fprintf (stderr, "\n    update_LR: RY <- " unsposFmt,                \
		                 (unspos) min ((sgnpos) RY, R));                     \
		}

#endif // not snoopAlgorithm and snoopSubprobs


static void update_LR_bounds
   (int			reversed,
	aliseg**	rightSeg,
	aliseg**	leftSeg,
	galign**	rightAlign,
	galign**	leftAlign,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2,
	sgnpos*		_L,
	sgnpos*		_R,
	unspos*		_LY,
	unspos*		_RY)
	{
	sgnpos		L  = *_L;
	sgnpos		R  = *_R;
	unspos		LY = *_LY;
	unspos		RY = *_RY;
#ifdef snoopAlgorithm
	int			snoop = true;
	int			havePrinted = false;
#endif // snoopAlgorithm

	dbg_timing_gapped_extend_sub (debugClockUpdateLrBounds);

#ifdef snoopAlgorithm
#if ((defined debugPosA1) && (defined debugPosB1))
	snoop = ((anchor1 == debugPosA1) && (anchor2 == debugPosA2))
	     || ((anchor1 == debugPosB1) && (anchor2 == debugPosB2));
#elif (defined debugPosA1)
	snoop = (anchor1 == debugPosA1) && (anchor2 == debugPosA2);
#elif (defined debugPosB1)
	snoop = (anchor1 == debugPosB1) && (anchor2 == debugPosB2);
#endif // debugPosA1,debugPosB1
#endif // snoopAlgorithm

	// handle forward search (DP row advances upward)

	if (!reversed)
		{
		if (*leftSeg != NULL)
			{
			if ((*leftSeg)->e1 >= row + anchor1)	// stay on same segment
				{
				if ((*leftSeg)->type == diagSeg) L++;
				snoopLR_forward_left_1;
				}
			else									// move to next segment
				{
				L = next_sweep_seg (/*right*/ false, leftSeg, leftAlign,
				                    row, anchor1, anchor2) + 1;
				snoopLR_forward_left_2;
				}
			}

		if (*leftSeg != NULL) // (next_sweep_seg may have changed *leftSeg)
			{
			snoopLR_forward_left_3;
			LY = (unspos) max ((sgnpos) LY, L);
			}

		if (*rightSeg != NULL)
			{
			if ((*rightSeg)->e1 >= row + anchor1)	// stay on same segment
				{
				if ((*rightSeg)->type == diagSeg) R++;
				snoopLR_forward_right_1;
				}
			else									// move to next segment
				{
				R = next_sweep_seg (/*right*/ true, rightSeg, rightAlign,
				                    row, anchor1, anchor2) - 1;
				snoopLR_forward_right_2;
				}
			}

		if (*rightSeg != NULL) // (next_sweep_seg may have changed *rightSeg)
			{
			snoopLR_forward_right_3;
			RY = (unspos) min ((sgnpos) RY, R);
			}
		}

	// handle reversed search (DP row advances downward)

	else
		{
		if (*rightSeg != NULL)
			{
			if ((*rightSeg)->b1 <= anchor1 - row)	// stay on same segment
				{
				if ((*rightSeg)->type == diagSeg) L++;
				snoopLR_reverse_left_1;
				}
			else									// move to next segment
				{
				L = prev_sweep_seg (/*right*/ true, rightSeg, rightAlign,
				                    row, anchor1, anchor2) + 1;
				snoopLR_reverse_left_2;
				}
			}

		if (*rightSeg != NULL) // (prev_sweep_seg may have changed *rightSeg)
			{
			snoopLR_reverse_left_3;
			LY = (unspos) max ((sgnpos) LY, L);
			}

		if (*leftSeg != NULL)
			{
			if ((*leftSeg)->b1 <= anchor1 - row)	// stay on same segment
				{
				if ((*leftSeg)->type == diagSeg) R++;
				snoopLR_reverse_right_1;
				}
			else									// move to next segment
				{
				R = prev_sweep_seg (/*right*/ false, leftSeg, leftAlign,
				                    row, anchor1, anchor2) - 1;
				snoopLR_reverse_right_2;
				}
			}

		if (*leftSeg != NULL) // (prev_sweep_seg may have changed *leftSeg)
			{
			snoopLR_reverse_right_3;
			RY = (unspos) min ((sgnpos) RY, R);
			}
		}

	*_L  = L;
	*_R  = R;
	*_LY = LY;
	*_RY = RY;

	dbg_timing_gapped_extend_add (debugClockUpdateLrBounds);
	}

//----------
//
// next_sweep_seg--
//	Move to the next leftSeg or rightSeg in a forward (upward) sweep.
// prev_sweep_seg--
//	Move to the previous leftSeg or rightSeg in a reverse (downward) sweep.
//
//----------
//
// Arguments:
//	int			lookRight:	true  => move to the next alignment to the right
//							         .. when we run off the end of an alignment
//							false => move to next alignment to left instead
//	aliseg**	bp:			Current constraining segment;  may be updated by
//							.. this function.
//	galign**	mp:			Alignment containing that segment;  may be updated
//							.. by this function.
//	unspos		row:		The sweep row.
//	unspos		anchor1:	The position at which the alignment
//	unspos		anchor2:	.. began.
//
// Returns:
//  The column position of the next (or previous) leftSeg or rightSeg if there
//	is one;  otherwise, the column position of the next (or previous) leftAlign
//	or rightAlign.  If there is also no such alignment, zero is returned.
//
//----------

static sgnpos next_sweep_seg
   (int			lookRight,
	aliseg**	bp,
	galign**	mp,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2)
	{
	sgnpos		col;

	dbg_timing_gapped_extend_sub (debugClockNextSweepSeg);

	// move to the next segment, if there is one

	*bp = (*bp)->nextSeg;
	if (*bp != NULL)
		{
		if (((*bp)->type == horzSeg) && ((*bp = (*bp)->nextSeg) == NULL))
			suicide ("Last alignment segment was horizontal");

		col = signed_difference ((*bp)->b2, anchor2);
		dbg_timing_gapped_extend_add (debugClockNextSweepSeg);
		return col;
		}

	// we've run off the end (top) of an alignment;  move to the next one

	if (lookRight) { *bp = (*mp)->rightSeg2;  *mp = (*mp)->rightAlign2; }
	          else { *bp = (*mp)->leftSeg2;   *mp = (*mp)->leftAlign2;  }

	if (*bp == NULL)
		{
		dbg_timing_gapped_extend_add (debugClockNextSweepSeg);
		return 0;	// no constraint; there was no "next alignment"
		}

	// figure out the column where the start of the new segment intersects the
	// line y=row

	if ((*bp)->type == diagSeg)
		col = (sgnpos) row
		    + signed_difference ((*bp)->b2, anchor2)
		    - signed_difference ((*bp)->b1, anchor1);
	else  // we jumped to a vertical segment
		col = signed_difference ((*bp)->b2, anchor2);

	dbg_timing_gapped_extend_add (debugClockNextSweepSeg);
	return col;
	}


static sgnpos prev_sweep_seg
   (int			lookRight,
	aliseg**	bp,
	galign**	mp,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2)
	{
	sgnpos		col;

	dbg_timing_gapped_extend_sub (debugClockPrevSweepSeg);

	// move to the previous segment, if there is one

	*bp = (*bp)->prevSeg;
	if (*bp != NULL)
		{
		if (((*bp)->type == horzSeg) && ((*bp = (*bp)->prevSeg) == NULL))
			suicide ("First alignment segment was horizontal");
		col = signed_difference (anchor2, (*bp)->e2);
		dbg_timing_gapped_extend_add (debugClockPrevSweepSeg);
		return col;
		}

	// we've run off the front (bottom) of an alignment;  move to the previous
	// one

	if (lookRight) { *bp = (*mp)->rightSeg1;  *mp = (*mp)->rightAlign1; }
	          else { *bp = (*mp)->leftSeg1;   *mp = (*mp)->leftAlign1;  }

	if (*bp == NULL)
		{
		dbg_timing_gapped_extend_add (debugClockPrevSweepSeg);
		return 0;	// no constraint; there was no "previous alignment"
		}

	// figure out the column where the end of the new segment intersects the
	// line y=row

	if ((*bp)->type == diagSeg)
		col = (sgnpos) row
		    + signed_difference (anchor2, (*bp)->e2)
		    - signed_difference (anchor1, (*bp)->e1);
	else  // we jumped to a vertical segment
		col = signed_difference (anchor2, (*bp)->e2);

	dbg_timing_gapped_extend_add (debugClockPrevSweepSeg);
	return col;
	}

//----------
//
// update_active_segs--
//	As we move to a new sweep row, update the list of segments that intersect
//	the sweep row within the feasible region.
//
//----------
//
// Arguments:
//	int			reversed:	true  => the DP row advances downward
//							false => it advances upward
//	activeseg**	active:		The list of active segments.  This may be updated
//							.. by this function.
//	galign**	alignList:	Alignments in advance of the sweep row.  This may
//							.. be updated by this function.
//	dpCell*		dp:			First DP cell in (conceptual) sweep row.  This is
//							.. indexed from LY to RY, inclusive.
//	unspos		row:		The sweep row.
//	unspos		anchor1:	The position at which the alignment
//	unspos		anchor2:	.. began.
//	unspos		LY, RY:		Current Y-drop constraints.
//
// Returns:
//  (nothing)
//
//----------

static aliseg* next_seg (aliseg* bp, int reversed)
	{ return (reversed)? bp->prevSeg : bp->nextSeg; }

static void update_active_segs
   (int			reversed,
	activeseg**	_active,
	galign**	_alignList,
	dpCell*		dp,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2,
	unspos		LY,
	unspos		RY)
	{
	activeseg*	active    = *_active;
	galign*		alignList = *_alignList;
	activeseg*	act;

	dbg_timing_gapped_extend_sub (debugClockUpdateActiveSegs);

	// process currently active segments (those that intersect the sweep row)

	for (act=active ; act!=NULL ; act=act->next)
		{
		if (act->type == horzSeg)
			suicide ("Impossible horizontal segment.");

		if (act->lastRow >= row)
			{ // sweep row still intersects this segment
			if (act->type == diagSeg)
				act->x++;
			if ((act->x >= LY) && (act->x <= RY))
				dp[act->x].mask = row;
			}
		else if ((act->seg = next_seg(act->seg, reversed)) != NULL)
		   	{ // sweep row intersects the next segment of this alignment;
		   	// move to the next segment and mask its intial DP cells
			build_active_seg (reversed, act, dp, row, anchor1, anchor2, LY, RY);
			if (act->type == horzSeg)
				{
				act->seg = next_seg (act->seg, reversed);
				build_active_seg (reversed, act, dp,
				                  row, anchor1, anchor2, LY, RY);
				}
			}
		else
			{ // sweep row has passed the end of this alignment
			act->filter = 1;  // (mark it for deletion)
			}
		}

	// add any other alignments the sweep row now intersects, adding the
	// first segment (from the appropriate end) to the active list and
	// masking its intial DP cells

	if (!reversed)
		{
		while ((alignList!=NULL) && (alignList->pos1 - anchor1 == row))
			{
			active = add_new_active (reversed, active, alignList,
			                         dp, row, anchor1, anchor2, LY, RY);
			alignList = alignList->next;
			}
		}
	else
		{
		while ((alignList != NULL) && (anchor1 - alignList->end1 == row))
			{
			active = add_new_active (reversed, active, alignList,
			                         dp, row, anchor1, anchor2, LY, RY);
			alignList = alignList->prev;
			}
		}

	filter_active_segs (&active, 0); // delete blocks whose filter is not 0

	*_active    = active;
	*_alignList = alignList;

	dbg_timing_gapped_extend_add (debugClockUpdateActiveSegs);
	}

//----------
//
// build_active_seg--
//	Create an active segment record for a given segment, and mask any DP
//	cells that it intersects.
//
//----------
//
// Arguments:
//	int			reversed:	true  => the DP row advances downward
//							false => it advances upward
//	activeseg*	act:		The active segment record, which already contains
//							.. the proper segment, but nothing else.
//	dpCell*		dp:			First DP cell in (conceptual) sweep row.  This is
//							.. indexed from LY to RY, inclusive.
//	unspos		row:		The sweep row.
//	unspos		anchor1:	The position at which the alignment
//	unspos		anchor2:	.. began.
//	unspos		LY, RY:		Current Y-drop constraints.
//
// Returns:
//  (nothing)
//
//----------

static void build_active_seg
   (int			reversed,
	activeseg*	act,
	dpCell*		dp,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2,
	unspos		LY,
	unspos		RY)
	{
	unspos		horzEnd, iMin, iMax, i;

	act->type = act->seg->type;

	// nota bene:  the following assigns to act->x and act->lastRow always
	//             result in non-neagtive values

	if (!reversed)
		{
		act->x       = act->seg->b2 - anchor2;
		act->lastRow = act->seg->e1 - anchor1;
		}
	else
		{
		act->x       = anchor2 - act->seg->e2;
		act->lastRow = anchor1 - act->seg->b1;
		}

	if (act->type != horzSeg)
		{
		if ((act->x >= LY) && (act->x <= RY))
			dp[act->x].mask = row;
		}
	else
		{
		horzEnd = (!reversed)? act->seg->e2 - anchor2
		                     : anchor2 - act->seg->b2;
		iMin = max (LY, act->x);
		iMax = min (RY, horzEnd);
		for (i=iMin ; i<=iMax ; i++)
			dp[i].mask = row;
		}
	}

//----------
//
// add_new_active--
//	Create a new active segment record and add it to the list, containing
//	the given alignment's terminal segment.
//
//----------
//
// Arguments:
//	int			reversed:	true  => the DP row advances downward
//							false => it advances upward
//	activeseg*	active:		The list of active segments.  Upon return, the
//							.. caller should assign this function's return
//							.. value to this.
//	galign*		alignList:	Alignments in advance of the sweep row.  This may
//							.. be updated by this function.
//	dpCell*		dp:			First DP cell in (conceptual) sweep row.  This is
//							.. indexed from LY to RY, inclusive.
//	unspos		row:		The sweep row.
//	unspos		anchor1:	The position at which the alignment
//	unspos		anchor2:	.. began.
//	unspos		LY, RY:		Current Y-drop constraints.
//
// Returns:
//	Pointer to the head of the active segment list, which may be the newly
//	created node.
//
//----------

static activeseg* add_new_active
   (int			reversed,
   	activeseg*	active,
   	galign*		alignList, 
	dpCell*		dp,
	unspos		row,
	unspos		anchor1,
	unspos		anchor2,
	unspos		LY,
	unspos		RY)
	{
	activeseg* act = malloc_or_die ("add_new_active", sizeof(activeseg));

	act->filter = 0;
	if (!reversed) act->seg = alignList->firstSeg;
	          else act->seg = alignList->lastSeg;
	act->next = active;
	build_active_seg (reversed, act, dp, row, anchor1, anchor2, LY, RY);

	return act;
	}

//----------
//
// filter_active_segs--
//	Remove (dispose of) active segments with filter values NOT equal to some
//	specified value.
//
//----------
//
// Arguments:
//	activeseg**	active:	List of segments.
//	int			filter:	The filter value of segments to KEEP.  The active
//						.. records for all other segments are removed from the
//						.. list and disposed of.
//
// Returns:
//  (nothing)
//
//----------

static void filter_active_segs
   (activeseg**	active,
	int			filter)
	{
	activeseg* prevAct, *act;

	dbg_timing_gapped_extend_sub (debugClockFilterActiveSegs);

	for (prevAct=NULL,act=(*active); act!=NULL ; )
		{
		if (act->filter == filter)
			{
			prevAct = act;
			act = act->next;
			}
		else if (prevAct != NULL)
			{
			prevAct->next = act->next;
			free_if_valid ("filter_active_segs", act);
			act = prevAct->next;
			}
		else
			{
			*active = act->next;
			free_if_valid ("filter_active_segs", act);
			act = *active;
			}
		}

	dbg_timing_gapped_extend_add (debugClockFilterActiveSegs);
	}

//----------
//
// format_alignment--
//	Process the edit script for a newly computed alignment by (1) storing a
//	linked-list form with the seeding MSP and (2) augmenting the edit script
//	into a form suitable for returning to the calling program.
//
//----------
//
// Arguments:
//	alignio* io:	The alignment to format.
//	galign* m:		The MSP that was used to anchor the alignment.
//
// Returns:
//	Pointer to the newly allocated alignment description.
//
//----------

static alignel* format_alignment
   (alignio*	io,
	galign*		m)
	{
	unspos		beg1, end1, beg2, end2;
	unspos		height, width, i, j, startI, startJ, run;
	u32			opIx;
	editscript*	script;
	u8*			seq1, *seq2;
	alignel*	a;

	beg1   = io->start1 + 1;
	end1   = io->stop1  + 1;
	beg2   = io->start2 + 1;
	end2   = io->stop2  + 1;
	script = io->script;
	seq1   = io->seq1;
	seq2   = io->seq2;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		startI = i;  startJ = j;
		run = edit_script_run_of_subs (script, &opIx);
		i += run; j += run;

		save_seg (m, beg1+startI-1, beg2+startJ-1, beg1+i-2, beg2+j-2);
		if (i < height || j < width)
			edit_script_indel_len (script, &opIx, &i, &j);
    	}

	a = malloc_or_die ("format_alignment", sizeof(alignel));
	a->script = script;
	a->beg1   = beg1;  a->beg2 = beg2;
	a->end1   = end1;  a->end2 = end2;
	a->seq1   = seq1;  a->seq2 = seq2;
	a->s      = io->s;
	a->next   = NULL;
	a->isTrivial = false;

	return a;
	}

//----------
//
// save_seg--
//	Add a gap-free segment to a seeding MSP, inserting a vertical or horizontal
//	piece before it if appropriate.
//
//----------
//
// Arguments:
//	galign*	m:		The MSP to add the segment to.
//	unspos	b1, b2:	The segment's starting position.
//	unspos	e1, e2:	The segment's ending position.
//
// Returns:
//	(nothing)
//
//----------

static void insert_seg_to_tail (galign* mp, aliseg* bp);

static void save_seg
   (galign*	m,
	unspos	b1,
	unspos	b2,
	unspos	e1,
	unspos	e2)
	{
	aliseg* bp = malloc_or_die ("save_seg bp", sizeof(aliseg));
	aliseg* bq;

	bp->b1 = b1;
	bp->b2 = b2;
	bp->e1 = e1;
	bp->e2 = e2;
	bp->type = diagSeg;

	// if the alignment is empty, create it with this as the first segment

	if (m->firstSeg == NULL)
		{
		m->firstSeg = bp->prevSeg = bp->nextSeg = bp;
		return;
		}

	// otherwise, we insert it at the tail, with a preceding vertical or
	// horizontal segment;  we assume the previous tail was a diagonal segment
	// (since they are the only type ever added to the tail);  further, we
	// assume that the previous tail ends on either the y=b1-1 or x=b2-1 line
	// (but not both) 

	bq = malloc_or_die ("save_seg bq", sizeof(aliseg));
	bq->type = ((b1 == m->firstSeg->prevSeg->e1+1)? horzSeg : vertSeg);
	bq->b1 = m->firstSeg->prevSeg->e1 + 1;
	bq->b2 = m->firstSeg->prevSeg->e2 + 1;
	bq->e1 = b1 - 1;
	bq->e2 = b2 - 1;

	insert_seg_to_tail (m, bq);
	insert_seg_to_tail (m, bp);
	}

static void insert_seg_to_tail (galign* mp, aliseg* bp)
	{
	bp->prevSeg                    = mp->firstSeg->prevSeg;
	bp->nextSeg                    = mp->firstSeg;
	mp->firstSeg->prevSeg->nextSeg = bp;
	mp->firstSeg->prevSeg          = bp;
	}

//----------
// [[-- a seed hit reporter function --]]
//
// gappily_extend_hsps--
//	Perform a gapped extension of a seed hit or HSP.
//
// Arguments and Return value: (see seed_search.h)
//
//----------

u32 gappily_extend_hsps
   (void*	_info,
	unspos	pos1,
	unspos	pos2,
	unspos	length,
	score	s)
	{
	hitrepgappily*	info = (hitrepgappily*) _info;
	seq*			seq1 = info->seq1;
	seq*			seq2 = info->seq2;
	seqpartition*	sp1  = &seq1->partition;
	seqpartition*	sp2  = &seq2->partition;
	partition*		p1, *p2;
	unspos			peak;
	alignio			io;
	galign			mp;
	aliseg*			bp, *bq;
	u32				returnVal;

	//fprintf (stderr, "working on segment " unsposSlashFmt " " unsposFmt "\n",
	//                 pos1-length, pos2-length, length);

	if (gapped_extend_dbgShowHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides
			 (stderr, seq1, pos1-length, seq2, pos2-length, length);
		}

	// (move pos1/pos2 from end of segment to start)

	pos1 -= length;
	pos2 -= length;

	// $$$ move this test to set_up_hit_processor()

	if (info->scoreThresh.t != 'S')
		suicidef ("gappily_extend_hsps can't handle score threshold %s",
		          score_thresh_to_string (&info->scoreThresh));

	// reduce the HSP to a single point

	peak = segment_peak (seq1->v+pos1, seq2->v+pos2, length, info->scoring);
	pos1 += peak;
	pos2 += peak;
	//length = 0; // (unnecessary)

#ifdef debugHspImmediate
	fprintf (stderr, "hsp: " unsposSlashFmt " -> " unsposSlashFmt "\n",
	                 pos1-peak, pos2-peak, pos1, pos2);
#endif

	if (gapped_extend_dbgShowAnchors)
		{
		segtable st;
		segment* seg = &st.seg[0];

		st.size          = 1;
		st.len           = 1;
		st.haveScores    = true;
		st.coverageLimit = 0;
		st.coverage      = length;
		st.lowScore      = s;
		seg->pos1        = pos1;
		seg->pos2        = pos2;
		seg->length      = 0;
		seg->s           = s;
		seg->id          = seq2->revCompFlags;
		seg->scoreCov    = length;
		seg->filter      = false;

		write_segments (stderr, &st, seq1, seq2, false);
		}

	// build the alignio record for ydrop_align()

	io.seq1 = seq1->v;
	io.seq2 = seq2->v;
	io.rev1 = info->rev1;
	io.rev2 = info->rev2;
	io.low1 = 0;   io.len1 = io.high1 = seq1->len;
	io.low2 = 0;   io.len2 = io.high2 = seq2->len;

	io.scoring    = info->scoring;
	io.yDrop      = info->yDrop;
	io.trimToPeak = info->trimToPeak;

	io.anchor1 = pos1;
	io.anchor2 = pos2;

	if (sp1->p != NULL)
		{
		p1 = lookup_partition (seq1, io.anchor1);
		io.low1  = p1->sepBefore + 1;
		io.high1 = p1->sepAfter;
		}

	if (sp2->p != NULL)
		{
		p2 = lookup_partition (seq2, io.anchor2);
		io.low2  = p2->sepBefore + 1;
		io.high2 = p2->sepAfter;
		}

	if (info->traceback == NULL)
		suicide ("gappily_extend_hsps was given a NULL traceback pointer.");
	io.tb = info->traceback;

	// perform alignment

	io.leftAlign  = NULL;		// (convince ydrop_align() that there are no
	io.rightAlign = NULL;		//  .. neighboring/bounding alignments
	io.leftSeg    = NULL;		//  .. to worry about)
	io.rightSeg   = NULL;
	io.aboveList  = NULL;
	io.belowList  = NULL;

	ydrop_align (&io);			// (find the gapped alignment)

#ifdef debugHspImmediate
	fprintf (stderr, "  gappily: " unsposSlashFmt " " unsposSlashFmt " " scoreFmt "\n",
	                 io.start1, io.start2, io.stop1, io.stop2, io.s);
#endif

    if (io.s < info->scoreThresh.s)
		{
#ifdef debugHspImmediate
		fprintf (stderr, "  gappily: (fails score thresh, " scoreFmt "<" scoreFmt ")\n",
		                 io.s, info->scoreThresh.s);
#endif
		free_if_valid ("gappily_extend_hsps io.script", io.script);
		mp.firstSeg = NULL;
		mp.align    = NULL;
		goto return_zero;		// (the alignment score is too low)
		}

	mp.firstSeg = NULL;			// (convert the alignment to a linked list)
	mp.align    = format_alignment (&io, &mp);
	mp.pos1     = io.start1;
	mp.pos2     = io.start2;
	mp.end1     = io.stop1;
	mp.end2     = io.stop2;

	if (mp.firstSeg == NULL)
		goto return_zero;		// (the alignment is empty)

	mp.lastSeg = mp.firstSeg->prevSeg;	// (record the alignment's tail and
	mp.firstSeg->prevSeg				//  .. detach the circular pointer)
	  = mp.lastSeg->nextSeg = NULL;

	if ((info->minIdentity > 0) || (info->maxIdentity < 1))
		{
		mp.align = filter_aligns_by_identity
		             (seq1, seq2, mp.align, info->minIdentity, info->maxIdentity);
		if (mp.align == NULL) goto return_zero;
		}

	if ((info->minCoverage > 0) || (info->maxCoverage < 1))
		{
		mp.align = filter_aligns_by_coverage
		             (seq1, seq2, mp.align, info->minCoverage, info->maxCoverage);
		if (mp.align == NULL) goto return_zero;
		}

	if ((info->minContinuity > 0) || (info->maxContinuity < 1))
		{
		mp.align = filter_aligns_by_continuity
		             (mp.align, info->minContinuity, info->maxContinuity);
		if (mp.align == NULL) goto return_zero;
		}

	if (info->minMatchCount > 0)
		{
		mp.align = filter_aligns_by_match_count
		             (seq1, seq2, mp.align, info->minMatchCount);
		if (mp.align == NULL) goto return_zero;
		}

	if (info->maxMismatchCount >= 0)
		{
		mp.align = filter_aligns_by_mismatch_count
					  (seq1, seq2, mp.align, info->maxMismatchCount);
		if (mp.align == NULL) goto return_zero;
		}

	if (info->maxSeparateGapsCount >= 0)
		{
		mp.align = filter_aligns_by_num_gaps
					  (mp.align, info->maxSeparateGapsCount);
		if (mp.align == NULL) goto return_zero;
		}

	if (info->maxGapColumnsCount >= 0)
		{
		mp.align = filter_aligns_by_num_gap_columns
					  (mp.align, info->maxGapColumnsCount);
		if (mp.align == NULL) goto return_zero;
		}

	if (info->deGapifyOutput) print_align_list_segments (mp.align);
	                     else print_align_list          (mp.align);

	free_align_list (mp.align);
	mp.align = NULL;

	returnVal = 1;
	goto cleanup;

return_zero:
	returnVal = 0;
	goto cleanup;

cleanup:
	if ((mp.align != NULL) && (mp.align->script != NULL))
		free_if_valid ("gappily_extend_hsps mp.script", mp.align->script);

	for (bp=mp.firstSeg ; bp!=NULL ; bp=bq)
		{ bq = bp->nextSeg;  free_if_valid ("gappily_extend_hsps seg", bp); }

	return returnVal;
	}

//----------
//
// dump_alignio_input, dump_alignio_output--
//	Dump the contents of ydrop_align's io record.
//
//----------

#ifdef snoopAlignioInput

static void dump_alignio_input
   (FILE*		f,
	alignio*	io)
	{
	fprintf (f, "=====\n");

	fprintf (f, "seq1       = %8p\n",           io->seq1);
	fprintf (f, "rev1       = %8p\n",           io->rev1);
	fprintf (f, "len1       = " unsposFmt "\n", io->len1);
	fprintf (f, "low1       = " unsposFmt "\n", io->low1);
	fprintf (f, "high1      = " unsposFmt "\n", io->high1);
	fprintf (f, "anchor1    = " unsposFmt "\n", io->anchor1);

	fprintf (f, "seq2       = %8p\n",           io->seq2);
	fprintf (f, "rev2       = %8p\n",           io->rev2);
	fprintf (f, "len2       = " unsposFmt "\n", io->len2);
	fprintf (f, "low2       = " unsposFmt "\n", io->low2);
	fprintf (f, "high2      = " unsposFmt "\n", io->high2);
	fprintf (f, "anchor2    = " unsposFmt "\n", io->anchor2);

	fprintf (f, "scoring    = %8p\n",           io->scoring);
	fprintf (f, "yDrop      = " scoreFmt "\n",  io->yDrop);
	fprintf (f, "trim       = %s\n",            (io->trimToPeak)? "yes" : "no");
	fprintf (f, "tb         = %8p\n",           io->tb);

	fprintf (f, "leftAlign  = %8p\n",           io->leftAlign);
	fprintf (f, "rightAlign = %8p\n",           io->rightAlign);
	fprintf (f, "leftSeg    = %8p\n",           io->leftSeg);
	fprintf (f, "rightSeg   = %8p\n",           io->rightSeg);
	fprintf (f, "aboveList  = %8p\n",           io->aboveList);
	fprintf (f, "belowList  = %8p\n",           io->belowList);
	}

#endif // snoopAlignioInput


#ifdef snoopAlignioOutput

static void dump_alignio_output
   (FILE*		f,
	alignio*	io)
	{
	fprintf (f, "s          = " scoreFmt "\n",  io->s);
	fprintf (f, "  start1   = " unsposFmt "\n", io->start1);
	fprintf (f, "  stop1    = " unsposFmt "\n", io->stop1);
	fprintf (f, "  start2   = " unsposFmt "\n", io->start2);
	fprintf (f, "  stop2    = " unsposFmt "\n", io->stop2);
//	fprintf (f, "  script   = %p\n",            io->script);
	}

#endif // snoopAlignioOutput

//----------
//
// score_alignment--
//	Determine the score of gapped alignment in two subsequences.
//
// Note that this is generally used only after an alignment has been modified.
// The gapped_extend() produces the same alignment score as the alignment is
// found.
//
//----------
//
// Arguments:
//	scoreset*	scoring:	The scoring scheme to use.
//	seq*		seq1:		The first sequence.
//	unspos 		pos1:		The subsequence start position in seq1 (origin-0).
//	seq*		seq2:		The second sequence.
//	unspos 		pos2:		The subsequence start position in seq2 (origin-0).
//	editscript* s:			The script describing the alignment.
//
// Returns:
//	The alignment's score.
//
//----------

score score_alignment
   (scoreset*	scoring,
	seq*		seq1,
	unspos		pos1,
	seq*		seq2,
	unspos		pos2,
	editscript* script)
	{
	u8*			s1 = seq1->v + pos1;
	u8*			s2 = seq2->v + pos2;
	u8*			stop;
	u32			opIx;
	editop		op;
	u32			rpt;
	score		similarity = 0;

	for (opIx=0 ; opIx<script->len ; opIx++)
		{
		// score s = similarity;

		op  = script->op[opIx];
		rpt = edit_op_repeat(op);
		if (rpt == 0) continue;
		op  = edit_op_operation(op);
		switch (op)
			{
			case editopSub:
				stop = s1 + rpt;
				while (s1 < stop)
					similarity += scoring->sub[*(s1++)][*(s2++)];
				pos1 += rpt;
				pos2 += rpt;
				//fprintf (stderr, "match " unsposFmt " -> " scoreFmt "\n", rpt, similarity - s);
				break;
			case editopIns:
				similarity -= scoring->gapOpen + (rpt * scoring->gapExtend);
				s2 += rpt;
				//fprintf (stderr, "insert " unsposFmt " -> " scoreFmt "\n", rpt, similarity - s);
				break;
			case editopDel:
				similarity -= scoring->gapOpen + (rpt * scoring->gapExtend);
				s1 += rpt;
				//fprintf (stderr, "delete " unsposFmt " -> " scoreFmt "\n", rpt, similarity - s);
				break;
			}

		}

	return similarity;
	}

//----------
//
// count_paired_bases--
//	Count the number of paired bases in an alignment.
//
//----------
//
// Arguments:
//	galign*	mp:	The MSP that was used to anchor the alignment.
//
// Returns:
//	The number of paired bases.
//
//----------

static u64 count_paired_bases
   (galign*		mp)
	{
	aliseg*		bp;
	u64			pairedBases;

	pairedBases = 0;
	for (bp=mp->firstSeg ; bp!=NULL ; bp=bp->nextSeg)
		{ if (bp->type == diagSeg) pairedBases += bp->e1+1 - bp->b1; }

	return pairedBases;
	}

//----------
//
// warn_for_paired_bases_limit--
//	Tell the user that this query exceeded the limit for paired bases.
//
//----------
//
// Arguments:
//	seq*	seq2:				The query sequence we were aligning.
//	u64		maxPairedBases:		(same meaning as for gapped_extend)
//	int		overlyPairedKeep:	(same meaning as for gapped_extend)
//
// Returns:
//	(nothing)
//
//----------

static void warn_for_paired_bases_limit
   (seq*			seq2,
	u64				maxPairedBases,
	int				overlyPairedKeep)
	{
	static int		firstReport = true;
	seqpartition*	sp2 = &seq2->partition;
	char*			name2;
	char			strand;

	if      (sp2->p != NULL)     name2 = "seq2"; // (seq2 is partitioned)
	else if (seq2->useFullNames) name2 = seq2->header;
	else                         name2 = seq2->shortHeader;

	strand = ((seq2->revCompFlags & rcf_rev) == 0)? '+' : '-';

	fprintf (stderr, "WARNING. Query %s (%c strand) contains more than %s paired bases.\n",
					 name2, strand, commatize(maxPairedBases));

	if (firstReport)
		{
		if (overlyPairedKeep)
			fprintf (stderr, "Any gapped alignments already found for this query/strand are reported but the\n"
                             "query/strand is not processed further.\n");
		else
			fprintf (stderr, "All gapped alignments for this query/strand are discarded and the query/strand\n"
                             "is not processed further.\n");
		firstReport = false;
		}
	}

//----------
//
// gapped_extend_zero_stats--
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

void gapped_extend_zero_stats
   (void)
	{
	dbg_timing_set_stat (numExtensions, 0);

#ifdef collect_stats

	// set 'em en masse to zero

	memset (&gappedExtendStats, 0, sizeof(gappedExtendStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	gapped_extend_set_stat (totalPeakScore, 0);

#endif // collect_stats
	}

//----------
//
// gapped_extend_show_stats--
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print the stats to.
//
// Returns:
//	(nothing)
//
//----------

void gapped_extend_show_stats
   (arg_dont_complain(FILE* f))
	{
	dbg_timing_report_stat (numExtensions, "gapped extensions");

#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, " number of anchors: %s\n", commatize(gappedExtendStats.numAnchors));
	fprintf (f, "  anchors >= %2d bp: %s\n", anchorPeakLen, commatize(gappedExtendStats.numPeaks));
	if (gappedExtendStats.numPeaks > 0)
		fprintf (f, "average peak score: %.1f\n", ((float)gappedExtendStats.totalPeakScore) / gappedExtendStats.numPeaks);
	fprintf (f, "  anchors extended: %s\n", commatize(gappedExtendStats.numAnchorsExtended));
	fprintf (f, " gapped extensions: %s\n", commatize(gappedExtendStats.numExtensions));
	fprintf (f, "  DP cells visited: %s\n", commatize(gappedExtendStats.dpCellsVisited));
	if (gappedExtendStats.numExtensions > 0)
		fprintf (f, "DP cells/extension: %s\n", commatize((2*gappedExtendStats.dpCellsVisited+gappedExtendStats.numExtensions)/(2*gappedExtendStats.numExtensions)));
	fprintf (f, "       max DP rows: %s\n", commatize(gappedExtendStats.maxDpRows));
	fprintf (f, "    max DP columns: %s\n", commatize(gappedExtendStats.maxDpColumns));

	if (gappedExtendStats.zallocCallsA != 0)
		fprintf (f, "    zalloc calls A: %s (%s bytes per)\n",
		            commatize(gappedExtendStats.zallocCallsA),
		            commatize((gappedExtendStats.zallocTotalA + gappedExtendStats.zallocCallsA/2) / gappedExtendStats.zallocCallsA));
		fprintf (f, "    zalloc total A: %s\n", commatize(gappedExtendStats.zallocTotalA));

	if (gappedExtendStats.zallocCallsB != 0)
		fprintf (f, "    zalloc calls B: %s (%s bytes per)\n",
		            commatize(gappedExtendStats.zallocCallsB),
		            commatize((gappedExtendStats.zallocTotalB + gappedExtendStats.zallocCallsB/2) / gappedExtendStats.zallocCallsB));
		fprintf (f, "    zalloc total B: %s\n", commatize(gappedExtendStats.zallocTotalB));

	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void gapped_extend_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "num_anchors=%" PRId64 "\n", gappedExtendStats.numAnchors);
#endif // collect_stats
	}

