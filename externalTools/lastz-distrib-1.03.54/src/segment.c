//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: segment.c
//
//----------
//
// segment--
//	Support for lists of matched segments of DNA sequences.
//
// This module implements the collection of a table of segments of a pair of
// DNA sequences.  A segment consists of an interval in each sequence (of the
// same length) and a score.  The actual DNA content of the interval is of
// no importance to this module-- all that is important is represented by a
// segment's length and score.
//
// The caller can specify a limit on the total length of segments in the table.
// The intent is to give the user the ability to collect alignments covering
// some meaningful portion of the sequences, even though she does not know what
// score threshold would accomplish this.  Without this capability, the user is
// forced to do some preliminary runs to try to estimate the score distribution,
// or set the score threshold so low that the process is overwhelmed by false
// positives.
//
// The limit (coverageLimit) is the number of bases of coverage that the caller
// wishes to allow.  However, we allow the sum of all the segment lengths to
// exceed this, by no more than one segment.  For example, if the caller
// specifies a limit of 10K bp, and if there is at least 10K bp of homologous
// segments in the two sequences, we will keep at least 10K bp in this table,
// but not much more.
//
//----------
//
// $$$ Caveats
//
// The coverageLimit stuff doesn't consider the possibility that the incoming
// segments could contain overlaps.  As of this writing (Mar/2008) we will only
// have overlaps if the caller is using a twin-hit seed.
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
#include "build_options.h"		// build options
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "diag_hash.h"			// diagonals hashing stuff

#define  segment_owner			// (make this the owner of its globals)
#include "segment.h"			// interface to this module

// miscellany

#define segtable_bytes(size) (sizeof(segtable) + (((size)-1)*sizeof(segment)))

//#define debugBinaryHeap

//----------
//
// prototypes for private functions
//
//----------

static void remove_root       (segtable* st);
static void record_tie_scores (segtable* st);
static int  record_tie_score  (segtable* st, int ix);

#ifdef debugBinaryHeap
static void validate_heap     (segtable* st, char* msg);
#endif // debugBinaryHeap

//----------
//
// new_segment_table--
//	Allocate an empty segment table.
//
//----------
//
// Arguments:
//	u32		size:			The number of entries to provide for.  This can be
//							.. increased later.
//	unspos	coverageLimit:	Limit on the total lengths of the segments in the
//							.. table (see discussion in file header above);
//							.. zero indicates no limit.
//
// Returns:
//	A pointer to a newly allocated segment table;  failures result in program
//	fatality.  The caller must eventually dispose of the table, with a call to
//	free_segment_table().
//
//----------

segtable* new_segment_table
   (u32			size,
	unspos		coverageLimit)
	{
	segtable*	st;
	size_t		bytesNeeded;

	// sanity check

	if (size < 1)
		suicidef ("in new_segment_table(), size can't be %d", size);

	// allocate

	bytesNeeded = segtable_bytes (size);
#if (SIZE_MAX > mallocLimit)
	if (bytesNeeded > mallocLimit) goto overflow;
#endif // overflow possible
	st = (segtable*) malloc_or_die ("new_segment_table", bytesNeeded);

	// initialize

	st->size          = size;
	st->len           = 0;
	st->haveScores    = false;
	st->coverageLimit = coverageLimit;
	st->coverage      = 0;
	st->lowScore      = worstPossibleScore;

	return st;

// failure exits

#if (SIZE_MAX > mallocLimit)
overflow:
	suicidef ("internal error, in new_segment_table()\n"
	          "table size (%s) exceeds allocation limit of %s;",
	          commatize(bytesNeeded),
	          commatize(mallocLimit));
	return NULL; // (doesn't get here)
#endif // overflow possible
	}

//----------
//
// empty_segment_table--
//	"Erase" all segments from a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:	The segment table to empty.
//
// Returns:
//	(nothing)
//
//----------

void empty_segment_table
   (segtable*	st)
	{
	st->len        = 0;
	st->haveScores = false;
	st->coverage   = 0;
	st->lowScore   = worstPossibleScore;
	}

//----------
//
// limit_segment_table--
//	Change a segment table's total length limit.
//
//----------
//
// Arguments:
//	segtable*	st:				The segment table to modify.
//	unspos		coverageLimit:	Limit on the total lengths of the segments in
//								.. the table;  zero indicates no limit.
//
// Returns:
//	(nothing)
//
//----------

void limit_segment_table
   (segtable*	st,
	unspos		coverageLimit)
	{
	int			ix, newLen;
	possum		cov, newCov;
	score		prevScore, newLow;
	segment*	seg, *tail;
	segment		tempSeg;

	st->coverageLimit = coverageLimit;

	// if this leaves us below the limit, we're done

	if ((st->coverageLimit == 0)
	 || (st->coverage < st->coverageLimit)
	 || (st->len < 2))
	 	return;

	// otherwise, may need to reduce;  begin by sorting segments by decreasing
	// score

	sort_segments (st, qSegmentsByDecreasingScore);

	// find the score that reduces coverage to no more than the limit

	newLen = st->len;
	newCov = st->coverage;
	newLow = st->lowScore;

	ix = st->len;
	seg = &st->seg[--ix];

	cov       = st->coverage - seg->length;
	prevScore = seg->s;

	while (ix > 0)
		{
		seg = &st->seg[--ix];

		// if this is the same as the segment after it, just keep looking (we
		// treat segments with the same score as one unbreakable entity)

		if (seg->s == prevScore)
			{ cov -= seg->length;  continue; }

		// if removing the segments above this one would be enough to get us
		// within the limit, we are finished;  note that what we remove will be
		// the segments above the (score-tied) segments above this one

		if (cov < st->coverageLimit)
			break;

		// potentially, we will remove the segments above this one, so record
		// that information and keep going

		newLen = ix+1;
		newCov = cov;
		newLow = seg->s;

		cov -= seg->length;
		prevScore = seg->s;
		}

	st->len      = newLen;
	st->coverage = newCov;
	st->lowScore = newLow;

	// now make the list a proper min-heap by simply reversing its order (which
	// orders it by increasing score), and add the tied-score information

	seg  = &st->seg[0];
	tail = &st->seg[st->len-1];
	while (seg < tail)
		{
		tempSeg   = *seg;
		*(seg++)  = *tail;
		*(tail--) = tempSeg;
		}

	record_tie_scores (st);
	}

//----------
//
// free_segment_table--
//	Deallocate a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:	The segment table to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_segment_table
   (segtable*	st)
	{
	free_if_valid ("free_segment_table", st);
	}

//----------
//
// read_segment_table--
//	Read a segment table from a file.
//
// The file consists of segments, one per line, like the ones below.  Note that
// this is NOT the same format as is written by write_segments or dump_segments
// (though it is similar to that written by write_segments).  There are three
// columns for the target-- name, start, end.  These are followed by three
// columns for the query (with the same meaning), the query strand, and an
// optional score.
// 
//		E88BQJZ01A3EQH 151 225  E86HODY01D81JM 14  88 +  6875
//		E88BQJZ01D4L6V  26 100  E86HODY01D81JM 10  84 +  6808
//		E88BQJZ01EVLNU  19  93  E86HODY01D81JM  7  81 +  6842
//		E88BQJZ01CEBPD   8  81  E86HODY01D81JM  9  82 +  7108
//		E88BQJZ01BLO6X 132 205  E86HODY01D81JM 11  84 -  7339
//		E88BQJZ01A2W3P 162 214  E86HODY01D81JM  2  54 -  5024
//		E88BQJZ01A9395  62 136  E86HODY01A323K 18  92 +  7231
//		E88BQJZ01DNC74  18  82  E86HODY01A323K  2  66 +  6418
//		E88BQJZ01CTR26  83 167  E86HODY01ASA7F 19 103 +  8034
//		E88BQJZ01C2TAC  95 181  E86HODY01ASA7F 15 101 +  8272
//
// Start and end are origin one, closed (thus the interval "154 228" has length
// 75 and is preceded by 153 bases in its sequence).  Negative strand intervals
// are measured from the 5' end of the query's *negative* strand (e.g from the
// opposite end as that used for the positive strand).  All target intervals are
// on the positive strand.  The query interval length *must* match the target
// interval length.  Segments without scores are given the score of zero.
//
// Sequence names for the query *must* appear in the same order as they do in
// the query file (but see note below about partitioned sequences).  For a given
// query, all positive strand intervals must appear before any negative strand
// intervals.  Sequence names for the target may appear in any order, and are
// only meaningful for partitioned sequences (see below);  otherwise they are
// ignored.  Intervals with names not found in the target or query are not
// allowed.
//
// A * can be used as a generic sequence name, in those cases where sequence
// names are either unknown or of no importance.
//
// A "#" is considered a comment.  Anything following a "#" (on the same line)
// is ignored.  Blank lines are also ignored.
//
// Note: Partitioned sequences (enabled at the command line by the "multi" file
//       action) are internally treated as a single sequence.  In this case,
//       query names can appear in any order.  However, *all* positive strands
//       must appear before any negative strands.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to read segments from.
//	char*		fName:	The name of the file (used only for reporting errors);
//						.. may be NULL.
//	segtable*	st:		The segment table to fill.  If this is NULL, it
//						.. indicates that this is the final call for the file
//						.. (see note about final call below).
//	seq*		target:	The sequence being searched.
//	seq*		query:	The sequence(s) being searched for.
//
// Returns:
//	A pointer to a the segment table;  if there was room in st for all the
//	segments, this is the same as st;  otherwise, this is a pointer to a new
//	copy of st, and the previous st has been deallocated;  reallocation
//	failures result in program fatality.
//
//----------
//
// Notes:
//	(1)	After the final query, this routine should be called with st == NULL,
//		to verify that all segments have been read.
//
//	(2)	WARNING: Since the line buffer is used to save information between
//		calls, this routine can not be used to handle multiple files in the
//		same run of a program.
//
//	(3)	[this note is also referenced from resolve_chore_query]
//		Example of position arithmetic for queries on negative strand in
//		partitions.  Suppose we have this partitioning of the query (note that
//		startLoc is 1-based):
//
//				 sepBefore sepAfter startLoc  trueLen contig header
//			[ 0]        0      481       51     2000      1  query17
//			[ 1]      481      951      532     2000      1  query17
//			[ 2]      951     1481     1002     2000      1  query17
//			[ 3]     1481     1851     1532     2000      1  query17
//			[ 4]     1851     2332       51     2000      2  query20
//			[ 5]     2332     2802      532     2000      2  query20
//			[ 6]     2802     3332     1002     2000      2  query20
//			[ 7]     3332     3702     1532     2000      2  query20
//			[ 8]     3702        0        0        0
//
//		and this segment:
//
//			target22 1392 1463  query17 272 343  -
//
//		Since the segment is counted along the reverse strand and the full
//		query is 2000 bp long, the interval 272..343 corresponds to the interval
//		1658..1729 along the forward strand (the segment file uses 1-based
//		intervals, so 2001-343 gives 1658, and 2001-272 gives 1729).  The reason
//		we need the forward strand is so we can identify the partition that
//		contains the interval.
//
//		Now, when the query is reverse-complemented, each partition is reverse-
//		complemented separately in place, so that the partition table remains
//		unchanged.  This means that (using 1-based intervals) partition 3
//		corresponds to forward-strand query17 1532..1900 and reverse-strand
//		query17 101..469.  It is in memory at 1482..1850, in reverse (where
//		"memory" means seq->v[]).
//
//		For this segment, we want the interval 1658..1729.  Along a number line,
//		this looks like this (1-based intervals):
//
//			reverse-strand sequence:  101............... 272..343 ......... 469
//			forward-strand sequence: 1900...............1729..1658.........1532
//			memory:                  1482...............1654..1725.........1850
//
//		So our interval is in memory at 1654..1725.  In zero-based indexes this
//		is 1653..1725.
//
//----------

segtable* read_segment_table
   (FILE*			f,
	char*			fName,
	segtable*		st,
	seq*			target,
	seq*			query)
	{
	// (parsing variables, preserve between calls)
	static char		line[1024];
	static int		pendingRewind = false;
	static int		pendingLine   = false;
	static int		pendingFirstAfterRewind = false;
	static int		lineNum = 0;
	static char*	tName, *qName, *tPrevName, *tPartName;
	static unspos	tStart, tEnd, qStart, qEnd;
	static char		qStrand;
	static score	s;
	static int		prevQueryWasPartitioned = false;
	// (normal local variables)
	int				firstAfterRewind;
	unspos			tSeqStart, tSeqEnd, qSeqStart, qSeqEnd;
	unspos			tSegStart, qSegStart, segLen;
	int				len, missingEol;
	char*			scan, *field;
	int				numItems, charsUsed;
	unspos			tOffset, tLen, qOffset, qLen, qTrue, qNegStart, qNegEnd;
	char*			queryName   = "";
	char			queryStrand = '+';
	seqpartition*	tSp = &target->partition;
	seqpartition*	qSp = &query->partition;
	partition*		tNamePart, *tPart, *qNamePart, *qPart;
	u32				tIx;
	int				err;

	if (fName == NULL) fName = "(filename not known)";

	if (st != NULL)
		{
		if (qSp->p != NULL)					// query is partitioned
			{
			queryName = "(partitioned query)";
			prevQueryWasPartitioned = true;
			}
		else								// query is not partitioned
			{
			queryName = (query->useFullNames)? query->header
											 : query->shortHeader;
			prevQueryWasPartitioned = false;
			}
		queryStrand = ((query->revCompFlags & rcf_rev) != 0)? '-' : '+';
		}

	// read the segments for this query/strand

	missingEol = false;

	firstAfterRewind = false;
	if (pendingRewind)
		{
		if (st != NULL)
			{
			err = fseek (f, 0, SEEK_SET);
			if (err != 0) goto rewind_failed;
			lineNum = 0;
			}
		pendingRewind = false;
		firstAfterRewind = true;
		}

	if ((st == NULL) && (pendingLine) && (pendingFirstAfterRewind))
		pendingLine = pendingFirstAfterRewind = false;

	while (true)
		{
		// get the next line, if we need one;  we also check for lines getting
		// split by fgets (the final line in the file might not have a newline,
		// but no internal lines can be that way)

		if (pendingLine)
			{
			pendingLine = false;
			firstAfterRewind = pendingFirstAfterRewind;
			goto parsing_finished;
			}
		else
			{
			if (fgets (line, sizeof(line), f) == NULL) break;
			lineNum++;

			if (missingEol)
				goto split_line;

			len = strlen(line);
			if (len == 0) continue;
			missingEol = (line[len-1] != '\n');
			}

		// trim blanks, end of line, and comments, and ignore blank lines

		if (line[len-1] == '\n') line[--len] = 0;

		field = strchr (line, '#');
		if (field != NULL) *field = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// see if this is a "rewind" command
		// $$$ we should make sure there's nothing left in the file

		if (strcmp (line, "rewind") == 0)
			{ pendingRewind = true;  break; }

		// parse the line

		if (segment_dbgAnchorParsing)
			fprintf (stderr, "segment line: \"%s\"\n", line);

		scan = line;

		if (*scan == 0) goto not_enough_fields;
		tName = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);

		if (*scan == 0) goto not_enough_fields;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		charsUsed = -1;
		numItems = sscanf (field, unsposFmtScanf "%n", &tStart, &charsUsed);
		if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;

		if (*scan == 0) goto not_enough_fields;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		charsUsed = -1;
		numItems = sscanf (field, unsposFmtScanf "%n", &tEnd, &charsUsed);
		if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
		if (tEnd < tStart) goto bad_target_interval;

		if (*scan == 0) goto not_enough_fields;
		qName = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);

		if (*scan == 0) goto not_enough_fields;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		charsUsed = -1;
		numItems = sscanf (field, unsposFmtScanf "%n", &qStart, &charsUsed);
		if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;

		if (*scan == 0) goto not_enough_fields;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		charsUsed = -1;
		numItems = sscanf (field, unsposFmtScanf "%n", &qEnd, &charsUsed);
		if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;
		if (qEnd < qStart) goto bad_query_interval;
		if (qEnd-qStart != tEnd - tStart) goto interval_length_mismatch;

		if (*scan == 0) goto not_enough_fields;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		if (strlen(field) != 1) goto bad_field;
		qStrand = *field;
		if ((qStrand != '+') && (qStrand != '-')) goto bad_strand;

		s = 0;
		if (*scan == 0) goto new_parse_finished;
		field = scan;
		scan = skip_darkspace (scan);
		*(scan++) = 0;
		scan = skip_whitespace (scan);
		charsUsed = -1;
		numItems = sscanf (field, scoreFmtScanf "%n", &s, &charsUsed);
		if ((numItems != 1) || (((u32)charsUsed) != strlen(field))) goto bad_field;

		if (*scan != 0) goto too_many_fields;

	new_parse_finished:

		if (segment_dbgAnchorParsing)
			fprintf (stderr, "  (line parsed)\n");

	parsing_finished:

		// it's a *syntactically* valid segment, but if we aren't accepting any
		// more segments, this is a failure

		if (st == NULL)
			goto extra_segments;

		// resolve query interval 
		// 	qSeqStart is in-file position of start of the resident piece-of-sequence (origin zero)
		//  qOffset   is index into query v[] of the start of the resident piece-of-sequence

		if (qStrand != queryStrand)
			{
			if (segment_dbgAnchorParsing)
				fprintf (stderr, "  (query strand mismatch-- %s%c vs %s%c)\n",
				                     qName, qStrand, queryName, queryStrand);
			goto query_name_mismatch;
			}

		if (qSp->p == NULL)					// query is not partitioned
			{
			if (strcmp (qName, "*") != 0)
				{
				if ((queryName != NULL)
				 && (queryName[0] != 0)
				 && (strcmp (qName, queryName) != 0))
				 	{
					if (segment_dbgAnchorParsing)
						fprintf (stderr, "  (query name mismatch--   %s%c vs %s%c)\n",
						                     qName, qStrand, queryName, queryStrand);
					goto query_name_mismatch;
					}
				}

			qSeqStart = query->startLoc - 1;
			qOffset   = 0;
			qLen      = query->len;
			qSeqEnd   = qSeqStart + qLen;
			if (qStrand == '-')
				{ // negative strand
				qTrue     = query->trueLen;
				qNegStart = qTrue - qSeqEnd;                    // (origin zero)
				qSeqEnd   = qTrue - qSeqStart;
				qSeqStart = qNegStart;
				}
			if (qStart <= qSeqStart) goto query_interval_before_start;
			if (qEnd   >  qSeqEnd)   goto query_interval_after_end;
			}
		else if (strcmp (qName, "*") == 0)	// query is partitioned and
			goto query_wild_card;			// .. name is wildcard
		else								// query is partitioned and
			{								// .. specific name is given
			qNamePart = lookup_named_partition (query, qName);
			if (qNamePart == NULL) goto bad_query_name;
			if (qStrand != '-')
				{ // positive strand
				qPart = lookup_partition_seq_pos (query, qNamePart, qStart);
				if (qPart == NULL) goto bad_query_position_forward;
				qSeqStart = qPart->startLoc - 1;
				qOffset   = qPart->sepBefore + 1;
				qLen      = qPart->sepAfter - qOffset;
				qSeqEnd   = qSeqStart + qLen;
				if (qEnd > qSeqEnd) goto query_interval_after_end;
				}
			else
				{ // negative strand, see note 3 above
				qTrue     = qNamePart->trueLen;
				qNegStart = qTrue+1 - qEnd;                     // (origin one)
				//иии (qNegEnd is never read, investigate)
				qNegEnd   = qNegStart + qEnd - qStart;
                                // for -Wunused-but-set-variable
                                (void) qNegEnd;
				qPart = lookup_partition_seq_pos (query, qNamePart, qNegStart);
				if (qPart == NULL) goto bad_query_position_reverse;
				qOffset   = qPart->sepBefore + 1;
				qLen      = qPart->sepAfter - qOffset;
				qSeqEnd   = qTrue - (qPart->startLoc-1);
				qSeqStart = qSeqEnd - qLen;                     // (origin zero)
				if (qStart <= qSeqStart) goto query_interval_before_reversed_end;
				}
			}

		// resolve target interval 

		if (tSp->p == NULL)					// target is not partitioned
			{
			tSeqStart = target->startLoc - 1;
			tOffset   = 0;
			tLen      = target->len;
			tSeqEnd   = tSeqStart + tLen;
			if (tStart <= tSeqStart) goto target_interval_before_start;
			if (tEnd   >  tSeqEnd)   goto target_interval_after_end;
			}
		else if (strcmp (tName, "*") == 0)	// target is partitioned and
			goto target_wild_card;			// .. name is wildcard
		else								// target is partitioned and
			{								// .. specific name is given
			tNamePart = lookup_named_partition (target, tName);
			if (tNamePart == NULL) goto bad_target_name;
			tPart = lookup_partition_seq_pos (target, tNamePart, tStart);
			if (tPart == NULL) goto bad_target_position;
			tSeqStart = tPart->startLoc - 1;
			tOffset   = tPart->sepBefore + 1;
			tLen      = tPart->sepAfter - tOffset;
			tSeqEnd   = tSeqStart + tLen;
			if (tEnd > tSeqEnd) goto target_interval_after_end;
			}

		// (phew!) it's a valid segment, add it to the table;  note that we use
		// the strand as an id

		tSegStart = tOffset + (tStart-1)-tSeqStart;
		qSegStart = qOffset + (qStart-1)-qSeqStart;
		segLen    = tEnd+1 - tStart;

		if (segment_dbgAnchorParsing)
			fprintf (stderr, "  --> adding segment " unsposFmt ".." unsposFmt
			                 "/" unsposFmt ".." unsposFmt "%c\n",
							  tSegStart, tSegStart + segLen,
							  qSegStart, qSegStart + segLen, qStrand);

		st = add_segment (st, tSegStart, qSegStart, segLen,
						  s, /*id*/ qStrand);
		continue;

		// the given target name is a wild card, and target is partitioned, so
		// we have to add a segment for every sequence in the target;  note that
		// some partitions may have the same name (e.g. if the user is using a
		// separator), so we have to look for the correct partition on a name-
		// by-name basis

	target_wild_card:

		tPrevName = "";
		for (tIx=0 ; tIx<tSp->len ; tIx++)
			{
			tNamePart = &tSp->p[tIx];
			tPartName = &tSp->pool[tNamePart->header];
			if (strcmp (tPartName, tPrevName) == 0)
				continue;  // (this partititon has same name as previous)
			tPrevName = tPartName;

			tPart = lookup_partition_seq_pos (target, tNamePart, tStart);
			if (tPart == NULL) goto bad_target_position;
			tSeqStart = tPart->startLoc - 1;
			tOffset   = tPart->sepBefore + 1;
			tLen      = tPart->sepAfter - tOffset;
			tSeqEnd   = tSeqStart + tLen;
			if (tEnd > tSeqEnd) goto target_interval_after_end;

			// add the segment to the table;  note that we use the strand as an
			// id

			tSegStart = tOffset + (tStart-1)-tSeqStart;
			qSegStart = qOffset + (qStart-1)-qSeqStart;
			segLen    = tEnd+1 - tStart;

			if (segment_dbgAnchorParsing)
				fprintf (stderr, "  --> adding segment " unsposFmt ".." unsposFmt
								 "/" unsposFmt ".." unsposFmt "%c\n",
								  tSegStart, tSegStart + segLen,
								  qSegStart, qSegStart + segLen, qStrand);

			st = add_segment (st, tSegStart, qSegStart, segLen,
							  s, /*id*/ qStrand);
			}
		continue;

		// interval name or strand did not match query;  this marks the end of
		// the list;  otherwise, we need to keep looking

	query_name_mismatch:
		pendingLine = true;
		pendingFirstAfterRewind = firstAfterRewind;
		break;
		}

	// success

	if (st != NULL) st->haveScores = true;
	return st;

	//////////
	// failure exits
	//////////

rewind_failed:
	suicidef ("failed to rewind segments file\n"
			  "in read_segment_table for %s, index fseek(0) returned %d",
			  fName, err);
	return NULL;

split_line:
	suicidef ("line is too long (%s: line %d)", fName, lineNum-1);
	return NULL;

not_enough_fields:
	suicidef ("line has too few fields (%s: line %d)", fName, lineNum);
	return NULL;

too_many_fields:
	suicidef ("line has too many fields (%s: line %d)", fName, lineNum);
	return NULL;

bad_field:
	suicidef ("bad field (%s: line %d, %s)", fName, lineNum, field);
	return NULL;

bad_target_interval:
	suicidef ("bad target interval (%s: line %d, " unsposFmt ">" unsposFmt ")",
	          fName, lineNum, tStart, tEnd);
	return NULL;

target_interval_before_start:
	suicidef ("target interval out of range (%s: line %d, " unsposFmt "<" unsposFmt ")",
	          fName, lineNum, tStart, tSeqStart+1);
	return NULL;

target_interval_after_end:
	// $$$ this should be more informative when a separator has been crossed
	suicidef ("target interval out of range (%s: line %d, " unsposFmt ">" unsposFmt ")",
	          fName, lineNum, tEnd, tSeqEnd);
	return NULL;

bad_query_interval:
	suicidef ("bad query interval (%s: line %d, " unsposFmt ">" unsposFmt ")",
	          fName, lineNum, qStart, qEnd);
	return NULL;

query_interval_before_start:
	if (qStrand == '-')
		suicidef ("query interval out of range (%s: line %d, " unsposFmt "<" unsposFmt")"
		          "\nminus strand subrange is " unsposFmt ".." unsposFmt,
		          fName, lineNum, qStart, qSeqStart+1, qSeqStart+1, qSeqEnd);
	else
		suicidef ("query interval out of range (%s: line %d, " unsposFmt "<" unsposFmt ")",
		          fName, lineNum, qStart, qSeqStart+1);
	return NULL;

query_interval_after_end:
	// $$$ this should be more informative when a separator has been crossed
	if (qStrand == '-')
		suicidef ("query interval out of range (%s: line %d, " unsposFmt ">" unsposFmt")"
		          "\nminus strand subrange is " unsposFmt ".." unsposFmt,
		          fName, lineNum, qEnd, qSeqEnd, qSeqStart+1, qSeqEnd);
	else
		suicidef ("query interval out of range (%s: line %d, " unsposFmt ">" unsposFmt ")",
		          fName, lineNum, qEnd, qSeqEnd);
	return NULL;

query_interval_before_reversed_end:
		suicidef ("query interval out of range (%s: line %d, " unsposFmt "<" unsposFmt ")"
		          "\nminus strand subrange is " unsposFmt ".." unsposFmt,
		          fName, lineNum, qStart, qSeqStart+1, qSeqStart+1, qSeqEnd);
	return NULL;

interval_length_mismatch:
	suicidef ("intervals have different lengths (%s: line %d, " unsposFmt " vs. " unsposFmt ")",
	          fName, lineNum, tEnd+1-tStart, qEnd+1-qStart);
	return NULL;

bad_strand:
	suicidef ("bad strand (%s: line %d, %c)",
	          fName, lineNum, qStrand);
	return NULL;

bad_target_name:
	suicidef ("bad target sequence name (%s: line %d, %s)",
	          fName, lineNum, tName);
	return NULL;

bad_target_position:
	suicidef ("bad target sequence name/position (%s: line %d, %s:" unsposFmt ")",
	          fName, lineNum, tName, tStart);
	return NULL;

bad_query_name:
	suicidef ("bad query sequence name (%s: line %d, %s)",
	          fName, lineNum, qName);
	return NULL;

bad_query_position_forward:
	suicidef ("bad query sequence name/position (%s: line %d, %s:" unsposFmt ")",
	          fName, lineNum, qName, qStart);
	return NULL;

bad_query_position_reverse:
	suicidef ("bad query sequence name/position (%s: line %d, %s:" unsposFmt ")"
			  "\n" unsposFmt " on forward strand",
			  fName, lineNum, qName, qStart, qNegStart);
	return NULL;

query_wild_card:
	suicidef ("bad query sequence name (%s: line %d, %s)\n"
	          "wildcard segment name (*) is not supported for queries with [multi]",
	          fName, lineNum, qName);
	return NULL;

extra_segments:
	if (prevQueryWasPartitioned)
		suicidef ("extra segments in file (%s: line %d, %s/%s%c)\n"
				  "(for this usage all + strand segments must appear before all - strand segments)",
				  fName, lineNum, tName, qName, qStrand);
	else
		suicidef ("extra segments in file (%s: line %d, %s/%s%c)\n"
				  "(for this usage segments must appear in the same order as the query file, with\n"
				  "all + strand segments before all - strand segments for each query)",
				  fName, lineNum, tName, qName, qStrand);
	return NULL;
	}

//----------
//
// add_segment--
//	Add a segment to a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:				The segment table to add to.
//	pos1, pos2, length, s, id:	The segment to add.  (Positions are origin-zero).
//
// Returns:
//	A pointer to a the segment table;  if there was room in st for the segment,
//	this is the same as st;  otherwise, this is a pointer to a new copy of st,
//	and the previous st has been deallocated;  reallocation failures result in
//	program fatality.
//
//----------

segtable* add_segment
   (segtable*	st,
	unspos		pos1,
	unspos		pos2,
	unspos		length,
	score		s,
	int			id)
	{
	u32			newSize;
	size_t		bytesNeeded;
	segment*	seg, *parent;
	segment		tempSeg;
	int			ix, pIx;
	int			tied, stopped;

//	fprintf (stderr, "add " unsposSlashSFmt " " unsposFmt " " scoreFmtSimple "; id %d\n",
//	               pos1+1, "+",
//	               pos2+1, ((id & rcf_rev) != 0)? "-" : "+",
//	               length, s, id);

	//////////
	// add the segment to the table, enlarging the table if needed, but
	// discarding the segment if it is low-scoring and the table has met its
	// coverage limit
	//////////

	// if the table is already full and this segment scores less than the
	// lowest score in the table, discard it

	if ((st->coverageLimit != 0)
	 && (st->coverage >= st->coverageLimit)
	 && (st->len > 0)
	 && (s < st->lowScore))
		return st;

	// if there's no room for the new segment, re-allocate

	if (st->len >= st->size)
		{
		newSize     = st->size + 100 + (st->size / 3);
		bytesNeeded = segtable_bytes (newSize);
#if (SIZE_MAX > mallocLimit)
		if (bytesNeeded > mallocLimit) goto overflow;
#endif // overflow possible
		st = (segtable*) realloc_or_die ("add_segment", st, bytesNeeded);
		st->size = newSize;
		}

	// add the segment, by appending it at the end

	seg = &st->seg[st->len++];
	seg->pos1     = pos1;
	seg->pos2     = pos2;
	seg->length   = length;
	seg->s        = s;
	seg->id       = id;
	seg->filter   = false;
	seg->scoreCov = (possum) length;

	st->coverage += length;
	if ((st->len == 1) || (s < st->lowScore)) st->lowScore = s;

	//////////
	// handle the transition between the two table states
	//	below-the-coverage-limit:  table is kept as a simple list
	//	met-the-coverage-limit:    table is kept as a proper min-heap
	//////////

	// if this segment leaves us below the limit, we're done

	if ((st->coverageLimit == 0)
	 || (st->coverage < st->coverageLimit))
		return st;

	// if this is the first time we've reached the limit, sort the segments to
	// create a proper min-heap, and add the tied-score information
	// nota bene:  if we reach here, st->coverageLimit > 0 and
	//             st->coverage >= st->coverageLimit

	if (st->coverage - length < st->coverageLimit)
		{
		sort_segments (st, qSegmentsByIncreasingScore);
		record_tie_scores (st);
		#ifdef debugBinaryHeap
		fprintf       (stderr, "\nafter sort:\n");
		dump_segments (stderr, st, NULL, NULL);
		validate_heap (st, "after sort");
		#endif // debugBinaryHeap
		goto prune;
		}

	//////////
	// maintain the min-heap property
	//////////

	#ifdef debugBinaryHeap
	//fprintf       (stderr, "\nbefore percolation:\n");
	//dump_segments (stderr, st, NULL, NULL);
	#endif // debugBinaryHeap

	// the rest of the list is a proper min-heap, so percolate the new segment
	// up the tree, while maintaining the tied-score information
	// nota bene:  if we reach here, length >= 2

	tied = false;
	for (ix=st->len-1 ; ix>0 ; )
		{
		pIx    = (ix-1) / 2;
		seg    = &st->seg[ix];
		parent = &st->seg[pIx];

		if (seg->s >= parent->s)
			{ tied = (seg->s == parent->s);  break; }

		// swap this segment with its parent, and adjust old parent's tied-score
		// subheap

		tempSeg = *seg;  *seg = *parent;  *parent = tempSeg;
		record_tie_score (st, ix);

		ix = pIx;
		}

	record_tie_score (st, ix);

	// if the new segment tied an existing score, we must continue to percolate
	// the tied-score info up the tree

	if (tied)
		{
		stopped = false;
		for (ix=(ix-1)/2 ; ix>0 ; ix=(ix-1)/2)
			{
			if (!record_tie_score (st, ix))
				{ stopped = true;  break; }
			}
		if (!stopped) record_tie_score (st, 0);
		}

	#ifdef debugBinaryHeap
	fprintf       (stderr, "\nafter percolation:\n");
	dump_segments (stderr, st, NULL, NULL);
	validate_heap (st, "after percolation");
	#endif // debugBinaryHeap

	//////////
	// remove low-scoring segments
	//////////

prune:

	// if removing the minimum scoring subheap would bring us below the
	// limit, no pruning is necessary

	if (st->coverage - st->seg[0].scoreCov < st->coverageLimit)
		return st;

	// otherwise, we must remove subheaps as long as doing so leaves us at or
	// above the limit

	while (st->coverage - st->seg[0].scoreCov >= st->coverageLimit)
		{
		s = st->seg[0].s;
		while (st->seg[0].s == s)
			{
			remove_root (st);
			#ifdef debugBinaryHeap
			fprintf       (stderr, "\nafter a pruning:\n");
			dump_segments (stderr, st, NULL, NULL);
			validate_heap (st, "after pruning");
			#endif // debugBinaryHeap
			}
		}

	st->lowScore = st->seg[0].s;

	#ifdef debugBinaryHeap
	fprintf       (stderr, "\nafter pruning:\n");
	dump_segments (stderr, st, NULL, NULL);
	validate_heap (st, "after pruning");
	#endif // debugBinaryHeap

	return st;

// failure exits

#if (SIZE_MAX > mallocLimit)
#define suggestions " consider using lastz_m40,"                              \
                    " or setting max_malloc_index for a special build,"       \
                    " or raising scoring threshold (--hspthresh or --exact)," \
                    " or breaking your target sequence into smaller pieces"

overflow:
	suicidef ("in add_segment()\n"
	          "table size (%s for %s segments) exceeds allocation limit of %s;\n"
	          suggestions,
	          commatize(bytesNeeded),
	          commatize(newSize),
	          commatize(mallocLimit));
	return NULL; // (doesn't get here)
#endif // overflow possible
	}


static void remove_root (segtable* st)
	{
	segment*	seg, *detached, *child, *rgtChild;
	u32			ix, childIx, rgtIx;

	// remove the root node's coverage

	seg = &st->seg[0];
	st->coverage -= seg->length;

	if (st->len <= 1)
		{ st->len = 0;  return; }

	// detach the final segment;  conceptually we will consider this as the new
	// root, then percolate it down to a new position

	detached = &st->seg[--st->len];

	if (st->len == 1)
		{ *seg = *detached;  return; }

	// recompute tied-score info up the tree from the detached node
	// nota bene:  st->len > 1, so the initial value of ix, (st->len-1)/2 >= 0

	for (ix=(st->len-1)/2 ; ix>0 ; ix=(ix-1)/2)
		{ if (!record_tie_score (st, ix)) break; }

	// the detached node may violate the min-heap property;  percolate it down
	// the tree until the property is re-established

	ix = 0;
	while (true)
		{
		childIx = 2*ix+1;
		if (childIx >= st->len) break;	// (ix is a leaf position)
		child = &st->seg[childIx];

		rgtIx = childIx + 1;
		if (rgtIx < st->len)
			{
			rgtChild = &st->seg[rgtIx];
			if (rgtChild->s < child->s)
				{ childIx = rgtIx;  child = rgtChild; }
			}

		if (detached->s <= child->s)// (if we place the detached node here
			break;					//  .. we'll satisfy the min-heap property)

		st->seg[ix] = *child;		// move the child up the tree

		ix = childIx;
		}

	st->seg[ix] = *detached;		// copy the detached node into the tree

	// now reverse that path, up the tree, to recompute tied-score info

	for ( ; ix>0 ; ix=(ix-1)/2)
		record_tie_score (st, ix);
	record_tie_score (st, 0);
	}

//----------
//
// record_tie_scores--
//	Record tied-score information in a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:		The segment table.
//
// Returns:
//	(nothing)
//
//----------

static void record_tie_scores
   (segtable*	st)
	{
	int			ix;

	// determine the coverage of all equal-scoring subheaps;  we start at the
	// tail and percolate coverage up the tree whenever we have a tie;  at the
	// end, seg->scoreCov will be equal to the sum of the lengths of the
	// subheap rooted at seg that consists of segments having the same score

	for (ix=st->len-1 ; ix>=0 ; ix--)
		record_tie_score (st, ix);
	}

//----------
//
// record_tie_score--
//	Record tied-score information for one node in a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:		The segment table.
//	int			ix:		The node index.
//
// Returns:
//	true  => score coverage for the node was changed by this operation;
//	false => score coverage for the node remained the same
//
//----------

static int record_tie_score
   (segtable*	st,
	int			ix)
	{
	possum		scoreCov;
	segment*	seg, *lftChild, *rgtChild;
	u32			lftIx, rgtIx;

	seg = &st->seg[ix];

	// scoreCov is the sum of the length of this segment, plus the left and
	// right subheaps *if* they exist and have the same score as this segment

	scoreCov = (possum) seg->length;

	lftIx = 2*ix+1;
	if (lftIx < st->len)
		{
		lftChild = &st->seg[lftIx];
		if (lftChild->s == seg->s)
			scoreCov += lftChild->scoreCov;

		rgtIx = lftIx + 1;
		if (rgtIx < st->len)
			{
			rgtChild = &st->seg[rgtIx];
			if (rgtChild->s == seg->s)
				scoreCov += rgtChild->scoreCov;
			}
		}

	if (scoreCov != seg->scoreCov)
		{ seg->scoreCov = scoreCov;  return true;  }

	return false;
	}

//----------
//
// split_segment_table--
//	Split a segment table into two tables.  All segments with a particular id
//	are left in the table, and all others are moved to a new table.
//
//----------
//
// Arguments:
//	segtable*	st:			The segment table to split.
//	int			id:			The id of segments to keep in the incoming table.
//	segtable**	leftovers:	The segment table in which to collect the leftovers.
//							.. This table must already exist.
//
// Returns:
//	(nothing)
//
//----------

void split_segment_table
   (segtable*	st,
	int			id,
	segtable**	_leftovers)
	{
	segtable*	leftovers = *_leftovers;
	possum		cov;
	score		lowScore;
	u32			srcIx, dstIx;
	segment*	seg;

	cov      = 0;
	lowScore = worstPossibleScore;

	//fprintf       (stderr, "incoming\n");
	//dump_segments (stderr, st, NULL, NULL);

	dstIx = 0;
	for (srcIx=0 ; srcIx<st->len ; srcIx++)
		{
		seg = &st->seg[srcIx];

		if (seg->id != id)
			{
			leftovers = add_segment (leftovers,
			                         seg->pos1, seg->pos2, seg->length,
			                         seg->s, seg->id);
			continue;
			}

		cov += seg->length;
		if ((dstIx == 0) || (seg->s < lowScore)) lowScore = seg->s;

		if (dstIx != srcIx) st->seg[dstIx] = *seg;
		dstIx++;
		}

	st->len      = dstIx;
	st->coverage = cov;
	st->lowScore = lowScore;

	//fprintf       (stderr, "\nid=%d\n", id);
	//dump_segments (stderr, st,        NULL, NULL);
	//fprintf       (stderr, "\nileftovers\n");
	//dump_segments (stderr, leftovers, NULL, NULL);
	//fprintf       (stderr, "\n");

	*_leftovers = leftovers;
	}

//----------
//
// score_segments--
//	Score every segment in a table, as the sum of the substitution scores along
//	the segment.
//
//----------
//
// Arguments:
//	segtable*	st:			The segment table to score.
//	seq*		seq1:		The sequence corresponding to the segments' pos1.
//	seq*		seq2:		The sequence corresponding to the segments' pos2.
//	scoreset*	scoring:	The scoring scheme;  usually this treats lowercase
//							.. letters as being 'bad'.
//
// Returns:
//	(nothing)
//
//----------

void score_segments
   (segtable*	st,
	seq*		seq1,
	seq*		seq2,
	scoreset*	scoring)
	{
	u32			ix;
	segment*	seg;
	u8*			s1, *s2;
	score		s;
	unspos		togo;

	for (ix=0,seg=st->seg ; ix<st->len ; ix++,seg++)
		{
		s1 = seq1->v + seg->pos1;
		s2 = seq2->v + seg->pos2;

		s = 0;
		for (togo=seg->length ; togo>0 ; togo--)
			s += scoring->sub[*(s1++)][*(s2++)];

		seg->s = s;
		}

	}

//----------
//
// sort_segments--
//	Sort a segment table.
// sort_some_segments--
//	Sort part of a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:		The segment table to add to.
//	u32			start:	(sort_some_segments only) Index into the table of the
//						.. first entry to be sorted.
//	u32			end:	(sort_some_segments only) Index into the table of the
//						.. first entry NOT to be sorted (i.e. the one after the
//						.. last index to be sorted).
//	int (*qCompare):	Comparison function to use (suitable for the standard
//						.. c function qsort).  For example, this could be
//						.. qSegmentsByPos2.
//
// Returns:
//	(nothing)
//
//----------

void sort_segments
   (segtable*	st,
	int (*qCompare) (const void* el1, const void* el2))
	{
	qsort (st->seg, st->len, sizeof(segment), qCompare);
	}


void sort_some_segments
   (segtable*	st,
	u32			start,
	u32			end,
	int (*qCompare) (const void* el1, const void* el2))
	{
	if (end <= start) return;    // (empty interval)
	if (end > st->len) return;   // (interval is beyond end of lest)
	qsort (st->seg+start, end-start, sizeof(segment), qCompare);
	}

//----------
//
// merge_segments--
//	Merge any overlapping segments in a segment table.
//
//----------
//
// Arguments:
//	segtable*	st:	The segment table to operate upon.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	Segments are considered to overlap if they are on the same diagonal and
//		share any positions.  Segments that adjoin, but do not overlap, are not
//		merged.
//	(2)	Any merged segment is given the maximum score of the segments being
//		merged.  If this is not appropriate, the caller will need to rescore
//		the segments afterwards.
//		$$$ The problem with making the caller rescore is that the caller has
//		$$$ .. no way to know which segments have been merged, thus *every*
//		$$$ .. segment would have to be rescored.  We could resolve this by
//		$$$ .. having the caller provide a callback routine to rescore a
//		$$$ .. segment.
//	(3) The min-heap property is NOT maintained.  We view this routine as a
//		post-processing step, and we do not expect additional segments to be
//		added to the table after this is called.
//
//----------

void merge_segments
   (segtable*	st)
	{
	u32			srcIx,   dstIx;
	segment*	srcSeg, *dstSeg;
	unspos		pos2, end2, srcPos2, srcEnd2;
	sgnpos		diag,       srcDiag;
	score		s,          srcS;

	// if we have fewer than two segments, there's nothing to merge

	if (st->len < 2) return;

	// sort segments by increasing pos2 along each diagonal

	sort_segments (st, qSegmentsByDiag);

	// start first segment

	srcSeg = st->seg;
	pos2   = srcSeg->pos2;
	diag   = diagNumber (srcSeg->pos1, pos2);
	end2   = pos2 + srcSeg->length;
	s      = srcSeg->s;
	srcSeg++;

	// scan segments, merging as needed;  note that any segment written to the
	// list is written to an earlier (or the same) index as the latest read

	dstIx = 0;
	for (srcIx=1 ; srcIx<st->len ; srcIx++,srcSeg++)
		{
		srcPos2 = srcSeg->pos2;
		srcDiag = diagNumber (srcSeg->pos1, srcPos2);
		srcEnd2 = srcPos2 + srcSeg->length;
		srcS    = srcSeg->s;

		if ((srcDiag == diag) && (srcPos2 < end2))
			{ // merge
			if (srcEnd2 > end2) end2 = srcEnd2;
			if (srcS    > s)    s    = srcS;
			continue;
			}

		// deposit the previous segment

		dstSeg = &st->seg[dstIx++];
		dstSeg->pos1   = (unspos) (diag + pos2);
		dstSeg->pos2   = pos2;
		dstSeg->length = end2 - pos2;
		dstSeg->s      = s;

		// start a new segment

		pos2 = srcPos2;
		diag = srcDiag;
		end2 = srcEnd2;
		s    = srcS;
		}

	// deposit the final segment

	dstSeg = &st->seg[dstIx++];
	dstSeg->pos1   = (unspos) (diag + pos2);
	dstSeg->pos2   = pos2;
	dstSeg->length = end2 - pos2;
	dstSeg->s      = s;

	// adjust the length of the list

	st->len = dstIx;
	}

//----------
//
// filter_marked_segments--
//	Remove from a table any segments marked for filtering.
//
//----------
//
// Arguments:
//	segtable*	st:	The segment table to operate upon.
//
// Returns:
//	(nothing)
//
//----------

void filter_marked_segments
   (segtable*	st)
	{
	segment*	srcSeg, *dstSeg;

	if (st      == NULL) return;
	if (st->seg == NULL) return;

	for (dstSeg=srcSeg=st->seg ; ((u32)(srcSeg-st->seg))<st->len ; srcSeg++)
		{
		if (srcSeg->filter) continue;
		if (srcSeg != dstSeg) *dstSeg = *srcSeg;
		dstSeg++;
		}

	st->len = dstSeg - st->seg;
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSegmentsByPos1--
//	Compare two segments by pos1, so that qsort will sort by increasing pos1.
//
//----------
//
// Arguments:
//	const void* (really segment*) segA:	Pointer to one segment.
//	const void* (really segment*) segB:	Pointer to another.
//
// Returns:
//	> 0 => segA is greater than segB.
//	= 0 => segA and segB are the same.
//	< 0 => segA is less than segB.
//
//----------

int qSegmentsByPos1
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;

	if      (segA->pos1   < segB->pos1)   return -1;
	else if (segA->pos1   > segB->pos1)   return  1;

	if      (segA->length < segB->length) return -1;
	else if (segA->length > segB->length) return  1;

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	if      (segA->s      < segB->s)      return -1;
	else if (segA->s      > segB->s)      return  1;

	return 0;
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSegmentsByPos2--
//	Compare two segments by pos2, so that qsort will sort by increasing pos2.
//
//----------
//
// Arguments:
//	const void* (really segment*) segA:	Pointer to one segment.
//	const void* (really segment*) segB:	Pointer to another.
//
// Returns:
//	> 0 => segA is greater than segB.
//	= 0 => segA and segB are the same.
//	< 0 => segA is less than segB.
//
//----------

int qSegmentsByPos2
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	if      (segA->length < segB->length) return -1;
	else if (segA->length > segB->length) return  1;

	if      (segA->pos1   < segB->pos1)   return -1;
	else if (segA->pos1   > segB->pos1)   return  1;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	if      (segA->s      < segB->s)      return -1;
	else if (segA->s      > segB->s)      return  1;

	return 0;
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSegmentsByDecreasingScore--
//	Compare two segments by score, so that qsort will sort by decreasing score.
// qSegmentsByIncreasingScore--
//	Compare two segments by score, so that qsort will sort by increasing score.
//
//----------
//
// Arguments:
//	const void* (really segment*) segA:	Pointer to one segment.
//	const void* (really segment*) segB:	Pointer to another.
//
// Returns:
//	(for qSegmentsByDecreasingScore)     (for qSegmentsByIncreasingScore)
//	> 0 => segA is greater than segB.    > 0 => segA is greater than segB.
//	= 0 => segA and segB are the same.   = 0 => segA and segB are the same.
//	< 0 => segA is less than segB.       < 0 => segA is less than segB.
//
//----------

int qSegmentsByDecreasingScore
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;

	if      (segA->s      < segB->s)      return  1;
	else if (segA->s      > segB->s)      return -1;

	if      (segA->length < segB->length) return -1; // if scores are equal we
	else if (segA->length > segB->length) return  1; // .. prefer shorter length

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	if      (segA->pos1   < segB->pos1)   return -1;
	else if (segA->pos1   > segB->pos1)   return  1;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	return 0;
	}


int qSegmentsByIncreasingScore
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;

	if      (segA->s      < segB->s)      return -1;
	else if (segA->s      > segB->s)      return  1;

	if      (segA->length < segB->length) return -1; // if scores are equal we
	else if (segA->length > segB->length) return  1; // .. prefer shorter length

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	if      (segA->pos1   < segB->pos1)   return -1;
	else if (segA->pos1   > segB->pos1)   return  1;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	return 0;
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSegmentsByDiag--
//	Compare two segments by diagonal, so that qsort will sort by increasing
//	diagonal, and by increasing pos2 along each diagonal. 
//
//----------
//
// Arguments:
//	const void* (really segment*) segA:	Pointer to one segment.
//	const void* (really segment*) segB:	Pointer to another.
//
// Returns:
//	> 0 => segA is greater than segB.
//	= 0 => segA and segB are the same.
//	< 0 => segA is less than segB.
//
//----------

int qSegmentsByDiag
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;
	sgnpos		diagA, diagB;

	// compare by diagonal

	diagA = diagNumber (segA->pos1, segA->pos2);
	diagB = diagNumber (segB->pos1, segB->pos2);

	if      (diagA < diagB) return -1;
	else if (diagA > diagB) return  1;

	// resort to tiebreakers

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	if      (segA->length < segB->length) return -1;
	else if (segA->length > segB->length) return  1;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	if      (segA->s      < segB->s)      return -1;
	else if (segA->s      > segB->s)      return  1;

	return 0;
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSegmentsById--
//	Compare two segments by id, so that qsort will sort by increasing id.
//
//----------
//
// Arguments:
//	const void* (really segment*) segA:	Pointer to one segment.
//	const void* (really segment*) segB:	Pointer to another.
//
// Returns:
//	> 0 => segA is greater than segB.
//	= 0 => segA and segB are the same.
//	< 0 => segA is less than segB.
//
//----------

int qSegmentsById
   (const void*	_segA,
	const void*	_segB)
	{
	segment*	segA = (segment*) _segA;
	segment*	segB = (segment*) _segB;

	if      (segA->id     < segB->id)     return -1;
	else if (segA->id     > segB->id)     return  1;

	if      (segA->s      < segB->s)      return -1;
	else if (segA->s      > segB->s)      return  1;

	if      (segA->pos1   < segB->pos1)   return -1;
	else if (segA->pos1   > segB->pos1)   return  1;

	if      (segA->length < segB->length) return -1;
	else if (segA->length > segB->length) return  1;

	if      (segA->pos2   < segB->pos2)   return -1;
	else if (segA->pos2   > segB->pos2)   return  1;

	return 0;
	}

//----------
//
// write_segments--
//	Write a segment table to a file, for debugging.
//
// The table will look something like what is shown below.  Intervals are
// origin-zero, half-open.
//
//	[0] apple 1126 1177 orange 411 461 + 0
//	[1] apple 1224 1267 orange 498 540 + 0
//	[2] apple 1262 1322 orange 562 621 + 0
//	[3] apple 1340 1377 orange 655 691 + 0
//	[4] apple 1471 1530 orange 807 865 + 0
//
//	[0] apple 576  624  orange 222 269 - 0
//	[1] apple 666  724  orange 326 383 - 0
//	[2] apple 856  930  orange 482 555 - 0
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	segtable*	st:			The segment table to dump.
//	seq*		target:		The target sequence the table relates to.
//	seq*		query:		The query sequence the table relates to.
//	int			withText:	true => add nucleotide text (for debug only)
//
// Returns:
//	(nothing)
//
//----------

void write_segments
   (FILE*			f,
	segtable*		st,
	seq*			target,
	seq*			query,
	int				withText)
	{
	u32				segIx;
	segment*		seg;
	static char*	tName, *qName;
	static unspos	tStart, qStart;
	seqpartition*	tSp = &target->partition;
	seqpartition*	qSp = &query->partition;
	partition*		tPart, *qPart;
	u8*				tV, *qV;
	u32				ix;

	for (segIx=0,seg=st->seg ; segIx<st->len ; segIx++,seg++)
		{
		tPart = NULL;
		if (tSp->p == NULL)		// target is not partitioned
			{
			tName = (target->useFullNames)? target->header : target->shortHeader;
			if ((tName == NULL) || (tName[0] == 0)) tName = "target";
			tStart = seg->pos1;
			}
		else					// target is partitioned
			{
			tPart = lookup_partition (target, seg->pos1);
			tName  = &tSp->pool[tPart->header];
			tStart = seg->pos1 - (tPart->sepBefore + 1);
			}

		qPart = NULL;
		if (qSp->p == NULL)		// query is not partitioned
			{
			qName = (query->useFullNames)? query->header : query->shortHeader;
			if ((qName == NULL) || (qName[0] == 0)) qName = "query";
			qStart = seg->pos2;
			}
		else					// query is partitioned
			{
			qPart = lookup_partition (query, seg->pos2);
			qName  = &qSp->pool[qPart->header];
			qStart = seg->pos2 - (qPart->sepBefore + 1);
			}

		fprintf (f, "[%d]", segIx);
		if (tPart == NULL)
			fprintf (f, " %s " unsposFmt " " unsposFmt,
			            tName, tStart, tStart+seg->length);
		else
			fprintf (f, " %s " unsposFmt "+" unsposFmt " " unsposFmt "+" unsposFmt,
			            tName, tPart->startLoc-1, tStart, tPart->startLoc-1, tStart+seg->length);
		if (qPart == NULL)
			fprintf (f, " %s " unsposFmt " " unsposFmt,
			            qName, qStart, qStart+seg->length);
		else
			fprintf (f, " %s " unsposFmt "+" unsposFmt " " unsposFmt "+" unsposFmt,
			            qName, qPart->startLoc-1, qStart, qPart->startLoc-1, qStart+seg->length);
		fprintf (f, " %c " scoreFmtSimple, seg->id, seg->s);

		if (withText)
			{
			tV = target->v + seg->pos1;
			qV = query->v  + seg->pos2;

			fprintf (f, " %lu:", (unsigned long) (tV - target->v));
			for (ix=0 ; ix<seg->length ; ix++)
				fprintf (f, "%c", dna_toprint(tV[ix]));
			fprintf (f, " %lu:", (unsigned long) (qV - query->v));
			for (ix=0 ; ix<seg->length ; ix++)
				fprintf (f, "%c", dna_toprint(qV[ix]));
			}

		fprintf (f, "\n");
		}
	}

//----------
//
// dump_segments--
//	Dump a segment table, for debugging.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	segtable*	st:			The segment table to dump.
//	char*		sym1,sym2:	Identifying strings to attach to pos1 and pos2,
//							.. respectively.  E.g. "+" or "-" for DNA strand.
//							.. These can be NULL, in which case we presume that
//							.. segment id fields are the second sequence's
//							.. revCompFlags, and we take the strand from that.
//
// Returns:
//	(nothing)
//
//----------

//#define show_separation

void dump_segments
   (FILE*		f,
	segtable*	st,
	char*		_sym1,
	char*		_sym2)
	{
	char*		sym1, *sym2;
	u32			ix;
	segment*	seg;
	sgnpos		diag;
#ifdef show_separation
	unspos		prevPos2 = 0;
	sgnpos		prevDiag = 0;
#endif // show_separation

	sym1 = _sym1;
	sym2 = _sym2;

	for (ix=0,seg=st->seg ; ix<st->len ; ix++,seg++)
		{
		if (_sym1 == NULL) sym1 = "+";
		if (_sym2 == NULL) sym2 = ((seg->id & rcf_rev) != 0)? "-" : "+";

		diag = diagNumber (seg->pos1, seg->pos2);
#ifdef show_separation
		if ((ix != 0) && (diag == prevDiag))
			fprintf (f, "[" unsposFmt "] " unsposSlashSFmt " %d " scoreFmtSimple
			            "; id %d, diag " sgnposFmt " * (" unsposFmt "), tied_length=" possumFmt "\n",
						ix, seg->pos1+1, sym1, seg->pos2+1, sym2, seg->length,
						seg->s, seg->id, diag, seg->pos2-prevPos2, seg->scoreCov);
		else
#endif // show_separation
			fprintf (f, "[%d] " unsposSlashSFmt " " unsposFmt " " scoreFmtSimple
			            "; id %d, diag " sgnposFmt ", tied_length=" possumFmt "\n",
						ix, seg->pos1+1, sym1, seg->pos2+1, sym2, seg->length,
						seg->s, seg->id, diag, seg->scoreCov);
#ifdef show_separation
		prevDiag = diag;
		prevPos2 = seg->pos2;
#endif // show_separation
		}

	}

//----------
//
// validate_heap--
//	Check whether a segment table is a valid min-heap, for debugging.
//
//----------
//
// Arguments:
//	segtable*	st:		The segment table to check.
//	char*		msg:	A message to output if validate fails.
//
// Returns:
//	(nothing;  failure causes program fatality)
//
//----------

#ifdef debugBinaryHeap

static void validate_heap
   (segtable*	st,
	char*		msg)
	{
	possum		scoreCov;
	segment*	seg, *lftChild, *rgtChild;
	int			ix, lftIx, rgtIx;

	if (st->coverage < st->coverageLimit)
		suicidef ("%s, below coverage limit", msg);

	for (ix=0,seg=st->seg ; ix<st->len ; ix++,seg++)
		{
		scoreCov = (possum) seg->length;

		lftIx = 2*ix+1;
		if (lftIx < st->len)
			{
			lftChild = &st->seg[lftIx];
			if (lftChild->s < seg->s)
				suicidef ("%s, node %d > node %d", msg, ix, lftIx);
			if (lftChild->s == seg->s)
				scoreCov += lftChild->scoreCov;

			rgtIx = lftIx + 1;
			if (rgtIx < st->len)
				{
				rgtChild = &st->seg[rgtIx];
				if (rgtChild->s < seg->s)
					suicidef ("%s, node %d > node %d", msg, ix, rgtIx);
				if (rgtChild->s == seg->s)
					scoreCov += rgtChild->scoreCov;
				}
			}

		if (scoreCov != seg->scoreCov)
			suicidef ("%s, node %d has bad score coverage", msg, ix);
		}

	}

#endif // debugBinaryHeap

