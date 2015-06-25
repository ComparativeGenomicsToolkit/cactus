//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: coverage_dist.c
//
//----------
//
// coverage_dist--
//	Support for collecting the query coverage distribution from alignments.
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
#include <ctype.h>				// standard C upper/lower stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

#define  coverage_dist_owner	// (make this the owner of its globals)
#include "coverage_dist.h"		// interface to this module

//----------
//
// filter_aligns_by_coverage--
//	Filter a list of alignments, removing any alignment that has coverage
//	percentage outside of a specified range.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	alignel*	alignList:		The list of alignments to operate upon.
//	float		minCoverage,	The range of query coverage that we will *keep*.
//				maxCoverage		.. These are values between 0 and 1.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
//
// Notes:
//	(1)	The denominator is the length of whichever sequence is shorter, and the
//		numerator is the length of the alignment in that sequence.
//	(2)	Memory for alignments that don't make the cut is deallocated here.
//	(3)	The returned list of alignments is in the same order as the incoming
//		.. list.
//	(4) Coverage is counted over all aligned bases, included those aligned to
//		.. gaps.
//
//----------

alignel* filter_aligns_by_coverage
   (seq*		seq1,
	seq*		seq2,
	alignel*	alignList,
	float		minCoverage,
	float		maxCoverage)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;
	float		minCov, maxCov;

	// process each alignment, collecting a list of those that are long enough

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next  = a->next;

		// if the alignment isn't long enough, skip it

		alignment_coverage (seq1, seq2, a, &numer, &denom);

		minCov = denom * minCoverage;
		maxCov = denom * maxCoverage;

		if ((numer < minCov) || (numer > maxCov))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_coverage a->script", a->script);
			free_if_valid ("filter_aligns_by_coverage a",         a);
			continue;
			}


		// this alignment is ok, add it to the end of the new list we're
		// building

		if (head == NULL) head = prev = a;
		             else { prev->next = a;  prev = a; }

		a->next = NULL;
		}

	return head;
	}

//----------
//
// alignment_coverage--
//	Measure an alignment's coverage percentage.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	alignel*	a:				The alignment to measure.
//	unspos*		numer, denom:	Place to return the coverage
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	The denominator is the length of whichever sequence is shorter, and the
//		numerator is the length of the segment.
//
//----------

void alignment_coverage
   (seq*			seq1,
	seq*			seq2,
	alignel*		a,
	unspos*			numer,
	unspos*			denom)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			seq1Len, seq2Len;

	if (!seq1->needTrueLen)
		suicidef ("internal error, in alignment_coverage, target does not have true lengths");

	if (!seq2->needTrueLen)
		suicidef ("internal error, in alignment_coverage, query does not have true lengths");

	// determine sequence lengths for this alignment

	if (sp1->p == NULL)
		seq1Len = seq1->trueLen;
	else
		{
		part    = lookup_partition (seq1, a->beg1-1);
		seq1Len = part->trueLen;
		}

	if (sp2->p == NULL)
		seq2Len = seq2->trueLen;
	else
		{
		part    = lookup_partition (seq2, a->beg2-1);
		seq2Len = part->trueLen;
		}

	// use shorter sequence as denominator

	if (seq1Len < seq2Len) { *numer = a->end1+1 - a->beg1;  *denom = seq1Len; }
	                  else { *numer = a->end2+1 - a->beg2;  *denom = seq2Len; }
	}

//----------
//
// filter_segments_by_coverage--
//	Filter a table of segments, removing any segment that has coverage
//	percentage outside of a specified range.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	segtable*	st:				The segment table to operate upon.
//	float		minCoverage,	The range of query coverage that we will *keep*.
//				maxCoverage		.. These are values between 0 and 1.
//
// Returns:
//	(nothing)
//
//----------

void filter_segments_by_coverage
   (seq*		seq1,
	seq*		seq2,
	segtable*	st,
	float		minCoverage,
	float		maxCoverage)
	{
	segment*	srcSeg, *dstSeg;
	unspos		numer, denom;
	float		minCov, maxCov;

	if (st      == NULL) return;
	if (st->seg == NULL) return;

	// process each segment, moving those that are long enough to the first
	// part of the list

	for (dstSeg=srcSeg=st->seg ; ((u32)(srcSeg-st->seg))<st->len ; srcSeg++)
		{
		// if the segment isn't long enough, skip it

		segment_coverage (seq1, seq2, srcSeg, &numer, &denom);

		minCov = denom * minCoverage;
		maxCov = denom * maxCoverage;

		if ((numer < minCov) || (numer > maxCov))
			continue; // (unwanted segment, skip it)

		// copy the segment

		if (srcSeg != dstSeg) *dstSeg = *srcSeg;
		dstSeg++;
		}

	// set the new length of the list

	st->len = dstSeg - st->seg;
	}

//----------
//
// filter_segment_by_coverage--
//	Filter a segment, reporting whether that segment has coverage percentage
//	outside of a specified range.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	segment*	seg:			The segment to consider.  Only pos1, pos2, and
//								.. length are required.
//	float		minCoverage,	The range of query coverage that we will *keep*.
//				maxCoverage		.. These are values between 0 and 1.
//
// Returns:
//	true if the segment if outside the specified range (i.e. that it fails to
//	pass the filter, and should be discarded)
//
//----------

int filter_segment_by_coverage
   (seq*		seq1,
	seq*		seq2,
	segment*	seg,
	float		minCoverage,
	float		maxCoverage)
	{
	unspos		numer, denom;
	float		minCov, maxCov;

	segment_coverage (seq1, seq2, seg, &numer, &denom);
	minCov = denom * minCoverage;
	maxCov = denom * maxCoverage;

	if ((numer < minCov) || (numer > maxCov))
		return true; // (unwanted segment, skip it)

	return false;
	}

//----------
//
// segment_coverage--
//	Measure a segment's coverage percentage.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	segment*	seg:			The segment to measure.
//	unspos*		numer, denom:	Place to return the coverage
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//	(1)	The denominator is the length of whichever sequence is shorter, and the
//		numerator is the length of the segment.
//
//----------

void segment_coverage
   (seq*			seq1,
	seq*			seq2,
	segment*		seg,
	unspos*			numer,
	unspos*			denom)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			seq1Len, seq2Len;

	if (!seq1->needTrueLen)
		suicidef ("internal error, in segment_coverage, target does not have true lengths");

	if (!seq2->needTrueLen)
		suicidef ("internal error, in segment_coverage, query does not have true lengths");

	// determine sequence lengths for this segment

	if (sp1->p == NULL)
		seq1Len = seq1->trueLen;
	else
		{
		part    = lookup_partition (seq1, seg->pos1);
		seq1Len = part->trueLen;
		}

	if (sp2->p == NULL)
		seq2Len = seq2->trueLen;
	else
		{
		part    = lookup_partition (seq2, seg->pos2);
		seq2Len = part->trueLen;
		}

	// use shorter sequence as denominator

	*numer = seg->length;
	*denom = (seq1Len < seq2Len)? seq1Len
	                            : seq2Len;
	}

