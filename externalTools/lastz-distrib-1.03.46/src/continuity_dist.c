//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: continuity_dist.c
//
//----------
//
// continuity_dist--
//	Support for collecting the query continuity (or gap rate) distribution from
//	alignments.
//
// Notes:
//	continuity = 1/(1+gaprate)
//  gaprate    = (1/continuity)-1
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

#define  continuity_dist_owner	// (make this the owner of its globals)
#include "continuity_dist.h"	// interface to this module

//----------
//
// filter_aligns_by_continuity--
//	Filter a list of alignments, removing any alignment that has continuity
//	percentage outside of a specified range.
//
//----------
//
// Arguments:
//	alignel*	alignList:		The list of alignments to operate upon.
//	float		minContinuity,	The range of query continuity in alignments that
//				maxContinuity:	.. we will *keep*.  These are values between 0
//								.. and 1.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
//
// Notes:
//	(1)	The numerator is the number of alignment columns *not* containing gaps
//      and the denominator is the number of total alignment columns.
//	(2)	Memory for alignments that don't make the cut is deallocated here.
//	(3)	The returned list of alignments is in the same order as the incoming
//		list.
//
//----------

// $$$ There's an inherent inefficiency here if the user asks us to filter
// $$$ .. alignments by several gap-based stats.  We'd like to compute the
// $$$ .. gap stats only once, then either carry that in the alignment record,
// $$$ .. or apply all the filters in one pass through the loop.

alignel* filter_aligns_by_continuity
   (alignel*	alignList,
	float		minContinuity,
	float		maxContinuity)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;
	float		minCont, maxCont;

	// process each alignment, collecting a list of those that are long enough

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		// if the alignment is too gappy, skip it

		alignment_continuity (a, &numer, &denom);

		minCont = denom * minContinuity;
		maxCont = denom * maxContinuity;

		if ((numer < minCont) || (numer > maxCont))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_continuity a->script", a->script);
			free_if_valid ("filter_aligns_by_continuity a",         a);
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
// filter_aligns_by_num_gaps--
//	Filter a list of alignments, removing any alignment that has more gaps than
//	a specified maximum (counting each run of gapped columns as a single gap).
//
//----------
//
// Arguments:
//	alignel*	alignList:				The list of alignments to operate upon.
//	s32			maxSeparateGapsCount:	The maximum number of gaps in alignments
//										.. that we will *keep*.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
//
// Notes:
//	(1)	Memory for alignments that don't make the cut is deallocated here.
//	(2)	The returned list of alignments is in the same order as the incoming
//		list.
//
//----------

alignel* filter_aligns_by_num_gaps
   (alignel*	alignList,
	s32			maxSeparateGapsCount)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		height, width, i, j, run;
	u32			opIx;
	s32			numGaps;

	// process each alignment, collecting a list of those that are long enough

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		// if the alignment is too gappy, skip it

		numGaps = 0;

		height = a->end1 - a->beg1 + 1;
		width  = a->end2 - a->beg2 + 1;
		opIx = 0;
		for (i=j=0 ; (i<height)||(j<width) ; )
			{
			// handle the next run

			run = edit_script_run_of_subs (a->script, &opIx);
			i += run; j += run;

			// handle the next indel

			if ((i < height) || (j < width))
				{
				edit_script_indel_len (a->script, &opIx, &i, &j);
				if (++numGaps > maxSeparateGapsCount)
					break;
				}
			}

		if (numGaps > maxSeparateGapsCount)
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_continuity a->script", a->script);
			free_if_valid ("filter_aligns_by_continuity a",         a);
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
// filter_aligns_by_num_gap_columns--
//	Filter a list of alignments, removing any alignment that has more gaps than
//	a specified maximum (counting each gapped column as a single gap).
//
//----------
//
// Arguments:
//	alignel*	alignList:			The list of alignments to operate upon.
//	s32			maxGapColumnsCount:	The maximum number of gaps in alignments
//									.. that we will *keep*.  Note that this
//									.. must not be negative.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
//
// Notes:
//	(1)	Memory for alignments that don't make the cut is deallocated here.
//	(2)	The returned list of alignments is in the same order as the incoming
//		list.
//
//----------

alignel* filter_aligns_by_num_gap_columns
   (alignel*	alignList,
	s32			maxGapColumnsCount)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;

	// process each alignment, collecting a list of those that are long enough

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		// if the alignment is too gappy, skip it

		alignment_continuity (a, &numer, &denom);

		if ((denom == 0) || (denom-numer > (u32) maxGapColumnsCount))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_continuity a->script", a->script);
			free_if_valid ("filter_aligns_by_continuity a",         a);
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
// alignment_continuity--
//	Compute the continuity of an gapped alignment block.  This is the number of
//	bases *not* aligned to gaps divided by the number of alignment columns.
//
//----------
//
// Arguments:
//	alignel*	a:				The alignment of interest.
//	unspos*		numer, denom:	Place to return the continuity fraction.  Note
//								.. that the returned denominator might be zero,
//								.. and that the fraction may be greater than 1.
//
// Returns:
//	(nothing)
//
//----------

void alignment_continuity
   (alignel*	a,
	unspos*		_numer,
	unspos*		_denom)
	{
	unspos		gapColumns, nonGapColumns;

	alignment_gap_rate (a, &gapColumns, &nonGapColumns);

	*_numer = nonGapColumns;
	*_denom = nonGapColumns + gapColumns;
	}

//----------
//
// alignment_gap_rate--
//	Compute the gap rate of an gapped alignment block.  This is the number of
//	bases aligned to gaps divided by the number of aligned bases.
//
//----------
//
// Arguments:
//	alignel*	a:				The alignment of interest.
//	unspos*		numer, denom:	Place to return the gap rate fraction.  Note
//								.. that the returned denominator might be zero,
//								.. and that the fraction may be greater than 1.
//
// Returns:
//	(nothing)
//
//----------

void alignment_gap_rate
   (alignel*	a,
	unspos*		_numer,
	unspos*		_denom)
	{
	unspos		beg1 = a->beg1;
	unspos		beg2 = a->beg2;
	unspos		height, width, i, j;
	u32			opIx;
    unspos		run;
	unspos		denom, gappedBases;

	height = a->end1 - beg1 + 1;
	width  = a->end2 - beg2 + 1;

	denom = 0;
	opIx  = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		run = edit_script_run_of_subs (a->script, &opIx);
		i += run; j += run;

		denom += run;

		if ((i < height) || (j < width))
			edit_script_indel_len (a->script, &opIx, &i, &j);
		}

	if (denom == 0)
		{ *_numer = *_denom = 0;  return; }

	gappedBases = (height - denom) + (width - denom);

	*_numer = gappedBases;
	*_denom = denom;
	}

