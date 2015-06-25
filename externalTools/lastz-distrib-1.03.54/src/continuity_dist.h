//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: continuity_dist.h
//
//----------

#ifndef continuity_dist_H		// (prevent multiple inclusion)
#define continuity_dist_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef continuity_dist_owner
#define global
#else
#define global extern
#endif

//----------
//
// prototypes for routines in continuity_dist.c
//
//----------

alignel* filter_aligns_by_continuity (alignel* alignList, 
		                              float minContinuity, float maxContinuity);
alignel* filter_aligns_by_num_gaps   (alignel* alignList, s32 maxSeparateGapsCount);
alignel* filter_aligns_by_num_gap_columns (alignel* alignList, s32 maxGapColumnsCount);
void     alignment_continuity        (alignel* a,
                                      unspos* numer, unspos* denom);
void     alignment_gap_rate          (alignel* a,
                                      unspos* numer, unspos* denom);

#undef global
#endif // continuity_dist_H
