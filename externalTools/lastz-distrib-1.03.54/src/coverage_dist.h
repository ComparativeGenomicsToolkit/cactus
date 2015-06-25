//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: coverage_dist.h
//
//----------

#ifndef coverage_dist_H			// (prevent multiple inclusion)
#define coverage_dist_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef coverage_dist_owner
#define global
#else
#define global extern
#endif

//----------
//
// prototypes for routines in coverage_dist.c
//
//----------

alignel* filter_aligns_by_coverage   (seq* seq1, seq* seq2, alignel* alignList, 
		                              float minCoverage, float maxCoverage);
void     alignment_coverage          (seq* seq1, seq* seq2, alignel* a,
                                      unspos* numer, unspos* denom);
void     filter_segments_by_coverage (seq* seq1, seq* seq2, segtable* st, 
		                              float minCoverage, float maxCoverage);
int      filter_segment_by_coverage  (seq* seq1, seq* seq2, segment* seg,
		                              float minCoverage, float maxCoverage);
void     segment_coverage            (seq* seq1, seq* seq2, segment* seg,
                                      unspos* numer, unspos* denom);

#undef global
#endif // coverage_dist_H
