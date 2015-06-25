//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: segment.h
//
//----------

#ifndef segment_H				// (prevent multiple inclusion)
#define segment_H

//----------
//
// other files
//
//----------

#include <stdio.h>				// standard C i/o stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff

// establish ownership of global variables

#ifdef segment_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef segment_owner
int segment_dbgAnchorParsing = false;	// true => debug anchor parsing
#else
global int segment_dbgAnchorParsing;
#endif

//----------
//
// data structures and types
//
//----------

// segment--
//	A segment is a pair of intervals with the same length.  This can represent
//	many things, for example an ungapped alignment between two sequences.

typedef struct segment
	{
	unspos	pos1;				// start of the interval, in one sequence (this
								// .. is origin-zero)
	unspos	pos2;				// start of the interval, in the other sequence
								// .. (this is origin-zero)
	unspos	length;				// length of the interval (e.g. number of
								// .. nucleotides
	score	s;					// segment score, (e.g. score of an ungapped
								// .. match)
	int		id;					// identifier to be used at caller's will
	possum	scoreCov;			// total lengths of the subheap (rooted at this
								// .. node) that has the same score as this
								// .. segment;  this is only valid when the
								// .. segtable is a proper min-heap
	int		filter;				// true => this segment should be discarded
	} segment;

// segment table--
//	A segment table is a list of segments.

typedef struct segtable
	{
	u32		size;				// the number of entries allocated for seq[]
	u32		len;				// the number of entries in seq[] that are
								// .. actually used
	int		haveScores;			// true  => the segments have been scored
								// false => they have not been scored
	unspos	coverageLimit;		// 'suggested' limit on the total lengths of
								// .. the segments in the table (see discussion
								// .. in file header of segment.c about how the
								// .. limit is honored);  zero indicates no
								// .. limit
	possum	coverage;			// total lengths of the segments in the table
	score	lowScore;			// score of lowest segment in the table;  if
								// .. there are no segments in the table this
								// .. is worstPossibleScore
	segment	seg[1];				// the segment table (variable-length array)
	} segtable;

//----------
//
// prototypes for routines in segment.c
//
//----------

segtable* new_segment_table          (u32 size, unspos coverageLimit);
void      empty_segment_table        (segtable* st);
void      limit_segment_table        (segtable* st, unspos coverageLimit);
void      free_segment_table         (segtable* st);
segtable* read_segment_table         (FILE* f, char* fName, segtable* st,
                                      seq* target, seq* query);
segtable* add_segment                (segtable* st,
                                      unspos pos1, unspos pos2, unspos length,
                                      score s, int id);
void      split_segment_table        (segtable* st, int id, segtable** leftovers);
void      score_segments             (segtable* st, seq* seq1, seq* seq2,
                                      scoreset* scoring);
void      sort_segments              (segtable* st,
                                      int (*qCompare) (const void* el1, const void* el2));
void      sort_some_segments         (segtable* st,
                                      u32 start, u32 end,
                                      int (*qCompare) (const void* el1, const void* el2));
void      merge_segments             (segtable* st);
void      filter_marked_segments     (segtable* st);
int       qSegmentsByPos1            (const void* segA, const void* segB);
int       qSegmentsByPos2            (const void* segA, const void* segB);
int       qSegmentsByDecreasingScore (const void* segA, const void* segB);
int       qSegmentsByIncreasingScore (const void* segA, const void* segB);
int       qSegmentsByDiag            (const void* segA, const void* segB);
int       qSegmentsById              (const void* segA, const void* segB);
void      write_segments             (FILE* f, segtable* st,
                                      seq* target, seq* query,
                                      int withText);
void      dump_segments              (FILE* f, segtable* st,
                                      char* sym1, char* sym2);

#undef global
#endif // segment_H
