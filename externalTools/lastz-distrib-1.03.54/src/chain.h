//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: chain.h
//
//----------

#ifndef chain_H					// (prevent multiple inclusion)
#define chain_H

// other files

#include "utilities.h"			// utility stuff
#include "segment.h"			// segment table management stuff

// establish ownership of global variables

#ifdef chain_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef chain_owner
int chain_dbgChaining = false;
int chain_dbgDumpTree = false;
#else
global int chain_dbgChaining;
global int chain_dbgDumpTree;
#endif

//----------
//
// data structures and types
//
//----------

// chain connection penalty functions--
//	Compute a chaining penalty between two segments.
//
//	Arguments:
//		segment*	seq1:	The first segment.
//		segment*	seq2:	The other segment.
//		int			scale:	Scaling factor (see note).
//
// Note:
//		The effect of the scale parameter is to permit integer arithmetic to be
//		used with very small gap penalties.  This is only useful if the scoring
//		type is an integer.

typedef score (*chainer) (segment* , segment* , int);

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int numAnchors;
	int numSegments;
	} chainStats;

// stats macros

#define chain_count_stat(field)   ++chainStats.field
#define chain_uncount_stat(field) --chainStats.field
#define chain_set_stat(field,val) (chainStats.field = val)
#define chain_add_stat(field,val) (chainStats.field += val)
#else
#define chain_count_stat(field)
#define chain_uncount_stat(field)
#define chain_set_stat(field,val)
#define chain_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void chain_zero_stats    (void);
void chain_show_stats    (FILE* f);
void chain_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in chain.c
//
//----------

score reduce_to_chain (segtable* st, score diagPen, score antiPen,
                       int scale, chainer connect);

#undef global
#endif // chain_H
