//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: tweener.h
//
//----------

#ifndef tweener_H				// (prevent multiple inclusion)
#define tweener_H

// other files

#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed matching stuff
#include "chain.h"				// segment chaining stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef tweener_owner
#define global
#else
#define global extern
#endif

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int   numTweeners;
	u64   totalArea;
	int   totalCoverage;
	int   tweenerCoverage;
	} tweenerStats;

// stats macros

#define tweener_count_stat(field)   ++tweenerStats.field
#define tweener_uncount_stat(field) --tweenerStats.field
#define tweener_set_stat(field,val) (tweenerStats.field = val)
#define tweener_add_stat(field,val) (tweenerStats.field += val)
#else
#define tweener_count_stat(field)
#define tweener_uncount_stat(field)
#define tweener_set_stat(field,val)
#define tweener_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void tweener_zero_stats    (void);
void tweener_show_stats    (FILE* f);
void tweener_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in tweener.c
//
//----------

alignel* tweener_interpolate (alignel* a, seq* seq1, seq* seq2,
                              int selfCompare, int inhibitTrivial,
                              const s8 charToBits[], seed* tweenSeed,
                              scoreset* scoring, scoreset* maskedScoring,
                              tback* tb, score xDrop, int gappedAllBounds,
                              score yDrop, int trimToPeak, score scoreThresh,
                              score diagPen, score antiPen,
                              int scale, chainer connect, u32 windowSize);

#undef global
#endif // tweener_H
