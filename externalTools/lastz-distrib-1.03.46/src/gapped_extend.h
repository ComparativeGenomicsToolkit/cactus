//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: gapped_extend.h
//
//----------

#ifndef gapped_extend_H			// (prevent multiple inclusion)
#define gapped_extend_H

// other files

#include "utilities.h"			// utility stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef gapped_extend_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef gapped_extend_owner
int gapped_extend_verbosity = 0;	// ranges from 0 (no info) to 10 (everything)
int gapped_extend_inhibitTruncationReport = false;
									// true => don't report alignment truncations
int gapped_extend_dbgShowIdentity = false;
int gapped_extend_dbgShowHsps     = false;
int gapped_extend_dbgShowAnchors  = false;
int gapped_extend_dbgAllowBatches = false;
#ifdef tryout
int gapped_extend_dbgTriviality   = false;
#endif // tryout
#else
global int gapped_extend_verbosity;
global int gapped_extend_inhibitTruncationReport;
global int gapped_extend_dbgShowIdentity;
global int gapped_extend_dbgShowHsps;
global int gapped_extend_dbgShowAnchors;
global int gapped_extend_dbgAllowBatches;
#ifdef tryout
global int gapped_extend_dbgTriviality;
#endif // tryout
#endif

//----------
//
// data structures and types
//
//----------

// traceback data structure

typedef struct tback
	{
	u32		size;				// the number of entries allocated for cell[]
	u8		space[1];			// the traceback cells
	} tback;

// special data structure for gappily_extend_hsps(info,...)

typedef struct hitrepgappily
	{
	seq*		seq1;
	seq*		seq2;
	u8*			rev1;
	u8*			rev2;
	scoreset*	scoring;
	score		yDrop;
	int			trimToPeak;
	sthresh		scoreThresh;
	tback*		traceback;
	float		minIdentity;
	float		maxIdentity;
	float		minCoverage;
	float		maxCoverage;
	float		minContinuity;
	float		maxContinuity;
	u32			minMatchCount;
	s32			maxMismatchCount;
	s32			maxSeparateGapsCount;
	s32			maxGapColumnsCount;
	int			deGapifyOutput;
	} hitrepgappily;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int64 numAnchors;
	int64 numAnchorsExtended;
	int64 numPeaks;
#if (scoreType == 'I')
	s64   totalPeakScore;
#else
	score totalPeakScore;
#endif // scoreType
	int64 numExtensions;
	int64 dpCellsVisited;
	int64 maxDpRows;
	int64 maxDpColumns;
	int64 zallocCallsA;
	int64 zallocTotalA;
	int64 zallocCallsB;
	int64 zallocTotalB;
	} gappedExtendStats;

// stats macros

#define gapped_extend_count_stat(field)   ++gappedExtendStats.field
#define gapped_extend_uncount_stat(field) --gappedExtendStats.field
#define gapped_extend_set_stat(field,val) (gappedExtendStats.field = val)
#define gapped_extend_add_stat(field,val) (gappedExtendStats.field += val)
#define gapped_extend_max_stat(field,val) if (val > gappedExtendStats.field) gappedExtendStats.field = val
#else
#define gapped_extend_count_stat(field)
#define gapped_extend_uncount_stat(field)
#define gapped_extend_set_stat(field,val)
#define gapped_extend_add_stat(field,val)
#define gapped_extend_max_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void gapped_extend_zero_stats    (void);
void gapped_extend_show_stats    (FILE* f);
void gapped_extend_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in gapped_extend.c
//
//----------

void     reduce_to_points     (seq* seq1, seq* seq2, scoreset* scoring,
                               segtable* anchors);
alignel* gapped_extend        (seq* seq1, u8* rev1, seq* seq2, u8* rev2,
                               int inhibitTrivial,
                               scoreset* scoring, segtable* anchors, tback* tb,
                               int allBounds, score yDrop, int trimToPeak,
                               sthresh scoreThresh,
                               u64 maxPairedBases,
                               int overlyPairedWarn, int overlyPairedKeep);
void     free_segment_batches (void);
tback*   new_traceback        (u32 size);
void     free_traceback       (tback* tb);
void     free_traceback_rows  (void);

u32      gappily_extend_hsps  (void* info,
                               unspos pos1, unspos pos2, unspos length,
                               score s);
score    score_alignment      (scoreset* scoring,
                               seq* seq1, unspos pos1,
                               seq* seq2, unspos pos2,
                               editscript* script);

#ifdef dbgTimingGappedExtend
void gapped_extend_timing_report (FILE* f);
#endif // dbgTimingGappedExtend

#undef global
#endif // gapped_extend_H
