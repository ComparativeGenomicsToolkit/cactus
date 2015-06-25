//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: seed_search.h
//
//----------

#ifndef seed_search_H			// (prevent multiple inclusion)
#define seed_search_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "pos_table.h"			// position table stuff
#include "segment.h"			// segment table management stuff

// establish ownership of global variables

#ifdef seed_search_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access
// nota bene: showProgress is a relic which is not currently used

#ifdef seed_search_owner
int seed_search_showProgress   = false; // true => make periodic progress reports
int seed_search_dbgDumpRawHits = false; // true => dump hits in find_table_matches()
int seed_search_dbgShowRawHits = false; // true => show hits in find_table_matches()
int seed_search_dbgShowHits    = false; // true => show hits in process_for_plain_hit(),
										//         .. process_for_simple_hit(),
										//         .. process_for_recoverable_hit(),
										//         .. or process_for_twin_hit()
int seed_search_dbgShowCoverage= false; // true => show bases hit count in
										//         .. seed_hit_search(),
int seed_search_dbgShowRejections= false;// true => report queries rejected in
										//         .. seed_hit_search(),


int seed_search_dbgSearchLimitExceeded = 0;
#ifdef snoopHspSubrange
int seed_search_dbgSubrangeHsps = false;// true => show all suboptimal HSPs
#endif // snoopHspSubrange
#else
global int seed_search_showProgress;
global int seed_search_dbgDumpRawHits;
global int seed_search_dbgShowRawHits;
global int seed_search_dbgShowHits;
global int seed_search_dbgShowCoverage;
global int seed_search_dbgShowRejections;
global int seed_search_dbgSearchLimitExceeded;
#ifdef snoopHspSubrange
global int seed_search_dbgSubrangeHsps;
#endif // snoopHspSubrange
#endif

//----------
//
// seed hit reporter functions--
//	Report the discovery of a seed hit or HSP.
//
//	Arguments:
//		void*	info:	Additional control/arguments specific to the seed hit
//						.. reporter being called.
//		unspos	pos1:	position, in the first sequence, of first character
//						.. *after* the last character in the match (origin-0)
//		unspos	pos2:	position, in the second sequence, of first character
//						.. *after* the last character in the match (origin-0)
//		unspos	length:	number of nucleotides
//		score	s:		the match's "score"
//
//	Returns:
//		The number of bases in the seed hit (or HSP);  0 if the hit is
//		rejected.
//
//----------

typedef u32 (*hitreporter) (void*, unspos, unspos, unspos, score);

//----------
//
// seed hit processor functions--
//	Process a simple seed hit for a given seed pattern.
//
//----------
//
// Arguments:
//		void*	info:	Additional control/arguments specific to the seed hit
//						.. processor being called.  This will usually be of a
//						.. type that includes at least the stuff in hitprocinfo.
//		unspos	pos1:	The hit position in sequence 1, relative to the
//						.. entire sequence (not to the interval).  This is
//						.. the first letter following the end of the match
//						.. (origin-0).
//		unspos	pos2:	The hit position in sequence 2 (with details the
//						.. same as for pos1).
//		unspos	length:	The length of the hit (number of nucleotides).
//
//	Returns:
//		The number of bases in the seed hit (or HSP);  0 if the hit is
//		rejected.
//
//----------

typedef u64 (*hitprocessor) (void*, unspos, unspos, unspos);

// basic data structure for all hit processors

typedef struct hitprocinfo
	{
	hitreporter	reporter;		// function to call to report each hit that is
								// .. 'good enough'
	void*		reporterInfo;	// value to pass thru with each call to reporter

	// filtering

	int			posFilter;		// true => discard seed hits outside of
	interval	targetInterval;	//         .. tStart,tEnd in target->v[] or
	interval	queryInterval;	//         .. qStart,qEnd in query->v[];  the
								//         .. intervals are origin-zero closed

	int			minMatches;		// filter criteria for each seed hit;  we
	int			maxTransversions;//.. require at least minMatches matches and no
	char*		filterPattern;	// .. more than maxTransversions transversions;
								// .. if minMatches<0 no filtering is performed;
								// .. if maxTransversions<0 there is no
								// .. transversion limit;  if filterPattern is
								// .. non-NULL the filter criteria only applies
								// .. to the pattern's "care" positions (and
								// .. not to any don't-care positions);
								// .. filterPattern points directly into the
								// .. active seed, and needn't be deallocated
	const s8*	charToBits;		// table to map sequence characters to two-bit
								// .. values, and illegal characters to -1;
								// .. indexed by a u8 value, 0..255;  normally
								// .. this will consider upper and lower case to
								// .. be the same

	// gap-free extension

	int			gfExtend;		// whether to extend seed hits into HSPs (one
								// .. of gfexXXX)
	seq*		seq1;
	seq*		seq2;
	scoreset*	scoring;		// the scoring scheme;  usually this treats
								// .. lowercase letters as being 'bad'
	score		xDrop;
	sthresh		hspThreshold;
	score		hspZeroThreshold; // max(0,hspThreshold.s)
	segtable**	anchors;
	int			entropicHsp;
	int			reportEntropy;
	} hitprocinfo;

// legal values for gfExtend;  note that gfexMismatch_min has to be 1, since
// values in the range gfexMismatch_min..gfexMismatch_max are identical to the
// number of mismatches they represent

enum
	{
	gfexNoExtend = -2,			// simple gfex, no hash collision detection
	gfexXDrop = -1,				// extend seed hits into HSPs, using xdrop
	gfexExact = 0,				// extend seed hits using exact match
	gfexMismatch_min = 1,		// 1..50 => extend seed hits using N-mismatch
	gfexMismatch_max = 50
	};

// special data structure for process_for_plain_hit(info,...) and
// process_for_simple_hit(info,...)

typedef struct hitprocsimple
	{
	hitprocinfo	hp;				// basic info
	} hitprocsimple;

// special data structure for process_for_twin_hit(info,...)

typedef struct hitproctwin
	{
	hitprocinfo	hp;				// basic info
	u32			minSpan;		// span threshold for hits to be considered
	u32			maxSpan;		// .. twins;  we require two (or more) hits
								// .. with end2-start1 between minSpan and
								// .. maxSpan, inclusive
	} hitproctwin;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

//#define maxHitsPerColumn 1000

global struct
	{
	int   withTrans;
	int   minMatches;
	int   maxTransversions;
	int   filterCaresOnly;
	int   isHspSearch;
	u32   searchLimit;

	int   wordsInSequence;
	int64 unresolvedSeedHits;
	int64 rawSeedHits;
	int64 hashCollisions;
	int64 hashFailures;
	int64 notEnoughMatches;
	int64 tooManyTransversions;
	int64 bpExtended;
	int64 lowScoringHsps;
	int64 hsps;

#ifdef snoopHspSubrange
	int64 suboptimalHsp;
	int64 suboptimalHspB;
#endif // snoopHspSubrange

#ifndef noSeedHitQueue
	int64 queueSeedsScanned;
	int64 queueSeedsExamined;
	int64 queueSeedsBlocked;
#endif // not noSeedHitQueue

#ifdef maxHitsPerColumn
	int64 hitsPerColumn[maxHitsPerColumn+2];
	u64   mostHitsInColumn;
#endif // maxHitsPerColumn
	} seedSearchStats;

// stats macros

#define seed_search_count_stat(field)   ++seedSearchStats.field
#define seed_search_uncount_stat(field) --seedSearchStats.field
#define seed_search_set_stat(field,val) (seedSearchStats.field = val)
#define seed_search_add_stat(field,val) (seedSearchStats.field += val)
#else
#define seed_search_count_stat(field)
#define seed_search_uncount_stat(field)
#define seed_search_set_stat(field,val)
#define seed_search_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void  seed_search_zero_stats       (void);
void  seed_search_show_stats       (FILE* f);
void  seed_search_generic_stats    (FILE* f, void (*func) (FILE*, const char*, ...));
int64 seed_search_hsps             (void);
int64 seed_search_low_scoring_hsps (void);
int64 seed_search_bp_extended      (void);

//----------
//
// prototypes for routines in seed_search.c
//
//----------

u64    seed_hit_search        (seq* seq1, postable* pt,
                               seq* seq2, unspos start, unspos end,
                               int selfCompare,
                               const s8 charToBits[], seed* hitSeed,
                               u32 searchLimit, u32 reportSearchLimit,
#ifdef densityFiltering
                               double maxDensity,
#endif // densityFiltering
                               hitprocessor processor, void* processorInfo);
void   free_seed_hit_search   (void);
u64    process_for_plain_hit  (void* info,
                               unspos pos1, unspos pos2, unspos length);
u64    process_for_simple_hit (void* info,
                               unspos pos1, unspos pos2, unspos length);
u64    process_for_recoverable_hit (void* info,
                               unspos pos1, unspos pos2, unspos length);
u64    process_for_twin_hit   (void* info,
                               unspos pos1, unspos pos2, unspos length);
float  discovery_probability  (seq* seq1, unspos pos1,
                               seq* seq2, unspos pos2, unspos length,
                               seed* hitSeed, u32 step);

#undef global
#endif // seed_search_H
