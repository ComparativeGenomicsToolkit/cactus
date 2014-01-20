//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: infer_scores.h
//
//----------

#ifndef infer_scores_H			// (prevent multiple inclusion)
#define infer_scores_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "sequences.h"			// sequence stuff
#include "pos_table.h"			// position table stuff
#include "seed_search.h"		// seed hit search stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef infer_scores_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef infer_scores_owner
int infer_scores_watchConverge   = false;	// true => report scores so we can watch convergence
int infer_scores_snoopConverge   = false;	// true => report stats so we can watch convergence
int infer_scores_showParams      = false;	// true => report the inference parameters  
int infer_scores_outputLav       = false;	// true => output LAV for inference alignments
int infer_scores_dbgShowIdentity = false;
#else
global int infer_scores_watchConverge;
global int infer_scores_snoopConverge;
global int infer_scores_showParams;
global int infer_scores_outputLav;
global int infer_scores_dbgShowIdentity;
#endif

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int64  coverageTotal;
	int64  coverageLow;
	int64  coverageHigh;
	int    lowIdentityBin;
	int    highIdentityBin;
	u32    subs[4][4];
	u32    n1[4];
	u32    n2[4];
	u32    m[4][4];
	double averageGapLength;
	double averageSegmentLength;
	double pExtend;
	double sExtend;
	double pOpen;
	double sOpen;
	double scaleBy;
	} inferScoresStats;

// stats macros

#define infer_scores_count_stat(field)   ++inferScoresStats.field
#define infer_scores_uncount_stat(field) --inferScoresStats.field
#define infer_scores_set_stat(field,val) (inferScoresStats.field = val)
#define infer_scores_add_stat(field,val) (inferScoresStats.field += val)
#define infer_scores_max_stat(field,val) if (val > inferScoresStats.field) inferScoresStats.field = val
#else
#define infer_scores_count_stat(field)
#define infer_scores_uncount_stat(field)
#define infer_scores_set_stat(field,val)
#define infer_scores_add_stat(field,val)
#define infer_scores_max_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void infer_scores_zero_stats        (void);
void infer_scores_show_stats_common (FILE* f);
void infer_scores_show_stats_subs   (FILE* f, int trial);
void infer_scores_show_stats_gaps   (FILE* f, int trial);
void infer_scores_show_stats        (FILE* f);
void infer_scores_generic_stats     (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in infer_scores.c
//
//----------

scoreset* drive_scoring_inference (void* params,
                                   seq* target, u8* targetRev,
                                   postable* targPositions,
                                   seq* query,
                                   tback* traceback);

void gather_stats_from_align_list (alignel* alignList, seq* seq1, seq* seq2);
void gather_stats_from_match      (seq* seq1, unspos pos1,
                                   seq* seq2, unspos pos2, unspos length);

void init_inference_stats_job     (seq* seq1, seq* seq2);
void print_inference_stats_job    (FILE* f);
void infer_stats_from_align_list  (alignel* alignList, seq* seq1, seq* seq2);
void infer_stats_from_match       (seq* seq1, unspos pos1,
                                   seq* seq2, unspos pos2, unspos length);

#undef global
#endif // infer_scores_H
