//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: infer_scores.c
//
//----------
//
// infer_scores--
//	Support for collecting alignment scoring inference stats from allignments.
//	"Alignment scoring stats" are background nucleotide counts, substitituion
//	counts, and gap length distribution.
//
// References:
//
//	[1]	"Scoring Pairwise Genomic Sequence Alignments"  F Chiaromonte, VB
//		Yap, W Miller.  Pacific Symposium on Biocomputing (2002), vol. 7, pp.
//		115-126
//
//	[2]	"Biological sequence analysis"  Durbin, Eddy, Krogh and Mitchison.
//		Cambridge University Press, 1998.  Pages 29-31.
//
//	[3] "Improved Pairwise Alignment of Genomic DNA". RS Harris, PhD Thesis,
//		Pennsylvania State University, 2007.
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
#include <math.h>				// standard C math stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "identity_dist.h"		// identity distribution "format" stuff
#include "output.h"				// alignment outout format stuff
#include "lastz.h"				// lastz program-wide stuff

#define  infer_scores_owner		// (make this the owner of its globals)
#include "infer_scores.h"		// interface to this module

//----------
//
// private global data
//
//----------

#define maxSubIterations 30		// max number of iterations for inferring
								// .. substitution scores
#define maxGapIterations 30		// max number of iterations for inferring gap
								// .. scores

#if (scoreType == 'I')
#define subCloseEnough 0
#define gapCloseEnough 0
#else
#define subCloseEnough .000001
#define gapCloseEnough .0001
#endif

// distribution--
//	A set of (length,count) pairs.
//
//	All memory is self-contained, so a direct call to free() can be used for
//	disposal.
//
//	$$$ Eventually we should improve this by using something like a balanced
//	    binary tree or a hash-table.  For now we just use the simple but slow
//	    method.
//
//	$$$ Moreover, we are collecting more stats than we need.  Currently we infer
//	    scores only from refGaps, secGaps, and segments.  And the two gap
//	    distributions can be combined.
//
//	$$$ Since we only need averages from the distributions, we could just store
//	    count and sum.  However, we will eventually try to make use of the
//	    actual shape of the distribution.  In particular, we need to make a
//	    correction for the inherent lack of short gaps.

typedef struct dpair
	{
	unspos	length;
	u64		count;
	} dpair;

typedef struct distn
	{
	u32		size;				// the number of entries allocated for items[]
	u32		len;				// the number of items used
	dpair*	items;				// the (length,count) pairs
	} distn;

// inference stats

typedef struct infstats
	{
	u64		count;
	possum	coverage;

	unspos	refBases;
	unspos	secBases;
	unspos	refBkgd[4];
	unspos	secBkgd[4];
	unspos	subs[4][4];

	distn*	refBlocks;			// alignment block length distributions
	distn*	secBlocks;
	distn*	refGaps;			// gap length distributions
	distn*	secGaps;
	distn*	refRuns;			// non-gap length distributions
	distn*	secRuns;
	distn*	segments;			// ungapped pair length distribution
	} infstats;

// short score arrays (for tracking score convergence)

typedef struct score2 { score s1,s2;             } score2;
typedef struct score6 { score s1,s2,s3,s4,s5,s6; } score6;

// private globals shared by all the routines herein

control*		params;
seq*			target;
u8*				targetRev;
postable*		targPositions;
seq*			query;
tback*			traceback;
float			minIdentity;
float			maxIdentity;
scoreset*		inferredScoring;

static int      statsActive = false;
static infstats infStats;

static infstats infStatsByPctId[numIdentityBins+1];

//----------
//
// prototypes for private functions
//
//----------

static double align_for_sub_scores        (score scaleTo);
static double align_for_gap_scores        (score scaleTo);
static void   align_for_stats             (hitprocessor hitProc,
                                           void* hitProcInfo);

static double infer_substitution_scores   (double pOpen, score scaleTo);
static double infer_gap_scores            (score scaleTo);

static void   copy_scores                 (scoreset* dst, scoreset* src);
static void   repair_scores               (scoreset* scoring, scoreset* masked);
static void   write_scores                (char* fileId, scoreset* ss,
                                           int withGapScores, int withExtras,
                                           int asInts);

static void   init_stats_for_inference    (seq* seq1, seq* seq2);
static void   erase_stats_for_inference   (void);
static void   free_stats_for_inference    (void);
static void   filter_stats_by_percentile  (void);
static void   combine_binned_stats        (int mergeSequences);

static void   init_stats                  (infstats* inf);
static void   erase_stats                 (infstats* inf);
static void   free_stats                  (infstats* inf);
static void   accumulate_stats_from_align (seq* seq1, unspos beg1, unspos end1,
                                           seq* seq2, unspos beg2, unspos end2,
                                           editscript* script,
                                           infstats* inf);
static void   accumulate_stats_from_match (seq* seq1, unspos pos1,
                                           seq* seq2, unspos pos2,
                                           unspos length,
                                           infstats* inf);

static void   print_bkgd_stats            (FILE* f, char* s, unspos bkgd[4]);
static void   print_subs_stats            (FILE* f, unspos subs[4][4]);
static void   print_blocks_stats          (FILE* f, char* s, distn* blocks);
static void   print_gaps_stats            (FILE* f, char* s, distn* gaps);
static void   print_runs_stats            (FILE* f, char* s, distn* runs);
static void   print_segments_stats        (FILE* f, distn* segments);

static distn* init_length_distribution    (u32 numEntries);
static void   erase_length_distribution   (distn* d);
static void   free_length_distribution    (distn* d);
static void   add_lengths_to_distribution (distn* src, distn** dst);
static void   add_length_to_distribution  (unspos length, distn** d);
static u64    number_of_instances         (distn* d);
static double average_length              (distn* d);
static void   print_length_distribution   (FILE* f, char* prefix, distn* d);
static int    qCompareByLength            (const void* pairA, const void* pairB);

//----------
//
// drive_scoring_inference--
//	Infer log odds alignment scores from the sequences.
//
// Scores for substitution and gaps are inferred by
//	(a) assuming some scoring set
//	(b) performing alignment
//	(c) filtering out alignment blocks that are too good or too bad
//	(d) inferring scores from the statistics of the resulting blocks
//
// Here we iterate that process until the scores converge.  In phase I we
// restrict ourselves to gap-free alignments and infer only substitution scores.
// Then in phase II we allow gapped alignments and infer gap scores (without
// changing substitution scores).  The reason for performing substitution score
// convergence in the absense of gaps is that gap-free alignment is much faster
// than gapped alignment.  However, there is no real justification for holding
// substitution scores constant during phase II;  this just matches the way the
// author of [3] had done some of his earlier experiments.
// 
// $$$ need to add a lot more here about convergence and the potential for it
//     .. to get into an orbit.
//
// $$$ describe percent-id filtering, percentiles, etc.
//
// $$$ we could allow substitution scores to change during phase II.
//
//----------
//
// Arguments:
//	control*		params:			Parameter set controlling the inference
//									.. searches.  This is actually declared as
//									.. void* to avoid a cyclic dependency in
//									.. the include files.
//	seq*			target:			The sequence being searched.
//	postable*		targPositions:	A table of positions of words in target.
//	u8*				targetRev:		The reverse (NOT reverse complement) of the
//									.. target sequence, as a zero-terminated
//									.. string;  this may be NULL if the caller
//									.. doesn't need/want to supply it.  It is
//									.. only needed if we will be inferring gap
//									.. scores.
//	seq*			query:			The sequence(s) to compare to the target.
//									.. Upon completion, this will be rewound
//									.. back to the starting position.
//	tback*			traceback:		Memory in which to track gapped alignment
//									.. traceback.
//
// Returns:
//	A pointer to the inferred scoring set.  The caller is responsible for
//	deallocating this.
//
//----------

static int close_enough_scores_2 (score2* u, score2* v);
static int close_enough_scores_6 (score6* u, score6* v);


scoreset* drive_scoring_inference
   (void*			_params,
	seq*			_target,
	u8*				_targetRev,
	postable*		_targPositions,
	seq*			_query,
	tback*			_traceback)
	{
	score6			subScores, pastSubScores[maxSubIterations+1];
	score2			gapScores, pastGapScores[maxGapIterations+1];
	char			scoreFileId[20];
	score			maxSubScore, minSubScore, scaleTo;
	double			oneOverMaxSubScore, minOverMaxSubScore;
	double			origHspThresholdRatio, origGappedThresholdRatio,
					origGapOpenRatio, origGapExtendRatio;
	double			hspThresholdRatio, gappedThresholdRatio,
					gapOpenRatio, gapExtendRatio;
	int				showAllScores, inOrbit;
	int				trial, oldTrial;
	scoreset*		temp;

	if (((control*) _params)->gappedThreshold.t != 'S')
		suicidef ("drive_scoring_inference can't handle score threshold %s",
		          score_thresh_to_string (&((control*) _params)->gappedThreshold));

	if ((((control*) _params)->minCoverage > 0)
	 || (((control*) _params)->maxCoverage < 1))
		suicidef ("drive_scoring_inference can't handle query coverage filtering");

#if (!defined infer_anything)
	if (((control*) _params)->ic.gapIterations > 0)
		suicide ("Gap scoring inference has not been shown to produce useful results and\n"
		         "is currently blocked.  To unblock gap scoring inference, contact the author.");
#endif

	// pass globals to the rest of this module
	// note: this makes this module non-threadsafe

	params        = (control*) _params;
	target        = _target;
	targetRev     = _targetRev;
	targPositions = _targPositions;
	query         = _query;
	traceback     = _traceback;

	if (params->ic.subIterations > maxSubIterations)
		params->ic.subIterations = maxSubIterations;

	if (params->ic.gapIterations > maxGapIterations)
		params->ic.gapIterations = maxGapIterations;

	if (params->ic.idIsPercentile)
		{
		minIdentity = params->minIdentity;  params->minIdentity = 0.0;
		maxIdentity = params->maxIdentity;  params->maxIdentity = 1.0;
		}

	origHspThresholdRatio    = params->hspThreshold.s;
	origGappedThresholdRatio = params->gappedThreshold.s;
	origGapOpenRatio         = params->scoring->gapOpen;
	origGapExtendRatio       = params->scoring->gapExtend;

	// determine limiting parameters, relative to the maximum substitution
	// score

	scaleTo = params->ic.inferScale;

	maxSubScore = max_in_score_matrix (params->scoring);
	minSubScore = min_in_score_matrix (params->scoring);
	oneOverMaxSubScore = 1.0 / maxSubScore;
	minOverMaxSubScore = (-minSubScore) / (double) maxSubScore;

	hspThresholdRatio = origHspThresholdRatio;

	switch (params->ic.hspThresholdIsRatio)
		{
		case ratioNone:        hspThresholdRatio    *= oneOverMaxSubScore;  break;
		case ratioMinSubScore: hspThresholdRatio    *= minOverMaxSubScore;  break;
		}

	// allocate a scoring set to receive inferred scores;  we will ping-pong
	// between this one and the params->scoring set;  and we will maintain the
	// params->maskedScoring set as the masked version of the inferred set

	inferredScoring = new_dna_score_set (NULL, 0, 0, 0, 0);

	// decide whether we are going to ouput every score set we try

	showAllScores = false;
	if (params->ic.inferFilename != NULL)
		showAllScores = (strstr (params->ic.inferFilename, "%s") != NULL);

	// report inference settings to the user

	if (infer_scores_showParams)
		{
		fprintf (stderr, "inference parameters\n");
		if (params->ic.inferFilename == NULL)
			fprintf (stderr, "  inferFilename      = (null)\n");
		else
			fprintf (stderr, "  inferFilename      = %s\n", params->ic.inferFilename);

		if (params->ic.idIsPercentile)
			{
			fprintf (stderr, "  min_identity       = %.2f%% (as percentile)\n", 100.0*params->minIdentity);
			fprintf (stderr, "  max_identity       = %.2f%% (as percentile)\n", 100.0*params->maxIdentity);
			}
		else
			{
			fprintf (stderr, "  min_identity       = %.2f (as identity)\n", 100.0*params->minIdentity);
			fprintf (stderr, "  max_identity       = %.2f (as identity)\n", 100.0*params->maxIdentity);
			}

		if ((params->ic.inferScale == 0) && (params->ic.writeAsInt))
			fprintf (stderr, "  inference_scale    = (no scaling) forced to int\n");
		else if ((params->ic.inferScale == 0) && (!params->ic.writeAsInt))
			fprintf (stderr, "  inference_scale    = (no scaling)\n");
		else if ((params->ic.inferScale != 0) && (params->ic.writeAsInt))
			fprintf (stderr, "  inference_scale    = " scoreFmtSimple " (forced to int)\n", params->ic.inferScale);
		else
			fprintf (stderr, "  inference_scale    = " scoreFmtSimple "\n", params->ic.inferScale);

		if (params->ic.hspThresholdIsRatio == ratioMaxSubScore)
			fprintf (stderr, "  hsp_threshold      = " scoreFmtSimple "*inference_scale\n", params->hspThreshold.s);
		else if (params->ic.hspThresholdIsRatio == ratioMinSubScore)
			fprintf (stderr, "  hsp_threshold      = " scoreFmtSimple "*worst_substitution\n", params->hspThreshold.s);
		else
			fprintf (stderr, "  hsp_threshold      = " scoreFmtSimple "\n", params->hspThreshold.s);

		if (params->ic.gappedThresholdIsRatio == ratioMaxSubScore)
			fprintf (stderr, "  gapped_threshold   = " scoreFmtSimple "*inference_scale\n", params->gappedThreshold.s);
		else if (params->ic.gappedThresholdIsRatio == ratioMinSubScore)
			fprintf (stderr, "  gapped_threshold   = " scoreFmtSimple "*worst_substitution\n", params->gappedThreshold.s);
		else
			fprintf (stderr, "  gapped_threshold   = " scoreFmtSimple "\n", params->gappedThreshold.s);

		fprintf (stderr, "  max_sub_iterations = %d\n", params->ic.subIterations);
		fprintf (stderr, "  max_gap_iterations = %d\n", params->ic.gapIterations);

		if (params->ic.gapOpenIsRatio == ratioMaxSubScore)
			fprintf (stderr, "  gap_open_penalty   = " scoreFmtSimple "*inference_scale\n", params->scoring->gapOpen);
		else if (params->ic.gapOpenIsRatio == ratioMinSubScore)
			fprintf (stderr, "  gap_open_penalty   = " scoreFmtSimple "*worst_substitution\n", params->scoring->gapOpen);
		else
			fprintf (stderr, "  gap_open_penalty   = " scoreFmtSimple "\n", params->scoring->gapOpen);

		if (params->ic.gapExtendIsRatio == ratioMaxSubScore)
			fprintf (stderr, "  gap_extend_penalty = " scoreFmtSimple "*inference_scale\n", params->scoring->gapExtend);
		else if (params->ic.gapExtendIsRatio == ratioMinSubScore)
			fprintf (stderr, "  gap_extend_penalty = " scoreFmtSimple "*worst_substitution\n", params->scoring->gapExtend);
		else
			fprintf (stderr, "  gap_extend_penalty = " scoreFmtSimple "\n", params->scoring->gapExtend);

		fprintf (stderr, "  step               = %d\n", params->step);
		fprintf (stderr, "  entropy            = %s\n", (params->entropicHsp)? "on" : "off");
		}

	//////////
	// Phase I-- iterate substitution scores inference
	//////////

	init_stats_for_inference (params->seq1, params->seq2);

	if (infer_scores_outputLav) params->outputFormat = fmtLavInfScores;
	                       else params->outputFormat = fmtInfScores;

	subScores.s1 = params->scoring->sub['A']['A'];
	subScores.s2 = params->scoring->sub['T']['T'];
	subScores.s3 = params->scoring->sub['A']['C'];
	subScores.s4 = params->scoring->sub['A']['G'];
	subScores.s5 = params->scoring->sub['A']['T'];
	subScores.s6 = params->scoring->sub['C']['G'];
	pastSubScores[0] = subScores;
	inOrbit = false;

	for (trial=1 ; (!inOrbit)&&(trial<=params->ic.subIterations) ; trial++)
		{
		// determine limiting parameters, relative to the current maximum
		// substitution score

		maxSubScore = max_in_score_matrix (params->scoring);
		params->hspThreshold.s = hspThresholdRatio * maxSubScore;
		params->xDrop          = 10                * maxSubScore;

		// output the scores we are trying

		print_job_header ();

		if (showAllScores)
			{
			sprintf (scoreFileId, "s%03d", trial-1);
			write_scores (scoreFileId, params->scoring,
			              /*withGapScores*/ false,
			              /*withExtras*/    true,
			              /*asInts*/        false);
			}

		if (infer_scores_dbgShowIdentity)
			printf ("===== starting iteration s%03d =====\n", trial-1);

		// perform alignment and infer scores from resulting alignments

		align_for_sub_scores (scaleTo);
		rewind_sequence_file (query);

		// if resulting score matrix is close enough to one we've seen before,
		// exit the main loop (by setting inOrbit = true)

		subScores.s1 = inferredScoring->sub['A']['A'];
		subScores.s2 = inferredScoring->sub['C']['C'];
		subScores.s3 = inferredScoring->sub['A']['C'];
		subScores.s4 = inferredScoring->sub['A']['G'];
		subScores.s5 = inferredScoring->sub['A']['T'];
		subScores.s6 = inferredScoring->sub['C']['G'];

		if (infer_scores_snoopConverge)
			{
			fprintf (stderr, "=== stats iteration %d ===\n", trial);
			fprintf (stderr, "alignments:%" PRIu64 " alignedBases:" possumFmt
			                 " refBases:" unsposFmt " secBases:" unsposFmt "\n",
			         infStats.count, infStats.coverage,
			         infStats.refBases, infStats.secBases);
			}

		if (infer_scores_watchConverge)
			{
			fprintf (stderr, "=== subs iteration %d ===\n", trial);
			fprintf (stderr, "AA:" scoreFmtSimple " CC:" scoreFmtSimple
			                " AC:" scoreFmtSimple " AG:" scoreFmtSimple
			                " AT:" scoreFmtSimple " CG:" scoreFmtSimple "\n",
			         subScores.s1, subScores.s2, subScores.s3,
			         subScores.s4, subScores.s5, subScores.s6);
			}

		for (oldTrial=trial-1 ; oldTrial>=0 ; oldTrial--)
			{
			if (close_enough_scores_6 (&subScores, &pastSubScores[oldTrial]))
				{ inOrbit = true;  break; }
			}

		pastSubScores[trial] = subScores;

		// ping-pong scoring sets

		temp            = inferredScoring;
		inferredScoring = params->scoring;
		params->scoring = temp;

		repair_scores (params->scoring, params->maskedScoring);

		if (params->showStats)
			infer_scores_show_stats_subs (params->statsFile, trial);
		}

	//////////
	// Phase II-- iterate gap scores inference
	//////////

	// copy the final substitution scores into both scoring sets

	copy_scores (/*to*/ inferredScoring, /*from*/ params->scoring);

	// determine limiting parameters, relative to the current maximum
	// substitution score

	maxSubScore = max_in_score_matrix (params->scoring);
	minSubScore = min_in_score_matrix (params->scoring);
	oneOverMaxSubScore = 1.0 / maxSubScore;
	minOverMaxSubScore = (-minSubScore) / (double) maxSubScore;

	hspThresholdRatio    = origHspThresholdRatio;
	gappedThresholdRatio = origGappedThresholdRatio;
	gapOpenRatio         = origGapOpenRatio;
	gapExtendRatio       = origGapExtendRatio;

	switch (params->ic.hspThresholdIsRatio)
		{
		case ratioNone:        hspThresholdRatio    *= oneOverMaxSubScore;  break;
		case ratioMinSubScore: hspThresholdRatio    *= minOverMaxSubScore;  break;
		}

	switch (params->ic.gappedThresholdIsRatio)
		{
		case ratioNone:        gappedThresholdRatio *= oneOverMaxSubScore;  break;
		case ratioMinSubScore: gappedThresholdRatio *= minOverMaxSubScore;  break;
		}

	switch (params->ic.gapOpenIsRatio)
		{
		case ratioNone:        gapOpenRatio         *= oneOverMaxSubScore;  break;
		case ratioMinSubScore: gapOpenRatio         *= minOverMaxSubScore;  break;
		}

	switch (params->ic.gapExtendIsRatio)
		{
		case ratioNone:        gapExtendRatio       *= oneOverMaxSubScore;  break;
		case ratioMinSubScore: gapExtendRatio       *= minOverMaxSubScore;  break;
		}

	params->hspThreshold.s     = hspThresholdRatio    * maxSubScore;
	params->gappedThreshold.s  = gappedThresholdRatio * maxSubScore;
	params->xDrop              = 10                   * maxSubScore;
	params->scoring->gapOpen   = gapOpenRatio         * maxSubScore;
	params->scoring->gapExtend = gapExtendRatio       * maxSubScore;

	gapScores.s1 = params->scoring->gapOpen;
	gapScores.s2 = params->scoring->gapExtend;
	pastGapScores[trial] = gapScores;
	inOrbit = false;

	for (trial=1 ; (!inOrbit)&&(trial<=params->ic.gapIterations) ; trial++)
		{
		// determine limiting parameters, relative to the current gap scores
		// (this setting of yDrop = O + 300E comes from BLASTZ defaults)

		params->yDrop = params->scoring->gapOpen
		              + 300 * params->scoring->gapExtend;

		// output the scores we are trying

		print_job_header ();

		if (showAllScores)
			{
			sprintf (scoreFileId, "g%03d", trial-1);
			write_scores (scoreFileId, params->scoring,
			              /*withGapScores*/ true,
			              /*withExtras*/    true,
			              /*asInts*/        false);
			}

		if (infer_scores_dbgShowIdentity)
			printf ("===== starting iteration g%03d =====\n", trial-1);

		// perform alignment and infer scores from resulting alignments

		align_for_gap_scores (scaleTo);
		rewind_sequence_file (query);

		// if resulting gap score pair is close enough to one we've seen before,
		// exit the main loop (by setting inOrbit = true)

		gapScores.s1 = inferredScoring->gapOpen;
		gapScores.s2 = inferredScoring->gapExtend;

		if (infer_scores_watchConverge)
			{
			fprintf (stderr, "=== gaps iteration %d ===\n", trial);
			fprintf (stderr, "O:" scoreFmtSimple " E:" scoreFmtSimple "\n",
			         gapScores.s1, gapScores.s2);
			}

		for (oldTrial=trial-1 ; oldTrial>=0 ; oldTrial--)
			{
			if (close_enough_scores_2 (&gapScores, &pastGapScores[oldTrial]))
				{ inOrbit = true;  break; }
			}

		pastGapScores[trial] = gapScores;

		// ping-pong scoring sets

		temp            = inferredScoring;
		inferredScoring = params->scoring;
		params->scoring = temp;

		repair_scores (params->scoring, params->maskedScoring);

		if (params->showStats)
			infer_scores_show_stats_gaps (params->statsFile, trial);
		}

	// ping-pong scoring sets (to compensate for the last ping-pong in the loop)

	temp            = inferredScoring;
	inferredScoring = params->scoring;
	params->scoring = temp;

	// output the resulting scores

	write_scores ("", inferredScoring,
	              /*withGapScores*/ (maxGapIterations>0),
	              /*withExtras*/    false,
	              /*asInts*/        params->ic.writeAsInt);

	// cleanup

	free_stats_for_inference ();

	if (params->ic.idIsPercentile)
		{
		params->minIdentity = minIdentity;
		params->maxIdentity = maxIdentity;
		}

	return inferredScoring;
	}


static int close_enough_scores_2
   (score2*	u,
	score2*	v)
	{
	score	diff;

	diff = u->s1 - v->s1;
	if (diff < -gapCloseEnough) return false;
	if (diff >  gapCloseEnough) return false;

	diff = u->s2 - v->s2;
	if (diff < -gapCloseEnough) return false;
	if (diff >  gapCloseEnough) return false;

	return true;
	}


static int close_enough_scores_6
   (score6*	u,
	score6*	v)
	{
	score	diff;

	diff = u->s1 - v->s1;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	diff = u->s2 - v->s2;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	diff = u->s3 - v->s3;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	diff = u->s4 - v->s4;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	diff = u->s5 - v->s5;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	diff = u->s6 - v->s6;
	if (diff < -subCloseEnough) return false;
	if (diff >  subCloseEnough) return false;

	return true;
	}

//----------
//
// align_for_sub_scores--
//	Perform ungapped alignment of the target sequence to every query, then
//	infer substitution scores from stats collected from those alignments.
//
//----------
//
// Arguments:
//	score	scaleTo:	The desired value for the maximum subsitution score.
//						.. If this is zero, no scaling is performed.
//
// Returns:
//	The scaling factor used to accomplish scaleTo.
//
//----------

static double align_for_sub_scores
   (score			scaleTo)
	{
	hitprocessor	hitProc;
	void*			hitProcInfo;
	double			scaleBy;

	params->chain = false;  params->gappedExtend = false; // equiv. to C=3
	set_up_hit_processor (params, false, &hitProc, &hitProcInfo);

	// align target vs all queries, collecting stats

	align_for_stats (hitProc, hitProcInfo);

	// combine stats into a single bin

	if (params->ic.idIsPercentile)
		filter_stats_by_percentile ();

	combine_binned_stats (/*merge sequences*/ true);

	// infer substitution scores (results are placed into inferredScoring)

	scaleBy = infer_substitution_scores (/*pOpen*/ 0.0, scaleTo);
	infer_scores_set_stat (scaleBy, scaleBy);

	return scaleBy;
	}

//----------
//
// align_for_gap_scores--
//	Perform gapped alignment of the target sequence to every query, then infer
//	infer gap scores from stats collected from those alignments.
//
//----------
//
// Arguments:
//	score	scaleTo:	The desired value for the maximum subsitution score.
//						.. If this is zero, no scaling is performed.
//	(additional input is implicit in pInf, sInf, q1Inf and q2Inf)
//
// Returns:
//	The scaling factor used to accomplish scaleTo.
//
//----------

static double align_for_gap_scores
   (score			scaleTo)
	{
	hitprocessor	hitProc;
	void*			hitProcInfo;
	double			scaleBy;

	params->chain = true;  params->gappedExtend = true; // equiv. to C=2
	set_up_hit_processor (params, false, &hitProc, &hitProcInfo);

	// align target vs all queries, collecting stats

	align_for_stats (hitProc, hitProcInfo);

	// combine stats into a single bin

	if (params->ic.idIsPercentile)
		filter_stats_by_percentile ();

	combine_binned_stats (/*merge sequences*/ true);

	// infer gap scores (results are placed into inferredScoring)

	scaleBy = infer_gap_scores (scaleTo);
	infer_scores_set_stat (scaleBy, scaleBy);

	return scaleBy;
	}

//----------
//
// align_for_stats--
//	Perform alignment (gapped or ungapped) of the target sequence to every
//	query, collecting stats from those alignments.
//
//----------
//
// Arguments:
//	hitprocessor processor:		Function to call for each hit to determine if
//								.. it is 'good enough'.
//	void*		processorInfo:	A value to pass thru with each call to
//								.. processor.
//
// Returns:
//	(nothing)
//
//----------

static void align_for_stats
   (hitprocessor	hitProc,
	void*			hitProcInfo)
	{
	int				collectHspsFromBoth;	// collect HSPs from both strands
											// .. before gapped stage
	int				hspsAreAdaptive;		// adaptive HSP scoring threshold
											// .. is being used
	int				abortQuery;


	hspsAreAdaptive = (params->hspThreshold.t != 'S');
	collectHspsFromBoth = hspsAreAdaptive;	// $$$ consider other conditions
											// $$$ .. used in mainline lastz

	// align target vs all queries, collecting stats

	erase_stats_for_inference ();

	while (load_sequence (query))
		{
		if (query->len == 0) continue;
		if (params->minMatchCountRatio != 0)
			params->minMatchCount = (u32) ceil (query->trueLen * params->minMatchCountRatio);
		if (params->whichStrand < 0)
			rev_comp_sequence (query, params->scoring->qToComplement);
		abortQuery = !start_one_strand (target, targPositions, query,
		                                /* empty anchors */ true,
		                                /* prev anchor count */ 0,
		                                hitProc, hitProcInfo);
		if (abortQuery) continue;
		if (!collectHspsFromBoth)
			finish_one_strand (target, targetRev, targPositions, query, NULL,
			                   traceback, NULL);

		if (params->whichStrand > 0)
			{
			rev_comp_sequence (query, params->scoring->qToComplement);
			abortQuery = !start_one_strand (target, targPositions, query,
			                                /* empty anchors */ !collectHspsFromBoth,
			                                /* prev anchor count */ 0,
			                                hitProc, hitProcInfo);
			if (!abortQuery)
				{
				rev_comp_sequence (query, params->scoring->qToComplement);
				continue;
				}
			if (collectHspsFromBoth) split_anchors (query->revCompFlags);
			finish_one_strand (target, targetRev, targPositions, query, NULL,
			                   traceback, NULL);
			if (collectHspsFromBoth)
				{
				swap_anchor_sets ();
				// we have to reverse query for subsequent call to finish_one_strand()
				rev_comp_sequence (query, params->scoring->qToComplement);
				}
			}

		if (collectHspsFromBoth)
			finish_one_strand (target, targetRev, targPositions, query, NULL,
			                   traceback, NULL);
		}

	}

//----------
//
// infer_substitution_scores--
//	Infer log odds substitution scoring matrix as per [1];  one addition we
//	make to [1] is that we involve pOpen in each score to properly account for
//	gaps in the three state pair FSA from [2];  for more information see the
//	description in infer_gap_scores().
//
//----------
//
// Arguments:
//	double	pOpen:		The gap-opening probability.
//	score	scaleTo:	The desired value for the maximum subsitution score.
//						.. All scores in the scoring matrix will be scaled by
//						.. the same amount so that the maximum is this value.
//						.. If this is zero, no scaling is performed.
//	(additional input is implicit in infStats.subs)
//
// Returns:
//	The scaling factor used to accomplish scaleTo.  Additional output is
//	implicit in inferredScoring, and internal state is saved in pInf, sInf,
//	q1Inf and q2Inf (so that subsequent calls to infer_gap_scores can make
//	use of it).
//
//----------

static double pInf[4][4], sInf[4][4];
static double qInf1[4], qInf2[4];

static void   sub_probs_to_log_scores   (double pOpen);
static double log_scores_to_scoring_set (score scaleTo);

//--- infer_substitution_scores

static double infer_substitution_scores
   (double	pOpen,
	score	scaleTo)
	{
	unspos	m[4][4];
	unspos	n1[4], n2[4];
	unspos	n;
	double	npairs;
	int		x, y, xx, yy;
	double	scaleBy;

	// collect column counts from the alignment stats

	for (x=0 ; x<4 ; x++)
		{
		n1[x] = 0;
		n2[x] = 0;
		for (y=0 ; y<4 ; y++)
			m[x][y] = 0;
		}

	for (x=0 ; x<4 ; x++)
			for (y=0 ; y<4 ; y++)
		{
		// "observe(n,x,y)"

		n = infStats.subs[x][y];

		xx = x;  yy = y;
		m[xx][yy] += n;
		n1[xx]    += n;
		n2[yy]    += n;

		xx = bits_to_complement[x];					// (for strand symmetry)
		yy = bits_to_complement[y];
		m[xx][yy] += n;
		n1[xx]    += n;
		n2[yy]    += n;

		xx = y;  yy = x;							// (for species symmetry)
		m[xx][yy] += n;
		n1[xx]    += n;
		n2[yy]    += n;

		xx = bits_to_complement[y];					// (for both strand and
		yy = bits_to_complement[x];					//  .. species symmetry)
		m[xx][yy] += n;
		n1[xx]    += n;
		n2[yy]    += n;
		}

#ifdef collect_stats
	for (x=0 ; x<4 ; x++) for (y=0 ; y<4 ; y++) infer_scores_set_stat (subs[x][y], infStats.subs[x][y]);
	for (x=0 ; x<4 ; x++)                       infer_scores_set_stat (n1[x],      n1[x]);
	for (y=0 ; y<4 ; y++)                       infer_scores_set_stat (n2[y],      n2[y]);
	for (x=0 ; x<4 ; x++) for (y=0 ; y<4 ; y++) infer_scores_set_stat (m[x][y],    m[x][y]);
#endif // collect_stats

	// validate the expected symmetry

	if ((n1[3] != n1[0])
	 || (n1[2] != n1[1])
	 || (n2[3] != n2[0])
	 || (n2[2] != n2[1])
	 || (m[3][3] != m[0][0])
	 || (m[2][2] != m[1][1])
	 || (m[1][0] != m[0][1])
	 || (m[2][3] != m[0][1])
	 || (m[3][2] != m[0][1])
	 || (m[2][0] != m[0][2])
	 || (m[1][3] != m[0][2])
	 || (m[3][1] != m[0][2])
	 || (m[3][0] != m[0][3])
	 || (m[2][1] != m[1][2]))
		suicidef ("internal error: non-symmetry in infer_substitution_scores\n"
		          "  n1:   %7d %7d %7d %7d\n"
		          "  n2:   %7d %7d %7d %7d\n"
		          "  m[0]: %7d %7d %7d %7d\n"
		          "  m[1]: %7d %7d %7d %7d\n"
		          "  m[2]: %7d %7d %7d %7d\n"
		          "  m[3]: %7d %7d %7d %7d\n",
		          n1[0],   n1[1],   n1[2],   n1[3],
		          n2[0],   n2[1],   n2[2],   n2[3],
		          m[0][0], m[0][1], m[0][2], m[0][3],
		          m[1][0], m[1][1], m[1][2], m[1][3],
		          m[2][0], m[2][1], m[2][2], m[2][3],
		          m[3][0], m[3][1], m[3][2], m[3][3]);

	// infer log odds scores from column counts
	//
	// nota bene: because of the symmetry folding performed above, we could
	//            simplify some of these counts, etc. (e.g. p(G) == p(T));  but
	//            the computational gain from doing so is inconsequential, and
	//            not doing so makes the correllation between the code and [1]
	//            more obvious

	npairs = (double) (n1[0] + n1[1] + n1[2] + n1[3]);

	for (x=0 ; x<4 ; x++)
		{
		if ((n1[x] == 0) || (n2[x] == 0))
			suicidef ("internal error in infer_substitution_scores"
			          " n1[%c] or n2[%c] is zero",
			          bits_to_nuc[x], bits_to_nuc[x]);
		qInf1[x] = n1[x] / npairs;
		qInf2[x] = n2[x] / npairs;
		for (y=0 ; y<4 ; y++)
			pInf[x][y] = m[x][y] / npairs;
		}

	sub_probs_to_log_scores (pOpen);

	// copy scores into inferredScoring, scaling if desired

	scaleBy = log_scores_to_scoring_set (scaleTo);

	inferredScoring->gapOpen   = 0;
	inferredScoring->gapExtend = 0;

	return scaleBy;
	}

//--- sub_probs_to_log_scores
//    (called by both infer_substitution_scores and infer_gap_scores)

static void sub_probs_to_log_scores (double pOpen)
	{
	double	overLog2 = 1 / log(2.0);
	int		x, y;

	for (x=0 ; x<4 ; x++)
			for (y=0 ; y<4 ; y++)
		{
		if (pInf[x][y] == 0)
			suicidef ("internal error in infer_substitution_scores"
			          " s[%c][%c] = -infinity",
			          bits_to_nuc[x], bits_to_nuc[y]);

		sInf[x][y] = log(pInf[x][y]/(qInf1[x]*qInf2[y])) * overLog2;
		if (pOpen != 0)
			sInf[x][y] += log(1-2*pOpen) * overLog2;
		}
	}

//--- log_scores_to_scoring_set
//    (called by both infer_substitution_scores and infer_gap_scores)

static double log_scores_to_scoring_set (score scaleTo)
	{
	int		x, y;
	u8		nuc1, nuc2;
	double	scaleBy;

	if (scaleTo <= 0)
		scaleBy = 1.0;
	else
		{
		double maxS = sInf[0][0];
		for (x=0 ; x<4 ; x++)
				for (y=0 ; y<4 ; y++)
				if (sInf[x][y] > maxS)
			{ maxS = sInf[x][y]; }

		scaleBy = ((double) scaleTo) / maxS;
		}

	for (x=0 ; x<4 ; x++)
		{
		nuc1 = (u8) bits_to_nuc[x];
		for (y=0 ; y<4 ; y++)
			{
			nuc2 = (u8) bits_to_nuc[y];
#if (scoreType == 'I')
			inferredScoring->sub[nuc1][nuc2] = round_score (scaleBy * sInf[x][y]);
#else
			inferredScoring->sub[nuc1][nuc2] = scaleBy * sInf[x][y];
#endif
			}
		}

	return scaleBy;
	}

//----------
//
// infer_gap_scores--
//	Infer log odds scores for gap open and extend.
//
// For the underlying gap open probability, we have:
//		avgSeg = 1/(2*p(stop))
//		=> p(stop) = 1/(2*avgSeg)
//		=> p(gap open) = 1 - 1/(2*avgSeg)
// The reason for 2 in the denominator is because a gap could occur in either
// sequence, so the real p(stopping a segment) = 2*p(stop).
//
// For the underlying gap extend probability, we have:
//		avgGap = 1/p(stop)
//		       = 1/(1-p(gap extend))
//		=> 1-p(gap extend) = 1/avgGap
//		=> p(gap extend) = 1-(1/avgGap)
//
// For the pair FSA in [2], the log odds scores are then
//		s'_xy     = log p_x,y + log(1-2p_open)
//		s'_open   = log p_open
//		s'_extend = log p_extend
//
// However there are two modifications that must be made.  First, the initial
// pair coming out of a gap whould be scored as log p_x,y + log(1-p_extend), so
// using s'_xy we have underscored it by log(1-p_extend) - log(1-2p_open).  We
// compensate for this by increasing the gap open score.  Second, an extra gap
// extend penalty is charged as the gap is being opened, so we subtract this
// from our gap open penalty.  Making these adjustments, the inferred scores
// are
//		s_xy     = log p_x,y + log(1-2p_open)
//		s_open   = log p_open + log(1-p_extend) - log(1-2p_open) - log p_extend
//		s_extend = log p_extend
//
//----------
//
// Arguments:
//	score	scaleTo:	The desired value for the maximum subsitution score.
//						.. See the description in infer_substitution_scores()
//						.. for further info.  Scaling is performed again as
//						.. part of this routine, since a non-zero pOpen changes
//						.. the substitution scores.
//	(other input is implicitly infStats.refGaps and infStats.segments)
//
// Returns:
//	The scaling factor used to accomplish scaleTo.  Additional output is
//	implicit in inferredScoring.
//
//----------

static double infer_gap_scores
   (score	scaleTo)
	{
	double	overLog2 = 1 / log(2.0);
	double	avgGap, avgSeg;
	double	pOpen, pExtend;
	double	sOpen, sExtend;
	double	scaleBy;

	pOpen = pExtend = sOpen = sExtend = 0;	// (placate compiler)

	if (number_of_instances (infStats.refGaps) == 0)
		suicide ("internal error in infer_gap_scores: no gaps");

	avgGap = average_length (infStats.refGaps);
	avgSeg = average_length (infStats.segments);

	infer_scores_set_stat (averageGapLength,     avgGap);
	infer_scores_set_stat (averageSegmentLength, avgSeg);

	// infer gap extend

	if (avgGap < 0)
		suicide ("internal error in infer_gap_scores: average gap doesn't exist");
	else if (avgGap == 1)
		suicide ("internal error in infer_gap_scores: average gap is 1");
	else
		{
		pExtend = 1 - (1/avgGap);
		sExtend = log(pExtend) * overLog2;
		infer_scores_set_stat (pExtend, pExtend);
		infer_scores_set_stat (sExtend, sExtend);
		}

	// infer gap open

	if (avgSeg < 0)
		suicide ("internal error in infer_gap_scores: average segment doesn't exist");
	else
		pOpen = 1 / (2*avgSeg);
		sOpen = (log(pOpen) - log(1-2*pOpen) + log(1-pExtend) - log(pExtend))
		      * overLog2;

		if (sOpen + sExtend >= 0)
			suicidef ("internal inconsistency, gap open \"reward\" in infer_gap_scores\n"
		              "(avgGap=%f pExtend=%f sExtend=%f avgSeg=%f pOpen=%f)\n"
		              "(+log(pOpen)=%f -log(1-2*pOpen)=%f +log(1-pExtend)=%f -log(pExtend)=%f)",
		              avgGap,pExtend,sExtend,avgSeg,pOpen,
		              log(pOpen)/log(2),-log(1-2*pOpen)/log(2),
		              log(1-pExtend)/log(2),-log(pExtend)/log(2));

		infer_scores_set_stat (pOpen, pOpen);
		infer_scores_set_stat (sOpen, sOpen);

	// recompute log odds substitution scores now that we have pOpen, and
	// rescale them

	sub_probs_to_log_scores (pOpen);
	scaleBy = log_scores_to_scoring_set (scaleTo);

	// scale and copy gap scores into inferredScoring

#if (scoreType == 'I')
	inferredScoring->gapOpen   = round_score (scaleBy * (-sOpen));
	inferredScoring->gapExtend = round_score (scaleBy * (-sExtend));
#else
	inferredScoring->gapOpen   = scaleBy * (-sOpen);
	inferredScoring->gapExtend = scaleBy * (-sExtend);
#endif

	return scaleBy;
	}

//----------
//
// copy_scores--
//	Copy one set of scores into the other (only uppercase substitution scores
//	are copied).
//
//----------
//
// Arguments:
//	scoreset*	dst:	The score set to copy into.
//	scoreset*	src:	The score set to copy from.
//
//----------
//
// Returns:
//	(nothing)
//
//----------

static void copy_scores
   (scoreset*	dst,
	scoreset*	src)
	{
	int		x, y;
	u8		nuc1, nuc2;

	for (x=0 ; x<4 ; x++)
		{
		nuc1 = (u8) bits_to_nuc[x];
		for (y=0 ; y<4 ; y++)
			{
			nuc2 = (u8) bits_to_nuc[y];
			dst->sub[nuc1][nuc2] = src->sub[nuc1][nuc2];
			}
		}

	}

//----------
//
// repair_scores--
//	Fix a set of inferred scores.
//
// The functions infer_substitution_scores() and infer_gap_scores() set values
// in a scores set, but do not propagate those scores to lower case.  Nor do
// they update the repeat-masked version of the score set.  This routine
// corrects for those shortcomings.
//
//----------
//
// Arguments:
//	scoreset*	scoring:		The score set, as created by either of the
//								.. infer_xxx_scores() functions.  This will be
//								.. modified by this routine.
//	scoreset*	maskedScoring:	The score set to contain a repeat-masked
//								.. version of the scoring matrix.
//
//----------
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free() should be used for this purpose.
//
//----------
//
// Notes:
//	(1)	In the resulting scoring matrix, upper and lower case characters are
//		are considered identical, so entries for lower case are copied from
//		upper case.
//
//----------

static void repair_scores
   (scoreset*	scoring,
	scoreset*	masked)
	{
	int		x, y;
	u8		nuc1, nuc2, nuc1low, nuc2low;
	int		c;
	score	sub, worstSub;

	worstSub = 0;

	for (x=0 ; x<4 ; x++)
		{
		nuc1    = (u8) bits_to_nuc[x];
		nuc1low = dna_tolower (nuc1);
		for (y=0 ; y<4 ; y++)
			{
			nuc2    = (u8) bits_to_nuc[y];
			nuc2low = dna_tolower (nuc2);
			sub = scoring->sub[nuc1][nuc2];
			scoring->sub[nuc1low][nuc2   ] = sub;
			scoring->sub[nuc1   ][nuc2low] = sub;
			scoring->sub[nuc1low][nuc2low] = sub;
			masked ->sub[nuc1   ][nuc2   ] = sub;
			if (sub < worstSub) worstSub = sub;
			}
		}

	for (x=0 ; x<4 ; x++)
		{
		nuc1    = (u8) bits_to_nuc[x];
		nuc1low = dna_tolower (nuc1);
		scoring->sub[nuc1   ]['N'    ] = worstSub;
		scoring->sub[nuc1low]['N'    ] = worstSub;
		scoring->sub[nuc1   ]['n'    ] = worstSub;
		scoring->sub[nuc1low]['n'    ] = worstSub;
		scoring->sub['N'    ][nuc1   ] = worstSub;
		scoring->sub['N'    ][nuc1low] = worstSub;
		scoring->sub['n'    ][nuc1   ] = worstSub;
		scoring->sub['n'    ][nuc1low] = worstSub;
		}

	scoring->sub['N']['N'] = worstSub;
	scoring->sub['N']['n'] = worstSub;
	scoring->sub['n']['N'] = worstSub;
	scoring->sub['n']['n'] = worstSub;

	// make sure scores for row and column zero are very very bad

	for (c=0 ; c<256 ; c++)
		scoring->sub[0][c] = scoring->sub[c][0] = veryBadScore;
	}

//----------
//
// write_scores--
//	Write the current scoring set to a file.
//
//----------
//
// Arguments:
//	char*		fileId:			A string to include in the file name, to
//								.. identify this scoring set.
//	scoreset*	ss:				The scoring set to write.
//	int			withGapScores:	true => write gap scores too
//	int			withExtras:		true => write extra values that are related to
//								        .. the scoring set, as comments
//	int			asInts:			write scores as integers regardless of scoreType
//
// Returns:
//  (nothing)
//
//----------

static void write_scores
   (char*		fileId,
	scoreset*	ss,
	int			withGapScores,
	int			withExtras,
	int			asInts)
	{
	char		name[201];
	int			replaced;
	FILE*		f;

	if (params->ic.inferFilename == NULL)
		f = stdout;
	else
		{
		strcpy (name, params->ic.inferFilename);

		replaced = false;
		if ((fileId == NULL) || (fileId[0] == 0))
			{
			if (!replaced) replaced = string_replace (name, sizeof(name), "_%s", fileId);
			if (!replaced) replaced = string_replace (name, sizeof(name), ".%s", fileId);
			}
		if (!replaced) replaced = string_replace (name, sizeof(name), "%s",  fileId);

		if ((!replaced) && (strstr (name, "%s") != NULL))
			suicidef ("unable to perform name substitution, try a shorter name than"
					  " %s", name);

		f = fopen_or_die (name, "wt");
		}

	// write the scores to the file

	if (asInts)
		write_score_set_as_ints (f, name, ss, withGapScores);
	else
		write_score_set         (f, name, ss, withGapScores);

	// write other useful info (as comments)

	if (withExtras)
		{
		                   fprintf (f, "\n");
		                   fprintf (f, "# hsp_threshold    = %s\n", score_thresh_to_string (&params->hspThreshold));
		if (withGapScores) fprintf (f, "# gapped_threshold = %s\n", score_thresh_to_string (&params->gappedThreshold));
		                   fprintf (f, "# x_drop           = " scoreFmtSimple "\n", params->xDrop);
		if (withGapScores) fprintf (f, "# y_drop           = " scoreFmtSimple "\n", params->yDrop);
		}

	fclose_if_valid (f);
	}

//----------
//
// init_stats_for_inference, erase_stats_for_inference, free_stats_for_inference--
//	Manage the by-percent-identity inference stats.
//
//----------

static void init_stats_for_inference
   (arg_dont_complain(seq* seq1),
	arg_dont_complain(seq* seq2))
	{
	u32	bin;

	init_stats (&infStats);
	for (bin=0 ; bin<=numIdentityBins ; bin++)
		init_stats (&infStatsByPctId[bin]);
	}


static void erase_stats_for_inference
   (void)
	{
	u32 bin;

	erase_stats (&infStats);	// (probably not necessary)
	for (bin=0 ; bin<=numIdentityBins ; bin++)
		erase_stats (&infStatsByPctId[bin]);
	}


static void free_stats_for_inference
   (void)
	{
	u32 bin;

	free_stats (&infStats);

	for (bin=0 ; bin<=numIdentityBins ; bin++)
		free_stats (&infStatsByPctId[bin]);
	}

//----------
//
// gather_stats_from_align_list--
//	Collect inference stats from a list of gapped alignments.
//
//----------
//
// Arguments:
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//
// Returns:
//	(nothing)
//
//----------

void gather_stats_from_align_list
   (alignel*	alignList,
	seq*		seq1,
	seq*		seq2)
	{
	alignel*	a;
	unspos		numer, denom;
	u32			bin;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		alignment_identity (seq1, seq2, a, &numer, &denom);
		bin = identity_bin (numer, denom);
		infStatsByPctId[bin].count++;
		infStatsByPctId[bin].coverage += denom;
		accumulate_stats_from_align (seq1, a->beg1-1, a->end1,
		                             seq2, a->beg2-1, a->end2,
		                             a->script, &infStatsByPctId[bin]);

		if (infer_scores_dbgShowIdentity)
			{
			// nota bene: positions written as 1-based
			printf (unsposSlashFmt " identity=" unsposSlashFmt
			        " (bin as " identityBinFormat ")\n",
			        a->beg1, a->beg2, numer, denom, bin_to_identity (bin));
			}
		}

	}

//----------
//
// gather_stats_from_match--
//	Collect inference stats from a single ungapped alignment.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos2:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the match.
//
// Returns:
//	(nothing)
//
//----------

void gather_stats_from_match
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	unspos	numer, denom;
	u32		bin;

	segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
	bin = identity_bin (numer, denom);
	infStatsByPctId[bin].count++;
	infStatsByPctId[bin].coverage += denom;
	accumulate_stats_from_match (seq1, pos1, seq2, pos2, length,
	                             &infStatsByPctId[bin]);
	}

//----------
//
// filter_stats_by_percentile--
//	Discard inference stats outside the desired percentile of percent-identity.
//
//----------
//
// Arguments:
//	(none;  input is implicitly infStatsByPctId[*])
//
// Returns:
//	(nothing;  output is implicitly infStatsByPctId[*])
//
//----------

static void filter_stats_by_percentile
   (void)
	{
	static const u32 noBin = (u32) -1;
	u32		bin, minBin;
	possum	cov, covTotal, covLo, covHi;

	// convert the percentiles to a range of coverage counts

	covTotal = 0;
	minBin   = noBin;
	for (bin=0 ; bin<=numIdentityBins ; bin++)
		{
		cov = infStatsByPctId[bin].coverage;
		if (cov == 0) continue;
		covTotal += cov;
		if (minBin == noBin) minBin = bin;
		}
	if (minBin == noBin) minBin = numIdentityBins;

	covLo = (covTotal * minIdentity) + 0.5;
	covHi = (covTotal * maxIdentity) + 0.5;

	infer_scores_set_stat (coverageTotal, covTotal);
	infer_scores_set_stat (coverageLow,   covLo);
	infer_scores_set_stat (coverageHigh,  covHi);

	// discard any bins outside the range

	for (bin=numIdentityBins+1 ; bin>0 ; )
		{
		bin--;
		cov = infStatsByPctId[bin].coverage;
		if (cov == 0) continue;
		erase_stats (&infStatsByPctId[bin]);
		covTotal -= cov;
		infer_scores_set_stat (highIdentityBin, bin);
		if (covTotal <= covHi) break;
		}

	covTotal = 0;
	for (bin=minBin ; bin<=numIdentityBins ; bin++)
		{
		cov = infStatsByPctId[bin].coverage;
		if (cov == 0) continue;
		erase_stats (&infStatsByPctId[bin]);
		covTotal += cov;
		infer_scores_set_stat (lowIdentityBin, bin);
		if (covTotal >= covLo) break;
		}

	if (infer_scores_dbgShowIdentity)
			for (bin=minBin ; bin<=numIdentityBins ; bin++)
		{
		cov = infStatsByPctId[bin].coverage;
		if (cov == 0) continue;
		printf ("bin: " identityBinFormat " cov=" possumFmt "\n",
				bin_to_identity (bin), cov);
		}
	}

//----------
//
// combine_binned_stats--
//	Combine inference stats into a single bin.
//
// Stats from all bins in infStatsByPctId are combined into infStats.
//
//----------
//
// Arguments:
//	int	mergeSequences:	true => merge stats for reference and secondary.
//	(other input is implicitly infStatsByPctId[*])
//
// Returns:
//	(nothing;  output is implicitly infStats)
//
//----------

static void combine_binned_stats
   (int			mergeSequences)
	{
	u32			bin;
	infstats*	inf;
	u8			c1, c2;

	erase_stats (&infStats);

	for (bin=0 ; bin<=numIdentityBins ; bin++)
		{
		inf = &infStatsByPctId[bin];
		if ((inf == NULL) || (inf->count == 0)) continue;

		infStats.count    += inf->count;
		infStats.coverage += inf->coverage;
		infStats.refBases += inf->refBases;
		infStats.secBases += inf->secBases;

		for (c1=0 ; c1<4 ; c1++)
			{
			infStats.refBkgd[c1] += inf->refBkgd[c1];
			infStats.secBkgd[c1] += inf->secBkgd[c1];
			for (c2=0 ; c2<4 ; c2++)
				infStats.subs[c1][c2] += inf->subs[c1][c2];
			}

		add_lengths_to_distribution (inf->refBlocks, &infStats.refBlocks);
		add_lengths_to_distribution (inf->refGaps,   &infStats.refGaps);
		add_lengths_to_distribution (inf->refRuns,   &infStats.refRuns);
		add_lengths_to_distribution (inf->segments,  &infStats.segments);

		if (mergeSequences)
			{
			add_lengths_to_distribution (inf->secBlocks, &infStats.refBlocks);
			add_lengths_to_distribution (inf->secGaps,   &infStats.refGaps);
			add_lengths_to_distribution (inf->secRuns,   &infStats.refRuns);
			}
		else
			{
			add_lengths_to_distribution (inf->secBlocks, &infStats.secBlocks);
			add_lengths_to_distribution (inf->secGaps,   &infStats.secGaps);
			add_lengths_to_distribution (inf->secRuns,   &infStats.secRuns);
			}
		}

	}

//----------
//
// infererence statistics sets routines--
//
//----------

static void init_stats
   (infstats*	inf)
	{
	inf->refBlocks = init_length_distribution (0);
	inf->secBlocks = init_length_distribution (0);
	inf->refGaps   = init_length_distribution (0);
	inf->secGaps   = init_length_distribution (0);
	inf->refRuns   = init_length_distribution (0);
	inf->secRuns   = init_length_distribution (0);
	inf->segments  = init_length_distribution (0);

	erase_stats (inf);
	}

static void erase_stats
   (infstats*	inf)
	{
	u8			c1, c2;

	inf->count    = 0;
	inf->coverage = 0;
	inf->refBases = 0;
	inf->secBases = 0;

	for (c1=0 ; c1<4 ; c1++)
		{
		inf->refBkgd[c1] = 0;
		inf->secBkgd[c1] = 0;
		for (c2=0 ; c2<4 ; c2++)
			inf->subs[c1][c2] = 0;
		}

	erase_length_distribution (inf->refBlocks);
	erase_length_distribution (inf->secBlocks);
	erase_length_distribution (inf->refGaps);
	erase_length_distribution (inf->secGaps);
	erase_length_distribution (inf->refRuns);
	erase_length_distribution (inf->secRuns);
	erase_length_distribution (inf->segments);
	}


static void free_stats
   (infstats*	inf)
	{
	free_length_distribution (inf->refBlocks);
	free_length_distribution (inf->secBlocks);
	free_length_distribution (inf->refGaps);
	free_length_distribution (inf->secGaps);
	free_length_distribution (inf->refRuns);
	free_length_distribution (inf->secRuns);
	free_length_distribution (inf->segments);
	}


static void accumulate_stats_from_align
   (seq*		seq1,
	unspos		beg1,
	unspos		end1,
	seq*		seq2,
	unspos		beg2,
	unspos		end2,
	editscript*	script,
	infstats*	inf)
	{
	unspos		height, width, i, j, prevI, prevJ;
	u32			opIx;
    unspos		run, refRun, secRun, indelLen, indelBases;
	u8*			s1, *s2;
	u8			c1,  c2;
	s8			cc1, cc2;
	unspos		denom;
	unspos		count, pairCount[4][4];
	unspos		ix;

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	add_length_to_distribution (height, &inf->refBlocks);
	add_length_to_distribution (width,  &inf->secBlocks);

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		pairCount[c1][c2] = 0;

	refRun = secRun = 0;
	opIx = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		prevI = i;  prevJ = j;
		run = edit_script_run_of_subs (script, &opIx);
		i      += run; j      += run;
		refRun += run; secRun += run;

		if (run > 0)
			{
			denom = count_substitutions (seq1, beg1-1+prevI,
			                             seq2, beg2-1+prevJ,
			                             run, pairCount);
			if (denom != 0)
				{
				inf->refBases += denom;
				inf->secBases += denom;
				add_length_to_distribution (denom, &inf->segments);
				}
			}

		if ((i < height) || (j < width))
			{
			prevI = i;  prevJ = j;
			edit_script_indel_len (script, &opIx, &i, &j);
			if (j != prevJ) // (deletion from reference sequence)
				{
				indelLen = j - prevJ;
				add_length_to_distribution (indelLen, &inf->refGaps);
				if (refRun > 0)
					{
					add_length_to_distribution (refRun, &inf->refRuns);
					refRun = 0;
					}
				indelBases = 0;
				s2 = seq2->v + beg2-1+prevJ;
				for (ix=0 ; ix<indelLen ; ix++)
					{
					cc2 = nuc_to_bits[*(s2++)];
					if (cc2 >= 0) { inf->secBkgd[(u8)cc2]++;  indelBases++; }
					}
				secRun        += indelBases;
				inf->secBases += indelBases;
				}
			if (i != prevI) // (deletion from second sequence)
				{
				indelLen = i - prevI;
				add_length_to_distribution (indelLen, &inf->secGaps);
				if (secRun > 0)
					{
					add_length_to_distribution (secRun, &inf->secRuns);
					secRun = 0;
					}
				indelBases = 0;
				s1 = seq1->v + beg1-1+prevI;
				for (ix=0 ; ix<indelLen ; ix++)
					{
					cc1 = nuc_to_bits[*(s1++)];
					if (cc1 >= 0) { inf->refBkgd[(u8)cc1]++;  indelBases++; }
					}
				refRun        += indelBases;
				inf->refBases += indelBases;
				}
			}
		}

	if (refRun > 0) add_length_to_distribution (refRun, &inf->refRuns);
	if (secRun > 0) add_length_to_distribution (secRun, &inf->secRuns);

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		{
		count = pairCount[c1][c2];
		inf->refBkgd[c1]  += count;
		inf->secBkgd[c2]  += count;
		inf->subs[c1][c2] += count;
		}
	}


static void accumulate_stats_from_match
   (seq*		seq1,
	unspos		pos1,
	seq*		seq2,
	unspos		pos2,
	unspos		length,
	infstats*	inf)
	{
	u8			c1, c2;
	unspos		denom;
	unspos		count, pairCount[4][4];

	// count substitutions

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		pairCount[c1][c2] = 0;

	denom = count_substitutions (seq1, pos1, seq2, pos2, length, pairCount);

	// collect stats

	inf->refBases += denom;
	inf->secBases += denom;
	add_length_to_distribution (denom, &inf->refBlocks);
	add_length_to_distribution (denom, &inf->secBlocks);
	add_length_to_distribution (denom, &inf->segments);

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		{
		count = pairCount[c1][c2];
		inf->refBkgd[c1]  += count;
		inf->secBkgd[c2]  += count;
		inf->subs[c1][c2] += count;
		}

	}

//----------
//
// miscellaneous printing routines--
//
//----------

static void print_bkgd_stats
   (FILE*	f,
	char*	s,
	unspos	bkgd[4])
	{
	u8		c;
	u8		nuc;

	fprintf (f, "    %-7s", s);

	for (c=0 ; c<4 ; c++)
		{
		nuc = (u8) bits_to_nuc[c];
		fprintf (f, " %c:" unsposFmt, nuc, bkgd[c]);
		}
	fprintf (f, "\n");
	}


static void print_subs_stats
   (FILE*	f,
	unspos	subs[4][4])
	{
	u8		c1,   c2;
	u8		nuc1, nuc2;

	for (c1=0 ; c1<4 ; c1++)
		{
		nuc1 = bits_to_nuc[c1];
		fprintf (f, "    ");
		for (c2=0 ; c2<4 ; c2++)
			{
			nuc2 = bits_to_nuc[c2];
			if (c2 != 0) fprintf (f, " ");
			fprintf (f, "%c%c:" unsposFmt, nuc1, nuc2, subs[c1][c2]);
			}
		fprintf (f, "\n");
		}
	}


static void print_blocks_stats
   (FILE*	f,
	char*	s,
	distn*	blocks)
	{
	fprintf (f, "    blocks in %s\n", s);
	print_length_distribution (f, "    ", blocks);
	}


static void print_gaps_stats
   (FILE*	f,
	char*	s,
	distn*	gaps)
	{
	fprintf (f, "    gaps in %s\n", s);
	print_length_distribution (f, "    ", gaps);
	}


static void print_runs_stats
   (FILE*	f,
	char*	s,
	distn*	runs)
	{
	fprintf (f, "    runs in %s\n", s);
	print_length_distribution (f, "    ", runs);
	}


static void print_segments_stats
   (FILE*	f,
	distn*	segments)
	{
	fprintf (f, "    segments\n");
	print_length_distribution (f, "    ", segments);
	}

//----------
//
// length distribution routines--
//
//----------

static distn* init_length_distribution
   (u32		numEntries)
	{
	u32		bytesMain, bytesHeap;
	distn*	d;

	// if there are no entries desired, we won't allocate any memory until
	// later (if and when anything is added to the distribution)

	if (numEntries == 0) return NULL;

	// allocate

	bytesMain = round_up_16 (sizeof(distn));
	bytesHeap = round_up_16 (numEntries * sizeof(dpair));
	d = malloc_or_die ("infer_stats distribution", bytesMain + bytesHeap);

	// hook up the internal array

	d->items = (dpair*) (((char*) d) + bytesMain);

	// initialize

	d->size = bytesHeap / sizeof(dpair);
	d->len  = 0;

	return d;
	}


static void erase_length_distribution
   (distn*	d)
	{
	if (d == NULL) return;
	d->len = 0;
	}


static void free_length_distribution
   (distn*	d)
	{
	free_if_valid ("infer_stats distribution", d);
	}


static void add_lengths_to_distribution
   (distn*	src,
	distn**	_dst)
	{
	distn*	dst = *_dst;
	u32		ix, iy;
	unspos	length;
	u64		count;
	u32		newEntries;
	u32		bytesMain, bytesHeap;

	if (src == NULL) return;

	// if the distribution hasn't been allocated yet, this amounts to a copy
	// operation

	if (dst == NULL)
		{
		newEntries = src->len;
		bytesMain  = round_up_16 (sizeof(distn));
		bytesHeap  = round_up_16 (newEntries * sizeof(dpair));
		(*_dst) = dst = malloc_or_die ("add_lengths_to_distribution",
		                               bytesMain + bytesHeap);

		dst->items = (dpair*) (((char*) dst) + bytesMain);
		dst->size  = bytesHeap / sizeof(dpair);
		dst->len   = newEntries;
		memcpy (/*to*/ dst->items, /*from*/ src->items, 
		        /*how much*/ newEntries * sizeof(dpair));
		return;
		}

	// otherwise, consider lengths one at a time

	for (iy=0 ; iy<src->len ; iy++)
		{
		length = src->items[iy].length;
		count  = src->items[iy].count;

		// locate this length;  if we've seen it before, just update the count

		for (ix=0 ; ix<dst->len ; ix++)
			{ if (dst->items[ix].length == length) break; }

		if (ix < dst->len)
			{ dst->items[ix].count += count;  continue; }

		// length wasn't found;  make sure there's enough room, then add an entry

		if (dst->len >= dst->size)
			{
			newEntries = 4*dst->size/3;
			if (src->size > newEntries) newEntries = src->size;

			bytesMain = round_up_16 (sizeof(distn));
			bytesHeap = round_up_16 (newEntries * sizeof(dpair));
			(*_dst) = dst = realloc_or_die ("add_lengths_to_distribution",
											dst, bytesMain + bytesHeap);

			dst->items = (dpair*) (((char*) dst) + bytesMain);
			dst->size  = bytesHeap / sizeof(dpair);
			}

		ix = dst->len++;
		dst->items[ix].length = length;
		dst->items[ix].count  = count;
		}

	}


static void add_length_to_distribution
   (unspos	length,
	distn**	_d)
	{
	distn*	d = *_d;
	u32		ix;
	u32		newEntries, len;
	u32		bytesMain, bytesHeap;

	// if the distribution hasn't been allocated yet, go do so

	if (d == NULL)
		{ newEntries = 1000;  len = 0;  goto alloc_distn; }

	// locate length

	for (ix=0 ; ix<d->len ; ix++)
		{ if (d->items[ix].length == length) break; }

	if (ix < d->len)
		{ d->items[ix].count++;  return; }

	// length wasn't found;  make sure there's enough room, then add an entry

	if (d->len >= d->size)
		{
		newEntries = 4*d->size/3;
		len        = d->len;

	alloc_distn:
		bytesMain = round_up_16 (sizeof(distn));
		bytesHeap = round_up_16 (newEntries * sizeof(dpair));
		(*_d) = d = realloc_or_die ("add_length_to_distribution",
									d, bytesMain + bytesHeap);

		d->items = (dpair*) (((char*) d) + bytesMain);
		d->size  = bytesHeap / sizeof(dpair);
		d->len   = len;
		}

	ix = d->len++;
	d->items[ix].length = length;
	d->items[ix].count  = 1;
	}


static u64 number_of_instances
   (distn*	d)
	{
	u32		ix;
	u64		count = 0;

	if (d == NULL) return 0;

	for (ix=0 ; ix<d->len ; ix++)
		count += d->items[ix].count;

	return count;
	}


static double average_length
   (distn*	d)
	{
	u32		ix;
	possum	sum   = 0;
	u64		count = 0;

	if (d == NULL) return -1;	// (since lengths are always strictly positive
								//  .. a negative average is impossible)

	for (ix=0 ; ix<d->len ; ix++)
		{
		count += d->items[ix].count;
		sum   += d->items[ix].count * d->items[ix].length;
		}

	if (count == 0) return -1;
	           else return ((double) sum) / count;
	}


static void print_length_distribution
   (FILE*	f,
	char*	prefix,
	distn*	d)
	{
	u32		ix;

	if ((d == NULL) || (d->len == 0))
		{ fprintf (f, "%s  (none)\n", prefix);  return; }

	qsort (d->items, d->len, sizeof(dpair), qCompareByLength);

	for (ix=0 ; ix<d->len ; ix++)
		fprintf (f, "%s  " unsposFmt ":" u64Fmt "\n",
		            prefix, d->items[ix].length, d->items[ix].count);
	}


static int qCompareByLength
   (const void*	_pairA,
	const void*	_pairB)
	{
	dpair*		pairA = (dpair*) _pairA;
	dpair*		pairB = (dpair*) _pairB;

	if      (pairA->length < pairB->length) return -1;
	else if (pairA->length > pairB->length) return  1;

	return 0;
	}

//==========
//
// The next four routines exist soley to support LASTZ's fmtInfStats output
// format.  This format allows collecting/reporting of stats upon which scoring
// inference can be performed by an external program.  Rather than binning the
// stats by pctid, they collect stats in a single bin (infStats instead of
// infStatsByPctId).
//
// This interface is supported only so that the early python version of INFERZ,
// which was used for the experiments discussed in [3], is still functional.
//
// The routines are
//	init_inference_stats_job
//	print_inference_stats_job
//	infer_stats_from_align_list
//	infer_stats_from_match
//
//==========

//----------
// 
// init_inference_stats_job--
//	Initialize inference stats.
//
//----------

void init_inference_stats_job
   (arg_dont_complain(seq* seq1),
	arg_dont_complain(seq* seq2))
	{
	if (statsActive)
		suicide ("attempt to open a second inference stats job");
	statsActive = true;

	init_stats (&infStats);
	}

//----------
//
// print_inference_stats_job--
//	Print inference stats.
//
//----------

static void private_print_inference_stats_job (FILE* f, infstats* inf);

void print_inference_stats_job
   (FILE*	f)
	{
	if (!statsActive)
		suicide ("attempt to close a non-existent inference stats job");

	private_print_inference_stats_job (f, &infStats);
	free_stats (&infStats);
	statsActive = false;
	}


static void private_print_inference_stats_job
   (FILE*		f,
	infstats*	inf)
	{
	char*		refSpecies = "seq1";
	char*		secSpecies = "seq2";

	if (!statsActive)
		suicide ("attempt to close a non-existent inference stats job");

	fprintf (f, "%s vs %s\n", refSpecies, secSpecies);
	fprintf (f, "  0%% < GC <= 100%%\n");

	fprintf (f, "    %-7s " unsposFmt " bases, " u64Fmt " gaps, " u64Fmt " runs\n",
	            refSpecies, inf->refBases,
		        number_of_instances (inf->refGaps),
		        number_of_instances (inf->refRuns));
	fprintf (f, "    %-7s " unsposFmt " bases, " u64Fmt " gaps, " u64Fmt " runs\n",
	            secSpecies, inf->secBases,
		        number_of_instances (inf->secGaps),
		        number_of_instances (inf->secRuns));

	print_bkgd_stats     (f, refSpecies, inf->refBkgd);
	print_bkgd_stats     (f, secSpecies, inf->secBkgd);
	print_subs_stats     (f,             inf->subs);
	print_blocks_stats   (f, refSpecies, inf->refBlocks);
	print_blocks_stats   (f, secSpecies, inf->secBlocks);
	print_gaps_stats     (f, refSpecies, inf->refGaps);
	print_gaps_stats     (f, secSpecies, inf->secGaps);
	print_runs_stats     (f, refSpecies, inf->refRuns);
	print_runs_stats     (f, secSpecies, inf->secRuns);
	print_segments_stats (f,             inf->segments);
	fprintf              (f, "\n");
	}

//----------
//
// infer_stats_from_align_list--
//	Collect inference stats from a list of gapped alignments.
//
//----------
//
// Arguments:
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//
// Returns:
//	(nothing)
//
//----------

void infer_stats_from_align_list
   (alignel*	alignList,
	seq*		seq1,
	seq*		seq2)
	{
	alignel*	a;

	for (a=alignList ; a!=NULL ; a=a->next)
		accumulate_stats_from_align (seq1, a->beg1-1, a->end1,
		                             seq2, a->beg2-1, a->end2,
		                             a->script, &infStats);
	}

//----------
//
// infer_stats_from_match--
//	Collect inference stats from a single ungapped alignment.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos2:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the match.
//
// Returns:
//	(nothing)
//
//----------

void infer_stats_from_match
   (seq*		seq1,
	unspos		pos1,
	seq*		seq2,
	unspos		pos2,
	unspos		length)
	{
	accumulate_stats_from_match (seq1, pos1, seq2, pos2, length, &infStats);
	}

//----------
//
// infer_scores_zero_stats--
//	Clear the statistics for this module.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

void infer_scores_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&inferScoresStats, 0, sizeof(inferScoresStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	inferScoresStats.averageGapLength     = 0.0;
	inferScoresStats.averageSegmentLength = 0.0;
	inferScoresStats.pExtend              = 0.0;
	inferScoresStats.sExtend              = 0.0;
	inferScoresStats.pOpen                = 0.0;
	inferScoresStats.sOpen                = 0.0;
	inferScoresStats.scaleBy              = 0.0;

#endif // collect_stats
	}

//----------
//
// infer_scores_show_stats_subs, infer_scores_show_stats_gaps,
// infer_scores_show_stats--
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print the stats to.
//
// Returns:
//	(nothing)
//
//----------

void infer_scores_show_stats_subs
   (arg_dont_complain(FILE* f),
	arg_dont_complain(int   trial))
	{
#ifdef collect_stats
	int		x, y;
	u32		totalSubs;

	if (f == NULL) return;

	fprintf (f, "(trial s%03d)\n", trial);
	infer_scores_show_stats_common (f);

	// inference counts

	totalSubs = 0;
	for (x=0 ; x<4 ; x++)
			for (y=0 ; y<4 ; y++)
		totalSubs += inferScoresStats.subs[x][y];

	fprintf (f, "%20s  %6s  ", "", "");
	for (y=0 ; y<4 ; y++)
		fprintf (f, " %6c", bits_to_nuc[y]);
	fprintf (f, "\n");

	for (x=0 ; x<4 ; x++)
		{
		if (x == 1) fprintf (f, "  inference counts: ");
		       else fprintf (f, "%20s", "");
		fprintf (f, "         %c", bits_to_nuc[x]);
		for (y=0 ; y<4 ; y++)
			fprintf (f, " %6u", inferScoresStats.subs[x][y]);
		fprintf (f, "\n");
		}

	fprintf (f, "           (total): %s\n", commatize(totalSubs));

	// inference observations

	fprintf (f, "%20s  %6s  ", "", "");
	for (y=0 ; y<4 ; y++)
		fprintf (f, " %6c", bits_to_nuc[y]);
	fprintf (f, "\n");

	fprintf (f, "%20s  %6s  ", "", "");
	for (y=0 ; y<4 ; y++)
		fprintf (f, " %6u", inferScoresStats.n2[y]);
	fprintf (f, "\n");

	fprintf (f, "%20s  %6s  ", "", "");
	for (y=0 ; y<4 ; y++)
		fprintf (f, " ------");
	fprintf (f, "\n");

	for (x=0 ; x<4 ; x++)
		{
		if (x == 0) fprintf (f, "      observations: ");
		       else fprintf (f, "%20s", "");
		fprintf (f, "%c", bits_to_nuc[x]);
		fprintf (f, " %6u |", inferScoresStats.n1[x]);
		for (y=0 ; y<4 ; y++)
			fprintf (f, " %6u", inferScoresStats.m[x][y]);
		fprintf (f, "\n");
		}

	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void infer_scores_show_stats_gaps
   (arg_dont_complain(FILE* f),
	arg_dont_complain(int   trial))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "(trial g%03d)\n", trial);
	infer_scores_show_stats_common (f);

	fprintf (f, "average gap length: %.13f\n", inferScoresStats.averageGapLength);
	fprintf (f, "average seg length: %.13f\n", inferScoresStats.averageSegmentLength);
	fprintf (f, "         p(extend): %.13f\n", inferScoresStats.pExtend);
	fprintf (f, "         s(extend): %.13f\n", inferScoresStats.sExtend);
	fprintf (f, "           p(open): %.13f\n", inferScoresStats.pOpen);
	fprintf (f, "           s(open): %.13f\n", inferScoresStats.sOpen);
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void infer_scores_show_stats_common
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "    total coverage: %s\n", commatize(inferScoresStats.coverageTotal));
	fprintf (f, "      low coverage: %s\n", commatize(inferScoresStats.coverageLow));
	fprintf (f, "     high coverage: %s\n", commatize(inferScoresStats.coverageHigh));
	fprintf (f, "      low identity: " identityBinLongFormat "\n", bin_top_to_identity   (inferScoresStats.lowIdentityBin));
	fprintf (f, "     high identity: " identityBinLongFormat "\n", bin_bottom_to_identity(inferScoresStats.highIdentityBin));
	fprintf (f, "          scale by: %.13f\n", inferScoresStats.scaleBy);
#endif // collect_stats
	}

void infer_scores_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "(no infer_scores stats)\n");
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void infer_scores_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "total_coverage: %" PRId64 "\n", inferScoresStats.coverageTotal);
	(*func) (f, "low_coverage: %" PRId64 "\n",   inferScoresStats.coverageLow);
	(*func) (f, "high_coverage: %" PRId64 "\n",  inferScoresStats.coverageHigh);
	(*func) (f, "low_identity: "  identityBinLongFormat "\n", bin_top_to_identity   (inferScoresStats.lowIdentityBin));
	(*func) (f, "high_identity: " identityBinLongFormat "\n", bin_bottom_to_identity(inferScoresStats.highIdentityBin));
#endif // collect_stats
	}

