//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: tweener.c
//
//----------
//
// tweener--
//	Given a list of alignments, search for alignments between them at a
//	higher-sensitivity ("weak" alignments).  We search for such "in-between
//	alignments" between properly ordered pairs of alignments and/or sequence
//	end points.
//
// We process the outer alignments in order and maintain list of alignments
// that have been processed but whose end point in seq1 is still potentially
// active (within window size of the start of the current alignment).
//
// For each alignment A we do the following:
//
//	(1)	Dismiss any now inactive alignments (alignments in the active list that
//		end too far before A starts).  If such an alignment is the end of a
//		chain of (1 or more) alignments, look for weak alignments to its
//		right.
//
//	(2)	Look for a current alignment, B, that overlaps A, i.e., B's end point is
//		within windowSize diagonals of where A starts, and follows A in one of
//		the sequences.  In the abnormal case that B ends after A ends (relative
//		to one of the sequences, mark A as "not the right end of a chain".
//
//	(3)	If no alignments overlap A, look for the alignment B that ends before A
//		(in both sequences) and is closest.  If it comes within windowSize bp in
//		both sequences, search for weak alignments between A and B.  If no such
//		B exists, think of A as being on the left end of a chain of (1 or more)
//		outer alignments, and search for weak alignments to its left.
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
#include "build_options.h"		// build options
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff
#include "seed_search.h"		// seed hit search stuff
#include "segment.h"			// segment table management stuff
#include "chain.h"				// segment chaining stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "edit_script.h"		// alignment edit script stuff

#define  tweener_owner			// (make this the owner of its globals)
#include "tweener.h"			// interface to this module

//----------
//
// private global data
//
//----------

int dbgShowInnerHsps = false;

#define verifyOrDie
//#define debugTweener

static alignel* innerList;

// set of alignments that have been checked, but may lie within windowSize bp
// of the start of a future alignment, relative to sequence #1

typedef struct active
	{
	alignel*	align;
	int			isRightEnd;	// right end of seen-so-far chain
	struct active* next;
	} active;

static active* activeList;
static active* activePool;

#undef  min
#define min(x,y) ((x)<(y)?(x):(y))
#undef  max
#define max(x,y) ((x)>(y)?(x):(y))

// private globals shared by all the routines under the umbrella of tweener()

static seq*		 seq1;
static seq*		 seq2;
static int       selfCompare;
static int       inhibitTrivial;
static const s8* upperCharToBits;
static seed*	 innerSeed;
static scoreset* scoring;
static scoreset* maskedScoring;
static tback*	 tb;
static score	 xDrop;
static int       gappedAllBounds;
static score	 yDrop;
static int		 trimToPeak;
static sthresh	 scoreThresh;
static score	 diagPen, antiPen;
static int		 scale;
static chainer	 connect;
static u32		 windowSize;

hitprocsimple	 simpleInfo;

#define numDefaultAnchors 100
static segtable* innerAnchors = NULL;

static seq* tweenSeq1;
static seq* tweenSeq2;
static seq* tweenRev1;
static seq* tweenRev2;

//----------
//
// prototypes for private functions
//
//----------

static void     try_bounded_align     (unspos b1, unspos e1,
                                       unspos b2, unspos e2);
static void     bounded_align         (unspos b1, unspos e1,
                                       unspos b2, unspos e2);
static u32      collect_inner_hsps    (void* info,
                                       unspos pos1, unspos pos2, unspos length,
                                       score s);
static alignel* merge_align           (alignel* a, alignel* b);
static void     activate              (alignel* a, int isRightEnd);
static active*  dismiss               (active* c);
static void     extract_subsequence   (seq* sf, unspos b, unspos e, seq* dst);
static void     copy_reverse_sequence (seq* sf, seq* dst);

//----------
//
// tweener_interpolate--
//	Interpolate "in-between" alignments in a chain of outer alignments.
//
//----------
//
// Arguments:
//	alignel*	alignList:	The list of 'outer' alignments.
//	 ...
//	u32			windowSize:	(see description in this file's header)
//
//	The other arguments are all pass-thru's to the following routines:
//		build_seed_position_table
//		seed_hit_search
//		reduce_to_chain
//		reduce_to_points
//		gapped_extend
//
// Returns:
//	A pointer to the updated alignment list.  This includes the outer
//	alignments (in the same memory as they were in in the input list) and any
//	inner alignments we find.
//
//----------

#define pairFmt  unsposFmt ".." unsposFmt
#define pairsFmt "(" unsposFmt "," unsposFmt ")..(" unsposFmt "," unsposFmt ")"

// macros for tweener_interpolate()

#ifndef verifyOrDie
#define verifyOrDie_1 ;
#define verifyOrDie_2 ;
#endif // not verifyOrDie

#ifdef verifyOrDie

#define verifyOrDie_1                                                        \
	for (a=alignList ; a!=NULL ; a=a->next)                                  \
		{                                                                    \
		if ((a->next != NULL) && (a->next->beg1 < a->beg1))                  \
			suicide ("outer alignments out of order");                       \
		}

#define verifyOrDie_2                                                        \
		for (c=activeList ; c!=NULL ; c=c->next)                             \
			{                                                                \
			b1 = c->align->end1;                                             \
			if (a1 > b1 + windowSize)                                        \
				suicidef ("tweener: impossible");                            \
			}

#endif // verifyOrDie


#ifndef debugTweener
#define debugTweener_1  ;
#define debugTweener_2  ;
#define debugTweener_3A ;
#define debugTweener_3B ;
#define debugTweener_3C ;
#endif // not debugTweener

#ifdef debugTweener

#define debugTweener_1                                                       \
	fprintf(stderr, "tweening these alignments:\n");                         \
	for (a=alignList ; a!=NULL ; a=a->next)                                  \
		fprintf(stderr, "  " pairsFmt "\n",                                  \
	                    a->beg1, a->beg2, a->end1, a->end2);                 \
	if (seq1->partition.p != NULL)                                           \
		print_partition_table (stderr, seq1);                                \
	if (seq2->partition.p != NULL)                                           \
		print_partition_table (stderr, seq2);                                \
	else                                                                     \
		fprintf(stderr, "seq2 = \"%s\"\n", seq2->header);

#define debugTweener_2                                                       \
	fprintf(stderr, "considering alignment " pairsFmt "\n",                  \
	                a1, a2, a->end1, a->end2);

#define debugTweener_3A                                                      \
	fprintf(stderr,"  overlap with " pairsFmt "?",                           \
	               b->beg1, b->beg2, b1, b2);

#define debugTweener_3B                                                      \
	fprintf(stderr, " (yes)\n");

#define debugTweener_3C                                                      \
	else fprintf(stderr, " (no)\n");

#endif // debugTweener

// tweener_interpolate--

alignel* tweener_interpolate
   (alignel*	alignList,
	seq*		_seq1,
	seq*		_seq2,
	int			_selfCompare,
	int			_inhibitTrivial,
	const s8	_upperCharToBits[],
	seed*		_innerSeed,
	scoreset*	_scoring,
	scoreset*	_maskedScoring,
	tback*		_tb,
	score		_xDrop,
	int			_gappedAllBounds,
	score		_yDrop,
	int			_trimToPeak,
	score		_scoreThresh,
	score		_diagPen,
	score		_antiPen,
	int			_scale,
	chainer		_connect,
	u32			_windowSize)
	{
	active*		c, *cNext;
	alignel*	a, *b;
	unspos		a1, a2, a1Lft;
	unspos		b1, b2;
	sgnpos		distToB, distToC;
	unspos		distD;
	int			isLeftEnd, isRightEnd, hasOverlap;

	// make sure we have something to do

	if (alignList == NULL)
		return NULL;

	//////////
	// setup
	//////////

	// copy parameters into globals

	seq1            = _seq1;
	seq2            = _seq2;
	selfCompare     = _selfCompare;
	inhibitTrivial  = _inhibitTrivial;
	upperCharToBits = _upperCharToBits;
	innerSeed       = _innerSeed;
	scoring         = _scoring;
	maskedScoring   = _maskedScoring;
	tb              = _tb;
	xDrop           = _xDrop;
	gappedAllBounds = _gappedAllBounds;
	yDrop           = _yDrop;
	trimToPeak      = _trimToPeak;
	scoreThresh.t   = 'S';
	scoreThresh.s   = _scoreThresh;
	diagPen         = _diagPen;
	antiPen         = _antiPen;
	scale           = _scale;
	connect         = _connect;
	windowSize      = round_up_2 (_windowSize); // (rounded up to an even number)

	// create sequence structures to use for each in-between block's
	// subsequences

	tweenSeq1 = new_sequence (windowSize);
	tweenSeq2 = new_sequence (windowSize);
	tweenRev1 = new_sequence (windowSize);
	tweenRev2 = new_sequence (windowSize);

	tweenSeq1->needTrueLen = tweenRev1->needTrueLen = seq1->needTrueLen;
	tweenSeq2->needTrueLen = tweenRev2->needTrueLen = seq2->needTrueLen;

	// set up part of the info record for process_for_simple_hit

	simpleInfo.hp.reporter         = collect_inner_hsps;
	simpleInfo.hp.reporterInfo     = NULL;
	simpleInfo.hp.minMatches       = -1;			// (no filtering)
	simpleInfo.hp.gfExtend         = gfexXDrop;
	simpleInfo.hp.seq1             = tweenSeq1;
	simpleInfo.hp.seq2             = tweenSeq2;
	simpleInfo.hp.scoring          = maskedScoring;
	simpleInfo.hp.xDrop            = xDrop;
	simpleInfo.hp.hspThreshold     = scoreThresh;
	simpleInfo.hp.hspZeroThreshold = (scoreThresh.t !='S')? 0
		                           : (scoreThresh.s >  0 )? scoreThresh.s
		                                                  : 0;
	simpleInfo.hp.anchors          = NULL;
	simpleInfo.hp.entropicHsp      = false;
	simpleInfo.hp.reportEntropy    = false;

	// create the table in which to collect the HSPs

	innerAnchors = new_segment_table (numDefaultAnchors, 0);

	// set up a pool of active elements (initially empty)

	activePool = NULL;

	//////////
	// process each viable in-between block (each inter-block gap, plus end
	// gaps)
	//////////

	activeList = NULL;
	innerList  = NULL;

	verifyOrDie_1
	debugTweener_1

	// process each outer alignment alignment in turn

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		a1 = a->beg1;
		a2 = a->beg2;
		a1Lft = (a1-1 < windowSize)? 0 : (a1 - windowSize);

		debugTweener_2

		// dismiss alignments that are too far from A

		while ((activeList != NULL) && (activeList->align->end1) < a1Lft)
			activeList = dismiss (activeList);

		c = activeList;
		while ((c != NULL) && ((cNext = c->next) != NULL))
			{
			if (a1Lft > cNext->align->end1) c->next = dismiss (cNext);
			                           else c = cNext;
			}

		verifyOrDie_2

		// look for an active alignment that overlaps A

		hasOverlap = false;
		for (c=activeList ; c!=NULL ; c=c->next)
			{
			b     = c->align;
			b1    = b->end1;
			b2    = b->end2;
			distD = abs((((sgnpos)b2) - ((sgnpos)b1)) 	// (distance between
                      - (((sgnpos)a2) - ((sgnpos)a1)));	//  .. diagonals) 

			debugTweener_3A

			if ((distD <= windowSize)
			 && ((b1 >= a1) || (b2 >= a2)))
				{
				hasOverlap = true;
				debugTweener_3B
				if ((b1 < a->end1) && (b2 < a->end2))
					c->isRightEnd = false;	// B ends properly -- before A ends
				else
					break;
				}
			debugTweener_3C
			}

		if (hasOverlap)
			{
			// if c is NULL, all overlaps were proper, so we're at the
			// right end of a chain

			isRightEnd = (c == NULL);
			activate (a, isRightEnd);
			continue; // don't try to tween on the left side of A
			}

		// find the closest alignment B to A such that B is active and ends
		// before the start of A in both sequences, but doesn't end more than 
		// windowSize bp before the start of A

		b = NULL;
		distToB = (sgnpos) (3*windowSize);
		isLeftEnd = true;		// A is the first outer alignment in a chain
		for (c=activeList ; c!=NULL ; c=c->next)
			{
			b1 = c->align->end1;
			b2 = c->align->end2;
			if ((b1 < a1) && (b2 < a2) && (a2 < b2 + windowSize))
				{
				isLeftEnd = false;
				if (c->isRightEnd)
					{
					distToC = ((sgnpos)a1) - ((sgnpos)b1)		// (manhattan
					        + ((sgnpos)a2) - ((sgnpos)b2);		//  .. distance)
					if (distToC < distToB)
						{ b = c->align;  distToB = distToC; }
					}
				c->isRightEnd = false;
				}
			}

		if (b != NULL)
			{
			b1 = b->end1;
			b2 = b->end2;
			try_bounded_align (b1, a1, b2, a2);
			}
		else if (isLeftEnd)
			{
			// A could be the first outer alignment in a chain;  look for
			// "in-between" alignments to the left

			b1 = (a1 <= windowSize/2)? 1 : (a1-windowSize/2);
			b2 = (a2 <= windowSize/2)? 1 : (a2-windowSize/2);
			try_bounded_align (b1, a1, b2, a2);
			}
		
		activate (a, true);
		}

	// align in windows after each chain-ending active alignment

	while (activeList != NULL)
		activeList = dismiss (activeList);

#ifdef collect_stats
	a1 = 0;
	for (a=innerList ; a!=NULL ; a=a->next)
		a1 += (a->end1 - a->beg1 + 1);
	tweener_add_stat (tweenerCoverage, a1);
#endif

	alignList = merge_align (alignList, innerList);

#ifdef collect_stats
	a1 = 0;
	for (a=alignList ; a!=NULL ; a=a->next)
		a1 += (a->end1 - a->beg1 + 1);
	tweener_add_stat (totalCoverage, a1);
#endif

	//////////
	// cleanup
	//////////

	free_sequence      (tweenSeq1);
	free_sequence      (tweenSeq2);
	free_sequence      (tweenRev1);
	free_sequence      (tweenRev2);
	free_segment_table (innerAnchors);

	for (c=activeList ; c!=NULL ; c=cNext) { cNext = c->next;  free(c); }
	for (c=activePool ; c!=NULL ; c=cNext) { cNext = c->next;  free(c); }

	return alignList;
	}

//----------
//
// try_bounded_align--
//	This routine is a wrapper for bounded_align(), to handle cases that arise
//	when either of sequence 1 or sequence 2 is partitioned.  In this case, we
//	must check whether the interval is split by a sequence boundary and break
//	apart the interval(s) passed to bounded_align().
//
//----------
//
// Arguments:
//	unspos	b1,e1:	Range of the rectangle in sequence 1 (origin 1, inclusive).
//	unspos	b2,e2:	Range of the rectangle in sequence 2.
//
// Returns:
//	(nothing)
//
//----------

// macros for bounded_align()

#ifndef debugTweener
#define debugTweener_E1 ;
#define debugTweener_E2 ;
#define debugTweener_E3 ;
#define debugTweener_E4 ;
#define debugTweener_E5 ;
#define debugTweener_E6 ;
#define debugTweener_E7 ;
#endif // not debugTweener

#ifdef debugTweener

#define debugTweener_E1                                                      \
	fprintf(stderr, "  try_bounded_align: " pairsFmt "\n",                   \
	                b1, b2, e1, e2);

#define debugTweener_E2                                                      \
	fprintf(stderr, "    seq1 partition1=" pairFmt " %s\n",                  \
				    part1->sepBefore, part1->sepAfter,                       \
	                &sp1->pool[part1->header]);                              \
	fprintf(stderr, "    seq1 partition2=" pairFmt " %s\n",                  \
				    part2->sepBefore, part2->sepAfter,                       \
	                &sp1->pool[part2->header]);                              \
	fprintf(stderr, "    part2-part1=%ld\n", part2-part1);

#define debugTweener_E3                                                      \
	fprintf(stderr, "    (E3) splitting " pairFmt                            \
	                " into " pairFmt "," pairFmt "\n",                       \
	                b1, e1, b1, e1left, b1right, e1);

#define debugTweener_E4                                                      \
	fprintf(stderr, "    seq2 partition1=" pairFmt " %s\n",                  \
				    part1->sepBefore, part1->sepAfter,                       \
	                &sp2->pool[part1->header]);                              \
	fprintf(stderr, "    seq2 partition2=" pairFmt " %s\n",                  \
				    part2->sepBefore, part2->sepAfter,                       \
	                &sp2->pool[part2->header]);                              \
	fprintf(stderr, "    part2-part1=%ld\n", part2-part1);

#define debugTweener_E5                                                      \
	fprintf(stderr, "    (E5) splitting " pairFmt                            \
	                " into " pairFmt "," pairFmt "\n",                       \
	                b2, e2, b2, e2left, b2right, e2);

#define debugTweener_E6                                                     \
	fprintf(stderr, "    seq1 partitionX=" pairFmt " %s\n",                  \
				    partX->sepBefore, partX->sepAfter,                       \
	                &sp1->pool[partX->header]);

#define debugTweener_E7                                                     \
	fprintf(stderr, "    seq2 partitionX=" pairFmt " %s\n",                  \
				    partY->sepBefore, partY->sepAfter,                       \
	                &sp2->pool[partY->header]);

#endif // debugTweener

// try_bounded_align--

static void try_bounded_align
   (unspos			b1,
	unspos			e1,
	unspos			b2,
	unspos			e2)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part1, *part2;
	partition*		partX1, *partX2, *partY1, *partY2, *partX, *partY;
	int				split1, split2;
	unspos			b1right, e1left, b2right, e2left;
	unspos			b1x, e1x, b2y, e2y;

	debugTweener_E1

	// if either interval is empty, let's not bother with it.

	if ((b1 == e1) || (b2 == e2)) return;

	// if neither sequence is partitioned, just pass the interval along to
	// bounded_align()

	if ((sp1->p == NULL)	// sequence 1 is not partitioned
	 && (sp2->p == NULL))	// sequence 2 is not partitioned
		{
		bounded_align (b1, e1, b2, e2);
		return;
		}

	// determine whether the interval is split in sequence 1

	split1  = false;
	e1left  = e1;
	b1right = b1;
	partX1  = partX2 = NULL;

	if (sp1->p != NULL)		// sequence 1 is partitioned
		{
		if      (seq1->v[b1-1] == 0) b1 += 1;
		else if (seq1->v[b1  ] == 0) b1 += 2;
		if      (seq1->v[e1-1] == 0) e1 -= 1;
		if (b1 >= e1) return;

		part1 = lookup_partition (seq1, b1-1);
		part2 = lookup_partition (seq1, e1-1);
		if (part1 != part2)
			{
			// split is b1 / part1->sepAfter / part2->sepBefore / e1
			//      --> b1 / e1left          /          b1right / e1
			debugTweener_E2
			e1left  = part1->sepAfter;
			b1right = part2->sepBefore+2;
			split1  = true;
			debugTweener_E3
			if (part2 - part1 > 1)
				{ partX1 = part1+1;  partX2 = part2-1; }
			}
		}

	// determine whether the interval is split in sequence 2

	split2  = false;
	e2left  = e2;
	b2right = b2;
	partY1  = partY2 = NULL;

	if (sp2->p != NULL)		// sequence 2 is partitioned
		{
		if      (seq2->v[b2-1] == 0) b2 += 1;
		else if (seq2->v[b2  ] == 0) b2 += 2;
		if      (seq2->v[e2-1] == 0) e2 -= 1;
		if (b2 >= e2) return;

		part1 = lookup_partition (seq2, b2-1);
		part2 = lookup_partition (seq2, e2-1);
		if (part1 != part2)
			{
			// split is b2 / part1->sepAfter / part2->sepBefore / e2
			//      --> b2 / e2left          /          b2right / e2
			debugTweener_E4
			e2left  = part1->sepAfter;
			b2right = part2->sepBefore+2;
			split2  = true;
			debugTweener_E5
			if (part2 - part1 > 1)
				{ partY1 = part1+1;  partY2 = part2-1; }
			}
		}

	// if neither sequence is split, just pass the interval along to
	// bounded_align()

	if ((!split1) && (!split2))
		{
		bounded_align (b1, e1, b2, e2);
		return;
		}

	// otherwise, send the split intervals to bounded_align(), plus intervals
	// for any intervening partitions

	bounded_align (b1, e1left, b2, e2left);
	bounded_align (b1right, e1, b2right, e2);

	if ((partX1 != NULL) && (partY1 == NULL))
		{
		// loop over intervening partitions in seq1
		for (partX=partX1 ; partX<=partX2 ; partX++)
			{
			debugTweener_E6
			b1x = partX->sepBefore+2;
			e1x = partX->sepAfter;
			bounded_align (b1x, e1x, b2, e2left);
			}
		}
	else if ((partX1 == NULL) && (partY1 != NULL))
		{
		// loop over intervening partitions in seq2
		for (partY=partY1 ; partY<=partY2 ; partY++)
			{
			debugTweener_E7
			b2y = partY->sepBefore+2;
			e2y = partY->sepAfter;
			bounded_align (b1, e1left, b2y, e2y);
			}
		}
	else if ((partX1 != NULL) && (partY1 != NULL))
		{
		// loop over intervening partitions in both sequences
		for (partX=partX1 ; partX<=partX2 ; partX++)
				for (partY=partY1 ; partY<=partY2 ; partY++)
			{
			b1x = partX->sepBefore+2;
			e1x = partX->sepAfter;
			b2y = partY->sepBefore+2;
			e2y = partY->sepAfter;
			bounded_align (b1x, e1x, b2y, e2y);
			}
		}
	}

//----------
//
// bounded_align--
//	Perform a high-sensitivity alignment within a specified rectangle.
//
//----------
//
// Arguments:
//	unspos	b1,e1:	Range of the rectangle in sequence 1 (origin 1,
//					.. inclusive).
//	unspos	b2,e2:	Range of the rectangle in sequence 2.
//
// Returns:
//	(nothing)
//
//----------

// macros for bounded_align()

#ifndef debugTweener
#define debugTweener_F1 ;
#define debugTweener_F2 ;
#define debugTweener_F3 ;
#define debugTweener_F4 ;
#endif // not debugTweener

#ifdef debugTweener

#define debugTweener_F1                                                      \
	fprintf(stderr, "  bounded_align: " pairsFmt "\n", b1, b2, e1, e2);

#define debugTweener_F2                                                      \
	{                                                                        \
	segtable* st = innerAnchors;                                             \
	segment* seg;                                                            \
	for (seg=st->seg ; (seg-st->seg)<st->len ; seg++)                        \
		fprintf(stderr, "    HSP " unsposSlashFmt " " unsposFmt " " scoreFmt "\n",  \
		                seg->pos1, seg->pos2, seg->length, seg->s);          \
	}

#define debugTweener_F3                                                      \
	{                                                                        \
	segtable* st = innerAnchors;                                             \
	segment* seg;                                                            \
	for (seg=st->seg ; (seg-st->seg)<st->len ; seg++)                        \
		fprintf(stderr, "    chained " unsposSlashFmt " " unsposFmt " " scoreFmt "\n", \
		                seg->pos1, seg->pos2, seg->length, seg->s);          \
	}

#define debugTweener_F4                                                      \
	for (aa=a ; aa!=NULL ; aa=aa->next)                                      \
		fprintf(stderr, "    new alignment " pairsFmt "\n",                  \
		                aa->beg1, aa->beg2, aa->end1, aa->end2);

#endif // debugTweener

// bounded_align--

static void bounded_align
   (unspos		b1,
	unspos		e1,
	unspos		b2,
	unspos		e2)
	{
	postable*	seq1Positions = NULL;
	alignel*	a = NULL;
	alignel*	aa;

	tweener_count_stat (numTweeners);
	tweener_add_stat   (totalArea, (e1-((u64)b1-1)) * (e2-((u64)b2-1)));

	debugTweener_F1

	// create the tiny subsequences

	extract_subsequence (seq1, b1-1, e1, tweenSeq1);
	extract_subsequence (seq2, b2-1, e2, tweenSeq2);

	// build word position table for the first sequence

	seq1Positions = build_seed_position_table
				       (tweenSeq1, 0, tweenSeq1->len, upperCharToBits,
				        innerSeed, /*step*/ 1);

	// find HSPs;  they get collected into innerAnchors[] by
	// collect_inner_hsps(), which is called from process_for_simple_hit()

	empty_segment_table (innerAnchors);
	seed_hit_search (tweenSeq1, seq1Positions,
	                 tweenSeq2, 0, tweenSeq2->len, /*selfCompare*/ false,
	                 upperCharToBits, innerSeed,
	                 /* searchLimit */ 0, 0,
#ifdef densityFiltering
	                 /*maxDensity*/ 0.0,
#endif // densityFiltering
	                 process_for_simple_hit, (void*) &simpleInfo);

	free_position_table (seq1Positions);

	debugTweener_F2

	// chain

	reduce_to_chain (innerAnchors, diagPen, antiPen, scale, connect);
	sort_segments   (innerAnchors, qSegmentsByPos1);

	debugTweener_F3

	// gapped extension

	if ((innerAnchors != NULL) && (innerAnchors->len != 0))
		{
		copy_reverse_sequence (tweenSeq1, tweenRev1);
		copy_reverse_sequence (tweenSeq2, tweenRev2);
		reduce_to_points (tweenSeq1, tweenSeq2, scoring, innerAnchors);
		a = gapped_extend (tweenSeq1, tweenRev1->v, tweenSeq2, tweenRev2->v,
		                   inhibitTrivial,
		                   scoring, innerAnchors, tb,
		                   gappedAllBounds, yDrop, trimToPeak, scoreThresh,
		                   /* no pairs limit */ 0, false, false);
		}

	// shift the positions, from subsequence back to sequence

	for (aa=a ; aa!=NULL ; aa=aa->next)
		{
		aa->seq1 =  seq1->v;
		aa->seq2 =  seq2->v;
		aa->beg1 += b1-1;
		aa->end1 += b1-1;
		aa->beg2 += b2-1;
		aa->end2 += b2-1;
		}

	debugTweener_F4

	innerList = merge_align (a, innerList);
	}

//----------
// [[-- a seed hit reporter function --]]
//
// collect_inner_hsps--
//	Collect a seed hit or HSP.
//
// Arguments and Return value: (see seed_search.h)
//
//----------

static u32 collect_inner_hsps
   (arg_dont_complain(void* info),
	unspos	pos1,
	unspos	pos2,
	unspos	length,
	score	s)
	{
	innerAnchors = add_segment (innerAnchors,
	                            pos1-length, pos2-length, length, s, /*id*/ 0);

	if (dbgShowInnerHsps)
		{
		fprintf (stderr, "\n");
		dump_aligned_nucleotides (stderr,
		                          tweenSeq1, pos1-length,
		                          tweenSeq2, pos2-length,
		                          length);
		}

	return 1;
	}

//----------
//
// merge_align--
//	Merge two beg1-ordered lists of alignments into a single beg1-ordered list.
//
//----------
//
// Arguments:
//	alignel*	a:	One input list.
//	alignel*	b:	The other.
//
// Returns:
//	A pointer to the merged list, which uses the same allocated memory as the
//	two input lists.
//
//----------

// macros for merge_align()

#ifndef verifyOrDie
#define verifyOrDie_3 ;
#endif // not verifyOrDie

#ifdef verifyOrDie

#define verifyOrDie_3                                                        \
	for (tail=a ; tail!=NULL ; tail=tail->next)                              \
		{                                                                    \
		if ((tail->next != NULL) && (tail->next->beg1 < tail->beg1))         \
			suicidef ("merge_align: first list out of order at " unsposFmt,  \
			          tail->next->beg1);                                     \
		}                                                                    \
	for (tail=b ; tail!=NULL ; tail=tail->next)                              \
		{                                                                    \
		if ((tail->next != NULL) && (tail->next->beg1 < tail->beg1))         \
			suicidef ("merge_align: second list out of order at " unsposFmt, \
			          tail->next->beg1);                                     \
		}

#endif // verifyOrDie

// merge_align--

static alignel* merge_align
   (alignel*	a,
	alignel*	b)
	{
	alignel* ret, *tail;

	verifyOrDie_3

	if (b == NULL) return a;
	if (a == NULL) return b;

	if (a->beg1 <= b->beg1) { ret = tail = a;  a = a->next; }
	                   else { ret = tail = b;  b = b->next; }

	while (a != NULL && b != NULL)
		{
		if (a->beg1 <= b->beg1) { tail = tail->next = a;  a = a->next; }
		                   else { tail = tail->next = b;  b = b->next; }
		}

	if (a == NULL) tail->next = b;
	          else tail->next = a;

	return ret;
	}

//----------
//
// activate--
//	Put an alignment at the front of the active list.
//
//----------
//
// Arguments:
//	alignel*	a:			The alignment.
//	int			isRightEnd:	true  => the alignment is the right end of a chain.
//							false => it isn't.
//
// Returns:
//	(nothing)
//
//----------

static void activate
   (alignel*	a,
	int			isRightEnd)
	{
	active*		c;

	// reclaim an element from the pool, otherwise allocate a new one

	if (activePool == NULL)   c = malloc_or_die ("activate", sizeof(*c));
	                   else { c = activePool;  activePool = c->next; }

	// create the element and add it to the head of the active list

	c->align      = a;
	c->isRightEnd = isRightEnd;
	c->next       = activeList;
	activeList    = c;
	}

//----------
//
// dismiss--
//	Remove an alignment from the active list.  If it is the outer alignment at
//	the right end of a chain, we look for "in-between" alignments beyond the
//	end.
//
//----------
//
// Arguments:
//	active*	c:		The alignment to dismiss.
//
// Returns:
//	c->next
//
//----------

static active* dismiss
   (active*	c)
   	{
   	active*	next;
	unspos	a1, a2, b1, b2;

	next = c->next;

	if (c->isRightEnd)
		{
		b1 = c->align->end1;
		b2 = c->align->end2;
		a1 = min (b1+windowSize/2, seq1->len);
		a2 = min (b2+windowSize/2, seq2->len);
		try_bounded_align (b1, a1, b2, a2);
		}

	// return c to the pool

	c->next = activePool;
	activePool = c;

	return next;
	}

//----------
//
// extract_subsequence--
//	Make a copy of a subsequence.
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to make a copy of.
//	unspos	b,e:	The index range of the subsequence to copy, origin-zero,
//					.. open end.
//	seq*	dst:	The sequence in which to place the copy.  Any previous
//					.. contents will be destroyed.
//
// Returns:
//	(nothing)
//
//----------

static void extract_subsequence
   (seq*	sf,
	unspos	b,
	unspos	e,
	seq*	dst)
	{
	u8*		s = sf->v;
	unspos	len, i;

	if (e <= b)
		suicidef ("internal error in extract_subsequence\n"
		          "  interval " unsposFmt ".." unsposFmt " is empty", b, e);

	len = e - b;
	sequence_long_enough (dst, len, false);

	for (i=0 ; i<len ; i++)
		dst->v[i] = s[b+i];
	dst->v[len] = 0;
	dst->len = len;

	dst->fileType = seq_type_nofile;
	dst->contig   = 1;
	dst->startLoc = 1;
	}

//----------
//
// copy_reverse_sequence--
//	Make a copy of a sequence, in reverse (*not* reverse complement).
//
//----------
//
// Arguments:
//	seq*	seq:	The sequence to make a copy of.
//	seq*	dst:	The sequence in which to place the copy.  Any previous
//					.. contents will be destroyed.
//
// Returns:
//	(nothing)
//
//----------

static void copy_reverse_sequence
   (seq*	sf,
	seq*	dst)
	{
	u32		len = sf->len;
	u8*		s, *d;

	sequence_long_enough (dst, len, false);

	for (s=sf->v+len,d=dst->v ; s>sf->v ; )
		*(d++) = *(--s);

	*d = 0;
	dst->len = len;
	}

//----------
//
// tweener_zero_stats--
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

void tweener_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&tweenerStats, 0, sizeof(tweenerStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// tweener_show_stats--
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

void tweener_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "number of tweeners: %s\n", commatize(tweenerStats.numTweeners));
	fprintf (f, "tweener total area: %s\n", commatize(tweenerStats.totalArea));
	if (tweenerStats.numTweeners > 0)
		fprintf (f, "  avg tweener area: %s\n", commatize((tweenerStats.totalArea / ((float) tweenerStats.numTweeners)) + .5));
	fprintf (f, "    total coverage: %s bp\n", commatize(tweenerStats.totalCoverage));
	fprintf (f, "  tweener coverage: %s bp\n", commatize(tweenerStats.tweenerCoverage));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void tweener_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "num_tweeners=%d\n",        tweenerStats.numTweeners);
	(*func) (f, "total_area=%" PRId64 "\n", tweenerStats.totalArea);
	(*func) (f, "total_coverage=%d\n",      tweenerStats.totalCoverage);
	(*func) (f, "tweener_coverage=%d\n",    tweenerStats.tweenerCoverage);
#endif // collect_stats
	}

