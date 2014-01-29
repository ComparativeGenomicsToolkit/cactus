//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: seed_search.c
//
//----------
//
// seed_search--
//	Support for finding "seed hits" in genomic sequences.
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
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed matching stuff
#include "pos_table.h"			// position table stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "segment.h"			// segment table management stuff

#define  seed_search_owner		// (make this the owner of its globals)
#include "seed_search.h"		// interface to this module

// debugging defines

//#define debugDiag 97-60724168	// if defined, breakdown what happens with
								// .. every seed hit on this pos1-pos2 diagonal

//#define debugResolvingSeeds	// if defined, show how partial seed hits are
//#define debugPos1 49615		// .. resolved for this seed hit (starting
//#define debugPos2 1345896		// .. positions of the hit)

#ifdef debugDiag				// if we are debugging a particular diagonal,
static int debugThisDiag;		// .. debugThisDiag will be set true only on
#else							// .. that diagonal;  otherwise it is set false
#define debugThisDiag false		// .. either at run-time or compile-time
#endif

#ifdef debugResolvingSeeds		// if we are debugging a particular seed hit,
static int debugThisHit;		// .. debugThisHit will be set true only on
#else							// .. that hit;  otherwise it is set false
#define debugThisHit false		// .. either at run-time or compile-time
#endif

//#define debugSearchPos2 2412	// if defined, breakdown what happens with this
								// .. position in sequence 2 (this is the right
								// .. end of the word;  the left end is at
								// .. X-(seedLength-1) counting from 1 at the
								// .. start of the sequence

// other debugging defines

//#define debugMismatchExtend	// if defined, show process of mismatch_extend_seed_hit

//#define snoopDiagHash			// if this is defined, extra code is added to
								// .. report on diag hash collisions

//#define snoopXDrop			// if this is defined, extra code is added to
								// .. report on x-drop extensions

//#define snoopEntropy			// if this is defined, extra code is added to
								// .. report entropy calculations

//#define snoopHspSubrange		// if this is defined, extra code is added to
								// .. report situations where an HSP subrange
								// .. scores higher than the HSP
//#define extendHspFromLeft		// if this is defined, HSP extension is from
								// .. the left edge of the hit, not the right

//#define densityFiltering		// if this is defined, we perform alignment
								// .. density filtering;  note that, as of this
								// .. writing, it has not been shown that
								// .. density filtering provides any benefit

//#define snoopPosFilter		// if this is defined, extra code is added to
								// .. filter_seed_hit_by_pos to report whether
								// .. or not each seed was discarded

//#define snoopBelowDiagonal	// if this is defined, extra code is added to
								// .. seed_hit_below_diagonal, for debugging

//#define snoopReporterCalls	// if this is defined, extra code is added to
								// .. to track calls to the seed hit reporter
								// .. function

#ifdef densityFiltering
#define densityCheckDepth2		// if this is defined, we check for density
								// .. filtering at depth 2

#define densityCheckDepth3		// if this is defined, we check for density
								// .. filtering at depth 3
#endif // densityFiltering

//----------
//
// stats to augment crude profiling
//
//----------

#ifndef dbgTiming
#define dbg_timing_set_stat(field,val)         ;
#define dbg_timing_count_stat(field)           ;
#define dbg_timing_report_stat(field,name)     ;
#define dbg_timing_report_big_stat(field,name) ;
#endif // not dbgTiming

#ifdef dbgTiming
struct
	{
	int64 ungappedExtensions;
	int   hsps;
	} seedSearchTimingStats;

#define dbg_timing_set_stat(field,val)         (seedSearchTimingStats.field = val)
#define dbg_timing_count_stat(field)           ++seedSearchTimingStats.field
#define dbg_timing_report_stat(field,name)     fprintf(stderr,"%-26s %d\n",         name":",seedSearchTimingStats.field)
#define dbg_timing_report_big_stat(field,name) fprintf(stderr,"%-26s %" PRId64 "\n",name":",seedSearchTimingStats.field)
#endif // dbgTiming

//----------
//
// private data
//
//----------

// private globals shared by all the routines under the umbrella of
// seed_hit_search()

static seq*			seq1;
static postable*	pt;
static seq*			seq2;
static unspos		start;
static unspos		end;
static int			selfCompare;
static int			sameStrand;  // (only meaningful if selfCompare is true)
static const s8*	upperCharToBits;
static seed*		hitSeed;
static u32			searchLimit;
static u32			reportSearchLimit;
static s32			searchToGo;
#ifdef densityFiltering
static u64			maxBasesAllowed;
#endif // densityFiltering
static hitprocessor	processor;
static void*		processorInfo;

static int          unblockedLeftExtension;	// true => seed extension routines
								// .. (e.g. xdrop_extend_seed_hit) shouldn't
								// .. block left-extension at previous diagEnd

// static data area for discovery_probability

static u32 foldedSize = 0;
static u8* foldedHits = NULL;

//----------
//
// prototypes for private functions
//
//----------

static u64   private_hit_search            (void);
static u64   private_hit_search_halfweight (void);
static u64   private_hit_search_resolve    (void);

static u64   find_table_matches         (u32 packed2, unspos pos2);
static u64   find_table_matches_resolve (u32 packed2, unspos pos2,
                                         u32 unpacked2, int transAllowed);

static int   seed_hit_below_diagonal (unspos pos1, unspos pos2);
static int   filter_seed_hit_by_pos  (hitprocinfo* hp,
                                      unspos pos1, unspos pos2, unspos length);
static int   filter_seed_hit_by_subs (hitprocinfo* hp,
                                      unspos pos1, unspos pos2, unspos length);
static score xdrop_extend_seed_hit   (hitprocinfo* hp,
                                      unspos* pos1, unspos* pos2, unspos* length);
static score match_extend_seed_hit   (hitprocinfo* hp,
                                      unspos* pos1, unspos* pos2, unspos* length);
static score mismatch_extend_seed_hit(hitprocinfo* hp,
                                      unspos* pos1, unspos* pos2, unspos* length);
static void  warn_for_search_limit   (void);
static void  dump_raw_hit            (FILE* f, unspos pos1, unspos pos2);

#ifdef debugDiag

static void  dump_extended_match     (FILE* f,
                                      seq* seq1, seq* seq2, sgnpos diag,
                                      u8* p1, u8* p2, u8* p3,
                                      u8* p4, u8* p5, u8* p6);
static char* pair_diagonal_as_text   (unspos pos1, unspos pos2);
static char* diagonal_as_text        (sgnpos diag);

#endif // debugDiag

#if ((defined snoopDiagHash) && (!defined debugDiag))
static char* pair_diagonal_as_text   (unspos pos1, unspos pos2);
static char* diagonal_as_text        (sgnpos diag);
#endif // snoopDiagHash && not debugDiag

#if ((defined snoopXDrop) && (!defined debugDiag))
static char* pair_diagonal_as_text   (unspos pos1, unspos pos2);
static char* diagonal_as_text        (sgnpos diag);
#endif // snoopXDrop && not debugDiag

#if (defined snoopXDrop)
static char* display_sequence_character (seq* _seq, u8 ch);
#endif // snoopXDrop

//----------
//
// seed_hit_search--
//	Search for seed hits between one sequence and another.
//
// The caller must already have built a table of seed-word positions in one of
// the sequences.
//
//----------
//
// Arguments:
//	seq*		seq1:			The sequence being searched.
//	postable*	pt:				A table of positions of words in seq1.
//	seq*		seq2:			The sequence being searched for.
//	unspos		start:			First sequence position to consider.  Zero is
//								.. the first possible position.
//	unspos		end:			One past the last sequence position to consider.
//								.. If this is zero, the sequence length is used.
//	int			selfCompare:	true => seq1 and seq2 are the same sequence.
//	s8			upperCharToBits[]: Table to map sequence characters to two-bit
//								.. values,and illegal characters to -1.
//	seed*		hitSeed:		The seed-word the table is based on.
//	u32			searchLimit:	The maximum number of "HSPs" allowed;  zero
//								.. indicates "no limit".  See note (3) below.
//	u32			reportSearchLimit:	The number to report (to the user) as the
//								.. search limit if searchLimit is reached.  Note
//								.. that searchLimit is per-search, whereas
//								.. reportSearchLimit is per-query.  The special
//								.. value of zero indicates that we should not
//								.. report this condition.
//	double		maxDensity:		(only if densityFiltering is #defined)
//								The maximum alignment density we will "allow"
//								.. before discarding a query sequence;  zero
//	                            .. means there is no limit (see note 2 below).
//	hitprocessor processor:		Function to call for each hit to determine if it
//								.. is 'good enough'.
//	void*		processorInfo:	A value to pass thru with each call to processor.
//
// Returns:
//	The number of bases in the seed hits (see note 2 below).
//
//----------
//
// Notes:
//
// (1)	This routine allocates and reuses memory via global pointers.  The
//		caller should make a call to free_seed_hit_search() to de-allocate this
//		memory, after all searches are complete.
//
// (2)	If the density limit is exceeded, the value u64max is returned.  It is
//		possible (in fact, a near certainty) that the processor (and associated
//		reporter) function will have been called.  It is up the the caller to
//		dispose of any seed hits already reported.
//
// (3)	searchLimit is not a hard limit, and there are many circumstances by
//		which we will report more seeds/hits/HSPs/alignments than the limit.
//		This allows us to avoid checking the limit after each and every hit.
//		The limit should be taken as giving us the permission to stop once we
//		have found that many hits.
//
//----------

#ifndef debugSearchPos2
#define debugSearchPos2_1 ;
#define debugSearchPos2_2 ;
#define debugSearchPos2_3 ;
#define debugSearchPos2_4 ;
#endif // not debugSearchPos2

#ifdef debugSearchPos2

#define debugSearchPos2_1                                                     \
	if (pos2 == debugSearchPos2)                                              \
		printf ("checking %s at seq 2 pos " unsposFmt " (matches)\n",         \
				seed_packed_to_string (hitSeed, packed), pos2);

#define debugSearchPos2_2                                                     \
	if (pos2 == debugSearchPos2)                                              \
		printf ("checking %s at seq 2 pos " unsposFmt " (one transition)\n",  \
				seed_packed_to_string (hitSeed, packedTrans), pos2);

#define debugSearchPos2_3                                                     \
	if (pos2 == debugSearchPos2)                                              \
		printf ("checking %s at seq 2 pos " unsposFmt " (one transition)\n",  \
				seed_packed_to_string (hitSeed, packedTrans), pos2);

#define debugSearchPos2_4                                                     \
	if (pos2 == debugSearchPos2)
		printf ("checking %s at seq 2 pos " unsposFmt " (two transitions)\n",
				seed_packed_to_string (hitSeed, packedTrans), pos2);

#endif // debugSearchPos2


u64 seed_hit_search
   (seq*			_seq1,
	postable*		_pt,
	seq*			_seq2,
	unspos			_start,
	unspos			_end,
	int				_selfCompare,
	const s8		_upperCharToBits[],
	seed*			_hitSeed,
	u32				_searchLimit,
	u32				_reportSearchLimit,
#ifdef densityFiltering
	double			_maxDensity,
#endif // densityFiltering
	hitprocessor	_processor,
	void*			_processorInfo)
	{
	u64				basesHit;
	seqpartition*	sp2;
	char*			name2;
	char			strand2;

	// sanity check

	if (_end == 0)
		_end = _seq2->len;

	if (_end <= _start)
		suicidef ("in seed_hit_search(), interval is void (%d-%d)",
		          _start, _end);

	if (_end > _seq2->len)
		suicidef ("in seed_hit_search(), interval end is bad (%d>%d)",
		          _end, _seq2->len);

	// allocate (or re-use) memory

	empty_diag_hash ();

	// pass globals to the rest of the search
	// note: this makes this module non-threadsafe

	seq1              = _seq1;
	pt                = _pt;
	seq2              = _seq2;
	start             = _start;
	end               = _end;
	selfCompare       = _selfCompare;
	sameStrand        = (selfCompare) && (seq1->revCompFlags == seq2->revCompFlags);
	upperCharToBits   = _upperCharToBits;
	hitSeed           = _hitSeed;
	searchLimit       = _searchLimit;
	reportSearchLimit = _reportSearchLimit;
	searchToGo        = _searchLimit;
#ifdef densityFiltering
	maxBasesAllowed   = _maxDensity * seq2->len;
#endif // densityFiltering
	processor	      = _processor;
	processorInfo     = _processorInfo;

	seed_search_set_stat (withTrans, hitSeed->withTrans);
	seed_search_set_stat (searchLimit, searchLimit);

	// perform the search

	if      (hitSeed->isHalfweight) basesHit = private_hit_search_halfweight ();
	else if (pt->asBits != NULL)    basesHit = private_hit_search_resolve    ();
	else                            basesHit = private_hit_search            ();

	// cleanup

	if (foldedHits != NULL)
		{ free_if_valid ("folded hits", foldedHits);  foldedHits = NULL; }

	if ((seed_search_dbgShowCoverage) && (basesHit > 0))
		{
		sp2 = &seq2->partition;
		if (sp2->p == NULL)		// sequence 2 is not partitioned
			{
			name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
			if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
			}
		else  					// sequence 2 is partitioned
			name2 = "seq2";

		strand2 = ((seq2->revCompFlags & rcf_rev) == 0)? '+' : '-';

		printf ("# seed bases hit in %s%c: " u64Fmt, name2, strand2, basesHit);
#ifdef densityFiltering
		if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
			printf (" (rejected)");
#endif // densityFiltering
		printf ("\n");
		}

#ifdef densityFiltering
	if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
		{
		if (seed_search_dbgShowRejections)
			{
			sp2 = &seq2->partition;
			if (sp2->p == NULL)		// sequence 2 is not partitioned
				{
				name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
				if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
				}
			else  					// sequence 2 is partitioned
				name2 = "seq2";

			strand2 = ((seq2->revCompFlags & rcf_rev) == 0)? '+' : '-';

			fprintf (stderr, "%s%c rejected due to excessive hsp density\n",
			                 name2, strand2);
			}

		return u64max;
		}
#endif // densityFiltering

	return basesHit;
	}


void free_seed_hit_search (void) { free_diag_hash (); }


// private_hit_search-- seed requires two bits per unpacked bp
//
// nota bene: I had expected that telling the compiler that the sequence is
//            not changing (by declaring vars as follows) would produce faster
//            code, but it was actually slightly slower.
//
//				const u8* const	qStart = seq2->v + start;
//				const u8* const	qStop  = seq2->v + end;
//				const u8* 		q;

static u64 private_hit_search (void)
	{
	int		seedLength;
	u8* 	qStart = seq2->v + start;
	u8* 	qStop  = seq2->v + end;
	u8* 	q;
	u64		w;
	s32		ww;
	u32		packed, packedTrans;
	int		nts;
	unspos	pos2;
	u32*	f1, *f2;
	u64		basesHit = 0;
#if ((defined collect_stats) && (defined maxHitsPerColumn))
	u64		prevRawHits, hitsInColumn;
#endif // collect_stats && maxHitsPerColumn

	seedLength = hitSeed->length;

	if (seedLength < 2)
		suicidef ("seed length must be at least two (yours is %d)", seedLength);

	if (seq2->len < (unsigned) seedLength)
		return 0; // (nothing to search for)

	// scan the sequence, processing each seed match

	for (q=qStart ; q<qStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(q<qStop) ; nts++)
			{
			ww = upperCharToBits[*(q++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; q<qStop ; q++)
			{
			ww = upperCharToBits[*q];				// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt

			pos2 = q-seq2->v + 1;
			packed = apply_seed (hitSeed, w);		// extract seed bits
			seed_search_count_stat (wordsInSequence);

			// generate seed hits for the complete seed match

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			prevRawHits = seedSearchStats.rawSeedHits;
#endif // collect_stats && maxHitsPerColumn
			debugSearchPos2_1;
			basesHit += find_table_matches (packed, pos2);

			// generate seed hits for all seed matches with 1 or 2 transitions

			if (hitSeed->withTrans == 1)
				{
				for (f1=hitSeed->transFlips ; *f1!=0 ; f1++)
					{
					packedTrans = packed ^ (*f1);
					debugSearchPos2_2;
					basesHit += find_table_matches (packedTrans, pos2);
					}
				}
			else if (hitSeed->withTrans >= 2)
				{
				for (f1=hitSeed->transFlips ; *f1!=0 ; f1++)
					{
					packedTrans = packed ^ (*f1);
					debugSearchPos2_3;
					basesHit += find_table_matches (packedTrans, pos2);
					for (f2=f1+1 ; *f2!=0 ; f2++)
						{
						packedTrans = packed ^ (*f1) ^ (*f2);
						debugSearchPos2_4;
						basesHit += find_table_matches (packedTrans, pos2);
						}
					}
				}

			if ((searchLimit > 0) && (searchToGo < 0))
				{ warn_for_search_limit ();  return basesHit; }

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			hitsInColumn = seedSearchStats.rawSeedHits - prevRawHits;
			if (hitsInColumn <= maxHitsPerColumn)
				seedSearchStats.hitsPerColumn[hitsInColumn]++;
			else
				{
				seedSearchStats.hitsPerColumn[maxHitsPerColumn+1]++;
				if (hitsInColumn > seedSearchStats.mostHitsInColumn)
					seedSearchStats.mostHitsInColumn = hitsInColumn;
				}
#endif // collect_stats && maxHitsPerColumn

#ifdef densityCheckDepth2
			if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
				return basesHit;
#endif // densityCheckDepth2
			}
		}

	return basesHit;
	}


// private_hit_search_halfweight-- seed requires one bit per unpacked bp

static u64 private_hit_search_halfweight (void)
	{
	int		seedLength;
	u8*		qStart = seq2->v + start;
	u8*		qStop  = seq2->v + end;
	u8*		q;
	u64		w;
	s32		ww;
	u32		packed;
	int		nts;
	unspos	pos2;
	u64		basesHit = 0;
#if ((defined collect_stats) && (defined maxHitsPerColumn))
	u64		prevRawHits, hitsInColumn;
#endif // collect_stats && maxHitsPerColumn

	seedLength = hitSeed->length;

	if (seedLength < 2)
		suicidef ("seed length must be at least two (yours is %d)", seedLength);

	if (seq2->len < (unsigned) seedLength)
		return 0; // (nothing to search for)

	// scan the sequence, processing each seed match

	for (q=qStart ; q<qStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(q<qStop) ; nts++)
			{
			ww = upperCharToBits[*(q++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y
			}

		// process each word of seedLength nucleotides

		for ( ; q<qStop ; q++)
			{
			ww = upperCharToBits[*q];				// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y

			pos2 = q-seq2->v + 1;
			packed = apply_seed (hitSeed, w);		// extract seed bits
			seed_search_count_stat (wordsInSequence);

			// generate seed hits for seed match

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			prevRawHits = seedSearchStats.rawSeedHits;
#endif // collect_stats && maxHitsPerColumn
			basesHit += find_table_matches (packed, pos2);
			if ((searchLimit > 0) && (searchToGo < 0))
				{ warn_for_search_limit ();  return basesHit; }

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			hitsInColumn = seedSearchStats.rawSeedHits - prevRawHits;
			if (hitsInColumn <= maxHitsPerColumn)
				seedSearchStats.hitsPerColumn[hitsInColumn]++;
			else
				{
				seedSearchStats.hitsPerColumn[maxHitsPerColumn+1]++;
				if (hitsInColumn > seedSearchStats.mostHitsInColumn)
					seedSearchStats.mostHitsInColumn = hitsInColumn;
				}
#endif // collect_stats && maxHitsPerColumn

#ifdef densityCheckDepth2
			if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
				return basesHit;
#endif // densityCheckDepth2
			}
		}

	return basesHit;
	}


// private_hit_search_resolve-- full seed requires two bits per unpacked bp,
//                              .. but not all seed bits are in the table, so
//                              .. the remaining bits must be resolved by direct
//                              .. comparison to sequence1

static u64 private_hit_search_resolve (void)
	{
	int		seedLength;
	int		transAllowed = hitSeed->withTrans;
	u8*		qStart = seq2->v + start;
	u8*		qStop  = seq2->v + end;
	u8*		q;
	u64		w;
	s32		ww;
	u32		packed, packedTrans;
	int		nts;
	unspos	pos2;
	u32*	f1, *f2;
	u64		basesHit = 0;
#if ((defined collect_stats) && (defined maxHitsPerColumn))
	u64		prevRawHits, hitsInColumn;
#endif // collect_stats && maxHitsPerColumn

	seedLength = hitSeed->length;

	if (seedLength < 2)
		suicidef ("seed length must be at least two (yours is %d)", seedLength);

	if (seq2->len < (unsigned) seedLength)
		return 0; // (nothing to search for)

	// scan the sequence, processing each seed match

	for (q=qStart ; q<qStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(q<qStop) ; nts++)
			{
			ww = upperCharToBits[*(q++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; q<qStop ; q++)
			{
			ww = upperCharToBits[*q];				// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt

			pos2 = q-seq2->v + 1;
			packed = apply_seed (hitSeed, w);		// extract seed bits
			seed_search_count_stat (wordsInSequence);

			// generate seed hits for the complete seed match

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			prevRawHits = seedSearchStats.rawSeedHits;
#endif // collect_stats && maxHitsPerColumn
			basesHit += find_table_matches_resolve (packed, pos2,
			                                        w, transAllowed);
			// generate seed hits for all seed matches with 1 or 2 transitions

			if (transAllowed == 1)
				{
				for (f1=hitSeed->transFlips ; *f1!=0 ; f1++)
					{
					packedTrans = packed ^ (*f1);
					basesHit += find_table_matches_resolve (packedTrans, pos2,
					                                        w, 0);
					}
				}
			else if (transAllowed >= 2)
				{
				for (f1=hitSeed->transFlips ; *f1!=0 ; f1++)
					{
					packedTrans = packed ^ (*f1);
					basesHit += find_table_matches_resolve (packedTrans, pos2,
					                                        w, 1);
					for (f2=f1+1 ; *f2!=0 ; f2++)
						{
						packedTrans = packed ^ (*f1) ^ (*f2);
						basesHit += find_table_matches_resolve (packedTrans, pos2,
						                                        w, 0);
						}
					}
				}

			if ((searchLimit > 0) && (searchToGo < 0))
				{ warn_for_search_limit ();  return basesHit; }

#if ((defined collect_stats) && (defined maxHitsPerColumn))
			hitsInColumn = seedSearchStats.rawSeedHits - prevRawHits;
			if (hitsInColumn <= maxHitsPerColumn)
				seedSearchStats.hitsPerColumn[hitsInColumn]++;
			else
				{
				seedSearchStats.hitsPerColumn[maxHitsPerColumn+1]++;
				if (hitsInColumn > seedSearchStats.mostHitsInColumn)
					seedSearchStats.mostHitsInColumn = hitsInColumn;
				}
#endif // collect_stats && maxHitsPerColumn

#ifdef densityCheckDepth2
			if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
				return basesHit;
#endif // densityCheckDepth2
			}
		}

	return basesHit;
	}

//----------
//
// find_table_matches, find_table_matches_resolve--
//	Given a packed word in sequence 2, find and process all its matches in the
//	table of sequence 1 word positions.
//
//	find_table_matches_resolve is used when the full seed is too big to use as
//	an index (into the position table), and we have to resolve any remaining
//	seed bits by comparing the sequences.
//
//----------
//
// Arguments:
//	u32		packed2:		The packed word, representing some seed/window of
//							.. nucleotides in sequence 2.
//	unspos	pos2:			The hit position in sequence 2.  This is the
//							.. position following the end of the hit.
//	u32		unpacked2:		(find_table_matches_resolve only)  The last 16
//							.. nucleotides from sequence 2, packed two bits per
//							.. nucleotide (oldest nt in most significant bits).
//	int		transAllowed:	(find_table_matches_resolve only)  The maximum
//							.. number of transitions allowed in the resolved
//							.. bits;  this is the same as the number of mis-
//							.. matches, because resolved bases can only be
//							.. matches or transitions (never transversions)
//
// Returns:
//	The number of bases in the seed hits.
//
//----------

static u64 find_table_matches
   (u32		packed2,
	unspos	pos2)
	{
	u32		seedLength, len1;
	unspos	adjStart = pt->adjStart;
	u32		step     = pt->step;
	unspos	pos, pos1;
	u64		basesHit = 0;

	seedLength = (unsigned) hitSeed->length;
	len1       = seedLength-1;

	if (pt->last[packed2] == 0)
		{
#ifdef debugSearchPos2
		if (pos2 == debugSearchPos2)
			printf (" no hits in sequence 1\n");
#endif
		return 0;
		}

	for (pos=pt->last[packed2] ; pos!=noPreviousPos ; pos=pt->prev[pos])
		{
		pos1 = adjStart + step*pos;

#ifdef debugSearchPos2
		if (pos2 == debugSearchPos2)
			printf ("  hit at pos " unsposFmt "\n", pos1);
#endif

		if ((selfCompare) && (seed_hit_below_diagonal (pos1, pos2)))
			continue;

		if (seed_search_dbgDumpRawHits)
			{
			if (seed_search_dbgShowRawHits)
				dump_raw_hit (stderr, pos1, pos2);
			else
				{
				printf ("\nraw seed hit " unsposSlashFmt "\n", pos1-len1, pos2-len1);
				dump_aligned_nucleotides (stdout,
				                          seq1, pos1-seedLength,
				                          seq2, pos2-seedLength,
				                          seedLength);
				}
			}

		// call the seed hit processor for this seed hit

		seed_search_count_stat (rawSeedHits);
		basesHit += (*processor) (processorInfo, pos1, pos2, seedLength);

#ifdef densityCheckDepth3
		if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
			return basesHit;
#endif // densityCheckDepth3
		}

	return basesHit;
	}


static u64 find_table_matches_resolve
   (u32		packed2,
	unspos	pos2,
	u32		unpacked2,
	int		transAllowed)
	{
	u32		seedLength, len1;
	unspos	adjStart = pt->adjStart;
	u32		step     = pt->step;
	unspos	pos, pos1, pos1Rel;
	u64		basesHit = 0;
	u32		unpacked1;
	int		mismatches;

	seedLength = (unsigned) hitSeed->length;
	len1       = seedLength-1;

	if (pt->last[packed2] == 0)
		return 0;

	for (pos=pt->last[packed2] ; pos!=noPreviousPos ; pos=pt->prev[pos])
		{
		pos1Rel = step*pos;
		pos1    = adjStart + pos1Rel;

		if ((selfCompare) && (seed_hit_below_diagonal (pos1, pos2)))
			continue;

		// resolve the remaining seed bits

#ifdef debugResolvingSeeds
		debugThisHit = (pos1 == debugPos1+len1) && (pos2 == debugPos2+len1);
#endif
		if (debugThisHit)
			mismatches = 0;		// (only to set a breakpoint here)

		unpacked1 = fetch_resolving_bits (pt, pos1Rel);

		if (debugThisHit)
			{
			printf ("\npartial seed hit " unsposSlashFmt "\n", pos1-len1, pos2-len1);
			dump_aligned_nucleotides (stdout,
			                          seq1, pos1-seedLength,
			                          seq2, pos2-seedLength,
			                          seedLength);
			printf ("  %08X %s\n", unpacked1, bits_to_nuc_string(unpacked1,16));
			printf ("  %08X %s\n", unpacked2, bits_to_nuc_string(unpacked2,16));
			}

		unpacked1 ^= unpacked2;					// combine bits into A-B-C-D-
		unpacked1 &= hitSeed->resolvingMask;	// ... where any 1 => mismatch
		unpacked1 += unpacked1 >> 17;			// shift bits to ----CADB
		mismatches = bit_count_16(unpacked1);	// count mismatches
		if (mismatches > transAllowed)
			{
			if (debugThisHit)
				printf ("  rejected (%d mismatches)\n", mismatches);
			seed_search_count_stat (unresolvedSeedHits);
			continue;
			}

		// seed is resolved, ship it

		if (debugThisHit)
			printf ("  accepted (%d mismatches)\n", mismatches);

		if (seed_search_dbgDumpRawHits)
			{
			if (seed_search_dbgShowRawHits)
				dump_raw_hit (stderr, pos1, pos2);
			else
				{
				printf ("\nraw seed hit " unsposSlashFmt "\n", pos1-len1, pos2-len1);
				dump_aligned_nucleotides (stdout,
				                          seq1, pos1-seedLength,
				                          seq2, pos2-seedLength,
				                          seedLength);
				}
			}

		// call the seed hit processor for this seed hit

		seed_search_count_stat (rawSeedHits);
		basesHit += (*processor) (processorInfo, pos1, pos2, seedLength);

#ifdef densityCheckDepth3
		if ((maxBasesAllowed > 0) && (basesHit > maxBasesAllowed))
			return basesHit;
#endif // densityCheckDepth3
		}

	return basesHit;
	}

//----------
// [[-- a seed hit processor function --]]
//
// process_for_plain_hit--
//	Process a seed hit for a given word, without bothering to check for
//	overlap with other hits, and without performing extension.
//
// Arguments and Return value: (see seed_search.h)
//
//----------
//
// Implementation:
//
// This is the simplest seed hit processor.  We do not use any of the diag
// hashing arrays (diagEnd, diagStart or diagActual).
//
//----------

u64 process_for_plain_hit
   (void*	_info,
	unspos	pos1,
	unspos	pos2,
	unspos	length)
	{
	hitprocsimple* info = (hitprocsimple*) _info;
	u32		basesHit;

	// filter by position (if specified)

	if ((info->hp.posFilter)
	 && (filter_seed_hit_by_pos (&info->hp, pos1, pos2, length)))
		return 0;

	// filter by match/transversion count (if specified)

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	if (seed_search_dbgShowHits)
		printf ("plain seed hit " unsposSlashFmt " (diag " sgnposFmt ")\n",
		        pos1-(length-1), pos2-(length-1), diagNumber (pos1, pos2));

	// report the hit

#ifdef snoopReporterCalls
	fprintf (stderr, "process_for_plain_hit reporting " unsposSlashFmt " #" unsposFmt " (to %p)\n",
	                 pos1, pos2, length, info->hp.reporter);
#endif
	basesHit = ((*info->hp.reporter) (info->hp.reporterInfo, pos1, pos2, length, 0));
	if (basesHit > 0) searchToGo--;
	return basesHit;
	}

//----------
// [[-- a seed hit processor function --]]
//
// process_for_simple_hit--
//	Process a seed hit for a given word, with one hit "good enough".
//
// Arguments and Return value: (see seed_search.h)
//
//----------
//
// Implementation:
//
// Note that seed hits arrive in increasing positions on sequence 2, thus the
// arrivals on any particular diagonal are also increasing.
//
// We record and report seed hits as we encounter them, but we do not report
// overlapping hits.  When a new seed hit arrives, if it overlaps the most
// recent on that hash-equivalent diagonal we record the new end but otherwise
// ignore the new seed hit.  This treatment thus "suffers" from undetected hash
// collisions.
//
// Note that we do not use the diagStart or diagActual arrays.
//
//----------

u64 process_for_simple_hit
   (void*	_info,
	unspos	pos1,
	unspos	pos2,
	unspos	length)
	{
	hitprocsimple* info = (hitprocsimple*) _info;
	u32		hDiag;
	score	s;
	u32		basesHit;
#ifdef snoopDiagHash
	unspos	start2 = pos2 - length;
#endif // snoopDiagHash

	// filter by position (if specified)

	if ((info->hp.posFilter)
	 && (filter_seed_hit_by_pos (&info->hp, pos1, pos2, length)))
		return 0;

	unblockedLeftExtension = false;

	// if we've already extended beyond this point on this hash-equivalent
	// diagonal, ignore this hit

	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		printf ("simp: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " " unsposDotsFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, pos2-length, pos2, diagEnd[hDiag]);
		if (diagEnd[hDiag] == hashInactiveEnd) printf ("simp: first hit on diagonal\n");
		else if (diagEnd[hDiag] > pos2-length) printf ("simp: hit discarded\n");
		}
#endif

	if (diagEnd[hDiag] == hashInactiveEnd)
		{
#ifdef snoopDiagHash
		fprintf (stderr, "  activating diag %9s"
		                 "                              "
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		activate_hashed_diag (hDiag);
		diagEnd[hDiag] = 0;
		}

	if (diagEnd[hDiag] > pos2-length)
		{
#ifdef snoopDiagHash
		fprintf (stderr, "  ignoring   diag %9s"
		                 "                              "
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		return 0;
		}

	// filter by match/transversion count (if specified)

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	// perform gap-free extension (if specified) and record the extent of seed
	// hits on this diagonal;  note that the extention routines (such as
	// match_extend_seed_hit) will record the extent of the extended hit

	if ((info->hp.gfExtend != gfexNoExtend) && (seed_search_dbgShowHits))
		{
		int isRev1 = ((seq1->revCompFlags & rcf_rev) != 0);
		int isRev2 = ((seq2->revCompFlags & rcf_rev) != 0);
		printf ("simple seed hit " unsposSlashCFmt " (diag " sgnposFmt ")\n",
		        pos1-(length-1), (isRev1)?'-':'+',
		        pos2-(length-1), (isRev2)?'-':'+',
		        diagNumber (pos1, pos2));
		}

	if (info->hp.gfExtend == gfexExact)
		{
		s = match_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if (info->hp.gfExtend == gfexXDrop)
		{
		s = xdrop_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if ((info->hp.gfExtend >= gfexMismatch_min)
	      && (info->hp.gfExtend <= gfexMismatch_max))
		{
		s = mismatch_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else // if (info->hp.gfExtend == gfexNoExtend)
		{
		diagEnd[hDiag] = pos2;
		s = 0;
#ifdef snoopDiagHash
	fprintf (stderr, "  setting    diag %9s"
	                 "                              "
	                 ", diagEnd[%04X] = " unsposFmt
	                 ", seed end = " unsposSlashFmt
	                 ", in seq 1: " unsposDotsFmt "\n",
	                 pair_diagonal_as_text(pos1,pos2),
	                 hDiag, diagEnd[hDiag],
	                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

	// report the hit

#ifdef snoopReporterCalls
	fprintf (stderr, "process_for_simple_hit reporting " unsposSlashFmt " #" unsposFmt " (to %p)\n",
	                 pos1, pos2, length, info->hp.reporter);
#endif
	basesHit = ((*info->hp.reporter) (info->hp.reporterInfo, pos1, pos2, length, s));
	if (basesHit > 0) searchToGo--;
	return basesHit;
	}

//----------
// [[-- a seed hit processor function --]]
//
// process_for_recoverable_hit--
//	Process a seed hit for a given word, with one hit "good enough", recovering
//	from hash collisions.
//
// Arguments and Return value: (see seed_search.h)
//
//----------
//
// Implementation:
//
// Note that seed hits arrive in increasing positions on sequence 2, thus the
// arrivals on any particular diagonal are also increasing.
//
// We record and report seed hits as we encounter them, but we may report hits
// that overlap.  When a new seed hit arrives, if it overlaps the most recent
// on that hash-equivalent diagonal AND is on a different actual diagonal, we
// treat it as a new hit.  This leads to potential reporting of duplicate HSPs
// (for example, if we get a hit on diag A, then B, then A again).  Overlaps
// that have the same actual diagonal are reported on the first hit.
//
// Note that we do not use the diagStart array.
//
//----------

u64 process_for_recoverable_hit
   (void*	_info,
	unspos	pos1,
	unspos	pos2,
	unspos	length)
	{
	hitprocsimple* info = (hitprocsimple*) _info;
	unspos	start2 = pos2 - length;
	sgnpos	diag;
	u32		hDiag;
	score	s;
	u32		basesHit;
#ifdef debugDiag
	static char prevSeq2Name[41] = "";
	static int  prevRevCompFlags;
#endif
#ifdef snoopDiagHash
	seq*	seq1 = info->hp.seq1;
	seq*	seq2 = info->hp.seq2;
#endif // snoopDiagHash

	// filter by position (if specified)

	if ((info->hp.posFilter)
	 && (filter_seed_hit_by_pos (&info->hp, pos1, pos2, length)))
		return 0;

	//////////
	// decide whether to discard this hit, based on the extent of previous hits
	// along a hash-equivalent diagonal
	//////////

	unblockedLeftExtension = true;

	// get the diagonal's hash value

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		if ((prevSeq2Name[0] == 0)
		 || (strcmp (seq2->header, prevSeq2Name) != 0)
		 || (seq2->revCompFlags != prevRevCompFlags))
		 	{
		 	strncpy (prevSeq2Name, seq2->header, sizeof(prevSeq2Name));
		 	prevRevCompFlags = seq2->revCompFlags;
			if      (prevRevCompFlags == rcf_forward) printf ("%s+\n", prevSeq2Name);
			else if (prevRevCompFlags == rcf_revcomp) printf ("%s-\n", prevSeq2Name);
			else                                      printf ("%s\n", prevSeq2Name);
			}

		printf ("sing: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " " unsposDotsFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, start2, pos2, diagEnd[hDiag]);
		}
#endif

	// if the diagonal was inactive, we treat it as a fresh hit

	if (diagEnd[hDiag] == hashInactiveEnd)
		{
		activate_hashed_diag (hDiag);
		diagEnd[hDiag] = 0;
		goto fresh_hit;
		}

	// if we have a collision, accept it as a fresh hit;  note that accepting
	// prevents us from being able to recognize a later hit on that same
	// (previous) diagonal later, in which case we would end up reporting it
	// more than once;  the premise of this routine is that multiple-extension
	// is preferable to missing a colliding hit

	if (diag != diagActual[hDiag])
		{
#ifdef snoopDiagHash
		fprintf (stderr, "%s:\n", seq2->header);
#endif // snoopDiagHash
		seed_search_count_stat (hashCollisions);
		if (start2 < diagEnd[hDiag])
			{
			seed_search_count_stat (hashFailures);
#ifdef debugDiag
			if (debugThisDiag)
				printf ("sing:    accepted in spite of hash failure\n");
#endif
#ifdef snoopDiagHash
			fprintf (stderr, "  recovery  on diag %9s"
			                 ", diagActual[%04X] = %9s"
			                 ", diagEnd[%04X] = " unsposFmt
			                 ", seed end = " unsposSlashFmt
			                 ", in seq 1: " unsposDotsFmt "\n",
			                 pair_diagonal_as_text(pos1,pos2),
			                 hDiag, diagonal_as_text(diagActual[hDiag]),
			                 hDiag, diagEnd[hDiag],
			                 pos1, pos2, start2, pos2);
			}
		else
			{
			fprintf (stderr, "  collision on diag %9s"
			                 ", diagActual[%04X] = %9s"
			                 ", diagEnd[%04X] = " unsposFmt
			                 ", seed end = " unsposSlashFmt
			                 ", in seq 1: " unsposDotsFmt "\n",
			                 pair_diagonal_as_text(pos1,pos2),
			                 hDiag, diagonal_as_text(diagActual[hDiag]),
			                 hDiag, diagEnd[hDiag],
			                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
			}
		goto fresh_hit;
		}

	// if this hit overlaps an earlier hit on the same diagonal, reject it;
	// note that we record the extent, but only if it increases it

	if (start2 < diagEnd[hDiag])
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("sing:    rejected (diagEnd[%04X] blocks at " unsposFmt ")\n",
			        hDiag, diagEnd[hDiag]);
#endif
		if (pos2 > diagEnd[hDiag])
			{
			diagEnd   [hDiag] = pos2;
			diagActual[hDiag] = diag;
#ifdef debugDiag
			if (debugThisDiag)
				printf ("sing: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
						pair_diagonal_as_text(pos1,pos2), pos1, pos2,
						hDiag, diagEnd[hDiag]);
#endif
			}
		return 0;
		}

	//////////
	// this hit is a keeper, as far as diagonal extent is concerned;  now
	// perform whatever other filtering is specified, extend it, and record
	// the extent
	//
	// note that we have to be careful and only record the extent if it would
	// increase it;  a previous hit on a hash-equivalent diagonal may have
	// already extended further
	//////////

fresh_hit:
	diagActual[hDiag] = diag;

	// filter by match/transversion count (if specified)

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		{
		if (pos2 > diagEnd[hDiag]) diagEnd[hDiag] = pos2;
		return 0;
		}

	// perform gap-free extension (if specified) and record the extent of seed
	// hits on this diagonal;  note that the extention routines (such as
	// match_extend_seed_hit) will record the extent of the extended hit

	if ((info->hp.gfExtend != gfexNoExtend) && (seed_search_dbgShowHits))
		{
		int isRev1 = ((seq1->revCompFlags & rcf_rev) != 0);
		int isRev2 = ((seq2->revCompFlags & rcf_rev) != 0);
		printf ("simple seed hit " unsposSlashCFmt " (diag " sgnposFmt ")\n",
		        pos1-(length-1), (isRev1)?'-':'+',
		        pos2-(length-1), (isRev2)?'-':'+',
		        diag);
		}

	if (info->hp.gfExtend == gfexExact)
		{
		s = match_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if (info->hp.gfExtend == gfexXDrop)
		{
		s = xdrop_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if ((info->hp.gfExtend >= gfexMismatch_min)
	      && (info->hp.gfExtend <= gfexMismatch_max))
		{
		s = mismatch_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else // if (info->hp.gfExtend == gfexNoExtend)
		{
		if (pos2 > diagEnd[hDiag]) diagEnd[hDiag] = pos2;
		s = 0;
		}

#ifdef debugDiag
	if (debugThisDiag)
		printf ("sing: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
				pair_diagonal_as_text(pos1,pos2), pos1, pos2,
				hDiag, diagEnd[hDiag]);
#endif

	if ((seed_search_dbgShowHits || debugThisDiag))
		dump_aligned_nucleotides (stdout,
		                          seq1, pos1-length,
		                          seq2, pos2-length,
		                          length);

#ifdef snoopReporterCalls
	fprintf (stderr, "process_for_recoverable_hit reporting " unsposSlashFmt " #" unsposFmt " (to %p)\n",
	                 pos1, pos2, length, info->hp.reporter);
#endif
	basesHit = ((*info->hp.reporter) (info->hp.reporterInfo, pos1, pos2, length, s));
	if (basesHit > 0) searchToGo--;
	return basesHit;
	}

//----------
// [[-- a seed hit processor function --]]
//
// process_for_twin_hit--
//	Process a seed hit for a given word, with two nearby hits required for a
//	hit to be "good enough".
//
// Arguments and Return value: (see seed_search.h)
//	The info argument points to a hitproctwin record specifying the twin's span
//	criteria.
//
//----------
//
// Implementation:
//
// Note that seed hits arrive in increasing positions on sequence 2, thus the
// arrivals on any particular diagonal are also increasing.
//
// We record seed hits as we encounter them, but delay reporting them until we
// get a second hit within the (minSpan,maxSpan) criteria.  The span of two hits
// is the length from the start of the first hit to the end of the second.
//
// When a new seed hit arrives, we check the distance between it and the most
// recent hit on the same diagonal.  If it is too far we record the new hit (and
// erase the previous one).  If it is too close we just extend the previous
// hit (record the new end).  If it is within (minSpan,maxSpan) we extend the
// previous hit and report it.  We will continue to extend the hit as long as
// new hits are close enough to the end, but will only report it once.
//
// Note that we expect that length < minSpan <= maxSpan.
//
//----------
//
// Gap computation example:
//
// Case 1: In the figure below, we have a seed of length 10, and we get simple
// hits at (end points) 11 and 36 (stars show hits, arrows show the position
// reported to this routine).  Suppose the criteria for twins is 20<=span<=30.
// The span of hits A and B is 35, so they are too far apart to qualify.
//
//		                         1         2         3         4
//		pos2:           1234567890123456789012345678901234567890
//		sequence:     ..XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
//		seed hit A:     **********^
//		seed hit B:                              **********^
//      AB span (35):   ===================================
//
// Case 2: Now suppose we had had an intervening simple hit at 18 (below).  The
// span of CD is 17, which is too short.  But the span of DE is 28, which sat-
// isfies the criteria.
//		                         1         2         3         4
//		pos2:           1234567890123456789012345678901234567890
//		sequence:     ..XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
//		seed hit C:     **********^
//		seed hit D:            **********^
//		seed hit E:                              **********^
//		CD span (17):   =================
//      DE span (28):          ============================
//
// The implementation here accumulates simple hits into an unresolved hit as
// long as the new hit is unsatisfactory when compared to the first hit and to
// the last hit.  It may be possible to contrive a 4-hit example in which hits
// 2 and 4 are satifactory but which this implementation will not report.  The
// author presently believes these are not a problem in practice.
//
//----------

//--- (implementation NOT using the seed hit queue) ---

#ifdef noSeedHitQueue

#error ***** non-seed-queue version of process_for_twin_hit() has a serious flaw *****

// The flaw is that diagEnd is used here to track the end of the most recent seed,
// while in xdrop_extend_seed_hit() it is used to block the left-extension of a
// seed hit.  Thus when we do find a valid twin seed, left-extension stops at the
// RIGHT end of the first seed of the pair.  So we miss some portion of the HSP
// and also incorrectly reject some HSPs because the right-extension is not enough,
// by itself, to meet the score threshold.  See the two "(PROBLEM)" notes below;
// these indicate where diagEnd is set in a way that defeats extension.

u64 process_for_twin_hit
   (void*		_info,
	unspos		pos1,
	unspos		pos2,
	unspos		length)
	{
	hitproctwin* info = (hitproctwin*) _info;
	unspos		start2 = pos2 - length;
	sgnpos		diag;
	u32			hDiag;
	u32			span;
	score		s;

	// filter by position (if specified)

	if ((info->hp.posFilter)
	 && (filter_seed_hit_by_pos (&info->hp, pos1, pos2, length)))
		return 0;

	//////////
	// decide whether to discard this hit, based on the extent of previous hits
	// and twinliness along a hash-equivalent diagonal
	//////////

	unblockedLeftExtension = false;

	// get the diagonal's hash value

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		printf ("twin: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " " unsposDotsFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, start2, pos2, diagEnd[hDiag]);
		}
#endif

	// if the diagonal was inactive, we treat it as a fresh hit;  note that we
	// might not really have to activate the diagonal here, because this hit
	// may get discarded;  but the penalty for unecessary activation is only
	// one extra pass through the loop in empty_diag_hash(), vs. having to do
	// this test again after we filter at fresh_hit

	if (diagEnd[hDiag] == hashInactiveEnd)
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    first hit on diagonal\n");
#endif
		activate_hashed_diag (hDiag);
		diagEnd[hDiag] = 0;
		goto fresh_hit;
		}

	// if we have a collision, reject/accept it based on whether we have
	// already reported a twin for the first hit;  once we have reported a twin,
	// it is safe to discard that hit in favor of the new one;  but *before*
	// that event, we must discard the new hit;  the reason is that if we
	// accept the new hit, we would be vulnerable to a situation in which two
	// orthologies colliding on hash-equivalent diagonals cause us to discard
	// the other before a twin is detected, resulting in us missing both of
	// them

	if (diag != diagActual[hDiag])
		{
		seed_search_count_stat (hashCollisions);
		if (pos2 >= diagEnd[hDiag] - length + info->maxSpan)
			goto fresh_hit; // (beyond maxSpan from end of previous hit)

		seed_search_count_stat (hashFailures);

		span = diagEnd[hDiag] - diagStart[hDiag]; // (span of old hit(s))
		if (span >= info->maxSpan)
			{
#ifdef debugDiag
			if (debugThisDiag)
				printf ("twin:    accepted in spite of hash failure"
				        " (span was " unsposFmt ", extent is " unsposFmt ")\n",
				        span, pos2 - (diagEnd[hDiag] - length));
#endif
			goto fresh_hit; // (old hit, on different diag, already reported)
			}

#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    rejected (hash failure on pending twin,"
			        " span was " unsposFmt ", extent is " unsposFmt ")\n",
			        span, pos2 - (diagEnd[hDiag] - length));
#endif

		return 0;
		}

	// this hit is on the same diagonal as the previous one;  if the last hit
	// already extends beyond this one, we can ignore this one;  this can happen
	// if it was extended (by xdrop_extend_seed_hit), in which case it has
	// already been reported

	if (pos2 <= diagEnd[hDiag])
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    rejected (diagEnd[%04X] blocks at " unsposFmt ")\n",
			        hDiag, diagEnd[hDiag]);
#endif
		return 0;
		}

	// if the span of the last hit with the new hit is too large, record this
	// as a fresh hit;  this is like seed hit B in case 1

	span = length + pos2 - diagEnd[hDiag];
	if (span > info->maxSpan)
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    span too long (" unsposFmt " > " unsposFmt ")\n", span, info->maxSpan);
#endif
		goto fresh_hit;
		}

	// otherwise, we will extend the previous hit, but we have to decide whether
	// to report it;  if that hit was previously reported, we just extend it;
	// this is like seed hits *after* E in case 2

	if (diagEnd[hDiag] - diagStart[hDiag] >= info->minSpan)
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    already reported (" unsposFmt " >= " unsposFmt ")\n",
			        diagEnd[hDiag]-diagStart[hDiag], info->minSpan);
#endif
		goto simple_extend_hit;
		}

	// if the combined length of the previous hit with this one added is too
	// short, we have not reached minSpan yet, so we just extend it;  this is
	// like seed hit D in case 2

	span = pos2 - diagStart[hDiag];
	if (span < info->minSpan)
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    not long enough yet (" unsposFmt " < " unsposFmt ")\n",
			        span, info->minSpan);
#endif
		goto simple_extend_hit;
		}

	// otherwise, the gap has met the (minSpan,maxSpan) citeria for the first
	// time;  this is like seed hit E in case 2

	goto fresh_twin_hit;

	//////////
	// this is a fresh hit, as far as diagonal extent is concerned;  perform
	// whatever other filtering is specified before recording it
	//////////

	// filter by match/transversion count (if specified)

fresh_hit:

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	// record it

	diagStart [hDiag] = start2;
	diagActual[hDiag] = diag;
#ifdef debugDiag
	if (debugThisDiag)
		printf ("twin: (diag %9s)      " unsposSlashFmt " diagStart[%04X] <-- " unsposFmt "\n",
		        pair_diagonal_as_text(pos1,pos2), pos1, pos2,
		        hDiag, diagStart[hDiag]);
#endif
	goto record_extent;

	//////////
	// this hit only requires that we record it's extent, as far as diagonal
	// extent is concerned;  perform whatever other filtering is specified
	// before recording it
	//////////

	// filter by match/transversion count (if specified)

simple_extend_hit:

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	// record it

record_extent:
	diagEnd[hDiag] = pos2;  // (PROBLEM)
#ifdef debugDiag
	if (debugThisDiag)
		printf ("twin: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
		        pair_diagonal_as_text(pos1,pos2), pos1, pos2,
		        hDiag, diagEnd[hDiag]);
#endif

	return 0;

	//////////
	// this hit is a keeper, as far as diagonal extent is concerned;  now
	// perform whatever other filtering is specified, extend it, and record
	// the extent
	//////////

fresh_twin_hit:

	// filter by match/transversion count (if specified)

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	// perform gap-free extension (if specified) and record the extent of seed
	// hits on this diagonal;  note that the extention routines (such as
	// match_extend_seed_hit) will record the extent of the extended hit

	length = span;

	if ((info->hp.gfExtend != gfexNoExtend) && (seed_search_dbgShowHits))
		printf ("twin seed hit " unsposSlashFmt "\n", pos1-(span-1), pos2-(span-1));

	if (info->hp.gfExtend == gfexExact)
		{
		s = match_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if (info->hp.gfExtend == gfexXDrop)
		{
		s = xdrop_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else if ((info->hp.gfExtend >= gfexMismatch_min)
	      && (info->hp.gfExtend <= gfexMismatch_max))
		{
		s = mismatch_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (s == noScore)
			return 0;
		}
	else // if (info->hp.gfExtend == gfexNoExtend)
		{
		diagEnd[hDiag] = pos2;  // (PROBLEM)
		s = 0;
		}

#ifdef debugDiag
	if (debugThisDiag)
		printf ("twin: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
		        pair_diagonal_as_text(pos1,pos2), pos1, pos2,
		        hDiag, diagEnd[hDiag]);
#endif

	if (seed_search_dbgShowHits || debugThisDiag)
		dump_aligned_nucleotides (stdout,
		                          seq1, pos1-length,
		                          seq2, pos2-length,
		                          length);

#ifdef snoopReporterCalls
	fprintf (stderr, "process_for_twin_hit reporting " unsposSlashFmt " #" unsposFmt " (to %p)\n",
	                 pos1, pos2, length, info->hp.reporter);
#endif
	return ((*info->hp.reporter) (info->hp.reporterInfo, pos1, pos2, length, s));
	}

#endif // noSeedHitQueue


//--- (implementation using the seed hit queue) ---

#ifndef noSeedHitQueue

u64 process_for_twin_hit
   (void*		_info,
	unspos		pos1,
	unspos		pos2,
	unspos		length)
	{
	hitproctwin* info = (hitproctwin*) _info;
	sgnpos		diag;
	u32			hDiag;
	unspos		oldDiagEnd, extent;
	u64			num;
	shqhit*		q;
	u32			span;
	unspos		start2;
	score		s;
#ifdef debugDiag
	u32			longestShortSpan = 0;
#endif

	// filter by position (if specified)

	if ((info->hp.posFilter)
	 && (filter_seed_hit_by_pos (&info->hp, pos1, pos2, length)))
		return 0;

	// filter by match/transversion count (if specified)

	if ((info->hp.minMatches >= 0)
	 && (filter_seed_hit_by_subs (&info->hp, pos1, pos2, length)))
		return 0;

	//////////
	// scan the seed hit queue for hits along this diagonal with valid span, or
	// for 'blocks' placed at the end of previously extended diagonals
	//////////

	unblockedLeftExtension = false;

	// get the diagonal's hash value

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		start2 = pos2 - length;
		printf ("\ntwin: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " " unsposDotsFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, start2, pos2, diagEnd[hDiag]);
		}
#endif

	// if the diagonal was inactive, activate it, then add this hit to the queue,
	// and exit

	if (diagEnd[hDiag] == hashInactiveEnd)
		{
#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin:    first hit on diagonal\n");
#endif
		activate_hashed_diag (hDiag);
		diagEnd[hDiag] = 0;

		enqueue_seed_hit (pos1, pos2, /*isBlock*/ false);

#ifdef debugDiag
		if (debugThisDiag)
			printf ("twin: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
					pair_diagonal_as_text(pos1,pos2), pos1, pos2,
					hDiag, diagEnd[hDiag]);
#endif

		return 0;
		}

	// scan the queue entries for this hashed diagonal, until one of the
	// following occurs:
	//	- we find a twin hit with an acceptable span
	//	- we reach the end of previously extended seed hits on this diagonal;
	//    any additional seed hits in the queue would be older than this block
	//    and thus needn't be considered
	//	- we reach the 'end' of the queue;  this is accomplished by checking
	//    the scanned record's 'number' to see if it would still be in the
	//    queue

	for (num = lastSeedHit[hDiag] ;
		 num > seedHitNum - seedHitQueueSize ;
		 num = q->prevHit)
		{
		q = &seedHitQueue[num % seedHitQueueSize];
		seed_search_count_stat (queueSeedsScanned);
		span = pos2 - (q->pos2 - length);
		if (span > info->maxSpan)
			{			// (span too long from previous hit)
#ifdef debugDiag
			if ((debugThisDiag) && (q->diag == diag))
				printf ("twin:    span too long (" unsposFmt " > " unsposFmt ")\n", span, info->maxSpan);
#endif
			break;
			}

		if (q->diag != diag)	// (not on correct diagonal)
			continue;

		seed_search_count_stat (queueSeedsExamined);

		if (q->isBlock)			// (we've already extended this far)
			{
			start2 = pos2 - length;
			if (start2 <= q->pos2)
				{
				seed_search_count_stat (queueSeedsBlocked);
				return 0;		// (this seed hit overlaps previous extension)
				}
			else
				break;			// (this seed hit is right of extension)
			}

		if (span < info->minSpan) // (span not long enough from previous hit)
			{
#ifdef debugDiag
			if ((debugThisDiag) && (span > longestShortSpan))
				longestShortSpan = span;
#endif
			continue;
			}

		goto twin_hit;	// (this hit is part of a twin with a desired span)
		}

#ifdef debugDiag
	if ((debugThisDiag) && (longestShortSpan > 0))
		printf ("twin:    not long enough yet (" unsposFmt " < " unsposFmt ")\n",
				longestShortSpan, info->minSpan);
#endif

	// we don't have a twin;  add this seed hit to the queue, and exit

	enqueue_seed_hit (pos1, pos2, /*isBlock*/ false);
	return 0;

	//////////
	// this hit is part of a twin with a desired span;  now extend it and
	// record the extent
	//////////

	// perform gap-free extension (if specified) and record the extent of seed
	// hits on this diagonal;  note that the extention routines (such as
	// match_extend_seed_hit) will record the extent of the extended hit

twin_hit:
#ifdef debugDiag
	if (debugThisDiag)
		printf ("twin:    twin hit, length = %u\n", span);
#endif

	length = span;

	if ((info->hp.gfExtend != gfexNoExtend) && (seed_search_dbgShowHits))
		printf ("twin seed hit " unsposSlashFmt "\n", pos1-(span-1), pos2-(span-1));

	if (info->hp.gfExtend == gfexExact)
		{
		oldDiagEnd = diagEnd[hDiag];
		s = match_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (diagEnd[hDiag] != oldDiagEnd)
			{
			extent = diagEnd[hDiag];
			enqueue_seed_hit (diagToPos1(diag,extent), extent, /*isBlock*/ true);
			if (s == noScore)
				enqueue_seed_hit (pos1, pos2, /*isBlock*/ false);
			}
		if (s == noScore)
			return 0;
		}
	else if (info->hp.gfExtend == gfexXDrop)
		{
		oldDiagEnd = diagEnd[hDiag];
		s = xdrop_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (diagEnd[hDiag] != oldDiagEnd)
			{
			extent = diagEnd[hDiag];
			enqueue_seed_hit (diagToPos1(diag,extent), extent, /*isBlock*/ true);
			}
		if (s == noScore)
			return 0;
		}
	else if ((info->hp.gfExtend >= gfexMismatch_min)
	      && (info->hp.gfExtend <= gfexMismatch_max))
		{
		oldDiagEnd = diagEnd[hDiag];
		s = mismatch_extend_seed_hit (&info->hp, &pos1, &pos2, &length);
		if (diagEnd[hDiag] != oldDiagEnd)
			{
			extent = diagEnd[hDiag];
			enqueue_seed_hit (diagToPos1(diag,extent), extent, /*isBlock*/ true);
			if (s == noScore)
				enqueue_seed_hit (pos1, pos2, /*isBlock*/ false);
			}
		if (s == noScore)
			return 0;
		}
	else // if (info->hp.gfExtend == gfexNoExtend)
		{
		diagEnd[hDiag] = pos2;
		enqueue_seed_hit (pos1, pos2, /*isBlock*/ true);
		s = 0;
		}

#ifdef debugDiag
	if (debugThisDiag)
		printf ("twin: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt "\n",
		        pair_diagonal_as_text(pos1,pos2), pos1, pos2,
		        hDiag, diagEnd[hDiag]);
#endif

	if (seed_search_dbgShowHits || debugThisDiag)
		dump_aligned_nucleotides (stdout,
		                          seq1, pos1-length,
		                          seq2, pos2-length,
		                          length);

#ifdef snoopReporterCalls
	fprintf (stderr, "process_for_twin_hit reporting " unsposSlashFmt " #" unsposFmt " (to %p)\n",
	                 pos1, pos2, length, info->hp.reporter);
#endif
	return ((*info->hp.reporter) (info->hp.reporterInfo, pos1, pos2, length, s));
	}

#endif // not noSeedHitQueue

//----------
//
// seed_hit_below_diagonal--
//	Determine whether a raw seed hit is below (or on) the diagonal.
//
// Arguments:
//	unspos	pos1, pos2:	The hit position in sequence 1 and 2.  This is the
//						.. position following the end of the hit (and origin
//						.. zero).
//
// Returns:
//	true if the hit is "below the digaonal", false if it is not.
//
//----------
//
// Notes:  [ similar notes appear in mirror_alignments() ]
//
//	(1) We assume, without checking, that seq1 and seq2 are essentially the
//		same.  I.e. that they have the same length, and if one is partitioned,
//		the other has the same partitions.
//
//	(2) The DP matrix is viewed as having sequence 1 along the x axis and
//		sequence 2 along the y axis, as in this diagram:
//
//			      +-------------+
//			  ^   | . . . . . / |
//			  |   | . . . . /   |
//			  |   | . . . /     |
//			seq 2 | . . /       |
//			  |   | . /         |
//			  |   | /           |
//			      +-------------+
//			       --- seq 1 -->
//
//	(3) The diagonal runs from lower-left to upper-right, shown as the slashed
//		line in the diagrom.
//
//	(4) Alignments "above the diagonal" have pos1 < pos2.  In the diagram, this
//		is the region filled with dots.  Alignments "below the diagonal" or "on
//		the diagonal" have pos1 >= pos2;  the region is empty in the diagram.
//
//	(5) When sequence 2 is on the minus strand, pos2 is actually in counted in
//		reverse.
//
//			       (conceptual)              (actual)
//			      +-------------+	      +-------------+
//			  ^   | . . . . . / |	  ^   | \ . . . . . |
//			  |   | . . . . /   |	  |   |   \ . . . . |
//			  |   | . . . /     |	  |   |     \ . . . |
//			seq 2 | . . /       |	seq 2 |       \ . . |
//			  |   | . /         |	  |   |         \ . |
//			  |   | /           |	  |   |           \ |
//			      +-------------+	      +-------------+
//			       --- seq 1 -->	       --- seq 1 -->
//
//	(6) When sequence 2 is partitioned, and on the minus strand, the situation
//		with positions in complicated by the fact that the partitions have been
//		reversed individually, not the sequence as a whole.
//
//			           (conceptual)	                     (actual)
//			      +-------------+-------+	      +---------------------+
//			  ^   | . . . . . . | . . / |	  ^   | . . . . . . | \ . . |
//			  |   | . . . . . . | . /   |	  |   | . . . . . . |   \ . |
//			  |   | . . . . . . | /     |	  |   | . . . . . . |     \ |
//			  |   +-------------+-------+	  |   +-------------+-------+
//			seq 2 | . . . . . / |       |	seq 2 | \ . . . . . |       |
//			  |   | . . . . /   |       |	  |   |   \ . . . . |       |
//			  |   | . . . /     |       |	  |   |     \ . . . |       |
//			  |   | . . /       |       |	  |   |       \ . . |       |
//			  |   | . /         |       |	  |   |         \ . |       |
//			  |   | /           |       |	  |   |           \ |       |
//			      +---------------------+	      +---------------------+
//			       --- seq 1 -->                   --- seq 1 -->
//
//----------

//=== stuff for snoopBelowDiagonal ===

#ifndef snoopBelowDiagonal
#define snoopBelowDiagonal_1 ;
#define snoopBelowDiagonal_2 ;
#define snoopBelowDiagonal_3 ;
#define snoopBelowDiagonal_4 ;
#define snoopBelowDiagonal_5 ;
#define snoopBelowDiagonal_6 ;
#define snoopBelowDiagonal_7 ;
#define snoopBelowDiagonal_8 ;
#endif // not snoopBelowDiagonal

#ifdef snoopBelowDiagonal

#define snoopBelowDiagonal_1                                                  \
	int nonTrivial = (pos1 != pos2);                                          \
	if (nonTrivial) fprintf (stderr, unsposSlashFmt,                          \
	                                 pos1-hitSeed->length,                    \
	                                 pos2-hitSeed->length);

#define snoopBelowDiagonal_2                                                  \
	if (nonTrivial) fprintf (stderr, " (same strand)");

#define snoopBelowDiagonal_3                                                  \
	if (nonTrivial)                                                           \
		fprintf (stderr, " --> %s\n",                                         \
		         (pos1 == pos2)? "on" : (pos1 >= pos2)? "below" : "above");

#define snoopBelowDiagonal_4                                                  \
	if (!nonTrivial) fprintf (stderr, unsposSlashFmt, pos1, pos2);            \
	nonTrivial = true;                                                        \
	fprintf (stderr, " (opposite strand)");

#define snoopBelowDiagonal_5                                                  \
	if (nonTrivial) fprintf (stderr, " (partitioned)");

#define snoopBelowDiagonal_6                                                  \
	if (nonTrivial) fprintf (stderr, " parts %ld/%ld",                        \
	                         part1 - sp1->p, part2 - sp2->p);

#define snoopBelowDiagonal_7                                                  \
	if (nonTrivial) fprintf (stderr, " --> %s (by part)\n",                   \
	                         (partIx1 >= partIx2)? "below" : "above");

#define snoopBelowDiagonal_8                                                  \
	if (nonTrivial)                                                           \
		fprintf (stderr, " " unsposSlashFmt " --> %s\n",                      \
		         pos1,pos2,                                                   \
		         (pos1 == pos2)? "on" : (pos1 >= pos2)? "below" : "above");

#endif // snoopBelowDiagonal


//=== seed_hit_below_diagonal() ===

static int seed_hit_below_diagonal
   (unspos			pos1,
	unspos			pos2)
	{
	seqpartition*	sp1, *sp2;
	partition*		part1, *part2;
	int				partIx1, partIx2;

	snoopBelowDiagonal_1;

	// same strand case;  see notes (2) thru (4)

	if (sameStrand)
		{
		snoopBelowDiagonal_2;
		snoopBelowDiagonal_3;
		return (pos1 >= pos2);
		}

	// opposite strand case;  see note (5)

	pos1 -= hitSeed->length;	// (note we don't bother to do this for same
	pos2 -= hitSeed->length;	//  .. strand since it doesn't affect the test)

	snoopBelowDiagonal_4;

	sp2 = &seq2->partition;
	if (sp2->p == NULL) // (seq2 is not partitioned)
		{
		pos2 = (seq2->len-1) - pos2;
		snoopBelowDiagonal_8;
		return (pos1 >= pos2);
		}

	// parititioned opposite strand case;  see note (6)

	snoopBelowDiagonal_5;

	sp1     = &seq1->partition;
	part1   = lookup_partition (seq1, pos1);
	part2   = lookup_partition (seq2, pos2);
	partIx1 = part1 - sp1->p;
	partIx2 = part2 - sp2->p;

	snoopBelowDiagonal_6;

	if (partIx1 != partIx2)
		{
		snoopBelowDiagonal_7;
		return (partIx1 >= partIx2);
		}

	pos2 = (part2->sepBefore + part2->sepAfter) - pos2;
	snoopBelowDiagonal_8;
	return (pos1 >= pos2);
	}

//----------
//
// filter_seed_hit_by_pos--
//	Determine whether a raw seed hit should be filtered, based on its position
//	in the target or query.
//
// Arguments:
//	hitprocinfo*	hp:		Pointer to record containing (among other things)
//							.. the filter criteria, tStart, tEnd, qStart and
//							.. qEnd.
//	unspos			pos1:	The hit position in sequence 1, relative to the
//							.. entire sequence (not to the interval).  This is
//							.. the first letter following the end of the match
//							.. (origin-0).
//	unspos			pos2:	The hit position in sequence 2 (with details the
//							.. same as for pos1).
//	unspos			length:	The length of the hit (number of nucleotides).
//
// Returns:
//	true if the hit should be filtered (discarded), false if it should be kept.
//
//----------
//
// Notes:
//
// (1)	The seed hit is discarded if it extends outside the allowable range on
//		.. either target (pos1) or query (pos2).
//
//----------

static int filter_seed_hit_by_pos
   (hitprocinfo*	hp,
	unspos			pos1,
	unspos			pos2,
	unspos			length)
	{
	unspos			tStart = hp->targetInterval.s;
	unspos			tEnd   = hp->targetInterval.e;
	unspos			qStart = hp->queryInterval.s;
	unspos			qEnd   = hp->queryInterval.e;

	pos1 -= length;		// (move from first location AFTER the
	pos2 -= length;		//  .. hit, to first location OF the hit)

#ifdef snoopPosFilter
	fprintf (stderr, "filter_seed_hit_by_pos(" unsposFmt "/" unsposFmt "#" unsposFmt ")"
	                 " vs [" unsposFmt "," unsposFmt "] / [" unsposFmt "," unsposFmt "]",
	                 pos1, pos2, length,
	                 tStart, tEnd, qStart, qEnd);

	if (pos1 < tStart)
		fprintf (stderr, " (discarded, before target start)\n");
	else if (pos1+length > tEnd)
		fprintf (stderr, " (discarded, beyond target end)\n");
	else if (pos2 < qStart)
		fprintf (stderr, " (discarded, before query start)\n");
	else if (pos2+length > qEnd)
		fprintf (stderr, " (discarded, beyond query end)\n");
	else
		fprintf (stderr, " (kept)\n");

#endif // snoopPosFilter

	// if the hit extends beyond target or query, discard it

	if ((pos1 < tStart) || (pos1+length > tEnd)) return true;
	if ((pos2 < qStart) || (pos2+length > qEnd)) return true;

	// otherwise, keep it

	return false;
	}

//----------
//
// filter_seed_hit_by_subs--
//	Determine whether a raw seed hit should be filtered, based on the number of
//	matches and transversions it contains.
//
// Arguments:
//	hitprocinfo*	hp:		Pointer to record containing (among other things)
//							.. the filter criteria, minMatches and
//							.. maxTransversions.
//	unspos			pos1:	The hit position in sequence 1.
//	unspos			pos2:	The hit position in sequence 2.
//	unspos			length:	The length of the hit.
//
// Returns:
//	true if the hit should be filtered (discarded), false if it should be kept.
//
//----------
//
// Notes:
//
// (1)	This is dependent on the specific 2-bit encoding of nucleotides, which
//		is defined (implicitly) in dna_utilities.c.  We assume that the least
//		significant of the two bits distinguishes between purines and
//		pyramidines, so if it is different in two nucleotides, they are a
//		transversion.
//
// (2)	This routine considers any substitution involving non-DNA characters to
//		be a transversion (anything outside of ACGTacgt).
//
//----------

#define is_transversion(bits1,bits2) ((((bits1)^(bits2))&1)==1)

static int filter_seed_hit_by_subs
   (hitprocinfo*	hp,
	unspos			pos1,
	unspos			pos2,
	unspos			length)
	{
	u8*				scan1 = hp->seq1->v + pos1;
	u8*				scan2 = hp->seq2->v + pos2;
	char*			pScan;
	unspos			remaining;
	int				matches, transversions;
	s8				bits1, bits2;

	// count the number of matches and transversions in this hit

	if (hp->filterPattern != NULL)
		{
		scan1 -= length;
		scan2 -= length;
		pScan =  hp->filterPattern;
		matches = transversions = 0;
		for (remaining=length ; remaining>0 ; remaining--)
			{
			if (*(pScan++) != '0')
				{
				bits1 = hp->charToBits[*scan1];
				bits2 = hp->charToBits[*scan2];

				if ((bits1 < 0) || (bits2 < 0))	// (negative => not ACGTacgt)
					transversions++;
				else if (bits1 == bits2)
					matches++;
				else if (is_transversion(bits1,bits2))
					transversions++;
				}
			scan1++;  scan2++;
			}
		}
	else
		{
		matches = transversions = 0;
		for (remaining=length ; remaining>0 ; remaining--)
			{
			bits1 = hp->charToBits[*(--scan1)];
			bits2 = hp->charToBits[*(--scan2)];

			if ((bits1 < 0) || (bits2 < 0))	// (negative => not ACGTacgt)
				transversions++;
			else if (bits1 == bits2)
				matches++;
			else if (is_transversion(bits1,bits2))
				transversions++;
			}
		}

	// if the counts don't meet the filter criteria, discard this hit

	if (matches < hp->minMatches)
		{ seed_search_count_stat (notEnoughMatches);      return true; }
	if ((hp->maxTransversions >= 0) && (transversions > hp->maxTransversions))
		{ seed_search_count_stat (tooManyTransversions);  return true; }

	// otherwise, keep it

	return false;
	}

//----------
//
// xdrop_extend_seed_hit--
//	Perform gap-free extension on a seed hit, and discard those that don't
//	score high enough.
//
//	Gap-free extension extends the hit/match in each direction, along the
//	diagonal, as long as the running score has not dropped too far below the
//	maximum score.
//
//	Low-scoring extensions are 'discarded' by the caller.  This routine only
//	makes the scoring decision.  Further, if an adpative scoring threshold is
//	being used, this routine will consider all extensions as being acceptible.
//
// Arguments:
//	hitprocinfo*	hp:		Pointer to record containing (among other things)
//							.. the extension controls and filtering criteria.
//	unspos*			pos1:	The hit position in sequence 1.  (see note 1 below)
//	unspos*			pos2:	The hit position in sequence 2.  (see note 1 below)
//	unspos*			length:	The length of the hit.           (see note 1 below)
//
// Returns:
//	The score of the extended match.  If the match should be filtered
//	(discarded), noScore is returned.
//
//----------
//
// Notes:
//
// (1)	Upon return, the values of pos1, pos2, and length will have been changed
//		to reflect the extended hit *IF* the match is to be kept (i.e. if the
//		return value is not noScore).
//
// (2)	We expect that the hits will arrive in increasing order on sequence 2.
//		See note 5.
//
// (3)	It may appear that we (incorrectly) assume that the entirety of
//		sequences 1 and 2 are fair game for the HSP.  However, we indirectly
//		halt processing when we encounter a NUL character.  This happens because
//		the scoring matrix contains veryBadScore for any substitution with NUL.
//		This in turn causes loop 1 or loop 2 to exit (because the running score
//		exceeds the xdrop threshold).  A NUL character indicates (1) the end of
//		a partition, either in a partitioned sequence or X-separated sequence,
//		or (2) the end of a chore.  Further, positional seed filtering
//		(performed by filter_seed_hit_by_pos) prevents us from being called for
//		hits outside the range of a chore.
//
// (4)	We assume hp->hspThreshold.t is either 'C' (count) or 'S' (score).
//
// (5)	Though hits arrive in increasing order on sequence 2, it is possible
//		that diagEnd[hDiag] > pos2.  This happens when a previous seed hit on a
//		hash-equivalent diagonal was extended (diagEnd records the righthand
//		limit of that extension).  In such case, left-extension of this seed hit
//		is prematurely halted (at pos2-length).  Right-extension is no affected.
//
//----------

#ifndef snoopXDrop
#define snoopXDrop_Pos   ;
#define snoopXDrop_Left  ;
#define snoopXDrop_Right ;
#define snoopXDrop_Score ;
#endif // not snoopXDrop

#ifdef snoopXDrop

#define snoopXDrop_Pos                                                       \
	fprintf (stderr,                                                         \
	         "xdrop, seed hit at " unsposSlashFmt "\tdiag=%s\n",             \
			 pos1, pos2, pair_diagonal_as_text(pos1,pos2));

#define snoopXDrop_Left                                                      \
	fprintf (stderr,                                                         \
			 "    left:  %s %s"                                              \
			 "  " scoreFmtStar                                               \
			 "  " scoreFmtStar "\n",                                         \
			 display_sequence_character (seq1, s1[-1]),                      \
			 display_sequence_character (seq2, s2[-1]),                      \
			 8, scoring->sub[s1[-1]][s2[-1]],                                \
			 8, runScore + scoring->sub[s1[-1]][s2[-1]]);

#define snoopXDrop_Right                                                     \
	fprintf (stderr,                                                         \
			 "    right: %s %s"                                              \
			 "  " scoreFmtStar                                               \
			 "  " scoreFmtStar "\n",                                         \
			 display_sequence_character (seq1, *s1),                         \
			 display_sequence_character (seq2, *s2),                         \
			 8, scoring->sub[*s1][*s2],                                      \
			 8, runScore + scoring->sub[*s1][*s2]);

#define snoopXDrop_Score                                                     \
	fprintf (stderr,                                                         \
	         "    score=" scoreFmtSimple "\n",                               \
			 similarity);

static char* display_sequence_character (seq* _seq, u8 ch)
	{
	static char	 s1[4];
	static char	 s2[4];
	static char* s = s2;

	s = (s == s1)? s2 : s1;	// (ping pong)

	if (_seq->fileType == seq_type_qdna) sprintf (s, "%02X", ch);
	                                else sprintf (s, "%c",   ch);

	return s;
	}

#endif // snoopXDrop


//--- xdrop_extend_seed_hit--

static score xdrop_extend_seed_hit
   (hitprocinfo*	hp,
	unspos*			_pos1,
	unspos*			_pos2,
	unspos*			_length)
	{
	unspos			pos1    = *_pos1;
	unspos			pos2    = *_pos2;
	unspos			length  = *_length;
	seq*			seq1    = hp->seq1;
	seq*			seq2    = hp->seq2;
	scoreset*		scoring = hp->scoring;
	score			xDrop   = hp->xDrop;
	sgnpos			diag, block2;
	unspos			oldDiagEnd, extent;
	u32				hDiag;
	u8*				s1, *s2, *stop, *leftStart, *rightStop, *rightBlock;
	score			similarity, runScore, leftScore, rightScore;
	int				adjustScore;
#ifdef debugDiag
	unspos			start2 = pos2 - length;
	u8*				p1, *p2, *p3, *p4;
#else
#ifdef collect_stats
	u8*				p1;
#endif
#endif
#ifdef snoopHspSubrange
	score			bestSubrangeScore;
	unspos			hspPos1, hspPos2, hspLen, subPos1, subPos2, subLen;
#endif
#ifdef snoopDiagHash
	unspos			start2 = pos2 - length;
#endif // snoopDiagHash

	//////////
	// get ready to extend the hit
	//////////

#ifdef debugDiag
	p1 = p2 = p3 = p4 = NULL;  // (satisfy the compiler)
#endif

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		printf ("gfex: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " " unsposDotsFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, start2, pos2, diagEnd[hDiag]);
		}
#endif

	snoopXDrop_Pos;

	//////////
	// extend to the left (loop 1)
	//
	// results:
	//	leftStart:  position of 1st nucleotide in extended match, in seq1
	//	leftScore:  score of bp added by left extension
	//	length:     possibly shortened (if the extension does not include the
	//	            .. entire hit)
	//////////

	s1 = seq1->v + pos1;	// start just past end of hit in both seq1 and seq2;
	s2 = seq2->v + pos2;	// .. we will pre-decrement before reads, so first
							// .. bp read are the ones at the tail of the hit

#ifdef extendHspFromLeft
	s1 = seq1->v + pos1 - length;
	s2 = seq2->v + pos2 - length;
#endif // extendHspFromLeft

	// determine stop location;  this is at the start of sequence 1, except
	// that if this diagonal ends (or is blocked) earlier in sequence 2, we
	// stop there
	// (see note 3;  instead of zero, we should use subsequence's start)

	if (unblockedLeftExtension) oldDiagEnd = 0;
	                       else oldDiagEnd = diagEnd[hDiag];
	block2 = (sgnpos) oldDiagEnd;
	if (block2 + diag > 0) stop = seq1->v + block2 + diag;
	                  else stop = seq1->v;

	// extend

	leftStart = s1;
	runScore  = leftScore = 0;

	while ((s1 > stop) && (runScore >= leftScore-xDrop))
		{
		snoopXDrop_Left;
		runScore += scoring->sub[*(--s1)][*(--s2)];
		if (runScore > leftScore)
			{
			leftStart = s1;
			leftScore = runScore;
			}
		}

	// adjust length if the extension is shorter than the hit

#ifndef extendHspFromLeft
	s2 = seq1->v + pos1 - length;	// (left end of hit)
	if (leftStart > s2)
		length -= leftStart - s2;
#endif // not extendHspFromLeft

#ifdef debugDiag
	p1 = s1;
	if (debugThisDiag)
		{
		p2 = leftStart;
		p3 = seq1->v+pos1-length;
		p4 = seq1->v+pos1;
		}
#endif
#ifdef collect_stats
	p1 = s1;
#endif

	//////////
	// extend to the right (loop 2)
	//
	// results:
	//	rightStop:  position of 1st nucleotide beyond extended match, in seq1
	//	similarity: increased by score of bp added by left and right extension
	//////////

	s1 = seq1->v + pos1;	// start just past end of hit in both seq1
	s2 = seq2->v + pos2;	// .. and seq2

#ifdef extendHspFromLeft
	s1 = seq1->v + pos1 - length;
	s2 = seq2->v + pos2 - length;
#endif // extendHspFromLeft

	// determine stop location;  this is at the end of sequence 1, except
	// that if this diagonal ends earlier in sequence 2, we stop there
	// (see note 3;  instead of sequence's end, we should use subsequence's end)

	block2 = (sgnpos) seq2->len;
	if ((sgnpos) seq1->len <= block2 + diag) stop = seq1->v + seq1->len;
	                                    else stop = seq1->v + block2 + diag;

	// extend

	rightStop = s1;
	runScore = rightScore = 0;

	while ((s1 < stop) && (runScore >= rightScore-xDrop))
		{
		snoopXDrop_Right;
		runScore += scoring->sub[*(s1++)][*(s2++)];
		if (runScore > rightScore)
			{
			rightStop  = s1;
			rightScore = runScore;
			}
		}
	rightBlock = s1;

	// adjust length if the extension is shorter than the hit

#ifdef extendHspFromLeft
	s2 = seq1->v + pos1;			// (past right end of hit)
	if (rightStop < s2)
		length -= s2 - rightStop;
#endif // extendHspFromLeft

#ifdef debugDiag
	if (debugThisDiag)
		dump_extended_match (stdout, seq1, seq2, diag,
		                     p1, p2, p3, p4, rightStop, s1);
#endif

	similarity = leftScore + rightScore;

	//////////
	// (for debugging only)
	// determine if some subrange of the HSP outscores the whole HSP
	//
	// We use the algorithm described in Bentley's "Programming Pearls" ("A
	// scanning Algorithm", section 8.4, page 81 in the second edition).
	//
	// bestSubrangeScore == Bentley's maxSoFar
	// subrangeScore     == Bentley's maxEndingHere
	//////////

#ifdef snoopHspSubrange

	{
	score subrangeScore;
	u8*   currLeft, *subLeft, *subRight;

	s1 = leftStart;
	s2 = seq2->v + diagToPos2 (diag, leftStart - seq1->v);
	currLeft = s1;
	subLeft = subRight = s1;

	subrangeScore = bestSubrangeScore = 0;
	runScore = 0;
	while (s1 < rightStop)
		{
		runScore      = runScore      + scoring->sub[*s1][*s2];
		subrangeScore = subrangeScore + scoring->sub[*(s1++)][*(s2++)];
		if (subrangeScore < 0)
			{ subrangeScore = 0;  currLeft = s1; }
		if (subrangeScore > bestSubrangeScore)
			{
			bestSubrangeScore = subrangeScore;
			subLeft  = currLeft;
			subRight = s1;
			}

		if (debugThisDiag)
			printf (unsposFmt ":"
			        " %c%c " scoreFmtSimple
			        " " scoreFmtSimple
			        " " scoreFmtSimple " " unsposFmt
			        " " scoreFmtSimple " " unsposFmt ".." unsposFmt
			        "\n",
			        (s1-1) - seq1->v,
			        s1[-1], s2[-1], scoring->sub[s1[-1]][s2[-1]],
			        runScore,
			        subrangeScore, currLeft - seq1->v,
			        bestSubrangeScore, subLeft - seq1->v, subRight - seq1->v);
		} 

	hspPos1 = hspPos2 = hspLen = subPos1 = subPos2 = subLen = 0;

	if (bestSubrangeScore > similarity)
		{
		hspPos1 = leftStart - seq1->v;
		hspPos2 = diagToPos2 (diag, leftStart - seq1->v);
		hspLen  = rightStop - leftStart;
		subPos1 = subLeft - seq1->v;
		subPos2 = diagToPos2 (diag, subLeft - seq1->v);
		subLen  = subRight - subLeft;
		seed_search_count_stat (suboptimalHsp);
		}
	}

#endif // snoopHspSubrange

	//////////
	// record the extent of HSP search on this diagonal
	//////////

	// record the extent

	extent = (unspos) (((sgnpos) (rightBlock-seq1->v)) - diag);
	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  setting    diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

#ifdef debugDiag
	if (debugThisDiag)
		{
		printf ("gfex: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt,
				pair_diagonal_as_text(pos1,pos2), pos1, pos2,
				hDiag, diagEnd[hDiag]);
		printf (" (" unsposSlashFmt " -> " unsposSlashFmt ")\n",
				(unspos) (rightStop-seq1->v),  (unspos) (rightStop-seq1->v-diag),
				(unspos) (rightBlock-seq1->v), diagEnd[hDiag]);
		}
#endif

	snoopXDrop_Score;

#ifdef collect_stats
	seed_search_add_stat (bpExtended, rightBlock-p1);
#endif

	//////////
	// update length of hit
	//////////

	pos1   = (unspos) (rightStop - seq1->v);
	pos2   = (unspos) (((sgnpos) pos1) - diag);
	length = (unspos) (rightStop - leftStart);

	//////////
	// if the extended hit's score is acceptable, but not very high, adjust
	// the score downward, based on the entropy of the sequences in the match;
	// note that we only adjust positive scores (since hspZeroThreshold ==
	// max(0,hspZeroThreshold)), otherwise the entropy adjustment would
	// *increase* the score when entropy is poor.
	//
	// When an adaptive scoring threshold is being used, we can't determine
	// what a reasonable "high enough" threshold is to not bother to perform
	// the entropy reduction, so we perform the reduction on any extended hit
	// that could potentially make the hsp table.
	//
	// $$$ Heuristically, we could still estimate some reasonable "high enough"
	//     .. threshold based on the current low score threshold, how many hits
	//     .. we have accepted/rejected so far, and how much room is left in
	//     .. the table.  Simpler schemes may also work well in practice.  This
	//     .. should be considered if the entropy calculation ends up being a
	//     .. significant time factor.
	//////////

	// decide whether to adjust the score

	if (!hp->entropicHsp)
		adjustScore = false;
	else if (hp->hspThreshold.t == 'S')	// (fixed score threshold)
		adjustScore = (similarity >=   hp->hspZeroThreshold)
	               && (similarity <= 3*hp->hspThreshold.s);
	else if (similarity <= 0)			// (adaptive score threshold, negative)
		adjustScore = false;
	else								// (adaptive score threshold, positive)
		{
		segtable* anchors = *(hp->anchors);
		adjustScore = (anchors->len > 0)
		           && (similarity >= anchors->lowScore);
		}

	// adjust it

	if (adjustScore)
		{
		double q = entropy (hp->seq1->v + pos1 - length,
		                    hp->seq2->v + pos2 - length,
		                    length);

		score rawS = similarity;
		similarity *= q;
		if ((similarity < hp->hspThreshold.s) && (hp->reportEntropy))
			fprintf(stderr, "hit of score " scoreFmtSimple
			                " at " unsposSlashFmt "#" unsposFmt " (diag " sgnposFmt " had block at " unsposFmt ")"
			                " fails entropy filter (%f)\n",
							rawS,
							pos1-length, pos2-length, length,
							diag, oldDiagEnd,
							q);
#ifdef snoopEntropy
		else
			fprintf(stderr, "hit of score " scoreFmtSimple
			                " at " unsposSlashFmt "#" unsposFmt " (diag " sgnposFmt " had block at " unsposFmt ")"
			                " passes entropy filter (%f)\n",
							rawS,
							pos1-length, pos2-length, length,
							diag, oldDiagEnd,
							q);
#endif // snoopEntropy

#ifdef snoopHspSubrange
		bestSubrangeScore *= q;
#endif // snoopHspSubrange
		}

	//////////
	// decide whether or not this extended seed hit is an hsp.
	//////////

	dbg_timing_count_stat (ungappedExtensions);

	// if it doesn't score high enough, discard it

	if ((hp->hspThreshold.t == 'S')				// (fixed score threshold)
	 && (similarity < hp->hspThreshold.s))
		{
		seed_search_count_stat (lowScoringHsps);
#ifdef snoopHspSubrange
		if (bestSubrangeScore >= hp->hspThreshold.s)
			{
			fprintf (stderr, "WARNING: discarded HSP " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple
							 " but subrange " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple "\n",
			                 hspPos1, hspPos2, hspLen, similarity,
			                 subPos1, subPos2, subLen, bestSubrangeScore);

			seed_search_count_stat (suboptimalHspB);
			}
		else if ((seed_search_dbgSubrangeHsps)
		      && (bestSubrangeScore > similarity))
			fprintf (stderr, "INFO: HSP " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple
							 " but subrange " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple "\n",
							 hspPos1, hspPos2, hspLen, similarity,
							 subPos1, subPos2, subLen, bestSubrangeScore);
#endif // snoopHspSubrange
		return noScore;
		}

#ifdef snoopHspSubrange
	if ((seed_search_dbgSubrangeHsps)
	 && (bestSubrangeScore > similarity))
		fprintf (stderr, "INFO: HSP " unsposSlashFmt "#" unsposFmt
						 " scores " scoreFmtSimple
						 " but subrange " unsposSlashFmt "#" unsposFmt
						 " scores " scoreFmtSimple "\n",
						 hspPos1, hspPos2, hspLen, similarity,
						 subPos1, subPos2, subLen, bestSubrangeScore);
#endif // snoopHspSubrange

	// it's a keeper

	*_pos1   = pos1;
	*_pos2   = pos2;
	*_length = length;

	if (hp->anchors != NULL)
		(*(hp->anchors))->haveScores = true;

	seed_search_count_stat (hsps);
	dbg_timing_count_stat  (hsps);

	return similarity;
	}

//----------
//
// match_extend_seed_hit--
//	Perform exact match extension on a seed hit, and discard those that aren't
//	long enough.
//
//	Exact match extension extends the hit/match in each direction, along the
//	diagonal, as long as it encounters matches.
//
//	Short extensions are 'discarded' by the caller.  This routine only makes the
//	decision.
//
// Arguments:
//	hitprocinfo*	hp:		Pointer to record containing (among other things)
//							.. the extension controls and filtering criteria.
//	unspos*			pos1:	The hit position in sequence 1.  (see note 1 below)
//	unspos*			pos2:	The hit position in sequence 2.  (see note 1 below)
//	unspos*			length:	The length of the hit.           (see note 1 below)
//
// Returns:
//	The "score" of the extended match;  this is actually the number of matching
//	bases (i.e. the length of the match).  If the match should be filtered
//	(discarded), noScore is returned.
//
//----------
//
// Notes:
//
// (1)	Upon return, the values of pos1, pos2, and length will have been changed
//		to reflect the extended hit *IF* the match is to be kept (i.e. if the
//		return value is not noScore).
//
// (2)	We expect that the hits will arrive in increasing order on sequence 2.  
//
// (3)	It may appear that we (incorrectly) assume that the entirety of
//		sequences 1 and 2 are fair game for the HSP.  However, we halt
//		processing when we encounter a NUL character, which indicates (1) the
//		end of a partition, either in a partitioned sequence or X-separated
//		sequence, or (2) the end of a chore.  Further, positional seed filtering
//		(performed by filter_seed_hit_by_pos) prevents us from being called for
//		hits outside the range of a chore.
//
// (4)	We assume hp->hspThreshold.t is 'S', but we treat it as the required
//		length of the match.
//
// (5)	Any non-ACGT is considered to be a mismatch (even if both sequences
//		have the same letter).
//
// (6)	Though hits arrive in increasing order on sequence 2, it is possible
//		that diagEnd[hDiag] > pos2.  This happens when a previous seed hit on a
//		hash-equivalent diagonal was extended (diagEnd records the righthand
//		limit of that extension).  In such case, left-extension of this seed
//		hit is prematurely halted (at pos2-length).  Right-extension is not
//		affected.
//
//----------

static score match_extend_seed_hit
   (hitprocinfo*	hp,
	unspos*			_pos1,
	unspos*			_pos2,
	unspos*			_length)
	{
	unspos			pos1   = *_pos1;
	unspos			pos2   = *_pos2;
	unspos			length = *_length;
	seq*			seq1   = hp->seq1;
	seq*			seq2   = hp->seq2;
	sgnpos			diag, block2;
	unspos			oldDiagEnd, extent;
	u32				hDiag;
	u8*				s1, *s2, *stop, *left, *right;
	u8				nuc1,  nuc2;
	s8				bits1, bits2;
#ifdef snoopDiagHash
	unspos			start2 = pos2 - length;
#endif // snoopDiagHash

	//////////
	// get ready to extend the hit
	//////////

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));

	if (debugThisDiag)
		{
		printf ("gfex: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, diagEnd[hDiag]);
		}
#endif

	//////////
	// validate that the hit is an exact match
	//////////

	s1   = seq1->v + pos1;
	s2   = seq2->v + pos2;
	stop = s1 - length;

	while (s1 > stop)
		{
		bits1 = hp->charToBits[*(--s1)];
		bits2 = hp->charToBits[*(--s2)];

		if ((bits1 != bits2)
		 || (bits1 < 0)				// (negative => not ACGTacgt)
		 || (bits2 < 0))
			{
			extent = s2 - seq2->v;  // (position of rightmost mismatch)
			goto hit_isnt_a_match;
			}
		}

	//////////
	// extend to the left
	//////////

	s1 = seq1->v + pos1 - length;	// start at start of hit in both seq1 and
	s2 = seq2->v + pos2 - length;	// .. seq2;  will pre-decrement before
									// .. reads, so first bp read are the ones
									// .. immediately in front of the hit

	// determine stop location;  this is at the start of sequence 1, except
	// that if this diagonal ends (or is blocked) earlier in sequence 2, we
	// stop there
	// (see note 3;  instead of zero, we should use subsequence's start)

	if (unblockedLeftExtension) oldDiagEnd = 0;
	                       else oldDiagEnd = diagEnd[hDiag];
	block2 = (sgnpos) oldDiagEnd;
	if (block2 + diag > 0) stop = seq1->v + block2 + diag;
	                  else stop = seq1->v;

	// extend

	if (s1 < stop)
		{
		s1--;	// if the new hit is left of the previous block (as can happen
		s2--;	// .. when called by process_for_recoverable_hit), the normal
				// .. loop (below) will fail to step, and so won't stop on a
				// .. mismatch;  so in this case we need to push the position
				// .. to one to the left of the start of the match
		}
	else
		{
		while (s1 >= stop)
			{
			if (s1 == stop)				// (this test is necessary since we
				{ s1--;  s2--;  break; }//  .. don't have a zero-terminator at
										//  .. the start of the sequence)
			nuc1  = *(--s1);
			bits1 = hp->charToBits[nuc1];
			nuc2  = *(--s2);
			bits2 = hp->charToBits[nuc2];

			if ((nuc1 == 0)				// (NUL => end of partition or chore)
			 || (nuc2 == 0)
			 || (bits1 != bits2)
			 || (bits1 < 0)				// (negative => not ACGTacgt)
			 || (bits2 < 0))
				break;
			}
		}

	left = s1;	// the first mismatch, or just to the left of the stop

	//////////
	// extend to the right
	//////////

	s1 = seq1->v + pos1 - 1;	// start at end of hit in both seq1 and seq2;
	s2 = seq2->v + pos2 - 1;	// .. will pre-increment before reads, so first
								// .. bp read are the ones immediately after the
								// .. hit

	// determine stop location;  this is at the end of sequence 1, except
	// that if this diagonal ends earlier in sequence 2, we stop there
	// (see note 3;  instead of sequence's end, we should use subsequence's end)

	block2 = (sgnpos) seq2->len;
	if ((sgnpos) seq1->len <= block2 + diag) stop = seq1->v + seq1->len;
	                                    else stop = seq1->v + block2 + diag;

	// extend

	while (s1 < stop)
		{
		nuc1  = *(++s1);
		bits1 = hp->charToBits[nuc1];
		nuc2  = *(++s2);
		bits2 = hp->charToBits[nuc2];

		if ((nuc1 == 0)				// (NUL => end of partition or chore)
		 || (nuc2 == 0)
		 || (bits1 != bits2)
		 || (bits1 < 0)				// (negative => not ACGTacgt)
		 || (bits2 < 0))
			break;
		}

	right = s1;	// the first mismatch, or at the stop

	//////////
	// record the extent of the search on this diagonal
	//////////

	// record the extent

	extent = (unspos) (((sgnpos) (right-seq1->v)) - diag);
	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  m setting  diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

	//////////
	// update length of hit
	//////////

	pos1   = (unspos) (right - seq1->v);
	pos2   = (unspos) (((sgnpos) pos1) - diag);
	length = (unspos) (right - (left + 1));

	//////////
	// decide whether or not this extended seed hit is long enough.
	//////////

	dbg_timing_count_stat (ungappedExtensions);

	// if it isn't long enough, discard it (see note 4)

	if (length < (unsigned) hp->hspThreshold.s)
		{
		seed_search_count_stat (lowScoringHsps);
		return noScore;
		}

	// it's a keeper

	*_pos1   = pos1;
	*_pos2   = pos2;
	*_length = length;

	seed_search_count_stat (hsps);
	dbg_timing_count_stat  (hsps);

	return (score) length;

	//////////
	// special exit for the case of the hit not being an exact match
	//////////

hit_isnt_a_match:

	// record the extent of the search on this diagonal

	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  nm setting diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

	seed_search_count_stat (lowScoringHsps);
	dbg_timing_count_stat  (ungappedExtensions);
	return noScore;
	}

//----------
//
// mismatch_extend_seed_hit--
//	Perform M-mismatch extension on a seed hit, and discard those that aren't
//	long enough.
//
//	Conceptually, M-mismatch extension extends the hit/match in each direction,
//	along the diagonal, as long as it encounters fewer than M mismatches.  In
//	actuality, it considers more then M mismatches because it scans both left
//	and right;  then picks the longest of several M-mismatch segments.
//
//	Short extensions are 'discarded' by the caller.  This routine only makes the
//	decision.
//
// Arguments:
//	hitprocinfo*	hp:		Pointer to record containing (among other things)
//							.. the extension controls and filtering criteria.
//	unspos*			pos1:	The hit position in sequence 1.  (see note 1 below)
//	unspos*			pos2:	The hit position in sequence 2.  (see note 1 below)
//	unspos*			length:	The length of the hit.           (see note 1 below)
//
// Returns:
//	The "score" of the extended match;  this is actually the length of the
// 	match (i.e. the number of matching bases plus the number of mismatches).
//	If the match should be filtered (discarded), noScore is returned.
//
//----------
//
// Notes:
//
// (1)	Upon return, the values of pos1, pos2, and length will have been changed
//		to reflect the extended hit *IF* the match is to be kept (i.e. if the
//		return value is not noScore).
//
// (2)	We expect that the hits will arrive in increasing order on sequence 2.  
//
// (3)	It may appear that we (incorrectly) assume that the entirety of
//		sequences 1 and 2 are fair game for the HSP.  However, we halt
//		processing when we encounter a NUL character, which indicates (1) the
//		end of a partition, either in a partitioned sequence or X-separated
//		sequence, or (2) the end of a chore.  Further, positional seed filtering
//		(performed by filter_seed_hit_by_pos) prevents us from being called for
//		hits outside the range of a chore.
//
// (4)	We assume hp->hspThreshold.t is 'S', but we treat it as the required
//		length of the match.
//
// (5)	Any non-ACGT is considered to be a mismatch (even if both sequences
//		have the same letter).
//
// (6)	Though hits arrive in increasing order on sequence 2, it is possible
//		that diagEnd[hDiag] > pos2.  This happens when a previous seed hit on a
//		hash-equivalent diagonal was extended (diagEnd records the righthand
//		limit of that extension).  In such case, left-extension of this seed
//		hit is prematurely halted (at pos2-length).  Right-extension is not
//		affected.
//
//----------
//
// Algorithm:
//
// The figure below shows an example of a 5 mismatch extension (M=5).  The
// asterisks represent a 19 bp seed hit.  Below that is the sequence of match
// (-) or mismatch (o) along the diagonal.  The seed contains two mismatches
// (E=2).
//
//	                         [--- seed hit ----]
//	                         *******************
//	sequence: o-----o-----oo---------o------o-----------o--------------o--o---o
//	(len=41)   [----------------- 5-mm ----------------]
//	(len=50)         [--------------------- 5-mm ---------------------]
//	(len=47)               [-------------------- 5-mm -------------------]
//	(len=50)                [--------------------- 5-mm ---------------------]
//
// We scan left until we find the 4th mismatch (M+1-E = 4).  Any of these can
// determine the start of a 5-mismatch segment containing the seed.  Note that
// the segment actually starts at the first base *after* that mismatch.  There
// are 4 possible starting points.  We scan right to find the end points for
// each of these 4 intervals, and choose the longest.  In the case of ties we
// prefer the one further to the left.
//
// The algorithm is complicated by the fact that we may hit either endpoint
// before seeing the requisite number of mismatches.  The example below shows
// a case where we fall 2 mismatches short during the left scan.  In this case
// we treat the stop point (x) as a mismatch.  Since this leaves us still one
// mismatch short, we skip the first mismatch during right-scanning.
//
//	                stop     [--- seed hit ----]
//	                  |      *******************
//	sequence:         x---oo---------o------o-----------o--------------o--o---o
//	(len=48)           [-------------------- 5-mm --------------------]
//	(len=47)               [-------------------- 5-mm -------------------]
//	(len=50)                [--------------------- 5-mm ---------------------]
//
// The extent (which limits subsequent processing on teh same diagonal) is set
// to the first of these that is true
//	(1) if seed hit had more than E mismatches, to the E+1st rightmost mismatch
//	(2) if the segment is long enough, just beyond the mismatch at the right
//	    end of the segment
//	(3) if seed hit had any mismatches, to the leftmost mismatch
//	(4) otherwise, to the leftmost mismatch to the right of the seed
//
// The justifciation for (3) is that when we reject all possible intervals for
// this seed hit, we have considered all intervals starting with anything to the
// left of the leftmost mismatch in the seed hit.  So we would like to consider
// intervals beginning from that point as part of a later seed hit. The
// justification for (1) and (4) is similar.
//
// The justification for (2) is that we don't want to 'report' segments that
// overlap.  Note that in a long segment of high identity, there will be many
// M-mismatch intervals that are long enough, all them overlapping.  Criteria
// (2) has the unfortunate side effect of often not reporting the longest
// acceptable interval in such a run.  However, when the segments are used as
// anchors for gapped extension, this side effect is unlikely to matter.  And
// the time savings gained from setting the extent as far to the right as
// possible, and from not reporting overlapping segments, can be substantial.
//
//----------

#ifndef debugMismatchExtend
#define debugMismatchExtend_1  ;
#define debugMismatchExtend_2  ;
#define debugMismatchExtend_3A ;
#define debugMismatchExtend_3B ;
#define debugMismatchExtend_3C ;
#endif // not debugMismatchExtend

#ifdef debugMismatchExtend

#ifdef debugDiag
static int debugMmExtendDiag;
#else
#define debugMmExtendDiag true
#endif

#define debugMismatchExtend_1                                                \
	if (debugMmExtendDiag)                                                   \
		{                                                                    \
		fprintf (stderr, "mmex: hit at " unsposSlashFmt " %d mm\n",          \
		                 pos1-length, pos2-length, E);                       \
		}

#define debugMismatchExtend_2                                                \
	if (debugMmExtendDiag)                                                   \
		{                                                                    \
		u8** mm;                                                             \
		fprintf (stderr, "mmex: lefties at");                                \
		for (mm=mmScan ; mm<mmStop ; mm++)                                   \
			{                                                                \
			if (mm == mmScan) fprintf (stderr, " ");                         \
			else              fprintf (stderr, ", ");                        \
			fprintf (stderr, unsposFmt, *mm - seq1->v);                      \
			}                                                                \
		if (mmShortfall > 0)                                                 \
			fprintf (stderr, " shortfall=%d", mmShortfall);                  \
		fprintf (stderr, "\n");                                              \
		}

#define debugMismatchExtend_3A                                               \
	if (debugMmExtendDiag)                                                   \
		{                                                                    \
		unspos p1 = (*mmScan)+1 - seq1->v;                                   \
		unspos p2 = diagToPos2 (diag,p1);                                    \
		fprintf (stderr, "mmex: %dmm interval"                               \
		                 " at " unsposSlashFmt " length %d (shorty)\n",      \
		                 M-mmShortfall, p1, p2, (s1-1)-*mmScan);             \
		}

#define debugMismatchExtend_3B                                               \
	if (debugMmExtendDiag)                                                   \
		{                                                                    \
		unspos p1 = (*mmScan)+1 - seq1->v;                                   \
		unspos p2 = diagToPos2 (diag,p1);                                    \
		fprintf (stderr, "mmex: %dmm interval"                               \
		                 " at " unsposSlashFmt " length %d\n",               \
		                 M, p1, p2, thisLength-1);                           \
		}

#define debugMismatchExtend_3C                                               \
	if (debugMmExtendDiag)                                                   \
		{                                                                    \
		unspos p1 = (*mmScan)+1 - seq1->v;                                   \
		unspos p2 = diagToPos2 (diag,p1);                                    \
		fprintf (stderr, "mmex: %dmm interval"                               \
		                 " at " unsposSlashFmt " length %d (right stop)\n",  \
		                 M, p1, p2, thisLength-1);                           \
		}


#endif // debugMismatchExtend


// mismatch_extend_seed_hit--

static score mismatch_extend_seed_hit
   (hitprocinfo*	hp,
	unspos*			_pos1,
	unspos*			_pos2,
	unspos*			_length)
	{
	unspos			pos1   = *_pos1;
	unspos			pos2   = *_pos2;
	unspos			length = *_length;
	seq*			seq1   = hp->seq1;
	seq*			seq2   = hp->seq2;
	sgnpos			diag, block2;
	unspos			oldDiagEnd, extent;
	u32				hDiag;
	u8*				s1, *s2, *stop, *left, *right;
	u8				nuc1,  nuc2;
	s8				bits1, bits2;
	int				M = hp->gfExtend;			// (max allowed mismatches)
	int				E;							// (number of mms in seed)
	u8*				mmLoc[gfexMismatch_max+1];	// (mm loc'ns as starters)
	u8**			mmScan, **mmStop;
	int				mmShortfall;
	unspos			thisLength, bestLength;
#ifdef snoopDiagHash
	unspos			start2 = pos2 - length;
#endif // snoopDiagHash

	//////////
	// get ready to extend the hit
	//////////

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

#ifdef debugDiag
	debugThisDiag = (hDiag == hashedDiag(debugDiag,0));
#ifdef debugMismatchExtend
	debugMmExtendDiag = debugThisDiag;
#endif

	if (debugThisDiag)
		{
		printf ("mmex: (diag %9s", pair_diagonal_as_text(pos1,pos2));
		printf ("|%9s|%04X) " unsposSlashFmt " end was " unsposFmt "\n",
		        diagonal_as_text(diagActual[hDiag]), hDiag,
		        pos1, pos2, diagEnd[hDiag]);
		}
#endif

	//////////
	// count the number of mismatches in the hit
	//////////

	s1   = seq1->v + pos1;
	s2   = seq2->v + pos2;
	stop = s1 - length;

	E = 0;
	extent = hashInactiveEnd;		// (this will remain unchanged through
									//  .. the following loop iff there are
									//  .. no mismatches in the hit)

	while (s1 > stop)
		{
		bits1 = hp->charToBits[*(--s1)];
		bits2 = hp->charToBits[*(--s2)];

		if ((bits1 != bits2)
		 || (bits1 < 0)				// (negative => not ACGTacgt)
		 || (bits2 < 0))
			{
			extent = s2 - seq2->v;  // (leftmost interesting mismatch in hit)
			if (++E > M)			// seed contains too many mismatches
				goto hit_isnt_a_match;
			}
		}

	debugMismatchExtend_1;

	//////////
	// extend left until the M+1-Eth mismatch, saving positions in an array;
	// note that we might not find that many mismatches, since we may hit the
	// stop first
	//////////

	s1 = seq1->v + pos1 - length;	// start at start of hit in both seq1 and
	s2 = seq2->v + pos2 - length;	// .. seq2;  will pre-decrement before
									// .. reads, so first bp read are the ones
									// .. immediately in front of the hit

	// determine stop location;  this is at the start of sequence 1, except
	// that if this diagonal ends (or is blocked) earlier in sequence 2, we
	// stop there
	// (see note 3;  instead of zero, we should use subsequence's start)

	if (unblockedLeftExtension) oldDiagEnd = 0;
	                       else oldDiagEnd = diagEnd[hDiag];
	block2 = (sgnpos) oldDiagEnd;
	if (block2 + diag > 0) stop = seq1->v + block2 + diag;
	                  else stop = seq1->v;

	// set up mmScan to the first location past the end of the mmLoc array;
	// we view the array as being of size M+1-E, since this is the number of
	// mismatches left of the seed that we have any interest in;  we will
	// pre-decrement this pointer as we collect mismatches into the array;
	// entries in the array will point to the position of a mismatch (or a fake
	// mismatch one left of the stop);  note that since M+1-E > 0, we will
	// always be looking for at least 1 mismatch

	mmScan = mmLoc + M+1-E;
	mmStop = mmScan; // (this is used in the right-scanning stage)

	// extend

	if (s1 < stop)
		{
		s1--;	// if the new hit is left of the previous block (as can happen
		s2--;	// .. when called by process_for_recoverable_hit), the normal
				// .. loop (below) will fail to step, and so won't stop on a
				// .. mismatch;  so in this case we need to push the position
				// .. to one to the left of the start of the match
		}
	else
		{
		while (s1 >= stop)
			{
			if (s1 == stop)				// (this test is necessary since we
				{ s1--;  s2--;  break; }//  .. don't have a zero-terminator at
										//  .. the start of the sequence)

			nuc1  = *(--s1);
			bits1 = hp->charToBits[nuc1];
			nuc2  = *(--s2);
			bits2 = hp->charToBits[nuc2];

			if ((nuc1 == 0)				// (NUL => end of partition or chore)
			 || (nuc2 == 0))
				break;

			if ((bits1 != bits2)
			 || (bits1 < 0)				// (negative => not ACGTacgt)
			 || (bits2 < 0))
				{
				*(--mmScan) = s1;		// save this as left endpoint
				if (mmScan == mmLoc)	// if we have enough, quit scanning
					break;
				}
			}
		}

	// if we did not get enough mismatches, add the point at which we stopped;
	// also record the number of mismatches that we *didn't* find, since we will
	// need to skip that many during the right-scanning stage;  note that this
	// guarantees that the array of collected interval starts is never empty

	if (mmScan > mmLoc)
		*(--mmScan) = s1;

	mmShortfall = mmScan - mmLoc;
	debugMismatchExtend_2;

	//////////
	// extend to the right, finding an ending mismatch for each of our
	// intervals
	//////////

	s1 = seq1->v + pos1 - 1;	// start at end of hit in both seq1 and seq2;
	s2 = seq2->v + pos2 - 1;	// .. will pre-increment before reads, so first
								// .. bp read are the ones immediately after the
								// .. hit

	// determine stop location;  this is at the end of sequence 1, except
	// that if this diagonal ends earlier in sequence 2, we stop there
	// (see note 3;  instead of sequence's end, we should use subsequence's end)

	block2 = (sgnpos) seq2->len;
	if ((sgnpos) seq1->len <= block2 + diag) stop = seq1->v + seq1->len;
	                                    else stop = seq1->v + block2 + diag;

	// extend

	bestLength = 0;
	left = right = NULL;

	while (s1 < stop)
		{
		nuc1  = *(++s1);
		bits1 = hp->charToBits[nuc1];
		nuc2  = *(++s2);
		bits2 = hp->charToBits[nuc2];

		if ((nuc1 == 0)					// (NUL => end of partition or chore)
		 || (nuc2 == 0))
			break;

		if ((bits1 != bits2)
		 || (bits1 < 0)					// (negative => not ACGTacgt)
		 || (bits2 < 0))
			{
			if (extent == hashInactiveEnd)
				extent = s2 - seq2->v;
			if (mmShortfall > 0)
				{
				debugMismatchExtend_3A;
				mmShortfall--;
				continue;
				}
			thisLength = s1 - *mmScan;	
			debugMismatchExtend_3B;
			if (thisLength > bestLength)// if this is the best interval so far
				{						// .. save the endpoints
				bestLength = thisLength;
				left       = *mmScan;
				right      = s1;
				}
			if (++mmScan == mmStop)		// if we have enough, quit scanning
				break;
			}
		}

	// if we did not get enough mismatches, treat the point at which we stopped
	// as another endpoint;  note that we don't honor mmShortfall here, because
	// the interval at this endpoint has less than M mismatches and if it is
	// long enough it will be an acceptible extension

	if (mmScan < mmStop) 
		{
		if (extent == hashInactiveEnd)
			extent = s2 - seq2->v;
		thisLength = s1 - *mmScan;	
		debugMismatchExtend_3C;
		if (thisLength > bestLength)	// if this is the best interval so far
			{							// .. save the endpoints
			left  = *mmScan;
			right = s1;
			}
		}

	if (left == NULL)
		suicide ("internal error (in mismatch_extend_seed_hit) found no interval");

	//////////
	// update length of hit
	//////////

	pos1   = (unspos) (right - seq1->v);
	pos2   = (unspos) (((sgnpos) pos1) - diag);
	length = (unspos) (right - (left + 1));

	//////////
	// record the extent of the search on this diagonal
	//////////

	if (length >= (unsigned) hp->hspThreshold.s)
		extent = (unspos) (((sgnpos) pos1+1) - diag);

	// record the extent

	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  m setting  diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

	//////////
	// decide whether or not this extended seed hit is long enough.
	//////////

	dbg_timing_count_stat (ungappedExtensions);

	// if it isn't long enough, discard it (see note 4)

	if (length < (unsigned) hp->hspThreshold.s)
		{
		seed_search_count_stat (lowScoringHsps);
		return noScore;
		}

	// it's a keeper

	*_pos1   = pos1;
	*_pos2   = pos2;
	*_length = length;

	seed_search_count_stat (hsps);
	dbg_timing_count_stat  (hsps);

	return (score) length;

	//////////
	// special exit for the case of the hit not being an exact match
	//////////

hit_isnt_a_match:

	// record the extent of the search on this diagonal

	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  nm setting diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

	seed_search_count_stat (lowScoringHsps);
	dbg_timing_count_stat  (ungappedExtensions);
	return noScore;
	}

//----------
//
// warn_for_search_limit--
//	Tell the user that this query exceeded the limit for HSPs.
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

static void warn_for_search_limit
   (void)
	{
	static int firstReport = true;
	char* name2 ;

	seed_search_dbgSearchLimitExceeded++;
	if (reportSearchLimit == 0) return;

	name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
	fprintf (stderr, "WARNING. Query \"%s\" contains more than %s HSPs.\n",
					 name2, commatize(reportSearchLimit));

	if (firstReport)
		{
		fprintf (stderr, "All HSPs for this query are discarded and the query is not processed further.\n");
		firstReport = false;
		}
	}

//----------
//
// discovery_probability--
//	Compute the probability that a particular HSP would be discovered by a given
//	(seed,step) search strategy, if it was equally likely to have occured at any
//	sequence position.
//
// For any step size Z greater than 1, our search process will miss some HSPs.
// For example, if an HSP contains only one seed hit, and that seed hit occurs
// in sequence 2 at an odd position, if Z=2 we will miss that seed hit (and thus
// the HSP).  More generally, if a seed hit occurs at position X, we will only
// find it if X == 0 modulo Z.
//
// This routine scans an HSP for the positions of all the seed hits it contains,
// and counts how many positional shifts would put at least one seed hit on a
// multiple of Z.  There are Z different positional shifts, so the probability
// that the HSP would have been discovered is C/Z, where C is that count.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The hit position in sequence 1.  This is the position
//					.. following the end of the hit.
//	seq*	seq2:	The other sequence.
//	unspos	pos2:	The hit position position in sequence 2 (same details as
//					.. pos1).
//	unspos	length:	The length of the alignment.
//	seed*	hitSeed: Seeding strategy for the hits that found this match.
//	u32		step:	Positional step size in the search for those hits.
//
// Returns:
//	(nothing)
//
//----------

float discovery_probability
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	seed*	hitSeed,
	u32		step)
	{
	u8*		aStart = seq1->v + pos1 - length;
	u8*		aStop  = seq1->v + pos1;
	u8*		bStart = seq2->v + pos2 - length;
	u8*		a, *b;
	int		aa, bb;
	u64		aUnpacked, bUnpacked;
	u32		aPacked,   bPacked;
	int		len;
	u32		*flip, flipBits, trans;
	u64		diffBits, transBits;
	u32		i;
	int		foundCount;

	// allocate a scratch array

	if ((foldedSize > step) && (foldedHits != NULL))
		{ free_if_valid ("folded hits", foldedHits);  foldedHits = NULL; }

	if (foldedHits == NULL)
		{
		foldedHits = (u8*) malloc_or_die ("", step);
		foldedSize = step;
		}

	// build the transition bits mask

	flipBits = 0;
	for (flip=hitSeed->transFlips ; (*flip)!=0 ; flip++)
		flipBits += *flip;

	transBits = seed_unpack (hitSeed, flipBits, NULL);

	// scan the alignment, checking for seed matches

	foundCount = 0;
	for (i=0 ; i<step ; i++)
		foldedHits[i] = 0;

	aUnpacked = bUnpacked = 0;

	for (a=aStart,b=bStart ; a<aStop ; )
		{
		// collect the first seedLength-1 columns

	empty:
		for (len=1 ; (len<hitSeed->length)&&(a<aStop) ; len++)
			{
			aa = nuc_to_bits[*(a++)];  if (aa < 0) goto empty;
			bb = nuc_to_bits[*(b++)];  if (bb < 0) goto empty;
			aUnpacked = (aUnpacked << 2) | aa;
			bUnpacked = (bUnpacked << 2) | bb;
			}

		// process each word of seedLength columns

		while (a<aStop)
			{
			aa = nuc_to_bits[*(a++)];  if (aa < 0) goto empty;
			bb = nuc_to_bits[*(b++)];  if (bb < 0) goto empty;
			aUnpacked = (aUnpacked << 2) | aa;
			bUnpacked = (bUnpacked << 2) | bb;
			aPacked = apply_seed (hitSeed, aUnpacked);
			bPacked = apply_seed (hitSeed, bUnpacked);

			if (aPacked == bPacked) goto got_a_hit;
			if (hitSeed->withTrans == 0) continue;

			diffBits = aUnpacked ^ bUnpacked;
			trans = (diffBits<<1) & transBits;               // (1=>transversion)
			if (trans != 0) continue;
			trans = (diffBits & ~(diffBits<<1)) & transBits; // (1=>transition)
			if (bit_count(trans) > hitSeed->withTrans) continue;

		got_a_hit:
			i = (a-(aStart+hitSeed->length)) % step;
			if (foldedHits[i] == 0)
				{ foldedHits[i] = 1;  foundCount++; }
			}
		}

	return foundCount / (float) step;
	}

//----------
//
// dump_raw_hit--
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	unspos	pos1, pos2:	The hit position in sequences 1 and 2.  This is the
//						.. position following the end of the hit.
//
// Returns:
//	(nothing)
//
//----------

static void dump_raw_hit
   (FILE*	f,
	unspos	pos1,
	unspos	pos2)
	{
	int		isRev1 = ((seq1->revCompFlags & rcf_rev) != 0);
	int		isRev2 = ((seq2->revCompFlags & rcf_rev) != 0);
	u32		seedLength, len1;

	seedLength = (unsigned) hitSeed->length;
	len1       = seedLength-1;

	fprintf      (f, "raw seed hit " unsposSlashCFmt " ",
	                 pos1-len1, (isRev1)?'-':'+',
	                 pos2-len1, (isRev2)?'-':'+');
	print_prefix (f, (char*) seq1->v + pos1-seedLength, (int) seedLength);
	fprintf      (f, "/");
	print_prefix (f, (char*) seq2->v + pos2-seedLength, (int) seedLength);
	fprintf      (f, "\n");
	}

//----------
//
// dump_extended_match--
//	Show an extended match, with its negative scoring flanks.
//
// Example:
//
//		4100: TTGCAAGAAGG ACAT[GGAAGGAA]GA ACGGATCTA
//		      GCTGTTATCAA ACAA[GGAAGGAA]GA CTTCTAGGT
//
//	- GGAAGGAA/GGAAGGAA is the seed match (in this case it was an 8-mer
//	  exact match)
//	- ACAT/ACAA improves the score on the left
//	- GA/GA improves the score on the right
//	- TTGCAAGAAGG/GCTGTTATCAA would drop the score on the left
//	- ACGGATCTA/CTTCTAGGT would drop the score on the right
//	- 4100 is the index of ACGGATCTA in sequence 1 (origin-0)
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	seq1, seq2:	The sequences.
//	sgnpos	diag:		The diagonal the match is on (pos1 - pos2).
//	u8*		p1..p6: 	Pointers into sequence 1, best described by the diagram
//						.. below.
//
//		4100: TTGCAAGAAGG ACAT[GGAAGGAA]GA ACGGATCTA
//		      ^           ^    ^        ^  ^        ^
//		      p1          p2   p3       p4 p5       p6
//
// Returns:
//	(nothing)
//
//----------

#ifdef debugDiag

static void dump_extended_match
   (FILE*	f,
	seq*	seq1,
	seq*	seq2,
	sgnpos	diag,
	u8*		p1,
	u8*		p2,
	u8*		p3,
	u8*		p4,
	u8*		p5,
	u8*		p6)
	{
	u8*		s1, *s2;

	fprintf (f, "\n");
	fprintf (f, "%9u: ", (unspos) (p1-seq1->v));
	for (s1=p1 ; s1<p2 ; s1++) fprintf (f, "%c", *s1);
	fprintf (f, " ");
	for ( ; s1<p3 ; s1++) fprintf (f, "%c", *s1);
	fprintf (f, "[");
	for ( ; s1<p4 ; s1++) fprintf (f, "%c", *s1);
	fprintf (f, "]");
	for ( ; s1<p5 ; s1++) fprintf (f, "%c", *s1);
	fprintf (f, " ");
	for ( ; s1<p6 ; s1++) fprintf (f, "%c", *s1);
	fprintf (f, "\n");

	fprintf (f, "%9u: ", (unspos) (p1-seq1->v - diag));
	s2 = seq2->v + ((sgnpos) (p1-seq1->v)) - diag;
	for (s1=p1 ; s1<p2 ; s1++) fprintf (f, "%c", *(s2++));
	fprintf (f, " ");
	for ( ; s1<p3 ; s1++) fprintf (f, "%c", *(s2++));
	fprintf (f, "[");
	for ( ; s1<p4 ; s1++) fprintf (f, "%c", *(s2++));
	fprintf (f, "]");
	for ( ; s1<p5 ; s1++) fprintf (f, "%c", *(s2++));
	fprintf (f, " ");
	for ( ; s1<p6 ; s1++) fprintf (f, "%c", *(s2++));
	fprintf (f, "\n");
	}

#endif // debugDiag

//----------
//
// diagonal_as_text--
//	Convert the diagonal of a hit to text, as either +diag, -diag, or 0.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	unspos	pos1:	The hit position in sequence 1.
//	unspos	pos2:	The hit position in sequence 2.
//
// Returns:
//	A string containing the text, static data local to this function.
//
//----------

#if ((defined debugDiag) || (defined snoopDiagHash) || (defined snoopXDrop))

static char* pair_diagonal_as_text (unspos pos1, unspos pos2)
	{ return diagonal_as_text(diagNumber(pos1,pos2)); }

static char* diagonal_as_text
   (sgnpos	diag)
	{
	static char s[25];  // (more than enough for 2^64 in decimal with sign)

	if      (diag == 0) sprintf (s, "0");
	else if (diag >  0) sprintf (s, "+" sgnposFmt, diag);
	else                sprintf (s,     sgnposFmt, diag);

	return s;
	}

#endif // debugDiag || snoopXDrop

//----------
//
// seed_search_zero_stats--
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

void seed_search_zero_stats
   (void)
	{
	dbg_timing_set_stat (ungappedExtensions, 0);
	dbg_timing_set_stat (hsps,               0);

#ifdef collect_stats

	// set 'em en masse to zero

	memset (&seedSearchStats, 0, sizeof(seedSearchStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// seed_search_show_stats--
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

void seed_search_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
#ifdef maxHitsPerColumn
	int		hits;
	int		haveAnyHits;
	char	scratch[19];
#endif // maxHitsPerColumn
#endif // collect_stats

	dbg_timing_report_big_stat (ungappedExtensions, "ungapped extensions");
	dbg_timing_report_stat     (hsps,               "HSPs");

#ifdef collect_stats

	if (f == NULL) return;

	fprintf (f, "  allow transition: %s\n", (seedSearchStats.withTrans==0)? "no" :
	                                        (seedSearchStats.withTrans==1)? "yes"
	                                                                      : "two");
	fprintf (f, "    words in seq 2: %s\n", commatize(seedSearchStats.wordsInSequence));
	fprintf (f, "   false seed hits: %s\n", commatize(seedSearchStats.unresolvedSeedHits));
	fprintf (f, "     raw seed hits: %s\n", commatize(seedSearchStats.rawSeedHits));
	if (seedSearchStats.rawSeedHits > 0)
		{
		fprintf (f, "   hash collisions: %s (%.2f%%)\n", commatize(seedSearchStats.hashCollisions),
		                                                 100*seedSearchStats.hashCollisions / (float) seedSearchStats.rawSeedHits);
		fprintf (f, "     hash failures: %s (%.2f%%)\n", commatize(seedSearchStats.hashFailures),
		                                                 100*seedSearchStats.hashFailures / (float) seedSearchStats.rawSeedHits);
		}
	fprintf (f, "       bp extended: %s\n", commatize(seedSearchStats.bpExtended));

#ifndef noSeedHitQueue
	fprintf (f, "     queue scanned: %s (%.1f)\n", commatize(seedSearchStats.queueSeedsScanned),
	                                               seedSearchStats.queueSeedsScanned / (float) seedSearchStats.rawSeedHits);
	fprintf (f, "    queue examined: %s (%.1f)\n", commatize(seedSearchStats.queueSeedsExamined),
	                                               seedSearchStats.queueSeedsExamined / (float) seedSearchStats.rawSeedHits);
	fprintf (f, "     queue blocked: %s (%.1f)\n", commatize(seedSearchStats.queueSeedsBlocked),
	                                               seedSearchStats.queueSeedsBlocked / (float) seedSearchStats.rawSeedHits);
#endif // not noSeedHitQueue
	fprintf (f, "-------------------\n");

	if (seedSearchStats.minMatches >= 0)
		{
		fprintf (f, "     matches >= %2d: %s\n", seedSearchStats.minMatches, commatize(seedSearchStats.notEnoughMatches));
		if (seedSearchStats.maxTransversions >= 0)
			fprintf (f, "transvers'ns <= %2d: %s\n", seedSearchStats.maxTransversions, commatize(seedSearchStats.tooManyTransversions));
		if (seedSearchStats.filterCaresOnly)
			fprintf (f, "     (cares only)\n");
		fprintf (f, "-------------------\n");
		}
	if (seedSearchStats.searchLimit > 0)
		fprintf (f, "      search limit: %s\n", commatize(seedSearchStats.searchLimit));
	if (seedSearchStats.isHspSearch)
		{
		int64 numExt = seedSearchStats.hsps+seedSearchStats.lowScoringHsps;
		fprintf (f, "     GF extensions: %s\n", commatize(numExt));
		fprintf (f, "     HSP wanna-bes: %s\n", commatize(seedSearchStats.lowScoringHsps));
		fprintf (f, "              HSPs: %s\n", commatize(seedSearchStats.hsps));
		if (numExt > 0)
			fprintf (f, "      bp/extension: %s\n", commatize((2*seedSearchStats.bpExtended+numExt)/(2*numExt)));
#ifdef snoopHspSubrange
		fprintf (f, "   suboptimal hsps: %s (%.2f%%)\n", commatize(seedSearchStats.suboptimalHsp),
		                                                 100*seedSearchStats.suboptimalHsp / (float) seedSearchStats.rawSeedHits);
		if (seedSearchStats.hsps + seedSearchStats.suboptimalHspB != 0)
			fprintf (f, "unjustly discarded: %s (%.2f%%)\n", commatize(seedSearchStats.suboptimalHspB),
			                                                 100*seedSearchStats.suboptimalHspB / (float) (seedSearchStats.hsps + seedSearchStats.suboptimalHspB));
		else
			fprintf (f, "unjustly discarded: 0\n");
#endif // snoopHspSubrange
		fprintf (f, "-------------------\n");
		}

#ifdef maxHitsPerColumn
	haveAnyHits = false;
	for (hits=1 ; hits<=maxHitsPerColumn+1 ; hits++)
		{
		if (seedSearchStats.hitsPerColumn[hits] != 0) { haveAnyHits = true;  break; }
		}

	if (haveAnyHits)
		{
		fprintf (f, "(seq 2 words with N raw seed hits)\n");
		for (hits=0 ; hits<=maxHitsPerColumn ; hits++)
			{
			if (seedSearchStats.hitsPerColumn[hits] == 0) continue;
			fprintf (f, "%18d: %s\n", hits, commatize(seedSearchStats.hitsPerColumn[hits]));
			}
		hits = maxHitsPerColumn + 1;
		if (seedSearchStats.hitsPerColumn[hits] != 0)
			{
			sprintf (scratch, "> %d", maxHitsPerColumn);
			fprintf (f, "%18s: %s\n", scratch, commatize(seedSearchStats.hitsPerColumn[hits]));
			fprintf (f, " max hits for word: %s\n", commatize(seedSearchStats.mostHitsInColumn));
			}
		fprintf (f, "-------------------\n");
		}
#endif // maxHitsPerColumn

#endif // collect_stats
	}

void seed_search_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "raw_seed_hits=%" PRId64 "\n",   seedSearchStats.rawSeedHits);
	(*func) (f, "hash_collisions=%" PRId64 "\n", seedSearchStats.hashCollisions);
	(*func) (f, "hash_failures=%" PRId64 "\n",   seedSearchStats.hashFailures);
	(*func) (f, "bp_extended=%" PRId64 "\n",     seedSearchStats.bpExtended);
#endif // collect_stats
	}

#ifdef collect_stats
int64 seed_search_hsps             (void) { return seedSearchStats.hsps; }
int64 seed_search_low_scoring_hsps (void) { return seedSearchStats.lowScoringHsps; }
int64 seed_search_bp_extended      (void) { return seedSearchStats.bpExtended; }
#endif // collect_stats

