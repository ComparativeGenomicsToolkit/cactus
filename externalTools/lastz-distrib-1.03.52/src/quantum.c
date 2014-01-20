//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: quantum.c
//
//----------
//
// quantum--
//	Support for finding "high scoring segment pairs" between a quantum DNA
//	sequence and a DNA sequence.
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
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "seed_search.h"		// seed hit search stuff

#define  quantum_owner			// (make this the owner of its globals)
#include "quantum.h"			// interface to this module

//----------
//
// data structures and types
//
//----------

//#define debugQuantumPartials	// if defined (and if quantum_dbgQuantumBall is
								// .. true), show the failed "partials" in the
								// .. computation of the ball of DNA words close
								// .. to each quantum word

// private globals shared by all the routines under the umbrella of
// quantum_seed_hit_search()

static seq*			seq1;
static postable*	pt;
static seq*			seq2;
static unspos		start;
static unspos		end;
static const s8*	charToBits;
static seed*		hitSeed;
static scoreset*	scoring;
static score		ballScore;
static hitprocessor	processor;
static void*		processorInfo;
static const u8*	bitsToSym;

//----------
//
// prototypes for private functions
//
//----------

static u32 private_quantum_word_hit_search (void);
static u32 private_quantum_seed_hit_search (void);

static u32 generate_dna_ball (u8* goal, u32 wordLen, u32 matchLen,
                              unspos goalEnd,
                              scoreset* scoring, score ballScore,
                              qdjudger judger, void* info);
static u32 judge_qd          (void* info, u8* qWord, u8* dWord,
                              u32 wordLen, u32 matchLen, unspos qEnd);

//----------
//
// quantum_seed_hit_search--
//	Search for high scoring segment pairs (HSPs) between a quantum sequence
//	and a DNA sequence.  HSPs are regions that align with high similarity.
//	Substitutions are allowed, but insertions and deletions are not.
//
// The caller must already have built a table of word positions in the DNA
// sequence.
//
//----------
//
// Arguments:
//	seq*		seq1:			The sequence being searched.
//	postable*	pt:				A table of positions of of words in seq1.
//	s8			charToBits[]:	Table to map DNA characters to two-bit
//								.. values, and illegal characters to -1.
//	seq*		seq2:			The quantum sequence being searched for.
//	unspos		start:			First sequence position to consider.  Zero is
//								.. the first possible position.
//	unspos		end:			One past the last sequence position to consider.
//								.. If this is zero, the sequence length is used.
//	s8			charToBits[]:	Table to map sequence characters to two-bit
//								.. values, and illegal characters to -1.
//	seed*		hitSeed:		The seed-word the table is based on.
//	scoreset*	scoring:		The alignment scoring parameters.  It is assumed
//								.. that a row corresponds to a DNA symbol (i.e.
//								.. sequence1 is DNA).
//	score		ballScore:		The minimum score required of a DNA word to be
//								.. considered 'in' a quantum word's ball.
//	hitprocessor processor:		Function to call for each hit to determine if it
//								.. is 'good enough'.
//	void*		processorInfo:	A value to pass thru with each call to processor.
//
// Returns:
//	The number of HSPs found.
//
// $$$ the return value should be changed to the number of bases covered, to be
// $$$ .. consistent with seed_hit_search()
//
//----------
//
// Notes:
//
// (1)	This routine allocates and reuses memory via global pointers.  The
//		caller should make a call to free_quantum_search() to de-allocate this
//		memory, after all searches are complete.
//
//----------

u32 quantum_seed_hit_search
   (seq*			_seq1,
	postable*		_pt,
	seq*			_seq2,
	unspos			_start,
	unspos			_end,
	const s8		_charToBits[],
	seed*			_hitSeed,
	scoreset*		_scoring,
	score			_ballScore,
	hitprocessor	_processor,
	void*			_processorInfo)
	{

	// sanity check

	if (_hitSeed->resolvingMask != 0)
		suicide ("quantum_seed_hit_search doesn't support overweight seeds");

	if (_hitSeed->type != 'S')
		suicide ("quantum_seed_hit_search only supports strict seeds"
		         " (1s and 0s only)");

	if (_hitSeed->withTrans != 0)
		suicide ("quantum_seed_hit_search doesn't support seeds with transitions");

	if (_end == 0)
		_end = _seq2->len;

	if (_end <= _start)
		suicidef ("in quantum_seed_hit_search(), interval is void (" unsposFmt "-" unsposFmt ")",
		          _start, _end);

	if (_end > _seq2->len)
		suicidef ("in quantum_seed_hit_search(), interval end is bad (" unsposFmt ">" unsposFmt ")",
		          _end, _seq2->len);

	// allocate (or re-use) memory

	empty_diag_hash ();

	// pass globals to the rest of the search
	// note: this makes this module non-threadsafe

	seq1          = _seq1;
	pt            = _pt;
	seq2          = _seq2;
	start         = _start;
	end           = _end;
	charToBits    = _charToBits;
	hitSeed       = _hitSeed;
	scoring       = _scoring;
	ballScore     = _ballScore;
	processor	  = _processor;
	processorInfo = _processorInfo;

	if (_scoring->rowsAreDna) bitsToSym = bits_to_nuc;
	                     else bitsToSym = _scoring->bottleneck;

	// perform search separately for match-seeds vs spaced-seeds

	if (hitSeed->weight == 2*hitSeed->length)
		return private_quantum_word_hit_search ();
	else
		return private_quantum_seed_hit_search ();
	}


void free_quantum_search (void) { free_diag_hash (); }


static u32 private_quantum_word_hit_search
   (void)
	{
	u32		wordLen = (unsigned) hitSeed->length;
	u8*		qStart  = seq2->v;
	u8*		qStop   = seq2->v + seq2->len;
	u8*		q;
	unspos	qPos;
	u32		numHsps = 0;

	if (seq2->len < wordLen)
		return 0; // (nothing to search for)

	// scan the sequence, finding the ball of DNA words 'close' to each quantum
	// word, and processing each seed match therein

	for (q=qStart,qPos=wordLen ; q<=qStop-wordLen ; q++,qPos++)
		numHsps += generate_dna_ball (q, wordLen, wordLen, qPos,
		                              scoring, ballScore, judge_qd, NULL);

	return numHsps;
	}


static u32 private_quantum_seed_hit_search
   (void)
	{
	u8		word[maxSeedLen+1];
	u32		matchLen = (unsigned)  hitSeed->length;
	u32		wordLen  = (unsigned) (hitSeed->weight / 2);
	u8*		qStart   = seq2->v;
	u8*		qStop    = seq2->v + seq2->len;
	u8*		q;
	unspos	qPos;
	u32		ix;
	u32*	shuffle = NULL;
	u32		numHsps = 0;

	if (seq2->len < matchLen)
		return 0; // (nothing to search for)

	// get shuffle list for this seed

	shuffle = seed_shuffle_list (hitSeed);
	if (shuffle[0] != wordLen)
		suicide ("in hsp_quantum_seed_search(), internal error");

	for (ix=0 ; ix<wordLen ; ix++)	// (remove length from shuffle)
		shuffle[ix] = shuffle[ix+1];

	// scan the sequence, finding the ball of DNA words 'close' to each
	// shuffled quantum word, and processing each seed match therein

	for (q=qStart,qPos=matchLen ; q<=qStop-matchLen ; q++,qPos++)
		{
		// fetch the shuffled quantum word at this location
		for (ix=0 ; ix<wordLen ; ix++) word[ix] = q[shuffle[ix]];
		// find matches to words close to it
		numHsps += generate_dna_ball (word, wordLen, matchLen, qPos,
		                              scoring, ballScore, judge_qd, NULL);
		}

	free_if_valid ("hsp_quantum_seed_search (shuffle)", shuffle);

	return numHsps;
	}

//----------
//
// generate_dna_ball--
//	Compute (and report) the 'ball' of DNA words 'close' to a quantum word.
//
// Formally, we compute
//
//		{ d | s(d,q) >= T }
//
// Where
//	d is a dna word
//	q is the goal quantum word
//	T is a similarity score threshold
//	s(.,.) is the similarity score between dna and quantum words, which is
//	       summed over the similarity scores for individual letters
//
//----------
//
// Arguments:
//	u8*			goal:		The word about which the ball will be generated.
//							.. This is a string of quantum characters, but need
//							.. not be zero-terminated.
//	u32			wordLen:	The word length (number of characters in the word).
//	u32			matchLen:	The length of the match the word represents.  This
//							.. can be longer than wordLen when a spaced seed is
//							.. being used.
//	unspos		goalEnd:	Position the goal represents in seq2.  This is
//							.. the index of the first position *after* the word.
//	scoreset*	scoring:	The alignment scoring parameters.  It is assumed
//							.. that a row corresponds to a DNA symbol (i.e.
//							.. sequence1 is DNA).
//	score		ballScore:	The minimum score required of a DNA word to be
//							.. considered 'in' the ball surrounding the
//							.. goal quantum word.
//	qdjudger	judger:		Function to call for each DNA word in the ball,
//							.. to determine if the DNA word results in any
//							.. high-scoring pairs.  This can be NULL if the
//							.. caller just wants to count the number of DNA
//							.. words in the ball.
//	void*		info:		Additional control/arguments specific to the judger
//							.. being called.
//
// Returns:
//	The number of 'good' DNA words in the ball.  Which words are good is
//	determined by the judger function.  If there is no judger all words are
//	considered good.
//
//----------

#define maxWord      16
#define alphabetSize 4

static u32 generate_dna_ball
   (u8*			goal,
	u32			wordLen,
	u32			matchLen,
	unspos		goalEnd,
	scoreset*	scoring,
	score		ballScore,
	qdjudger	judger,
	void*		info)
	{
	score		minNeeded [maxWord];	// minNeeded[i] is min score needed from
										// .. bases 0..i to have a chance to
										// .. reach ballScore
	s8			citizenVal[maxWord];	// (each location is -1, 0, 1, 2, or 3)
	u8			citizenDna[maxWord+1];	// (each location is A, C, G, T)
	u32			dnaWords, goodWords;
	int			ix, sym;
	score		symScore, maxScore, wordScore, bestScore;

	if (wordLen == 0)
		suicidef ("wordLen is zero in generate_dna_ball");
	else if (wordLen > 16)
		suicidef ("wordLen=%u is too large in generate_dna_ball", wordLen);

	quantum_count_stat (qWordsExamined);

	//////////
	// precompute running minimum requirement
	//////////

	// the following #if clause is a workaround for gcc bug 37861, which
	// .. erroneously reports:  "array subscript is above array bounds";  see
	// .. http://gcc.gnu.org/bugzilla/show_bug.cgi?id=37861
	// update dec/2009: since various versions of gcc complain about this, and
	// .. since wordLen>0 due to the fact that it is unsigned and has failed
	// .. the wordLen==0 test above, I've commented-out all the #if stuff, and
	// .. hope that the optimizer is smart enough to remove the wordLen>0 check

//#if ((defined(__GNUC__)) && (GCC_VERSION == 40302))
	if (wordLen > 0) minNeeded[wordLen-1] = ballScore;
//#else
//	minNeeded[wordLen-1] = ballScore;
//#endif

	maxScore = 0;

	for (ix=(int)wordLen-1 ; ix>=0 ; ix--)
		{
		bestScore = scoring->sub[bitsToSym[0]][goal[ix]];
		for (sym=1 ; sym<alphabetSize ; sym++)
			{
			symScore = scoring->sub[bitsToSym[sym]][goal[ix]];
			if (symScore > bestScore) bestScore = symScore;
			}

		if (ix > 0) minNeeded[ix-1] = minNeeded[ix] - bestScore;
		maxScore += bestScore;
		}

	if (quantum_dbgQuantumBall)
		{
		if (maxScore >= ballScore)
			fprintf (stderr, "candidate: %s  maxScore=" scoreFmtSimple "\n",
			                 quantum_word_string (goal, wordLen, 3), maxScore);
#ifdef debugQuantumPartials
		else
			fprintf (stderr, "  no ball: %s"
			                 "  maxScore=" scoreFmtSimple
			                 "<" scoreFmtSimple "\n",
			                 quantum_word_string (goal, wordLen, 3),
			                 maxScore, ballScore);
#endif // debugQuantumPartials
		}

	if (maxScore < ballScore)
		{
		quantum_count_stat (dWordsInBallDistrib[0]);
		return 0;
		}

	//////////
	// generate the ball
	//////////

	dnaWords = goodWords = 0;

	citizenDna[wordLen] = 0;
	citizenVal[0] = -1;
	wordScore = 0;

	ix = 0;
	while (ix >= 0)
		{
		// subtract the score for the symbol in this position

		if (citizenVal[ix] >= 0)
			wordScore -= scoring->sub[citizenDna[ix]][goal[ix]];

		// if we've tried all symbols in this position, backtrack

		if (citizenVal[ix] == alphabetSize-1)
			{ ix--;  continue; }

		// try the next symbol in this position

		citizenVal[ix]++;
		citizenDna[ix] = bitsToSym[(u8)citizenVal[ix]];

		// add score for this symbol, and if it's not enough, prune (and go try
		// the next symbol)

		wordScore += scoring->sub[citizenDna[ix]][goal[ix]];
		if (wordScore < minNeeded[ix])
			{
#ifdef debugQuantumPartials
			if (quantum_dbgQuantumBall)
				fprintf (stderr, "  partial: %s  score=" scoreFmtSimple "\n",
				                 quantum_word_string (citizenDna, ix+1, 3), wordScore);
#endif // debugQuantumPartials
			continue;
			}

		// if we don't have a full word yet, advance to the next position

		if (ix < (int)wordLen-1)
			{ citizenVal[++ix] = -1;  continue; }

		// we have a word that occupies the 'ball'-- report it (and then go try
		// the next symbol)

		dnaWords++;

		if (judger == NULL)
			goodWords++;
		else
			{
			quantum_count_stat (dWordsInBall);
			goodWords += (*judger) (info, goal, (u8*) citizenVal,
			                        wordLen, matchLen, goalEnd);
			}
		}

#ifdef collect_stats
	if (dnaWords < qstatMaxDnaWords-1)
		quantum_count_stat (dWordsInBallDistrib[dnaWords]);
	else
		quantum_count_stat (dWordsInBallDistrib[qstatMaxDnaWords-1]);
#endif // collect_stats

	return goodWords;
	}

//----------
// [[-- qdjudger function --]]
//
// judge_qd--
//	Determine whether a quantum word and a dna word lead to any high-scoring
//	pairs.
//
//----------
//
// Arguments (as per qdjudger functions):
//	void*	info:		(not used)
//	u8*		qWord:		The quantum word.
//	u8*		dWord:		The DNA word (2 bits per letter).
//	u32		wordLen:	The word length.
//	u32		matchLen:	The length of the match the word represents.
//	unspos	qEnd:		Position the quantum word represents in seq2.  This is
//						.. the index of the first position *after* the word.
//
// Returns:
//	The number of seed hits found.
//
// Notes:
//	We assume the scoring parameters are such that a row corresponds to a DNA
//	symbol (i.e. sequence1 is DNA).
//
//----------

static u32 judge_qd
   (arg_dont_complain(void* info),
	u8*			qWord,
	u8*			dWord,
	u32			wordLen,
	u32			matchLen,
	unspos		qEnd)
	{
	unspos		adjStart = pt->adjStart;
	u32			step     = pt->step;
	u32			dnaPacked;
	unspos		pos, pos1;
	u32			ix;
	u32			numHits = 0;

	// pack the dna word into a table index

	dnaPacked = 0;
	for (ix=0 ; ix<wordLen ; ix++)
		dnaPacked = (dnaPacked << 2) + dWord[ix];

	// show the programmer what's going on

	if (quantum_dbgQuantumBall)
		{
		score dnaScore;
		u8    dna[maxWord+1];	// (each location is a bottleneck
								//  .. character, e.g. A, C, G or T)

		dnaScore = 0;
		for (ix=0 ; ix<wordLen ; ix++)
			{
			dna[ix]  =  bitsToSym[dWord[ix]];
			dnaScore += scoring->sub[dna[ix]][qWord[ix]];
			}
		dna[wordLen] = 0;

		fprintf (stderr, "  in ball: %s  score=" scoreFmtSimple,
		                 quantum_word_string (dna, wordLen, 3), dnaScore);

		if (pt->last[dnaPacked] == 0)
			fprintf (stderr, "  no hits\n");
		}

	// process

	if (pt->last[dnaPacked] == 0)
		return 0;

	for (pos=pt->last[dnaPacked] ; pos!=noPreviousPos ; pos=pt->prev[pos])
		{
		pos1 = adjStart + step*pos;
		quantum_count_stat (wordHits);
		numHits += (*processor) (processorInfo, pos1, qEnd, matchLen);
		}

	if (quantum_dbgQuantumBall)
		fprintf (stderr, "  %u hits\n", numHits);

#ifdef collect_stats
	if (numHits < qstatMaxHits-1)
		quantum_count_stat (wordLookupDistrib[numHits]);
	else
		quantum_count_stat (wordLookupDistrib[qstatMaxHits-1]);
#endif // collect_stats

	return numHits;
	}

//----------
//
// quantum_zero_stats--
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

void quantum_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&quantumStats, 0, sizeof(quantumStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// quantum_show_stats--
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

void quantum_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	int		numEntries;
	int		ix, lastIx;

	if (f == NULL) return;

	fprintf (f, "     max sym score: " scoreFmtSimple "\n", quantumStats.maxSymScore);
	fprintf (f, "     min sym score: " scoreFmtSimple "\n", quantumStats.minSymScore);
	fprintf (f, "    max word score: " scoreFmtSimple "\n", quantumStats.maxWordScore);
	fprintf (f, "        ball score: " scoreFmtSimple "\n", quantumStats.ballScore);
	fprintf (f, "            x drop: " scoreFmtSimple "\n", quantumStats.xDrop);
	fprintf (f, "     hsp threshold: %s\n",                 score_thresh_to_string (&quantumStats.hspThreshold));
	fprintf (f, "-------------------\n");
	fprintf (f, "   qWords examined: %s\n",                 commatize (quantumStats.qWordsExamined));
	fprintf (f, "    dWords in ball: %s\n",                 commatize (quantumStats.dWordsInBall));
	fprintf (f, "  dWords per qWord: %.2f\n",               quantumStats.dWordsInBall / (float) quantumStats.qWordsExamined);
	fprintf (f, "         word hits: %s\n",                 commatize (quantumStats.wordHits));
	fprintf (f, "    hits per qWord: %.2f\n",               quantumStats.wordHits / (float) quantumStats.qWordsExamined);
	fprintf (f, "    hits per dWord: %.2f\n",               quantumStats.wordHits / (float) quantumStats.dWordsInBall);
	fprintf (f, "word hits extended: %s\n",                 commatize (quantumStats.wordHitsExtended));
	fprintf (f, "     extended/hits: %.2f%%\n",             100*quantumStats.wordHitsExtended / (float) quantumStats.wordHits);
	fprintf (f, "        HSPs found: %s\n",                 commatize (quantumStats.hspsFound));
	fprintf (f, "         HSPs/hits: %.2f%%\n",             100*quantumStats.hspsFound / (float) quantumStats.wordHits);
	fprintf (f, "-------------------\n");

	fprintf (f, "dWords in ball distribution (by actual usage)\n");
	lastIx = 0;
	numEntries = qstatMaxDnaWords;
	for (ix=0 ; ix<numEntries-1 ; ix++)
		{
		if (quantumStats.dWordsInBallDistrib[ix] == 0) continue;
		lastIx = ix;
		fprintf (f, "              [%2d]: %s\n", ix, commatize (quantumStats.dWordsInBallDistrib[ix]));
		}
	if (quantumStats.dWordsInBallDistrib[numEntries-1] != 0)
		fprintf (f, "             [>%2d]: %s\n", numEntries-2, commatize (quantumStats.dWordsInBallDistrib[numEntries-1]));

	fprintf (f, "-------------------\n");
	fprintf (f, "dWords per lookup distribution (by actual usage)\n");
	lastIx = 0;
	numEntries = qstatMaxHits;
	for (ix=0 ; ix<numEntries-1 ; ix++)
		{
		if (quantumStats.wordLookupDistrib[ix] == 0) continue;
		lastIx = ix;
		fprintf (f, "              [%2d]: %s\n", ix, commatize (quantumStats.wordLookupDistrib[ix]));
		}
	if (quantumStats.wordLookupDistrib[numEntries-1] != 0)
		fprintf (f, "             [>%2d]: %s\n", numEntries-2, commatize (quantumStats.wordLookupDistrib[numEntries-1]));

#endif // collect_stats
	}

void quantum_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
	;  // add these later if desired
	}

