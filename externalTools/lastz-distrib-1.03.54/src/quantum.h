//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: quantum.h
//
//----------

#ifndef quantum_H				// (prevent multiple inclusion)
#define quantum_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff
#include "seed_search.h"		// seed hit search stuff

// establish ownership of global variables

#ifdef quantum_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef quantum_owner
int quantum_dbgQuantumBall = false;	// true => show details of the ball of DNA
									// .. words close to each quantum word
#else
global int quantum_dbgQuantumBall;
#endif


//----------
//
// data structures and types
//
//----------

// quantum-to-dna judger functions--
//	Judge whether a quantum-word and dna-word comprise a high-scoring pair.
//
//	Arguments:
//		void*	info:		(pass-thru argument)
//		u8*		qword:		The quantum word.  This is a string of characters,
//							.. but is *not* zero-terminated.
//		u8*		dword:		The DNA word, with nucleotides encoded as two bits.
//		u32		wordLen:	The word length (number of characters in the words).
//		u32		matchLen:	The length of the match the word represents;  this
//							.. can be longer than wordLen when a spaced seed is
//							.. being used.
//		unspos	qEnd:		Position the quantum word represents in some
//							sequence.  This is the index of the first position
//							*after* the word.
//
//	Returns:
//		The number of HSP's derived from this pair.

typedef u32 (*qdjudger) (void*, u8*, u8*, u32, u32, unspos);

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	score	maxSymScore;
	score	minSymScore;
	score	maxWordScore;
	score	ballScore;
	score	xDrop;
	sthresh	hspThreshold;
	u32		qWordsExamined;
    u32		dWordsInBall;
    u32		wordHits;
    u32		wordHitsExtended;
    u32		hspsFound;
    u32		dWordsInBallDistrib[22];
    u32		wordLookupDistrib[152];
	} quantumStats;

#define qstatMaxDnaWords entriesof(quantumStats.dWordsInBallDistrib)
#define qstatMaxHits     entriesof(quantumStats.wordLookupDistrib)

// stats macros

#define quantum_count_stat(field)   ++quantumStats.field
#define quantum_uncount_stat(field) --quantumStats.field
#define quantum_set_stat(field,val) (quantumStats.field = val)
#define quantum_add_stat(field,val) (quantumStats.field += val)
#else
#define quantum_count_stat(field)
#define quantum_uncount_stat(field)
#define quantum_set_stat(field,val)
#define quantum_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void quantum_zero_stats    (void);
void quantum_show_stats    (FILE* f);
void quantum_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in quantum.c
//
//----------

u32  quantum_seed_hit_search (seq* seq1, postable* pt,
                              seq* seq2, unspos start, unspos end,
                              const s8 charToBits[], seed* hitSeed,
                              scoreset* scoring, score ballScore,
                              hitprocessor processor, void* processorInfo);
void free_quantum_search     (void);

#undef global
#endif // quantum_H
