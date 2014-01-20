//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: pos_table.h
//
//----------

#ifndef pos_table_H				// (prevent multiple inclusion)
#define pos_table_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "seeds.h"				// seed matching stuff

// establish ownership of global variables

#ifdef pos_table_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access
// nota bene: showProgress is a relic which is not currently used

#ifdef pos_table_owner
int   pos_table_showProgress    = false; // true => make periodic progress reports
int   pos_table_dbgShowWords    = false; // true => show 'words' in add_word()
int   pos_table_dbgShowDiscards = false; // true => show 'words' discarded in
seed* pos_table_dbgSeed         = NULL;  //         .. limit_position_table()
#else
global int   pos_table_showProgress;
global int   pos_table_dbgShowWords;
global int   pos_table_dbgShowDiscards;
global seed* pos_table_dbgSeed;
#endif

//----------
//
// data structures and types
//
//----------

// position dumper functions--
//	Dump one position (from a position table) to a file.
//
//	Arguments:
//		FILE*		f:			The file to print to.
//		postable*	pt:			The table containing the position.
//		int			field:		Which field to dump (one of posdump_xxx).
//		u64			fieldVal:	The value of the field being dumped.

typedef void (*posdumper) (FILE*, void*, int, u64);

enum
	{
	posdump_index = 0,		// dump index
	posdump_index_space,	// dump as much space as index
	posdump_word,			// dump word
	posdump_word_space,		// dump as much space as word
	posdump_position		// dump position
	};

// position table--
//	A position table maps a word to a list of positions of that word in the
//	associated sequence.  Think of a word as a W-mer with two bits per
//	nucleotide (though it may be bits selected from a larger window and then
//	packed).
//
//	Positions are recorded as the index (into a subinterval of the sequence) of
//	the first character *after* the end of the word.  For example, "abra" would
//	be recorded as positions 4 and 11 in "abracadabra".  Positions are relative
//	to a (start,end) subinterval of the sequence.  If the sequence is
//	"elabrabracadabrazil" and the subinterval is (5,16(, the positions of "abra"
//	are still 4 and 11.
//
//		index:    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
//		sequence: e  l  a  b  r  a  b  r  a  c  a  d  a  b  r  a  z  i  l
//		interval:               (-------------------------------(
//		position:                0  1  2  3  4  5  6  7  8  9 10
//		"abra":                              X                    X
//
//	This implementation consists of two arrays-- last and prev.  Last has an
//	entry for every possible word and gives the rightmost position in the
//	sequence where that word exists.  Prev has an entry for each position in
//	the sequence and gives the position of the first duplicate of the word, to
//	its left.  Thus all the positions of a given word can be found by traversing
//	a chain of integers (the equivalent of a linked list), with a special value
//	(noPrevious) indicating the end of the list.  In contrast, an empty list has
//	last set to zero.  Any position for which a seed is *not* stored in the
//	table (i.e. is not part of a linked list) has prev set to zero.
//
//	The reason for the incongruity for the end marker (0 for empty list,
//	noPrevious for end of list) is that it easier to allocate memory filled
//	with zeros, thus the entries in last and prev that we never touch will
//	contain zero.  For dynamic masking, in order to remove seeds from the table,
//	we need to be able to quickly determine whether a position (to be masked) is
//	in the table or not.  This is distinguishable by prev zero (not in table) or
//	non-zero (in table).
//
//	The positions in both last and prev are indexes into prev, and thus are
//	positions (relative to granularity) in the *subinterval* of the sequence
//	defined by start and end.  Granularity (G) specifies a subset of the
//	positions that can be stored.  If G=1, every position can be stored;  if
//	G=2 only even positions can be stored; if G=5 only every 5th position can
//	be stored.  Since the values in last and prev are indexes into prev, the
//	following conversions are helpful:
//		P = S + G*I
//		I = (P-S)/G
//	where
//		P = the position in the sequence
//		S = the adjusted start of the subinvterval (adjStart)
//		I = the index into prev
//		G = the granularity
//
//	note (1):	The memory for last, prev, and (if needed) asBits, is usually
//				allocated as part of the same block as the postable struct.
//				Thus the whole kit and kaboodle can be deallocated with a
//				single call to free (but you should use free_position_table).
//				In some cases, however, the caller will install pointers for
//				last, prev and asBits to point elsewhere.  In this case the
//				caller is responsible for making sure they get cleaned up.

typedef struct postable
	{
	u32				allocLast;	// actually number of entries allocated for the
	unspos			allocPrev;	// .. last and prev arrays;  some may be unused

	int				wordBits;	// number of *bits* in a word (for a simple
								// .. W-mer match, this would be 2W)
	u32				wordEntries;// number of entries in last[]

	posdumper		dump;		// function to dump a position to the console
								// .. (only used for debugging; can be NULL)
	void*			dumpInfo;	// custom argument for dump()

	unspos			start;		// first sequence position of interest
	unspos			end;		// sequence position just beyond the last
								// .. position of interest
	unspos			adjStart;	// sequence position corresponding to first
								// .. entry in prev[];  this is guaranteed to
								// .. be a multiple of postable.step
	u32				step;		// granularity of positions (see note above)

	// (see note (1) about whether last[], prev[], and asBits[] arrays are
	// allocated *within* this same malloc block, or not

	unspos*			last;		// array giving the rightmost position of each
								// .. word (see text above) (also see note 1)
	unspos*			prev;		// array giving the position of the first
								// .. duplicate to the left of each position
								// .. (see text above) (also see note 1);  the
								// .. first entry in this array corresponds to
								// .. sequence postion postable.adjStart
	u32*			asBits;		// a packed version of the sequence the table
								// .. catalogs;  these are two bits per bp
								// .. (encoded using whatever upperCharToBits[]
								// .. array is in use);  the first word
								// .. corresponds to sequence locations adjStart
								// .. thru adjStart+15, with adjStart in the two
								// .. most significant bits;  this field can be
								// .. NULL; (see note 1)
	} postable;

#define noPreviousPos ((unspos) -1)

// position count distribution--

typedef struct poscount
	{
	unspos			count;		// the number of positions some particular seed
								// .. word occurs at in the sequence
	unspos			occurrences;// the number of words which have that count; a
								// .. zero indicates the end of a list
	} poscount;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int   wordWeight;
	int   wordSpace;
	int   wordsInTable;
    int   wordsPresent;
    int   singletonWords;
    int   basesParsed;
    u32   wordCountLimit;
    u32   maxWordCountChasm;
	int   discardedWords;
	int   protectedWords;
	int   intervalsMasked;
	int   maskedIntervalBases;
	int   wordsRemovedFromTable;
    int   maskBasesParsed;
	} posTableStats;

// stats macros

#define pos_table_count_stat(field)   ++posTableStats.field
#define pos_table_uncount_stat(field) --posTableStats.field
#define pos_table_set_stat(field,val) (posTableStats.field = val)
#define pos_table_add_stat(field,val) (posTableStats.field += val)
#else
#define pos_table_count_stat(field)
#define pos_table_uncount_stat(field)
#define pos_table_set_stat(field,val)
#define pos_table_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void pos_table_zero_stats       (void);
void pos_table_show_stats       (FILE* f, postable* pt);
void pos_table_show_stats_after (FILE* f);

//----------
//
// prototypes for routines in pos_table.c
//
//----------

postable* build_seed_position_table (seq* seq, unspos start, unspos end,
                                     const s8 upperCharToBits[],
                                     seed* seed, u32 step);
postable* build_quantum_seed_position_table
                                    (seq* seq, unspos start, unspos end,
                                     u8* bottleneck, const charvec qToBest[],
                                     seed* seed, u32 step);
void      mask_seed_position_table  (postable* pt,
                                     seq* seq, unspos start, unspos end,
                                     const s8 upperCharToBits[], seed* hitSeed);
postable* new_position_table        (int wordBits, unspos start, unspos end,
                                     u32 step, int allocLast, int allocPrev,
                                     int allocBits);
void      free_position_table       (postable* pt);
u32       fetch_resolving_bits      (postable* pt, unspos pos1);
void      dump_position_table       (FILE* f, postable* pt, seed* hitSeed,
                                     int showPositions, int showCounts);
unspos    count_position_table      (postable* pt);
void      limit_position_table      (postable* pt, u32 limit, u32 maxChasm);
poscount* position_table_count_distribution (postable* pt);
u32       find_position_table_limit (postable* pt, float keep);

#undef global
#endif // pos_table_H
