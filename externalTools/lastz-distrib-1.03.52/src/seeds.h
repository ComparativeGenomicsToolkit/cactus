//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: seeds.h
//
//----------

#ifndef seeds_H					// (prevent multiple inclusion)
#define seeds_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff

// debugging defines

//#define hardCodedSeed			// if this is defined, a hard-coded seed packing
								// .. routine replaces apply_seed; see pack_seed
								// .. below for more details

//----------
//
// data structures and types
//
//----------

// seeds--
//	A structure describing the implementation of a seed;  note that "shift"
//	and mask point to arrays included within the structure, so that only the
//	structure itself, and any seeds down the "next" link, need be freed.  To
//	apply the seed to a sequence s (stored as as two bits per location), the
//	following is computed (see apply_seed):
//		sum over i=0..numParts-1 of (s>>shift[i] & mask[i])

// $$$ some of these should be unsigned

typedef struct seed
	{
	struct seed* next;			// next seed in a linked list

	char	type;				// the type of seed pattern
								//   'S' => strict      (only 1's and 0's)
								//   'H' => half-weight (only T's and 0's)
								//   'R' => overweight  (use resolvingMask)
								//   '_' => nothing special

	int		length;				// seed length, measured in locations
	int		weight;				// seed weight, measured in bits (half locs)
	char*	pattern;			// the seed, represented as a string of 1s
								// .. (matches), 0s (don't cares) and Ts
								// .. (transition-allows);  this may be NULL
	int		numParts;			// length of shift[] and mask[] arrays
	int*	shift;				// array of shift counts, indexed by
								// .. 0..numParts-1
	u32*	mask;				// array of masks, indexed by 0..numParts-1

	u32		resolvingMask; 		// mask (in unpacked bit order) to pick out
								// .. the seed bits which are *not* accounted
								// .. for by shift[] and mask[], and which have
								// .. to be resolved by looking at the sequences

	int		revComp;			// true => pack such that k-mers are identical
								//         ..  to their reverse-complements

	int		isHalfweight;		// false => each unpacked bp is 2-bit nucleotide
								// true  => each unpacked bp is 1-bit R/Y
	int		withTrans;			// non-zero => we allow 1 or 2 transitions in
								//     .. any of the 'match' positions;  assumes
								//     .. that the seed pattern is "strict"
	u32*	transFlips;			// array of words (in packed bit order) for each
								// .. bit that can be flipped to match a
								// .. transition;  each word has only one bit
								// .. set;  terminated by an empty word (zero);
								// .. this can be NULL;  space for this is part
								// .. of this block
	} seed;

#define maxSeedLen       31		// maximum locations allowed in a seed
#define maxHwSeedLen     63		// maximum locations allowed in a half-weight
								// .. seed
#define maxSeedBitWeight 31		// maximum bit weight allowed in a seed;  bit
								// .. weight is 2*matches+transitions
#define maxResolvedBits  16		// maximum number of "extra" bits in a resolving
								// .. seed

#define seed_12of19  "1110100110010101111"
#define seed_14of22  "1110101100110010101111"

//----------
//
// hard-coded seed packing routine
//
// The following gives us a way to plug in a specific seed-packing routine.
// This is provided only as a means to compare timing performance of seed
// packing optimized at compile-time vs the general packing routine.  It is up
// to the user to make sure she specifies the same seed on the command line, or
// all bets are off.
//
// The program seed_function can be used to create the packing routine.
//
//----------

#ifndef hardCodedSeed
u32 apply_seed (seed* seed, u64 word);
#endif

#ifdef hardCodedSeed
#define apply_seed(seed,word) pack_seed(word)
#ifdef straightforwardSeed
static inline u32 pack_seed (u64 word) // pack_1110100110010101111
	{
    return ( word        & 0x000000FF)
         | ((word >>  2) & 0x00000300)
         | ((word >>  4) & 0x00000C00)
         | ((word >>  8) & 0x0000F000)
         | ((word >> 12) & 0x00030000)
         | ((word >> 14) & 0x00FC0000);
	}
#else
static inline u32 pack_seed (u64 word) // pack_1110100110010101111
	{
    return ( word        & 0x00F0CCFF)
         | ((word >> 16) & 0x000F3000)
         | ((word >> 28) & 0x00000300);
	}
#endif
#endif

//----------
//
// prototypes for routines in seeds.c
//
//----------

int    parse_seeds_string        (char* s, seed** seed, int maxIndexBits);
int    parse_strict_seeds_string (char* s, seed** seed, int maxIndexBits);
seed*  reconstruct_seed          (char type, int length, int weight,
                                  char* pattern,
                                  u32 resolvingMask, int revComp,
                                  int isHalfweight, int numParts,
                                  int* shift, u32* mask, u32* transFlips);
seed*  copy_seeds                (seed* seed);
void   free_seeds                (seed* seed);
int    is_same_seed              (seed* seed1, seed* seed2);
char*  seed_pattern              (seed* seed);
u32*   seed_shuffle_list         (seed* seed);
void   print_seeds               (FILE* f, seed* seed);
char*  seed_packed_to_string     (seed* seed, u32 word);
char*  seed_packed_to_string2    (seed* seed, u32 word,
                                  const u8* bitToChar, const u8* bitsToChar);
u64    seed_unpack               (seed* seed, u32 word, u64* seedBits);

#endif // seeds_H
