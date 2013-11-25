//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: seeds.c
//
//----------
//
// seeds--
//	Support for generalized spaced seeds, including transitions.
//
// A "seed pattern" (usually shortened here to "seed") describes a criteria
// for local matches (called "seed hits" or just "hits") between two DNA
// sequences.  The seed specifies which positions in a small local window must
// match and which are allowed to differ.  This is the basis of an anchor-
// finding heuristic for discovering orthology between the two sequences.
//
// A seed is described by a string over the symbols {1,0,T}, with 0 prohibited
// from each end of the string.  A 1 indicates a position which must match, a 0
// allows any mismatch, and a T allows transition mismatches.  Seeds with only
// 1 and 0 bits are called "strict seeds".  For the seed 1101100100101111, the
// two sequences below have a hit at the locations shown, because the sequences
// match at every 1 position in the seed.
//
//		sequence 1:	..GGACCTCTTCTCGCTCTATATAAGCGGTGG..
//		match/trans:      | || ||  | || |||| | :   |
//		sequence 2:	    ..CACTACTGTCGCTCTATATGAACGTGATGT..
//		seed:		        1101100100101111
//
// In addition to the basic criteria above, a seed can further specify that one
// or two transitions are allowed among its match bits (the 1's).  Thus if the
// same seed allowed one transition, we would also have the hit shown below.  If
// the seed did not allow a transition, this would *not* be a hit, in spite of
// the fact that the sequences have more matches in this window that in the
// window above;  matches in the spaces (the 0s) are irrelevant.
//
//		sequence 1:	..GGACCTCTTCTCGCTCTATATAAGCGGTGG..
//		match/trans:      || :| |||||||||||:|:||  :
//		sequence 2:	  ..CACTACTGTCGCTCTATATGAACGTGATGT..
//		                        1101100100101111
//		transitions:                         *
//
// Another extension of the hit criteria is to allow T positions.  A T specifies
// that the position must contain either a match or a transion, but not a
// transversion.  So if our seed were 1101T00100101T11 we would have both of
// those hits.  Note that the transitions allowed by a T location are separate
// from the one or two transitions that can be allowed at match positions.
//
//		sequence 1:	..GGACCTCTTCTCGCTCTATATAAGCGGTGG..
//		match/trans:      | || ||  | || |||| | :   |
//		sequence 2:	    ..CACTACTGTCGCTCTATATGAACGTGATGT..
//		seed:		        1101T00100101T11
//
//		sequence 1:	..GGACCTCTTCTCGCTCTATATAAGCGGTGG..
//		match/trans:      || :| |||||||||||:|:||  :
//		sequence 2:	  ..CACTACTGTCGCTCTATATGAACGTGATGT..
//		                        1101T00100101T11
//
// The "length" (L) of a seed is the number of locations in its string.  This
// corresponds to the length of a corresponding hit.  The "weight" (W) of a
// strict seed is the number of 1s;  for seeds with Ts we count each T as half
// of a 1.  1101100100101111 is called a "10 of 16" seed (length 16, weight 10).
// The same terminology doesn't apply to non-strict seeds;  1101T00100101T11
// has weight 9 and length 16 but is not rightfully a 9 of 16 seed.  The "bit
// weight" of a seed is twice its weight.
//
// Yet another extension is the concept of "half-weight" seeds.  These are seeds
// consisting entirely of Ts and 0s.  This facilitates a heuristic process in
// which the seed is used to identify hits with transitions or matches in the
// prescribed locations, and then the hits are further qualified by requiring
// a minimum number of matches over the length of the seed (including in the
// spaces).
//
// Seed hits are usually found by the following process (this is implemented
// in some other module, but is included in this discussion to give motivation
// for the discussion of "overweight seeds" below).  A window of L nucleotides
// slides across sequence 1.  The nucleotides in the window are converted to a
// two-bit-per-base word, the bits relevant to the seed are packed into a
// smaller word, and a list of all the locations at which that word occurs is
// kept.  Then the second sequence is scanned in a similar manner, with each
// packed word used to locate the list of matching positions.  In practice the
// packed word is used as an index into a table of lists, and so we often call
// it the index.
//
// Here's how a packing for 1101100100101111 might work.  Since we have two bits
// per base, we really need to collect bits in pairs.  The general idea is shown
// below, but we won't go into great detail here about how this is accomplished.
// The key concept is that we extract 2W bits from a word of size 2L, and
// combine the relevant bits to form a unique word of size 2W, the index.
//
//		seed string:	1 1 0 1 1 0 0 1 0 0 1 0 1 1 1 1
//		seed bits:		abcd--efgh----ij----kl--mnopqrst
//		packed bits:	            ghijabcdklefmnopqrst
//
// In practice, the size of the index is limited by the amount of memory
// available for the sequence position table.  On 32-bit machines with 1G of
// memory the practical limit is around 26 index bits.  To allow for heavier
// seeds than this, we handle "overweight seeds".  The actual index is limited
// to the practical maximum, and the remaining seed bits are resolved by
// comparison to the sequence.
//
// For example, consider the weight-14 seed 1110101100110010101111, and suppose
// the maximum index size is 23 bits.  We need to reduce the seed's 28 bits to
// 23, so we change the last 5 1s to Ts: 1110101100110010T0TTTT (why we choose
// this particular reduction will become clear in a moment).  Any hit for the
// full seed will be a hit for this partial seed.  But (in random sequences)
// only 1 in 32 hits to the partial seed are true hits for the full seed.  Upon
// detection of a hit (to the partial seed) we resolve whether it is a true hit
// by comparing the missing bits (the "resolving" bits).  These are the most-
// signifcant bits from each of the five induced T positions.  If these all
// match (as they do in the example below), then the hit is a true hit for the
// full seed.  Due to the way we chose the resolving bits, if they don't match,
// the number of bits that don't match gives the number of transitions over
// those seed positions (mismatches can't be transversions since the least-
// signifcant bits would not have matched).  
//
//		sequence 1:	.. C T C T T C T C G C T C T A T A T A A G C G ..
//		match/trans:   | | |   | | | | : | | | : : |   | : | | | |
//		sequence 2:	.. C T C G T C T C A C T C C G T C T G A G C G ..
//		-------------------------------------------------------------
//		full seed:     1 1 1 0 1 0 1 1 0 0 1 1 0 0 1 0 1 0 1 1 1 1
//		partial seed:  1 1 1 0 1 0 1 1 0 0 1 1 0 0 1 0 T 0 T T T T
//		-------------------------------------------------------------
//		seq1 bits:    01110111110111011001110111001100110000100110
//		seq2 bits:    01110110110111010001110101101101111000100110
//		seq1 r-bits:  --------------------------------1---0-1-0-1-
//		seq2 r-bits:  --------------------------------1---0-1-0-1-
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
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff

#define  seeds_owner			// (make this the owner of its globals)
#include "seeds.h"				// interface to this module

// complain if someone has tried to set limits higher than we can support

#if maxSeedLen > 31
#error ***** maxSeedLen is too large (this module only supports maxSeedLen <= 31) *****
#endif

#if maxHwSeedLen > 63
#error ***** maxHwSeedLen is too large (this module only supports maxHwSeedLen <= 63) *****
#endif

#if maxSeedBitWeight > 31
#error ***** maxSeedBitWeight is too large (this module only supports maxSeedBitWeight <= 31) *****
#endif

#if maxResolvedBits > 16
#error ***** maxResolvedBits is too large (this module only supports maxResolvedBits <= 16) *****
#endif

#define maintainFlippedBitOrder	// if defined, the transition flip order is
								// .. the same for overweight seeds as for the
								// .. equivalent full seed

//----------
//
// prototypes for private functions
//
//----------

static seed* parse_one_seed (char* s, char* e, int transitionsOk,
                             int maxIndexBits);
static seed* new_seed       (int numParts, int patternLen, int numFlips);
static int   best_shift     (u32 uncoveredBits, u64 seedBits);

//----------
//
// parse_seeds_string, parse_strict_seeds_string--
//	Convert a seed(s) string into a collection of bits 'implementing' the
//	seed(s).  parse_strict_seeds_string() doesn't allow transitions ('T').
//
//----------
//
// Arguments:
//	char*	s:				The seeds, represented as a comma-separated string.
//							.. See parse_one_seed() for the format of each seed.
//	seed**	seed:			A place to store the resulting implementation of the
//							.. seeds.  This will be a linked list of seed
//							.. structures that have been allocated from the
//							.. heap, and the caller must eventually dispose of
//							.. it, with a call to free_seeds().
//	int		maxIndexBits:	The maximum weight of the seed (in bits) that can
//							.. be directly supported.  If the seed is heavier
//							.. than this, it will be implemented as an
//							.. overweight seed.  Zero indicates that overweight
//							.. seeds should not be created.
//
// Returns:
//	The *bit* weight of the seed (the maximum of the weights if there are
//	multiple seeds).
//
//----------
//
// Notes:
//
// (1)	Failures result in program fatality.
//
// (2)	No seed can be longer than 31 locations.
//
// (3)	Seed weights are limited by maxIndexBits.  The weight of the seed's
//		index will not exceed maxIndexBits.  Seeds heavier than that will be
//		made as overweight seeds.  maxIndexBits is limited to 31.
//
// (4)	The internal representation of the seed is dependent on the specific
//		2-bit encoding of nucleotides, which is defined (implicitly) in
//		dna_utilities.c.  It requires that the least significant of the two bits
//		distinguishes between purines and pyramidines.
//
//----------

static int _parse_seeds_string (char* s, seed** seed, int transitionsOk,
                                int maxIndexBits);

int parse_seeds_string (char* s, seed** seed, int maxIndexBits)
	{ return _parse_seeds_string (s, seed, /*transitionsOk*/ true, maxIndexBits); }

int parse_strict_seeds_string (char* s, seed** seed, int maxIndexBits)
	{ return _parse_seeds_string (s, seed, /*transitionsOk*/ false, maxIndexBits); }

static int _parse_seeds_string
   (char*	s,
	seed**	_seed,
	int		transitionsOk,
	int		maxIndexBits)
	{
	seed*	tail = NULL;
	seed*	newSeed;
	int		maxWeight;
	char*	terminator;

	//////////
	// convert each comma-separated piece of the string into a seed,
	// collecting them into a linked list
	//////////

	*_seed = NULL;
	maxWeight = -1;

	while (true)
		{
		terminator = strchr (s, ',');
		if (terminator == NULL) terminator = s + strlen (s);

		newSeed = parse_one_seed (s, terminator-1, transitionsOk, maxIndexBits);
		if (*_seed == NULL)
			{
			tail = *_seed = newSeed;
			maxWeight = newSeed->weight;
			}
		else
			{
			tail->next = newSeed;
			tail       = newSeed;
			if (newSeed->weight > maxWeight)
				maxWeight = newSeed->weight;
			}

		if (*terminator == 0)
			break;
		s = terminator+1;
		}

	// return weight, counted in bits

	return maxWeight;
	}

//----------
//
// parse_one_seed--
//	Convert a single seed string into a collection of bits 'implementing' the
//	seed.
//
//----------
//
// Arguments:
//	char*	s:				The seed, represented as a string.  This is a string
//							.. of 1s (representing matches), 0s (don't cares),
//							.. and Ts (transition allows).  X may be substituted
//							.. for 0, and spaces may be used (they are ignored).
//							.. Any leading or trailing don't cares will be
//							.. removed.  The string need not be terminated.
//	char*	e:				The end of the seed (a pointer to the last character
//							.. in the string). 
//	int		transitionsOk:	true  => seed is allowed to contain 'T'.
//							false => seed may only contain match and don't-care.
//	int		maxIndexBits:	The maximum weight of the seed (in bits) that can
//							.. be directly supported.  If the seed is heavier
//							.. than this, it will be implemented as an
//							.. overweight seed.  Zero indicates that overweight
//							.. seeds should not be created.
//
// Returns:
//	The resulting implementation of the seed.
//
//----------
//
// Notes:
//
// (1)	The same restrictions from parse_seed_string() are in force.
//
// (2)	The seed 'implementation' (as a collection of shits and masks) is not
//		optimal in general.
//
//----------

static seed* parse_one_seed
   (char*	_s,
	char*	e,
	int		transitionsOk,
	int		maxIndexBits)
	{
	char	pattern[maxHwSeedLen+1];
	char*	s, *ss, *p;
	char	type;
	int		isStrict, isHalfweight;
	int		length;				// seed length, measured in locations
	int		weight;				// seed weight, measured in bits
	u64		seedBits;			// the seed, as two bits per location
	u64		flipBits;			// the transition flip bits
	u32		resolveBits;		// bits removed from seed, which must be
								// .. resolved by looking at the sequences
	int		shift, bitsPer;		// number of bits to shift the seed
	int		matches;			// number of matches in the seed, and how many
	int		matchesToKeep;		// .. to keep for an overweight seed
	u32		mask;				// mask to apply to the shifted seed
	u32		wBits;				// the W least significant bits, where W=weight
	u32		covered;			// bits (of wBits) which we have covered so far
	u64		remBits;			// bits (of seedBits) which we haven't taken
								// .. care of yet
	int		numParts;			// number of masked-shifts needed
	seed*	seed;

	if (maxIndexBits > maxSeedBitWeight)
		suicidef ("max index bits cannot exceed %d (it's %d).",
				  maxSeedBitWeight, maxIndexBits);

	//////////
	// determine the length and weight of the seed
	//////////

	// skip leading don't cares

	s = _s;

	while ((s <= e) && ((*s == '0') || (*s == 'X') || (*s == 'x')))
		s++;

	if (s > e)
		suicide ("seed string is empty!");

	// skip trailing don't cares

	while ((*e == '0') || (*e == 'X') || (*e == 'x'))
		e--;

	// scan string, to determine if seed is "strict", "half-weight", or a mixture

	isStrict     = true;
	isHalfweight = true;
	matches      = 0;
	weight       = 0;

	for (ss=s ; ss<=e ; ss++)
		{
		switch (*ss)
			{
			case '1':
				isHalfweight = false;
				matches++;
				weight += 2;
				break;

			case 'T':
			case 't':
				isStrict = false;
				weight++;
				break;

			case '0':
			case 'X':
			case 'x':
				break;
			}
		}

	if      (isStrict)     type = 'S';
	else if (isHalfweight) type = 'H';
	else                   type = '_';

	//////////
	// if the seed will be too heavy, turn it into an overweight seed 
	//////////

	matchesToKeep = matches;
	if ((maxIndexBits > 0) && (weight > maxIndexBits))
		{
		int toResolve = weight - maxIndexBits;
		if (toResolve > matches)
			suicidef ("seed (%s) requires more resolving bits (%d) than it has matches (%d).",
			          _s, toResolve, matches);
		if (toResolve > maxResolvedBits)
			suicidef ("seed (%s) requires more resolving bits (%d) than are allowed (%d).",
			          _s, toResolve, maxResolvedBits);
		type = 'R';
		matchesToKeep -= toResolve;
		}

	//////////
	// scan the string, converting each location into one bit (for half-weight
	// seeds) or two bits (for seeds with matches)
	//////////

	resolveBits = 0;
	seedBits    = 0;
	flipBits    = 0;
	bitsPer     = (type=='H')? 1 : 2;
	matches     = 0;
	length      = 0;
	weight      = 0;

	for (ss=s,p=pattern ; ss<=e ; ss++)
		{
		switch (*ss)
			{
			default:
			bad_character:
				if (isprint (*ss))
					suicidef ("seed string %s contains illegal character %c",
					          _s, *ss);
				else
					suicidef ("seed string %s contains illegal character %02X",
					          _s, *ss);

			case ' ':
			case '\t':
			case '\n':
				break;

			case '1':
				if (matches >= matchesToKeep)
					{
					if ((resolveBits << 2) < resolveBits) // (overflow)
						suicidef ("resolving bits in seed string %s are spread too widely",
						          _s);
					resolveBits = (resolveBits << bitsPer) + 2;
					goto transition;
					}
				resolveBits <<= bitsPer;
				seedBits = (seedBits << bitsPer) + 3;
				flipBits = (flipBits << bitsPer) + 2;
				matches++;
				length++;
				weight += 2;
				*(p++) = '1';
				break;

			case 'T':
			case 't':
				if (!transitionsOk) goto bad_character;
				resolveBits <<= bitsPer;
			transition:
				seedBits = (seedBits << bitsPer) + 1;
				flipBits <<= bitsPer;
				length++;
				weight++;
				*(p++) = 'T';
				break;

			case '0':
			case 'X':
			case 'x':
				resolveBits <<= bitsPer;
				seedBits = (seedBits << bitsPer) + 0;
				flipBits <<= bitsPer;
				length++;
				*(p++) = '0';
				break;
			}
		}

	*p = 0; // terminae pattern string

	// sanity check on sizes

	if (type == 'H')
		{
		if (length > maxHwSeedLen)
			suicidef ("half-weight seed string (%s) cannot have length exceeding %d (it's %d).",
			          _s, maxHwSeedLen, length);
		}
	else
		{
		if (length > maxSeedLen)
			suicidef ("seed string (%s) cannot have length exceeding %d (it's %d).",
			          _s, maxSeedLen, length);
		}

	if (weight > maxSeedBitWeight)
		suicidef ("seed string (%s) cannot have bit weight exceeding %d (it's %d).",
		          _s, maxSeedBitWeight, weight);

	if (weight == 0)
		suicidef ("seed string (%s) cannot have zero weight.", _s);

	//////////
	// figure out how to implement the seed
	//
	// we want to find the minimum set of masked-shifts that will bring all
	// the seed bits from scattered positions 0..length-1 into a covering of
	// the positions 0..weight-1;  rather than try to find an optimal set, we
	// use a greedy algorithm to find a good set
	//
	// note:  for many seeds it is possible to find smaller sets of masked-
	//        shifts using more sophistocated (and time-consuming) algorithms
	//////////

	wBits = (1L << weight) - 1;

	// first masked-shift in the set will be shift-zero

	covered  = seedBits & wBits;
	remBits  = seedBits - covered;
	numParts = 1;

	// take whatever masked-shift will cover the most bits;  keep doing so
	// until all bits are covered

	while (covered != wBits)
		{
		shift   = best_shift ((~covered) & wBits, remBits);
		mask    = (remBits >> shift) & (~covered) & wBits;
		covered = covered + mask;
		remBits = remBits - (((u64)mask) << shift);
		numParts++;
		}

	//////////
	// record the seed implementation
	//
	// we run the same algorithm again, depositing the masked-shifts into a
	// seed structure
	//////////

	// allocate the seed structure

	if (type == 'H')
		seed = new_seed (numParts, length, 0);
	else
		seed = new_seed (numParts, length, bit_count_64(flipBits));

	seed->next          = NULL;
	seed->type          = type;
	seed->length        = length;
	seed->weight        = weight;
	seed->isHalfweight  = (type == 'H');
	seed->withTrans     = 0;
	seed->resolvingMask = resolveBits;
	seed->revComp       = false;

	strcpy (seed->pattern, pattern);

	// first masked-shift in the set is shift-zero

	covered        = seedBits & wBits;
	remBits        = seedBits - covered;
	numParts       = 1;
	seed->shift[0] = 0;
	seed->mask [0] = covered;

	// take whatever masked-shift will cover the most bits;  keep doing so
	// until all bits are covered

	while (covered != wBits)
		{
		shift                 = best_shift ((~covered) & wBits, remBits);
		mask                  = (remBits >> shift) & (~covered) & wBits;
		covered               = covered + mask;
		remBits               = remBits - (((u64)mask) << shift);
		seed->shift[numParts] = shift;
		seed->mask [numParts] = mask;
		numParts++;
		}

	// separate the transition-flip bits into a list of single-bit values

	if (seed->transFlips != NULL)
		{
#ifdef maintainFlippedBitOrder
		u64  rightBit;
		u32* f=seed->transFlips;

		while (flipBits != 0)
			{
			rightBit = flipBits-(flipBits&(flipBits-1));// isolate rightmost 1
			flipBits -= rightBit;						// remove it
			*(f++) = apply_seed (seed, rightBit);		// add it to the list
			}
		*f = 0; // terminate the transFlips array
#else // not maintainFlippedBitOrder
		u32  packed = apply_seed (seed, flipBits);
		u32  rightBit;
		u32* f=seed->transFlips;

		while (packed != 0)
			{
			rightBit = packed - (packed & (packed-1));	// isolate rightmost 1
			packed -= rightBit;							// remove it
			*(f++) =  rightBit;							// add it to the list
			}
		*f = 0; // terminate the transFlips array
#endif // maintainFlippedBitOrder
		}

	// return the seed

	return seed;
	}

//----------
//
// new_seed--
//	Allocate a new seed structure.
//
//----------
//
// Arguments:
//	int		numParts:	The number of masked-shifts that will be needed to
//						.. implement the seed.
//	int		patternLen:	The number of bytes to allow for a pattern (we allocate
//						.. one additional byte to allow for a terminating
//						.. zero).
//	int		numFlips:	The number of transition flips to allow for (usually
//						.. this should be the number of full match positions in
//						.. seed).
//
// Returns:
//	A pointer to the newly allocated seed, which the caller will have to
//	dispose of eventually.  The routine free_seeds() should be used for this
//	purpose.
//
//----------

static seed* new_seed
   (int		numParts,
	int		patternLen,
	int		numFlips)
	{
	seed*	s;
	int		bytesNeeded, bytesMain, bytesShift, bytesMask, bytesFlips;

	// figger out how many bytes we need

	bytesMain   = round_up_8 (sizeof(seed));
	bytesShift  = round_up_8 (numParts * sizeof(s->shift[0]));
	bytesMask   = round_up_8 (numParts * sizeof(s->mask [0]));
	bytesFlips  = 0;
	if (numFlips > 0) bytesFlips = (numFlips+1) * sizeof(s->transFlips[0]);
	if      (patternLen > 0) patternLen += 1;
	else if (patternLen < 0) patternLen = 0;
	bytesNeeded = bytesMain + bytesShift + bytesMask + bytesFlips + patternLen;

	// allocate

	s = (seed*) zalloc_or_die ("new_seed", bytesNeeded);

	// hook up the internal arrays

	s->shift = (int*) (((char*) s)        + bytesMain);
	s->mask  = (u32*) (((char*) s->shift) + bytesShift);
	s->transFlips = NULL;
	if (numFlips > 0)
		s->transFlips = (u32*) (((char*) s->mask) + bytesMask);
	if (patternLen > 0)
		s->pattern = (char*) (((char*) s->mask) + bytesMask + bytesFlips);

	// initialize

	s->numParts = numParts;

	return s;
	}

//----------
//
// reconstruct_seed--
//	Build a single seed string from information about a previously-constructed
//	seed.  (Usually this information would come from a file).
//
//----------
//
// Arguments:
//	(the arguemnts have the same meaning as in the seed structure definition)
//
// Returns:
//	The resulting implementation of the seed.
//
//----------

seed* reconstruct_seed
   (char	type,
	int		length,
	int		weight,
	char*	pattern,
	u32		resolvingMask,
	int		revComp,
	int		isHalfweight,
	int		numParts,
	int*	shift,
	u32*	mask,
	u32*	transFlips)
	{
	seed*	s;
	int		numFlips, ix;

	// count flips

	for (numFlips=0 ; transFlips[numFlips]!=0 ; numFlips++)
		;

	// allocate (this links up arrays and sets numParts)

	if (pattern == NULL)
		s = new_seed (numParts, 0, numFlips);
	else
		s = new_seed (numParts, strlen(pattern), numFlips);

	// copy fields

	s->type          = type;
	s->length        = length;
	s->weight        = weight;
	s->resolvingMask = resolvingMask;
	s->revComp       = revComp;
	s->isHalfweight  = isHalfweight;

	if (pattern != NULL)
		strcpy (s->pattern, pattern);

	// copy arrays

	for (ix=0 ; ix<numParts  ; ix++) s->shift[ix]      = shift[ix];
	for (ix=0 ; ix<numParts  ; ix++) s->mask[ix]       = mask[ix];
	for (ix=0 ; ix<=numFlips ; ix++) s->transFlips[ix] = transFlips[ix];

	return s;
	}

//----------
//
// copy_seeds--
//	Make a copy of a list of seed structures.
//
//----------
//
// Arguments:
//	seed*	seed:	The linked list of seeds to copy.  Note that all seeds in
//					.. the list are copied.
//
// Returns:
//	A pointer to the newly allocated seed, which the caller will have to
//	dispose of eventually.  The routine free_seeds() should be used for this
//	purpose.
//
//----------

static seed* copy_seed (seed* _seed);


seed* copy_seeds
   (seed*	_seed)
	{
	seed*	head = NULL;
	seed*	prev = NULL;
	seed*	s;

	for ( ; _seed!=NULL ; _seed=_seed->next)
		{
		s = copy_seed (_seed);
		s->next = NULL;
		if (prev == NULL) head       = s;
		             else prev->next = s;
		prev = s;
		}

	return head;
	}


static seed* copy_seed
   (seed*	_seed)
	{
	u32*	f;
	int		numFlips = 0;
	seed*	s;
	int		ix;

	// allocate a seed with enough room

	numFlips = 0;
	if (_seed->transFlips != NULL)
		{ for (f=_seed->transFlips ; *f!=0 ; f++) numFlips++; }

	if (_seed->pattern == NULL)
		s = new_seed (_seed->numParts, 0, numFlips);
	else
		s = new_seed (_seed->numParts, strlen(_seed->pattern), numFlips);

	// copy the simple fields

	s->next         = NULL;
	s->type         = _seed->type;
	s->length       = _seed->length;
	s->weight       = _seed->weight;
	s->revComp      = _seed->revComp;
	s->isHalfweight = _seed->isHalfweight;

	if (_seed->pattern != NULL)
		strcpy (s->pattern, _seed->pattern);

	// copy parts

	s->numParts      = _seed->numParts;
	s->resolvingMask = _seed->resolvingMask;

	for (ix=0 ; ix<_seed->numParts ; ix++)
		{
		s->shift[ix] = _seed->shift[ix];
		s->mask [ix] = _seed->mask [ix];
		}

	// copy transition flips

	s->withTrans = _seed->withTrans;
	if (_seed->transFlips != NULL)
		{
		for (ix=0 ; ix<numFlips ; ix++)
			s->transFlips[ix] = _seed->transFlips[ix];
		s->transFlips[numFlips] = 0;
		}

	return s;
	}

//----------
//
// free_seeds--
//	De-allocate a list of seed structures.
//
//----------
//
// Arguments:
//	seed*	seed:	The linked list of seeds to de-allocate.
//
// Returns:
//	(nothing)
//
//----------

void free_seeds
   (seed*	_seed)
	{
	seed*	next;

	for ( ; _seed!=NULL ; _seed=next)
		{ next = _seed->next;  free_if_valid ("free_seeds", _seed); }
	}

//----------
//
// is_same_seed--
//	Determine whether two seeds are identical, including having the same index
//	encoding.
//
//----------
//
// Arguments:
//	seed*	seed1, seed2:	The seeds to compare.  Only the single seed pointed
//							.. to is compared;  any remaining seeds in the
//							.. linked list are ignored.
//
// Returns:
//	(nothing)
//
//----------

int is_same_seed
   (seed*	seed1,
	seed*	seed2)
	{
	int		part;

	if (seed1 == seed2) return true;
	if (seed1 == NULL)  return false;
	if (seed2 == NULL)  return false;

	if (seed1->type         != seed2->type)         return false;
	if (seed1->length       != seed2->length)       return false;
	if (seed1->weight       != seed2->weight)       return false;
	if (seed1->revComp      != seed2->revComp)      return false;
	if (seed1->isHalfweight != seed2->isHalfweight) return false;
	if (seed1->withTrans    != seed2->withTrans)    return false;

	if (seed1->numParts     != seed2->numParts)     return false;
	for (part=0 ; part<seed1->numParts ; part++)
		{
		if (seed1->mask [part] != seed2->mask [part]) return false;
		if (seed1->shift[part] != seed2->shift[part]) return false;
		}

	if (seed1->type == 'R')
		{ if (seed1->resolvingMask != seed2->resolvingMask) return false; }

	return true;
	}

//----------
//
// seed_pattern--
//	Create a string describing a list of seed structures.
//
//----------
//
// Arguments:
//	seed*	seed:	The linked list of seeds.
//
// Returns:
//  A string containing the representation of the seed.  This string is
//  actually static data belonging to this routine, so the caller must copy
//  it if more than one such string is to be used simultaneously.
//
//----------

char* seed_pattern
   (seed*	_seed)
	{
	static	char s[70];
	char*	ss;
	seed*	seed;
	int		firstInList;
	int		part;
	u64		seedBits;
	int		bitsPer;
	u32		mask;
	int		loc;
	char	ch;

	// convert each seed in the list to a string of 1TX's, and collect them
	// in the pattern string

	ss = s;
	firstInList = true;

	for (seed=_seed ; seed!=NULL ; seed=seed->next)
		{
		// recover this seed's bits

		seedBits = 0;

		for (part=0 ; part<seed->numParts ; part++)
			seedBits |= ((u64) seed->mask[part]) << seed->shift[part];

		// convert it to a pattern string

		if (!firstInList)
			{
			if (ss-s >= (int) sizeof(s)-1) goto full;
			*(ss++) = ',';
			}

		bitsPer = (seed->type=='H')? 1 : 2;
		mask    = (seed->type=='H')? 1 : 3;

		for (loc=seed->length-1 ; loc>=0 ; loc--)
			{
			switch ((seedBits >> (bitsPer*loc)) & mask)
				{
				default: // (to placate compiler, can't happen)
				case 3: ch = '1'; break;
				case 2: ch = '?'; break;
				case 1: ch = 'T'; break;
				case 0: ch = '0'; break;
				}

			if (ss-s >= (int) sizeof(s)-1) goto full;
			*(ss++) = ch;
			}

		firstInList = false;
		}

	// add resolving bits

	seed = _seed;
	if (seed->type == 'R')
		{
		for (loc=0 ; loc<16 ; loc++)
			{ if (seed->resolvingMask >> (2*loc) == 0) break; }

		if (loc > 0)
			{
			if (ss-s >= (int) sizeof(s)-1) goto full;
			*(ss++) = '/';

			for (loc-- ; loc>=0 ; loc--)
				{
				switch ((seed->resolvingMask >> (2*loc)) & 3)
					{
					default: // (to placate compiler, can't happen)
					case 3: ch = '?'; break;
					case 2: ch = 'R'; break;
					case 1: ch = '?'; break;
					case 0: ch = '0'; break;
					}

				if (ss-s >= (int) sizeof(s)-1) goto full;
				*(ss++) = ch;
				}
			}
		}

	// terminate and return

	*ss = 0;
	return s;

full:
	*ss = 0;
	ss[-3] = '.';
	ss[-2] = '.';
	ss[-1] = '.';

	return s;
	}

//----------
//
// seed_shuffle_list--
//	Create a list of indexes describing how a seed shuffles locations when it
//	is applied.
//
// Example:
//	The seed 11101100101111 is implemented as a list of masks and shifts like
//	this:
//
//		       seed |  1 1 1 0 1 1 0 0 1 0 1 1 1 1 
//		       bits | 1111110011110000110011111111 --> 11111111111111111111
//		------------+------------------------------------------------------
//		mask  shift |
//		F0CFF    0  |         11110000110011111111 --> 11110000110011111111
//		0F000   10  |   1111000000000000           -->     1111000000000000
//		00300   18  | 1100000000                   -->           1100000000
//
//	From the standpoint of rearrangements of a character string, this looks
//	like this:
//
//		       seed | 11101100101111 
//		 characters | ABCDEFGHIJKLMN -->  E  F  B  C  I  A  K  L  M  N
//		------------+-------------------------------------------------
//		F0CFF    0  |     EF  I KLMN -->  E  F        I     K  L  M  N
//		0F000   10  |  BC            -->        B  C
//		00300   18  | A              -->                 A
//		------------+-------------------------------------------------
//		                                  4  5  1  2  8  0 10 11 12 13
//
//	If we consider a pointer pointing at the leftmost character of the original
//	string (the A in this example), we can create the packed word by taking
//	characters at indexes 4 (E), 5 (F), 1 (B), 2 (C), 8 (I), 0 (A), 10 (K),
//	11 (L), 12 (M), 13 (N).  This is the list of indexes returned by this
//	function.
//
//----------
//
// Arguments:
//	seed*	seed:	The seed (if it happens to be a list, only the first seed
//					.. is processed).  This must not contain transitions (only
//					.. match and don't care are allowed).
//
// Returns:
//  A list of indexes, in <length,list> form (see note below).  This is
//	allocated from the heap, and the caller is responsible for disposing of it.
//	Failure results in fatality.
//
// Note:
//	<length,list> form means the first entry in the list is the number of items
//	in the list, NOT including this entry.  So, for example, the list of primes
//	less than ten would be stored as 4,2,3,5,7.  4 is the length and the array
//	has total length 5.
//
//----------

u32* seed_shuffle_list
   (seed*	seed)
	{
	int		length = seed->length;
	int		weight = seed->weight/2;
	u32*	indexes;
	u32		mask, siteBits;
	int		part, origPos, listPos;

	// allocate the index list

	indexes = (u32*) malloc_or_die ("seed_shuffle_list", (1+weight)*sizeof(u32));

	indexes[0] = (unsigned) weight;

	// deposit indexes

	for (part=0 ; part<seed->numParts ; part++)
		{
		mask    = seed->mask[part];
		origPos = (length-1) - (seed->shift[part]/2);
		listPos = weight;

		for ( ; mask!=0 ; mask>>=2,origPos--,listPos--)
			{
			siteBits = mask & 3;
			if (siteBits == 0) continue;
			if (siteBits != 3)
				suicide ("seed contains things other than don't-care and match");

			if (listPos < 1)
				suicide ("internal error, seed weight and masks conflict");

			indexes[listPos] = (unsigned) origPos;
			}
		}

	return indexes;
	}

//----------
//
// print_seeds--
//	Print a list of seed structures.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seed*	seed:	The linked list of seeds to print.
//
// Returns:
//	(nothing)
//
//----------

void print_seeds
   (FILE*	f,
	seed*	seed)
	{
	int		part;
	u64		seedBits;

	for ( ; seed!=NULL ; seed=seed->next)
		{
		// recover the seed

		seedBits = 0;

		for (part=0 ; part<seed->numParts ; part++)
			seedBits |= ((u64) seed->mask[part]) << seed->shift[part];

		fprintf (f, "%016llX\n", (unsigned long long) seedBits);

		// print the masked-shifts

		for (part=0 ; part<seed->numParts ; part++)
			fprintf (f, "  ( >> %2d) & %08X\n",
			            seed->shift[part], seed->mask[part]);

		// print the resolving mask

		if (seed->resolvingMask != 0)
			fprintf (f, "  resolve:   %08X\n", seed->resolvingMask);
		}
	}

//----------
//
// seed_packed_to_string, seed_packed_to_string2--
//	Convert bits, packed as per a seed, to a character string.
//
//----------
//
// Arguments:
//	seed*	seed:		The seed.  This describes how the bits were packed.
//	u32		word:		The nucleotides, packed as per the seed.
//	u8*		bitToChar:	Mapping from a 0/1 value to the character for that
//						.. value (e.g. "RY").
//	u8*		bitsToChar:	Mapping from a 0/1/2/3 value to the character for that
//						.. value (e.g. "ACGT").
//
// Returns:
//  A string containing the nucleotide characters.  This string is actually
//  static data belonging to this routine, so the caller must copy it if more
//  than one such string is to be used simultaneously.
//
//----------

char* seed_packed_to_string (seed* seed, u32 word)
	{ return seed_packed_to_string2 (seed, word, bit_to_pur_pyr, bits_to_nuc); }

char* seed_packed_to_string2
   (seed*		seed,
	u32			word,
	const u8*	bitToChar,
	const u8*	bitsToChar)
	{
	static char	s[maxHwSeedLen+1];
	u64			unpackedWord, unpackedSeed;
	int			numChars;
	int			bitsPer;
	char*		ss;
	u32			twoWordBits,  twoSeedBits, mask;

	// unpack the bits

	unpackedWord = seed_unpack (seed, word, &unpackedSeed);

	// convert to characters

	numChars = seed->length;
	if (numChars > (int) sizeof(s)-1) numChars = sizeof(s)-1;

	// convert each bit pair to a character

	ss = s;
	bitsPer = (seed->type=='H')? 1 : 2;
	mask    = (seed->type=='H')? 1 : 3;

	while (numChars-- > 0)
		{
		twoWordBits = (unpackedWord >> (bitsPer*numChars)) & mask;
		twoSeedBits = (unpackedSeed >> (bitsPer*numChars)) & mask;

		switch (twoSeedBits)
			{
			case 0:     *(ss++) = 'x';                      break;
			case 1: if (twoWordBits < 2)
						*(ss++) = bitToChar[twoWordBits];
					else
						*(ss++) = '?';                      break;
			case 2:     *(ss++) = '?';                      break;
			case 3:     *(ss++) = bitsToChar[twoWordBits];  break;
			}
		}

	*ss = 0;
	return s;
	}

//----------
//
// seed_unpack--
//	Convert bits, packed as per a seed, to a character string.
//
//----------
//
// Arguments:
//	seed*	seed:		The seed.  This describes how the bits were packed.
//	u32		word:		The nucleotides, packed as per the seed.
//	u64*	seedBits:	Place to return the bits that are 'active' in the seed.
//						.. This may be NULL.
//
// Returns:
//  The nucleotides, unpacked to their original bit format.  Any bits that are
//	not 'active' in the seed are zero.
//
//----------

u64 seed_unpack
   (seed*	seed,
	u32		word,
	u64*	seedBits)
	{
	u64		unpackedWord, unpackedSeed, partMask;
	int		part;

	//////////
	// unpack the bits
	//////////

	unpackedWord = 0; // the bits we have
	unpackedSeed = 0; // the bits we could have had

	for (part=0 ; part<seed->numParts ; part++)
		{
		partMask     =  (u64) seed->mask[part];
		unpackedWord |= (word & partMask) << seed->shift[part];
		unpackedSeed |=         partMask  << seed->shift[part];
		}

	if (seedBits != NULL) *seedBits = unpackedSeed;
	return unpackedWord;
	}

//----------
//
// apply_seed--
//	Apply a seed to a word, extracting and packing the seed bits.
//
//----------
//
// Arguments:
//	seed*	seed:	The seed.  This describes how to extract bits from the
//					.. word and pack them.
//	u64		word:	A word of consecutive nts.  This must have at least as
//					.. many nts as the seed needs (seed->length).  If it has
//					.. more nts than that, the extras are ignored.
//
// Returns:
//	The packed word.
//
//----------

#ifndef hardCodedSeed

u32 apply_seed
   (seed*	seed,
	u64		word)
	{
	int		part;
	u64		rcWord = 0;
	u32		packedWord = 0;
	int		seedBits = 0;

	// perform reverse-complement if necessary;  for complementing seeds both
	// a k-mer and its reverse complement are represented by whichever is
	// numerically lowest

	if (seed->revComp)
		{
		if (seed->type == 'H')			// half-weight seed
			{
			seedBits = seed->length;
			rcWord   = rev_comp_by_bits (word, seed->length);
			}
		else if (seed->type == 'R')		// overweight seed
			suicide ("internal error: overweight seeds cannot be complementing");
		else
			{
			seedBits = 2*seed->length;
			rcWord   = rev_comp_by_pairs (word, seed->length);
			}

		word &= (((u64)1)<<seedBits) - 1;
		//printf ("%0*llX/%s",
		//        (seedBits+3)/4, word,   seed_packed_to_string(seed,word));
		//printf (" -> %0*llX/%s\n",
		//        (seedBits+3)/4, rcWord, seed_packed_to_string(seed,rcWord));

		if (rcWord < word) word = rcWord;
		}

	// apply the seed by combining the masked-shifts of the word

	for (part=0 ; part<seed->numParts ; part++)
		packedWord |= (word >> seed->shift[part]) & seed->mask[part];

	return packedWord;
	}

#endif // not hardCodedSeed

//----------
//
// best_shift--
//	Determine the amount to shift a seed to fill in the most uncovered bits.
//
//----------
//
// Arguments:
//	u32		uncoveredBits:	The bits we want to cover.
//	u64		seedBits:		The bits in the seed that we can shift into any
//							.. coverage position.
//
// Returns:
//	The best right shift count.
//
//----------

static int best_shift
   (u32		uncoveredBits,
	u64		seedBits)
	{
	int		coverage, bestCoverage;
	int		shift,    bestShift;

	bestCoverage = -1;
	bestShift    = -1;

	for (shift=0 ; seedBits!=0 ; seedBits>>=1,shift++)
		{
		coverage = bit_count (seedBits & uncoveredBits);
		if (coverage > bestCoverage)
			{ bestCoverage = coverage;  bestShift = shift; }
		}

	return bestShift;
	}

// Scraps no longer need, saved for future use.
//
// unused routine that determines a seed's corresponding bit pattern
//
//static u64 seed_bits
//   (seed*	seed)
//	{
//	u64		seedBits;
//	int		part;
//
//	if (seed->type == 'R') seedBits = seed->resolvingMask;
//	                  else seedBits = 0;
//
//	for (part=0 ; part<seed->numParts ; part++)
//		seedBits |= ((u64) seed->mask[part]) << seed->shift[part];
//
//	return seedBits;
//	}

