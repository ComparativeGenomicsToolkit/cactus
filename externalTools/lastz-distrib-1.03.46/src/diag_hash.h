//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: diag_hash.h
//
//----------

#ifndef diag_hash_H				// (prevent multiple inclusion)
#define diag_hash_H

// other files

#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff

// establish ownership of global variables

#ifdef diag_hash_owner
#define global
#else
#define global extern
#endif

// debugging defines

//#define noSeedHitQueue		// if this is defined, the seed hit queue is
								// .. NOT used for finding twin seed hits (they
								// .. are found using other techniques)

//----------
//
// diagonal hash--
//	(see diag_hash.c for details)
//
// Helpful formulas:
//	diag = pos1 - pos2			diagNumber(pos1,pos2)
//	pos1 = diag + pos2			diagToPos1(diag,pos2)
//	pos2 = pos1 - diag			diagToPos2(diag,pos1)
//
//----------

// note that diagHashSize must be a power of two

#ifdef puny_diag_hash			// test case; very small hash forces collisions
#define diagHashSize			((u32) 16)
#else							// normal case
#ifdef tiny_diag_hash			// test case; small hash forces collisions
#define diagHashSize			((u32) 1024)
#else							// normal case
#ifdef huge_diag_hash			// test case; no collisions on 1M vs 1M
#define diagHashSize			((u32) 4194304)
#else							// normal case
#ifdef diag_hash_size			// override;  beware: diag_hash_size must be a power of 2
#define diagHashSize			((u32) diag_hash_size)
#else							// normal case
#define diagHashSize			((u32) 65536)
#endif
#endif
#endif
#endif

#define diagNumber(pos1,pos2)	(((sgnpos)(pos1))-((sgnpos)(pos2)))
#define hashedDiag(pos1,pos2)	((diagNumber(pos1,pos2)) & (diagHashSize-1))
#define diagToPos1(diag,pos2)	(((unspos)(diag))+((unspos)(pos2)))
#define diagToPos2(diag,pos1)	(((unspos)(pos1))-((unspos)(diag)))

//--- initialized variables for this module ---

#ifdef diag_hash_owner
unspos*  diagEnd    = NULL;		// arrays (indexed by 0..diagHashSize-1) to
unspos*  diagStart  = NULL;		// .. track the extent of discovered hits on a
sgnpos*  diagActual = NULL;		// .. given diagonal (modulo diagHashSize);
								// .. no values in any of these arrays are valid
								// .. if the corresponding value in diagEnd is
								// .. hashInactiveEnd;  diagEnd (in addition to
								// .. indicating validity) records the end of
								// .. of the most recent hit or extension, as a
								// .. position in sequence 2;  diagStart records
								// .. the start of the most recent hit (or
								// .. series of overlapping hits in some cases);
								// .. diagActual records the actual diagonal
								// .. represented (this is used only to
								// .. detect/resolve collisions)

u32*     diagActive = NULL;		// array (indexed by 0..numdiagActive-1) of
int      numDiagActive;			// .. positions in diagEnd and diagStart that
								// ..  have a non-zero value

//--- external access to the variables, for other modules ---

#else
extern unspos* diagEnd;
extern unspos* diagStart;
extern sgnpos* diagActual;
extern u32*    diagActive;
extern int     numDiagActive;
#endif

// macros for accessing the hash

#define activate_hashed_diag(h) diagActive[numDiagActive++] = h;
#define hashInactiveEnd ((unspos) -1)

//----------
//
// seed hit queue--
//	(see diag_hash.c for details)
//
//----------

#ifndef noSeedHitQueue

#define defaultSeedHitQueueSize	256*1024

typedef struct shqhit			// seed hit queue element
	{
	u64			prevHit;		// the number of the most recent seed hit on
								// .. this hashed diagonal;  0 indicates no hit
	int			isBlock;		// true  => this is an end-of-extension
								// false => this is a seed hit
	unspos		pos2;			// the position of the seed hit in sequence 2;
								// .. this is the position following the end of
								// .. the hit (or extension)
	sgnpos		diag;			// the diagonal containing the seed hit
	} shqhit;

#ifdef diag_hash_owner
int		seedHitQueueSize = 0;	// (N) number of entries in seedHitQueue[]
shqhit*	seedHitQueue = NULL; 	// the most recent N hits, as a cyclic queue
u64*	lastSeedHit = NULL;		// the number of the most recent seed hit on
								// .. each hashed diagonal;  lastSeedHit mod N
								// .. is the index into seedHitQueue[] of that
								// .. seed hit;  0 indicates no hit
u32		seedHitQueueColumns = 0;// the desired number of sequence 2 positions
								// .. covered by the queue;  this is only used
								// .. to generate a warning if we drop below it
u64		seedHitNum;				// identifying number of the most recent seed
								// .. hit (the first seed hit is number N+1)
#else
extern int     seedHitQueueSize;
extern shqhit* seedHitQueue;
extern u64*    lastSeedHit;
extern u32     seedHitQueueColumns;
extern u64     seedHitNum;
#endif

#endif // not noSeedHitQueue

//----------
//
// prototypes for routines in diag_hash.c
//
//----------

void empty_diag_hash  (void);
void free_diag_hash   (void);

#ifndef noSeedHitQueue
void _enqueue_seed_hit (unspos pos1, unspos pos2, int isBlock);
#ifndef diag_hash_owner
#define enqueue_seed_hit(pos1,pos2,isBlock) \
               if (seedHitQueueSize > 0) _enqueue_seed_hit(pos1,pos2,isBlock)
#endif
#endif // not noSeedHitQueue

#undef global
#endif // diag_hash_H
