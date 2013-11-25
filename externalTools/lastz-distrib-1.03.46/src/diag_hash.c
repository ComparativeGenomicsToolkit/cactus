//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: diag_hash.c
//
//----------
//
// diag_hash--
//	Support for managing a hash table of values associated with alignment
//	diagonals.
//
// Our search algorithms need to track hits on each diagonal to detect
// overlapping hits and nearby multiple hits.  For each diagonal, we store the
// position in the second sequence of the end of most recent hit (pos2), and in
// some cases we also store the start position of a run of overlapping or
// nearly-overlapping hits.
//
// Ideally we would have an entry for each possible diagonal.  But with
// sequence lengths approaching a quarter billion, such an array would require
// too much memory (memory usage is nominally 16*(L1+L2) bytes, so two .25Gbase
// sequences would require 8Gbytes).  Instead, we use a hash array much smaller
// than the number of hits likely to be "active" at any one time.  It is
// critical that we store pos2 for this purpose rather than pos1, because pos2
// always increases during the search while pos1 hops around the sequence.  By
// guaranteeing the order that we (effectively) scan each diagonal, we limit
// our exposure to collisions in the hash.  We also store the number of the
// actual diagonal associated with each hash location, so that we can detect
// collisions.
//
// The effect of a collision is that a hit that should be combined with a
// nearby one to its left won't be.  Because the second sequence is being
// scanned from left to right, a hit stays "dangerous" (i.e., it can cover up
// a subsequent hit on a different hash-equivalent diagonal) only until the
// sweep along sequence 2 reaches pos2 for the end of that hit.  If we were
// storing pos1, a hit could stay dangerous for the entire sweep.
//
// The hash funtion used simply uses some number of the least significant bits
// of the diagonal.  For this reason the hash size H (diagHashSize) must be a
// power of 2.  Hash-equivalent diagonals are not random but are separated by H
// intervening diagonals.
//
// Helpful formulae:
//	diag = pos1 - pos2
//	pos1 = diag + pos2
//	pos2 = pos1 - diag
//
// We could reduce memory use by using two (or more) hashes with different
// moduli (relatively prime).  The tradeoff is that we'd be spending more time
// computing the hash values and reading/writing them from/to the hash arrays.
//
// Allocation size... By default we allocate 65K entries for these tables.  In
// a 31-bit sequence index enviroment, this works out to 1M bytes (4 tables
// with 32-bit entries => 4*4*65K = 1M).  This can be overridden at compile
// time by defining diag_hash_size (to some power of 2).  For example, setting
// diag_hash_size to 4M in a 32-bit sequence index enviroment works out to 84M
// (3 tables with-32 bit entries and one with 64-bit entries => 5*4*4M = 84M).
//
//----------
//
// Seed hit queue
// ==============
//
// The seed hit queue keeps track of the N most recent seed hits, searchable by
// hashed diagonal.  It is assumed that N is large enough that even in the worst
// case we have all the hits for the most recent W words in seqeunce 2, where W
// is the span allowed for twin seed hits.
//
// The queue is implemented as a simple cyclic array with a write head that
// wraps around.  This saves us the effort of disposing of stale hit records;
// they are just overwritten.  Within the queue, we keep a linked list of the
// hits belonging to each hashed diagonal.
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
#include <limits.h>				// standard C value limit stuff
#include "build_options.h"		// build options

#define  diag_hash_owner		// (make this the owner of its globals)
#include "diag_hash.h"			// interface to this module

// debugging defines

//#define debugDiag 23420-14467	// if defined, breakdown what happens on this
								// .. pos1-pos2 diagonal

//----------
//
// prototypes for private functions
//
//----------

//static void dump_seed_hit_queue (void);

//----------
//
// empty_diag_hash--
//	Create or re-create an empty diagonals hash and seed hit queue.
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
//
// Notes:
//
// (1)	This routine allocates and reuses memory via global pointers.  The
//		caller should make a call to free_diag_hash() to de-allocate this
//		memory when it is no longer needed.
//
//----------

void empty_diag_hash (void)
	{
	u32 hDiag;

#ifndef noSeedHitQueue
	//printf ("\n(erasing seed hit queue)\n");
#endif // not noSeedHitQueue

	// allocate (or re-use) memory;  note that only diagEnd requires that we
	// fill it with zeros

	if (diagEnd == NULL)
		{
		// allocate new array
		// caveat: sometimes we don't need all of this (e.g. process_for_simple_hit)
		//         .. so allocating all these together is wasteful

		diagEnd    = malloc_or_die ("empty_diag_hash (diagEnd)",
		                            diagHashSize * sizeof(unspos));
		diagStart  = malloc_or_die ("empty_diag_hash (diagStart)",
		                            diagHashSize * sizeof(unspos));
		diagActual = malloc_or_die ("empty_diag_hash (diagActual)",
		                            diagHashSize * sizeof(sgnpos));
		diagActive = malloc_or_die ("empty_diag_hash (diagActive)",
		                            diagHashSize * sizeof(u32));

		for (hDiag=0 ; hDiag<diagHashSize ; hDiag++)
			diagEnd[hDiag] = hashInactiveEnd;

		numDiagActive = 0;

#ifdef debugDiag
		printf ("init: (diag %9s|%9s|%04X) end is %d\n",
		        "","", debugDiag, hashInactiveEnd);
#endif

#ifndef noSeedHitQueue
		if (seedHitQueueSize > 0)
			{
			seedHitQueue = zalloc_or_die ("empty_diag_hash (seedHitQueue)",
			                              seedHitQueueSize * sizeof(shqhit));
			lastSeedHit  = zalloc_or_die ("empty_diag_hash (lastSeedHit)",
			                              diagHashSize * sizeof(u64));
			seedHitNum   = seedHitQueueSize;
			}
		else
			{
			seedHitQueue = NULL;
			lastSeedHit  = NULL;
			}
#endif // not noSeedHitQueue
		}
	else
		{
		// clear the previously-allocated data

#ifndef noSeedHitQueue
		if (seedHitQueueSize > 0)
			{
			while (numDiagActive > 0)
				{
				hDiag = diagActive[--numDiagActive];
				diagEnd[hDiag] = hashInactiveEnd;
				lastSeedHit[hDiag] = 0;
#ifdef debugDiag
				if (hDiag == hashedDiag (debugDiag,0))
					printf ("insq: (diag %9s|%9s|%04X) end is %d\n",
					        "","", hDiag, hashInactiveEnd);
#endif
				}
			seedHitNum = seedHitQueueSize;
			}
		else
			{
#endif // not noSeedHitQueue
			while (numDiagActive > 0)
				{
				hDiag = diagActive[--numDiagActive];
				diagEnd[hDiag] = hashInactiveEnd;
#ifdef debugDiag
				if (hDiag == hashedDiag (debugDiag,0))
					printf ("isq:  (diag %9s|%9s|%04X) end is %d\n",
					        "","", hDiag, hashInactiveEnd);
#endif
				}
#ifndef noSeedHitQueue
			}
#endif // not noSeedHitQueue
		}

	}

//----------
//
// free_diag_hash--
//	Dispose of the diagonals hash and seed hit queue.
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

void free_diag_hash
   (void)
	{
	free_if_valid ("free_diag_hash (diagEnd)",      diagEnd);       diagEnd      = NULL;
	free_if_valid ("free_diag_hash (diagStart)",    diagStart);     diagStart    = NULL;
	free_if_valid ("free_diag_hash (diagActual)",   diagActual);    diagActual   = NULL;
	free_if_valid ("free_diag_hash (diagActive)",   diagActive);    diagActive   = NULL;

#ifndef noSeedHitQueue
	free_if_valid ("free_diag_hash (seedHitQueue)", seedHitQueue);  seedHitQueue = NULL;
	free_if_valid ("free_diag_hash (lastSeedHit)",  lastSeedHit);   lastSeedHit  = NULL;
#endif // not noSeedHitQueue
	}

//----------
//
// enqueue_seed_hit--
//	Record a seed hit in the seed hit queue.
//
// Nota Bene: we name this _enqueue_seed_hit here so that that a macro in
//            .. diaghash.h can make this a conditional call
//
//----------
//
// Arguments:
//	unspos		pos1:		The hit position in sequence 1.
//	unspos		pos2:		The hit position in sequence 2.
//	int			isBlock:	true  => this is an end-of-extension
//							false => this is a seed hit
//
// Returns:
//	(nothing)
//
//----------
//
// Assumptions:
//	seedHitQueueSize > 0
//	lastSeedHit  != NULL
//	seedHitQueue != NULL
//
//----------

#ifndef noSeedHitQueue

void _enqueue_seed_hit
   (unspos		pos1,
	unspos		pos2,
	int			isBlock)
	{
	static int haveShortfall = false;
	sgnpos		diag;
	u32			hDiag;
	shqhit*		q;

	diag  = diagNumber (pos1, pos2);
	hDiag = hashedDiag (pos1, pos2);

	//printf ("\n[%04X] adding " sgnposFmt "/" unsposFmt "\n",
	//        hDiag, diag, pos2);

	seedHitNum++;
	q = &seedHitQueue[seedHitNum % seedHitQueueSize];

	if (seedHitNum > (u32) (2*seedHitQueueSize))
		{
		//printf ("[%04X] removing " sgnposFmt "/" unsposFmt " (" unsposFmt " columns)\n",
		//        hashedDiag (q->diag,0), q->diag, q->pos2, pos2-q->pos2);
		if ((!haveShortfall)
		 && (!q->isBlock)
		 && (pos2 - q->pos2 <= seedHitQueueColumns))
			{
			// $$$ with some work we could realloc the queue here;  we'd need
			//     .. to rotate the entries so that the current pointer is at
			//     .. the old end (so new entries go into the virgin area and
			//     .. don't overwrite old entries);  this would require updating
			//     .. every prevHit link
			haveShortfall = true;
			fprintf (stderr,
			         "seed hit queue shortfall at " unsposSlashFmt "\n",
			         (unspos) (diag+pos2), pos2);
			}
		}

	if (lastSeedHit[hDiag] <= seedHitNum - seedHitQueueSize)
		q->prevHit = 0; // (last seed hit is stale, no longer in queue)
	else
		q->prevHit = lastSeedHit[hDiag];
	lastSeedHit[hDiag] = seedHitNum;

	q->isBlock = isBlock;
	q->pos2    = pos2;
	q->diag    = diag;

	//dump_seed_hit_queue ();
	}

#endif // not noSeedHitQueue

//----------
//
// dump_seed_hit_queue--
//	Dump the contents of the seed hit queue.
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

#ifndef noSeedHitQueue

//static void dump_seed_hit_queue
//   (void)
//	{
//	u32		hDiag;
//	u64		ix;
//	shqhit*	q;
//
//	for (hDiag=0 ; hDiag<diagHashSize ; hDiag++)
//		{
//		ix = lastSeedHit[hDiag];
//		if (ix == 0) continue;
//		if (ix <= seedHitNum - seedHitQueueSize)
//			continue;	// (the index is stale)
//
//		printf ("[%04X]", hDiag);
//		while (ix > seedHitNum - seedHitQueueSize)
//			{
//			q = &seedHitQueue[ix % seedHitQueueSize];
//			printf (" " sgnposFmt "/" unsposFmt, q->diag, q->pos2);
//			ix = q->prevHit;
//			}
//		printf ("\n");
//		}
//	}

#endif // not noSeedHitQueue

