//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: pos_table.c
//
//----------
//
// pos_table--
//	Support for creating a table of positions of words in genomic sequences.
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

#define  pos_table_owner		// (make this the owner of its globals)
#include "pos_table.h"			// interface to this module

// debugging defines

//#define debugTablePos1 111	// if defined, breakdown what happens with this
								// .. position in sequence 1 (this is the right
								// .. end of the word;  the left end is at
								// .. X-(seedLength-1) counting from 1 at the
								// .. start of the sequence

//----------
//
// private data
//
//----------

typedef struct dumpinfo
	{
	seed*	seed;				// the seed used to pack words
	u8*		bitsToAlphabet;		// mapping from 2-bit values to characters (e.g.
								// .. "ACGT")
	} dumpinfo;

//----------
//
// stats to augment crude profiling
//
//----------

#ifndef dbgTiming
#define dbg_timing_set_stat(field,val)     ;
#define dbg_timing_count_stat(field)       ;
#define dbg_timing_report_stat(field,name) ;
#endif // not dbgTiming

#ifdef dbgTiming
struct
	{
	int   wordsInTable;
	} posTableTimingStats;

#define dbg_timing_set_stat(field,val)     (posTableTimingStats.field = val)
#define dbg_timing_count_stat(field)       ++posTableTimingStats.field
#define dbg_timing_report_stat(field,name) fprintf(stderr,"%-26s %d\n", \
                                              name":",posTableTimingStats.field)
#endif // dbgTiming

//----------
//
// prototypes for private functions
//
//----------

static void record_seed_positions
                         (postable* pt, seq* seq,
                          const s8 upperCharToBits[], seed* hitSeed);
static void record_seed_positions_halfweight
                         (postable* pt, seq* seq,
                          const s8 upperCharToBits[], seed* hitSeed);
static void record_seed_positions_bits
                         (postable* pt, seq* seq,
                          const s8 upperCharToBits[], seed* hitSeed);
static void record_seed_positions_quantum
                         (postable* pt, seq* seq,
                          const charvec qToBest[], seed* hitSeed);
static void mask_seed_positions
                         (postable* pt, seq* seq, unspos start, unspos end,
                          const s8 upperCharToBits[], seed* hitSeed);
static void mask_seed_positions_halfweight
                         (postable* pt, seq* seq, unspos start, unspos end,
                          const s8 upperCharToBits[], seed* hitSeed);
static int  position_is_in_table (postable* pt, unspos position);
static void add_word     (postable* pt, u32 word, unspos position);
static void remove_word  (postable* pt, u32 word, unspos position);
static void dump_word_position
                         (FILE* f, postable* pt, int field, u64 fieldVal);
static void dump_seed_position
                         (FILE* f, postable* pt, int field, u64 fieldVal);
static void dump_quantum_seed_position
                         (FILE* f, postable* pt, int field, u64 fieldVal);

//----------
//
// build_seed_position_table--
//	Create a table of the positions of all seed-words in an interval of a
//	sequence.  The basic idea of a word is a series of W consecutive bases, but
//	it is generalized to that of a seed containing W (specific) bits over L
//	bases.
//
//----------
//
// Arguments:
//	seq*	seq:				The sequence to build the position table of.
//	unspos	start:				First sequence position to consider.  Zero is
//								.. the first possible position.
//	unspos	end:				One past the last sequence position to consider.
//								.. If this is zero, the sequence length is used.
//	s8		upperCharToBits[]:	Table to map sequence characters to two-bit
//								.. values, and illegal characters to -1.
//	seed*	hitSeed:			The seed-word to base the table on.
//	u32		step:				Positional step size indicating the granularity
//								.. of the positions stored.  For example, step=5
//								.. means only every 5th position will be
//								.. stored.  Step=1 means all positions are
//								.. stored.
//
// Returns:
//	A pointer to a newly allocated table of positions;  failures result in
//	program fatality.  The caller must eventually dispose of the table, with a
//	call to free_position_table().
//
//----------

postable* build_seed_position_table
   (seq*		seq,
	unspos		start,
	unspos		end,
	const s8	upperCharToBits[],
	seed*		hitSeed,
	u32			step)
	{
	postable*	pt;
	dumpinfo	di;

	// sanity check

	if (step < 1)
		suicidef ("in build_seed_position_table(), step can't be %u", step);

	if (end == 0)
		end = seq->len;

	if (end <= start)
		suicidef ("in build_seed_position_table(), interval is void (" unsposDotsFmt ")",
		          start, end);

	if (end > seq->len)
		suicidef ("in build_seed_position_table(), interval end is bad (" unsposFmt ">" unsposFmt ")",
		          end, seq->len);

	// create an empty table

	pt = new_position_table (hitSeed->weight, start, end, step,
	                         true, true, (hitSeed->type == 'R'));

	// install dumper

	pt->dump     = (posdumper) dump_seed_position;
	pt->dumpInfo = &di;
	di.seed           = hitSeed;
	di.bitsToAlphabet = NULL;

	// fill the table

	if (hitSeed->isHalfweight)
		record_seed_positions_halfweight (pt, seq, upperCharToBits, hitSeed);
	else if (pt->asBits != NULL)
		record_seed_positions_bits       (pt, seq, upperCharToBits, hitSeed);
	else
		record_seed_positions            (pt, seq, upperCharToBits, hitSeed);

	return pt;
	}

//----------
//
// build_quantum_seed_position_table--
//	Create a table of the positions of all seed-words in an interval of a
//	sequence (similar to build_seed_position_table).  Here a word consists of
//	W quantum bases, which is reduced to the closest (highest scoring) word in
//	the bottleneck alphabet (which can be thought of as A, C, G, T;  see note
//	below).
//
//----------
//
// Arguments:
//	seq*		seq:			The sequence to build the position table of.
//	unspos		start:			First sequence position to consider.  Zero is
//								.. the first possible position.
//	unspos		end:			One past the last sequence position to consider.
//								.. If this is zero, the sequence length is used.
//	u8*			bottleneck:		The bottleneck alphabet.
//	charvec		qToBest[]:		Table to map a quantum character to the two-bit
//								.. code(s) for the 'closest' bottleneck
//								.. character(s).
//	seed*		hitSeed:		The seed-word to base the table on.
//	u32			step:			Positional step size indicating the granularity
//								.. of the positions stored.  For example, step=5
//								.. means only every 5th position will be
//								.. stored.  Step=1 means all positions are
//								.. stored.
//
// Returns:
//	A pointer to a newly allocated table of positions;  failures result in
//	program fatality.  The caller must eventually dispose of the table, with a
//	call to free_position_table().
//
//----------
//
// (1)	The bottleneck alphabet is invisible to this routine.  Its effect is
//		completely described by the qToBits[] table.
//
//----------

postable* build_quantum_seed_position_table
   (seq*			seq,
	unspos			start,
	unspos			end,
	u8*				bottleneck,
	const charvec	qToBest[],
	seed*			hitSeed,
	u32				step)
	{
	postable*		pt;
	dumpinfo		di;

	// sanity check

	if (step < 1)
		suicidef ("in build_quantum_seed_position_table(), step can't be %u", step);

	if (end == 0)
		end = seq->len;

	if (end <= start)
		suicidef ("in build_quantum_seed_position_table(), interval is void (" unsposDotsFmt ")",
		          start, end);

	if (end > seq->len)
		suicidef ("in build_quantum_seed_position_table(), interval end is bad (" unsposFmt ">" unsposFmt ")",
		          end, seq->len);

	if (hitSeed->type != 'S')
		suicide ("(internal error in build_quantum_seed_position_table: strict seeds only)\n");

	// create an empty table

	pt = new_position_table (hitSeed->weight, start, end, step,
	                         true, true, (hitSeed->type == 'R'));

	// install dumper

	pt->dump     = (posdumper) dump_quantum_seed_position;
	pt->dumpInfo = &di;
	di.seed           = hitSeed;
	di.bitsToAlphabet = bottleneck;

	// fill the table

	record_seed_positions_quantum (pt, seq, qToBest, hitSeed);

	return pt;
	}

//----------
//
// record_seed_positions,
// record_seed_positions_halfweight,
// record_seed_positions_bits--
// record_seed_positions_quantum--
//	Record the positions of all spaced-seed words in (a subinterval of) a
//	sequence.  The subinterval is defined by (pt->start,pt->end).  The only
//	difference between these versions is
//	  - the normal version encodes at two bits per nucleotide (before packing)
//	  - the half-weight version encodes only one bit per nucleotide
//	  - the bits version encodes as two bits and also makes a copy of the
//	    the sequence as a bit stream
//	  - the quantum version is just like the normal version, but breaks ties
//	    when mapping quantum characters to bit pairs
//
//----------
//
// Arguments:
//	postable* pt:				The position table in which to record the
//								.. positions.
//	seq*	seq:				The sequence to build the position table of.
//	s8		upperCharToBits[]:	(see note 2) Table to map sequence characters
//								.. to two-bit values,and illegal characters to
//								.. -1.
//	seed*	hitSeed:			The seed-word to base the table on.
//
// Returns:
//	(nothing)
//
//----------
//
// Notes:
//
// (1)	record_seed_positions_halfweight is dependent on the specific 2-bit
//		encoding of nucleotides, which is defined (implicitly) in
//		dna_utilities.c.  We assume that the least significant of the two bits
//		distinguishes between purines and pyramidines.
//
// (2)	record_seed_positions_quantum replaces the upperCharToBits argument with
//		qToBest:
//
//			charvec qToBest[]:	Table to map a quantum character to the two-bit
//								.. code(s) for the 'closest' bottleneck
//								.. character(s).
//
//----------
//
// Normally we slide a 64-bit window along the sequence and collect a
// bit-packed version of the nucleotides in that window.  64 bits corresponds
// to 32 bases, or, for record_seed_positions_halfweight, 64 bases.  When we
// have accumulated enough bits to satsify the seed, we pack them according to
// the seed and record the word/position pair in the table.
//
// Words that contain 'illegal' bases are excluded from the table.  This is
// accomplished by restarting the collection of bits whenever an illegal base
// is encountered.  The corresponding positions are never recorded in the table,
// thus their prev values remain zero.  All positions recorded in the table have
// prev non-zero.
//
// A step size can be used to limit the number of locations stored in the
// table.  This reduces memory needs and also increases overall speed (since
// later processing will have fewer matches to deal with).  Only positions that
// are multiples of the step size are stored.  For 'short' step sizes (step no
// longer than the seed) all bases are collected, but packing and recording are
// only performed at such positions.
//
// For 'long' step sizes (step longer than seed), we only collect bases that can
// possibly be part of the word for such a location.  This is accomplished by
// "skip-ahead", moving the sequence pointer past several useless bases.  There
// are two places we should skip ahead.  The first is when we have recorded a
// word.  The other is when we hit a bad base and restart collection.
//
// In the following examples, step size Z=15 and seed length L=10.  oo indicates
// a good base, xx a bad base, and -- a base we don't care about.
//
// In the first case, after we have just recorded a word, we have something like
// this:
//
//	oo oo -- -- -- -- --[oo oo oo oo oo oo oo oo oo oo]-- -- -- -- -- oo oo
//	28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51
//	       s              *                             *
//
// When s is 30 we report the word from 20..29.  At this point the next word we
// would possibly report would occur when s=45.  So we should skip ahead to
// 35.  The formula is s' = s + L-Z.
//
// In the second case, we have something like this:
//
//	-- --[oo oo xx -- -- -- -- -- -- --]-- -- -- -- -- oo oo oo oo oo oo oo
//	33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56
//	                s                                   *
//
// We were scanning along toward 45 and encounted a bad nucleotide at 37 (and s
// has already been incremented to 38).  This kills any chance of having a word
// to record at 45.  The next possible word will come when s=60, so we need to
// skip ahead to s=50.  The formula is s' = s + Z-1 - (s+L-1 mod Z).
//
// Since the result of mod will be in the range 0..Z-1, s' will be in the range
// of s to s+Z-1.  We want the result to be -L modulo Z so that the next L bases
// form a word when we reach Z.  Below we prove this formula produces s' == -L,
// and since s <= s' < Z, s' is the desired position.
//
//	s' == s + Z-1 - (s+L-1 mod Z)  modulo Z
//	   == s - 1 - (s+L-1)
//	   == s - 1 - s - L + 1
//	   == s-s - L + 1-1
//	   == -L
//
//----------

static void record_seed_positions
   (postable*	pt,
	seq*		seq,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + pt->start;
	u8*			seqStop  = seq->v + pt->end;
	u8*			s;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;

	seedLength = (unsigned) hitSeed->length;

	if (seq->len < seedLength)
		return; // (nothing to search for)

	// scan the sequence, adding each packed word to the table

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		if (step > seedLength)						// for large steps, skip
			{										// .. ahead to the next
			pos = s - seq->v;						// .. viable start position
			s = s + (step-1) - ((pos+seedLength-1) % step);
			}

	empty_skipped:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0)					// .. on step boundary
				{
#ifdef debugTablePos1
				if (pos == debugTablePos1)
					printf ("seq 1 pos " unsposFmt " not on z-step boundary\n", pos);
#endif
				continue;
				}

			packed = apply_seed (hitSeed, w);		// extract seed bits
			add_word (pt, packed, pos);				// add it to the table
#ifdef debugTablePos1
			if (pos == debugTablePos1)
				printf ("recording %s at seq 1 pos " unsposFmt "\n",
				        seed_packed_to_string (hitSeed, packed), pos);
#endif

			if (step > seedLength)					// for large steps, skip
				{									// .. directly to the start
				s += step - seedLength;				// .. of the next possible
				goto empty_skipped;					// .. word we'll record
				}
			}
		}

	}


static void record_seed_positions_halfweight
   (postable*	pt,
	seq*		seq,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + pt->start;
	u8*			seqStop  = seq->v + pt->end;
	u8*			s;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;

	seedLength = (unsigned) hitSeed->length;

	if (seq->len < seedLength)
		return; // (nothing to search for)

	// scan the sequence, adding each packed word to the table

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		if (step > seedLength)						// for large steps, skip
			{										// .. ahead to the next
			pos = s - seq->v;						// .. viable start position
			s = s + (step-1) - ((pos+seedLength-1) % step);
			}

	empty_skipped:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0) continue;		// .. on step boundary

			packed = apply_seed (hitSeed, w);		// extract seed bits
			add_word (pt, packed, pos);				// add it to the table

			if (step > seedLength)					// for large steps, skip
				{									// .. directly to the start
				s += step - seedLength;				// .. of the next possible
				goto empty_skipped;					// .. word we'll record
				}
			}
		}

	}


static void record_seed_positions_bits
   (postable*	pt,
	seq*		seq,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + pt->start;
	u8*			seqStop  = seq->v + pt->end;
	u8*			s;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;
	u32*		tp;
	u32			asBits;
	int			numNts;

	seedLength = (unsigned) hitSeed->length;

	if (seq->len < seedLength)
		return; // (nothing to search for)

	// scan the sequence, adding each packed word to the table, and accum-
	// ulating the 'top bits' for each nucleotide
	//
	// notes: (1) 'bad' characters are encoded the same as some good character
	//            (probably a T);  this does not cause problems because no
	//            position is recorded (in the table) that would use any of
	//            those bad bits
	//        (2) unlike the other two record_seed_positions routines, here we
	//            cannot skip ahead when we have long step sizes, because we
	//            we miss encoding the intervening nucleotides

	tp = pt->asBits;
	asBits = 0;
	numNts = (int) (pt->start - pt->adjStart);

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			asBits = (asBits << 2) | (ww & 3);
			if (++numNts == 16)
				{ *(tp++) = asBits;  numNts = 0; }
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (basesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			asBits = (asBits << 2) | (ww & 3);
			if (++numNts == 16)
				{ *(tp++) = asBits;  numNts = 0; }
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0) continue;		// .. on step boundary

			packed = apply_seed (hitSeed, w);		// extract seed bits
			add_word (pt, packed, pos);				// add it to the table
			}
		}

	if (numNts > 0)
		*tp = asBits << (2*(16-numNts));
	}


static void record_seed_positions_quantum
   (postable*	pt,
	seq*		seq,
	const charvec qToBest[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + pt->start;
	u8*			seqStop  = seq->v + pt->end;
	u8*			s;
	u8			ch;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;
	int			numTied;

	seedLength = (unsigned) hitSeed->length;

	if (seq->len < seedLength)
		return; // (nothing to search for)

	// scan the sequence, adding each packed word to the table

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		if (step > seedLength)						// for large steps, skip
			{										// .. ahead to the next
			pos = s - seq->v;						// .. viable start position
			s = s + (step-1) - ((pos+seedLength-1) % step);
			}

	empty_skipped:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (basesParsed);
			ch = *(s++);							// fetch next char
			numTied = qToBest[ch].len;
			if (numTied == -1) goto empty;			// bad char => start over
			if (numTied == 1)
				ww = qToBest[ch].v[0];				// map next char
			else
				ww = qToBest[ch].v[(s-seq->v)%numTied];	// map next char
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (basesParsed);
			ch = *(s++);							// fetch next char
			numTied = qToBest[ch].len;
			if (numTied == -1) goto empty;			// bad char => start over
			if (numTied == 1)
				ww = qToBest[ch].v[0];				// map next char
			else
				ww = qToBest[ch].v[(s-seq->v)%numTied];	// map next char
			w = (w << 2) | ww;						// append next nt

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0)					// .. on step boundary
				{
#ifdef debugTablePos1
				if (pos == debugTablePos1)
					printf ("seq 1 pos " unsposFmt " not on z-step boundary\n", pos);
#endif
				continue;
				}

			packed = apply_seed (hitSeed, w);		// extract seed bits
			add_word (pt, packed, pos);				// add it to the table
#ifdef debugTablePos1
			if (pos == debugTablePos1)
				printf ("recording %s at seq 1 pos " unsposFmt "\n",
				        seed_packed_to_string (hitSeed, packed), pos);
#endif

			if (step > seedLength)					// for large steps, skip
				{									// .. directly to the start
				s += step - seedLength;				// .. of the next possible
				goto empty_skipped;					// .. word we'll record
				}
			}
		}

	}

//----------
//
// mask_seed_position_table--
//	Remove masked seeds from a position table.  A masked seed is one that
//	contains a masked base.
//
//----------
//
// Arguments:
//	postable*	pt:			The position table to operate on.
//	unspos		start,end:	The range of sequence positions to consider.  Any
//							.. seed enclosed in this range is removed from the
//							.. table.  Origin-0, end-exclusive.  If end==0, the
//							.. sequence length is used.
//	(all other arguments are as per build_seed_position_table)
//
// Returns:
//	(nothing)
//
//----------

void mask_seed_position_table
   (postable*	pt,
	seq*		seq,
	unspos		start,
	unspos		end,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	// sanity check

	if (end == 0)
		end = seq->len;

	if (end <= start)
		suicidef ("in mask_seed_position_table(), interval is void (" unsposFmt "-" unsposFmt ")",
		          start, end);

	if (end > seq->len)
		suicidef ("in mask_seed_position_table(), interval end is bad (" unsposFmt ">" unsposFmt ")",
		          end, seq->len);

	pos_table_count_stat (intervalsMasked);
	pos_table_add_stat   (maskedIntervalBases, end-start);

	// mask the table
	//
	// note that if the table contains a copy of the the sequence as a bit
	// stream (when pt->asBits is non NULL), we are unable to mask the bits in
	// that bit stream (because it has just two bits per base and provides no
	// way to encode a masked base);  this causes no problem, since we remove
	// the corresponding seed the bits we don't clear will never be used

	if (hitSeed->isHalfweight)
		mask_seed_positions_halfweight (pt, seq, start, end, upperCharToBits, hitSeed);
	else
		mask_seed_positions            (pt, seq, start, end, upperCharToBits, hitSeed);
	}

//----------
//
// mask_seed_positions,
// mask_seed_positions_halfweight,
//	Remove the positions of any masked spaced-seed words in (a subinterval of)
//	a sequence.  The subinterval is defined by the (start,end) arguments, not
//	the start,end values in the position table structure.  The only
//	difference between these versions is
//		- the normal version encodes at two bits per nucleotide (before packing)
//		- the half-weight version encodes only one bit per nucleotide
//
//----------
//
// Arguments:
//	postable* pt:				The position table in which to un-record the
//								.. positions.
//	seq*	seq:				The sequence the position table was built for.
//	unspos	start, end:			The interval to check (same meaning as for
//								.. mask_seed_position_table).  Origin-0, end-
//								.. exclusive.
//	s8		upperCharToBits[]:	Character to bit mapping that was used to build
//								.. the table.
//	seed*	hitSeed:			The seed-word to table is based on.
//
// Returns:
//	(nothing)
//
//----------

static void mask_seed_positions
   (postable*	pt,
	seq*		seq,
	unspos		start,
	unspos		end,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + start;
	u8*			seqStop  = seq->v + end;
	u8*			s;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;

	seedLength = (unsigned) hitSeed->length;

	if (end-start < seedLength)
		return; // (nothing to search for)

	// scan the sequence, removing each packed word from the table

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		if (step > seedLength)						// for large steps, skip
			{										// .. ahead to the next
			pos = s - seq->v;						// .. viable start position
			s = s + (step-1) - ((pos+seedLength-1) % step);
			}

	empty_skipped:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (maskBasesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (maskBasesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 2) | ww;						// append next nt

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0) continue;		// .. on step boundary

			if (!position_is_in_table(pt,pos)) 		// make sure position is
				continue;							// .. currently in the table

			packed = apply_seed (hitSeed, w);		// extract seed bits
			remove_word (pt, packed, pos);			// remove it from the table

			if (step > seedLength)					// for large steps, skip
				{									// .. directly to the start
				s += step - seedLength;				// .. of the next possible
				goto empty_skipped;					// .. word we'll record
				}
			}
		}

	}


static void mask_seed_positions_halfweight
   (postable*	pt,
	seq*		seq,
	unspos		start,
	unspos		end,
	const s8	upperCharToBits[],
	seed*		hitSeed)
	{
	u32			step = pt->step;
	u32			seedLength;
	u8*			seqStart = seq->v + start;
	u8*			seqStop  = seq->v + end;
	u8*			s;
	u64			w;
	s32			ww;
	u32			packed;
	u32			nts;
	unspos		pos;

	seedLength = (unsigned) hitSeed->length;

	if (end-start < seedLength)
		return; // (nothing to search for)

	// scan the sequence, removing each packed word to the table

	for (s=seqStart ; s<seqStop ; )
		{
		// collect the first seedLength-1 nucleotides

	empty:
		if (step > seedLength)						// for large steps, skip
			{										// .. ahead to the next
			pos = s - seq->v;						// .. viable start position
			s = s + (step-1) - ((pos+seedLength-1) % step);
			}

	empty_skipped:
		w = 0L;
		for (nts=1 ; (nts<seedLength)&&(s<seqStop) ; nts++)
			{
			pos_table_count_stat (maskBasesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y
			}

		// process each word of seedLength nucleotides

		for ( ; s<seqStop ; )
			{
			pos_table_count_stat (maskBasesParsed);
			ww = upperCharToBits[*(s++)];			// map next char
			if (ww < 0) goto empty;					// bad char => start over
			w = (w << 1) | (ww & 1);				// append next R/Y

			pos = s - seq->v;						// make sure position is
			if ((pos % step) != 0) continue;		// .. on step boundary

			if (!position_is_in_table(pt,pos)) 		// make sure position is
				continue;							// .. currently in the table

			packed = apply_seed (hitSeed, w);		// extract seed bits
			remove_word (pt, packed, pos);			// remove it from the table

			if (step > seedLength)					// for large steps, skip
				{									// .. directly to the start
				s += step - seedLength;				// .. of the next possible
				goto empty_skipped;					// .. word we'll record
				}
			}
		}

	}

//----------
//
// new_position_table--
//	Allocate a new, empty, position table structure.
//
//----------
//
// Arguments:
//	int		wordBits:	The number of *bits* in a word (roughly speaking, twice
//						.. the number of nucleotides).
//	unspos	start, end:	Range of sequence positions that will be used.
//	u32		step:		The granularity of the positions that will be stored.
//	int		allocLast:	true => allocate space for the last[] array.
//	int		allocPrev:	true => allocate space for the prev[] array.
//	int 	allocBits:	true => allocate space in which to save 2 bits per bp.
//
// Returns:
//	A pointer to the newly allocated position table, which the caller will
//	have to dispose of eventually.  The routine free_position_table() should
//	be used for this purpose.
//
//----------
//
// total memory used (give or take a few bytes) is
//		4 * (2^W + L/G)         =   2^(W+2) + 4L/G
// where
//		W = wordBits
//		L = sequence length (end-start)
//		G = granularity (step)
//
// When allocBits is true an additional L/4 bytes is used.  This is usually
// insignificant relative to the rest.  For example, this is equal to 4L/G if G
// is 16 (an absurdly large value for G).
//
// some examples:
//
//		 W |    L | G |  mem  |
//		---+------+---+-------+
//		20 | 250M | 1 |  .94G | chr1, 10-of-L seed
//		20 | 250M | 2 |  .47G |
//		---+------+---+-------+
//		24 | 250M | 1 |  .99G | chr1, 12-of-L seed
//		24 | 250M | 2 |  .53G |
//		---+------+---+-------+
//		24 |  50M | 1 |  .25G | chr21, 12-of-L seed
//		24 |  50M | 2 |  .16G |
//		---+------+---+-------+
//		28 |   2M | 1 | 1.01G | ENCODE region, 14-of-L seed
//		28 |   2M | 2 | 1.00G |
//		---+------+---+-------+
//
//----------
//
// Relationship of start,end,adjStart,step
//
//		start = 33        adjStart = start - (start%step)
//		end   = 47           "     = 33 - 3
//		step  = 5            "     = 30
//
//
//	sequence:  ..|28|29|30|31|32|33|34|35|36|37|38|39|40|41|42|43|44|45|46|47|..
//	interval:                   (-----------------------------------------(
//	prev:              [ 0]           [ 1]           [ 2]           [ 3]
//	                                                ^  ^  ^  ^  ^  ^  ^  ^  ^
//	for word length = 6                             |  |  |  |  |  |  |  |  |
//	window end = 39, discarded  (-----------------( +  |  |  |  |  |  |  |  |
//	window end = 40, saved as 2    (-----------------( +  |  |  |  |  |  |  |
//	window end = 41, discarded        (-----------------( +  |  |  |  |  |  |
//	window end = 42, discarded           (-----------------( +  |  |  |  |  |
//	window end = 43, discarded              (-----------------( +  |  |  |  |
//	window end = 44, discarded                 (-----------------( +  |  |  |
//	window end = 45, saved as 3                   (-----------------( +  |  |
//	window end = 46, discarded                       (-----------------( +  |
//	window end = 47, discarded                          (-----------------( +
//
//----------

postable* new_position_table
   (int			wordBits,
	unspos		start,
	unspos		end,
	u32			step,
	int			allocLast,
	int			allocPrev,
	int			allocBits)
	{
	postable*	pt;
	unspos		adjStart;
	u32			wordEntries;
	unspos		prevEntries;
	u64			bytesNeeded, bytesStruct, bytesLast, bytesPrev, bytesAsBits;

	if (wordBits > 28)
		suicidef ("new_position_table can't support >28 seed bits (%d requested)", wordBits);

	// figger out how many bytes we need

	adjStart     = start - (start % step);	// (force adjStart down to a
											//  multiple of the granularity)

	wordEntries = ((u32) 1) << wordBits;
	prevEntries = 1 + ((end-adjStart) / step);

	bytesStruct = round_up_16 (sizeof(postable));
	bytesNeeded = bytesStruct;

	bytesLast   = 0;
	bytesPrev   = 0;
	bytesAsBits = 0;

	if (allocLast)
		{
		bytesLast   =  round_up_16 (((u64) wordEntries) * sizeof(pt->last[0]));
		bytesNeeded += bytesLast;
		if (bytesLast > mallocLimit) goto overflow_last;
		}

	if (allocPrev)
		{
		bytesPrev   =  round_up_16 (((u64) prevEntries) * sizeof(pt->prev[0]));
		bytesNeeded += bytesPrev;
		if (bytesPrev > mallocLimit) goto overflow_prev;
		}

	if (allocBits)
		{
		bytesAsBits =  round_up_16((end-adjStart+3) / 4);
		bytesNeeded += bytesAsBits;
		if (bytesAsBits > mallocLimit) goto overflow_as_bits;
		}

	if (bytesNeeded > mallocLimit) goto overflow;

	//fprintf (stderr, "wordBits    = %d\n", wordBits);
	//fprintf (stderr, "wordEntries = %s\n", commatize(wordEntries));
	//fprintf (stderr, "bytesLast   = %s\n", commatize(bytesLast));
	//fprintf (stderr, "\n");
	//fprintf (stderr, "start       = %s\n", commatize(start));
	//fprintf (stderr, "end         = %s\n", commatize(end));
	//fprintf (stderr, "prevEntries = %s\n", commatize(prevEntries));
	//fprintf (stderr, "bytesPrev   = %s\n", commatize(bytesPrev));
	//fprintf (stderr, "\n");
	//fprintf (stderr, "bytesAsBits = %s\n", commatize(bytesAsBits));
	//fprintf (stderr, "\n");
	//fprintf (stderr, "bytesNeeded = %s\n", commatize(bytesNeeded));

	// allocate

	pt = (postable*) zalloc_or_die ("new_position_table", bytesNeeded);

	// initialize control fields

	pt->allocLast   = wordEntries;
	pt->allocPrev   = prevEntries;
	pt->wordBits    = wordBits;
	pt->wordEntries = wordEntries;
	pt->start       = start;
	pt->adjStart    = adjStart;
	pt->end         = end;
	pt->step        = step;
	pt->dump        = NULL;
	pt->dumpInfo    = NULL;

	// hook up the internal arrays;  note that we do not need to initialize
	// their contents, since allocation filled them with zeros

	pt->last   = (unspos*) (((char*) pt)       + bytesStruct);
	pt->prev   = (unspos*) (((char*) pt->last) + bytesLast);
	pt->asBits = (u32*)    (((char*) pt->prev) + bytesPrev);

	if (bytesLast   == 0) pt->last   = NULL;
	if (bytesPrev   == 0) pt->prev   = NULL;
	if (bytesAsBits == 0) pt->asBits = NULL;

	return pt;

// failure exits

#define suggestions " consider using lastz_m40,"                            \
                    " or setting max_malloc_index for a special build,"     \
                    " or breaking your target sequence into smaller pieces"


overflow:
	{
	char* tempStruct = commatize(bytesStruct);
	char* tempLast   = commatize(bytesLast);
	char* tempPrev   = commatize(bytesPrev);
	char* tempAsBits = commatize(bytesAsBits);
	suicidef ("in new_position_table(), structure size (%s+%s+%s+%s = %s) exceeds allocation limit of %s;"
	          suggestions,
	          tempStruct, tempLast, tempPrev, tempAsBits, commatize(bytesNeeded),
	          commatize(mallocLimit));
	return NULL; // (doesn't get here)
	}

overflow_last:
	suicidef ("in new_position_table(), last[] array size (%s) exceeds allocation limit of %s;"
	          suggestions,
	          commatize(bytesStruct), commatize(mallocLimit));
	return NULL; // (doesn't get here)

overflow_prev:
	suicidef ("in new_position_table(), prev[] array size (%s) exceeds allocation limit of %s;"
	          suggestions,
	          commatize(bytesPrev), commatize(mallocLimit));
	return NULL; // (doesn't get here)

overflow_as_bits:
	suicidef ("in new_position_table(), asBits[] array size (%s) exceeds allocation limit of %s;"
	          suggestions,
	          commatize(bytesAsBits), commatize(mallocLimit));
	return NULL; // (doesn't get here)
	}

//----------
//
// free_position_table--
//	De-allocate a position table.
//
//----------
//
// Arguments:
//	postable*	pt:	The position table to de-allocate.
//
// Returns:
//	(nothing)
//
//----------

void free_position_table (postable*	pt)
	{ free_if_valid ("free_position_table (table)", pt); }

//----------
//
// fetch_resolving_bits--
//	Fetch a 16-nucleotide word from a position table's packed representation of
//	sequence 1.
//
//----------
//
// Arguments:
//	postable*	pt:		The position table.  We assume pt->asBits is non-NULL.
//	unspos		pos1:	The position following the end of that word in the
//						.. sequence (origin-0, relative to pt->adjStart).  Note
//						.. that this is relative to the adjusted subinterval.
//
// Returns:
//	16 consecutive nucleotides from the sequence, as *32* bits.  From most
//	significant to least, bit pairs represent positions P-16 to P-1 (where P is
//	pos2).
//
//----------
//
// Notes:
//
// (1)	Schematic of the fetch.  Each x is a bit, and the bars show the boundary
//		of the 32-bit words in the asBits array.  ix is an index into the array,
//		while iy is an index into the sequence.
//
//		ix:  ..          5                                               6  ..
//		seq: .. xx xx xx|xx xx xx xx xx xx xx xx xx xx xx xx xx xx xx xx|xx ..
//		iy:  .. 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 ..
//		pos1=94:                                                    *
//		result:   [xx xx|xx xx xx xx xx xx xx xx xx xx xx xx xx xx]
//
//----------

#define wordSize (8*sizeof(pt->asBits[0]))	// (must be 32 or the return type
#define halfSize (wordSize/2)				//  .. is wrong!)

u32 fetch_resolving_bits
   (postable*	pt,
	unspos		pos1)
	{
	unspos		ix;
	int			shift;
	u32			seqBits;

	// split pos1 into array index and bit position

	ix   =  pos1 / halfSize;
	pos1 %= halfSize;

	// if bit position is zero we just fetch and return

	if (ix == 0) seqBits = 0;
	        else seqBits = pt->asBits[ix-1];

	if (pos1 == 0) return seqBits;

	// otherwise we have to shift and bring in more bits from the next
	// array word
	//
	// when wordSize=32, halfSize=16, and
	//		pos1=15 gives a shift of 2
	//		pos1=1  gives a shift of 30

	shift = (int) (2*(halfSize-pos1));

	return (seqBits        << (wordSize-shift))
	     + (pt->asBits[ix] >> shift);
	}

//----------
//
// position_is_in_table--
//	Determine if a position has a word stored in a position table.
//
//----------
//
// Arguments:
//	postable*	pt:			The position table.
//	unspos		position:	The position following the end of the word in the
//							.. sequence (origin-0, and divided by step).  Note
//							.. that this is relative to the *sequence*, and not
//							.. to the subinterval defined by pt->start,pt->end.
//							.. We expect (but do not check) that this is an
//							.. exact multiple of pt->step.
//
// Returns:
//	(nothing)
//
//----------

static int position_is_in_table
   (postable*	pt,
	unspos		position)
	{
	// convert the position to a prev[] index;  note that we expect (but do not
	// check) that position and start are both exact multiples of step

	position = (position - pt->adjStart) / pt->step;

	// see if that position is part of any list

	return (pt->prev[position] != 0);
	}

//----------
//
// add_word--
//	Add a word/position pair to a position table.
//
//----------
//
// Arguments:
//	postable*	pt:			The position table to add to.
//	u32			word:		The word to add.
//	unspos		position:	The position following the end of that word in the
//							.. sequence (origin-0, and divided by step).  Note
//							.. that this is relative to the *sequence*, and not
//							.. to the subinterval defined by pt->start,pt->end.
//							.. We expect (but do not check) that this is an
//							.. exact multiple of pt->step.
//
// Returns:
//	(nothing)
//
//----------

static void add_word
   (postable*	pt,
	u32			word,
	unspos		position)
	{
	u32			step = pt->step;
	unspos		oldLast;

	// convert the position to a prev[] index;  note that we expect (but do not
	// check) that position and start are both exact multiples of step

	position = (position - pt->adjStart) / step;

	// add the node to the front of the appropriate list

	oldLast = pt->last[word];
	if (oldLast == 0) pt->prev[position] = noPreviousPos;// was empty => end-of-list
	             else pt->prev[position] = oldLast;      // not empty => prepend
	pt->last[word] = position;

	// track some stats

	pos_table_count_stat  (wordsInTable);
	dbg_timing_count_stat (wordsInTable);

	if (oldLast == 0)						// (first occurence of this word)
		{
		pos_table_count_stat (wordsPresent);
		pos_table_count_stat (singletonWords);
		}
	else if (pt->prev[oldLast] == noPreviousPos)// (second occurence of this word)
		{
		pos_table_uncount_stat (singletonWords);
		}

	// debug

	if (pos_table_dbgShowWords)
		{
		posdumper dump = (pt->dump != NULL)? pt->dump
		                                   : (posdumper) dump_word_position;

		printf  ("adding ");
		(*dump) (stdout, pt, posdump_word,     word);
		printf  ("/");
		(*dump) (stdout, pt, posdump_position, position);
		printf  (" to table, prev is " unsposFmt "\n", pt->prev[position]);
		}

	}

//----------
//
// remove_word--
//	Remove a word/position pair from a position table.
//
//----------
//
// Arguments:
//	postable*	pt:			The position table to remove from.
//	u32			word:		The word to remove.
//	unspos		position:	The position following the end of that word in the
//							.. sequence (origin-0, and divided by step).  Note
//							.. that this is relative to the *sequence*, and not
//							.. to the subinterval defined by pt->start,pt->end.
//							.. We expect (but do not check) that this is an
//							.. exact multiple of pt->step.
//
// Returns:
//	(nothing)
//
//----------

static void remove_word
   (postable*	pt,
	u32			word,
	unspos		position)
	{
	u32			step = pt->step;
	unspos		pos, predPos;
	posdumper	dump = NULL;

	if (pos_table_dbgShowWords)
		{
		dump = (pt->dump != NULL)? pt->dump
		                         : (posdumper) dump_word_position;

		printf  ("removing ");
		(*dump) (stdout, pt, posdump_word,     word);
		printf  ("/");
		(*dump) (stdout, pt, posdump_position, position);
		printf  (" from table");
		}

	// convert the position to a prev[] index;  note that we expect (but do not
	// check) that position and start are both exact multiples of step

	position = (position - pt->adjStart) / step;

	// make sure the list for this word isn't empty

	if (pt->last[word] == 0)
		{
		if (pos_table_dbgShowWords) printf (" (list was empty)\n");
		return;
		}

	// find this position in the list

	predPos = noPreviousPos;
	for (pos=pt->last[word] ; pos!=noPreviousPos ; pos=pt->prev[pos])
		{
		if (pos == position) break;
		predPos	= pos;
		}

	if (pos != position)			// this position wasn't in the list
		{
		if (pos_table_dbgShowWords) printf (" (not found in list)\n");
		return;
		}

	// remove this position from the list

	if (predPos != noPreviousPos) 	// this position was *not* first in list
		pt->prev[predPos] = pt->prev[pos];
	else							// this position *was* the first in list
		{
		if (pt->prev[pos] == noPreviousPos) // the list is now empty
			pt->last[word] = 0;
		else
			pt->last[word] = pt->prev[pos];
		}

	pt->prev[pos] = 0;				// indicate position is no longer in table

	// track some stats

	pos_table_count_stat (wordsRemovedFromTable);

	if (pos_table_dbgShowWords)
		{
		if (predPos != noPreviousPos)	// this position was *not* first in list
			printf (", prev[" unsposFmt "] <- " unsposFmt "\n",
			        predPos, pt->prev[predPos]);
		else
			{
			printf (", last[");
			(*dump) (stdout, pt, posdump_word, word);
			printf ("] <- " unsposFmt "\n", pt->last[word]);
			}
		}
	}

//----------
//
// dump_position_table--
//	Dump the contents of a single position table.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	postable*	pt:				The position table to print.
//	seed*		hitSeed:		The seed-word the table was built for.  This is used
//								.. solely to unpack the table indexes.  This can be
//								.. NULL.
//	int			showPositions:	true => show a list of positions for each entry.
//	int			showCounts:		true => show count for each entry
//
// Returns:
//	(nothing)
//
//----------

// $$$ this needs to be updated to use the bottleneck alphabet rather than
// $$$ .. assuming ACGT

void dump_position_table
   (FILE*		f,
	postable*	pt,
	seed*		hitSeed,
	int			showPositions,
	int			showCounts)
	{
	posdumper	dump = (pt->dump != NULL)? pt->dump
				                         : (posdumper) dump_word_position;
	unspos		adjStart = pt->adjStart;
	u32			step     = pt->step;
	u32			w;
	unspos		pos;
	unspos		count;
	char*		s;

	for (w=0 ; w<pt->wordEntries ; w++)
		{
		if (pt->last[w] == 0) continue;
		if (pt->last[w] == noPreviousPos) continue;
		(*dump) (f, pt, posdump_index, w);

		if (hitSeed == NULL)
			fprintf (f, ":");
		else
			{
			s = seed_packed_to_string (hitSeed, w);
			fprintf (f, "/%s:", s);
			}

		if (showCounts)
			{
			count = 0;
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				count++;
			fprintf (f, " " unsposFmt, count);
			}

		if (showPositions)
			{
			fprintf (f, " ");
			pos = pt->last[w];
			(*dump) (f, pt, posdump_position, adjStart+step*pos);
			for (pos=pt->prev[pos] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				{
				fprintf (f, ",");
				(*dump) (f, pt, posdump_position, adjStart+step*pos);
				}
			}

		fprintf (f, "\n");
		}

	}

//----------
//
// dump_word_position--
//	Dump a word/position pair.
//
//----------
//
// Arguments:
//	(see posdumper description in pos_table.h)
//	The pt->dumpInfo field is expected to contain an int* pointing to the
//	.. word length
//
// Returns:
//	(nothing)
//
//----------

static void dump_word_position
   (FILE*		f,
	postable*	pt,
	int			field,
	u64			fieldVal)
	{
	int			wordLen = *(int*) pt->dumpInfo;

	switch (field)
		{
		default:
			break;
		case posdump_index:
			fprintf (f, "%0*X", (pt->wordBits+3)/4, (u32) fieldVal);
			break;
		case posdump_index_space:
			fprintf (f, "%*s", ((pt->wordBits+3)/4), "");
			break;
		case posdump_word:
			fprintf (f, "%s", bits_to_nuc_string (fieldVal, wordLen));
			break;
		case posdump_word_space:
			fprintf (f, "%*s", wordLen, "");
			break;
		case posdump_position:
			fprintf (f, unsposFmt, (unspos) fieldVal);
			break;
		}
	}

//----------
//
// dump_seed_position--
//	Dump a seed/position pair.
//
//----------
//
// Arguments:
//	(see posdumper description in pos_table.h)
//
// Returns:
//	(nothing)
//
//----------

static void dump_seed_position
   (FILE*		f,
	postable*	pt,
	int			field,
	u64			fieldVal)
	{
	dumpinfo*	di = pt->dumpInfo;
	seed*		hitSeed = di->seed;

	switch (field)
		{
		default:
			break;
		case posdump_index:
			fprintf (f, "%0*X", (pt->wordBits+3)/4, (u32) fieldVal);
			break;
		case posdump_index_space:
			fprintf (f, "%*s", (pt->wordBits+3)/4, "");
			break;
		case posdump_word:
			fprintf (f, "%s", seed_packed_to_string (hitSeed, fieldVal));
			break;
		case posdump_word_space:
			fprintf (f, "%*s", hitSeed->length, "");
			break;
		case posdump_position:
			fprintf (f, unsposFmt, (unspos) fieldVal);
			break;
		}
	}

//----------
//
// dump_quantum_seed_position--
//	Dump a seed/position pair for a quantum sequence.
//
//----------
//
// Arguments:
//	(see posdumper description in pos_table.h)
//	The pt->dumpInfo field is expected to contain a seed* pointing to the seed
//
// Returns:
//	(nothing)
//
//----------

static void dump_quantum_seed_position
   (FILE*		f,
	postable*	pt,
	int			field,
	u64			fieldVal)
	{
	dumpinfo*	di = pt->dumpInfo;
	seed*		hitSeed        = di->seed;
	u8*			bitsToAlphabet = di->bitsToAlphabet;
	char*		s;

	switch (field)
		{
		default:
			break;
		case posdump_index:
			fprintf (f, "%0*X", (pt->wordBits+3)/4, (u32) fieldVal);
			break;
		case posdump_index_space:
			fprintf (f, "%*s", (pt->wordBits+3)/4, "");
			break;
		case posdump_word:
			s = seed_packed_to_string2 (hitSeed, fieldVal, NULL, bitsToAlphabet);
			if    (*s != 0) fprintf (f,  "%02X", *(s++));
			while (*s != 0) fprintf (f, " %02X", *(s++));
			break;
		case posdump_word_space:
			fprintf (f, "%*s", hitSeed->length, "");
			break;
		case posdump_position:
			fprintf (f, unsposFmt, (unspos) fieldVal);
			break;
		}
	}


//----------
//
// count_position_table--
//	Count the number of positions in a table.
//
//----------
//
// Arguments:
//	postable*	pt:		The position table to count.
//
// Returns:
//	The number of words in the table.
//
//----------

unspos count_position_table
   (postable*	pt)
	{
	u32			w;
	unspos		pos;
	unspos		count;

	count = 0;
	for (w=0 ; w<pt->wordEntries ; w++)
		{
		if (pt->last[w] == 0) continue;
		for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
			count++;
		}

	return count;
	}

//----------
//
// limit_position_table--
//	Remove any words from a table that occur too frequently.
//
//----------
//
// Arguments:
//	postable*	pt:			The position table to modify.
//	u32			limit:		Words occurring more often than this are removed
//							.. from the table.
//	u32			maxChasm:	The maximum length of an interval of discarded
//							.. seed word positions that will be tolerated.
//							.. Some seed word positions may be preotected from
//							.. removal so that no interval will exceed this.
//							.. The value zero indicates that there is no such
//							.. limit.
//
// Returns:
//	(nothing)
//
//----------

static void breakup_chasm (char* protected, unspos startPos, unspos endPos,
                           unspos maxChasm);

void limit_position_table
   (postable*	pt,
	u32			_limit,
	unspos		maxChasm)
	{
	u32			w;
	unspos		pos, next;
	unspos		count;
	unspos		limit = _limit;
	char*		protected = NULL;

	pos_table_set_stat (wordCountLimit,    _limit);
	pos_table_set_stat (maxWordCountChasm, maxChasm);

	maxChasm /= pt->step;		// (convert maxChasm into the step realm)

	//////////
	// if we have a limit to the length of discard intervals, create a list of
	// "protected" positions
	//////////

	if (maxChasm > 0)
		{
		size_t	bytesNeeded;
		int		inChasm;
		unspos	chasmStart = 0; // (placation assignment)

		// create a list of position marks;  at this point all positions are
		// marked as "unprotected"

		bytesNeeded = pt->allocPrev * sizeof(char);
		protected = (char*) zalloc_or_die ("protected seed word positions", bytesNeeded);

		// scan position table and mark any positions that we intend to discard;
		// such positions will *all* temporarily be marked as "protected"

		for (w=0 ; w<pt->wordEntries ; w++)
			{
			if (pt->last[w] == 0) continue;
	
			count = 0;
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				count++;
			if (count <= limit) continue;
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				protected[pos] = true;
			}

		// scan marks and mark some positions in long intervals as protected

		inChasm = false;
		for (pos=0 ; pos<pt->allocPrev ; pos++)
			{
			if (protected[pos])
				{
				if (!inChasm)
					{ chasmStart = pos;  inChasm = true; }
				protected[pos] = false; // (breakup_chasm will set it back to
				continue;				//  .. true if necessary)
				}
			if (!inChasm)
				continue;
			inChasm = false;
			if (pos - chasmStart > maxChasm)
				breakup_chasm (protected, chasmStart, pos, maxChasm);
			}

		if ((inChasm) && (pos - chasmStart >= maxChasm))
			breakup_chasm (protected, chasmStart, pos, maxChasm);
		}

	//////////
	// dump the positions that will be limited
	//////////

	if (pos_table_dbgShowDiscards)
		{
		unspos    adjStart = pt->adjStart;
		u32       step     = pt->step;
		posdumper dump     = (pt->dump != NULL)? pt->dump
											   : (posdumper) dump_word_position;
		char*     s;
		unspos    numWords, numDiscarded;

		numWords = numDiscarded = 0;

		for (w=0 ; w<pt->wordEntries ; w++)
			{
			if (pt->last[w] == 0) continue;

			count = 0;
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				count++;
			numWords += count;
			if (count <= limit) continue;
			if (maxChasm > 0)
				{
				for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
					{ if (protected[pos]) count--; }
				}
			numDiscarded += count;
			}

		fprintf (stderr, "discarding %s/%s (%.2f%%) for maxwordcount=%d\n",
		                 commatize(numDiscarded), commatize(numWords),
		                 100.0*numDiscarded/numWords, _limit);

		for (w=0 ; w<pt->wordEntries ; w++)
			{
			if (pt->last[w] == 0) continue;

			count = 0;
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				count++;
			if (count <= limit) continue;

			(*dump) (stderr, pt, posdump_index, w);

			if (pos_table_dbgSeed == NULL)
				fprintf (stderr, ":");
			else
				{
				s = seed_packed_to_string (pos_table_dbgSeed, w);
				fprintf (stderr, "/%s:", s);
				}

			fprintf (stderr, " ");
			pos = pt->last[w];
			(*dump) (stderr, pt, posdump_position, adjStart+step*pos);
			for (pos=pt->prev[pos] ; pos!=noPreviousPos ; pos=pt->prev[pos])
				{
				fprintf (stderr, ",");
				(*dump) (stderr, pt, posdump_position, adjStart+step*pos);
				if ((maxChasm > 0) && (protected[pos]))
					fprintf (stderr, "*");
				}
			fprintf (stderr, "\n");
			}
		}

	//////////
	// discard positions from the table
	//////////

	for (w=0 ; w<pt->wordEntries ; w++)
		{
		if (pt->last[w] == 0) continue;

		count = 0;
		for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
			count++;
		if (count <= limit) continue;

		pos_table_add_stat (discardedWords, count);

		if (maxChasm == 0)
			{
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=next)
				{ next = pt->prev[pos];  pt->prev[pos] = noPreviousPos; }
			pt->last[w] = noPreviousPos;
			}
		else
			{
			unspos* pred;

			pred = &pt->last[w];
			for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=next)
				{
				next = pt->prev[pos];
				if (protected[pos])
					pred = &pt->prev[pos];
				else
					{
					*pred = next;
					pt->prev[pos] = noPreviousPos;
					}
				}
			}
		}

	//////////
	// erase any marks we made for the purpose of limiting discard intervals
	//////////

	free_if_valid ("protected seed word positions", protected);
	}


// breakup_chasm-- mark enough points in an interval to meet the maximum-chasm
//   criterion.  The algorithm used is similar to Brensenham's line drawing
//   algorithm, starting at position -1/2, stepping along the interval in steps
//   of N/D, truncating to an integer, and marking that position.

static void breakup_chasm
   (char*	protected,
	unspos	startPos,
	unspos	endPos,
	unspos	maxChasm)
	{
	unspos	pos, len, markNum;
	s64		numer;
	u64		denom;

	len = endPos - startPos;
	denom = 1 + (len / (maxChasm+1));	// (number of sub-intervals)
	numer = (denom/2) - denom;			// (intentionally wraps to 'negative')
	for (markNum=1 ; markNum<denom ; markNum++)
		{
		numer += len + 1;				// (numer is no longer 'negative')
		if (numer < 0)
			suicide ("internal error, in breakup_chasm()");
		pos   =  numer / denom;
		protected[startPos + pos] = true;
		}

	}

//----------
//
// find_position_table_limit--
//	Determine the word count limit that corresponds to a specified fraction of
//	low-occurring seed word positions.
//
//----------
//
// Arguments:
//	postable*	pt:		The position table.
//	float		keep:	A lower bound on the fraction of seed word positions
//						.. that are to be kept.  This value is between 0 and 1.
//
// Returns:
//	A limit, suitable for limit_position_table.
//
//----------

#define maxPossibleLimit ((u32) -1)

u32 find_position_table_limit
   (postable*	pt,
	float		keep)
	{
	poscount*	posDist, *pd;
	unspos		numPositions, minToKeep;
	unspos		limit;

	posDist = position_table_count_distribution (pt);

	// determine how many seed word positions we are required to keep

	numPositions = 0;
	for (pd=posDist ; pd->occurrences!=0 ; pd++)
		numPositions += pd->count * pd->occurrences;

	minToKeep = (unspos) ceil (numPositions * keep);

	// scan the list, counting up positions until we meet or exceed the
	// requirement (note that the posDist entries are in order of increasing
	// count)

	limit = 0;
	for (pd=posDist ; pd->occurrences!=0 ; pd++)
		{
		if (pd->count * pd->occurrences >= minToKeep)
			{ limit = pd->count;  break; }
		minToKeep -= pd->count * pd->occurrences;
		}

	free_if_valid ("seed word position counts distribution",  posDist);

	if (limit > maxPossibleLimit) return maxPossibleLimit;  // (bloody unlikely)
	                         else return limit;
	}

//----------
//
// position_table_count_distribution--
//	Determine the distribution of occurence counts in a table.
//
//----------
//
// Arguments:
//	postable*	pt:	The position table.
//
// Returns:
//	A pointer to the newly allocated count distribution, which the caller will
//	have to dispose of eventually (using free()).  This is an array ordered by
//	increasing count, and terminated by an entry with zero occurrences.
//
//----------

static int qIncreasingCount (const void* _a, const void* _b);
static int qIncreasingCount (const void* _a, const void* _b)
	{
	poscount* a = (poscount*) _a;
	poscount* b = (poscount*) _b;
	if      (a->count < b->count) return -1;
	else if (a->count > b->count) return  1;
	else                          return  0;
	}


poscount* position_table_count_distribution
   (postable*	pt)
	{
#ifndef noMemoryWrappers
	char*		idString = "position_table_count_distribution";
#endif // not noMemoryWrappers
	poscount*	posDist, *pd;
	int			countsAllocated, countsUsed;
	size_t		bytesNeeded;
	u32			w;
	unspos		pos, count;
	int			ix;

	// allocate an array to hold the distribution;  for hg18.chr1 with 13-mers,
	// the number of distict counts is 3,259

	countsAllocated = 3500;
	countsUsed      = 0;

	bytesNeeded = countsAllocated * sizeof(poscount);
	posDist = (poscount*) malloc_or_die (idString, bytesNeeded);

	posDist[0].occurrences = 0;	// (terminate list)

	// scan the table, count the locations for each word, and add to the
	// distibution
	// $$$ re-implement the list/search using heapsort

	for (w=0 ; w<pt->wordEntries ; w++)
		{
		if (pt->last[w] == 0) continue;

		count = 0;
		for (pos=pt->last[w] ; pos!=noPreviousPos ; pos=pt->prev[pos])
			count++;

		pd = NULL;
		for (ix=0 ; ix<countsUsed ; ix++)
			{ if (posDist[ix].count == count) { pd = &posDist[ix];  break; }}

		if (pd == NULL)
			{
			if (countsUsed+1 == countsAllocated)
				{
				countsAllocated += 1000;
				bytesNeeded = countsAllocated * sizeof(poscount);
				posDist = (poscount*) realloc_or_die (idString, posDist, bytesNeeded);
				}
			pd = &posDist[countsUsed++];
			pd->count = count;			// (note that pd->occurrences == 0)
			(pd+1)->occurrences = 0;	// (terminate list)
			}

		pd->occurrences++;
		}

	// sort by decreasing count

	qsort (posDist, countsUsed, sizeof(poscount), qIncreasingCount);

	return posDist;
	}

//----------
//
// pos_table_zero_stats--
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

void pos_table_zero_stats
   (void)
	{
	dbg_timing_set_stat (wordsInTable, 0);

#ifdef collect_stats

	// set 'em en masse to zero

	memset (&posTableStats, 0, sizeof(posTableStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// pos_table_show_stats,
// pos_table_show_stats_after--
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print the stats to.
//	postable*	pt:			The relevant position table (this can be NULL).
//
// Returns:
//	(nothing)
//
//----------

void pos_table_show_stats
   (arg_dont_complain(FILE*     f),
	arg_dont_complain(postable* pt))
	{
#ifdef collect_stats
	char		weight[10];
#endif // collect_stats

	dbg_timing_report_stat (wordsInTable, "DNA words in table");

#ifdef collect_stats

	if (f == NULL) return;

	if ((posTableStats.wordWeight % 2) == 0)
		sprintf (weight, "%d",   posTableStats.wordWeight / 2);
	else
		sprintf (weight, "%d.5", posTableStats.wordWeight / 2);

	fprintf (f, "word len or weight: %s\n",     weight);
	fprintf (f, "    possible words: %s\n",     commatize (posTableStats.wordSpace));
	fprintf (f, "-------------------\n");
	fprintf (f, "DNA words in table: %s\n",     commatize (posTableStats.wordsInTable));
//	if (pt != NULL)
//		fprintf (f, "          (actual): %s\n", commatize (tableCount));
	fprintf (f, "distinct DNA words: %s\n",     commatize (posTableStats.wordsPresent));
	fprintf (f, "   singleton words: %s\n",     commatize (posTableStats.singletonWords));
	if (posTableStats.discardedWords > 0)
		{
		fprintf (f, "        kept words: %s\n", commatize (posTableStats.wordsInTable - posTableStats.discardedWords));
		fprintf (f, "   discarded words: %s\n", commatize (posTableStats.discardedWords));
		}
	if (posTableStats.wordSpace > 0)
		fprintf (f, " distinct/possible: %.2f%%\n", (((u64) 100)*posTableStats.wordsPresent) / (float) posTableStats.wordSpace);
	if (posTableStats.wordsInTable > 0)
		{
		fprintf (f, "    distinct/words: %.2f%%\n", (((u64) 100)*posTableStats.wordsPresent)   / (float) posTableStats.wordsInTable);
		fprintf (f, "   singleton/words: %.2f%%\n", (((u64) 100)*posTableStats.singletonWords) / (float) posTableStats.wordsInTable);
		if (posTableStats.discardedWords > 0)
			{
			fprintf (f, "        kept/words: %.2f%%", (((u64) 100)*(posTableStats.wordsInTable-posTableStats.discardedWords)) / (float) posTableStats.wordsInTable);
			fprintf (f, " (word count limit = %u)\n", posTableStats.wordCountLimit);
			fprintf (f, "        (max chasm = %u)\n", posTableStats.maxWordCountChasm);
			fprintf (f, "   discarded/words: %.2f%%\n", (((u64) 100)*posTableStats.discardedWords) / (float) posTableStats.wordsInTable);
			fprintf (f, "   protected/words: %.2f%%\n", (((u64) 100)*posTableStats.protectedWords) / (float) posTableStats.wordsInTable);
			}
		}
	fprintf (f, "      bases parsed: %s\n",     commatize (posTableStats.basesParsed));
	fprintf (f, "-------------------\n");

#endif // collect_stats
	}


void pos_table_show_stats_after
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats

	if (f == NULL) return;

	fprintf (f, "  intervals masked: %s\n", commatize (posTableStats.intervalsMasked));
	fprintf (f, "masked i'val bases: %s\n", commatize (posTableStats.maskedIntervalBases));
	if (posTableStats.intervalsMasked > 0)
		fprintf (f, "bases/masked i'val: %.1f\n", ((float) posTableStats.maskedIntervalBases) / posTableStats.intervalsMasked);
	fprintf (f, " DNA words removed: %s\n", commatize (posTableStats.wordsRemovedFromTable));
	fprintf (f, " mask bases parsed: %s\n", commatize (posTableStats.maskBasesParsed));
	fprintf (f, "-------------------\n");

#endif // collect_stats
	}

