//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: identity_dist.c
//
//----------
//
// identity_dist--
//	Support for collecting the percent identity distribution from alignments.
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
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

#define  identity_dist_owner	// (make this the owner of its globals)
#include "identity_dist.h"		// interface to this module

// debugging defines

//#define snoopCountSubs		// if this is defined, extra code is added to
								// .. count_substitutions to track its operation

//#define snoopIdentityCalc		// if this is defined, extra code is added to
//								// .. to aid debugging of alignment_identity

//----------
//
// private global data
//
//----------
static int statsActive = false;

static unspos identityCount   [numIdentityBins+1];
static possum identityCoverage[numIdentityBins+1];

//----------
//
// filter_aligns_by_identity--
//	Filter a list of alignments, removing any alignment that has percent
//	identity outside of a specified range.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	alignel*	alignList:		The list of alignments to operate upon.
//	float		minIdentity,	The range of percent identity in alignments that
//				maxIdentity		.. we will *keep*.  These are values between 0
//								.. and 1.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
// Notes:
//	(1)	Memory for alignments that don't make the cut is deallocated here.
//	(2)	The returned list of alignments is in the same order as the incoming
//		.. list.
//	(3) Identity is counted over all aligned bases;  gaps are not counted.
//		$$$ consider whether gaps should be included (perhaps make this a
//		$$$ .. user option so inferz can experiment with it).
//----------

// $$$ There's an inherent inefficiency here if the user asks us to filter
// $$$ .. alignments by several idenitity-based stats.  We'd like to compute
// $$$ .. the identity stats only once, then either carry that in the alignment
// $$$ .. record, or apply all the filters in one pass through the loop.

alignel* filter_aligns_by_identity
   (seq*		seq1,
	seq*		seq2,
	alignel*	alignList,
	float		minIdentity,
	float		maxIdentity)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		alignment_identity (seq1, seq2, a, &numer, &denom);

		if ((denom == 0)
		 || (numer < denom * minIdentity)
		 || (numer > denom * maxIdentity))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_identity a->script", a->script);
			free_if_valid ("filter_aligns_by_identity a",         a);
			if (identity_dist_dbgShowIdentity)
				{
				// nota bene: positions written as 1-based
				printf ("discarding " unsposSlashSFmt " identity=" unsposSlashFmt "\n",
				        a->beg1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+",
				        a->beg2, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+",
				        numer, denom);
				}
			continue;
			}

		if (identity_dist_dbgShowIdentity)
			{
			// nota bene: positions written as 1-based
			printf ("keeping " unsposSlashSFmt " identity=" unsposSlashFmt "\n",
			        a->beg1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+",
			        a->beg2, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+",
			        numer, denom);
			}

		// this alignment is ok, add it to the end of the new list we're
		// building

		if (head == NULL) head = prev = a;
		             else { prev->next = a;  prev = a; }

		a->next = NULL;
		}

	return head;
	}

//----------
//
// alignment_identity--
//	Compute the identity fraction of an gapped alignment block.
//
//----------
//
// Arguments:
//	seq*		seq1:			The first sequence.
//	seq*		seq2:			The second sequence.
//	alignel*	a:				The alignment of interest.
//	unspos*		numer, denom:	Place to return the identity fraction.  Note
//								.. that the returned denominator might be zero.
//
// Returns:
//	(nothing)
//
//----------

//=== stuff for snoopIdentityCalc ===

#ifndef snoopIdentityCalc
#define debugSnoopIdentityCalc_1 ;
#define debugSnoopIdentityCalc_2 ;
#endif // not snoopIdentityCalc

#ifdef snoopIdentityCalc

#define debugSnoopIdentityCalc_1                                              \
	fprintf (stderr, "alignment_identity("                                    \
	                 unsposSlashFmt ".." unsposSlashFmt " "                   \
	                 unsposFmt "x" unsposFmt ")\n",                           \
	                 a->beg1, a->beg2, a->end1, a->end2, width, height);

#define debugSnoopIdentityCalc_2                                              \
	fprintf (stderr, "  " unsposFmt " bp run from "                           \
	                 unsposSlashFmt ".." unsposSlashFmt "\n",                 \
	                 run, prevI, prevJ, i, j);

#endif // snoopIdentityCalc


//=== stuff for alignment_identity ===

void alignment_identity
   (seq*		seq1,
	seq*		seq2,
	alignel*	a,
	unspos*		_numer,
	unspos*		_denom)
	{
	unspos		beg1 = a->beg1;
	unspos		beg2 = a->beg2;
	unspos		height, width, i, j, prevI, prevJ;
	u32			opIx;
    unspos		run;
	u8			c1, c2;
	unspos		denom, matches;
	unspos		pairCount[4][4];

	height = a->end1 - beg1 + 1;
	width  = a->end2 - beg2 + 1;
	debugSnoopIdentityCalc_1

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		pairCount[c1][c2] = 0;

	denom = 0;
	opIx  = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		prevI = i;  prevJ = j;
		run = edit_script_run_of_subs (a->script, &opIx);
		i += run; j += run;
		debugSnoopIdentityCalc_2

		denom += count_substitutions (seq1, beg1-1+prevI,
		                              seq2, beg2-1+prevJ,
		                              run, pairCount);

		if ((i < height) || (j < width))
			edit_script_indel_len (a->script, &opIx, &i, &j);
		}

	if (denom == 0)
		{ *_numer = *_denom = 0;  return; }

	matches = 0;
	for (c1=0 ; c1<4 ; c1++)
		matches += pairCount[c1][c1];

	*_numer = matches;
	*_denom = denom;
	}

//----------
//
// filter_segments_by_identity--
//	Filter a table of segments, removing any segment that has percent identity
//	outside of a specified range.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	segtable*	st:				The segment table to operate upon.
//	float		minIdentity,	The range of percent identity in segments that
//				maxIdentity		.. we will *keep*.  These are values between 0
//								.. and 1.
//
// Returns:
//	(nothing)
//
//----------

void filter_segments_by_identity
   (seq*		seq1,
	seq*		seq2,
	segtable*	st,
	float		minIdentity,
	float		maxIdentity)
	{
	segment*	srcSeg, *dstSeg;
	unspos		numer, denom;

	if (st      == NULL) return;
	if (st->seg == NULL) return;

	for (dstSeg=srcSeg=st->seg ; ((u32)(srcSeg-st->seg))<st->len ; srcSeg++)
		{
		segment_identity (seq1, srcSeg->pos1,
		                  seq2, srcSeg->pos2, srcSeg->length,
		                  &numer, &denom);
		if ((denom == 0)
		 || (numer < denom * minIdentity)
		 || (numer > denom * maxIdentity))
			continue; // (unwanted segment, skip it)
		if (srcSeg != dstSeg) *dstSeg = *srcSeg;
		dstSeg++;
		}

	st->len = dstSeg - st->seg;
	}

//----------
//
// filter_segment_by_identity--
//	Filter a segment, reporting whether that segment has percent identity
//	outside of a specified range.
//
//----------
//
// Arguments:
//	seq*	seq1:			The first sequence.
//	unspos 	pos1:			The subsequence start position in seq1 (origin-0).
//	seq*	seq2:			The second sequence.
//	unspos 	pos2:			The subsequence start position in seq2 (origin-0).
//	unspos 	length:			The length of the subsequence.
//	float	minIdentity,	The range of percent identity in segments that
//			maxIdentity		.. we will *keep*.  These are values between 0
//							.. and 1.
//
// Returns:
//	true if the segment if outside the specified range (i.e. that it fails to
//	pass the filter, and should be discarded)
//
//----------

int filter_segment_by_identity
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	float	minIdentity,
	float	maxIdentity)
	{
	unspos	numer, denom;

	segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
	if ((denom == 0)
	 || (numer < denom * minIdentity)
	 || (numer > denom * maxIdentity))
		return true; // (unwanted segment, skip it)

	return false;
	}

//----------
//
// segment_identity--
//	Compute the identity fraction of an ungapped alignment segment.
//
//----------
//
// Arguments:
//	seq*	seq1:			The first sequence.
//	unspos 	pos1:			The subsequence start position in seq1 (origin-0).
//	seq*	seq2:			The second sequence.
//	unspos 	pos2:			The subsequence start position in seq2 (origin-0).
//	unspos 	length:			The length of the subsequence.
//	unspos*	numer, denom:	Place to return the identity fraction.  Note that
//							.. the returned denominator might be zero.
//
// Returns:
//	(nothing)
//
//----------

void segment_identity
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	unspos*	_numer,
	unspos*	_denom)
	{
	u8		c1, c2;
	unspos	denom, matches;
	unspos	pairCount[4][4];

	// count substitutions and see if we pass our percent identity filter

	for (c1=0 ; c1<4 ; c1++)
			for (c2=0 ; c2<4 ; c2++)
		pairCount[c1][c2] = 0;

	denom = count_substitutions (seq1, pos1, seq2, pos2, length, pairCount);
	if (denom == 0)
		{ *_numer = *_denom = 0;  return; }

	matches = 0;
	for (c1=0 ; c1<4 ; c1++)
		matches += pairCount[c1][c1];

	*_numer = matches;
	*_denom = denom;
	}

//----------
//
// count_substitutions--
//	Count the number of each  type of base substitution in two subsequences.
//
//----------
//
// Arguments:
//	seq*	seq1:			The first sequence.
//	unspos 	pos1:			The subsequence start position in seq1 (origin-0).
//	seq*	seq2:			The second sequence.
//	unspos 	pos2:			The subsequence start position in seq2 (origin-0).
//	unspos 	length:			The length of the subsequence.  Note that this may
//							.. be zero.
//	unspos	count[4][4]:	Array in which to accumulate the substitutions.
//							.. This is indexed by [c1][c2], where c1 and c2 are
//							.. fromsequence 1 and 2, respectively, and code for
//							.. nucleotides as per bits_to_nuc[].  Note that
//							.. wew *accumulate* into this array;  we don't clear
//							.. it first.
//
// Returns:
//	The number of new substitutions or matches we've counted.
//
//----------
//
// Note:  Masked (lowercase) bp are counted the same as unmasked (uppercase),
//        but illegal values like 'N' or '-' are completely ignored.
//
//----------

//=== stuff for snoopCountSubs ===

#ifndef snoopCountSubs
#define debugSnoopCountSubs_1 ;
#endif // not snoopCountSubs

#ifdef snoopCountSubs

#define debugSnoopCountSubs_1                                                 \
	{                                                                         \
	fprintf (stderr, "count_substitutions(" unsposSlashFmt " #" unsposFmt "): ", \
	                 pos1, pos2, length);                                     \
	for (ix=0 ; ix<length ; ix++) fprintf (stderr, "%c", s1[ix]);             \
	fprintf (stderr, " ");                                                    \
	for (ix=0 ; ix<length ; ix++) fprintf (stderr, "%c", s2[ix]);             \
	fprintf (stderr, "\n");                                                   \
	}

#endif // snoopCountSubs


//=== count_substitutions ===

unspos count_substitutions
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	unspos	count[4][4])
	{
	u8*		s1 = seq1->v + pos1;
	u8*		s2 = seq2->v + pos2;
	s8		c1, c2;
	unspos	ix;
	unspos	denom = 0;

	if (length == 0)
		return 0;

	debugSnoopCountSubs_1

	for (ix=0 ; ix<length ; ix++)
		{
		c1 = nuc_to_bits[*(s1++)];
		c2 = nuc_to_bits[*(s2++)];
		if ((c1 >= 0) && (c2 >= 0))
			{ count[(u8)c1][(u8)c2]++;  denom++; }
		}

	return denom;
	}

//----------
//
// filter_aligns_by_match_count--
//	Filter a list of alignments, removing any alignment that has fewer matched
//	bases than a specified minimum.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	alignel*	alignList:		The list of alignments to operate upon.
//	u32			minMatchCount:	The minimum number of matched bases in
//								.. alignments that we will *keep*.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
// Notes:
//	(1)	Memory for alignments that don't make the cut is deallocated here.
//	(2)	The returned list of alignments is in the same order as the incoming
//		.. list.
//	(3) Match-count is counted over all aligned bases, and counts only matches
//		.. (not substitutions or gaps).
//----------

alignel* filter_aligns_by_match_count
   (seq*		seq1,
	seq*		seq2,
	alignel*	alignList,
	u32			minMatchCount)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		alignment_identity (seq1, seq2, a, &numer, &denom);

		if ((denom == 0) || (numer < minMatchCount))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_match_count a->script", a->script);
			free_if_valid ("filter_aligns_by_match_count a",         a);
			continue;
			}

		// this alignment is ok, add it to the end of the new list we're
		// building

		if (head == NULL) head = prev = a;
		             else { prev->next = a;  prev = a; }

		a->next = NULL;
		}

	return head;
	}

//----------
//
// filter_segments_by_match_count--
//	Filter a table of segments, removing any segment that has fewer matched
//	bases than a specified minimum.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:		The sequences.
//	segtable*	st:				The segment table to operate upon.
//	u32			minMatchCount:	The minimum number of matched bases in segments
//								.. that we will *keep*.
//
// Returns:
//	(nothing)
//
//----------

void filter_segments_by_match_count
   (seq*		seq1,
	seq*		seq2,
	segtable*	st,
	u32			minMatchCount)
	{
	segment*	srcSeg, *dstSeg;
	unspos		numer, denom;

	if (st      == NULL) return;
	if (st->seg == NULL) return;

	for (dstSeg=srcSeg=st->seg ; ((u32)(srcSeg-st->seg))<st->len ; srcSeg++)
		{
		segment_identity (seq1, srcSeg->pos1,
		                  seq2, srcSeg->pos2, srcSeg->length,
		                  &numer, &denom);
		if ((denom == 0) || (numer < minMatchCount))
			continue; // (unwanted segment, skip it)
		if (srcSeg != dstSeg) *dstSeg = *srcSeg;
		dstSeg++;
		}

	st->len = dstSeg - st->seg;
	}

//----------
//
// filter_segment_by_match_count--
//	Filter a segment, reporting whether that segment has fewer matched bases
//	than a specified minimum.
//
//----------
//
// Arguments:
//	seq*	seq1:			The first sequence.
//	unspos 	pos1:			The subsequence start position in seq1 (origin-0).
//	seq*	seq2:			The second sequence.
//	unspos 	pos2:			The subsequence start position in seq2 (origin-0).
//	unspos 	length:			The length of the subsequence.
//	u32		minMatchCount:	The minimum number of matched bases in segments
//							.. that we will *keep*.
//
// Returns:
//	true if the segment if outside the specified range (i.e. that it fails to
//	pass the filter, and should be discarded)
//
//----------

int filter_segment_by_match_count
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	u32		minMatchCount)
	{
	unspos	numer, denom;

	segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
	if ((denom == 0) || (numer < minMatchCount))
		return true; // (unwanted segment, skip it)

	return false;
	}

//----------
//
// filter_aligns_by_mismatch_count--
//	Filter a list of alignments, removing any alignment that has more mismatched
//	bases than a specified maximum.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:			The sequences.
//	alignel*	alignList:			The list of alignments to operate upon.
//	s32			maxMismatchCount:	The maximum number of mismatched bases in
//									.. alignments that we will *keep*.  Note
//									.. that this must not be negative.
//
// Returns:
//	A pointer to the list of remaining alignments.
//
//----------
// Notes:
//	(1)	Memory for alignments that don't make the cut is deallocated here.
//	(2)	The returned list of alignments is in the same order as the incoming
//		.. list.
//	(3) Mismatch-count is counted over all aligned bases, and counts only
//		.. mismatches (not matches or gaps).
//----------

alignel* filter_aligns_by_mismatch_count
   (seq*		seq1,
	seq*		seq2,
	alignel*	alignList,
	s32			maxMismatchCount)
	{
	alignel*	a, *next;
	alignel*	head, *prev;
	unspos		numer, denom;

	head = prev = NULL;
	for (a=alignList ; a!=NULL ; a=next)
		{
		next = a->next;

		alignment_identity (seq1, seq2, a, &numer, &denom);

		if ((denom == 0) || (denom-numer > (u32) maxMismatchCount))
			{ // (unwanted alignment, discard it)
			free_if_valid ("filter_aligns_by_mismatch_count a->script", a->script);
			free_if_valid ("filter_aligns_by_mismatch_count a",         a);
			continue;
			}

		// this alignment is ok, add it to the end of the new list we're
		// building

		if (head == NULL) head = prev = a;
		             else { prev->next = a;  prev = a; }

		a->next = NULL;
		}

	return head;
	}

//----------
//
// filter_segments_by_mismatch_count--
//	Filter a table of segments, removing any segment that has more mismatched
//	bases than a specified maximum.
//
//----------
//
// Arguments:
//	seq*		seq1, seq2:			The sequences.
//	segtable*	st:					The segment table to operate upon.
//	s32			maxMismatchCount:	The maximum number of matched bases in
//									.. segments that we will *keep*.  Note that
//									.. this must not be negative.
//
// Returns:
//	(nothing)
//
//----------

void filter_segments_by_mismatch_count
   (seq*		seq1,
	seq*		seq2,
	segtable*	st,
	s32			maxMismatchCount)
	{
	segment*	srcSeg, *dstSeg;
	unspos		numer, denom;

	if (st      == NULL) return;
	if (st->seg == NULL) return;

	for (dstSeg=srcSeg=st->seg ; ((u32)(srcSeg-st->seg))<st->len ; srcSeg++)
		{
		segment_identity (seq1, srcSeg->pos1,
		                  seq2, srcSeg->pos2, srcSeg->length,
		                  &numer, &denom);
		if ((denom == 0) || (denom-numer > (u32) maxMismatchCount))
			continue; // (unwanted segment, skip it)
		if (srcSeg != dstSeg) *dstSeg = *srcSeg;
		dstSeg++;
		}

	st->len = dstSeg - st->seg;
	}

//----------
//
// filter_segment_by_mismatch_count--
//	Filter a segment, reporting whether that segment has more mismatched bases
//	than a specified maximum.
//
//----------
//
// Arguments:
//	seq*	seq1:				The first sequence.
//	unspos 	pos1:				The subsequence start position in seq1 (origin-0).
//	seq*	seq2:				The second sequence.
//	unspos 	pos2:				The subsequence start position in seq2 (origin-0).
//	unspos 	length:				The length of the subsequence.
//	s32		maxMismatchCount:	The maximum number of matched bases in
//								.. segments that we will *keep*.  Note that
//								.. this must not be negative.
//
// Returns:
//	true if the segment if outside the specified range (i.e. that it fails to
//	pass the filter, and should be discarded)
//
//----------

int filter_segment_by_mismatch_count
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	s32		maxMismatchCount)
	{
	unspos	numer, denom;

	segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
	if ((denom == 0) || (denom-numer > (u32) maxMismatchCount))
		return true; // (unwanted segment, skip it)

	return false;
	}

//----------
//
// init_identity_dist_job--
//	Initialize percent identity distribution.
//
//----------

void init_identity_dist_job
   (arg_dont_complain(seq* seq1),
	arg_dont_complain(seq* seq2))
	{
	u32		bin;

	if (statsActive)
		suicide ("attempt to open a second identity distribution job");
	statsActive = true;

	for (bin=0 ; bin<=numIdentityBins ; bin++)
		{
		identityCount   [bin] = 0;
		identityCoverage[bin] = 0;
		}
	}

//----------
//
// print_identity_dist_job--
//	Print the percent identity distribution.
//
//----------

void print_identity_dist_job
   (FILE*	f)
	{
	static const u32 noBin = (u32) -1;
	u32		bin, minBin, maxBin;
	float	binCenter;

	if (!statsActive)
		suicide ("attempt to close a non-existent identity distribution job");

	minBin = maxBin = noBin;
	for (bin=0 ; bin<=numIdentityBins ; bin++)
		{
		if (identityCount[bin] == 0) continue;
		maxBin = bin;
		if (minBin == noBin) minBin = bin;
		}
	if (minBin == noBin) minBin = maxBin = numIdentityBins;

	if (minBin > 0)               minBin--;		// inferz likes to have an empty
	if (maxBin < numIdentityBins) maxBin++;		// .. bin before and after the
												// .. table

	for (bin=minBin ; bin<=maxBin ; bin++)
		{
		binCenter = bin / ((float) numIdentityBins);
		fprintf (f, identityBinFormat "\t" unsposFmt "\t" possumFmt "\n",
		            binCenter, identityCount[bin], identityCoverage[bin]);
		}

	statsActive = false;
	}

//----------
//
// identity_dist_from_align_list--
//	Collect percent identity distribution from a list of gapped alignments.
//
//----------
//
// Arguments:
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//
// Returns:
//	(nothing)
//
//----------

void identity_dist_from_align_list
   (alignel*	alignList,
	seq*		seq1,
	seq*		seq2)
	{
	alignel*	a;
	unspos		numer, denom;
	u32			bin;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		alignment_identity (seq1, seq2, a, &numer, &denom);
		bin = identity_bin (numer, denom);

		identityCount   [bin]++;
		identityCoverage[bin] += denom;
		}
	}

//----------
//
// identity_dist_from_match--
//	Collect percent identity distribution from a single ungapped alignment.
//
//----------
//
// Arguments:
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos2:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the match.
//
// Returns:
//	(nothing)
//
//----------

void identity_dist_from_match
   (seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	unspos	numer, denom;
	u32		bin;

	segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
	bin = identity_bin (numer, denom);

	identityCount   [bin]++;
	identityCoverage[bin] += denom;
	}

