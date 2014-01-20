//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: masking.c
//
//----------
//
// masking--
//	Provide dynamic masking for sequences that do not have their repeats
//	masked.
//
// The basic premise is that if anything in the first sequence aligns to a
// repeat in the second sequence, it will align to all copies of the repeat.
// Thus, positions that align frequently are probably repeat elements.
//
// For every position in sequence 1, we keep a count of how many times it has
// been part of an alignment (even if it's part of an indel).  When the count
// reaches the given threshold, we mark that position with an 'x' so that
// subsequent alignments will give it a large scoring penalty.
//
// This is most useful if the caller is comparing sequence 1 to many other
// sequences.  The first comparison(s) will get no advantage-- nothing will be
// masked.  But if enough sequences are compared, the later sequences will
// benefit from the masking.
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
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

#define  masking_owner			// (make this the owner of its globals)
#include "masking.h"			// interface to this module

//----------
//
// miscellaneous data
//
//----------

//#define noCallback			// if this is defined, mask_interval doesn't
								// bother to call the callback routine for
								// masked intervals;  this is to test whether
								// there is any speed advantage to removing
								// masked seeds from the position table;  note
								// that there are other good reasons to remove
								// masked seeds, regardless of speed

//----------
//
// prototypes for private functions
//
//----------

static unspos mask_interval (u8* fwd, u8* rev, unspos beg, unspos end,
                             census* cen,
                             void(*func) (unspos beg, unspos end, void* info),
                             void* info);

//----------
//
// new_census--
//	Create a new census of a given length.
//
//----------
//
// Arguments:
//	new_census	len:		The length of the sequence to be censused.
//	char		kind:		The type of count[] array to allocate (see the
//							.. definition of the census type in masking.h).
//	u32			maskThresh:	The count threshold.  Note that we don't range
//							.. check this to compare its validity vs. kind.
//
// Returns:
//	Pointer to the new census (failures are fatal);  caller must eventually
//	dispose of this memory with a call to free().
//
//----------

census* new_census
   (unspos	len,
	char	kind,
	u32		maskThresh)
	{
	census*	cen;
	size_t	bytesNeeded;

	if (strchr ("BWL", kind) == NULL) kind = 'B';
	if (len == 0) len = 1;

	if (kind == 'B')
		{
		bytesNeeded = sizeof(census) + (len*sizeof(u8));
		cen = (census*) zalloc_or_die ("new_census", bytesNeeded);
		cen->count8 = (u8*) (cen + 1);
		}
	else if (kind == 'W')
		{
		bytesNeeded = sizeof(census) + (len*sizeof(u16));
		cen = (census*) zalloc_or_die ("new_census", bytesNeeded);
		cen->count16 = (u16*) (cen + 1);
		}
	else // if (kind == 'L')
		{
		bytesNeeded = sizeof(census) + (len*sizeof(u32));
		cen = (census*) zalloc_or_die ("new_census", bytesNeeded);
		cen->count32 = (u32*) (cen + 1);
		}

	cen->len        = len;
	cen->kind       = kind;
	cen->maskThresh = maskThresh;

	return cen;
	}

//----------
//
// census_mask_segments--
//	Given a list of ungapped alignment segments, increment the count of all
//	aligned locations and for any that meet or exceed the threshold, mask them
//	from the sequence.
//
//----------
//
// Arguments:
//	segtable*	st:			The segments.
//	u8*			fwd, rev:	The sequence and its reverse.  rev may be NULL.
//	census*		cen:		The census to update.
//	void (*func):			A callback function to report masked intervals to.
//	        (unspos beg,	.. beg,end is the interval that is about to be
//	         unspos end,	.. masked (origin 1, inclusive).  info is passed
//	         void* info)	.. through from our caller.
//	void*		info:		Pass-thru for func.
//
// Returns:
//	A count of the number of bases masked.  This does not include any bases that
//	were previously masked.
//
//----------

unspos census_mask_segments
   (segtable*	st,
	u8*			fwd,
	u8*			rev,
	census*		cen,
	void		(*func)(unspos beg, unspos end, void* info),
	void*		info)
	{
	segment*	seg;
	u32			ix;
	unspos		beg, end, pos;
	unspos		count = 0;

	if (cen == NULL)
		return 0;

	for (ix=0,seg=st->seg ; ix<st->len ; ix++,seg++)
		{
		if ((seg->length > cen->len) || (seg->pos1 > cen->len - seg->length))
			suicide ("census_mask_segments, internal error");
		beg = seg->pos1;
		end = beg + seg->length;
		switch (cen->kind)
			{
			case 'B':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count8[pos]  < u8max)  cen->count8 [pos]++; }
				break;
			case 'W':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count16[pos] < u16max) cen->count16[pos]++; }
				break;
			case 'L':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count32[pos] < u32max) cen->count32[pos]++; }
				break;
			}
		if (cen->maskThresh > 0)
			count += mask_interval (fwd, rev, beg+1, end, cen, func, info);
		}

	masking_add_stat (maskedBases, count);
	return count;
	}

//----------
//
// census_mask_aligns--
//	Given a list of alignments, increment the count of all aligned locations
//	and for any that meet or exceed the threshold, mask them from the sequence.
//
// Note that the entire aligned interval is counted/masked, even indels.
//
//----------
//
// Arguments:
//	alignel*	alignList:	The alignments.
//	u8*			fwd, rev:	The sequence and its reverse.  rev may be NULL.
//	census*		cen:		The census to update.
//	void (*func):			A callback function to report masked intervals to.
//	        (unspos beg,	.. beg,end is the interval that is about to be
//	         unspos end,	.. masked (origin 1, inclusive).  info is passed
//	         void* info)	.. through from our caller.
//	void*		info:		Pass-thru for func.
//
// Returns:
//	A count of the number of bases masked.  This does not include any bases that
//	were previously masked.
//
//----------

unspos census_mask_aligns
   (alignel*	alignList,
	u8*			fwd,
	u8*			rev,
	census*		cen,
	void		(*func)(unspos beg, unspos end, void* info),
	void*		info)
	{
	alignel*	a;
	unspos		beg, end, pos;
	unspos		count = 0;

	if (cen == NULL)
		return 0;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		if (a->beg1 < 1)
			suicide ("census_mask_aligns, internal error");
		if (a->end1 > cen->len)
			suicide ("census_mask_aligns, internal error");
		beg = a->beg1 - 1;
		end = a->end1;
		switch (cen->kind)
			{
			case 'B':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count8[pos]  < u8max)  cen->count8 [pos]++; }
				break;
			case 'W':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count16[pos] < u16max) cen->count16[pos]++; }
				break;
			case 'L':
				for (pos=beg ; pos<end ; pos++)
					{ if (cen->count32[pos] < u32max) cen->count32[pos]++; }
				break;
			}
		if (cen->maskThresh > 0)
			count += mask_interval (fwd, rev, beg+1, end, cen, func, info);
		}

	masking_add_stat (maskedBases, count);
	return count;
	}

//----------
//
// mask_interval--
//	Given an interval, mask any position that meets or exceeds the count
//	threshold.
//
//----------
//
// Arguments:
//	u8*			fwd, rev:	The sequence and its reverse.  rev may be NULL.
//	unspos		beg, end:	The interval to check.  Origin-1, inclusive.
//	census*		cen:		The census describing which bases to mask.
//	void (*func):			A callback function to report masked intervals to.
//	        (unspos beg,	.. beg,end is the interval that is about to be
//	         unspos end,	.. masked, see note 1 (origin 1, inclusive).  info
//	         void* info)	.. is passed through from our caller.
//	void*		info:		Pass-thru for func.
//
// Returns:
//	A count of the number of bases masked.  This does not include any bases that
//	were previously masked.
//
//----------
//
// Notes:
//	(1)	The intervals passed to the callback function are not the same as the
//		interval passed to this routine.  Only sub-intervals that reach the
//		masking threshold for the first time are passed to the callback.
//
//----------

#ifndef noCallback		// normal version, callback used

static unspos mask_interval
   (u8*		fwd,
	u8*		rev,
	unspos	beg,
	unspos	end,
	census*	cen,
	void	(*func)(unspos beg, unspos end, void* info),
	void*	info)
	{
	static const u32 noRun = (u32) -1;
	unspos	revLen, runBeg, pos, j;
	unspos	basesMasked = 0;
	u32		count;

	if (cen == NULL)
		return 0;

	if ((beg < 1) || (end > cen->len))
		suicide ("mask_interval, internal error");

	revLen = cen->len - 1;

	runBeg = noRun;
	for (pos=beg-1 ; pos<end ; pos++)
		{
		count = (cen->kind == 'B')? cen->count8 [pos]
			  : (cen->kind == 'W')? cen->count16[pos]
								  : cen->count32[pos];
		if ((cen->maskThresh > 0)
		 && (count >= cen->maskThresh)
		 && (dna_isupper(fwd[pos])))
			{ if (runBeg == noRun) runBeg = pos; }
		else if (runBeg != noRun)
			{
			func (runBeg+1, pos, info);
			for (j=runBeg ; j<pos ; j++)
				{
				fwd[j] = 'x';  basesMasked++;
				if (rev != NULL) rev[revLen-j] = 'x';
				}
			runBeg = noRun;
			}
		}

	if (runBeg != noRun)
		{
		func (runBeg+1, end, info);
		for (j=runBeg ; j<end ; j++)
			{
			fwd[j] = 'x';  basesMasked++;
			if (rev != NULL) rev[revLen-j] = 'x';
			}
		}

	return basesMasked;
	}

#endif // not noCallback


#ifdef noCallback		// test version, callback not used

static int mask_interval
   (u8*		fwd,
	u8*		rev,
	unspos	beg,
	unspos	end,
	census*	cen,
	void	(*func)(unspos beg, unspos end, void* info),
	void*	info)
	{
	unspos	revLen, pos;
	unspos	basesMasked = 0;
	u32		count;

	if ((cen == NULL)
	 || (cen->maskThresh == 0))
		return 0;

	if ((beg < 1) || (end > cen->len))
		suicide ("mask_interval, internal error");

	revLen = cen->len - 1;
	for (pos=beg-1 ; pos<end ; pos++)
		{
		count = (cen->kind == 'B')? cen->count8 [pos]
			  : (cen->kind == 'W')? cen->count16[pos]
								  : cen->count32[pos];
		if (count >= cen->maskThresh)
			{
			if (dna_isupper(fwd[pos])) { fwd[pos] = 'x';  basesMasked++; }
			if (rev != NULL)             rev[revLen-pos] = 'x';
			}
		}

	return basesMasked;
	}

#endif // noCallback

//----------
//
// count_masked_bases--
//	Report masked bases in a sequence.  The masked bases could have been in the
//	sequence when it was input, or they could have been added during running of
//	the program (e.g. by dynamic masking).
//
//----------
//
// Arguments:
//	seq*	seq:			The sequence that has been masked.
//	int		maskChar:		The character to consider as a mask.  Normally this
//							.. is a character (e.g. 'X' or 'N').  However, the
//							.. value -1 means that we should mask by changing
//							.. to lowercase.
//
// Returns:
//	The number of bases masked.
//
//----------

unspos count_masked_bases
   (seq*	_seq,
	int		maskChar)
	{
	unspos	ix;
	unspos	numMasked = 0;
	int		isMasked;

	for (ix=0 ; ix<=_seq->len ; ix++)
		{
		if (maskChar >= 0)
			isMasked = (_seq->v[ix] == maskChar);
		else
			isMasked = ((_seq->v[ix] >= 'a') && (_seq->v[ix] <= 'z'));

		if (isMasked) numMasked++;
		}

	return numMasked;
	}

//----------
//
// report_census_intervals--
//	Report runs in a census of positions that meet or exceed the threshold.
//
//----------
//
// Arguments:
//	census*		cen:		The census to report.
//	void (*func):			The reporting function.  beg,end is the interval
//	        (unspos beg,	.. (origin 1, inclusive).  info is passed through
//	         unspos end,	.. from our caller.
//	         void* info)
//	void*		info:		Pass-thru for func.
//
// Returns:
//	The number of intervals reported.
//
//----------

unspos report_census_intervals
   (census*	cen,
	void	(*func)(unspos beg, unspos end, void* info),
	void*	info)
	{
	static const u32 noRun = (u32) -1;
	unspos	runBeg, pos;
	unspos	numIntervals = 0;
	u32		count;

	if (cen == NULL)
		return 0;

	runBeg = noRun;
	for (pos=0 ; pos<cen->len ; pos++)
		{
		count = (cen->kind == 'B')? cen->count8 [pos]
			  : (cen->kind == 'W')? cen->count16[pos]
								  : cen->count32[pos];
		if (count >= cen->maskThresh)
			{ if (runBeg == noRun) runBeg = pos; }
		else if (runBeg != noRun)
			{
			if (func != NULL) func (runBeg+1, pos, info);
			numIntervals++;
			runBeg = noRun;
			}
		}

	if (runBeg != noRun)
		{
		if (func != NULL) func (runBeg+1, cen->len, info);
		numIntervals++;
		}

	return numIntervals;
	}

//----------
//
// report_masked_intervals--
//	Report masked runs in a sequence.  The masked bases could have been in the
//	sequence when it was input, or they could have been added during running of
//	the program (e.g. by dynamic masking).
//
//----------
//
// Arguments:
//	seq*	seq:			The sequence that has been masked.
//	int		maskChar:		The character to consider as a mask.  Normally this
//							.. is a character (e.g. 'X' or 'N').  However, the
//							.. value -1 means that we should mask by changing
//							.. to lowercase.
//	void (*func):			The reporting function.  beg,end is the interval
//	        (unspos beg,	.. (origin 1, inclusive).  info is passed through
//	         unspos end,	.. from our caller.
//	         void* info)
//	void*		info:		Pass-thru for func.
//
// Returns:
//	The number of intervals reported.
//
//----------

unspos report_masked_intervals
   (seq*	_seq,
	int		maskChar,
	void	(*func)(unspos beg, unspos end, void* info),
	void*	info)
	{
	static const u32 noRun = (u32) -1;
	unspos	runBeg, ix;
	unspos	numIntervals = 0;
	int		isMasked;

	runBeg = noRun;
	for (ix=0 ; ix<=_seq->len ; ix++)
		{
		if (maskChar >= 0)
			isMasked = (_seq->v[ix] == maskChar);
		else
			isMasked = ((_seq->v[ix] >= 'a') && (_seq->v[ix] <= 'z'));

		if (isMasked)
			{ if (runBeg == noRun) runBeg = ix; }
		else if (runBeg != noRun)
			{
			if (func != NULL) func (runBeg+1, ix, info);
			numIntervals++;
			runBeg = noRun;
			}
		}

	if (runBeg != noRun)
		{
		if (func != NULL) func (runBeg+1, _seq->len, info);
		numIntervals++;
		}

	return numIntervals;
	}

//----------
// [[-- a report_census_intervals or report_masked_intervals callback function --]]
//
// print_masking_interval--
//	Write a masking interval to a file, in a format compatible with being read
//	by mask_sequence() or mask_sequence_keep().
// print_masking_interval_3--
//	Write a masking interval to a file, including sequence names.  This format
//	is NOT compatible with mask_sequence() or mask_sequence_keep().
//
// A typical masking file (resulting from several calls to this routine) looks
// like this (for print_masking_interval):
//
//	1527933 3184039
//	4165389 6877343
//	7374477 7902860
//
// or like this (for print_masking_interval_3):
//
//	chr1 1527933 3184039
//	chr1 4165389 6877343
//	chr3 7374477 7902860
//
// Each line describes a region to be masked.  Indexes are one-based, and
// inclusive on both ends.
//
//----------
//
// Arguments:
//	unspos beg, end:	The interval, in the target sequence, that is to be
//						.. printed.  Origin-1, inclusive.
//	void* info:			(really pmiInfo*)
//							info->f:	The file to write to.
//							info->seq:	The sequence that has been masked.
//
// Returns:
//  (nothing)
//
//----------

void print_masking_interval
   (unspos	beg,
	unspos	end,
	void*	info)
	{
	FILE*	f    = ((pmiInfo*) info)->f;
	seq*	_seq = ((pmiInfo*) info)->seq;

	beg += _seq->startLoc - 1;
	end += _seq->startLoc - 1;

	fprintf (f, unsposFmt " " unsposFmt "\n", beg, end);
	}

void print_masking_interval_3
   (unspos	beg,
	unspos	end,
	void*	info)
	{
	FILE*			f    = ((pmiInfo*) info)->f;
	seq*			_seq = ((pmiInfo*) info)->seq;
	seqpartition*	sp   = &_seq->partition;
	partition*		part;
	char*			name;
	unspos			offset;

	// figure out position offsets and names

	if (sp->p == NULL)		// sequence 1 is not partitioned
		{
		name = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
		if ((name == NULL) || (name[0] == 0)) name = "seq1";
		offset = 0;
		}
	else					// sequence 1 is partitioned
	 	{
		part   = lookup_partition (_seq, beg-1);
		name   = &sp->pool[part->header];
		offset = part->sepBefore + 1;
		}

	beg += _seq->startLoc - offset - 1;
	end += _seq->startLoc - offset - 1;

	// print the interval

	fprintf (f, "%s " unsposFmt " " unsposFmt "\n", name, beg, end);
	}

//----------
//
// print_census--
//	Print all locations in a census that meet or exceed the threshold.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	_seq:		The sequence being counted.  If this is NULL, census
//						.. positions are not named.  If it is not NULL, the
//						.. sequence name is shown with each census position.
//	census*	cen:		The census to print.
//	char	delimiter:	The character to use between fields (e.g. tab or space).
//
// Returns:
//	(nothing)
//
//----------

void print_census
   (FILE*			f,
	seq*			_seq,
	census*			cen,
	char			delimiter)
	{
	seqpartition*	sp;
	partition*		nextPart;
	u32				nextIx;
	char*			name;
	unspos			offset;
	unspos			pos;
	u32				count;

	if (cen == NULL) return;

	// simple print with no sequence names

	if (_seq == NULL)
		{
		for (pos=0 ; pos<cen->len ; pos++)
			{
			count = (cen->kind == 'B')? cen->count8 [pos]
			      : (cen->kind == 'W')? cen->count16[pos]
			                          : cen->count32[pos];
			if (count >= cen->maskThresh)
				fprintf(f, unsposFmt "%c%u\n", pos+1, delimiter, count);
			}
		return;
		}

	// print with same sequence name, for non-partitioned sequence

	sp = &_seq->partition;
	if (sp->p == NULL)
		{
		name = (_seq->useFullNames)? _seq->header : _seq->shortHeader;
		if ((name == NULL) || (name[0] == 0)) name = "seq1";
		for (pos=0 ; pos<cen->len ; pos++)
			{
			count = (cen->kind == 'B')? cen->count8 [pos]
			      : (cen->kind == 'W')? cen->count16[pos]
			                          : cen->count32[pos];
			if (count >= cen->maskThresh)
				fprintf(f, "%s%c" unsposFmt "%c%u\n",
				           name, delimiter, pos+1, delimiter, count);
			}
		return;
		}

	// print with sequence names, for partitioned sequence

	nextPart = sp->p;
	nextIx   = 0;
	name     = NULL;
	offset   = 0;

	for (pos=0 ; pos<cen->len ; pos++)
		{
		if (pos == nextPart->sepBefore)
			{
			if (nextIx < sp->len)
				{
				name    = &sp->pool[nextPart->header];
				offset  = nextPart->sepBefore + 1;
				nextPart++;  nextIx++;
				}
			else
				name = NULL;
			}
		else if (name != NULL)
			{
			count = (cen->kind == 'B')? cen->count8 [pos]
			      : (cen->kind == 'W')? cen->count16[pos]
			                          : cen->count32[pos];
			if (count >= cen->maskThresh)
				fprintf(f, "%s%c" unsposFmt "%c%u\n",
				           name, delimiter, pos+1-offset, delimiter, count);
			}
		else
			suicidef ("internal error in print_census\n");
		}

	}

//----------
//
// masking_zero_stats--
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

void masking_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&maskingStats, 0, sizeof(maskingStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// masking_show_stats--
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

void masking_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "      masked bases: %s\n", commatize(maskingStats.maskedBases));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void masking_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "masked_bases=%d\n", maskingStats.maskedBases);
#endif // collect_stats
	}

