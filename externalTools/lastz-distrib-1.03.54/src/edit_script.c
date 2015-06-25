//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: edit_script.c
//
//----------
//
// edit_script--
//	Support for representing alignments as a series of substitute, insert,
//	and delete.
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

#define  edit_script_owner		// (make this the owner of its globals)
#include "edit_script.h"		// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static void edit_script_make_room (editscript** s, u32 entries);
static void edit_script_put       (editscript** s, u32 op, u32 rpt);

//----------
//
// free_align_list--
//	Dispose of a list of alignments.
//
//----------
//
// Arguments:
//	alignel*	a:	The list of alignments to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_align_list
   (alignel*	a)
	{
	alignel*	b;

	while (a != NULL)
		{
		b = a->next;
		free_if_valid ("free_align_list a->script", a->script);
		free_if_valid ("free_align_list a",         a);
		a = b;
		}
	}

//----------
//
// edit_script_new--
//	Allocate an alignment edit script.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A pointer to an empty edit script.  The caller is responsible for disposing
//	of this memory, for which purpose free() can be used.
//
//----------
//
// Notes:  NULL is never returned-- failure to allocate is a fatal error.
//
//----------

editscript* edit_script_new
   (void)
	{
	u32			entries = 12;
	editscript*	s;

	// allocate

	s = zalloc_or_die ("edit_script_new", edit_script_bytes(entries));

	// initialize;  note that by use of zalloc we already have
	//	s->len    = 0;
	//	s->tailOp = 0;

	s->size = entries;

	return s;
	}

//----------
//
// edit_script_make_room--
//	Make sure an alignment edit script has enough unused entries available, and
//	enlarge if it doesn't.
//
//----------
//
// Arguments:
//	editscript**	s:			(pointer to) The script to check/enlarge.  If
//								.. reallocation is required we may alter this.
//	u32				entries:	The number of unused entries required.
//	(none)
//
// Returns:
//	Nothing.  Note that a reallocation failure is a fatal error.
//
//----------

static void edit_script_make_room
   (editscript**	_s,
	u32				entries)
	{
	editscript*		s = *_s;

	// do we have enough space already?

	entries += s->len;
	if (s->size >= entries)
		return;  // (yes)

	// reallocate

	entries += entries/2;	// (anticipate 50% future growth)

	*_s = s = realloc_or_die ("edit_script_has_room",
	                          s, edit_script_bytes(entries));

	s->size = entries;
	}

//----------
//
// edit_script_copy--
//	Allocate a copy of an alignment edit script.
//
//----------
//
// Arguments:
//	editscript*	s:	The script to copy.
//
// Returns:
//	A pointer to a copy of the edit script.  The caller is responsible for
//	disposing of this memory, for which purpose free() can be used.
//
//----------
//
// Notes:  NULL is never returned-- failure to allocate is a fatal error.
//
//----------

editscript* edit_script_copy
   (editscript*	s)
	{
	u32			entries;
	size_t		bytesNeeded;
	editscript*	newS;

	// allocate

	entries     = s->len;
	bytesNeeded = edit_script_bytes(entries);
	newS = zalloc_or_die ("edit_script_copy", bytesNeeded);

	// initialize

	memcpy (/*to*/ newS, /*from*/ s, /*how much*/ bytesNeeded);
	newS->size = entries;

	return newS;
	}

//----------
//
// edit_script_add--
//	Add an repeated operation to an alignment edit script, merging it with the
//	tail of the script if possible.
//
//----------
//
// Arguments:
//	editscript** s:		The script to add to.  If the script has to be
//						.. enlarged, this value may change upon return.
//	u32			 op:	The operation to add.
//	unspos		 rpt:	The repeat count (i.e. how many copies of op to add).
//
// Returns:
//	(nothing).
//
//----------

void edit_script_add
   (editscript**	_s,
	u32				op,
	unspos			rpt)
	{
	editscript*		s = *_s;
	editop*			tail;
	u32				tailRpt;

	// if this operation matches the one currently in the tail, increase the
	// repeat count on the tail;  if the repeat count doesn't have enough count
	// left, fall through to the loop

	if (edit_op_operation(s->tailOp) == op)
		{
		tail    = s->op + s->len-1;
		tailRpt = edit_op_repeat (*tail);

		if (tailRpt + rpt <= maxEditopRepeat)
			{
			*tail = edit_op_add_repeat (*tail, rpt);
			return;
			}
		else
			{
			*tail = edit_op (op, maxEditopRepeat);
			rpt = tailRpt + rpt - maxEditopRepeat;
			}
		}

	// loop, adding new operation(s) to the end of the script

	while (rpt > maxEditopRepeat)
		{
		edit_script_put (_s, op, maxEditopRepeat);
		rpt -= maxEditopRepeat;
		}

	edit_script_put (_s, op, rpt);
	}

//----------
//
// edit_script_put--
//	Add an repeated operation to an alignment edit script, tacking a new entry
//	onto the script.
//
//----------
//
// Arguments:
//	editscript** s:			The script to add to.  If the script has to be
//							.. enlarged, this value may change upon return.
//	u32			 op, rpt:	The operation to add, including a repeat count.
//
// Returns:
//	(nothing).
//
//----------

static void edit_script_put
   (editscript**	_s,
	u32				op,
	u32				rpt)
	{
	editscript*		s;

	edit_script_make_room (_s, 1); // (make sure we have room for one operation)

	s = *_s;
	s->op[s->len++] = edit_op (op, rpt);
	s->tailOp = op;
	}

//----------
//
// edit_script_append--
//	Copy one alignment edit script to the end of another.
//
//----------
//
// Arguments:
//	editscript**	dst:	(pointer to) The script to copy to.  If
//							.. reallocation is required we may alter this.
//	editscript*		src:	(pointer to) The script to copy from.
//
// Returns:
//	A pointer to an empty edit script.  The caller is responsible for disposing
//	of this memory, for which purpose free() can be used.
//
//----------

void edit_script_append
   (editscript**	_dst,
	editscript*		src)
	{
	editscript*		dst;
	editop*			s, *d;
	u32				toCopy;
	u32				sOp;
	u32				sRpt, dRpt;

	if (src->len == 0) return;

	// make sure we have enough room

	edit_script_make_room (_dst, src->len);
	dst = *_dst;

	// copy dst to src

	s  = src->op;
	d  = dst->op + dst->len-1;
	toCopy = src->len;

	sOp = edit_op_operation (*s);
	if (sOp == dst->tailOp)
		{
		dRpt = edit_op_repeat (*d);
		sRpt = edit_op_repeat (*s);
		if (dRpt + sRpt <= maxEditopRepeat)
			*d = edit_op_add_repeat (*d, sRpt);
		else
			{
			*(d++) = edit_op (sOp, maxEditopRepeat);
			*d     = edit_op (sOp, dRpt + sRpt - maxEditopRepeat);
			dst->len++;
			}
		s++;  toCopy--;
		}
	d++;

	memcpy (d, s, toCopy*sizeof(editop));

	dst->len    += toCopy;
	dst->tailOp =  src->tailOp;
	}

//----------
//
// edit_script_reverse--
//	Reverse the items in an alignment edit script, in place.
//
//----------
//
// Arguments:
//	editscript*	s:	The script to modify.
//
// Returns:
//	(nothing)
//
//----------

void edit_script_reverse
   (editscript*	s)
	{
	u32			i, j;
	editop		t;

	if (s->len < 2) return;

	for (i=0,j=s->len-1 ; i<j ; i++,j--)
		{ t = s->op[i];  s->op[i] = s->op[j];  s->op[j] = t; }
	}

//----------
//
// edit_script_mirror--
//	Flip an alignment edit script across the main diagonal, in place.
//
// This amounts to changing deletions to insertions, and vice-versa.
//
//----------
//
// Arguments:
//	editscript*	s:	The script to modify.
//
// Returns:
//	(nothing)
//
//----------

void edit_script_mirror
   (editscript*	s)
	{
	u32			i;
	editop		op;
	u32			rpt;

	for (i=0 ; i<s->len ; i++)
		{
		op = s->op[i];
		switch (edit_op_operation(op))
			{
			case editopIns:
				rpt = edit_op_repeat(op);
				s->op[i] = edit_op (editopDel,rpt);
				break;
			case editopDel:
				rpt = edit_op_repeat(op);
				s->op[i] = edit_op (editopIns,rpt);
				break;
			default:
				// do nothing
				break;
			}
		}
	}

//----------
//
// edit_script_trim_head--
//	Trim some number of steps off the head of an alignment edit script, in
//	place.
//
//----------
//
// Arguments:
//	editscript*	s:		The script to modify.
//	unspos		len:	The number of steps to trim.  A step is a one-base
//						.. step, in any direction.
//
// Returns:
//	(nothing)
//
//----------

void edit_script_trim_head
   (editscript*	s,
	unspos		len)
	{
	u32			i, j;
	editop		op;
	u32			rpt;
	int			shortScript;

	if (s->len == 0) return; // the alignment is empty
	if (len == 0)    return; // nothing to trim

	// scan to find the first segment that we won't completely skip

	shortScript = true;
	for (i=0 ; i<s->len ; i++)
		{
		op = s->op[i];
		rpt = edit_op_repeat(op);
		if (rpt > len)
			{ shortScript = false;  break; }
		len -= rpt;
		}

	if (shortScript) // the alignment didn't have enough steps
		{ s->len = 0;  return; }

	// if we skipped whole segments, shift the remaining segments forward

	if (i > 0)
		{
		for (j=i ; j<s->len ; j++)
			s->op[j-i] = s->op[j];
		s->len -= i;
		}

	// if we have anything else to trim (other than those whole segments),
	// trim the first segment

	if (len > 0)
		{
		op = edit_op_operation(op);
		s->op[0] = edit_op (op,rpt-len);
		}
	}

//----------
//
// edit_script_upper_truncate--
//	Truncate an alignment edit script at the main diagonal, as it crosses from
//	upper triangle to lower triangle.
//
// We assume that the first sequence is along the positive strand, and the
// second is along the negative strand.
//
//----------
//
// Arguments:
//	editscript*	s:			The script to modify.
//	unspos*		pos1, pos2:	(pointer to) The starting position of the script,
//							.. in the DP matrix, both relative to the
//							.. corresponding positive strand.  If we return
//							.. true, the position of the truncated ending is
//							.. written to these values.  If the truncated
//							.. alignment is empty, seqposInfinity is written to
//							.. both values.
//
// Returns:
//	true if truncation has occured, in which case the new ending values for
//	pos1 and pos2 have been written;  false otherwise.
//
//----------
//
// notes:
//	(1) Positions with respect to the diagonal are as follows:
//	      pos1 <  pos2    =>  above the diagonal
//	      pos1 == pos2    =>  on the diagonal
//	      pos1 >  pos2    =>  below the diagonal
//	    However, we will truncate at a point with either of these cases:
//	      pos1 == pos2    =>  on the diagonal
//	      pos1 == pos2+1  =>  adjacent to the diagonal (below it)
//	    The latter case is necessary because an alignment can cross the diagonal
//	    without, technically, containing any positions on the diagonal.
//
//----------

int edit_script_upper_truncate
   (editscript*	s,
	unspos*		_pos1,
	unspos*		_pos2)
	{
	u32			i;
	editop		op;
	u32			rpt;
	unspos		pos1, pos2, prevPos1, prevPos2;
	int			reachesDiagonal;
	unspos		limit = 0; // (to placate compiler)

	//fprintf (stderr, "s=%08X l=%08X\n", s->size, s->len);

	// handle special cases

	if (s->len == 0)
		{
		// the alignment is empty
		return false;
		}

	pos1 = (*_pos1);
	pos2 = (*_pos2);

	//fprintf (stderr, unsposSlashFmt "\n", pos1, pos2);

	if (pos1 > pos2)
		{
		// the alignment starts below the diagonal, discard all of it
		s->len = 0;
		(*_pos1) = seqposInfinity;
		(*_pos2) = seqposInfinity;
		return true;
		}

	// scan for the first segment that touches or crosses the main diagonal

	reachesDiagonal = false;
	for (i=0 ; i<s->len ; i++)
		{
		prevPos1 = pos1;
		prevPos2 = pos2;

		op  = s->op[i];
		rpt = edit_op_repeat(op);
		op  = edit_op_operation(op);
		switch (op)
			{
			case editopSub: pos1 += rpt;  pos2 -= rpt;  limit = pos2+1; break;
			case editopIns:               pos2 -= rpt;  limit = pos2;   break;
			case editopDel: pos1 += rpt;                limit = pos2;   break;
			}

		if (pos1 >= limit)
			{ reachesDiagonal = true;  break; }
		}

	if (!reachesDiagonal) return false;

	//fprintf (stderr, "diagonal reached\n");
	//fprintf (stderr, unsposCommaFmt " -> " unsposCommaFmt "\n",
	//                 prevPos1, prevPos2, pos1, pos2);

	// truncate the list at the crossing segment

	s->len = i+1;

	// split the crossing segment (unless we're lucky enough that it ended
	// at the diagonal)

	if (pos1 > pos2)
		{
		switch (op)
			{
			case editopSub:
				// prevPos    pos      new pos    rpt
				// (90,110)  (105,95)  (100,100)  10
				// (90,111)  (104,95)  (101,100)  11
				rpt = (prevPos2+1 - prevPos1) / 2;
				s->op[i] = edit_op (editopSub,rpt);
				pos1 = prevPos1 + rpt;
				pos2 = prevPos2 - rpt;
				break;
			case editopIns:
				// prevPos    pos       new pos    rpt
				// (100,110)  (100,95)  (100,100)  10
				rpt = prevPos2 - prevPos1;
				s->op[i] = edit_op (editopIns,rpt);
				pos1 = prevPos1;
				pos2 = prevPos2 - rpt;
				break;
			case editopDel:
				// prevPos    pos       new pos    rpt
				// (90,100)  (105,100)  (100,100)  10
				rpt = prevPos2 - prevPos1;
				s->op[i] = edit_op (editopDel,rpt);
				pos1 = prevPos1 + rpt;
				pos2 = prevPos2;
				break;
			}

		//fprintf (stderr, "truncated to %u -> " unsposCommaFmt "\n",
		//                 rpt, pos1, pos2);
		}

	(*_pos1) = pos1;
	(*_pos2) = pos2;
	return true;
	}

//----------
//
// edit_script_run_of_subs, edit_script_run_of_subs_match--
//	Find the length of the current run of substitutions in an alignment edit
//	script.
//
//----------
//
// Arguments:
//	editscript*	s:		The script being parsed.
//	u32*		opIx:	Current parse location in the script.  Upon return this
//						.. is updated to point to the next operation beyond the
//						.. run.
//	const u8*	p, q:	Pointer to the sequences' nucleotides (corresponding to
//						.. the parse location).  These are used only to provide
//						.. a match count, and can be NULL if match is NULL.
//	unspos*		match:	Place to return the number of nucleotide matches in the
//						.. run (including upper/lower mismatches).
//
// Returns:
//	The length of the run.  Note that this could be zero.
//
//----------

u32 edit_script_run_of_subs
   (editscript*	s,
	u32*		_opIx)
	{
	u32			opIx = (u32) *_opIx;
	u32			rpt, run;

	run = 0;
	while ((opIx < s->len) && (edit_op_operation(s->op[opIx]) == editopSub))
		{
		rpt = edit_op_repeat(s->op[opIx]);  opIx++;
		run += rpt;
		}

	*_opIx = opIx;
	return run;
	}


u32 edit_script_run_of_subs_match
   (editscript*	s,
	u32*		_opIx,
	const u8*	p,
	const u8*	q,
	unspos*		_match)
	{
	u32			opIx  = *_opIx;
	unspos		match = *_match;
	u32			rpt, run;
	u8			pCh, qCh;

	run = 0;
	match = 0;
	while ((opIx < s->len)
	    && (edit_op_operation(s->op[opIx]) == editopSub))
		{
		rpt = edit_op_repeat(s->op[opIx]);  opIx++;
		run += rpt;
		while (rpt-- > 0)
			{
            pCh = *(p++);  qCh = *(q++);
            if (dna_toupper(pCh) == dna_toupper(qCh)) match++;
			}
		}

	*_opIx  = opIx;
	*_match = match;
	return run;
	}

//----------
//
// edit_script_indel_len--
//	Find the length of the current "run" of indels in an alignment edit script.
//
//----------
//
// Arguments:
//	editscript*	s:		The script being parsed.
//	u32*		opIx:	Current parse location in the script.  Upon return this
//						.. is updated to point to the next operation beyond the
//						.. indel.
//	unspos*		i, j:	Current parse location in the sequences.  Upon return
//						.. these are updated.
//
// Returns:
//	The length of the run.
//
//----------

u32 edit_script_indel_len
   (editscript*	s,
	u32*		opIx,
	unspos*		i,
	unspos*		j)
	{
	editop		op;
	u32			rpt;

	if (s->len <= (u32) *opIx)
		return 0;

	op = s->op[*opIx];
	rpt = edit_op_repeat(op);

	switch (edit_op_operation(op))
		{
		case editopIns: *j += rpt;  break;
		case editopDel: *i += rpt;  break;
		}

	(*opIx)++;
	return rpt;
	}

//----------
//
// edit_script_overall_len--
//	Find the length of an alignment edit script, along both sequences.
//
//----------
//
// Arguments:
//	editscript*	s:		The script to modify.
//	unspos*		i, j:	Place to return the lengths.  i is along the first
//						.. sequence;  j is along the second.
//
// Returns:
//	(nothing)
//
//----------

void edit_script_overall_len
   (editscript*	s,
	unspos*		_i,
	unspos*		_j)
	{
	u32			opIx;
	editop		op;
	u32			rpt;
	unspos		i, j;

	i = j = 0;

	for (opIx=0 ; opIx<s->len ; opIx++)
		{
		op  = s->op[opIx];
		rpt = edit_op_repeat(op);
		switch (edit_op_operation(op))
			{
			case editopSub: i += rpt;  j += rpt;  break;
			case editopIns:            j += rpt;  break;
			case editopDel: i += rpt;             break;
			}
		}

	(*_i) = i;
	(*_j) = j;
	}

//----------
//
// dump_edit_script--
//	Print the raw contents of an alignment edit script.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print to.
//	editscript* s:	The script to print.
//
// Returns:
//	(nothing).
//
//----------

void dump_edit_script
   (FILE*		f,
	editscript*	s)
	{
	char*		opName[4] = { "???", "INS", "DEL", "SUB" };
	u32			op, rpt;
	u32			ix;

	for (ix=0 ; ix<s->len ; ix++)
		{
		op  = edit_op_operation (s->op[ix]);
		rpt = edit_op_repeat    (s->op[ix]);
		fprintf (f, "%dx%s\n", rpt, opName[op]);
		}
	}

