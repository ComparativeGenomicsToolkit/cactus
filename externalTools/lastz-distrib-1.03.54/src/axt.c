//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: axt.c
//
//----------
//
// axt--
//	Support for printing alignments in AXT format.
//
// AXT format is a well-established pairwise alignment format.  As of Jan/2009,
// a spec for AXT files can be found at
//	http://genome.ucsc.edu/goldenPath/help/axt.html
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
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "identity_dist.h"		// identity distribution "format" stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "output.h"				// alignment outout format stuff
#include "genpaf.h"				// genpaf alignment format stuff

#define  axt_owner				// (make this the owner of its globals)
#include "axt.h"				// interface to this module

// alignment counter

static u64 axtAlignmentNumber;

//----------
//
// print_axt_job_header--
//	Print axt format job header.
//
//----------

void print_axt_job_header
   (FILE*		f,
	char*		_programName,
	char*		_args,
	scoreset*	scoring,
	sthresh*	hspThreshold,
	sthresh*	gappedThreshold,
	score		xDrop,
	score		yDrop)
	{
	char*		programName = _programName;
	char*		args        = _args;

	if (programName == NULL) programName = "(no name)";
	if (args        == NULL) args        = "";

	fprintf (f, "# %s %s\n", programName, args);
	fprintf (f, "#\n");
	fprintf (f, "# hsp_threshold      = %s\n", score_thresh_to_string (hspThreshold));
	fprintf (f, "# gapped_threshold   = %s\n", score_thresh_to_string (gappedThreshold));
	fprintf (f, "# x_drop             = " scoreFmtSimple "\n", xDrop);
	fprintf (f, "# y_drop             = " scoreFmtSimple "\n", yDrop);
	print_score_matrix_prefix (f, scoring, true, "# ");

	axtAlignmentNumber = (u64) -1;     // caveat:  this only works properly if
	                                   // .. we only write one axt file at a
	                                   // .. time, and write it completely
	}

//----------
//
// print_axt_job_footer--
//	Print axt format job footer.
//
//----------

void print_axt_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_axt_header--
//	Print axt format query header.
//
//----------

void print_axt_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2))
	{
	// (do nothing)
	}

//----------
//
// print_axt_align_list--
//	Print a list of gapped alignments in axt format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	alignel*	alignList:		The list of alignments to print.
//	seq*		seq1:			One sequence.
//	seq*		seq2:			Another sequence.
//	int			withComments:	true => print comments as well
//	char*		extras:			Extra fields to print on the summary line.
//								.. These are genpaf keys, defined in genpaf.h
//								.. (genpafXXX values).  Currently only
//								.. genpafSize2 is supported.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_axt_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			withComments,
	char*		extras)
	{
	alignel*	a;
	unspos		numer, denom;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		if (withComments)
			{
			// report identity
			alignment_identity (seq1, seq2, a, &numer, &denom);
			fprintf (f, "# identity=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");

			// report coverage
			alignment_coverage (seq1, seq2, a, &numer, &denom);
			fprintf (f, "# coverage=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");
			}

		print_axt_align (f,
		                 seq1, a->beg1-1, a->end1,
		                 seq2, a->beg2-1, a->end2,
		                 a->script, a->s, extras);
		}

	}

//----------
//
// print_axt_align--
//	Print a single gapped alignment in axt format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	seq*		seq1:			One sequence.
//	unspos		beg1, end1:		Range of positions in sequence 1 (origin 0).
//	seq*		seq2:			Another sequence.
//	unspos		beg2, end2:		Range of positions in sequence 2 (origin 0).
//	editscript*	script:			The script describing the path the alignment
//								.. takes in the DP matrix.
//	score		s:				The alignment's score.
//	char*		extras:			Extra fields to print on the summary line.
//								.. These are genpaf keys, defined in genpaf.h
//								.. (genpafXXX values).  Currently only
//								.. genpafSize2 is supported.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_axt_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	score			s,
	char*			extras)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, startI, startJ, run;
	u32				opIx;
	u8*				p, *q;
	unspos			ix, len1, len2;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq2Len, seq2True;
	char			strand2;

	if ((extras != NULL) && (strlen(extras) != 1) && (extras[0] != genpafSize2))
		suicide ("internal error: print_axt_align doesn't support extras");

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_axt_align");

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	len1 = height = end1 - beg1 + 1;
	len2 = width  = end2 - beg2 + 1;

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1-1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2Len   = seq2->len;
		seq2True  = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, beg2-1);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		}

	//////////
	// print summary line
	//////////

	axtAlignmentNumber++;

	start1 = beg1-1 - offset1 + startLoc1;

	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = beg2-1 - offset2 + startLoc2;
		strand2 = '+';
		}
	else
		{
		start2  = beg2-1 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		strand2 = '-';
		}

	fprintf (f, u64Fmt " %s " unsposFmt " " unsposFmt
	            " %s "    unsposFmt " " unsposFmt " %c " scoreFmt,
	            axtAlignmentNumber,
	            name1, start1, start1+len1-1,
	            name2, start2, start2+len2-1, strand2,
	            s);
	if ((extras != NULL) && (strlen(extras) == 1) && (extras[0] == genpafSize2))
		fprintf (f, " " unsposFmt, seq2Len);
	fprintf (f, "\n");

	//////////
	// print aligning path in sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)
	//////////

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i-1;
		q = seq2->v+beg2+j-1;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*p));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i-1;
			startJ = j;  q = seq2->v+beg2+j-1;

			edit_script_indel_len (script, &opIx, &i, &j);

			if (i != startI)
				{
				for ( ; startI<i ; startI++)
					{ fprintf (f, "%c", dna_toprint(*p));  p++; }
				}

			if (j != startJ)
				{
				for ( ; startJ<j ; startJ++)
					{ fprintf (f, "-");  q++; }
				}
			}
		}

	fprintf (f, "\n");

	//////////
	// print aligning path in sequence 2
	//////////

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i-1;
		q = seq2->v+beg2+j-1;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*q));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i-1;
			startJ = j;  q = seq2->v+beg2+j-1;

			edit_script_indel_len (script, &opIx, &i, &j);

			if (i != startI)
				{
				for ( ; startI<i ; startI++)
					{ fprintf (f, "-");  p++; }
				}

			if (j != startJ)
				{
				for ( ; startJ<j ; startJ++)
					{ fprintf (f, "%c", dna_toprint(*q));  q++; }
				}
			}
		}

	fprintf (f, "\n\n");
	}

//----------
//
// print_axt_match--
//	Print an hsp in axt format.
//
//----------
//
// Arguments:
//	FILE*	f:				The file to print to.
//	seq*	seq1:			One sequence.
//	unspos	pos1:			The position, in seq1, of first character in the
//							.. match (origin-0).
//	seq*	seq2:			Another sequence.
//	unspos	pos1:			The position, in seq2, of first character in the
//							.. match (origin-0).
//	unspos	length:			The number of nucleotides in the HSP.
//	score	s:				The HSP's score.
//	int		withComments:	true => print comments as well
//	char*	extras:			Extra fields to print on the summary line.  These
//							.. are genpaf keys, defined in genpaf.h (genpafXXX
//							.. values).  Currently only genpafSize2 is
//							.. supported.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_axt_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	int				withComments,
	char*			extras)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, end1, start2, end2;
	unspos			startLoc1, startLoc2;
	unspos			seq2Len, seq2True;
	char			strand2;
	unspos			ix;
	segment			seg;
	unspos			numer, denom;

	if ((extras != NULL) && (strlen(extras) != 1) && (extras[0] != genpafSize2))
		suicide ("internal error: print_axt_match doesn't support extras");

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_axt_match");

	if (withComments)
		{
		// report identity
		segment_identity (seq1, pos1, seq2, pos2, length, &numer, &denom);
		fprintf (f, "# identity=" unsposSlashFmt, numer, denom);
		if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
		fprintf (f, "\n");

		// report coverage
		seg.pos1   = pos1;
		seg.pos2   = pos2;
		seg.length = length;
		segment_coverage (seq1, seq2, &seg, &numer, &denom);
		fprintf (f, "# coverage=" unsposSlashFmt, numer, denom);
		if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
		fprintf (f, "\n");
		}

	// figure out position offsets and names

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, pos1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2Len   = seq2->len;
		seq2True  = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, pos2);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		}

	// print summary line

	axtAlignmentNumber++;

	start1 = pos1 - offset1 + startLoc1;
	end1   = start1 + length;

	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = pos2 - offset2 + startLoc2;
		end2    = start2 + length;
		strand2 = '+';
		}
	else
		{
		start2  = pos2 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		end2    = start2 + length;
		strand2 = '-';
		}

	fprintf (f, u64Fmt " %s " unsposFmt " " unsposFmt
	            " %s "    unsposFmt " " unsposFmt " %c " scoreFmt,
	            axtAlignmentNumber,
	            name1, start1, end1-1,
	            name2, start2, end2-1, strand2,
	            s);
	if ((extras != NULL) && (strlen(extras) == 1) && (extras[0] == genpafSize2))
		fprintf (f, " " unsposFmt, seq2Len);
	fprintf (f, "\n");

	// print aligning segment of sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s1[ix]));
	fprintf (f, "\n");

	// print aligning segment of sequence 2

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s2[ix]));
	fprintf (f, "\n\n");
	}

//----------
//
// print_axt_comment--
//	Print a comment in a axt file.
//
// Note that the axt format does *not* support comments, so these records will
// invariably spoil the axt file for most downstream tools.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to print to.
//	const char*	format:	A format string, as per printf.
//	...:				(same as for printf)
//
// Returns:
//	(nothing)
//
//----------

void print_axt_comment
   (FILE*		f,
	const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);
	vprint_axt_comment (f, format, args);
	va_end (args);
	}

void vprint_axt_comment
   (FILE*		f,
	const char*	format,
	va_list		args)
	{
	fprintf (f, "# ");
	if (format != NULL)
		vfprintf (f, format, args);
	fprintf  (f, "\n");
	}

