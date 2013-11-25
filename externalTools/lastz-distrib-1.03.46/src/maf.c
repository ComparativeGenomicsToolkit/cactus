//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: maf.c
//
//----------
//
// maf--
//	Support for printing alignments in MAF format.
//
// MAF format is a well-established multiple alignment format.  As of Jan/2009,
// a spec for MAF files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format5
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
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "identity_dist.h"		// identity distribution stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "continuity_dist.h"	// query continuity distribution stuff

#define  maf_owner				// (make this the owner of its globals)
#include "maf.h"				// interface to this module

static int max_digits (s64 num1, s64 num2);

// debugging defines

//#define debugSeq1Beg 858		// if defined, only alignments entirely within
//#define debugSeq1End 1153		// .. this range in sequence 1 are output;  note
								// .. that these positions are origin-zero

//#define snoopMafAlignList		// if this is defined, extra code is added to
								// .. track alignment lists

//----------
//
// print_maf_job_header--
//	Print maf format job header.
//
//----------

void print_maf_job_header
   (FILE*		f,
	char*		_programName,
	char*		_args,
	scoreset*	scoring,
	sthresh*	hspThreshold,
	sthresh*	gappedThreshold,
	score		xDrop,
	score		yDrop,
	int			withComments)
	{
	char*		programName = _programName;
	char*		args        = _args;

	if (!withComments) return;

	if (programName == NULL) programName = "(no name)";
	if (args        == NULL) args        = "";

	fprintf (f, "##maf version=1 scoring=%s\n", programName);
	fprintf (f, "# %s %s\n", programName, args);
	fprintf (f, "#\n");
	fprintf (f, "# hsp_threshold      = %s\n", score_thresh_to_string (hspThreshold));
	if ((gappedThreshold->t == 'S') || (hspThreshold->t == 'S'))
		fprintf (f, "# gapped_threshold   = %s\n", score_thresh_to_string (gappedThreshold));
	else
		fprintf (f, "# gapped_threshold   = (derived from hsp_threshold)\n");
	fprintf (f, "# x_drop             = " scoreFmtSimple "\n", xDrop);
	fprintf (f, "# y_drop             = " scoreFmtSimple "\n", yDrop);
	print_score_matrix_prefix (f, scoring, true, "# ");
	}

//----------
//
// print_maf_job_footer--
//	Print maf format job footer.
//
//----------

void print_maf_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_maf_header--
//	Print maf format query header.
//
//----------

void print_maf_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2))
	{
	// (do nothing)
	}

//----------
//
// print_maf_align_list--
//	Print a list of gapped alignments in maf format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	alignel*	alignList:		The list of alignments to print.
//	seq*		seq1:			One sequence.
//	seq*		seq2:			Another sequence.
//	int			withComments:	true => print comments as well
//
// Returns:
//	(nothing)
//
//----------

//=== stuff for snoopMafAlignList ===

#ifndef snoopMafAlignList
#define snoopMafAlignList_1 ;
#endif // not snoopMafAlignList

#ifdef snoopMafAlignList

#define snoopMafAlignList_1                                                    \
	fprintf (stderr, "print_maf_align_list  a=%08lX"                           \
	                 "  a->seq1=%08lX  a->seq2=%08lX\n",                       \
	                 (long) a, (long) a->seq1, (long) a->seq2);

#endif // snoopMafAlignList


// print_maf_align_list--

void print_maf_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			withComments)
	{
	alignel*	a;
	unspos		numer, denom;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		snoopMafAlignList_1;
		if (withComments)
			{
			unspos height, width, i, j, prevI, prevJ, run;
			u32    opIx;

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

			// report continuity
			alignment_continuity (a, &numer, &denom);
			fprintf (f, "# continuity=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");

			// report alignment path

			fprintf (f, "# cigar=");

			height = a->end1 - a->beg1 + 1;
			width  = a->end2 - a->beg2 + 1;

			opIx = 0;
			for (i=j=0 ; (i< height)||(j<width) ; )
				{
				run = edit_script_run_of_subs (a->script, &opIx);
				if (run > 0)
					{
					fprintf (f, unsposFmt "m", run);
					i += run; j += run;
					}
		
				if ((i < height) || (j < width))
					{
					prevI = i;  prevJ = j;
					edit_script_indel_len (a->script, &opIx, &i, &j);
					if (i > prevI)
						fprintf (f, unsposFmt "d", i - prevI);
					if (j > prevJ)
						fprintf (f, unsposFmt "i", j - prevJ);
					}
				}
			fprintf (f, "\n");
			}

		print_maf_align (f,
		                 seq1, a->beg1-1, a->end1,
		                 seq2, a->beg2-1, a->end2,
		                 a->script, a->s);
		}
	}

//----------
//
// print_maf_align--
//	Print a single gapped alignment in maf format.
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
//
// Returns:
//	(nothing)
//
//----------

static char* rcfSuffix[4] = { "", "~", "~", "" };

void print_maf_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	score			s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, run;
	u32				opIx;
	u8*				p, *q;
	unspos			ix;
	char*			name1, *name2, *pref2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			startI, startJ;
	int				len1, len2, nameW, startW, endW, lenW;

#ifdef debugSeq1Beg
	if ((beg1 < debugSeq1Beg) || (end1 > debugSeq1End)) return;
#endif // debugSeq1Beg

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	// report diagonal

	if (maf_dbgReportDiag)
		fprintf (f, "# diagonal=" sgnposFmt "\n", diagNumber(beg1,beg2));

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		seq1Len   = seq1->len;
		seq1True  = seq1->trueLen;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1-1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		seq1Len   = part->sepAfter - offset1;
		seq1True  = part->trueLen;
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
		startLoc2 = part->startLoc;
		offset2   = part->sepBefore + 1;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		}

	//////////
	// print summary line
	//////////

	fprintf (f, "a score=" scoreFmt "\n", s);

	//////////
	// print aligning path in sequence 1
	//////////

	// figure out fields and widths

	pref2 = ((maf_distinguishNames) && (strcmp (name1, name2) == 0))? "~" : "";
	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = beg1-1 - offset1 + startLoc1;
		strand1 = '+';
		}
	else
		{
		start1  = beg1-1 - offset1 + seq1True+2 - (startLoc1 + seq1Len);
		strand1 = '-';
		}
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

	len1  =                  strlen (name1) + strlen (suff1);
	len2  = strlen (pref2) + strlen (name2) + strlen (suff2);
	nameW = (len1 >= len2)? len1 : len2;

	startW = max_digits (start1, start2);
	endW   = max_digits (end1+1-beg1, end2+1-beg2);
	lenW   = max_digits (seq1True, seq2True);

	// print aligning path in sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)

	fprintf (f, "s %s%s%*s" unsposStarFmt " " unsposStarFmt " %c " unsposStarFmt " ",
	            name1, suff1, nameW+1-len1, " ",
	            startW, start1-1, endW, end1+1-beg1, strand1, lenW, seq1True);

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

	fprintf (f, "s %s%s%s%*s" unsposStarFmt " " unsposStarFmt " %c " unsposStarFmt " ",
	            pref2, name2, suff2, nameW+1-len2, " ",
	            startW, start2-1, endW, end2+1-beg2, strand2, lenW, seq2True);

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
// print_maf_match--
//	Print an hsp in maf format.
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
//
// Returns:
//	(nothing)
//
//----------

void print_maf_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	int				withComments)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	char*			name1, *name2, *pref2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	int				len1, len2, nameW, startW, lenW;
	unspos			ix;
	segment			seg;
	unspos			numer, denom;

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_maf_match");

#ifdef debugSeq1Beg
	if ((pos1 < debugSeq1Beg) || (pos1+length > debugSeq1End)) return;
#endif // debugSeq1Beg

	// report diagonal

	if (maf_dbgReportDiag)
		fprintf (f, "# diagonal=" sgnposFmt "\n", diagNumber(pos1,pos2));

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

		// report alignment path
		fprintf (f, "# cigar=" unsposFmt "m\n", length);
		}

	// figure out position offsets and names

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		seq1Len   = seq1->len;
		seq1True  = seq1->trueLen;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, pos1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		seq1Len   = part->sepAfter - offset1;
		seq1True  = part->trueLen;
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

	fprintf (f, "a score=" scoreFmt "\n", s);

	// figure out fields and widths

	pref2 = ((maf_distinguishNames) && (strcmp (name1, name2) == 0))? "~" : "";
	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = pos1 - offset1 + startLoc1;
		strand1 = '+';
		}
	else
		{
		start1  = pos1 - offset1 + seq1True+2 - (startLoc1 + seq1Len);
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = pos2 - offset2 + startLoc2;
		strand2 = '+';
		}
	else
		{
		start2  = pos2 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		strand2 = '-';
		}

	len1  =                  strlen (name1) + strlen (suff1);
	len2  = strlen (pref2) + strlen (name2) + strlen (suff2);
	nameW = (len1 >= len2)? len1 : len2;

	startW = max_digits (start1, start2);
	lenW   = max_digits (seq1True, seq2True);

	// print aligning segment of sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)

	fprintf (f, "s %s%s%*s" unsposStarFmt " " unsposFmt " %c " unsposStarFmt " ",
	            name1, suff1, nameW+1-len1, " ",
	            startW, start1-1, length, strand1, lenW, seq1True);

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s1[ix]));
	fprintf (f, "\n");

	// print aligning segment of sequence 2

	fprintf (f, "s %s%s%s%*s" unsposStarFmt " " unsposFmt " %c " unsposStarFmt " ",
	            pref2, name2, suff2, nameW+1-len2, " ",
	            startW, start2-1, length, strand2, lenW, seq2True);

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s2[ix]));
	fprintf (f, "\n\n");
	}

//----------
//
// print_maf_comment--
//	Print a comment in a maf file.
//
// Note that the maf format does *not* support comments, so these records will
// invariably spoil the maf file for most downstream tools.
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

void print_maf_comment
   (FILE*		f,
	const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);
	vprint_maf_comment (f, format, args);
	va_end (args);
	}

void vprint_maf_comment
   (FILE*		f,
	const char*	format,
	va_list		args)
	{
	fprintf (f, "# ");
	if (format != NULL)
		vfprintf (f, format, args);
	fprintf  (f, "\n");
	}

//----------
//
// max_digits--
//	Figure out the number of digits required to print either of two numbers.
//
//----------
//
// Arguments:
//	s64 num1,num2: The two numbers.
//
// Returns:
//	The number of characters required for printing.
//
//----------

static int max_digits
   (s64		num1,
	s64		num2)
	{
	int w1 = snprintf (NULL, 0, s64Fmt, num1);
	int w2 = snprintf (NULL, 0, s64Fmt, num2);
	return (w1 >= w2)? w1 : w2;
	}
