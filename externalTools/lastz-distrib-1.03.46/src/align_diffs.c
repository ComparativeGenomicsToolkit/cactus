//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: align_diffs.c
//
//----------
//
// align_diffs--
//	Support for printing alignments in a textual differences format.
//
// For an alignment like this:
//
//	s phiX                  4294 35 + 5386 CCCCCAACTTGATATTAATAACACTATAGACCACC
//	s HWI-EAS91_1_306UPAAXX    1 35 -   36 CCCCCATCTTGATATTAATAACACTATAGACCACC
//
// The output will be something like this (but all on one line):
//
//	phiX                  4300 4301 + 5386
//  HWI-EAS91_1_306UPAAXX    7    8 -   36
//  A T
//  CCCCCAACTTGATATTAATAACACTATAGACCACC CCCCCATCTTGATATTAATAACACTATAGACCACC
//
// Note that we use the same position conventions are per MAF format, zero-based,
// half-open, and negative strand positions are counted from the 5' end of the
// *negative* strand.
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
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

#define  align_diffs_owner		// (make this the owner of its globals)
#include "align_diffs.h"		// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static void print_align_difference (FILE* f,
                                    seq* seq1, unspos beg1, unspos end1,
                                    seq* seq2, unspos beg2, unspos end2,
                                    editscript* script,
                                    unspos diffPos1, u8* diffText1,
                                    unspos diffPos2, u8* diffText2,
                                    unspos diffLength,
                                    int withBlocks);

static void print_match_difference (FILE* f,
                                    seq* seq1, unspos pos1, unspos diffPos1,
                                    seq* seq2, unspos pos2, unspos diffPos2,
                                    unspos length, unspos diffLength,
                                    int withBlocks);

//----------
//
// print_align_diffs_job_header--
//	Print a alignment differences job header.
//
//----------

void print_align_diffs_job_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(char* programName),
	arg_dont_complain(char* name1),
	arg_dont_complain(char* name2))
	{
	// (do nothing)
	}

//----------
//
// print_align_diffs_job_footer--
//	Print a alignment differences job footer.
//
//----------

void print_align_diffs_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_align_diffs_header--
//	Print a alignment differences query header.
//
//----------

void print_align_diffs_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2))
	{
	// (do nothing)
	}

//----------
//
// print_align_diffs_align_list--
//	Print a list of gapped alignments, textually.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//	int			withBlocks:	true => include full alignment block in output.
//	int			inhibitN:	true => don't report mismatch-with-N differences.
//
// Returns:
//	(nothing)
//
//----------

void print_align_diffs_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			withBlocks,
	int			inhibitN)
	{
	alignel* a;

	for (a=alignList ; a!=NULL ; a=a->next)
		print_align_diffs_align (f,
		                        seq1, a->beg1-1, a->end1,
		                        seq2, a->beg2-1, a->end2,
		                        a->script,
		                        withBlocks, inhibitN);
	}

//----------
//
// print_align_diffs_align--
//	Print a single gapped alignment, textually.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	seq*		seq1:		One sequence.
//	unspos		beg1, end1:	Range of positions in sequence 1 (origin 0).
//	seq*		seq2:		Another sequence.
//	unspos		beg2, end2:	Range of positions in sequence 2 (origin 0).
//	editscript*	script:		The script describing the path the alignment takes
//							.. in the DP matrix.
//	int			withBlocks:	true => include full alignment block in output.
//	int			inhibitN:	true => don't report mismatch-with-N differences.
//
// Returns:
//	(nothing)
//
//----------

void print_align_diffs_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	int				withBlocks,
	int				inhibitN)
	{
	unspos			height, width, i, j, run;
	u32				opIx;
	u8*				p, *q;
	s8				b1,  b2;
	unspos			ix;
	int				isMatch, mismatchRun, gapLen;

	height = end1 - beg1;
	width  = end2 - beg2;

	// find and report each difference

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		u32 startI = i;
		u32 startJ = j;

		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i;
		q = seq2->v+beg2+j;
		mismatchRun = 0;
		for (ix=0 ; ix<run ; ix++)
			{
			b1 = nuc_to_bits[*(p++)];
			b2 = nuc_to_bits[*(q++)];

			if (inhibitN) isMatch = (b1 < 0) || (b2 < 0) || (b1 == b2);
			         else isMatch = (b1 == b2);

			if (!isMatch)
				mismatchRun++;
			else if (mismatchRun != 0)
				{
				print_align_difference (f,
				                        seq1, beg1, end1,
				                        seq2, beg2, end2,
				                        script,
				                        i+ix-mismatchRun, p-1-mismatchRun,
				                        j+ix-mismatchRun, q-1-mismatchRun,
				                        mismatchRun,
				                        withBlocks);
				mismatchRun = 0;
				}
			}

		if (mismatchRun != 0)
			print_align_difference (f,
			                        seq1, beg1, end1,
			                        seq2, beg2, end2,
			                        script,
			                        i+ix-mismatchRun, p-mismatchRun,
			                        j+ix-mismatchRun, q-mismatchRun,
			                        mismatchRun,
			                        withBlocks);


		i += run;
		j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i;
			startJ = j;  q = seq2->v+beg2+j;

			edit_script_indel_len (script, &opIx, &i, &j);

			if (i != startI)
				{
				gapLen = i - startI;
				print_align_difference (f,
				                        seq1, beg1, end1,
				                        seq2, beg2, end2,
				                        script,
										i-gapLen, p,
										j, NULL,
										gapLen,
				                        withBlocks);
				p += gapLen;
				}

			if (j != startJ)
				{
				gapLen = j - startJ;
				print_align_difference (f,
				                        seq1, beg1, end1,
				                        seq2, beg2, end2,
				                        script,
										i, NULL,
										j-gapLen, q,
										gapLen,
				                        withBlocks);
				q += gapLen;
				}
			}
		}

	}


static void print_align_difference
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	unspos			diffPos1,
	u8*				diffText1,
	unspos			diffPos2,
	u8*				diffText2,
	unspos			diffLength,
	int				withBlocks)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, run;
	u32				opIx;
	u8*				p, *q;
	unspos			ix;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			startI, startJ;
	unspos			diffLength1, diffLength2;

	height = end1 - beg1;
	width  = end2 - beg2;

	//////////
	// figure out the alignment's length
	//////////

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);
		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			edit_script_indel_len (script, &opIx, &i, &j);
		}

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
		part = lookup_partition (seq1, beg1);
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
		part = lookup_partition (seq2, beg2);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		}

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = beg1 + diffPos1 - offset1 + startLoc1;
		strand1 = '+';
		}
	else
		{
		start1  = beg1 + diffPos1 - offset1 + seq1True+2 - (startLoc1 + seq1Len);
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = beg2 + diffPos2 - offset2 + startLoc2;
		strand2 = '+';
		}
	else
		{
		start2  = beg2 + diffPos2 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		strand2 = '-';
		}

	diffLength1 = (diffText1 != NULL)? diffLength : 0;
	diffLength2 = (diffText2 != NULL)? diffLength : 0;

	//////////
	// print positional information
	//////////

	fprintf (f, "%s\t" unsposFmt "\t" unsposFmt "\t%c\t" unsposFmt "\t",
	            name1, start1-1, start1-1+diffLength1, strand1, seq1True);

	fprintf (f, "%s\t" unsposFmt "\t" unsposFmt "\t%c\t" unsposFmt "\t",
	            name2, start2-1, start2-1+diffLength2, strand2, seq2True);

	// print the aligned difference

	if (diffText1 != NULL)
		{
		for (ix=0 ; ix<diffLength ; ix++)
			fprintf (f, "%c", dna_toprint(diffText1[ix]));
		}
	else
		{
		for (ix=0 ; ix<diffLength ; ix++)
			fprintf (f, "-");
		}

	fprintf (f, "\t");
	if (diffText2 != NULL)
		{
		for (ix=0 ; ix<diffLength ; ix++)
			fprintf (f, "%c", dna_toprint(diffText2[ix]));
		}
	else
		{
		for (ix=0 ; ix<diffLength ; ix++)
			fprintf (f, "-");
		}

	if (!withBlocks)
		goto skip_alignment_block;

	//////////
	// print aligning path in sequence 1
	//////////

	fprintf (f, "\t");

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i;
		q = seq2->v+beg2+j;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*p));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i;
			startJ = j;  q = seq2->v+beg2+j;

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

	//////////
	// print aligning path in sequence 2
	//////////

	fprintf (f, "\t");

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i;
		q = seq2->v+beg2+j;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*q));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i;
			startJ = j;  q = seq2->v+beg2+j;

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

skip_alignment_block:
	fprintf (f, "\n");
	}

//----------
//
// print_align_diffs_match--
//	Print an hsp in a alignment differences.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	seq1:		One sequence.
//	unspos	pos1:		The position, in seq1, of first character in the match
//						.. (origin-0).
//	seq*	seq2:		Another sequence.
//	unspos	pos2:		The position, in seq2, of first character in the match
//						.. (origin-0).
//	unspos	length:		The number of nucleotides in the HSP.
//	int		withBlocks:	true => include full alignment block in output.
//	int		inhibitN:	true => don't report mismatch-with-N differences.
//
// Returns:
//	(nothing)
//
//----------

void print_align_diffs_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	int				withBlocks,
	int				inhibitN)
	{
	u8*				s1, *s2;
	s8				b1,  b2;
	unspos			ix;
	int				isMatch, mismatchRun;

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_align_diffs_match");

	s1 = seq1->v + pos1;
	s2 = seq2->v + pos2;

	// find and report each difference

	mismatchRun = 0;
	for (ix=0 ; ix<length ; ix++)
		{
		b1 = nuc_to_bits[s1[ix]];
		b2 = nuc_to_bits[s2[ix]];

		if (inhibitN) isMatch = (b1 < 0) || (b2 < 0) || (b1 == b2);
		         else isMatch = (b1 == b2);

		if (!isMatch) { mismatchRun++;  continue; }
		if (mismatchRun == 0) continue;

		print_match_difference (f,
		                        seq1, pos1, pos1+ix-mismatchRun,
		                        seq2, pos2, pos2+ix-mismatchRun,
		                        length, mismatchRun,
		                        withBlocks);
		mismatchRun = 0;
		}

	if (mismatchRun != 0)
		print_match_difference (f,
		                        seq1, pos1, pos1+length-mismatchRun,
		                        seq2, pos2, pos2+length-mismatchRun,
		                        length, mismatchRun,
		                        withBlocks);
	}


static void print_match_difference
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	unspos			diffPos1,
	seq*			seq2,
	unspos			pos2,
	unspos			diffPos2,
	unspos			length,
	unspos			diffLength,
	int				withBlocks)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			ix;

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

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = diffPos1 - offset1 + startLoc1;
		strand1 = '+';
		}
	else
		{
		start1  = diffPos1 - offset1 + seq1True+2 - (startLoc1 + seq1Len);
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = diffPos2 - offset2 + startLoc2;
		strand2 = '+';
		}
	else
		{
		start2  = diffPos2 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		strand2 = '-';
		}

	// print positional information

	fprintf (f, "%s\t" unsposFmt "\t" unsposFmt "\t%c\t" unsposFmt "\t",
	            name1, start1-1, start1-1+diffLength, strand1, seq1True);

	fprintf (f, "%s\t" unsposFmt "\t" unsposFmt "\t%c\t" unsposFmt "\t",
	            name2, start2-1, start2-1+diffLength, strand2, seq2True);

	// print the aligned difference

	for (ix=0 ; ix<diffLength ; ix++)
		fprintf (f, "%c", dna_toprint(s1[diffPos1-pos1+ix]));

	fprintf (f, "\t");
	for (ix=0 ; ix<diffLength ; ix++)
		fprintf (f, "%c", dna_toprint(s2[diffPos2-pos2+ix]));

	// print aligning segments

	if (withBlocks)
		{
		fprintf (f, "\t");
		for (ix=0 ; ix<length ; ix++)
			fprintf (f, "%c", dna_toprint(s1[ix]));

		fprintf (f, "\t");
		for (ix=0 ; ix<length ; ix++)
			fprintf (f, "%c", dna_toprint(s2[ix]));
		}
		
	fprintf (f, "\n");
	}

