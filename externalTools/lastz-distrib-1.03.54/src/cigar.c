//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: cigar.c
//
//----------
//
// cigar--
//	Support for printing alignments in CIGAR format.
//
// CIGAR format is a pairwise alignment format that describes alignment blocks
// in a run-length format.  As of Jan/2009, a spec for CIGAR files could be
// found at
//	http://may2005.archive.ensembl.org/Docs/wiki/html/EnsemblDocs/CigarFormat.html
// However, as of Jan/2012, that page no longer exists, and it is not known
// where anything like a "spec" can be found.
//
// The treatment of intervals on the - strand is not addressed at the above
// link.  CIGAR is also produced by exonerate (which may be where it origin-
// ated).  The following alignment output from exonerate shows how - strand is
// treated.  It is counted along the + strand, and then listed in reverse order.
//
//	cigar: CAT 11407 11062 - PIG 13828 14153 + 892 M 11 I 1 M 176 I 2 ...
//
//	11407 : TGAGTGTTGAAGTAAACTTGCCAAGTTATCTTTATAGGTATCAGTCCATCGTTAGATTTG : 11348
//	        | | | ||| | || ||||||||||||| ||||||| |||||| |||||| || ||| ||
//	13829 : TTAATTTTGTA-TAGACTTGCCAAGTTAACTTTATATGTATCATTCCATCATTGGATGTG : 13887
//
//	11347 : TTATAGCACACATGCACATTGCTTAGCTAACTGAAACATATCAGAAGAATTTATTATAAT : 11288
//	         |||| ||||| ||| |||| | |||| |||  ||| |||||| |||| |||||||||||
//	13888 : GTATAACACACCTGCTCATTCCCTAGCAAACCCAAAAATATCAAAAGATTTTATTATAAT : 13947
//
//	11287 : GCTTAGATAACTTAGGGAATTGCCTACCAGGAAGTAGTGAATATCCGGAACAAGCTCCTC : 11228
//	        ||||||||||||| | | | |||||| |  | ||||||||| ||  || || ||   | |
//	13948 : GCTTAGATAACTTGGAGCAGTGCCTATCCTGTAGTAGTGAAGATGTGGGACGAGAGTCAC : 14007
//
//	11227 : TATGGACTCCCATAGAGTCAAGGTCAGCTTGATGCTATTATGCTTGCTATGGGCCATCAG : 11168
//	         ||||| |  ||     ||| ||||||||||||||||| |||| | |  |  || |||| 
//	14008 : AATGGATT--CA-----TCAGGGTCAGCTTGATGCTATAATGCCTACCCTAAGCAATCAA : 14060
//
//	11167 : AGAGGGGCCTTTGTAACTTTTAGAGTTTGTGAAGGGACGCCTAAATTTGACCATTATCCT : 11108
//	        | ||            ||||||||||||||||||| || |   ||||  |||||||| | 
//	14061 : ACAG------------CTTTTAGAGTTTGTGAAGGTACACTAGAATTCTACCATTATGCA : 14108
//
//	11107 : GCTTTTCAATAATTTCACATGTAAATCAAAAATAGGGCATATTTA : 11063
//	        || |||||||||| ||| |||| | ||||| || |||||||||||
//	14109 : GCATTTCAATAATATCATATGTGACTCAAAGATTGGGCATATTTA : 14153
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

#define  cigar_owner			// (make this the owner of its globals)
#include "cigar.h"				// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static void print_cigar_mismatchy_run (FILE* f, u8* s1, u8* s2, unspos length,
                                       int letterAfter, int hideSingles, int lowercase);

//----------
//
// print_cigar_job_header--
//	Print cigar format job header.
//
//----------

void print_cigar_job_header
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_cigar_job_footer--
//	Print cigar format job footer.
//
//----------

void print_cigar_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_cigar_header--
//	Print cigar format query header.
//
//----------

void print_cigar_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2))
	{
	// (do nothing)
	}

//----------
//
// print_cigar_align_list--
//	Print a list of gapped alignments in cigar format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	alignel*	alignList:		The list of alignments to print.
//	seq*		seq1:			One sequence.
//	seq*		seq2:			Another sequence.
//	int			withInfo:		true  => include info before cigar path string
//	int			markMismatches:	true  => use =/X syntax for non-indel runs
//								false => use M syntax instead
//	int			letterAfter:	true  => letters after count (defying spec)
//	int			hideSingles:	true  => don't both to print count if count==1
//	int			lowercase:		true  => use lower case letters (defying spec)
//	int			withNewLine:	true  => write new-lines after each alignment
//
// Returns:
//	(nothing)
//
//----------

void print_cigar_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			withInfo,
	int			markMismatches,
	int			letterAfter,
	int			hideSingles,
	int			lowercase,
	int			withNewLine)
	{
	alignel*	a;

	for (a=alignList ; a!=NULL ; a=a->next)
		print_cigar_align (f,
		                   seq1, a->beg1-1, a->end1,
		                   seq2, a->beg2-1, a->end2,
		                   a->script, a->s,
		                   withInfo, markMismatches, letterAfter,
		                   hideSingles, lowercase, withNewLine);
	}

//----------
//
// print_cigar_align--
//	Print a single gapped alignment in cigar format.
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
//	int			withInfo:		true  => include info before cigar path string
//	int			markMismatches:	true  => use =/X syntax for non-indel runs
//								false => use M syntax instead
//	int			letterAfter:	true  => letters after count (defying spec)
//	int			hideSingles:	true  => don't both to print count if count==1
//	int			lowercase:		true  => use lower case letters (defying spec)
//	int			withNewLine:	true  => write new-line after the alignment
//
// Returns:
//	(nothing)
//
//----------

static char* rcfSuffix[4] = { "", "~", "~", "" };

void print_cigar_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			_end1,
	seq*			seq2,
	unspos			beg2,
	unspos			_end2,
	editscript*		script,
	score			s,
	int				withInfo,
	int				markMismatches,
	int				letterAfter,
	int				hideSingles,
	int				lowercase,
	int				withNewLine)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + beg1;
	u8*				s2 = seq2->v + beg2;
	unspos			height, width, i, j, prevI, prevJ, run;
	u32				opIx;
	char*			name1, *name2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2, end1, end2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len;
	char			strand1, strand2;
	char			chM = (lowercase)? 'm' : 'M';
	char			chD = (lowercase)? 'd' : 'D';
	char			chI = (lowercase)? 'i' : 'I';

	height = _end1 - beg1;
	width  = _end2 - beg2;

	// figure out position offsets and names

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		seq1Len   = seq1->len;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		seq1Len   = part->sepAfter - offset1;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2Len   = seq2->len;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, beg2);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2Len   = part->sepAfter - offset2;
		}

	// figure out strandedness

	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = beg1-1 - offset1 + startLoc1;
		end1    = start1 + height;
		strand1 = '+';
		}
	else
		{
		start1  = startLoc1 + seq1Len + offset1 - (beg1+1);
		end1    = start1 - height;
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = beg2-1 - offset2 + startLoc2;
		end2    = start2 + width;
		strand2 = '+';
		}
	else
		{
		start2  = startLoc2 + seq2Len + offset2 - (beg2+1);
		end2    = start2 - width;
		strand2 = '-';
		}

	// print the alignment

	if (withInfo)
		fprintf (f, "cigar:"
		            " %s%s " unsposFmt " " unsposFmt " %c"
		            " %s%s " unsposFmt " " unsposFmt " %c"
		            " " scoreFmt,
		            name2, suff2, start2, end2,  strand2,
		            name1, suff1, start1, end1, strand1,
		            s);

	opIx = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		run = edit_script_run_of_subs (script, &opIx);
		if (run > 0)
			{
			if (markMismatches)
				print_cigar_mismatchy_run (f, s1+i, s2+j, run,
				                           letterAfter, hideSingles, lowercase);
			else
				{
				if (letterAfter) fprintf (f, unsposFmt "%c",   run, chM);
				            else fprintf (f, " %c " unsposFmt, chM, run);
				}
			i += run; j += run;
			}

		if ((i < height) || (j < width))
			{
			prevI = i;  prevJ = j;
			edit_script_indel_len (script, &opIx, &i, &j);
			if (i > prevI)
				{
				if (!letterAfter)
					fprintf (f, " %c " unsposFmt, chD, i - prevI);
				else if ((hideSingles) && (i - prevI == 1))
					fprintf (f, "%c", chD);
				else
					fprintf (f, unsposFmt "%c", i - prevI, chD);
				}
			if (j > prevJ)
				{
				if (!letterAfter)
					fprintf (f, " %c " unsposFmt, chI, j - prevJ);
				else if ((hideSingles) && (j - prevJ == 1))
					fprintf (f, "%c", chI);
				else
					fprintf (f, unsposFmt "%c", j - prevJ, chI);
				}
			}
		}

	if (withNewLine)
		fprintf (f, "\n");
	}

//----------
//
// print_cigar_match--
//	Print an hsp in cigar format.
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
//	int		withInfo:		true  => include info before cigar path string
//	int		markMismatches:	true  => use =/X syntax for non-indel runs
//							false => use M syntax instead
//	int		letterAfter:	true  => letters after count (defying spec)
//	int		hideSingles:	true  => don't both to print count if count==1
//	int		lowercase:		true  => use lower case letters (defying spec)
//	int		withNewLine:	true  => write new-line after the alignment
//
// Returns:
//	(nothing)
//
//----------

void print_cigar_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	int				withInfo,
	int				markMismatches,
	int				letterAfter,
	int				hideSingles,
	int				lowercase,
	int				withNewLine)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	char*			name1, *name2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2, end1, end2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len;
	char			strand1, strand2;
	char			chM = (lowercase)? 'm' : 'M';

	// figure out position offsets and names

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		seq1Len   = seq1->len;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, pos1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		seq1Len   = part->sepAfter - offset1;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2Len   = seq2->len;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, pos2);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2Len   = part->sepAfter - offset2;
		}

	// figure out strandedness

	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = pos1-1 - offset1 + startLoc1;
		end1    = start1 + length;
		strand1 = '+';
		}
	else
		{
		start1  = startLoc1 + seq1Len + offset1 - (pos1+1);
		end1    = start1 - length;
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = pos2-1 - offset2 + startLoc2;
		end2    = start2 + length;
		strand2 = '+';
		}
	else
		{
		start2  = startLoc2 + seq2Len + offset2 - (pos2+1);
		end2    = start2 - length;
		strand2 = '-';
		}

	// print the alignment

	if (withInfo)
		fprintf (f, "cigar:"
		            " %s%s " unsposFmt " " unsposFmt " %c"
		            " %s%s " unsposFmt " " unsposFmt " %c"
		            " " scoreFmt,
		            name2, suff2, start2, end2, strand2,
		            name1, suff1, start1, end1, strand1,
		            s);

	if (markMismatches)
		print_cigar_mismatchy_run (f, s1, s2, length,
		                           letterAfter, hideSingles, lowercase);
	else
		{
		if (!letterAfter)
			fprintf (f, " %c " unsposFmt, chM, length);
		else if ((hideSingles) && (length == 1))
			fprintf (f, "%c", chM);
		else
			fprintf (f, unsposFmt "%c", length, chM);
		}

	if (withNewLine)
		fprintf (f, "\n");
	}

//----------
//
// print_cigar_mismatchy_run--
//	Print a run of match/mismatch in (new) cigar format, using =/X syntax.
//
//----------
//
// Arguments:
//	FILE*	f:				The file to print to.
//	u8*  	s1:				Start of the run in one sequence.
//	u8*  	s2:				Start of the run in the other sequence.
//	unspos	length:			The number of nucleotides in the run.
//	int		letterAfter:	true => letters after count (defying spec)
//	int		hideSingles:	true => don't both to print count if count==1
//	int		lowercase:		true => use lower case letters (defying spec)
//
// Returns:
//	(nothing)
//
//----------

static void print_cigar_mismatchy_run
   (FILE*	f,
	u8*  	s1,
	u8*  	s2,
	unspos	length,
	int		letterAfter,
	int		hideSingles,
	int		lowercase)
	{
	char	chX = (lowercase)? 'x' : 'X';

	unspos	ix;
	int		runIsMm, runLen;
	s8		b1, b2;
	char	ch;

	//for (ix=0 ; ix<length ; ix++) fprintf (f, "%c", s1[ix]);  fprintf (f, "~");
	//for (ix=0 ; ix<length ; ix++) fprintf (f, "%c", s2[ix]);  fprintf (f, "~");

	runIsMm = false;
	runLen  = 0;
	for (ix=0 ; ix<length ; ix++)
		{
		b1 = nuc_to_bits[s1[ix]];
		b2 = nuc_to_bits[s2[ix]];
		if ((b1 == b2) && (b1 >= 0)) // match
			{
			if (!runIsMm) { runLen++;  continue; }
			if (runLen > 0)
				{
				if (!letterAfter)
					fprintf (f, " %c " unsposFmt, chX, runLen);
				else if ((hideSingles) && (runLen == 1))
					fprintf (f, "%c", chX);
				else
					fprintf (f, unsposFmt "%c", runLen, chX);
				}
			runIsMm = false;
			runLen  = 1;
			}
		else						// mismatch
			{
			if (runIsMm) { runLen++;  continue; }
			if (runLen > 0)
				{
				if (!letterAfter)
					fprintf (f, " = " unsposFmt, runLen);
				else if ((hideSingles) && (runLen == 1))
					fprintf (f, "=");
				else
					fprintf (f, unsposFmt "=", runLen);
				}
			runIsMm = true;
			runLen  = 1;
			}
		}

	if (runLen > 0)
		{
		ch = (runIsMm)? chX : '=';
		if (!letterAfter)
			fprintf (f, " %c " unsposFmt, ch, runLen);
		else if ((hideSingles) && (runLen == 1))
			fprintf (f, "%c", ch);
		else
			fprintf (f, unsposFmt "%c", runLen, ch);
		}

	}

