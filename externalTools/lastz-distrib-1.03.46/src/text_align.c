//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: text_align.c
//
//----------
//
// text_align--
//	Support for printing alignments in a textual alignment format.
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
#include "diag_hash.h"			// diagonals hashing stuff

#define  text_align_owner		// (make this the owner of its globals)
#include "text_align.h"			// interface to this module

//----------
//
// private data
//
//----------

#define maxDisplayPerRow     50
#define dnaDisplayPerRow     50
#define quantumDisplayPerRow 20

typedef struct aligndisplay
	{
	FILE*	f;						// the file to print to
	int		displayPerRow;			// number of characters in a displayed row
	int		rev1, rev2;				// true => the corresponding sequence is
									//         .. a reverse-complement
	unspos	beg1, beg2;				// index into seq1 and seq2 of start of
									// .. current line (origin1, inclusive)
	unspos	loc1, loc2;				// current index into seq1 and seq2
	int		ix;						// current index into row1 and row2
	int		quantum1;				// sequence 1 is quantum DNA
	qcode*	qCoding1;				// table to map sequence 1 quantum symbols
									// .. to probabilities (this may be NULL)
	int		quantum2;				// sequence 2 is quantum DNA
	qcode*	qCoding2;				// table to map sequence 2 quantum symbols
									// .. to probabilities (this may be NULL)
	u8		gap1, gap2;				// character for gap in seq1 and seq2
	u8		row1[maxDisplayPerRow+1];	// current row of seq1
	u8		row2[maxDisplayPerRow+1];	// current row of seq2
	} aligndisplay;

//----------
//
// prototypes for private functions
//
//----------

static void          expand_segment (seq* seq1, unspos* pos1,
                                     seq* seq2, unspos* pos2, unspos* length,
                                     u32 expandLeft, u32 expandRight);
static aligndisplay* display_init   (FILE* f,
                                     unspos beg1, int rev1,
                                     unspos beg2, int rev2,
                                     int quantum1, qcode* qCoding1,
                                     int quantum2, qcode* qCoding2);
static void          display_finish (aligndisplay* disp);
static void          display_add    (aligndisplay* disp, u8 ch1, u8 ch2);
static void          display_print  (aligndisplay* disp);

//----------
//
// print_text_align_job_header--
//	Print a textual alignment job header.
//
//----------

void print_text_align_job_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(char* programName),
	arg_dont_complain(char* name1),
	arg_dont_complain(char* name2),
	arg_dont_complain(int	oneBased))
	{
	}

//----------
//
// print_text_align_job_footer--
//	Print a textual alignment job footer.
//
//----------

void print_text_align_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_text_align_header--
//	Print a textual alignment query header.
//
//----------

void print_text_align_header
   (arg_dont_complain(FILE* f),
	arg_dont_complain(seq*  seq1),
	arg_dont_complain(seq*  seq2),
	arg_dont_complain(int	oneBased))
	{
	// (do nothing)
	}

//----------
//
// print_text_align_align_list--
//	Print a list of gapped alignments, textually.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//	int			oneBased:	true  => show positions as origin 1
//							false => show positions as origin 0
//	u32			expand:		Number of extra bp to print at the ends of matches,
//							.. to provide context.
//
// Returns:
//	(nothing)
//
//----------

void print_text_align_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			oneBased,
	u32			expand)
	{
	alignel* a;

	for (a=alignList ; a!=NULL ; a=a->next)
		print_text_align_align (f,
		                        seq1, a->beg1-1, a->end1,
		                        seq2, a->beg2-1, a->end2,
		                        a->script, a->s, oneBased, expand);
	}

//----------
//
// print_text_align_align--
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
//	score		s:			The alignment's score.
//	int			oneBased:	true  => show positions as origin 1
//							false => show positions as origin 0
//	u32			expand:		Number of extra bp to print at the ends of matches,
//							.. to provide context.
//
// Returns:
//	(nothing)
//
//----------
//
// Typical output:
//
//		 61 TCTATCGGTAACCTAATAGA--GACTGAAGCTTACCCCTATGATCTTTGA
//		    |||||||||||||||||||   |||    |  || |  | | |   ||
//		 41 TCTATCGGTAACCTAATAGTTTGACGACGGTGTATCTATGTTAGTGTTAT
//
//		109 CTGGAGTTGTTACGCGATAT-CTTTACCTGTTATCTGGCAC
//		    |||    | | | |  || | ||||||||||||||||||||
//		 91 CTG---CTATAAAGACATCTACTTTACCTGTTATCTGGCAC
//
//----------

void print_text_align_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	score			s,
	int				oneBased,
	u32				expand)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	const int		rev1 = ((seq1->revCompFlags & rcf_rev) != 0);
	const int		rev2 = ((seq2->revCompFlags & rcf_rev) != 0);
	unspos			extra1, extra2;
	u32				expandLeft, expandRight;
	unspos			offset1, offset2;
	unspos			seq1Len, seq2Len;
	unspos			dispBeg1, dispBeg2;
	aligndisplay*	disp;
	unspos			height, width, i, j, run;
	u32				opIx;
	u8*				p, *q;
	u32				ix;
	int				bo = (oneBased)? 0 : -1;

	expandLeft = expandRight = 0;
	if (expand > 0)
		{
		expandLeft = (u32) beg1;
		if (((u32) beg2) < expandLeft) expandLeft = (u32) beg2;
		if (expand       < expandLeft) expandLeft = expand;
		beg1 -= expandLeft;
		beg2 -= expandLeft;

		extra1 = seq1->len - end1;
		extra2 = seq2->len - end2;
		expandRight = (u32) extra1;
		if (((u32) extra2) < expandRight) expandRight = (u32) extra2;
		if (expand         < expandRight) expandRight = expand;
		end1 += expandRight;
		end2 += expandRight;
		}

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	// report diagonal

	if (text_align_dbgReportDiag)
		fprintf (f, "# diagonal=" sgnposFmt "\n", diagNumber(beg1,beg2));

	//////////
	// figure out the alignment's length
	//////////

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);
		if ((i==0) && (j==0))    run += expandLeft;
		if (opIx == script->len) run += expandRight;
		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			edit_script_indel_len (script, &opIx, &i, &j);
		}

	fprintf (f, "score:" scoreFmt " length:(" unsposFmt " " unsposFmt ")\n",
	            s, i, j);

	//////////
	// figure out position offsets
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		offset1 = 0;
		seq1Len = seq1->len;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1);
		offset1 = part->sepBefore + 1;
		seq1Len = part->sepAfter - offset1;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		offset2 = 0;
		seq2Len = seq2->len;
		}
	else					// sequence 2 is partitioned
	 	{
		part    = lookup_partition (seq2, beg2);
		offset2 = part->sepBefore + 1;
		seq2Len = part->sepAfter - offset2;
		}

	//////////
	// draw the alignment (non-printables are printed as '*' but such should
	// never be seen unless there is a problem elsewhere)
	//////////

	dispBeg1 = (rev1)? (seq1Len+1 + bo - beg1) : (beg1 + bo - offset1);
	dispBeg2 = (rev2)? (seq2Len+1 + bo - beg2) : (beg2 + bo - offset2);

	disp = display_init (f, dispBeg1, rev1, dispBeg2, rev2,
	                     (seq1->fileType == seq_type_qdna), seq1->qCoding,
	                     (seq2->fileType == seq_type_qdna), seq2->qCoding);
	if (disp == NULL)
		return;

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		u32 startI = i;
		u32 startJ = j;

		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);
		if ((i==0) && (j==0))    run += expandLeft;
		if (opIx == script->len) run += expandRight;

		p = seq1->v+beg1+i-1;
		q = seq2->v+beg2+j-1;
		for (ix=0 ; ix<run ; ix++)
			{ display_add (disp, dna_toprint(*p), dna_toprint(*q));  p++;  q++; }

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
					{ display_add (disp, dna_toprint(*p), disp->gap2);  p++; }
				}

			if (j != startJ)
				{
				for ( ; startJ<j ; startJ++)
					{ display_add (disp, disp->gap1, dna_toprint(*q));  q++; }
				}
			}
		}

	display_finish (disp);
	}

//----------
//
// print_text_align_match--
//	Print an hsp in a textual alignment.
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
//	score	s:			The HSP's score.
//	int		oneBased:	true  => show positions as origin 1
//						false => show positions as origin 0
//	u32		expand:		Number of extra bp to print at the ends of matches, to
//						.. provide context.
//
// Returns:
//	(nothing)
//
//----------

static void print_quantum_match (FILE* f,
                                 seq* seq1, unspos pos1,
                                 seq* seq2, unspos pos2, unspos length,
                                 score s, int oneBased);

static char quantum_match_char  (qcode* qCoding1, u8 ch1,
                                 qcode* qCoding2, u8 ch2);


void print_text_align_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	int				oneBased,
	u32				expand)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s1, *s2;
	s8				b1,  b2;
	unspos			offset1, offset2;
	unspos			startLoc1, startLoc2;
	char			c;
	unspos			ix;
	int				bo = (oneBased)? 0 : -1;

	if (expand > 0)
		expand_segment (seq1, &pos1, seq2, &pos2, &length, expand, expand);
	s1 = seq1->v + pos1;
	s2 = seq2->v + pos2;

	if ((seq1->fileType == seq_type_qdna)
	 || (seq2->fileType == seq_type_qdna))
		{
		print_quantum_match (f, seq1, pos1, seq2, pos2, length, s, oneBased);
		return;
		}

	// report diagonal

	if (text_align_dbgReportDiag)
		fprintf (f, "# diagonal=" sgnposFmt "\n", diagNumber(pos1,pos2));

	fprintf (f, "score:" scoreFmt " length:" unsposFmt "\n", s, length);

	// figure out position offsets

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part      = lookup_partition (seq1, pos1);
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		}
	else					// sequence 2 is partitioned
	 	{
		part      = lookup_partition (seq2, pos2);
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		}

	// print aligning segment of sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)

	fprintf (f, unsposStarFmt ": ", 10, pos1 + bo - offset1 + startLoc1);
	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s1[ix]));
	fprintf (f, "\n");

	// print match bars

	fprintf (f, "%10s  ", "");
	for (ix=0 ; ix<length ; ix++)
		{
		b1 = nuc_to_bits[s1[ix]];
		b2 = nuc_to_bits[s2[ix]];

		if ((b1 < 0) || (b2 < 0))
			c = ' ';
		else if (b1 == b2)
			c = '|';
		else if (bits_to_pur_pyr[(u8)b1] == bits_to_pur_pyr[(u8)b2])
			c = ':';
		else
			c = ' ';

		fprintf (f, "%c", c);
		}
	fprintf (f, "\n");

	// print aligning segment of sequence 2

	fprintf (f, unsposStarFmt ": ", 10, pos2 + bo - offset2 + startLoc2);
	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toprint(s2[ix]));
	fprintf (f, "\n\n");
	}


void print_quantum_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s,
	int				oneBased)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	int				quantum1 = (seq1->fileType == seq_type_qdna);
	int				quantum2 = (seq2->fileType == seq_type_qdna);
	qcode*			qCoding1 = seq1->qCoding;
	qcode*			qCoding2 = seq2->qCoding;
	partition*		part;
	u8*				s1 = seq1->v + pos1;
	u8*				s2 = seq2->v + pos2;
	unspos			offset1, offset2;
	unspos			startLoc1, startLoc2;
	unspos			ix;
	char3			pField;
	u8				nucIx, ch1, ch2;
	int				bo = (oneBased)? 0 : -1;

	fprintf (f, "score:" scoreFmt " length:" unsposFmt "\n", s, length);

	// figure out position offsets

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part      = lookup_partition (seq1, pos1);
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		}
	else					// sequence 2 is partitioned
	 	{
		part      = lookup_partition (seq2, pos2);
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		}

	// print sequence 1 probabilities

	if (qCoding1 != NULL)
		{
		for (nucIx=0 ; nucIx<sizeof(qCoding1->dna) ; nucIx++)
			{
			fprintf (f, "%10c:", qCoding1->dna[nucIx]);
			for (ix=0 ; ix<length ; ix++)
				{
				ch1 = s1[ix];
				pField = prob_to_string(qCoding1->p[ch1][nucIx]);
				fprintf (f, " %s", pField.s);
				}
			fprintf (f, "\n");
			}
		}

	// print aligning segment of sequence 1

	fprintf (f, unsposStarFmt ":", 10, pos1 + bo - offset1 + startLoc1);
	if (seq1->fileType == seq_type_qdna)
		{ for (ix=0 ; ix<length ; ix++) fprintf (f, "  %02X", s1[ix]); }
	else
		{ for (ix=0 ; ix<length ; ix++) fprintf (f, "  %c ", s1[ix]); }
	fprintf (f, "\n");

	// print the match/mismatch row

	if ((( quantum1) && (qCoding1 != NULL) && ( quantum2) && (qCoding2 != NULL))
	 || ((!quantum1) && (qCoding1 == NULL) && ( quantum2) && (qCoding2 != NULL))
	 || (( quantum1) && (qCoding1 != NULL) && (!quantum2) && (qCoding2 == NULL)))
		{
		fprintf (f, "%10s ", "");
		for (ix=0 ; ix<length ; ix++)
			fprintf (f, "  %c ", quantum_match_char (qCoding1, s1[ix],
			                                         qCoding2, s2[ix]));
		fprintf (f, "\n");
		}

	// print aligning segment of sequence 2

	fprintf (f, unsposStarFmt ":", 10, pos2 + bo - offset2 + startLoc2);
	if (seq2->fileType == seq_type_qdna)
		{ for (ix=0 ; ix<length ; ix++) fprintf (f, "  %02X", s2[ix]); }
	else
		{ for (ix=0 ; ix<length ; ix++) fprintf (f, "  %c ", s2[ix]); }
	fprintf (f, "\n");

	// print sequence 2 probabilities

	if (qCoding2 != NULL)
		{
		for (nucIx=0 ; nucIx<sizeof(qCoding2->dna) ; nucIx++)
			{
			fprintf (f, "%10c:", qCoding2->dna[nucIx]);
			for (ix=0 ; ix<length ; ix++)
				{
				ch2 = s2[ix];
				pField = prob_to_string(qCoding2->p[ch2][nucIx]);
				fprintf (f, " %s", pField.s);
				}
			fprintf (f, "\n");
			}
		}

	fprintf (f, "\n");
	}


static char quantum_match_char
   (qcode*	qCoding1,
	u8		ch1,
	qcode*	qCoding2,
	u8		ch2)
	{
	double	pDiff, pDiffSum;
	int		chIx;
	char*	lookup1, *lookup2;
	u8		nucIx1, nucIx2, ch;

	// if we have no coding, just return a blank

	if ((qCoding1 == NULL) && (qCoding2 == NULL))
		return ' ';

	// if one of the codings is absent, make sure it is qCoding2

	if (qCoding1 == NULL)
		{
		qcode* tqc;
		u8     tch;
		tqc = qCoding1;  qCoding1 = qCoding2;  qCoding2 = tqc;
		tch = ch1;       ch1      = ch2;       ch1      = tch;
		}

	// handle the one coding case

	if (qCoding2 == NULL)
		{
		lookup1 = strchr(qCoding1->dna,(char)ch2);
		if (lookup1 != NULL)
			{
			nucIx1 = lookup1 - qCoding1->dna;
			if      (qCoding1->p[ch1][nucIx1] >= .75) return '|';
			else if (qCoding1->p[ch1][nucIx1] >= .40) return ':';
			}
		return ' ';
		}

	// handle the two coding case

	pDiffSum = 0.0;
	for (chIx=0 ; chIx<4 ; chIx++)
		{
		ch = bits_to_nuc[chIx];
		lookup1 = strchr(qCoding1->dna,(char)ch);
		lookup2 = strchr(qCoding2->dna,(char)ch);
		if ((lookup1 != NULL) && (lookup2 != NULL))
			{
			nucIx1 = lookup1 - qCoding1->dna;
			nucIx2 = lookup2 - qCoding2->dna;
			pDiff = qCoding1->p[ch1][nucIx1] - qCoding2->p[ch2][nucIx2];
			if (pDiff < 0) pDiff = -pDiff;
			pDiffSum += pDiff;
			}
		}

	if      (1-pDiffSum >= .75) return '|';
	else if (1-pDiffSum >= .40) return ':';

	return ' ';
	}

//----------
//
// expand_segment--
//	Expand a segment by adding bases to its ends.
//
//----------
//
// Arguments:
//	seq*	seq1:			One sequence.
//	unspos*	pos1:			The position, in seq1, of first character in the
//							.. match (origin-0).
//	seq*	seq2:			Another sequence.
//	unspos*	pos1:			The position, in seq2, of first character in the
//							.. match (origin-0).
//	unspos*	length:			The number of nucleotides in the HSP.
//	u32		expandLeft:		Number of extra bp to print add at the left end.
//	u32		expandRight:	Number of extra bp to print add at the right end.
//
// Returns:
//  (nothing)
//
//----------

static void expand_segment
   (seq*	seq1,
	unspos*	pos1,
	seq*	seq2,
	unspos*	pos2,
	unspos*	length,
	u32		expandLeft,
	u32		expandRight)
	{
	unspos	beg1 = *pos1;
	unspos	beg2 = *pos2;
	unspos	end1 = beg1 + *length;
	unspos	end2 = beg2 + *length;
	unspos	extra1, extra2;

	if (expandLeft > 0)
		{
		if (beg1 < (unspos) expandLeft) expandLeft = beg1;
		if (beg2 < (unspos) expandLeft) expandLeft = beg2;
		beg1 -= expandLeft;
		beg2 -= expandLeft;
		}

	if (expandRight > 0)
		{
		extra1 = seq1->len - end1;
		extra2 = seq2->len - end2;
		if (extra1 < (unspos) expandRight) expandRight = extra1;
		if (extra2 < (unspos) expandRight) expandRight = extra2;
		end1 += expandRight;
		// end2 += expandRight; (not needed)
		}

	*pos1   = beg1;
	*pos2   = beg2;
	*length = end1 - beg1;
	}

//----------
//
// display_init--
//	Initialize an alignment display.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	unspos	beg1:   	Location of the start of the display, for sequence 1
//						(origin1, inclusive).
//  int		rev1:   	true => sequence 1 is a reverse-complement.
//	unspos	beg2:   	Location of the start of the display, for sequence 2.
//  int		rev2:   	true => sequence 2 is a reverse-complement.
//	int		quantum1:	true => sequence 1 is quantum DNA
//	qcode*	qCoding1:	Table to map sequence 1 quantum symbols to
//						.. probabilities (this may be NULL)
//	int		quantum2:	true => sequence 2 is quantum DNA
//	qcode*	qCoding2:	Table to map sequence 2 quantum symbols to
//						.. probabilities (this may be NULL)
//
// Returns:
//  A pointer to the newly allocated display.  If there is a failure, it is
//  reported to the user and NULL is returned.
//
//----------

static aligndisplay* display_init
   (FILE*	f,
	unspos	beg1,
	int		rev1,
	unspos	beg2,
	int		rev2,
	int		quantum1,
	qcode*	qCoding1,
	int		quantum2,
	qcode*	qCoding2)
	{
	aligndisplay*	disp;

	// allocate the display
	// note that we don't call malloc_or_die here, because we don't want to kill
	// an alignment in this case.

	disp = (aligndisplay*) malloc (sizeof(aligndisplay));
	if (disp == NULL)
		{
		fprintf (stderr, "unable to allocate alignment display for " unsposSlashFmt "\n",
		                 beg1, beg2);
		return NULL;
		}

	// initialize it

	disp->f    = f;
	disp->beg1 = disp->loc1 = beg1;  disp->rev1 = rev1;
	disp->beg2 = disp->loc2 = beg2;  disp->rev2 = rev2;
	disp->ix   = 0;

	disp->gap1 = disp->gap2 = '-';
	disp->displayPerRow = dnaDisplayPerRow;

	disp->quantum1 = quantum1;
	disp->qCoding1 = NULL;
	disp->quantum2 = quantum2;
	disp->qCoding2 = NULL;

	if (quantum1)
		{
		disp->qCoding1 = qCoding1;
		disp->gap1     = 0;
		disp->displayPerRow = quantumDisplayPerRow;
		}

	if (quantum2)
		{
		disp->qCoding2 = qCoding2;
		disp->gap2     = 0;
		disp->displayPerRow = quantumDisplayPerRow;
		}

	return disp;
	}

//----------
//
// display_finish--
//	Finish an alignment display.  The final pair of lines is printed and the
//  display is de-allocated.
//
//----------
//
// Arguments:
//	aligndisplay* disp: The display to finish.
//
// Returns:
//  (nothing)
//
//----------

static void display_finish
   (aligndisplay* disp)
	{
	if (disp->ix > 0)
		{ display_print (disp);  printf ("\n"); }

	free (disp);
	}

//----------
//
// display_add--
//	Add an aligned pair of characters to an alignment display.
//
//----------
//
// Arguments:
//	aligndisplay*	disp:		The display to add to.
//	u8				ch1, ch2:	The aligned characters.  If either of these is
//								.. '-', this is considered an indel.
//
// Returns:
//  (nothing)
//
//----------

static void display_add
   (aligndisplay*	disp,
	u8				ch1,
	u8				ch2)
	{

	// if there's no more room, push stuff out to the console

	if (disp->ix >= disp->displayPerRow)
		display_print (disp);

	// add these characters

	disp->row1[disp->ix] = ch1;
	disp->row2[disp->ix] = ch2;
	disp->ix++;

	// update the sequence positions

	if (ch1 != disp->gap1) { if (disp->rev1) disp->loc1--;
	                                    else disp->loc1++; }
	if (ch2 != disp->gap2) { if (disp->rev2) disp->loc2--;
	                                    else disp->loc2++; }
	}

//----------
//
// display_print--
//	Print one pair of lines of an alignment display.  This should only be
//  called by display_add() or display_finish().
//
//----------
//
// Arguments:
//	aligndisplay* disp: The display to print.
//
// Returns:
//  (nothing)
//
//----------

static void quantum_display_print (aligndisplay* disp);


static void display_print
   (aligndisplay* disp)
	{
	FILE*	f = disp->f;
	int		digits = 10;
	s8		b1, b2;
	char	c;
	int		ix;

	if ((disp->quantum1) || (disp->quantum2))
		{ quantum_display_print (disp);  return; }

	// terminate the lines

	disp->row1[disp->ix] = disp->row2[disp->ix] = 0;

	// print a top spacer

	fprintf (f, "\n");

	// print the top (first) sequence

	fprintf (f, unsposStarFmt " %s\n", digits, disp->beg1, disp->row1);

	// print the match/mismatch row

	fprintf (f, "%*s ", digits, "");
	for (ix=0 ; ix<disp->ix ; ix++)
		{
		b1 = nuc_to_bits[disp->row1[ix]];
		b2 = nuc_to_bits[disp->row2[ix]];

		if ((disp->row1[ix] == disp->gap1)
		 || (disp->row2[ix] == disp->gap2))
			c = '-';
		else if ((b1 < 0) || (b2 < 0))
			c = ' ';
		else if (b1 == b2)
			c = '|';
		else if (bits_to_pur_pyr[(u8)b1] == bits_to_pur_pyr[(u8)b2])
			c = ':';
		else
			c = ' ';

		fprintf (f, "%c", c);
		}
	fprintf (f, "\n");

	// print the bottom (second) sequence

	fprintf (f, unsposStarFmt " %s\n", digits, disp->beg2, disp->row2);

	// prepare for the next line

	disp->beg1 = disp->loc1;
	disp->beg2 = disp->loc2;
	disp->ix = 0;
	}


static void quantum_display_print
   (aligndisplay* disp)
	{
	FILE*	f = disp->f;
	int		quantum1 = disp->quantum1;
	int		quantum2 = disp->quantum2;
	qcode*	qCoding1 = disp->qCoding1;
	qcode*	qCoding2 = disp->qCoding2;
	int		digits = 10;
	int		ix;
	char3	pField;
	u8		nuc, ch1, ch2;

	// print a top spacer

	fprintf (f, "\n");

	// print sequence 1 probabilities

	if (qCoding1 != NULL)
		{
		for (nuc=0 ; nuc<sizeof(qCoding1->dna) ; nuc++)
			{
			fprintf (f, "%*c:", digits, qCoding1->dna[nuc]);
			for (ix=0 ; ix<disp->ix ; ix++)
				{
				ch1 = disp->row1[ix];
				if (ch1 == disp->gap1)
					{ fprintf (f, "  ..");  continue; }
				ch1 = disp->row1[ix];
				if (ch1 == disp->gap1)
					{ fprintf (f, "  ,,");  continue; }
				pField = prob_to_string(qCoding1->p[ch1][nuc]);
				fprintf (f, " %s", pField.s);
				}
			fprintf (f, "\n");
			}
		}

	// print aligning text for sequence 1

	fprintf (f, unsposStarFmt " ", digits, disp->beg1);
	for (ix=0 ; ix<disp->ix ; ix++)
		{
		if (disp->row1[ix] == disp->gap1) fprintf (f, "  --");
		else if (quantum1)                fprintf (f, "  %02X", disp->row1[ix]);
		                             else fprintf (f, "  %c ",  disp->row1[ix]);
		}
	fprintf (f, "\n");

	// print the match/mismatch row

	if ((( quantum1) && (qCoding1 != NULL) && ( quantum2) && (qCoding2 != NULL))
	 || ((!quantum1) && (qCoding1 == NULL) && ( quantum2) && (qCoding2 != NULL))
	 || (( quantum1) && (qCoding1 != NULL) && (!quantum2) && (qCoding2 == NULL)))
		{
		fprintf (f, "%*s ", digits, "");
		for (ix=0 ; ix<disp->ix ; ix++)
			fprintf (f, "  %c ", quantum_match_char (qCoding1, disp->row1[ix],
			                                         qCoding2, disp->row2[ix]));
		fprintf (f, "\n");
		}

	// print aligning text for sequence 2

	fprintf (f, unsposStarFmt " ", digits, disp->beg2);
	for (ix=0 ; ix<disp->ix ; ix++)
		{
		if (disp->row2[ix] == disp->gap2) fprintf (f, "  --");
		else if (quantum2)                fprintf (f, "  %02X", disp->row2[ix]);
		                             else fprintf (f, "  %c ",  disp->row2[ix]);
		}
	fprintf (f, "\n");

	// print sequence 2 probabilities

	if (qCoding2 != NULL)
		{
		for (nuc=0 ; nuc<sizeof(qCoding2->dna) ; nuc++)
			{
			fprintf (f, "%*c:", digits, qCoding2->dna[nuc]);
			for (ix=0 ; ix<disp->ix ; ix++)
				{
				ch1 = disp->row1[ix];
				if (ch1 == disp->gap1)
					{ fprintf (f, "  ..");  continue; }
				ch2 = disp->row2[ix];
				if (ch2 == disp->gap2)
					{ fprintf (f, "  ,,");  continue; }
				pField = prob_to_string(qCoding2->p[ch2][nuc]);
				fprintf (f, " %s", pField.s);
				}
			fprintf (f, "\n");
			}
		}

	// prepare for the next line

	disp->beg1 = disp->loc1;
	disp->beg2 = disp->loc2;
	disp->ix = 0;
	}

