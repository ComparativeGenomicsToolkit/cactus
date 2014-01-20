//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: lav.c
//
//----------
//
// lav--
//	Support for printing alignments in LAV format.
//
// LAV format is the well-established pairwise alignment format produced by
// blastz.
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
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "masking.h"			// dynamic masking stuff

#define  lav_owner				// (make this the owner of its globals)
#include "lav.h"				// interface to this module

// debugging defines

//#define snoopLavAlignList		// if this is defined, extra code is added to
								// .. track alignment lists

//----------
//
// print_lav_job_header--
//	Print lav format job header.
//
//----------

void print_lav_job_header
   (FILE*		f,
	char*		_programName,
	char*		_name1,
	char*		_name2,
	char*		_args,
	scoreset*	scoring,
	sthresh*	hspThreshold,
	sthresh*	gappedThreshold,
	u8			dynamicMasking,
	int			withExtras,
	score		xDrop,
	score		yDrop)
	{
	char*		programName = _programName;
	char*		name1       = _name1;
	char*		name2       = _name2;
	char*		args        = _args;

	if (programName == NULL) programName = "(no name)";
	if (name1       == NULL) name1       = "(no name)";
	if (name2       == NULL) name2       = "(no name)";
	if (args        == NULL) args        = "";

	fprintf            (f, "#:lav\n");
	fprintf            (f, "d {\n");
	fprintf            (f, "  \"%s %s %s %s\n", programName, name1, name2, args);
	print_score_matrix (f, scoring, false);
	fprintf            (f, "  O = " scoreFmtSimple
	                       ", E = " scoreFmtSimple
	                       ", K = %s"
	                       ", L = %s"
	                       ", M = %d",
	                       scoring->gapOpen, scoring->gapExtend,
	                       score_thresh_to_string (hspThreshold),
	                       score_thresh_to_string (gappedThreshold),
	                       dynamicMasking);
	if (withExtras)
		fprintf        (f, ", X = " scoreFmtSimple
		                   ", Y = " scoreFmtSimple,
		                   xDrop, yDrop);
	fprintf            (f, "\"\n}\n");
	}

//----------
//
// print_lav_job_footer--
//	Print lav format job footer.
//
//----------

void print_lav_job_footer
   (FILE*	f)
	{
	fprintf (f, "#:eof\n");
	}

//----------
//
// print_lav_header--
//	Print lav format query header.
//
//----------

void print_lav_header
   (FILE*	f,
	seq*	seq1,
	seq*	seq2)
	{
	char*	rcfShortSuffix[4] = { "", "~", "~-", "-" };
	char*	rcfLongSuffix [4] = { "",                        // forward
	                              "~",                       // complement
	                              "~ (reverse complement)",  // reverse
	                              " (reverse complement)" }; // rev-comp
	char*	name1   = seq1->filename;
	char*	name2   = seq2->filename;
	char*	header1 = seq1->header;
	char*	header2 = seq2->header;
	u32		contig1 = seq1->contig;
	u32		contig2 = seq2->contig;

	if (name1   == NULL) name1   = "(no name)";
	if (name2   == NULL) name2   = "(no name)";
	if (header1 == NULL) header1 = "(no header)";
	if (header2 == NULL) header2 = "(no header)";

	fprintf (f, "#:lav\n");
	fprintf (f, "s {\n");
	fprintf (f, "  \"%s%s\" " unsposFmt " " unsposFmt " %d %u\n",
	            name1, rcfShortSuffix[seq1->revCompFlags],
	            seq1->startLoc, seq1->startLoc+seq1->len-1,
	            ((seq1->revCompFlags & rcf_rev) != 0)?1:0, contig1);
	fprintf (f, "  \"%s%s\" " unsposFmt " " unsposFmt " %d %u\n",
	            name2, rcfShortSuffix[seq2->revCompFlags],
	            seq2->startLoc, seq2->startLoc+seq2->len-1,
	            ((seq2->revCompFlags & rcf_rev) != 0)?1:0, contig2);
	fprintf (f, "}\n");

	fprintf (f, "h {\n");
	fprintf (f, "   \"%s%s\"\n", header1, rcfLongSuffix[seq1->revCompFlags]);
	fprintf (f, "   \"%s%s\"\n", header2, rcfLongSuffix[seq2->revCompFlags]);
	fprintf (f, "}\n");
	}

//----------
//
// print_lav_align_list--
//	Print a list of gapped alignments in lav format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//
// Returns:
//	(nothing)
//
//----------

//=== stuff for snoopLavAlignList ===

#ifndef snoopLavAlignList
#define snoopLavAlignList_1 ;
#endif // not snoopLavAlignList

#ifdef snoopLavAlignList

#define snoopLavAlignList_1                                                    \
	fprintf (stderr, "print_lav_align_list  a=%08lX"                           \
	                 "  a->seq1=%08lX  a->seq2=%08lX\n",                       \
	                 (long) a, (long) a->seq1, (long) a->seq2);

#endif // snoopLavAlignList


// print_lav_align_list--

void print_lav_align_list
   (FILE*			f,
	alignel*		alignList,
	seq*			seq1,
	seq*			seq2)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	alignel*		a;

	if ((sp1->p != NULL) || (sp2->p != NULL))
		suicide ("lav format can't handle multi-sequences");
		// the issue is that we'd have to check if the partition changed
		// since the previous alignment, and generate an s/h-stanza pair

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		snoopLavAlignList_1;
		print_lav_align (f,
		                 a->seq1, a->beg1-1, a->end1,
		                 a->seq2, a->beg2-1, a->end2,
		                 a->script, a->s);
		}
	}

//----------
//
// print_lav_align--
//	Print a single gapped alignment in lav format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	const u8*	seq1:		One sequence.
//	unspos		beg1, end1:	Range of positions in sequence 1 (origin 0).
//	const u8*	seq2:		Another sequence.
//	unspos		beg2, end2:	Range of positions in sequence 2 (origin 0).
//	editscript*	script:		The script describing the path the alignment takes
//							.. in the DP matrix.
//	score		s:			The alignment's score.
//
// Returns:
//	(nothing)
//
//----------

static int align_match_percent (unspos run, unspos match);

void print_lav_align
   (FILE*		f,
	const u8*	seq1,
	unspos		beg1,
	unspos		end1,
	const u8*	seq2,
	unspos		beg2,
	unspos		end2,
	editscript*	script,
	score		s)
	{
	unspos		height, width, i, j, prevI, prevJ;
	unspos		run, match;
	u32			opIx;

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	fprintf (f, "a {\n  s " scoreFmtSimple "\n"
	            "  b " unsposFmt " " unsposFmt "\n"
	            "  e " unsposFmt " " unsposFmt "\n",
	            s, beg1, beg2, end1, end2);

	opIx = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		prevI = i;  prevJ = j;
		run = edit_script_run_of_subs_match (script, &opIx,
		                                     seq1+beg1+i-1, seq2+beg2+j-1,
		                                     &match);

		if (run > 0)
			{
			i += run; j += run;
			fprintf (f, "  l " unsposFmt " " unsposFmt
			            " "    unsposFmt " " unsposFmt " %d\n",
			            beg1+prevI, beg2+prevJ, beg1+i-1, beg2+j-1,
			            align_match_percent (run, match));
			}

		if ((i < height) || (j < width))
			edit_script_indel_len (script, &opIx, &i, &j);
		}

	fprintf (f, "}\n");
	}


static int align_match_percent (unspos run, unspos match)
	{
	possum numer, denom;

	if (run == 0) return 0;  // (not clear what should be returned in this case)

	numer = 200 * ((possum) match) + ((possum) run);
	denom =   2 * ((possum) run);

	return numer / denom;	// 100*match/run, rounded
	}

//----------
//
// print_lav_match--
//	Print an hsp in lav format.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos1:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the HSP.
//	score	s:		The HSP's score.
//
// Returns:
//	(nothing)
//
//----------

void print_lav_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	unspos			end1 = pos1 + length;
	unspos			end2 = pos2 + length;
	int				pctId;

	if ((sp1->p != NULL) || (sp2->p != NULL))
		suicide ("lav format can't handle multi-sequences");

	// compute percent identity

	pctId = percent_identical (seq1, pos1, seq2, pos2, length);

	// print it

	fprintf (f, "a {\n");
	fprintf (f, "  s " scoreFmtSimple "\n", s);
	fprintf (f, "  b " unsposFmt " " unsposFmt "\n",
	            pos1+1, pos2+1);
	fprintf (f, "  e " unsposFmt " " unsposFmt "\n",
	            end1,   end2);
	fprintf (f, "  l " unsposFmt " " unsposFmt
	            " "    unsposFmt " " unsposFmt " %d\n",
	            pos1+1, pos2+1, end1, end2, pctId);
	fprintf (f, "}\n");
	}


void print_lavscore_match	// same as regular lav except we output the score
   (FILE*			f,		// .. wherever the pctid would normally go;  this
	seq*			seq1,	// .. is to allow compatibility with some very old
	unspos			pos1,	// .. programs (Dblast, chain)
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	score			s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	unspos			end1 = pos1 + length;
	unspos			end2 = pos2 + length;

	if ((sp1->p != NULL) || (sp2->p != NULL))
		suicide ("lav format can't handle multi-sequences");

	// print it

	fprintf (f, "a {\n");
	fprintf (f, "  s " scoreFmtSimple "\n", s);
	fprintf (f, "  b " unsposFmt " " unsposFmt "\n",              pos1+1, pos2+1);
	fprintf (f, "  e " unsposFmt " " unsposFmt "\n",              end1,   end2);
	fprintf (f, "  l " unsposFmt " " unsposFmt " " unsposFmt " " unsposFmt " " scoreFmtSimple "\n",
	                                        pos1+1, pos2+1, end1, end2, s);
	fprintf (f, "}\n");
	}

//----------
//
// print_lav_m_stanza--
// print_lav_census_stanza--
// print_lav_x_stanza--
//
//----------

static void print_lav_m_interval (unspos b, unspos e, void* info);
static FILE* plmiFile = NULL;

void print_lav_m_stanza (FILE* f, census* cen)
	{
	unspos n;

	plmiFile = f;

	fprintf (f, "m {\n");
	n = 0;
	if (cen != NULL)
		n = report_census_intervals (cen, print_lav_m_interval, NULL);
	fprintf (f, "  n " unsposFmt "\n", n);
	fprintf (f, "}\n");
	}

static void print_lav_m_interval
   (unspos b, unspos e, arg_dont_complain(void* info))
	{ fprintf (plmiFile, "  x " unsposFmt " " unsposFmt "\n", b, e); }


void print_lav_census_stanza (FILE* f, census* cen)
	{
	fprintf      (f, "Census {\n");
	print_census (f, NULL, cen, ' ');
	fprintf      (f, "}\n");
	}


void print_lav_x_stanza (FILE* f, unspos numMasked)
	{ fprintf (f, "x {\n  n " unsposFmt "\n}\n", numMasked); }

//----------
//
// print_lav_comment_open, print_lav_comment_close--
//	Support general comment printing in a lav file.
//
// Note that the lav format does *not* support comments.  Here we use a
// d-stanza, the presence of which may spoil the lav file for some downstream
// tools.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to print to.
//
// Returns:
//	(print_lav_comment_open only) A string which the caller should use as a
//	prefix on all lines within the comment (this may be NULL).
//
//----------

char* print_lav_comment_open
   (FILE* f)
	{
	fprintf (f, "#:lav\n");
	fprintf (f, "d {\n");

	return NULL;
	}


void print_lav_comment_close
   (FILE* f)
	{
	fprintf (f, "}\n");
	}


//----------
//
// print_lav_comment--
//	Print a comment in a lav file.
//
// Note that the lav format does *not* support comments, so these records will
// invariably spoil the lav file for most downstream tools.
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

void print_lav_comment
   (FILE*		f,
	const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);
	vprint_lav_comment (f, format, args);
	va_end (args);
	}

void vprint_lav_comment
   (FILE*		f,
	const char*	format,
	va_list		args)
	{
	fprintf (f, "# ");
	if (format != NULL)
		vfprintf (f, format, args);
	fprintf  (f, "\n");
	}

