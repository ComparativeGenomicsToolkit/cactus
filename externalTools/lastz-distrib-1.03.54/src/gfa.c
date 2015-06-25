//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: gfa.c
//
//----------
//
// gfa--
//	Support for printing alignments in GFA format.
//
// GFA format is for gap-free alignments, with one per line, like the one
// below.  This implemention does not provide enough information that it can be
// converted to a LAV file (full GFA is intended to do so).  This line
// corresponds to an "a" stanza for an alignment starting at 10825 of sequence
// one's + strand and at 8530 on sequence two's - strand (positions are origin
// one), of length 74, score 4137.  GFA's optional percent-identity field is
// not reported.
//
//		a 10825+/8530- 74 4137
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
#include "diag_hash.h"			// diagonals hashing stuff

#define  gfa_owner				// (make this the owner of its globals)
#include "gfa.h"				// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static int is_comment_string (const char* s);

//----------
//
// print_gfa_job_header--
//	Print gfa format job header.
//
//----------

void print_gfa_job_header
   (FILE*		f,
	char*		_programName,
	char*		_name1,
	char*		_name2)
	{
	char*		programName = _programName;
	char*		name1       = _name1;
	char*		name2       = _name2;

	if (programName == NULL) programName = "(no name)";
	if (name1       == NULL) name1       = "(no name)";
	if (name2       == NULL) name2       = "(no name)";

	fprintf (f, "d");
	fprintf (f, " %s %s %s", programName, name1, name2);
	fprintf (f, "\n");
	}

//----------
//
// print_gfa_job_footer--
//	Print gfa format job footer.
//
//----------

void print_gfa_job_footer
   (arg_dont_complain(FILE* f))
	{
	// (do nothing)
	}

//----------
//
// print_gfa_header--
//	Print gfa format header.
//
//----------

void print_gfa_header
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

	fprintf (f, "s ");
	fprintf (f, "\"%s%s\" " unsposFmt " " unsposFmt " %d %u ",
	            name1, rcfShortSuffix[seq1->revCompFlags],
	            seq1->startLoc, seq1->startLoc+seq1->len-1,
	            ((seq1->revCompFlags & rcf_rev) != 0)?1:0, contig1);
	fprintf (f, "\"%s%s\" " unsposFmt " " unsposFmt " %d %u",
	            name2, rcfShortSuffix[seq2->revCompFlags],
	            seq2->startLoc, seq2->startLoc+seq2->len-1,
	            ((seq2->revCompFlags & rcf_rev) != 0)?1:0, contig2);
	fprintf (f, "\n");

	fprintf (f, "h ");
	fprintf (f, "\"%s%s\"",  header1, rcfLongSuffix[seq1->revCompFlags]);
	fprintf (f, " \"%s%s\"", header2, rcfLongSuffix[seq2->revCompFlags]);
	fprintf (f, "\n");
	}

//----------
//
// print_gfa_align_list--
//	Print a list of gapped alignments in gfa format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	scoreset*	scoring:	The scoring scheme to use.  This may be NULL.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//
// Returns:
//	(nothing)
//
//----------

void print_gfa_align_list
   (FILE*		f,
	scoreset*	scoring,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2)
	{
	alignel*	a;

	for (a=alignList ; a!=NULL ; a=a->next)
		print_gfa_align (f, scoring,
		                 seq1, a->beg1-1, a->end1,
		                 seq2, a->beg2-1, a->end2,
		                 a->script);
	}

//----------
//
// print_gfa_align--
//	Print a single gapped alignment in gfa format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	scoreset*	scoring:	The scoring scheme to use.  This may be NULL.
//	seq*		seq1:		One sequence.
//	unspos		beg1, end1:	Range of positions in sequence 1 (origin 0).
//	seq*		seq2:		Another sequence.
//	unspos		beg2, end2:	Range of positions in sequence 2 (origin 0).
//	editscript*	script:		The script describing the path the alignment takes
//							.. in the DP matrix.
//
// Returns:
//	(nothing)
//
//----------

void print_gfa_align
   (FILE*			f,
	scoreset*		scoring,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	unspos			height, width, i, j, prevI, prevJ, run;
	u32				opIx;
	score			s;

	if ((sp1->p != NULL) || (sp2->p != NULL))
		suicide ("gfa format can't handle multi-sequences");
		// the issue is that we'd have to check if the partition changed
		// since the previous alignment, and generate an s/h-stanza pair

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	// print the overall alignment info

	s = 0;
	if (scoring != NULL)
		{
		opIx = 0;
		for (i=j=0 ; (i< height)||(j<width) ; )
			{
			prevI = i;  prevJ = j;
			run = edit_script_run_of_subs (script, &opIx);
			i += run; j += run;
			s += score_match (scoring, seq1, beg1-1+prevI, seq2, beg2-1+prevJ, run);

			if ((i < height) || (j < width))
				{
				run = edit_script_indel_len (script, &opIx, &i, &j);
				if (run > 0)
					s -= scoring->gapOpen + run*scoring->gapExtend;
				}
			}
		}

	fprintf (f, "A " unsposSlashSFmt " " unsposSlashFmt " " scoreFmtSimple "\n",
	            beg1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            beg2, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            height, width, s);

	// print the alignment's segments

	opIx = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		prevI = i;  prevJ = j;
		run = edit_script_run_of_subs (script, &opIx);
		i += run; j += run;

		s = 0;
		if (scoring != NULL)
			s = score_match (scoring, seq1, beg1-1+prevI, seq2, beg2-1+prevJ, run);
		print_gfa_match (f, seq1, beg1-1+prevI, seq2, beg2-1+prevJ, run, s);

		if ((i < height) || (j < width))
			edit_script_indel_len (script, &opIx, &i, &j);
		}

	}

//----------
//
// print_gfa_match--
//	Print a match in gfa format.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The position, in seq1, of first character in the match
//					.. (origin-0).
//	seq*	seq2:	Another sequence.
//	unspos	pos2:	The position, in seq2, of first character in the match
//					.. (origin-0).
//	unspos	length:	The number of nucleotides in the match.
//	score	s:		The match's score.
//
// Returns:
//	(nothing)
//
//----------

void print_gfa_match
   (FILE*	f,
	seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	score	s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	int		pctId;
	sgnpos	diag = diagNumber (pos1, pos2);

	if ((sp1->p != NULL) || (sp2->p != NULL))
		suicide ("gfa format can't handle multi-sequences");

	// compute percent identity

	pctId = percent_identical (seq1, pos1, seq2, pos2, length);

	// print it

	fprintf (f, "a " unsposSlashSFmt " " unsposFmt " " scoreFmtSimple " %d ; diag " sgnposFmt "\n",
	            pos1+1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            pos2+1, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            length, s, pctId, diag);
	}

//----------
//
// print_gfa_generic--
//	Print a generic record in gfa format.
//
// Generic records allow the caller to play fast and loose with the gfa file,
// adding records with stanza-types that most gfa readers will (hopefully)
// ignore.  It can be a good way to pass additional data along to a downstream
// program.
//
//----------
//
// Arguments:
//	FILE*		f:		The file to print to.
//	char		stanza:	The single character identifying some stanza type.
//	const char*	format:	A format string, as per printf.
//	...:				(same as for printf)
//
// Returns:
//	(nothing)
//
//----------

void print_gfa_generic
   (FILE*		f,
	char		stanza,
	const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);
	vprint_gfa_generic (f, stanza, format, args);
	va_end (args);
	}

void vprint_gfa_generic
   (FILE*		f,
	char		stanza,
	const char*	format,
	va_list		args)
	{
	fprintf (f, "%c", stanza);
	if (format != NULL)
		{
		fprintf  (f, " ");
		vfprintf (f, format, args);
		}
	fprintf  (f, "\n");
	}

//----------
//
// parse_gfa_s_record--
//	Parse an s-record from a GFA file.
//
// A-records look something like this.
//
//		s "td/human.fa" 1 1877426 0 1 "td/mouse.fa" 1 1736124 0 1
//
//----------
//
// Arguments:
//	char*	rec:		The record to parse (a zero-terminated string).
//	(remaining arguments are self-explanatory;  each can be NULL if the caller
//   has no interest in that field)
//
// Returns:
//	true if successful, false if not.
//
//----------

int parse_gfa_s_record
   (char*	_rec,
	char**	_name1,
	char**	_name2)
	{
	int		scanned;
	char*	rec = copy_string (_rec);
	char*	s, *n1, *n2;
	unspos	start1, stop1, start2, stop2;
	int		rc1, contig1, rc2, contig2;
	char*	name1, *name2;
	int		len;

	// skip 's'

	s = rec;
	if (*s != 's') goto abort;
	s = skip_whitespace(s+1);

	// parse the first filename

	if (*s != '"') goto abort;
	s++;
	n1 = s;
	s = skip_til (s, "\"");
	if (*s != '"') goto abort;
	*(s++) = 0; // terminate n1
	if (s[-2] == '-') s[-2] = 0;
	s = skip_whitespace(s);

	// parse the four int fields for sequence 1

	scanned = -1;
	sscanf (s, unsposFmtScanf " " unsposFmtScanf " %d %d%n",
	        &start1, &stop1, &rc1, &contig1, &scanned);

	if (scanned == -1)
		goto abort;

	s = skip_whitespace(s+scanned);

	// parse the second filename

	if (*s != '"') goto abort;
	s++;
	n2 = s;
	s = skip_til (s, "\"");
	if (*s != '"') goto abort;
	*(s++) = 0; // terminate n2
	if (s[-2] == '-') s[-2] = 0;
	s = skip_whitespace(s);

	// parse the four int fields for sequence 2

	scanned = -1;
	sscanf (s, unsposFmtScanf " " unsposFmtScanf " %d %d%n",
	        &start2, &stop2, &rc2, &contig2, &scanned);

	if ((scanned == -1) || (!is_blank_string (s+scanned)))
		goto abort;

	// build the file names;  the format is name[s..e][-] where s..e defines a
	// subrange of the file and - indicates reverse-complement

	if (_name1 != NULL)
		{
		len = snprintf (NULL, 0, "%s[" unsposDotsFmt "]%s", n1, start1, stop1, (rc1==1)?"-":"");
		name1 = malloc_or_die ("parse_s_record (name1)", len+1);
		sprintf (name1, "%s[" unsposDotsFmt "]%s", n1, start1, stop1, (rc1==1)?"-":"");
		*_name1 = name1;
		}

	if (_name2 != NULL)
		{
		len = snprintf (NULL, 0, "%s[" unsposDotsFmt "]%s", n2, start2, stop2, (rc2==1)?"-":"");
		name2 = malloc_or_die ("parse_s_record (name2)", len+1);
		sprintf (name2, "%s[" unsposDotsFmt "]%s", n2, start2, stop2, (rc2==1)?"-":"");
		*_name2 = name2;
		}

	// success

	free (rec);
	return true;

	// failure

abort:
	free (rec);
	return false;
	}

//----------
//
// parse_gfa_a_record--
//	Parse an a-record from a GFA file.
//
// A-records look something like this.  The score and pctid fields are optional
// (though of course if you have pctid you have to also have score).
//
//		a start1+/start2+ length score pctid ; comment
//
// - Start1 and start2 are origin-1, relative to the 5's end of their strand.
//   However, the *values* returned by this routine are origin-zero.
//
// - When score is not present, it is assigned noScore.
//
// - When pctid is not present, it is assigned -1.
//
//----------
//
// Arguments:
//	char*	rec:		The record to parse (a zero-terminated string).
//	(remaining arguments are self-explanatory;  each can be NULL if the caller
//   has no interest in that field)
//
// Returns:
//	true if successful, false if not.
//
//----------

int parse_gfa_a_record
   (char*	rec,
	unspos*	_start1,	// (origin-zero)
	char*	_strand1,
	unspos*	_start2,	// (origin-zero)
	char*	_strand2,
	unspos*	_length,
	score*	_s,
	int*	_pctId)
	{
	int		scanned;
	unspos	start1, start2, length;
	char	strand1, strand2;
	int		pctId;
	score	s;

	scanned = -1;
	sscanf (rec, "a " unsposSlashCFmtScanf " " unsposFmtScanf " " scoreFmtScanf " %d%n",
	             &start1, &strand1, &start2, &strand2,
	             &length, &s, &pctId, &scanned);

	if ((scanned != -1) && (!is_comment_string (rec+scanned)))
		return false;

	if (scanned == -1)
		{
		pctId = -1;
		sscanf (rec, "a " unsposSlashCFmtScanf " " unsposFmtScanf " " scoreFmtScanf "%n",
	             &start1, &strand1, &start2, &strand2,
	             &length, &s, &scanned);
		if ((scanned != -1) && (!is_comment_string (rec+scanned)))
			return false;
		}

	if (scanned == -1)
		{
		s = noScore;
		sscanf (rec, "a " unsposSlashCFmtScanf " " unsposFmtScanf "%n",
	             &start1, &strand1, &start2, &strand2,
	             &length, &scanned);
		if ((scanned != -1) && (!is_comment_string (rec+scanned)))
			return false;
		}

	if (scanned == -1)
		return false;

	if ((length <= 0) || (start1 <= 0) || (start2 <= 0)
	 || ((strand1 != '+') && (strand1 != '-'))
	 || ((strand2 != '+') && (strand2 != '-')))
		return false;

	if (_start1  != NULL) *_start1  = start1-1;
	if (_strand1 != NULL) *_strand1 = strand1;
	if (_start2  != NULL) *_start2  = start2-1;
	if (_strand2 != NULL) *_strand2 = strand2;
	if (_length  != NULL) *_length  = length;
	if (_s       != NULL) *_s       = s;
	if (_pctId   != NULL) *_pctId   = pctId;

	return true;
	}

//----------
//
// is_comment_string--
//	Determine if a string contains a comment, or only blank characters.
//
//----------

static int is_comment_string
   (const char*	s)
	{
	char* ss = (char*) s;

	for ( ; *ss!=0 ; ss++)
		{
		if (isspace(*ss)) continue;
		return (*ss == ';');
		}

	return true;
	}

