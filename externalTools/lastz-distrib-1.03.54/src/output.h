//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: output.h
//
//----------

#ifndef output_H				// (prevent multiple inclusion)
#define output_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "masking.h"			// dynamic masking stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// data structures and types
//
//----------

// internal codes for output formats
//
// nota bene: The entries in formatNames[] must 'line up with' the entries in
// this enum.  Thus, if additional output types are added, they must be added
// to both the enum and to formatNames.  fmt_max is not a real format, but can
// be used in the code to determine the largest index into formatNames.
//
// Unfortunately, I haven't figured out how to get the compiler to check for
// any inconsistencies here.

enum
	{
	fmtBad  = -1,
	fmt_min = 0,
	fmtGfa = fmt_min,			// output alignments in GFA format
	fmtGfaNoScore,				// .. in GFA format but without scores
	fmtLav,						// .. in LAV format
	fmtLavComment,				// .. in LAV format with comments
	fmtLavScore,				// .. in LAV format with scores where pctid is
	fmtLavText,					// .. in LAV format with as-text
	fmtAxt,						// .. in AXT format
	fmtAxtComment,				// .. in AXT format with comments
	fmtAxtGeneral,				// .. in AXT format with extra (general) fields
	fmtMaf,						// .. in MAF format
	fmtMafComment,				// .. in MAF format with comments
	fmtMafNoComment,			// .. in MAF format with no comments at all
	fmtSoftSam,					// .. in SAM format, soft masking
	fmtSoftSamNoHeader,			// .. in SAM format, soft masking, no header
	fmtHardSam,					// .. in SAM format, hard masking
	fmtHardSamNoHeader,			// .. in SAM format, hard masking, no header
	fmtCigar,					// .. in standard CIGAR format
	fmtGenpaf,					// .. in 'standard' GENPAF format
	fmtGenpafNoHeader,			// .. in 'standard' GENPAF format, no header
	fmtGenpafNameHeader,		// .. in GENPAF format with names header line
	fmtGenpafBlast,				// .. in 'standard' BLASTN format
	fmtGenpafBlastNoHeader,		// .. in 'standard' BLASTN format, no header
	fmtText,					// .. as text
	fmtZeroText,				// .. as text (zero-based)
	fmtHspComp,					// .. as text, showing composition of each HSP
	fmtDiffs,					// .. as alignment differences
	fmtDiffsNoBlocks,			// .. as alignment differences, without blocks
	fmtInfStats,				// .. as scoring inference stats
	fmtIdDist,					// .. as identity distribution
	fmtDeseed,					// .. as text for deseed program
	fmtInfScores,				// .. collect scoring inference stats (this
								// .. cannot be directly chosen by the user)
	fmtLavInfScores,			// .. fmtLav + fmtInfScores (debugging only)
	fmtNone,					// don't bother to output
	fmt_max = fmtNone
	};

#ifdef output_owner
char* formatNames[] = {"GFA","GFANOSCORE",
                       "LAV","lav+","LAVSCORE","lav+text",
                       "AXT","axt+",NULL,"MAF","maf+","maf-",
                       "sam","sam-","hardsam","hardsam-","cigar",
                       "general","general-",NULL,"blastn","blastn-",
                       "text", "ztext", "comp", "diffs", "diffs-",
                       "infstats","iddist","deseed",
                       "infscores","lav+infscores",
                       "none" };
#else
extern char* formatNames[];
#endif

//----------
//
// prototypes for routines in output.c
//
//----------

void  init_output_for_query     (void);
void  init_output_for_strand    (void);
void  print_align_list_segments (alignel* alignList);
void  print_job_header          (void);
void  print_job_footer          (void);
void  print_header              (void);
void  print_align_list          (alignel* alignList);
char* print_comment_open        (void);
void  print_comment_close       (void);
void  print_eof_comment         (void);
void  print_match               (unspos pos1, unspos pos2, unspos length,
                                 score s);
void  print_m_stanza            (census* cen);
void  print_census_stanza       (census* cen);
void  print_x_stanza            (unspos numMasked);
void  print_generic             (FILE* f, const char* format, ...);

#undef global
#endif // output_H
