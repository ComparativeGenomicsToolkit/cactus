//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: maf.h
//
//----------

#ifndef maf_H					// (prevent multiple inclusion)
#define maf_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef maf_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef maf_owner
int maf_distinguishNames = false;	// true => add a "~" prefix to the second
									// sequence name when names are identical
int maf_dbgReportDiag = false;		// true => report diagonal as a maf comment
#else
global int maf_distinguishNames;
global int maf_dbgReportDiag;
#endif

//----------
//
// prototypes for routines in maf.c
//
//----------

void print_maf_job_header (FILE* f,
                           char* programName, char* args, scoreset* scoring,
                           sthresh* hspThreshold, sthresh* gappedThreshold,
                           score xDrop, score yDrop,
                           int withComments);
void print_maf_job_footer (FILE* f);
void print_maf_header     (FILE* f, seq* seq1, seq* seq2);
void print_maf_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                           int withComments);
void print_maf_align      (FILE* f,
                           seq* seq1, unspos beg1, unspos end1,
                           seq* seq2, unspos beg2, unspos end2,
                           editscript* script, score s);
void print_maf_match      (FILE* f,
                           seq* seq1, unspos pos1,
                           seq* seq2, unspos pos2, unspos length,
                           score s, int withComments);
void print_maf_comment    (FILE* f, const char* format, ...);
void vprint_maf_comment   (FILE* f, const char* format, va_list args);

#undef global
#endif // maf_H
