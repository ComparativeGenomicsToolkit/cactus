//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: axt.h
//
//----------

#ifndef axt_H					// (prevent multiple inclusion)
#define axt_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// prototypes for routines in axt.c
//
//----------

void print_axt_job_header (FILE* f,
                           char* programName, char* args, scoreset* scoring,
                           sthresh* hspThreshold, sthresh* gappedThreshold,
                           score xDrop, score yDrop);
void print_axt_job_footer (FILE* f);
void print_axt_header     (FILE* f, seq* seq1, seq* seq2);
void print_axt_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                           int withComments, char* extras);
void print_axt_align      (FILE* f,
                           seq* seq1, unspos beg1, unspos end1,
                           seq* seq2, unspos beg2, unspos end2,
                           editscript* script, score s, char* extras);
void print_axt_match      (FILE* f,
                           seq* seq1, unspos pos1,
                           seq* seq2, unspos pos2, unspos length,
                           score s, int withComments, char* extras);
void print_axt_comment    (FILE* f, const char* format, ...);
void vprint_axt_comment   (FILE* f, const char* format, va_list args);

#endif // axt_H
