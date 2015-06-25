//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: maf.h
//
//----------

#ifndef cigar_H					// (prevent multiple inclusion)
#define cigar_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef cigar_owner
#define global
#else
#define global extern
#endif

//----------
//
// prototypes for routines in maf.c
//
//----------

void print_cigar_job_header (FILE* f);
void print_cigar_job_footer (FILE* f);
void print_cigar_header     (FILE* f, seq* seq1, seq* seq2);
void print_cigar_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                             int withInfo, int markMismatches, int letterAfter,
                             int hideSingles, int lowercase, int withNewLine);
void print_cigar_align      (FILE* f,
                             seq* seq1, unspos beg1, unspos end1,
                             seq* seq2, unspos beg2, unspos end2,
                             editscript* script, score s,
                             int withInfo, int markMismatches, int letterAfter,
                             int hideSingles, int lowercase, int withNewLine);
void print_cigar_match      (FILE* f,
                             seq* seq1, unspos pos1,
                             seq* seq2, unspos pos2, unspos length,
                             score s,
                             int withInfo, int markMismatches, int letterAfter,
                             int hideSingles, int lowercase, int withNewLine);

#undef global
#endif // maf_H
