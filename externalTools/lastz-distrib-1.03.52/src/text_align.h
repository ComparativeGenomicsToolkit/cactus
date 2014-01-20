//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: text_align.h
//
//----------

#ifndef text_align_H			// (prevent multiple inclusion)
#define text_align_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef text_align_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef text_align_owner
int text_align_dbgReportDiag = false; // true => report diagonal
#else
global int text_align_dbgReportDiag;
#endif

//----------
//
// prototypes for routines in text_align.c
//
//----------

void print_text_align_job_header (FILE* f,
                                  char* programName, char* name1, char* name2,
                                  int oneBased);
void print_text_align_job_footer (FILE* f);
void print_text_align_header     (FILE* f, seq* seq1, seq* seq2,
                                  int oneBased);
void print_text_align_align_list (FILE* f,
                                  alignel* alignList, seq* seq1, seq* seq2,
                                  int oneBased, u32 expand);
void print_text_align_align      (FILE* f,
                                  seq* seq1, unspos beg1, unspos end1,
                                  seq* seq2, unspos beg2, unspos end2,
                                  editscript* script, score s,
                                  int oneBased, u32 expand);
void print_text_align_match      (FILE* f,
                                  seq* seq1, unspos pos1,
                                  seq* seq2, unspos pos2, unspos length,
                                  score s, int oneBased, u32 expand);

#undef global
#endif // text_align_H
