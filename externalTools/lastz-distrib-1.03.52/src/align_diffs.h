//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: align_diffs.h
//
//----------

#ifndef align_diffs_H			// (prevent multiple inclusion)
#define align_diffs_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// prototypes for routines in align_diffs.c
//
//----------

void print_align_diffs_job_header (FILE* f,
                                   char* programName, char* name1, char* name2);
void print_align_diffs_job_footer (FILE* f);
void print_align_diffs_header     (FILE* f, seq* seq1, seq* seq2);
void print_align_diffs_align_list (FILE* f,
                                   alignel* alignList, seq* seq1, seq* seq2,
                                   int withBlocks, int inhibitN);
void print_align_diffs_align      (FILE* f,
                                   seq* seq1, unspos beg1, unspos end1,
                                   seq* seq2, unspos beg2, unspos end2,
                                   editscript* script,
                                   int withBlocks, int inhibitN);
void print_align_diffs_match      (FILE* f,
                                   seq* seq1, unspos pos1,
                                   seq* seq2, unspos pos2, unspos length,
                                   int withBlocks, int inhibitN);

#endif // align_diffs_H
