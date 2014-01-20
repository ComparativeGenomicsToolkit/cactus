//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sam.h
//
//----------

#ifndef sam_H					// (prevent multiple inclusion)
#define sam_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// prototypes for routines in sam.c
//
//----------

char* sam_rg_tags          (char* readGroup, char** errorText);
void  print_sam_job_header (FILE* f, char* readGroup);
void  print_sam_header     (FILE* f, seq* seq1, seq* seq2);
void  print_sam_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                            int softMasked, char* rgTags);
void  print_sam_align      (FILE* f,
                            seq* seq1, unspos beg1, unspos end1,
                            seq* seq2, unspos beg2, unspos end2,
                            editscript* script, score s,
                            int softMasked, char* rgTags);
void  print_sam_match      (FILE* f,
                            seq* seq1, unspos pos1,
                            seq* seq2, unspos pos2, unspos length,
                            score s, int softMasked, char* rgTags);

#endif // sam_H
