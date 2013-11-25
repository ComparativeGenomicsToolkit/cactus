//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: gfa.h
//
//----------

#ifndef gfa_H					// (prevent multiple inclusion)
#define gfa_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// prototypes for routines in gfa.c
//
//----------

void print_gfa_job_header (FILE* f,
                           char* programName, char* name1, char* name2);
void print_gfa_job_footer (FILE* f);
void print_gfa_header     (FILE* f, seq* seq1, seq* seq2);
void print_gfa_align_list (FILE* f, scoreset* scoring, alignel* alignList,
                           seq* seq1, seq* seq2);
void print_gfa_align      (FILE* f, scoreset* scoring,
                           seq* seq1, unspos beg1, unspos end1,
                           seq* seq2, unspos beg2, unspos end2,
                           editscript* script);
void print_gfa_match      (FILE* f,
                           seq* seq1, unspos pos1,
                           seq* seq2, unspos pos2, unspos length,
                           score s);
void print_gfa_generic    (FILE* f, char stanza, const char* format, ...);
void vprint_gfa_generic   (FILE* f, char stanza, const char* format,
                           va_list args);
int  parse_gfa_s_record   (char* rec, char** name1, char** name2);
int  parse_gfa_a_record   (char* rec,
                           unspos* start1, char* strand1,
                           unspos* start2, char* strand2, unspos* length,
                           score* s, int* pctid);

#endif // gfa_H
