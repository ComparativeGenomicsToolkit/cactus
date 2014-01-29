//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: lav.h
//
//----------

#ifndef lav_H					// (prevent multiple inclusion)
#define lav_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "masking.h"			// dynamic masking stuff
#include "edit_script.h"		// alignment edit script stuff

//----------
//
// prototypes for routines in lav.c
//
//----------

void  print_lav_job_header    (FILE* f,
                               char* programName, char* name1, char* name2,
                               char* args, scoreset* scoring,
                               sthresh* hspThreshold, sthresh* gappedThreshold,
                               u8 dynamicMasking,
                               int withExtras, score xDrop, score yDrop);
void  print_lav_job_footer    (FILE* f);
void  print_lav_header        (FILE* f, seq* seq1, seq* seq2);
void  print_lav_align_list    (FILE* f, alignel* alignList, seq* seq1, seq* seq2);
void  print_lav_align         (FILE* f,
                               const u8* seq1, unspos beg1, unspos end1,
                               const u8* seq2, unspos beg2, unspos end2,
                               editscript* script, score s);
void  print_lav_match         (FILE* f,
                               seq* seq1, unspos pos1,
                               seq* seq2, unspos pos2, unspos length,
                               score s);
void  print_lavscore_match    (FILE* f,
                               seq* seq1, unspos pos1,
                               seq* seq2, unspos pos2, unspos length,
                               score s);
void  print_lav_m_stanza      (FILE* f, census* cen);
void  print_lav_census_stanza (FILE* f, census* cen);
void  print_lav_x_stanza      (FILE* f, unspos numMasked);
char* print_lav_comment_open  (FILE* f);
void  print_lav_comment_close (FILE* f);
void  print_lav_comment       (FILE* f, const char* format, ...);
void  vprint_lav_comment      (FILE* f, const char* format, va_list args);

#endif // lav_H
