//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: masking.h
//
//----------

#ifndef masking_H				// (prevent multiple inclusion)
#define masking_H

// other files

#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "segment.h"			// segment table management stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef masking_owner
#define global
#else
#define global extern
#endif

//----------
//
// data structures and types
//
//----------

typedef struct census
	{
	unspos	len;			// the length of the sequence being censused
	char	kind;			// the type of counting array
							// .. 'B' => count8  is used
							// .. 'W' => count16 is used
							// .. 'L' => count32 is used
							// .. (anthing else => count8)
	u32		maskThresh;		// the count threshold;  any position with a count
							// .. this high is to be masked
	u8*		count8;			// how many times each base has been part of an
	u16*	count16;		// .. alignment (variable-length array);  the length
	u32*	count32;		// .. of this array is len, and the entries
							// .. correspond to sequence positions 1..len;
							// .. only one of these is active, as indicated by
							// .. census.kind;  the memory for this array is
							// .. allocated as part of this structure
	} census;


// info struct for print_masking_interval (used for the "info" argument)

typedef struct pmiInfo
	{
	FILE*	f;				// the file to write to
	seq*	seq;			// the sequence that has been masked
	} pmiInfo;

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	int   maskedBases;
	} maskingStats;

// stats macros

#define masking_count_stat(field)   ++maskingStats.field
#define masking_uncount_stat(field) --maskingStats.field
#define masking_set_stat(field,val) (maskingStats.field = val)
#define masking_add_stat(field,val) (maskingStats.field += val)
#else
#define masking_count_stat(field)
#define masking_uncount_stat(field)
#define masking_set_stat(field,val)
#define masking_add_stat(field,val)
#endif // collect_stats

// prototypes for stats routines

void masking_zero_stats    (void);
void masking_show_stats    (FILE* f);
void masking_generic_stats (FILE* f, void (*func) (FILE*, const char*, ...));

//----------
//
// prototypes for routines in masking.c
//
//----------

census* new_census               (unspos len, char kind, u32 maskThresh);
unspos  census_mask_segments     (segtable* st, u8* fwd, u8* rev, census* cen,
                                  void(*func) (unspos beg, unspos end, void* info),
                                  void* info);
unspos  census_mask_aligns       (alignel* a, u8* fwd, u8* rev, census* cen,
                                  void(*func) (unspos beg, unspos end, void* info),
                                  void* info);
unspos  count_masked_bases       (seq* _seq, int maskChar);
unspos  report_census_intervals  (census* cen,
                                  void(*func) (unspos beg, unspos end, void* info),
                                  void* info);
unspos  report_masked_intervals  (seq* seq, int maskChar,
                                  void(*func) (unspos beg, unspos end, void* info),
                                  void* info);
void    print_masking_interval   (unspos beg, unspos end, void* info);
void    print_masking_interval_3 (unspos beg, unspos end, void* info);
void    print_census             (FILE* f, seq* _seq, census* cen, char delimiter);

#undef global
#endif // masking_H
