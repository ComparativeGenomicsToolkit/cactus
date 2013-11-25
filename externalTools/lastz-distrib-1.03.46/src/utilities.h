//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: utilities.h
//
//----------

#ifndef utilities_H				// (prevent multiple inclusion)
#define utilities_H

#include <stdio.h>				// standard C i/o stuff

// GNU compiler version

#ifdef __GNUC__
#define GCC_VERSION (10000 * __GNUC__       \
				   +   100 * __GNUC_MINOR__ \
				   +         __GNUC_PATCHLEVEL__)
#endif

// establish ownership of global variables

#ifdef utilities_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef utilities_owner
int utilities_dbgDumpFilePointers = false; // true => dump file pointers in, e.g, fopen_or_die() and getc_or_die()
#else
global int utilities_dbgDumpFilePointers;
#endif

//----------
//
// data structures and types
//
//----------

// sized data types;  these generally come from stdint.h but on some older
// platforms it may not exist, in which case override_stdint can be enabled
// and the proper types set up here

#ifdef override_stdint

typedef signed   char       int8;
typedef signed   char       s8;
typedef unsigned char       u8;
typedef          short int  int16;
typedef          short int  s16;
typedef unsigned short int  u16;
typedef          long  int  int32;
typedef          long  int  s32;
typedef unsigned long  int  u32;
typedef          long  long int64;
typedef          long  long s64;
typedef unsigned long  long u64;

#else

#include <stdint.h>
typedef int8_t   int8;
typedef int8_t   s8;
typedef uint8_t  u8;
typedef int16_t  int16;
typedef int16_t  s16;
typedef uint16_t u16;
typedef int32_t  int32;
typedef int32_t  s32;
typedef uint32_t u32;
typedef int64_t  int64;
typedef int64_t  s64;
typedef uint64_t u64;

#endif // override_stdint

#define u8max  255U
#define u16max 65535U
#define u32max 4294967295U
#define u64max 18446744073709551615LLU

#define s8max  127
#define s16max 32767
#define s32max 2147483647
#define s64max 9223372036854775807LL

// short strings

typedef struct char3 { char s[4]; } char3;

// macro to round data structure sizes to the next larger multiple of 8, 16,
// or m

#define round_up_2(b)   ((((u64) (b))+1)&(~1))
#define round_up_4(b)   ((((u64) (b))+3)&(~3))
#define round_up_8(b)   ((((u64) (b))+7)&(~7))
#define round_up_16(b)  ((((u64) (b))+15)&(~15))
#define round_up_32(b)  ((((u64) (b))+31)&(~31))
#define round_up_64(b)  ((((u64) (b))+63)&(~63))
#define round_up_128(b) ((((u64) (b))+127)&(~127))
#define round_up_256(b) ((((u64) (b))+255)&(~255))
#define round_up_512(b) ((((u64) (b))+511)&(~511))
#define round_up_1K(b)  ((((u64) (b))+1023)&(~1023))
#define round_up_2K(b)  ((((u64) (b))+2047)&(~2047))
#define round_up_4K(b)  ((((u64) (b))+4095)&(~4095))
#define round_up_8K(b)  ((((u64) (b))+8191)&(~8191))
#define round_up_16K(b) ((((u64) (b))+16383)&(~16383))
#define round_up(b,m)   (((((u64) (b))+((m)-1))/(m))*(m))

// macro to count the number of entries in a staticly declared array

#define entriesof(array) (sizeof(array)/sizeof((array)[0]))

// silly type check defeaters

#define ustrlen(s)     (strlen((char*)(s)))
#define ustrcmp(s1,s2) (strcmp((char*)(s1),(char*)(s2)))
#define ustrcpy(s1,s2) (strcpy((char*)(s1),(char*)(s2)))
#define ustrchr(s,c)   (strchr((char*)(s),(char)(c)))

#define strleni(s)     ((int)(strlen(s)))

// macro to convince gnu c compiler not to complain about unusued function
// arguments

#ifdef __GNUC__
#define arg_dont_complain(arg) arg __attribute__ ((unused))
#else
#define arg_dont_complain(arg) arg
#endif // __GNUC__

// printf macros for sized integers

#ifdef override_inttypes
#define s64Fmt  "%jd"
#define u64Fmt  "%ju"
#define u64xFmt "%jX"
#else
#include <inttypes.h>
#define s64Fmt  "%" PRId64
#define u64Fmt  "%" PRIu64
#define u64xFmt "%" PRIX64
#endif // override_inttypes

//----------
//
// memory allocation routines in utilities.c
//
// These routines wrap the memory allcoation routines in the standard library,
// and permit us to use some post-processing tools to inspect memory usage.
// They are affected by two compile-time #defines:
//
//		trackMemoryUsage:	If defined, the routines will write a detailed
//							memory de/re/allocation history to stderr, which
//							can be processed by memory_sniffer (part of the
//							lastz tools).
//
//		noMemoryWrappers:	If defined, calls to the memory wrappers are
//							replaced by calls to the standard library routines.
//							This is useful in conjunction with Valgrind's
//							heap profiler, Massif.
//
//----------

#ifdef noMemoryWrappers
#define malloc_or_die(id,size)    malloc (size)
#define zalloc_or_die(id,size)    calloc (1,size)
#define realloc_or_die(id,p,size) realloc (p,size)
#define free_if_valid(id,p)       free (p)
#endif // noMemoryWrappers

#ifndef noMemoryWrappers
void* malloc_or_die  (char* id, size_t size);
void* zalloc_or_die  (char* id, size_t size);
void* realloc_or_die (char* id, void* _p, size_t size);
void  free_if_valid  (char* id, void* p);
#endif // not noMemoryWrappers


#ifdef trackMemoryUsage
#define memory_checkpoint(fmt)       fprintf(stderr,fmt)
#define memory_checkpoint_1(fmt,i)   fprintf(stderr,fmt,i)
#define memory_checkpoint_2(fmt,i,s) fprintf(stderr,fmt,i,s)
#endif // trackMemoryUsage

#ifndef trackMemoryUsage
#define memory_checkpoint(fmt)       ;
#define memory_checkpoint_1(fmt,i)   ;
#define memory_checkpoint_2(fmt,i,s) ;
#endif // not trackMemoryUsage

//----------
//
// malloc sizing range
//	(see also "sequence sizing types" in sequences.h)
//	These types control the range of dynamically allocated block sizes we can
//	handle.
//
//	Allocation lengths are normally assumed to be small enough to fit into a
//	32-bit integer.  This gives a maximum length of about 4 billion bytes.
//	Since the biggest allocation expected is four bytes per each base in a
//	sequence, this limits the maximum sequence to about 1 billion bp (long
//	even for possum chromosomes).  The programmer can override this at compile
//	time by defining max_malloc_index as 31 or 40 (we also allow 20, but that
//	is only to test whether the mechanism actually works).
//
//----------

#define mallocOverhead 16

#if defined(max_malloc_index)
#define maxMallocIndex max_malloc_index
#else
#define maxMallocIndex 32
#endif

#if (maxMallocIndex == 31)
#define mallocLimit ((u32max/2)-mallocOverhead)
#elif (maxMallocIndex == 32)
#define mallocLimit (u32max-mallocOverhead)
#elif (maxMallocIndex == 40)
#define mallocLimit (1099511627776LLU-mallocOverhead)
#elif (maxMallocIndex == 20) // for debug only
#define mallocLimit ((u32max/(1<<12))-mallocOverhead)
#else
#error ***** undecipherable max malloc length definition *****
#endif

//----------
//
// prototypes for routines in utilities.c
//
//----------

FILE*  fopen_or_die              (const char* name, const char* mode);
int    fclose_if_valid           (FILE* f);
int    getc_or_die               (FILE* f, char* filename);
int    print_prefix              (FILE* f, const char* s, int n);
char*  copy_string               (const char* s);
char*  copy_prefix               (const char* s, int n);
char*  concatenate_strings       (const char* s1, const char* s2);
void   append_char               (char** s, u32* size, u32* len, char ch);
void   append_u8                 (u8**   s, u32* size, u32* len, u8   ch);
int    strcmp_prefix             (const char* str1, const char* str2);
int    strcmp_suffix             (const char* str1, const char* str2);
int    strncmp_suffix            (const char* str1, const char* str2, size_t n);
int    is_blank_string           (const char* s);
int    string_to_int             (const char* s);
int    string_to_unitized_int    (const char* s, int byThousands);
int64  string_to_unitized_int64  (const char* s, int byThousands);
int    hex_string_to_int         (const char* s);
double string_to_double          (const char* s);
double string_to_unitized_double (const char* s, int byThousands);
double pct_string_to_double      (const char* s);
char3  prob_to_string            (double p);
int    string_replace            (char* s, int len, char* sub, char* rep);
char*  trim_string               (char* s);
char*  skip_whitespace           (char* s);
char*  skip_darkspace            (char* s);
char*  skip_til                  (char* s, char* chars);
char*  skip_while                (char* s, char* chars);
char*  find_tabbed_tag           (char* s, char* tag);
int    tabbed_tag_length         (char* tag);
int    is_valid_lastz_version    (char* s);
int    is_later_lastz_version    (char* s1, char* s2);
char*  commatize                 (const int64 v);
char*  ucommatize                (const u64 v);
char*  unitize                   (const int64 v, int byThousands);
char*  hex_64_string             (const int64 v);
u64    swap_64_halves            (const u64 v);
u64    swap_two32_endian         (const u64 v);
u32    swap_32_endian            (const u32 v);
int    bit_count                 (u32 bits);
int    bit_count_64              (u64 bits);
int    bit_count_16              (u32 bits);
int    bit_count_8               (u8 bits);
u32    hassock_hash              (const void* key, u32 len);
void   suicide                   (const char* message);
void   suicidef                  (const char* format, ...);
void   suicide_with_perror       (const char* message);
void   suicidef_with_perror      (const char* format, ...);

#undef global
#endif // utilities_H
