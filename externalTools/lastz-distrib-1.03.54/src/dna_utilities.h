//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: dna_utilities.h
//
//----------

#ifndef dna_utilities_H			// (prevent multiple inclusion)
#define dna_utilities_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <float.h>				// (for FLT_MIN)
#include "utilities.h"			// utility stuff

// establish ownership of global variables

#ifdef dna_utilities_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef dna_utilities_owner
int dna_utilities_dbgShowQToBest = false;	// true => report quantum scores
											//         .. qToBest to stderr
#else
global int dna_utilities_dbgShowQToBest;
#endif

//----------
//
// score values--
//	Scores used for sequence comparisons are normally signed 32-bit integers,
//	but the programmer can override this at compile time by defining score_type
//	as one of 'F', 'D', or 'I'.  Note that some effort must be taken to get
//	the apostrophes into the definition from the compiler command line, such as
//	-Dscore_type=\'F\' .
//
//----------
//
// Notes:
//
//	(1)	For 32-bit integers, W (the worst possible score) is about minus 2.1
//		billion (-2,100,000,000).
//
//	(2)	We set negInfinity to a value higher than the worst possible (i.e. less
//		negative.  This permits values *less* than negative infinity, so that
//		we can subtract values from negative infinity without running the risk
//		that the result underflows (and becomes positive).  By setting this to
//		90% of W, we are safe as long as valid scores for pieces of alignments
//		do not get worse than -W/10.  For 32-bit integers, negInfinity is about
//		-1,900,000,000 and -W/10 is about -200,000,000.  We expect gap and
//		substitution scores to be much smaller than this.  As long as xdrop and
//		ydrop are smaller than this, we should be OK.  Typically we expect
//		substitution scores in the range of hundreds.
//
//	(3)	We set veryBadScore to a negative value such that, when added to
//		negative infinity, we will not get underflow.  We also want veryBadScore
//		to be so bad that no alignment can possibly include a substition with
//		that score.  We set this to half the difference between W and negative
//		infinity.  On 32 bit machines this is about -100,000,000
//
//----------

#if defined(score_type)
#define scoreType score_type
#else
#define scoreType 'I'
#endif

#if (scoreType == 'I')
typedef s32 score;
#elif (scoreType == 'F')
typedef float score;
#elif (scoreType == 'D')
typedef double score;
#else
#error ***** undecipherable score type definition *****
#endif

#ifdef dna_utilities_owner
char dna_utilities_scoreType = scoreType;
#else
extern char dna_utilities_scoreType;
#endif

#if (scoreType == 'I')
#ifdef override_inttypes
#define scoreFmt       "%d"
#define scoreFmtStar   "%*d"
#define scoreFmtSimple "%d"
#define scoreFmtScanf  "%d"
#else
#include <inttypes.h>
#define scoreFmt       "%"  PRId32
#define scoreFmtStar   "%*" PRId32
#define scoreFmtSimple "%"  PRId32
#define scoreFmtScanf  "%"  SCNd32
#endif // override_inttypes
#elif (scoreType == 'F')
#define scoreFmt       "%f"
#define scoreFmtStar   "%*f"
#define scoreFmtSimple "%f"
#define scoreFmtScanf  "%f"
#elif (scoreType == 'D')
#define scoreFmt       "%le"
#define scoreFmtStar   "%*.6le"
#define scoreFmtSimple "%lf"
#define scoreFmtScanf  "%le"
#endif


#if (scoreType == 'I')
#define worstPossibleScore (-0x7FFFFFFF-1)
#define bestPossibleScore    0x7FFFFFFF
#elif (scoreType == 'F') || (scoreType == 'D')
#define worstPossibleScore -FLT_MAX		// nota bene: FLT_MIN wasn't what I
#define bestPossibleScore  FLT_MAX		// .. thought it was
#endif

#define noScore      worstPossibleScore							// see note (1)
#define negInfinity  ((score) (0.9*worstPossibleScore))			// see note (2)
#define veryBadScore (-((negInfinity-worstPossibleScore)/2))	// see note (3)

// $$$ The following is problematic.  Some modules need to know the size of
//     the score type, at *compile* time.  Standard C doesn't seem to provide
//     this capability, as sizeof seems to not work in the preprocessor

#if (scoreType == 'I')
#define score_sz 4
#elif (scoreType == 'F')
#define score_sz 4
#elif (scoreType == 'D')
#define score_sz 8
#endif

//----------
//
// data structures and types
//
//----------

// scoring sets--
//	A score set defines how alignments are scored.  It includes a scoring
//	matrix for substitutions and a linear gap penalty.
//
//	The matrix is implemented as a fixed 256x256 array.  Row and column zero of
//	are reserved for scores so bad that nothing will match to a NUL character.

typedef score scorerow[256];

typedef struct charvec
	{
	s8			len;			// number of characters in the vector;  -1;
								// .. indicates an empty vector
	u8			v[4];			// two-bit codes of some of the characters in
								// .. the bottleneck alphabet
	} charvec;

typedef struct scoreset
	{
	u8			rowChars[256];	// the characters we expect to see used as
	u8			colChars[256];	// .. row (sequence1) and column (sequence2)
								// .. indexes;  these are zero-terminated
								// .. strings;  all other indexes lead to
								// .. 'background' scores
	int			rowsAreDna;		// true => row indexes are A,C,G,T
	int			colsAreDna;		// true => column indexes are A,C,G,T
	u8			badRow, badCol;	// the characters used for the bad scoring
								// .. row and column
	int			gapOpenSet;		// true => gapOpen was set explicitly (as
								// .. opposed to default value)
	score		gapOpen;		// (non-negative) penalty for opening an
								// .. alignment gap;  note that we also apply
								// .. gapExtend on open
	int			gapExtendSet;	// true => gapExtend was set explicitly
	score		gapExtend;		// (non-negative) penalty for extending an
								// .. alignment gap
	u8			bottleneck[5];	// bottleneck alphabet for quantum DNA (includes
								// a zero terminator);  only valid if rowsAreDna
								// .. is false
	charvec		qToBest[256];	// array to map a quantum base in the row
								// .. alphabet to the two-bit codes for the
								// .. 'closest' character(s) in the bottleneck
								// .. alphabet;  indexes that are not a valid
								// .. row character have qToBest[].len == -1;
								// .. only valid if the row alphabet is not
								// .. ACGT
	u8*			qToComplement;	// (similar to nuc_to_complement) array to
								// .. map a quantum base in the column alphabet
								// .. to its complement;  this may be NULL
	scorerow	sub[256];		// maps a nucleotide pair to a score;  indexed
								// .. by [c1][c2], where c1 and c2 are
								// characters (in the range 0..255) from
								// sequence 1 and 2, respectively
	} scoreset;

// extended score set, extended to include related command-line parameters;
// note that ss field must be first, so that we can safely typecast this
// structure to a scoreset

typedef struct exscoreset
	{							//  .. command-line parameters)
	scoreset	ss;
	int			xDropSet;
	score		xDrop;
	int			yDropSet;
	score		yDrop;
	int			stepSet;
	u32			step;
	int			hspThresholdSet;
	score		hspThreshold;
	int			gappedThresholdSet;
	score		gappedThreshold;
	int			ballScoreSet;
	float		ballScoreFactor;
	score		ballScore;
	int			seedSet;
	char*		seed;
	} exscoreset;

// quantum dna codes--
//	A qcode describes a mapping from a quantum sequence character to its
//	probability vector in {A, C, G, T}

typedef struct qcode
	{
	char		dna[4];			// (usually "ACGT")
	double		p[256][4];		// p[sym][*] = dna probability vector for sym
								// sum of p[sym][*] is 1.0 if sym is meaningful
								// p[0][*] is usually meaningless
	} qcode;

// adapative score thresholds

typedef struct sthresh			// (could do as a C 'union', but why bother?)
	{
	char		t;				// type of threshold ('S', 'P', 'C')
	score		s;				// threshold as a score
	double		p;				// threshold as a percentage of bases in target
	u32			c;				// threshold as a count of bases in target
	} sthresh;

//----------
//
// globally available data in dna_utilities.c
//
//----------

#ifndef dna_utilities_owner
extern const s8    nuc_to_bits[256];
extern const s8    upper_nuc_to_bits[256];
extern const u8*   bits_to_nuc;
extern const u8*   bit_to_pur_pyr;
extern const u8*   bits_to_pur_pyr;
extern const u8    nuc_to_complement[256];
extern const u8    bits_to_complement[4];
extern       score HOXD70[4][4];
extern const score HOXD70_open;
extern const score HOXD70_extend;
extern const score HOXD70_X;
extern const score HOXD70_fill;
extern       score unitScores[4][4];
extern const double unitScores_open;
extern const double unitScores_extend;
extern const double unitScores_X;
extern const double unitScores_fill;
extern const double unitScores_thresh;
#endif

//----------
//
// prototypes for routines in dna_utilities.c
//
//----------

// macros to replace ctype.h (the standard C routines isupper, toupper, etc.
// behave strangely on 0x80..0xFF in some implementations)
// BEWARE: these all will fail if ch has side effects (such as ch = *s++)  

#define dna_isupper(ch)  (((ch)>='A')&&((ch)<='Z'))
#define dna_islower(ch)  (((ch)>='a')&&((ch)<='z'))
#define dna_isalpha(ch)  ((dna_isupper(c))||(dna_islower(c)))
#define dna_isprint(ch)  (((ch)>=0x20)&&((ch)<=0x7E))
#define dna_isxdigit(ch) ((((ch)>='0')&&((ch)<='9'))||(((ch)>='A')&&((ch)<='F'))||(((ch)>='a')&&((ch)<='f')))

#define dna_toupper(ch)  (dna_islower(ch)?((ch)-'a'+'A'):(ch))
#define dna_tolower(ch)  (dna_isupper(ch)?((ch)-'A'+'a'):(ch))
#define dna_toprint(ch)  (dna_isprint(ch)?(ch):'*')

// prototypes for real routines

scoreset*   new_dna_score_set         (score template[4][4],
                                       score xScore, score fillScore,
                                       score gapOpen, score gapExtend);
void        free_score_set            (char* id, scoreset* ss);
scoreset*   copy_score_set            (scoreset* ss);
scoreset*   masked_score_set          (scoreset* ss);
exscoreset* read_score_set_by_name    (char* name);
exscoreset* read_score_set            (FILE* f, char* name);
void        ambiguate_n               (scoreset* ss,
                                       score nVsN, score nVsNonN);
void        ambiguate_iupac           (scoreset* ss,
                                       score nVsN, score nVsNonN);
void        write_score_set_by_name   (char* name, scoreset* ss,
                                       int withGapScores);
void        write_score_set           (FILE* f, char* name, scoreset* ss,
                                       int withGapScores);
void        write_score_set_as_ints   (FILE* f, char* name, scoreset* ss,
                                       int withGapScores);
void        dump_dna_score_set        (FILE* f, scoreset* ss);
score       string_to_score           (const char* s);
void        scale_score_set           (scoreset* ss, double scale);
int         round_score               (double v);
void        print_score_matrix        (FILE* f, scoreset* ss, int withExtras);
void        print_score_matrix_lf     (FILE* f, scoreset* ss, int withExtras,
                                       char lineFeedCh);
void        print_score_matrix_prefix (FILE* f, scoreset* ss, int withExtras,
                                       char* prefix);
void        dump_score_set            (FILE* f, scoreset* ss,
                                       u8* rowChars, u8* colChars);
void        dump_lower_score_set      (FILE* f, scoreset* ss);
void        dump_full_score_set       (FILE* f, scoreset* ss);
void        print_quantum_word        (FILE* f, qcode* coding, u8* q,
                                       u32 wordLen);
void        print_quantum_dna_match   (FILE* f, qcode* coding, u8* q, u8* d,
                                       u32 wordLen);
char*       quantum_word_string       (u8* q, u32 wordLen, int symWidth); 
score       max_in_score_matrix       (scoreset* ss);
score       min_in_score_matrix       (scoreset* ss);
void        resolve_score_thresh      (sthresh* threshold, u32 denom);
sthresh     string_to_score_thresh    (const char* s);
char*       score_thresh_to_string    (const sthresh* threshold);
double      blastz_score_to_ncbi_bits (score bzScore);
double      blastz_score_to_ncbi_expectation (score bzScore);
qcode*      new_quantum_code          (void);
qcode*      read_quantum_code_by_name (char* name);
qcode*      read_quantum_code         (FILE* f, char* name);
int         print_dna_similarities    (FILE* f, const char* s1, const char* s2,
                                       int n);
char*       bits_to_nuc_string        (u64 word, int numChars);
double      entropy                   (u8* s, u8* t, int len);
double      entropy_lower_ok          (u8* s, u8* t, int len);
u64         rev_comp_by_pairs         (u64 word, int length);
u64         rev_comp_by_bits          (u64 word, int length);
char*       char_to_description       (char ch);

#undef global
#endif // dna_utilities_H
