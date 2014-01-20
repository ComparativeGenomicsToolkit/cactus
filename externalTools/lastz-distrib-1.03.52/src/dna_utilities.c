//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: dna_utilities.c
//
//----------
//
// dna_utilities--
//	Utility functions relating to DNA.
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <ctype.h>				// standard C upper/lower stuff
#include <math.h>				// standard C math stuff
#include "build_options.h"		// build options

#define  dna_utilities_owner	// (make this the owner of its globals)
#include "dna_utilities.h"		// interface to this module

// debugging defines

//#define bottleneckBiasOK		// if this is defined, the mapping from quantum
								// .. symbols to best scoring bottleneck char
								// .. will be biased toward earlier chars in the
								// .. alphabet

//----------
//
// globally available data
//
//----------

// nucleotide encoding--
//	nuc_to_bits maps an ascii character to a 2 bit nucleotide code;  the 2-bit
//	coding is designed so that the following are true:
//    nuc_to_bits['A'] xor nuc_to_bits['G'] is 2
//    nuc_to_bits['C'] xor nuc_to_bits['T'] is 2
//	Do not change this code without maintaining that relationship. 

#define __ -1
#define A_ 0
#define C_ 1
#define G_ 2
#define T_ 3

const s8 nuc_to_bits[256] =
	{
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,A_,__,C_,__,__,__,G_,__,__,__,__,__,__,__,__, // 4x (upper case)
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,A_,__,C_,__,__,__,G_,__,__,__,__,__,__,__,__, // 6x (lower case)
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__	 // Fx
	};

const s8 upper_nuc_to_bits[256] =
	{
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 0x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 1x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 2x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 3x (numbers)
	__,A_,__,C_,__,__,__,G_,__,__,__,__,__,__,__,__, // 4x (upper case)
	__,__,__,__,T_,__,__,__,__,__,__,__,__,__,__,__, // 5x (upper case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 6x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 7x (lower case)
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 8x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // 9x
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ax
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Bx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Cx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Dx
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__, // Ex
	__,__,__,__,__,__,__,__,__,__,__,__,__,__,__,__	 // Fx
	};

const u8* bits_to_nuc     = (u8*) "ACGT";
const u8* bit_to_pur_pyr  = (u8*) "RY";		// purine (AG) or pyramidine (CT)
const u8* bits_to_pur_pyr = (u8*) "RYRY";	// purine (AG) or pyramidine (CT)


const u8 nuc_to_complement[256] = // assumes upper/lower iupac code
	{
	0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F, // 0x
	0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F, // 1x
	0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x28,0x29,0x2A,0x2B,0x2C,0x2D,0x2E,0x2F, // 2x
	0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,0x3A,0x3B,0x3C,0x3D,0x3E,0x3F, // 3x (numbers)
	0x40,'T', 'V', 'G', 'H', 0x45,0x46,'C', 'D', 0x49,0x4A,'M', 0x4C,'K', 'N' ,0x4F, // 4x (upper case)
	0x50,0x51,'Y', 'S', 'A', 0x55,'B', 'W', 0x58,'R', 0x5A,0x5B,0x5C,0x5D,0x5E,0x5F, // 5x (upper case)
	0x60,'t', 'v', 'g', 'h', 0x65,0x66,'c', 'd', 0x69,0x6a,'m', 0x6c,'k', 'n' ,0x6f, // 6x (lower case)
	0x70,0x71,'y', 's', 'a', 0x75,'b', 'w', 0x78,'r', 0x7a,0x7b,0x7c,0x7d,0x7e,0x7f, // 7x (lower case)
	0x80,0x81,0x82,0x83,0x84,0x85,0x86,0x87,0x88,0x89,0x8A,0x8B,0x8C,0x8D,0x8E,0x8F, // 8x
	0x90,0x91,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9A,0x9B,0x9C,0x9D,0x9E,0x9F, // 9x
	0xA0,0xA1,0xA2,0xA3,0xA4,0xA5,0xA6,0xA7,0xA8,0xA9,0xAA,0xAB,0xAC,0xAD,0xAE,0xAF, // Ax
	0xB0,0xB1,0xB2,0xB3,0xB4,0xB5,0xB6,0xB7,0xB8,0xB9,0xBA,0xBB,0xBC,0xBD,0xBE,0xBF, // Bx
	0xC0,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,0xC8,0xC9,0xCA,0xCB,0xCC,0xCD,0xCE,0xCF, // Cx
	0xD0,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,0xD7,0xD8,0xD9,0xDA,0xDB,0xDC,0xDD,0xDE,0xDF, // Dx
	0xE0,0xE1,0xE2,0xE3,0xE4,0xE5,0xE6,0xE7,0xE8,0xE9,0xEA,0xEB,0xEC,0xED,0xEE,0xEF, // Ex
	0xF0,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,0xF8,0xF9,0xFA,0xFB,0xFC,0xFD,0xFE,0xFF  // Fx
	};

                                // A   C   G   T
const u8 bits_to_complement[4] = { T_, G_, C_, A_ };
#define A_ 0
#define C_ 1
#define G_ 2
#define T_ 3

#undef __
#undef A_
#undef C_
#undef G_
#undef T_
#undef ___

// default substitution scores

score HOXD70[4][4] =
	{//         A     C     G     T
	/* A */	{  91, -114,  -31, -123 },
	/* C */	{-114,  100, -125,  -31 },
	/* G */	{ -31, -125,  100, -114 },
	/* T */	{-123,  -31, -114,   91 },
	};

const score HOXD70_open   = 400;
const score HOXD70_extend = 30;
const score HOXD70_X      = -1000;
const score HOXD70_fill   = -100;

score unitScores[4][4] =
	{//       A   C   G   T
	/* A */	{ 1, -1, -1, -1 },
	/* C */	{-1,  1, -1, -1 },
	/* G */	{-1, -1,  1, -1 },
	/* T */	{-1, -1, -1,  1 },
	};

const double unitScores_open   =   3.25;		// 400/123
const double unitScores_extend =   0.24375;		//  30/123
const double unitScores_X      = -10.0;
const double unitScores_fill   =  -1.0;
const double unitScores_thresh =  30.0;

//----------
//
// prototypes for private functions
//
//----------

// macro to test membership in a string

#define in_string(ch,str)     (((ch)!=0) && (strchr(((char*)(str)),(ch))!=NULL))
#define not_in_string(ch,str) (((ch)==0) || (strchr(((char*)(str)),(ch))==NULL))

// real functions

static exscoreset* create_extended_score_set (void);
static scoreset*   create_score_set          (void);
static char*       quantum_visual            (int ch);

//----------
//
// new_dna_score_set--
//	Create a new score set.
//
//----------
//
// Arguments:
//	score	template[4][4]:		The template containing the scores, with rows
//								.. and columns corresponding to bits_to_nuc[].
//								.. A row corresponds to a character in sequence
//								.. 1 and a column corresponds to a character in
//								.. sequence 2.  This can be NULL if the caller
//								.. doesn't care about initial scores.
//	score	badScore:			The score to use for row and column 'X'.
//	score	fillScore:			The score to use for all other rows and
//								.. columns.
//	score	gapOpen, gapExtend:	Gap scoring parameters.
//
//----------
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free() should be used for this purpose.
//
//----------
//
// Notes:
//	(1)	In the resulting scoring matrix, upper and lower case characters are
//		are considered identical, so entries for lower case are copied from
//		upper case.
//
//----------

scoreset* new_dna_score_set
   (score		template[4][4],
	score		badScore,
	score		fillScore,
	score		gapOpen,
	score		gapExtend)
	{
	scoreset*	ss;
	u8*			s, *d;
	int			r, c, rowCh, colCh, rowLower, colLower;
	int			len;

	//////////
	// allocate the score set
	//////////

	// allocate

	ss = create_score_set ();

	// set character set

	ustrcpy (ss->rowChars, bits_to_nuc);
	s = ss->rowChars;  len = ustrlen(s);  d = s + len;
	for ( ; len>0 ; s++,len--)
		{ *(d++) = dna_tolower(*s);  *d = 0; }

	ustrcpy (ss->colChars, bits_to_nuc);
	s = ss->colChars;  len = ustrlen(s);  d = s + len;
	for ( ; len>0 ; s++,len--)
		{ *(d++) = dna_tolower(*s);  *d = 0; }

	ss->badRow = ss->badCol = 'X';

	ss->rowsAreDna = true;
	ss->colsAreDna = true;

	// set gap scoring parameters

	ss->gapOpen   = gapOpen;
	ss->gapExtend = gapExtend;

	//////////
	// fill the array with a filler score and make sure scores for row and
	// column zero are very very bad
	//////////

	// fill row 0

	for (c=0 ; c<256 ; c++)
		ss->sub[0][c] = veryBadScore;

	// fill row 1

	ss->sub[1][0] = veryBadScore;
	for (c=1 ; c<256 ; c++)
		ss->sub[1][c] = fillScore;

	// copy row 1 to the remaining rows

	for (r=2 ; r<256 ; r++)
		memcpy (/*to*/ ss->sub[r], /*from*/ ss->sub[1],
		        /*how much*/ sizeof(scorerow));

	//////////
	// set up the remaining rows
	//////////

	// fill in X scores

	for (c=0 ; c<256 ; c++)
		{
		  ss->sub['X'][ c ]
		= ss->sub['x'][ c ]
		= ss->sub[ c ]['X']
		= ss->sub[ c ]['x'] = badScore;
		}

	// copy scores from the template

	if (template != NULL)
		{
		for (r=0 ; r<4 ; r++)
				for (c=0 ; c<4 ; c++)
			{
			rowCh    = bits_to_nuc[r];
			colCh    = bits_to_nuc[c];
			rowLower = dna_tolower(rowCh);
			colLower = dna_tolower(colCh);

			ss->sub[rowCh   ][colCh   ] = template[r][c];
			ss->sub[rowCh   ][colLower] = template[r][c];
			ss->sub[rowLower][colCh   ] = template[r][c];
			ss->sub[rowLower][colLower] = template[r][c];
			}
		}

	return ss;
	}

//----------
//
// create_score_set, create_extended_score_set--
//	Create a new score set but don't initialize it.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free_score_set() should be used for
//	this purpose.
//
// Notes:
//	- The internal rows of the scoring matrix are *not* initialized.
//
//----------

static scoreset* create_either_score_set (int isExtended);


static exscoreset* create_extended_score_set ()
	{ return (exscoreset*) create_either_score_set (true); }
static scoreset* create_score_set ()
	{ return create_either_score_set (false); }


static scoreset* create_either_score_set
   (int isExtended)
	{
	scoreset*	ss;
	u32			bytesNeeded;
	int			r;

	// allocate

	if (isExtended) bytesNeeded = round_up_8 (sizeof(exscoreset));
	           else bytesNeeded = round_up_8 (sizeof(scoreset));
	ss = (scoreset*) zalloc_or_die ("create_score_set", bytesNeeded);

	// initialize

	ss->rowChars[0] = 0;
	ss->colChars[0] = 0;
	ss->badRow = ss->badCol = -1;

	ss->gapOpenSet   = false;
	ss->gapExtendSet = false;

	ss->bottleneck[0] = 0;
	for (r=0 ; r<256 ; r++) ss->qToBest[r].len = -1;

	ss->qToComplement = NULL;

	// return the score set

	return ss;
	}

//----------
//
// free_score_set--
//	Deallocate a score set, along with any associated memory.
//
//----------
//
// Arguments:
//	char* 		id:	an identifying string to be used when trackMemoryUsage is
//					.. turned on;  this can be NULL.
//	scoreset*	ss:	The score set to dispose of.
//
// Returns:
//	(nothing)
//
//----------

void free_score_set
   (char*		id,
	scoreset*	ss)
	{
	char*		message = " (qToComplement)";
	char		temp[200];

	if (ss == NULL) return;
	if (id == NULL) id = "free_score_set";

	if (strlen(id) + strlen(message) + 1 > sizeof(temp))
		free_if_valid (id, ss->qToComplement);
	else
		{
		strcpy (temp, id);
		strcpy (temp + strlen(id), message);
		free_if_valid (temp, ss->qToComplement);
		}

	free_if_valid (id, ss);
	}

//----------
//
// copy_score_set--
//	Create a copy of a score set.
//
//----------
//
// Arguments:
//	scoreset*	ss:	The score set to copy.
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free() should be used for this purpose.
//
//----------

scoreset* copy_score_set
   (scoreset*	ss)
	{
	scoreset*	ssNew;
	int			r, c;

	// alloacte new score set

	ssNew = create_score_set ();

	if (ss->qToComplement != NULL)
		ssNew->qToComplement = (u8*) malloc_or_die ("copy_score_set (qToComplement)", 256);

	// copy gap scoring parameters

	ssNew->gapOpen      = ss->gapOpen;
	ssNew->gapOpenSet   = ss->gapOpenSet;

	ssNew->gapExtend    = ss->gapExtend;
	ssNew->gapExtendSet = ss->gapExtendSet;

	// copy character set

	ssNew->badRow = ss->badRow;
	ssNew->badCol = ss->badCol;

	ustrcpy (ssNew->rowChars, ss->rowChars);
	ustrcpy (ssNew->colChars, ss->colChars);

	ssNew->rowsAreDna = ss->rowsAreDna;
	ssNew->colsAreDna = ss->colsAreDna;

	ustrcpy (ssNew->bottleneck, ss->bottleneck);
	for (r=0 ; r<256 ; r++) ssNew->qToBest[r] = ss->qToBest[r];

	if (ss->qToComplement != NULL)
		{
		for (c=0 ; c<256 ; c++)
			ssNew->qToComplement[c] = ss->qToComplement[c];
		}

	// copy substitution scores

	memcpy (/*to*/ ssNew->sub, /*from*/ ss->sub, /*how much*/ sizeof(ss->sub));

	return ssNew;
	}

//----------
//
// masked_score_set--
//	Create a copy of a score set, with all lower case entries given 'bad'
//	scores.
//
//----------
//
// Arguments:
//	scoreset*	ss:	The score set to copy.  Columns can be either DNA or
//					.. quantum DNA, but rows must be DNA.
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free() should be used for this purpose.
//
//----------

scoreset* masked_score_set
   (scoreset*	ss)
	{
	score		badScore;
	scoreset*	ssNew;
	int			r, c, goodRow;
	u8*			rr, *cc, *d;
	int			nIsARow, nIsACol;

	// $$$ the tests for rowsAreDna and colsAreDna should be replaced by
	// $$$ .. a new field indicating whether rows/cols are maskable, and
	// $$$ .. which characters survive masking

	// copy the score set and reduce copy's good characters to upper case 

	ssNew = copy_score_set (ss);

	if (ss->rowsAreDna)
		{
		d = ssNew->rowChars;
		for (rr=ss->rowChars ; *rr!=0 ; rr++)
			{ if (dna_isupper (*rr)) { *(d++) = *rr;  *d = 0; } }
		}

	if (ss->colsAreDna)
		{
		d = ssNew->colChars;
		for (cc=ss->colChars ; *cc!=0 ; cc++)
			{ if (dna_isupper (*cc)) { *(d++) = *cc;  *d = 0; } }
		}

	// mask the copy;  fill each lower case row or column with a bad score

	goodRow = ss->rowChars[0];
	badScore = ss->sub[goodRow][ss->badCol];

	if (ss->rowsAreDna)
		{
		nIsARow = (ustrchr (ssNew->rowChars, 'N') != NULL);

		for (rr=ss->rowChars ; *rr!=0 ; rr++)
				if (!dna_isupper (*rr))
			{ for (c=1 ; c<256 ; c++) ssNew->sub[*rr][c] = badScore; }
		if (!nIsARow) for (c=1 ; c<256 ; c++) ssNew->sub['N'][c] = badScore;
		for (c=1 ; c<256 ; c++) ssNew->sub['n'][c] = badScore;
		for (c=1 ; c<256 ; c++) ssNew->sub['X'][c] = badScore;
		}

	if (ss->colsAreDna)
		{
		nIsACol = (ustrchr (ssNew->colChars, 'N') != NULL);

		for (cc=ss->colChars ; *cc!=0 ; cc++)
				if (!dna_isupper (*cc))
			{ for (r=1 ; r<256 ; r++) ssNew->sub[r][*cc] = badScore; }
		if (!nIsACol) for (r=1 ; r<256 ; r++) ssNew->sub[r]['N'] = badScore;
		for (r=1 ; r<256 ; r++) ssNew->sub[r]['n'] = badScore;
		for (r=1 ; r<256 ; r++) ssNew->sub[r]['X'] = badScore;
		}

	return ssNew;
	}

//----------
//
// read_score_set_by_name, read_score_set--
//	Read a new score set from a file (see format description below).
//
//----------
//
// Arguments:
//	FILE*	f:				(read_score_set only) The file that scoring data is
//							.. to be read from.  This should already be open for
//							.. text read.
//	char*	name:			The name of the file that scoring data is to be read
//							.. from.  For read_score_set this is only used for
//							.. reporting problems to the user (and may be NULL).
//
// Returns:
//	A pointer to the newly allocated score set, which the caller will have to
//	dispose of eventually.  The routine free() should be used for this purpose.
//
//----------
//
// Score Set File Format
// =====================
//
// Here's an example:
//
//		# This matches the default scoring set for blastz
//
//		bad_score          = X:-1000  # used for sub['X'][*] and sub[*]['X']
//		fill_score         = -100     # used when sub[*][*] not defined
//		gap_open_penalty   =   30
//		gap_extend_penalty =  400
//
//			 A     C     G     T
//		A   91  -114   -31  -123
//		C -114   100  -125   -31
//		G  -31  -125   100  -114
//		T -123   -31  -114    91
//
// The score set consists of a substitution matrix and other settings.  The
// other settings come first.  Any line may contain a comment, # is the comment
// character.
//
// Labels can either be single characters, or two-digit hexadecimal character
// values (the value 00 is not allowed).  Rows and columns of the matrix need
// not have the same labels or range, so, for example, a matrix might describe
// scoring between the 15-letter ambiguity code and the 4-letter DNA code.
//
// Row labels are optional, and if absent it is presumed that they are the same
// as the column labels (and in the same order).  This allows us to read blastz
// score files.
//
// For quantum alphabets, column labels may indicate reverse complement pairing.
// This is done by specifying each column as, e.g. A~T.  If any labels indicate
// complement they all must.
//
// The bad_score setting is optional, and the X shown above is the character
// for which all scores will be marked as bad.  A separate character can be
// specified for rows and columns by using <row>:<col>:<score>.  Both <row>
// and <col> are optional.
//
// The other settings, fill_score, gap_open_penalty and gap_extend_penalty, are
// also optional, and defaults compatible with blastz are used.
//
// Rows correspond to characters in sequence 1 and columns correspond to
// characters in sequence 2.
//
// Score values can be floating-point if the library is built with an
// appropriate scoreType.
//
//----------

static int parse_char_code         (char** s, int* comp, char terminator);
static int parse_char_code_zero_ok (char** s, int* comp, char terminator);
static int parse_char_code_common  (char** s, int* comp, char terminator, int zeroOk);
static int is_dna_alphabet         (u8* alphabet);
static u8  two_char_as_hex         (u8 ch1, u8 ch2);
static int parse_bottleneck        (char* s, u8 bottleneck[5]);


exscoreset* read_score_set_by_name
   (char*		name)
	{
	FILE*		f;
	exscoreset*	xss;

	if (name == NULL)
		suicide ("can't open NULL file in read_score_set_by_name()");

	f = fopen_or_die (name, "rt");
	xss = read_score_set (f, name);
	fclose_if_valid (f);

	return xss;
	}


exscoreset* read_score_set
   (FILE*		f,
	char*		_name)
	{
	// $$$ why is the line buffer declared as static?
	static char	line[256*25+1];	// (must hold 256 fields, up to 25 chars each)
	char*		name = _name;
	u8			rowChars[256], colChars[256], colComps[256];
	u8*			scanRowChars, *scanColChars;
	int			lineNum, len, missingEol;
	int			badRow, badCol, r, c, compC;
	u8			bottleneck[5];
	score		badScore, fillScore, gapOpen, gapExtend,
	            hspThreshold, gappedThreshold, xDrop, yDrop, ballScore;
	float		ballScoreFactor;
	u32			step;
	char*		seed;
	scorerow	fillRowData;
	exscoreset*	xss;
	char*		s, *prevS, *waffle, *end;
	int			numRows, numCols, numFields, fieldCount, ix, iy;
	int			haveBottleneck, haveFillScore, haveGapOpen, haveGapExtend,
	            haveHspThreshold, haveGappedThreshold,
				haveXDrop, haveYDrop, haveStep, haveBallScore, haveSeed;
	char*		valString, *scan, *colon;
	int			valLength, finalField, haveComps;
	u8*			src, *dst, ch;

	if (name == NULL)
		name = "(unnamed file)";

	//////////
	// read assignments
	//////////

	badScore  = -1000;	// (deafults from blastz)
	fillScore =  -100;
	gapOpen   = HOXD70_open;
	gapExtend = HOXD70_extend;

	badCol = badRow = -1;
	haveBottleneck = haveFillScore = haveGapOpen = haveGapExtend = false;
	haveHspThreshold = haveGappedThreshold = false;
	haveXDrop = haveYDrop = haveStep = haveBallScore = haveSeed = false;
	hspThreshold = gappedThreshold = xDrop = yDrop = ballScore = 0;
	ballScoreFactor = -1;
	step = 0;
	seed = NULL;

	lineNum    = 0;
	missingEol = false;

	while (fgets (line, sizeof(line), f) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (missingEol)
			suicidef ("line is too long (%s: line %d)", name, lineNum-1);

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		// trim blanks, end of line, and comments, and ignore blank lines

		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// if it doesn't contain an assignment, more on to the next phase

		valString = strchr (line, '=');
		if (valString == NULL) break;

		// parse the assignment

		*(valString++) = 0;
		trim_string (line);
		trim_string (valString);

		if ((!haveBottleneck)
		 && (strcmp (line, "bottleneck") == 0))
			{
			if (!parse_bottleneck (valString, bottleneck))
				suicidef ("invalid bottleneck alphabet (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			haveBottleneck = true;
			}
		else if ((badCol == -1)
		      && ((strcmp (line, "bad") == 0)
		       || (strcmp (line, "bad_score") == 0)))
			{
			// parse [<row>:[<col>:]]<score>, no whitespace allowed

			scan = valString;

			colon = strchr (scan, ':');
			if (colon != NULL)
				{
				badCol = badRow = parse_char_code_zero_ok (&scan, NULL, ':');
				scan   = colon+1;
				if (badCol < 0)
					suicidef ("invalid bad_score character code (%s: line %d) %s=%s",
					          name, lineNum, line, valString);
				}

			colon = strchr (scan, ':');
			if (colon != NULL)
				{
				badRow = parse_char_code_zero_ok (&scan, NULL, ':');
				scan   = colon+1;
				if (badRow < 0)
					suicidef ("invalid bad_score character code (%s: line %d) %s=%s",
					          name, lineNum, line, valString);
				}

			badScore = string_to_score (scan);
			}
		else if ((!haveFillScore)
		      && ((strcmp (line, "fill") == 0)
		       || (strcmp (line, "fill_score") == 0)))
			{
			fillScore     = string_to_score (valString);
			haveFillScore = true;
			}
		else if ((!haveGapOpen)
		      && ((strcmp (line, "O") == 0)
		       || (strcmp (line, "open") == 0)
		       || (strcmp (line, "gap_open") == 0)
		       || (strcmp (line, "gap_open_penalty") == 0)))
			{
			gapOpen     = string_to_score (valString);
			haveGapOpen = true;
			}
		else if ((!haveGapExtend)
		      && ((strcmp (line, "E") == 0)
		       || (strcmp (line, "extend") == 0)
		       || (strcmp (line, "gap_extend") == 0)
		       || (strcmp (line, "gap_extend_penalty") == 0)))
			{
			gapExtend     = string_to_score (valString);
			haveGapExtend = true;
			}
		else if ((!haveHspThreshold)
		      && ((strcmp (line, "K") == 0)
		       || (strcmp (line, "hsp_thresh") == 0)
		       || (strcmp (line, "hsp_threshold") == 0)))
			{
			hspThreshold     = string_to_score (valString);
			haveHspThreshold = true;
			}
		else if ((!haveGappedThreshold)
		      && ((strcmp (line, "L") == 0)
		       || (strcmp (line, "gapped_thresh") == 0)
		       || (strcmp (line, "gapped_threshold") == 0)))
			{
			gappedThreshold     = string_to_score (valString);
			haveGappedThreshold = true;
			}
		else if ((!haveXDrop)
		      && ((strcmp (line, "X") == 0)
		       || (strcmp (line, "x_drop") == 0)))
			{
			xDrop     = string_to_score (valString);
			haveXDrop = true;
			if (xDrop <= 0)
				suicidef ("invalid x-drop threshold (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			}
		else if ((!haveYDrop)
		      && ((strcmp (line, "Y") == 0)
		       || (strcmp (line, "y_drop") == 0)))
			{
			yDrop     = string_to_score (valString);
			haveYDrop = true;
			if (yDrop <= 0)
				suicidef ("invalid y-drop threshold (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			}
		else if ((!haveStep)
		      && ((strcmp (line, "Z") == 0)
		       || (strcmp (line, "step") == 0)))
			{
			step = string_to_int (valString);
			haveStep = true;
			if (step <= 0)
				suicidef ("invalid step (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			}
		else if ((!haveBallScore)
		      && (strcmp (line, "ball") == 0))
			{
			valLength = strlen(valString);
			if ((valLength > 0) && (valString[valLength-1] == '%'))
				{
				ballScoreFactor = pct_string_to_double (valString);
				haveBallScore   = true;
				if ((ballScoreFactor <= 0) || (ballScoreFactor > 1))
					suicidef ("invalid quantum ball score (%s: line %d) %s=%s",
					          name, lineNum, line, valString);
				}
			else
				{
				ballScore     = string_to_score (valString);
				haveBallScore = true;
				if (ballScore <= 0)
					suicidef ("invalid quantum ball score (%s: line %d) %s=%s",
					          name, lineNum, line, valString);
				}
			}
		else if ((!haveSeed)
		      && (strcmp (line, "T") == 0))
			{
			if (strcmp (valString, "1") == 0)
				seed = copy_string ("T=1");
			else if (strcmp (valString, "2") == 0)
				seed = copy_string ("T=2");
			else if (strcmp (valString, "3") == 0)
				seed = copy_string ("T=3");
			else if (strcmp (valString, "4") == 0)
				seed = copy_string ("T=4");
			else
				suicidef ("invalid seed (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			haveSeed = true;
			}
		else if ((!haveSeed)
		      && (strcmp (line, "seed") == 0))
			{
			if ((strcmp (valString, "12of19,transition") == 0)
			 || (strcmp (valString, "12_of_19,transition") == 0))
				seed = copy_string ("T=1");
			else if ((strcmp (valString, "12of19,notransition") == 0)
			      || (strcmp (valString, "12_of_19,no_transition") == 0))
				seed = copy_string ("T=2");
			else if ((strcmp (valString, "14of22,transition") == 0)
			      || (strcmp (valString, "14_of_22,transition") == 0))
				seed = copy_string ("T=3");
			else if ((strcmp (valString, "14of22,notransition") == 0)
			      || (strcmp (valString, "14_of_22,no_transition") == 0))
				seed = copy_string ("T=4");
			else
				suicidef ("invalid seed (%s: line %d) %s=%s",
				          name, lineNum, line, valString);
			haveSeed = true;
			}
		else
			suicidef ("invalid name in assignment (%s: line %d) %s=%s",
			          name, lineNum, line, valString);
		}

	//////////
	// read column characters
	//////////

	// current line caused us to exit assignment stage, so it must contain
	// the column headers

	for (c=0 ; c<256 ; c++)
		colComps[c] = 0;

	haveComps = -1;

	colChars[0] = 0;
	for (s=line,scanColChars=colChars ; *s!=0 ; )
		{
		prevS = s;
		c = parse_char_code (&s, &compC, ' ');

		if (c <= 0)
			suicidef ("invalid character code in %s:line %d at \"%s\"",
			          name, lineNum, s);

		if (compC < 0)
			suicidef ("invalid complement in %s:line %d at \"%s\"",
			          name, lineNum, s);

		if (in_string (c, colChars))
			suicidef ("duplicate character code in %s:line %d at \"%s\"",
			          name, lineNum, s);

		if (haveComps == -1)
			haveComps = (compC != 0);
		else if (haveComps)
			{
			if (compC == 0)
				suicidef ("missing complement in %s:line %d at \"%s\"",
				          name, lineNum, prevS);
			}
		else // if (!haveComps)
			{
			if (compC != 0)
				suicidef ("missing complement(s) in %s:line %d before \"%s\"",
				          name, lineNum, prevS);
			}

		*(scanColChars++) = c;  *scanColChars = 0;
		colComps[c] = compC;
		}

	numCols = scanColChars - colChars;

	if ((badCol >= 0) && (in_string (badCol, colChars)))
		suicidef ("character code for bad_score can't also be a matrix column\n"
		          "(%s: line %d)",
		          name, lineNum);

	if (numCols == 0)
		suicidef ("matrix has no column headers (%s: line %d)",
		          name, lineNum);

	// validate complements

	if (haveComps)
		{
		for (ix=0 ; ix<numCols ; ix++)
			{
			c     = colChars[ix];
			compC = colComps[c];

			if (!in_string (compC, colChars))
				suicidef ("complement (%s~%s) not in column alphabet in %s:line %d",
				          quantum_visual(c), quantum_visual(compC), name, lineNum);

			if (colComps[compC] != c)
				suicidef ("complement (%s~%s~%s) is not symmetric in %s:line %d",
				          quantum_visual(c), quantum_visual(compC),
				          quantum_visual(colComps[compC]), name, lineNum);
			}
		}

	//////////
	// create the scoring matrix
	//////////

	xss = create_extended_score_set ();

	// fill it with a filler score

	for (c=0 ; c<256 ; c++)
		fillRowData[c] = fillScore;

	for (r=0 ; r<256 ; r++)
		memcpy (/*to*/ xss->ss.sub[r], /*from*/ fillRowData,
		        /*how much*/ sizeof(scorerow));

	// disable the extra parameters

	xss->hspThresholdSet    = false;
	xss->gappedThresholdSet = false;
	xss->xDropSet           = false;
	xss->yDropSet           = false;
	xss->stepSet            = false;
	xss->ballScoreSet       = false;
	xss->ballScoreFactor    = -1;
	xss->seedSet            = false;

	xss->seed               = NULL;

	//////////
	// read the scoring matrix, filling in over the scores we've already filled
	//////////

	scanRowChars = rowChars;
	*scanRowChars = 0;

	numFields = -1;
	iy = 0;

	while (fgets (line, sizeof(line), f) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (missingEol)
			suicidef ("line is too long (%s: line %d)", name, lineNum-1);

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		// trim blanks, end of line, and comments, and ignore blank lines

		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// count the number of fields

		fieldCount = 0;
		for (s=line ; *s!=0 ; )
			{
			s = skip_darkspace  (s);
			s = skip_whitespace (s);
			fieldCount++;
			}

		if (numFields < 0)
			{
			numFields = fieldCount;
			if ((numFields != numCols)
			 && (numFields != numCols+1))
				suicidef ("wrong number of score columns (%s: line %d)",
				          name, lineNum);
			}
		else if (fieldCount != numFields)
			suicidef ("inconsistent number of score columns (%s: line %d)",
			          name, lineNum);

		// first field is character code for the row;  for blastz compatibility
		// we just assign the next DNA character

		s = line;

		if (numFields == numCols)
			{
			if (iy >= numCols)
				suicidef ("too many score rows (%s: line %d): \"%s\"",
				          name, lineNum, line);
			r = colChars[iy++];
			*(scanRowChars++) = r;  *scanRowChars = 0;
			}
		else
			{
			r = parse_char_code (&s, NULL, ' ');
			if (r <= 0)
				suicidef ("invalid row character code (%s: line %d) %s=%s",
				          name, lineNum, line, s);
			if (in_string (r, rowChars))
				suicidef ("duplicate row character code (%s: line %d): \"%s\"",
				          name, lineNum, line);

			*(scanRowChars++) = r;  *scanRowChars = 0;
			}

		// remaining fields are the rows in this column

		for (ix=0 ; ix<numCols ; ix++)
			{
			if (*s == 0) // (can't happen due to earlier test)
				suicidef ("not enough score columns (%s: line %d)",
				          name, lineNum);

			c = colChars[ix];

			end = skip_darkspace (s);
			finalField = (*end == 0);
			*end = 0;

			xss->ss.sub[r][c] = string_to_score (s);
			if (finalField) s = end;
			           else s = skip_whitespace (end+1);
			}
		}

	numRows = scanRowChars - rowChars;

	if (numFields < 0)
		suicidef ("scores file %s contains no score rows", name);
	if ((numFields == numCols) && (numRows != numCols))
		suicidef ("not enough score rows, line (%s: line %d): \"%s\"",
				  name, lineNum, line);

	if ((badRow >= 0) && (in_string (badRow, rowChars)))
		suicide ("character code for bad_score can't also be a matrix row");

	//////////
	// finish off the scoring matrix
	//////////

	ustrcpy (xss->ss.colChars, colChars);
	ustrcpy (xss->ss.rowChars, rowChars);

	xss->ss.gapOpen      = gapOpen;
	xss->ss.gapOpenSet   = haveGapOpen;
	xss->ss.gapExtend    = gapExtend;
	xss->ss.gapExtendSet = haveGapExtend;

	if ((xss->ss.gapOpenSet) && (xss->ss.gapOpen + xss->ss.gapExtend <= 0))
		suicidef (" (in %s) "scoreFmt " is not a valid gap open penalty with extension penalty " scoreFmt "\n"
		          "(open can be negative but the sum has to be postive)\n",
		          name, xss->ss.gapOpen, xss->ss.gapExtend);
	if ((xss->ss.gapExtendSet) && (xss->ss.gapExtend < 0))
		suicidef (scoreFmt " is not a valid gap extension penalty (in %s)\n",
		          xss->ss.gapExtend, name);

	// set any extra parameters we have

	if (haveHspThreshold)
		{ xss->hspThresholdSet    = true;  xss->hspThreshold    = hspThreshold; }
	if (haveGappedThreshold)
		{ xss->gappedThresholdSet = true;  xss->gappedThreshold = gappedThreshold; }
	if (haveXDrop)
		{ xss->xDropSet           = true;  xss->xDrop           = xDrop; }
	if (haveYDrop)
		{ xss->yDropSet           = true;  xss->yDrop           = yDrop; }
	if (haveStep)
		{ xss->stepSet            = true;  xss->step            = step; }
	if (haveBallScore)
		{ xss->ballScoreSet       = true;  xss->ballScore       = ballScore;
		                                   xss->ballScoreFactor = ballScoreFactor; }
	if (haveSeed)
		{ xss->seedSet            = true;  xss->seed            = seed; }

	// if the columns are DNA, make lower case columns equivalent to upper case

	xss->ss.colsAreDna = is_dna_alphabet (colChars);

	if (xss->ss.colsAreDna)
		{
		if (badCol < 0) badCol = 'X';
		for (ix=0 ; ix<numCols ; ix++)
			{
			c = colChars[ix];
			for (iy=0 ; iy<numRows ; iy++)
				{
				r = rowChars[iy];
				xss->ss.sub[r][c+'a'-'A'] = xss->ss.sub[r][c];
				}
			}

		src = xss->ss.colChars;
		dst = src + ustrlen(src);
		for ( ; *src!=0 ; src++)
			{
			ch = dna_tolower (*src);
			if (ustrchr (xss->ss.colChars, ch) == NULL)
				{ *(dst++) = ch;  *dst = 0; }
			}
		}

	// if the rows are DNA, make lower case rows equivalent to upper case

	xss->ss.rowsAreDna = is_dna_alphabet (rowChars);

	if (xss->ss.rowsAreDna)
		{
		if (badRow < 0) badRow = 'X';
		for (ix=0 ; ix<numRows ; ix++)
			{
			r = rowChars[ix];
			memcpy (/*to*/ xss->ss.sub[r+'a'-'A'], /*from*/ xss->ss.sub[r],
			        /*how much*/ sizeof(scorerow));
			}

		src = xss->ss.rowChars;
		dst = src + ustrlen(src);
		for ( ; *src!=0 ; src++)
			{
			ch = dna_tolower (*src);
			if (ustrchr (xss->ss.rowChars, ch) == NULL)
				{ *(dst++) = ch;  *dst = 0; }
			}
		}

	// fill the bad row and column

	if (badCol == -1) badCol = 0;	// (if rows and/or cols were DNA, these
	if (badRow == -1) badRow = 0;	//  .. would already be set to 'X')

	xss->ss.badRow = badRow;
	xss->ss.badCol = badCol;

	for (c=0 ; c<256 ; c++) xss->ss.sub[badRow][c] = badScore;
	for (r=0 ; r<256 ; r++) xss->ss.sub[r][badCol] = badScore;

	// make sure scores for row and column zero are very very bad

	for (c=0 ; c<256 ; c++)
		xss->ss.sub[0][c] = xss->ss.sub[c][0] = veryBadScore;

	//////////
	// create complement-mapping table
	//////////

	if (haveComps)
		{
		xss->ss.qToComplement = (u8*) malloc_or_die ("read_score_set (qToComplement)", 256);
		for (c=0 ; c<256 ; c++)
			xss->ss.qToComplement[c] = colComps[c];
		}

	//////////
	// create bottleneck-related fields
	//////////

	// set neutered defaults

	xss->ss.bottleneck[0] = 0;
	for (r=0 ; r<256 ; r++) xss->ss.qToBest[r].len = -1;

	// if we don't have a quantum row alphabet, we can't have a bottleneck

	if ((haveBottleneck) && (xss->ss.rowsAreDna))
		suicidef ("invalid bottleneck alphabet (%s in %s), rows are DNA",
				  bottleneck, name);

	// if we don't have a quantum column alphabet, the bottleneck has to be DNA

	if ((haveBottleneck) && (xss->ss.colsAreDna)
	 && (ustrcmp (bottleneck, (u8*) "ACGT") != 0))
		suicidef ("invalid bottleneck alphabet (%s in %s), columns are DNA",
				  bottleneck, name);

	// if we have quantum rows and DNA columns but no bottleneck, assign
	// a DNA bottleneck

	if ((!haveBottleneck) && (!xss->ss.rowsAreDna) && (xss->ss.colsAreDna))
		{ ustrcpy (bottleneck, (u8*) "ACGT");  haveBottleneck = true; }

	// if we have quantum rows and quantum columns, we gotta have a bottleneck

	if ((!haveBottleneck) && (!xss->ss.rowsAreDna) && (!xss->ss.colsAreDna))
		suicidef ("missing bottleneck alphabet (in %s)", name);

	// ok, let's fill in the fields

	if (haveBottleneck)
		{
		u8*		rr;
		u8		bits;
		charvec	bestBits;
		score	bestScore, thisScore;

#ifdef bottleneckBiasOK
		bestBits.v[0]									// (placate compiler)
		   = bestBits.v[1]
		   = bestBits.v[2]
		   = bestBits.v[3]
		   = 0;
#endif // no bottleneckBiasOK

		// all bottleneck chars must be in column alphabet

		if ((ustrchr (xss->ss.colChars, bottleneck[0]) == NULL)
		 || (ustrchr (xss->ss.colChars, bottleneck[1]) == NULL)
		 || (ustrchr (xss->ss.colChars, bottleneck[2]) == NULL)
		 || (ustrchr (xss->ss.colChars, bottleneck[3]) == NULL))
			suicidef ("invalid bottleneck alphabet (%s in %s)"
			          ", not contained in column alphabet",
			          bottleneck, name);

		ustrcpy (xss->ss.bottleneck, bottleneck);

		// find 'closest' match for each row character

		for (rr=xss->ss.rowChars ; *rr!=0 ; rr++)
			{
			r = *rr;
			c = bottleneck[0];
			bestBits.len = 0;
			bestScore    = worstPossibleScore;
			for (bits=0 ; bits<4 ; bits++)
				{
				c         = bottleneck[bits];
				thisScore = xss->ss.sub[r][c];
				if (thisScore > bestScore)
					{ // (this character is 'closest' so far)
					bestBits.len  = 1;
					bestBits.v[0] = bits;
					bestScore     = thisScore;
					}
#ifndef bottleneckBiasOK
				else if (thisScore == bestScore)
					{ // (this character is tied for 'closest' so far)
					bestBits.v[(u8)(bestBits.len++)] = bits;
					}
#endif // not bottleneckBiasOK
				}
			if (bestBits.len == 0) // (can't happen, but compiler frets)
				bestBits.len = -1;
			xss->ss.qToBest[r] = bestBits;
			}

		if (dna_utilities_dbgShowQToBest)
			{
			for (rr=xss->ss.rowChars ; *rr!=0 ; rr++)
				{
				r = *rr;
				bestBits = xss->ss.qToBest[r];
				fprintf (stderr, "qToBest[%02X]:", r);
				if (bestBits.len == -1)
					fprintf (stderr, " (none)");
				for (ix=0 ; ix<bestBits.len ; ix++)
					{
					bits = bestBits.v[ix];
					c    = bottleneck[bits];
					fprintf (stderr, " %02X", c);
					}
				fprintf (stderr, "\n");
				}
			}
		}

	return xss;
	}


static int parse_char_code (char** s, int* comp, char terminator)
	{ return parse_char_code_common (s, comp, terminator, false); }

static int parse_char_code_zero_ok (char** s, int* comp, char terminator)
	{ return parse_char_code_common (s, comp, terminator, true); }

static int parse_char_code_common
   (char**	_s,
	int*	comp,
   	char	terminator,
   	int		zeroOk)
	{
	char*	s = *_s;
	int		cc, cc2;
	char	follower;

	// parse the (first) character code

	cc       = *(s++);
	follower = *s;

	if (isxdigit(follower))
		{
		s++;
		if (isxdigit(cc)) cc = two_char_as_hex(cc,follower);
		             else cc = -1;
		if ((!zeroOk) && (cc == 0)) cc = -1;
		}

	// parse the (second) character code

	cc2 = 0;

	if ((comp != NULL) && (*s != '~'))
		; // (cc2 = 0 is correct)
	else if ((comp != NULL) && (*s == '~'))
		{
		s++;
		cc2      = *(s++);
		follower = *s;

		if (isxdigit(follower))
			{
			s++;
			if (isxdigit(cc2)) cc2 = two_char_as_hex(cc2,follower);
			              else cc2 = -1;
			if ((!zeroOk) && (cc2 == 0)) cc2 = -1;
			}
		}

	// eat up trailing whitespace

	if (terminator == ' ')
		{
		if ((*s != 0) && (!isspace (*s))) cc = 0;
		                             else s = skip_whitespace (s);
		}
	else if (terminator != 0)
		{
		if (*s != terminator) cc = -1;
		                 else s++;
		}

	if ((cc >= 0) && (cc2 >= 0)) *_s = s;
	if (comp != NULL) *comp = cc2;
	return cc;
	}


static int is_dna_alphabet
   (u8*	alphabet)
	{
	int match = 0;

	if (ustrchr (alphabet,'A') != NULL) match++;
	if (ustrchr (alphabet,'C') != NULL) match++;
	if (ustrchr (alphabet,'G') != NULL) match++;
	if (ustrchr (alphabet,'T') != NULL) match++;

	if (ustrlen (alphabet) == 4)
		return (match == 4);

	if (ustrlen (alphabet) == 5)
		return ((match == 4) && (ustrchr (alphabet,'N') != NULL));

	if (ustrchr (alphabet,'a') != NULL) match++;
	if (ustrchr (alphabet,'c') != NULL) match++;
	if (ustrchr (alphabet,'g') != NULL) match++;
	if (ustrchr (alphabet,'t') != NULL) match++;

	if (ustrlen (alphabet) == 8)
		return (match == 8);

	if (ustrlen (alphabet) == 9)
		return ((match == 8) && (ustrchr (alphabet,'N') != NULL));

	return false;
	}

static u8 two_char_as_hex	// assumes both characters are valid hex digits
   (u8	ch1,
	u8	ch2)
	{
	return 16 * ((ch1<='9')? (ch1-'0') : (ch1<='F')? (10+ch1-'A') : (10+ch1-'a'))
	          + ((ch2<='9')? (ch2-'0') : (ch2<='F')? (10+ch2-'A') : (10+ch2-'a'));
	}

static int parse_bottleneck  // (returns true => success, false => failure)
   (char*	_s,
	u8		bottleneck[5])
	{
	char*	s = _s;
	int		cc;
	char	follower;
	int		i;

	// parse the four symbols, separated by spaces

	for (i=0 ; i<4 ; i++)
		{
		cc = *(s++);
		if (cc == 0) return false;
		follower = *s;

		if ((follower != 0) && (!isspace(follower)))
			{
			s++;
			if (isxdigit(cc)) cc = two_char_as_hex(cc,follower);
			             else return false;
			if (cc == 0) return false;
			}

		bottleneck[i] = cc;

		// eat up trailing whitespace

		if (*s != 0) s = skip_whitespace (s);
		}

	if (*s != 0) return false;

	return true;
	}

//----------
//
// ambiguate_n, ambiguate_iupac--
//	Change the substitution scores for N so that alignments to N treat N as an
//	ambiguous base, rather than as a sequence-splicing character.
//
//	For ambiguate_iupac, we consider any of the IUPAC ambiguous-but-not-ACGT
//	characters to be the same as an N.
//
// We expect that nVsN is zero and nVsNonN is no worse than -2*E (i.e. twice
// the gap extend penalty).  This value is chosen to give preference to this
// alignment instead of the second, gapped one (regardless of the number of
// consecutive Ns):
//
//   TTCTCttcttacttcttcttcttcttcttcttcttcttctTC
//   TTCTCTTCTNNNNNNNNNNNNNNNNNNNNNNNNNTCTTNTTC
//
//   TTCTCttct-------------------------tacttcttcttcttcttcttcttcttcttctTC
//   TTCTCTTCTNNNNNNNNNNNNNNNNNNNNNNNNNT-------------------------CTTNTTC
//
//----------
//
// Arguments:
//	scoreset*	ss:			The score set to modify.
//	score		nVsN:		Score for N-to-N substitutions.
//	score		nVsNonN:	Score for N-to-non-N substitutions.
//
// Returns:
//	(nothing)
//
//----------

void ambiguate_n
   (scoreset*	ss,
	score		nVsN,
	score		nVsNonN)
	{
	u8*			rr, *cc;
	int			ch, chLow;

	ss->sub['N']['N'] = nVsN;
	ss->sub['N']['n'] = nVsN;
	ss->sub['n']['N'] = nVsN;
	ss->sub['n']['n'] = nVsN;

	if (ss->colsAreDna)
		{
		for (rr=ss->rowChars ; *rr!=0 ; rr++)
			{
			ch = *rr;
			if (ch == 'N') continue;
			chLow = dna_tolower(ch);
			ss->sub[ch]   ['N'] = nVsNonN;
			ss->sub[ch]   ['n'] = nVsNonN;
			ss->sub[chLow]['N'] = nVsNonN;
			ss->sub[chLow]['n'] = nVsNonN;
			}
		}

	if (ss->rowsAreDna)
		{
		for (cc=ss->colChars ; *cc!=0 ; cc++)
			{
			ch = *cc;
			if (ch == 'N') continue;
			chLow = dna_tolower(ch);
			ss->sub['N'][ch]    = nVsNonN;
			ss->sub['n'][ch]    = nVsNonN;
			ss->sub['N'][chLow] = nVsNonN;
			ss->sub['n'][chLow] = nVsNonN;
			}
		}

	}


void ambiguate_iupac
   (scoreset*	ss,
	score		nVsN,
	score		nVsNonN)
	{
	static u8*	ambiggies = (u8*) "NnBDHKMRSVWYbdhkmrsvwy";
	u8*			rr, *cc;
	int			ch, chLow, rrLow, ccLow;

	// set all ambi-vs-ambi scores as N-vs-N or N-vs-nonN, as appropriate

	for (rr=ambiggies ; *rr!=0 ; rr++)
			for (cc=ambiggies ; *cc!=0 ; cc++)
		{
		rrLow = dna_tolower(*rr);
		ccLow = dna_tolower(*cc);
		if (rrLow == ccLow) ss->sub[*rr][*cc] = nVsN;
		else                ss->sub[*rr][*cc] = nVsNonN;
		}

	// set all non-ambi rows as N-vs-nonN

	if (ss->rowsAreDna)
		{
		for (rr=ss->rowChars ; *rr!=0 ; rr++)
			{
			ch    = *rr;
			chLow = dna_tolower(ch);
			for (cc=ambiggies ; *cc!=0 ; cc++)
				{
				if ((ch == 'N') && ((*cc == 'N') || (*cc == 'n'))) continue;
				ss->sub[ch]   [*cc] = nVsNonN;
				ss->sub[chLow][*cc] = nVsNonN;
				}
			}
		}

	// set all non-ambi columns as N-vs-nonN

	if (ss->colsAreDna)
		{
		for (cc=ss->colChars ; *cc!=0 ; cc++)
			{
			ch    = *cc;
			chLow = dna_tolower(ch);
			for (rr=ambiggies ; *rr!=0 ; rr++)
				{
				if ((ch == 'N') && ((*rr == 'N') || (*rr == 'n'))) continue;
				ss->sub[*rr][ch]    = nVsNonN;
				ss->sub[*rr][chLow] = nVsNonN;
				}
			}
		}

	}

//----------
//
// write_score_set_by_name, write_score_set, write_score_set_as_ints--
//	Write a new score set to a file (see format description above, in the header
//	for read_score_set_by_name).  The write_score_set_as_ints() version will
//	write the scores as though (scoreType == 'I').  This allows score sets
//	written by a floating point version of the program to be read by an integer
//	version.
//
//----------
//
// Arguments:
//	FILE*		f:		(write_score_set only) The file that scoring data is to
//						.. be written to.  This should already be open for text
//						.. write.
//	char*		name:	The name of the file that scoring data is to be written
//						.. to.  For write_score_set this is only used for
//						.. reporting problems to the user (and may be NULL).
//	scoreset*	ss:		The scoring set to write.
//	int			withGapScores: true => write gap scores too
//
// Returns:
//	(nothing)
//
//----------

static void private_write_score_set (FILE* f, scoreset* ss,
                                     int withGapScores, int asInts);


void write_score_set_by_name
   (char*		name,
	scoreset*	ss,
	int			withGapScores)
	{
	FILE*		f;

	if (name == NULL)
		suicide ("can't open NULL file in write_score_set_by_name()");

	f = fopen_or_die (name, "wt");
	write_score_set (f, name, ss, withGapScores);
	fclose_if_valid (f);
	}


void write_score_set
   (FILE* f, arg_dont_complain(char* name), scoreset* ss, int withGapScores)
	{ private_write_score_set (f, ss, withGapScores, false); }


void write_score_set_as_ints
   (FILE* f, arg_dont_complain(char* name), scoreset* ss, int withGapScores)
	{ private_write_score_set (f, ss, withGapScores, true); }


static void private_write_score_set
   (FILE*		f,
	scoreset*	ss,
	int			withGapScores,
	int			asInts)
	{
	char		s[101];
	u8*			rr, *cc;
	score		minSub;
	int			vWidth, w;
	char*		wssScoreFmt, *wssScoreFmtStar;
	score		v;

	if (asInts)
		{
		wssScoreFmt     = "%d";
		wssScoreFmtStar = "%*d";
		}
	else
		{
#if (scoreType == 'I')
		wssScoreFmt     = "%d";
		wssScoreFmtStar = "%*d";
#elif (scoreType == 'F')
		wssScoreFmt     = "%.6f";
		wssScoreFmtStar = "%*.6f";
#elif (scoreType == 'D')
		wssScoreFmt     = "%.6f";
		wssScoreFmtStar = "%*.6f";
#endif
		}

	if ((!ss->rowsAreDna) || (!ss->colsAreDna))
		suicide ("write_score_set only handles DNA scoring matrices");

	//////////
	// write non-matrix fields
	//////////

	// determine minimum substitution score

	minSub = 0;
	for (rr=ss->rowChars ; *rr!=0 ; rr++)
			for (cc=ss->colChars ; *cc!=0 ; cc++)
		{ if (ss->sub[*rr][*cc] < minSub) minSub = ss->sub[*rr][*cc]; }

	// write the fields

	if (withGapScores) vWidth = 18;
	              else vWidth = 10;

	fprintf (f, "# (a LASTZ scoring set, created by \"LASTZ --infer\")\n");
	fprintf (f, "\n");

	v = 10 * minSub;
	fprintf (f, "%-*s = %c:", vWidth, "bad_score", ss->badRow);
	if (asInts) fprintf (f, wssScoreFmt, round_score (v));
	       else fprintf (f, wssScoreFmt, v);
	fprintf (f, " # used for sub[%c][*] and sub[*][%c]\n", ss->badRow, ss->badRow);

	v = minSub;
	fprintf (f, "%-*s = ", vWidth, "fill_score");
	if (asInts) fprintf (f, wssScoreFmt, round_score (v));
	       else fprintf (f, wssScoreFmt, v);
	fprintf (f, "    # used when sub[*][*] not otherwise defined\n");

	if (withGapScores)
		{
		v = ss->gapOpen;
		fprintf (f, "%-*s = ", vWidth, "gap_open_penalty");
		if (asInts) fprintf (f, wssScoreFmt, round_score (v));
		       else fprintf (f, wssScoreFmt, v);
		fprintf (f, "\n");

		v = ss->gapExtend;
		fprintf (f, "%-*s = ", vWidth, "gap_extend_penalty");
		if (asInts) fprintf (f, wssScoreFmt, round_score (v));
		       else fprintf (f, wssScoreFmt, v);
		fprintf (f, "\n");
		}

	fprintf (f, "\n");

	//////////
	// write subsitution scores
	//////////

	// determine field width

	w = 3;
	for (rr=ss->rowChars ; *rr!=0 ; rr++)
			if ((!ss->rowsAreDna) || (dna_isupper (*rr)))
			for (cc=ss->colChars ; *cc!=0 ; cc++)
			if ((!ss->colsAreDna) || (dna_isupper (*cc)))
		{
		v = ss->sub[*rr][*cc];
		if (asInts) sprintf (s, wssScoreFmt, round_score (v));
		       else sprintf (s, wssScoreFmt, v);
		if (strleni(s)+1 > w) w = strlen(s)+1;
		}

	// write them

	fprintf (f, " ");
	for (cc=ss->colChars ; *cc!=0 ; cc++)
		{
		if ((ss->colsAreDna) && (!dna_isupper (*cc))) continue;
		fprintf (f, " %*c", w, *cc);
		}
	fprintf (f, "\n");

	for (rr=ss->rowChars ; *rr!=0 ; rr++)
		{
		if ((ss->rowsAreDna) && (!dna_isupper (*rr))) continue;
		fprintf (f, "%c", *rr);
		for (cc=ss->colChars ; *cc!=0 ; cc++)
			{
			if ((ss->colsAreDna) && (!dna_isupper (*cc))) continue;
			v = ss->sub[*rr][*cc];
			fprintf (f, " ");
			if (asInts) fprintf (f, wssScoreFmtStar, w, round_score (v));
			       else fprintf (f, wssScoreFmtStar, w, v);
			}
		fprintf (f, "\n");
		}

	}

//----------
//
// dump_dna_score_set--
//	Dump the substitution matrix of a score set that is expected to contain
//	scores for any DNA characters.  (Intended for debugging)
//
//----------
//
// Arguments:
//	FILE*		f:	The file to dumpt to.
//	scoreset*	ss:	The score set.
//
// Returns:
//	(nothing)
//
//----------

void dump_dna_score_set
   (FILE*		f,
	scoreset*	ss)
	{
	u8* alphabet = (u8*) "ACGTacgt";
	u8* ch1, *ch2;

	fprintf (stderr, "bad:  " scoreFmtSimple "\n",
	         ss->sub[ss->badRow][ss->badCol]);
	fprintf (stderr, "fill: " scoreFmtSimple "\n",
	         ss->sub['B']['B']);

	fprintf (f, "   ");
	for (ch2=alphabet ; *ch2!=0 ; ch2++)
		fprintf (f, " %5c", *ch2);
	fprintf (f, "\n");

	for (ch1=alphabet ; *ch1!=0 ; ch1++)
		{
		fprintf (f, "%c: ", *ch1);
		for (ch2=alphabet ; *ch2!=0 ; ch2++)
#if (scoreType == 'I')
			fprintf (f, " %5d",    ss->sub[*ch1][*ch2]);
#else
			fprintf (f, " %12.6f", ss->sub[*ch1][*ch2]);
#endif
		fprintf (f, "\n");
		}

	}

//----------
//
// string_to_score--
//	Parse a string for the score value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The score value of the string.  Note that the string *must not* contain
//	anything other than a valid score-- failures result in fatality.
//
//----------

score string_to_score
   (const char*	s)
	{
#if (scoreType == 'I')
	return string_to_unitized_int (s, /* byThousands */ true);
#elif (scoreType == 'F')
	return (float) string_to_double (s);
#elif (scoreType == 'D')
	return string_to_double (s);
#endif
	}

//----------
//
// scale_score_set--
//	Multiply each match/substitution score in a set by a constant.
//
//----------
//
// Arguments:
//	scoreset*	ss:		The score set to modify.  Note that *only* the match/
//						.. substitution scores are modified, and *only* for
//						.. legitimate characters (those in ss->rowChars and
//						.. ss->colChars).
//	double		scale:	How much to scale the scores (e.g. 1.0 means scores
//						.. are not changed).
//
// Returns:
//	(nothing)
//
//----------

void scale_score_set
   (scoreset*	ss,
	double		scale)
	{
	int			r, c;

	for (r=0 ; r<256 ; r++)
			for (c=0 ; c<256 ; c++)
		ss->sub[r][c] *= scale;
	}

//----------
//
// round_score--
//	Round a score value to the nearest integer.
//
//----------
//
// Arguments:
//	double	v:		The score to round.  We declare this as 'double' rather
//					.. than 'score', because the routine is most useful when
//					.. the caller has computed a score as a real, and needs to
//					.. convert it to an int.
//
// Returns:
//	The integer nearest to v.
//
//----------

int round_score
   (double	v)
	{
	if (v >= 0) return (int) (v + .5);
	       else return (int) (v - .5);
	}

//----------
//
// print_score_matrix, print_score_matrix_lf,print_score_matrix_prefix--
//	Print the meaningful contents of a score set's matrix.
// dump_score_set--
//	Dump the contents of a score set.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	scoreset*	ss:			The score set to print.
//	int			withExtras:	true => show extra stuff, such as gap_open_penalty
//	char		lineFeedCh:	The character to use to separate rows of the
//							.. matrix.
//	char*		prefix:		A string to print before each row of the matrix.
//
// Returns:
//	(nothing)
//
//----------

static void private_print_score_matrix (FILE* f, scoreset* ss, int withExtras,
                                        char lineFeedCh, char* prefix);


void print_score_matrix (FILE* f, scoreset* ss, int withExtras)
	{ private_print_score_matrix (f, ss, withExtras, '\n', NULL); }

void print_score_matrix_lf (FILE* f, scoreset* ss, int withExtras, char lineFeedCh)
	{ private_print_score_matrix (f, ss, withExtras, lineFeedCh, NULL); }

void print_score_matrix_prefix (FILE* f, scoreset* ss, int withExtras, char* prefix)
	{ private_print_score_matrix (f, ss, withExtras, '\n', prefix); }


static void private_print_score_matrix
   (FILE*		f,
	scoreset*	ss,
	int			withExtras,
	char		lineFeedCh,
	char*		prefix)
	{
	int			width;
	int			rowsHidden, rowsAsHex, colsAsHex;
	u8*			r, *c;
	char		s[3];

	if (prefix == NULL) prefix = "";

	// determine how character codes should be printed;  rowsHidden is for
	// old-style blastz compatibility, so we only set it false if any post-
	// blastz features are active

	rowsAsHex  = false;
	for (r=ss->rowChars ; *r!=0 ; r++)
		{
		if ((!isprint (*r)) || (isspace (*r)))
			{ rowsAsHex = true;  break; }
		}

	colsAsHex = false;
	for (c=ss->colChars ; *c!=0 ; c++)
		{
		if ((!isprint (*c)) || (isspace (*c)))
			{ colsAsHex = true;  break; }
		}

	rowsHidden = ((!rowsAsHex) && (!colsAsHex) && (!withExtras));

	// print assignments

	if (withExtras)
		{
		fprintf (f, "%sgap_open_penalty   = " scoreFmt "%c",
		            prefix, ss->gapOpen,   lineFeedCh);
		fprintf (f, "%sgap_extend_penalty = " scoreFmt "%c",
		            prefix, ss->gapExtend, lineFeedCh);
		}

	// print the matrix

	if (lineFeedCh != '\n')
		width = 1;
	else if ((dna_utilities_scoreType == 'F') || (dna_utilities_scoreType == 'D'))
		width = 13;
	else
		width = 4;

	fprintf (f, "%s", prefix);
	if (lineFeedCh == '\n')
		{
		fprintf (f,   (rowsHidden)? " "
		            : (rowsAsHex)?  "    "
		            :               "   ");
		}
	for (c=ss->colChars ; *c!=0 ; c++)
		{
		if ((ss->colsAreDna) && (!dna_isupper (*c))) continue;
		if (colsAsHex) sprintf (s, "%02X", *c);
		          else sprintf (s, "%c",   *c);
		fprintf (f, " %*s", width, s);
		}
	fprintf (f, "%c", lineFeedCh);

	for (r=ss->rowChars ; *r!=0 ; r++)
		{
		if ((ss->rowsAreDna) && (!dna_isupper (*r))) continue;
		fprintf (f, "%s", prefix);
		if (lineFeedCh == '\n')
			fprintf (f, (rowsAsHex)? "  " : " ");
		if (!rowsHidden)
			{
			if (rowsAsHex) sprintf (s, "%02X", *r);
			          else sprintf (s, "%c",   *r);
			fprintf (f, "%2s", s);
			}
		for (c=ss->colChars ; *c!=0 ; c++)
			{
			if ((ss->colsAreDna) && (!dna_isupper (*c))) continue;
			fprintf (f, " " scoreFmtStar, width, ss->sub[*r][*c]);
			}
		fprintf (f, "%c", lineFeedCh);
		}

	}


void dump_score_set
   (FILE*		f,
	scoreset*	ss,
	u8*			rowChars,
	u8*			colChars)
	{
	int			width;
	u8*			r, *c;

	if (rowChars == NULL) rowChars = ss->rowChars;
	if (colChars == NULL) rowChars = ss->colChars;

	if ((dna_utilities_scoreType == 'F') || (dna_utilities_scoreType == 'D'))
		width = 13;
	else
		width = 5;

	fprintf (f, "rowChars = %s\n", ss->rowChars);
	fprintf (f, "colChars = %s\n", ss->colChars);

	fprintf (f, "%2s %8s ", "", "");
	for (c=colChars ; *c!=0 ; c++)
		fprintf (f, " %*s", width, quantum_visual(*c));
	fprintf (f, "\n");

	for (r=rowChars ; *r!=0 ; r++)
		{
		fprintf (f, "%2s %8p:", quantum_visual(*r), ss->sub[*r]);
		for (c=colChars ; *c!=0 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[*r][*c]);
		fprintf (f, "\n");
		}
	}

void dump_lower_score_set
   (FILE*		f,
	scoreset*	ss)
	{
	int			width;
	u8*			r, *c;
	u8			rr;

	if ((dna_utilities_scoreType == 'F') || (dna_utilities_scoreType == 'D'))
		width = 13;
	else
		width = 5;

	fprintf (f, "%2s %8s ", "", "");
	for (c=ss->colChars ; *c!=0 ; c++)
		fprintf (f, " %*s", width, quantum_visual(*c));
	for (c=ss->colChars ; *c!=0 ; c++)
		fprintf (f, " %*s", width, quantum_visual(*c+'a'-'A'));
	fprintf (f, "\n");

	for (r=ss->rowChars ; *r!=0 ; r++)
		{
		fprintf (f, "%2s %8p:", quantum_visual(*r), ss->sub[*r]);
		for (c=ss->colChars ; *c!=0 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[*r][*c]);
		for (c=ss->colChars ; *c!=0 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[*r][*c+'a'-'A']);
		fprintf (f, "\n");
		}

	for (r=ss->rowChars ; *r!=0 ; r++)
		{
		rr = *r+'a'-'A';
		fprintf (f, "%2s %8p:", quantum_visual(rr), ss->sub[rr]);
		for (c=ss->colChars ; *c!=0 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[rr][*c]);
		for (c=ss->colChars ; *c!=0 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[rr][*c+'a'-'A']);
		fprintf (f, "\n");
		}
	}

void dump_full_score_set
   (FILE*		f,
	scoreset*	ss)
	{
	int			width;
	int			r, c;

	if ((dna_utilities_scoreType == 'F') || (dna_utilities_scoreType == 'D'))
		width = 13;
	else
		width = 4;

	fprintf (f, "%2s %8s ", "", "");
	for (c=0 ; c<256 ; c++)
		fprintf (f, " %*s", width, quantum_visual(c));
	fprintf (f, "\n");

	for (r=0 ; r<256 ; r++)
		{
		fprintf (f, "%2s %8p:", quantum_visual(r), ss->sub[r]);
		for (c=0 ; c<256 ; c++)
			fprintf (f, " " scoreFmtStar, width, ss->sub[r][c]);
		fprintf (f, "\n");
		}
	}

static char* quantum_visual (int ch)
	{
	static char	 s1[10], s2[10], s3[10];
	static char* s = s2;

	s = (s == s1)? s2 : (s == s2)? s3 : s1;	// (ping pong pung)

	if ((isprint (ch)) && (!isspace (ch))) sprintf (s, "%c",   ch);
	                                  else sprintf (s, "%02X", ch);

	return s;
	}

//----------
//
// resolve_score_thresh--
//	If an adaptive score threshold is a percentage, convert it to a count.
//
//----------
//
// Arguments:
//	sthresh*	threshold:	The threshold to resolve.
//	u32			denom:		The denominator that the threshold's percentage
//							.. would be relative to.
//
// Returns:
//	(nothing)
//
//----------

void resolve_score_thresh
   (sthresh*	threshold,
	u32			denom)
	{
	if (threshold->t != 'P') return;

	threshold->c = (u32) (threshold->p * denom + 0.5);
	threshold->t = 'C';
	}

//----------
//
// string_to_score_thresh--
//	Parse a string for the adaptive score threshold value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The score value of the string.  Note that the string *must not* contain
//	anything other than a valid adaptive score threshold-- failures result in
//	fatality.
//
//----------

sthresh string_to_score_thresh
   (const char*	s)
	{
	sthresh		threshold;
	int			len;

	memset (&threshold, 0, sizeof(threshold));	// placate compilter

	if (strcmp_prefix (s, "top") != 0)
		{
		threshold.t = 'S';
		threshold.s = string_to_score (s);
		return threshold;
		}

	len = strlen(s);
	if ((len > 3) && (s[len-1] == '%'))
		{
		threshold.t = 'P';
		threshold.p = pct_string_to_double (s+3);
		return threshold;
		}

	threshold.t = 'C';
	threshold.c = string_to_unitized_int (s+3, /* byThousands */ true);
	return threshold;
	}

//----------
//
// score_thresh_to_string--
//	Convert an adaptive score threshold to a string.
//
//----------
//
// Arguments:
//	sthresh* threshold:	The threshold to convert.
//
// Returns:
//  A string containing the threshold, as text.  This string is actually static
//  data belonging to this routine, so the caller must copy it if more than one
//  such string is to be used simultaneously.
//
//----------

char* score_thresh_to_string
   (const sthresh*	threshold)
	{
	static char	 s1[41];
	static char	 s2[41];
	static char* s = s2;

	s = (s == s1)? s2 : s1;	// (ping pong)

	if      (threshold->t == 'S') sprintf (s, scoreFmtSimple, threshold->s);
	else if (threshold->t == 'P') sprintf (s, "top%.1f%%",    100*threshold->p);
	else if (threshold->t == 'C') sprintf (s, "top%d",        threshold->c);
	else                          sprintf (s, "(unrecognized)");

	return s;
	}

//----------
//
// blastz_score_to_ncbi_bits--
//	Convert (b)lastz score to bit score in NCBI sense. 
// blastz_score_to_ncbi_expectation--
//	Convert (b)lastz score to expectation in NCBI sense.
//
// Convert (b)lastz-style scores to NCBI BLAST scores.  The conversion is not
// exact, and only provides a quick-and-dirty estimate of the corresponding
// score that would be reported by NCBI BLAST.
//
// Note:  These are borrowed with permission from the UCSC genome browser's
//        source code tree, from src/lib/blastOut.c .  Per communication with
//        Jim Kent, these are from a part of that code which is considered to
//        be in the public domain.
//
//        The routines have been syntactically modified here, to fit the style
//        of lastz's source code.  At UCSC they are called blastzScoreToNcbiBits
//        and blastzScoreToNcbiExpectation.
//
//----------
//
// Arguments:
//	score	bzScore:	The (b)lastz score to convert.
//	(none)
//
// Returns:
//	A converted score.
//
//----------

double blastz_score_to_ncbi_bits
   (score	bzScore)
	{
	return bzScore * 0.0205;
	}

double blastz_score_to_ncbi_expectation
   (score	bzScore)
	{
	double bits = bzScore * 0.0205;
	double logProb = -bits * log(2);
	return 3.0e9 * exp(logProb);
	}

//----------
//
// new_quantum_code--
//	Create a new quantum dna code.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A pointer to the newly allocated quantum dna code, which the caller will
//	have to dispose of eventually.  The routine free() should be used for this
//	purpose.
//
//----------

qcode* new_quantum_code
   (void)
	{
	return (qcode*) zalloc_or_die ("new_q_code", sizeof(qcode));
	}

//----------
//
// read_quantum_code_by_name, read_quantum_code--
//	Read a new quantum dna code from a file (see format description below).
//
//----------
//
// Arguments:
//	FILE*	f:		(read_quantum_code only) The file that code is to be read
//					.. from.  This should already be open for text read.
//	char*	name:	The name of the file that code is to be read from.  For
//					.. For read_quantum_code this is only used for reporting
//					.. problems to the user (and may be NULL).
//
// Returns:
//	A pointer to the newly allocated quantum dna code, which the caller will
//	have to dispose of eventually.  The routine free() should be used for this
//	purpose.
//
//----------
//
// Quantum Code File Format
// ========================
//
// Here's an example.  Note that blanks lines and # comments are ignored.  Rows
// begin the quantum symbol (as either a single character or a 2-digit
// hexadecimal representation).  Columns represent probabilities of A, C, G,
// and T, in that order.
//
//		01	0.125041	0.080147	0.100723	0.694088
//		02	0.111162	0.053299	0.025790	0.809749
//		03	0.065313	0.007030	0.004978	0.922679
//		 ... more rows here ...
//		FF	0.209476	0.014365	0.755682	0.020477
//
//----------

static int parse_quantum_profile (char* s, int* sym,
                                  double* pa, double* pc, double* pg, double* pt);

qcode* read_quantum_code_by_name
   (char*	name)
	{
	FILE*	f;
	qcode*	cc;

	if (name == NULL)
		suicide ("can't open NULL file in read_quantum_code_by_name()");

	f = fopen_or_die (name, "rt");
	cc = read_quantum_code (f, name);
	fclose_if_valid (f);

	return cc;
	}


qcode* read_quantum_code
   (FILE*		f,
	char*		_name)
	{
	static char	line[5*25+1];	// (must hold 5 fields, up to 25 chars each)
	char*		name = _name;
	qcode*		qc;
	int			seen[256];
	int			lineNum, len, missingEol;
	char*		waffle;
	int			sym;
	double		pa, pc, pg, pt;

	pa = pc = pg = pt = 0; // (placate compiler)

	if (name == NULL)
		name = "(unnamed file)";

	// allocate

	qc = new_quantum_code ();

	for (sym=0 ; sym<256 ; sym++)
		{
		seen[sym] = false;
		qc->p[sym][0] = qc->p[sym][1] = qc->p[sym][2] = qc->p[sym][3] = 0.0;
		}

	qc->dna[0] = bits_to_nuc[0];
	qc->dna[1] = bits_to_nuc[1];
	qc->dna[2] = bits_to_nuc[2];
	qc->dna[3] = bits_to_nuc[3];

	//////////
	// read it
	//////////

	lineNum    = 0;
	missingEol = false;

	while (fgets (line, sizeof(line), f) != NULL)
		{
		lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (missingEol)
			suicidef ("line is too long (%s: line %d)", name, lineNum-1);

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');

		// trim blanks, end of line, and comments, and ignore blank lines

		if (line[len-1] == '\n') line[--len] = 0;

		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;

		trim_string (line);
		if (line[0] == 0) continue;

		// parse it

		if (!parse_quantum_profile (line, &sym, &pa, &pc, &pg, &pt))
			suicidef ("invalid quantum code (%s: line %d) %s",
			          name, lineNum, line);

		if (seen[sym])
			suicidef ("quantum code %02X occurs more than once in %s",
			          sym, name);

		seen[sym] = true;
		qc->p[sym][0] = pa;
		qc->p[sym][1] = pc;
		qc->p[sym][2] = pg;
		qc->p[sym][3] = pt;
		}

	return qc;
	}


static int parse_quantum_profile  // (returns true => success, false => failure)
   (char*	_s,
	int*	sym,
	double*	pa,
	double*	pc,
	double*	pg,
	double*	pt)
	{
	char*	s = _s;
	u8		ch;
	double	prob, numer, denom;
	int		items, charsUsed;
	int		i;

	// parse the symbol

	charsUsed = -1;
	items = sscanf (s, "%c%n", &ch, &charsUsed);
	if ((items == 1) && (ch != 0))
		{ *sym = ch;  s += charsUsed; }
	else
		{
		charsUsed = -1;
		items = sscanf (s, "%x%n", sym, &charsUsed);
		if ((items != 1) || (*sym < 1) || (*sym > 255)) return false;
		s += charsUsed;
		}

	// parse the four probabilities

	for (i=0 ; i<4 ; i++)
		{
		charsUsed = -1;
		items = sscanf (s, " %lf/%lf%n", &numer, &denom, &charsUsed);
		if (items == 2)
			{ prob = numer / denom;  s += charsUsed; }
		else
			{
			items = sscanf (s, " %lf%n", &prob, &charsUsed);
			if (items != 1) return false;
			s += charsUsed;
			}

		switch (i)
			{
			case 0: *pa = prob;  break;
			case 1: *pc = prob;  break;
			case 2: *pg = prob;  break;
			case 3: *pt = prob;  break;
			}
		}

	if (*s != 0) return false;

	return true;
	}

//----------
//
// print_quantum_word--
//	Print a quantum word, as something like a position weight matrix.
// print_quantum_dna_match--
//	Print a quantum word with a dna word lined up with it.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	qcode*	coding:		The code mapping quantum character to probaility
//						.. vector.  This may be NULL, in which case only
//						.. a hexadecimal representation of the word is
//						.. printed.
//	u8*		q:			The quantum word to print.  This is a string of
//						characters, but is *not* zero-terminated.
//	u8*		d:			(print_quantum_dna_match only) The dna word to
//						.. print.  Like q, *not* zero-terminated.
//	u32		wordLen:	The word length (number of quantum symbols in the word).
//
// Returns:
//	(nothing).
//
//----------

void print_quantum_word (FILE* f, qcode* coding, u8* q, u32 wordLen)
	{ print_quantum_dna_match (f, coding, q, NULL, wordLen); }

void print_quantum_dna_match
   (FILE*	f,
	qcode*	coding,
	u8*		q,
	u8*		d,
	u32		wordLen)
	{
	char3	field;
	int		width = 4;	// (must be at least 4)
	u32		ix;
	u32		nuc;

	// print dna sequence

	if (d != NULL)
		{
		fprintf (f, "  ");
		for (ix=0 ; ix<wordLen ; ix++)
			fprintf (f, "%*s%c ", width-2, "", d[ix]);
		fprintf (f, "\n");
		}

	// print quantum sequence in hex

	fprintf (f, "  %s\n", quantum_word_string (q, wordLen, width));

	if (coding == NULL) return;

	// print position weight matrix

	for (nuc=0 ; nuc<sizeof(coding->dna) ; nuc++)
		{
		fprintf (f, "%c:", coding->dna[nuc]);
		for (ix=0 ; ix<wordLen ; ix++)
			{
			field = prob_to_string(coding->p[q[ix]][nuc]);
			fprintf (f, " %*s", width-1, field.s);
			}
		fprintf (f, "\n");
		}

	}

//----------
//
// quantum_word_string--
//	Convert a quantum word to a string, showing each symbol in hex.
//
//----------
//
// Arguments:
//	qcode*	coding:		The code mapping quantum character to probaility
//						.. vector.  This may be NULL, in which case only
//						.. a hexadecimal representation of the word is
//						.. printed.
//	u8*		q:			The quantum word to print.  This is a string of
//						characters, but is *not* zero-terminated.
//	u32		wordLen:	The word length (number of quantum symbols in the word).
//	int		symWidth:	Number of characters to devote to each quantum symbol.
//						.. This must be at least 2.
//
// Returns:
//	A string representing that quantum word.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The memory containing the returned string belongs to this routine, as
//		static memory.  There are only two such memory blocks, and they are
//		used on alternate calls.  So when you make more than two calls, the
//		results of previous calls are clobbered.
//
//----------

char* quantum_word_string
   (u8*		q,
	u32 	wordLen,
	int		symWidth)
	{
	static char	 s1[200];
	static char	 s2[200];
	static char* s = s2;
	char*	ss;
	u32		ix;

	s = (s == s1)? s2 : s1;	// (ping pong)

	// sanity check

	if (symWidth < 2) symWidth = 2;
	if (wordLen * symWidth + 1 > sizeof(s1))
		suicide ("internal error in quantum_word_string()");

	// print quantum sequence to the string

	ss = s;
	for (ix=0 ; ix<wordLen ; ix++)
		{ sprintf (ss, "%*s%02X", symWidth-2, "", q[ix]);  ss += symWidth; }

	return s;
	}

//----------
//
// max_in_score_matrix, min_in_score_matrix--
//	Report the maximum or minimum value in a score set's matrix, over the
//	meaningful rows and columns.
//
//----------
//
// Arguments:
//	scoreset*	ss:			The score set to examine
//
// Returns:
//	The maximum value in ss->sub[r][c], where r is in ss->rowChars and c is in
//	ss->colChars.
//
//----------

score max_in_score_matrix
   (scoreset*	ss)
	{
	u8*			r, *c;
	score		maxScore;

	maxScore = worstPossibleScore;

	for (r=ss->rowChars ; *r!=0 ; r++)
			for (c=ss->colChars ; *c!=0 ; c++)
		{ if (ss->sub[*r][*c] > maxScore) maxScore = ss->sub[*r][*c]; }

	return maxScore;
	}


score min_in_score_matrix
   (scoreset*	ss)
	{
	u8*			r, *c;
	score		minScore;

	minScore = bestPossibleScore;

	for (r=ss->rowChars ; *r!=0 ; r++)
			for (c=ss->colChars ; *c!=0 ; c++)
		{ if (ss->sub[*r][*c] < minScore) minScore = ss->sub[*r][*c]; }

	return minScore;
	}

//----------
//
// print_dna_similarities--
//	Print the similarities of two nucleotide strings (or their prefixes).
//	This prints a | when the strings contain the same nucleotide, and a ':'
//	when they contain a transition.  Spaces are printed otherwise.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print to.
//	const char*	s1:	One string.
//	const char*	s2:	The other string.
//	int			n:	The number of characters to compare (and print).
//
// Returns:
//	The number of characters printed, similar to print_prefix, except that
//	the printing is terminated at the shorter of the two strings.
//
//----------

int print_dna_similarities
   (FILE*		f,
	const char*	s1,
	const char*	s2,
	int			n)
	{
	int			ix;
	u8			ch1, ch2;

	if (n < 1) return 0;

	for (ix=0 ; ix<n ; ix++)
		{
		ch1 = (u8) s1[ix];
		ch2 = (u8) s2[ix];

		if ((ch1 == 0) || (ch2 == 0)) break;

		if (ch1 == ch2)
			fprintf (f, "|");
		else if ((nuc_to_bits[ch1] < 0) || (nuc_to_bits[ch2] < 0))
			fprintf (f, "?");
		else if ((nuc_to_bits[ch1]) == (nuc_to_bits[ch2]))
			fprintf (f, "|");
		else if ((nuc_to_bits[ch1]&1) == (nuc_to_bits[ch2]&1))
			fprintf (f, ":");
		else
			fprintf (f, " ");
		}

	return ix;
	}

//----------
//
// bits_to_nuc_string--
//	Convert a nucleotide bit string to a character string.
//
//----------
//
// Arguments:
//	u32		word:		The nucleotides, packed as two bits each.
//	int		numChars:	The number of characters that should be produced.  The
//						.. least significant 2*numChars bits of word are read.
//
// Returns:
//  A string containing the nucleotide characters.  This string is actually
//  static data belonging to this routine, so the caller must copy it if more
//  than one such string is to be used simultaneously.
//
//----------

char* bits_to_nuc_string
   (u64			word,
	int			numChars)
	{
	static char s[33];
	char*		ss;
	u32			twoBits;

	if (numChars > (int) sizeof(s)-1) numChars = sizeof(s)-1;

	// convert each bit pair to a character

	ss = s;

	while (numChars-- > 0)
		{
		twoBits = (word >> (2*numChars)) & 3;
		*(ss++) = bits_to_nuc[twoBits];
		}

	*ss = 0;
	return s;
	}

//----------
//
// entropy--
//	Compute the entropy of a sequence pair.
// entropy_lower_ok--
//	Compute the entropy of a sequence pair, with upper/lower case seen as
//	equivalent.
//
// WARNING:  These functions are only valid for DNA sequences.
//
//----------
//
// Arguments:
//	u8*	s:		One sequence.
//	u8*	t:		The other sequence.
//	int	len:	The length of the sequences.
//
// Returns:
//	The entropy of the sequence pair.
//
//----------

static double compute_entropy (u8* s, u8* t, int len, int lowerOk);

double entropy (u8* s, u8* t, int len)
	{ return compute_entropy (s, t, len, false); }

double entropy_lower_ok (u8* s, u8* t, int len)
	{ return compute_entropy (s, t, len, true); }

static double compute_entropy
   (u8*		s,
	u8*		t,
	int		len,
	int		lowerOk)
	{
#ifndef disallowEntropy
	int		count[256];
	double	pA, pC, pG, pT, qA, qC, qG, qT;
	int		ix, cA, cC, cG, cT;

	count['A'] = count['C'] = count['G'] = count['T'] = 0;
	if (lowerOk)
		count['a'] = count['c'] = count['g'] = count['t'] = 0;

	if (lowerOk)
		{
		for (ix=0; ix<len; ix++)
			{ if (dna_toupper(s[ix]) == dna_toupper(t[ix])) count[s[ix]]++; }
		count['A'] += count['a'];
		count['C'] += count['c'];
		count['G'] += count['g'];
		count['T'] += count['t'];
		}
	else
		{
		for (ix=0; ix<len; ix++)
			{ if (s[ix] == t[ix]) count[s[ix]]++; }
		}

	cA = count['A'];
	cC = count['C'];
	cG = count['G'];
	cT = count['T'];

	if ((cA+cC+cG+cT) < 20)
		return (double) 1.0;

	pA = ((double) cA) / ((double) len);
	pC = ((double) cC) / ((double) len);
	pG = ((double) cG) / ((double) len);
	pT = ((double) cT) / ((double) len);

	qA = (cA != 0)? log(pA) : 0.0;
	qC = (cC != 0)? log(pC) : 0.0;
	qG = (cG != 0)? log(pG) : 0.0;
	qT = (cT != 0)? log(pT) : 0.0;

	return -(pA*qA + pC*qC + pG*qG + pT*qT)/log(4.0);
#else
	return 0.0;
#endif // disallowEntropy
	}

//----------
//
// rev_comp_by_pairs, rev_comp_by_bits--
//	Compute the reverse complement a word, either as 2-bits per base (standard
//	nucleotide coding) or as 1 bit per base.
//
//----------
//
// Arguments:
//	u64 word:	The bits to reverse-complement (stored in the least significant
//				.. bits)
//	int length:	The number of bases in the word
//
// Returns:
//	The reverse-complement of the word, stored in the least significant bits.
//
//----------

static u8 revCompByteByPairs[256] =
	{
	0xFF,0xBF,0x7F,0x3F,0xEF,0xAF,0x6F,0x2F,0xDF,0x9F,0x5F,0x1F,0xCF,0x8F,0x4F,0x0F,
	0xFB,0xBB,0x7B,0x3B,0xEB,0xAB,0x6B,0x2B,0xDB,0x9B,0x5B,0x1B,0xCB,0x8B,0x4B,0x0B,
	0xF7,0xB7,0x77,0x37,0xE7,0xA7,0x67,0x27,0xD7,0x97,0x57,0x17,0xC7,0x87,0x47,0x07,
	0xF3,0xB3,0x73,0x33,0xE3,0xA3,0x63,0x23,0xD3,0x93,0x53,0x13,0xC3,0x83,0x43,0x03,
	0xFE,0xBE,0x7E,0x3E,0xEE,0xAE,0x6E,0x2E,0xDE,0x9E,0x5E,0x1E,0xCE,0x8E,0x4E,0x0E,
	0xFA,0xBA,0x7A,0x3A,0xEA,0xAA,0x6A,0x2A,0xDA,0x9A,0x5A,0x1A,0xCA,0x8A,0x4A,0x0A,
	0xF6,0xB6,0x76,0x36,0xE6,0xA6,0x66,0x26,0xD6,0x96,0x56,0x16,0xC6,0x86,0x46,0x06,
	0xF2,0xB2,0x72,0x32,0xE2,0xA2,0x62,0x22,0xD2,0x92,0x52,0x12,0xC2,0x82,0x42,0x02,
	0xFD,0xBD,0x7D,0x3D,0xED,0xAD,0x6D,0x2D,0xDD,0x9D,0x5D,0x1D,0xCD,0x8D,0x4D,0x0D,
	0xF9,0xB9,0x79,0x39,0xE9,0xA9,0x69,0x29,0xD9,0x99,0x59,0x19,0xC9,0x89,0x49,0x09,
	0xF5,0xB5,0x75,0x35,0xE5,0xA5,0x65,0x25,0xD5,0x95,0x55,0x15,0xC5,0x85,0x45,0x05,
	0xF1,0xB1,0x71,0x31,0xE1,0xA1,0x61,0x21,0xD1,0x91,0x51,0x11,0xC1,0x81,0x41,0x01,
	0xFC,0xBC,0x7C,0x3C,0xEC,0xAC,0x6C,0x2C,0xDC,0x9C,0x5C,0x1C,0xCC,0x8C,0x4C,0x0C,
	0xF8,0xB8,0x78,0x38,0xE8,0xA8,0x68,0x28,0xD8,0x98,0x58,0x18,0xC8,0x88,0x48,0x08,
	0xF4,0xB4,0x74,0x34,0xE4,0xA4,0x64,0x24,0xD4,0x94,0x54,0x14,0xC4,0x84,0x44,0x04,
	0xF0,0xB0,0x70,0x30,0xE0,0xA0,0x60,0x20,0xD0,0x90,0x50,0x10,0xC0,0x80,0x40,0x00
	};

static u8 revCompByteByBits[256] =
	{
	0xFF,0x7F,0xBF,0x3F,0xDF,0x5F,0x9F,0x1F,0xEF,0x6F,0xAF,0x2F,0xCF,0x4F,0x8F,0x0F,
	0xF7,0x77,0xB7,0x37,0xD7,0x57,0x97,0x17,0xE7,0x67,0xA7,0x27,0xC7,0x47,0x87,0x07,
	0xFB,0x7B,0xBB,0x3B,0xDB,0x5B,0x9B,0x1B,0xEB,0x6B,0xAB,0x2B,0xCB,0x4B,0x8B,0x0B,
	0xF3,0x73,0xB3,0x33,0xD3,0x53,0x93,0x13,0xE3,0x63,0xA3,0x23,0xC3,0x43,0x83,0x03,
	0xFD,0x7D,0xBD,0x3D,0xDD,0x5D,0x9D,0x1D,0xED,0x6D,0xAD,0x2D,0xCD,0x4D,0x8D,0x0D,
	0xF5,0x75,0xB5,0x35,0xD5,0x55,0x95,0x15,0xE5,0x65,0xA5,0x25,0xC5,0x45,0x85,0x05,
	0xF9,0x79,0xB9,0x39,0xD9,0x59,0x99,0x19,0xE9,0x69,0xA9,0x29,0xC9,0x49,0x89,0x09,
	0xF1,0x71,0xB1,0x31,0xD1,0x51,0x91,0x11,0xE1,0x61,0xA1,0x21,0xC1,0x41,0x81,0x01,
	0xFE,0x7E,0xBE,0x3E,0xDE,0x5E,0x9E,0x1E,0xEE,0x6E,0xAE,0x2E,0xCE,0x4E,0x8E,0x0E,
	0xF6,0x76,0xB6,0x36,0xD6,0x56,0x96,0x16,0xE6,0x66,0xA6,0x26,0xC6,0x46,0x86,0x06,
	0xFA,0x7A,0xBA,0x3A,0xDA,0x5A,0x9A,0x1A,0xEA,0x6A,0xAA,0x2A,0xCA,0x4A,0x8A,0x0A,
	0xF2,0x72,0xB2,0x32,0xD2,0x52,0x92,0x12,0xE2,0x62,0xA2,0x22,0xC2,0x42,0x82,0x02,
	0xFC,0x7C,0xBC,0x3C,0xDC,0x5C,0x9C,0x1C,0xEC,0x6C,0xAC,0x2C,0xCC,0x4C,0x8C,0x0C,
	0xF4,0x74,0xB4,0x34,0xD4,0x54,0x94,0x14,0xE4,0x64,0xA4,0x24,0xC4,0x44,0x84,0x04,
	0xF8,0x78,0xB8,0x38,0xD8,0x58,0x98,0x18,0xE8,0x68,0xA8,0x28,0xC8,0x48,0x88,0x08,
	0xF0,0x70,0xB0,0x30,0xD0,0x50,0x90,0x10,0xE0,0x60,0xA0,0x20,0xC0,0x40,0x80,0x00
	};

u64 rev_comp_by_pairs (u64 word, int length)
	{
	int bytes  = (2*length+7) / 8;
	int adjust = 8*bytes - 2*length;
	u64 revComp = 0;

	while (bytes-- > 0)
		{
		revComp = (revComp << 8) + revCompByteByPairs[word & 0xFF];
		word >>= 8;
		}

	return revComp >> adjust;
	}

u64 rev_comp_by_bits (u64 word, int length)
	{
	int bytes  = (length+7) / 8;
	int adjust = 8*bytes - length;
	u64 revComp = 0;

	while (bytes-- > 0)
		{
		revComp = (revComp << 8) + revCompByteByBits[word & 0xFF];
		word >>= 8;
		}

	return revComp >> adjust;
	}

//----------
//
// char_to_description--
//	Convert a character to an english description.  Only "oddball" characters
//	are described-- those characters which would not be understandable if
//	printed in an error message.
//
//----------
//
// Arguments:
//	char ch: The character to describe.
//
// Returns:
//	Pointer to a string describing the character.
//
//----------

typedef struct c2d_el
	{
	u8		ch;
	char*	description;
	} c2d_el;

static c2d_el c2dLookup[] =
	{{'!',  "exclamation point \"!\""},
	 {'"',  "double quote"},
	 {'#',  "waffle/number sign/pound \"#\""},
	 {'$',  "dollar sign \"$\""},
	 {'%',  "percent sign \"%\""},
	 {'&',  "ampersand \"&\""},
	 {'\'', "single quote/apostrophe \"'\""},
	 {'(',  "open parenthesis \"(\""},
	 {')',  "closing parenthesis \")\""},
	 {'*',  "asterisk \"*\""},
	 {'+',  "plus sign \"+\""},
	 {',',  "comma \",\""},
	 {'-',  "minus sign \"-\""},
	 {'.',  "period/dot/stop \".\""},
	 {'/',  "slash \"/\""},
	 {':',  "colon \":\""},
	 {';',  "semicolon \";\""},
	 {'<',  "less than sign \"<\""},
	 {'=',  "equals sign \"=\""},
	 {'>',  "greater than sign \">\""},
	 {'?',  "question mark \"?\""},
	 {'@',  "at sign \"@\""},
	 {'[',  "opening bracket \"[\""},
	 {'\\', "backslash \"\\\""},
	 {']',  "closing bracket \"]\""},
	 {'^',  "caret/circumflex \"^\""},
	 {'_',  "underscore \"_\""},
	 {'{',  "opening brace \"{\""},
	 {'|',  "vertical bar \"|\""},
	 {'}',  "closing brace \"}\""},
	 {'~',  "tilde/squiggle sign \"~\""},
	 {0,    NULL}};


char* char_to_description
   (char	ch)
	{
	char*			desc;
	static	char	_desc[50];
	c2d_el*			scan;

	// describe punctuation as per table

	desc = NULL;
	for (scan=c2dLookup ; scan->description!=NULL ; scan++)
		{ if (scan->ch == ch) { desc=scan->description; break; } }

	if (desc != NULL)
		return desc;

	// describe digits as "the digit X"

	if (('0' <= ch) && (ch <= '9'))
		{
		sprintf (_desc, "the digit %c", ch);
		return _desc;
		}

	// describe uppercase as "uppercase X"

	if (('A' <= ch) && (ch <= 'Z'))
		{
		sprintf (_desc, "uppercase %c", ch);
		return _desc;
		}

	// describe lowercase as "lowercase x"

	if (('a' <= ch) && (ch <= 'z'))
		{
		sprintf (_desc, "lowercase %c", ch);
		return _desc;
		}

	// describe anything else as "ascii XX"

	sprintf (_desc, "ascii %02X", ch);
	return _desc;
	}

