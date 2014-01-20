//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: genpaf.h
//
//----------

#ifndef genpaf_H				// (prevent multiple inclusion)
#define genpaf_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef genpaf_owner
#define global
#else
#define global extern
#endif

// type codes for printable fields;  note that these single character codes
// are not something the user ever sees, so the only real requirements are that
// they be unique, and that they agree with definitions of genpafStandardKeys
// and other genpafXXXKeys;  they needn't be printable characters

enum
	{
	genpafAlign1           = 'A',
	genpafStart1           = 'B',
	genpafCoverage         = 'C',
	genpafStrand1          = 'D',
	genpafEnd1             = 'E',
	genpafQualsAlign1      = 'F',
	genpafGapRate          = 'G',  // deprecated;  users should use continuity
	genpafIdentity         = 'I',
	genpafTargetNucs       = 'J',
	genpafTargetQuals      = 'K',
	genpafLength1          = 'L',
	genpafName1            = 'N',
	genpafNumber1          = 'O',
	genpafStart1Blast      = 'Q',  // (we don't expect the user to set this directly)
	genpafEnd1Blast        = 'R',  // (we don't expect the user to set this directly)
	genpafSize1            = 'S',
	genpafText1            = 'T',
	genpafAlignmentColumns = 'W',
	genpafNA               = 'X',
	genpafStart1Zero       = 'Z',
	genpafAlign2           = 'a',
	genpafStart2           = 'b',
	genpafContinuity       = 'c',
	genpafStrand2          = 'd',
	genpafEnd2             = 'e',
	genpafQualsAlign2      = 'f',
	genpafIdentityPct      = 'i',
	genpafIdentityFrac     = 'k',
	genpafLength2          = 'l',
	genpafBlastIdentityPct = 'm',
	genpafName2            = 'n',
	genpafNumber2          = 'o',
	genpafQueryNucs        = 'p',
	genpafQueryQuals       = 'q',
	genpafSize2            = 's',
	genpafText2            = 't',
	genpafMatch            = 'u',
	genpafMismatch         = 'v',
	genpafAlignedPairs     = 'w',
	genpafGapColumns       = 'x',
	genpafSeparateGaps     = 'y',
	genpafStart2Zero       = 'z',
	genpafStart1DotPlot    = '0',
	genpafEnd1DotPlot      = '1',
	genpafStart2DotPlot    = '2',
	genpafEnd2DotPlot      = '3',
	genpafCoveragePct      = '6',
	genpafCoverageFrac     = '7',
	genpafContinuityPct    = '8',
	genpafContinuityFrac   = '9',
	genpafCigarLower       = '@',
	genpafCR               = '!',
	genpafScore            = '#',
	genpafBlastBitScore    = '$',  // (we don't expect the user to set this directly)
	genpafBlastEValue      = '%',  // (we don't expect the user to set this directly)
	genpafCigar            = '&',
	genpafChoreId          = '*',
	genpafEnd2OnPlus       = ',',
	genpafDiagonal         = '/',
	genpafInfoSeparator    = ';',
	genpafStart2OnPlus     = '<',
	genpafTextDiff         = '=',
	genpafStart2ZeroOnPlus = '>',
	genpafAlignmentNum     = '[',
	genpafShingle          = '\\',
	genpafAlignmentNumZero = ']',
	genpafCigarX           = '_',
	genpafCigarXLower      = '^',
	genpafMarker           = '~'
	};


#define genpafStandardKeys      "#NDSZEndszeIC"
#define genpafMappingKeys       "NZEnd>,IC^"
#define genpafSegmentKeys       "NBEnbed#"
#define genpafBlastKeys         "nNmWvybeQR%$"
#define genpafRDotplotKeys      "02!13!XX"
#define genpafRDotplotScoreKeys "02#!13#!XXX"


#define genpafTDName        "diff"
#define genpafTDInfoDefault ".:x--X"	// (indexed by these next definitons)
#define genpafTDInfoMatch        0
#define genpafTDInfoTransition   1
#define genpafTDInfoTransversion 2
#define genpafTDInfoInsert1      3
#define genpafTDInfoInsert2      4
#define genpafTDInfoOther        5
#define genpafTDInfoSize         6

#define genpafTNucsName     "nucs1"
#define genpafTQualsName    "quals1"
#define genpafQNucsName     "nucs2"
#define genpafQQualsName    "quals2"

typedef struct stringtokey
	{
	char*	name;
	char	key;
	} stringtokey;

#ifdef genpaf_owner
global stringtokey genpafName[] =
	{
	{ "name1",          genpafName1            },
	{ "number1",        genpafNumber1          },
	{ "strand1",        genpafStrand1          },
	{ "size1",          genpafSize1            },
	{ "start1",         genpafStart1           },
	{ "zstart1",        genpafStart1Zero       },
	{ "end1",           genpafEnd1             },
	{ "length1",        genpafLength1          },
	{ "align1",         genpafAlign1           },
	{ "text1",          genpafText1            },
	{ "qalign1",        genpafQualsAlign1      },
	{ "name2",          genpafName2            },
	{ "number2",        genpafNumber2          },
	{ "strand2",        genpafStrand2          },
	{ "size2",          genpafSize2            },
	{ "start2",         genpafStart2           },
	{ "zstart2",        genpafStart2Zero       },
	{ "start2+",        genpafStart2OnPlus     },
	{ "zstart2+",       genpafStart2ZeroOnPlus },
	{ "end2",           genpafEnd2             },
	{ "end2+",          genpafEnd2OnPlus       },
	{ "length2",        genpafLength2          },
	{ "align2",         genpafAlign2           },
	{ "text2",          genpafText2            },
	{ "qalign2",        genpafQualsAlign2      },
	{ "nmatch",         genpafMatch            },
	{ "nmismatch",      genpafMismatch         },
	{ "npair",          genpafAlignedPairs     },
	{ "ncolumn",        genpafAlignmentColumns },
	{ "ngap",           genpafSeparateGaps     },
	{ "cgap",           genpafGapColumns       },
	{ genpafTDName,     genpafTextDiff         },
	{ "cigar",          genpafCigar            },
	{ "cigar-",         genpafCigarLower       },
	{ "cigarx",         genpafCigarX           },
	{ "cigarx-",        genpafCigarXLower      },
	{ "diagonal",       genpafDiagonal         },
	{ "shingle",        genpafShingle          },
	{ "score",          genpafScore            },
	{ "identity",       genpafIdentity         },
	{ "idfrac",         genpafIdentityFrac     },
	{ "id%",            genpafIdentityPct      },
	{ "blastid%",       genpafBlastIdentityPct },
	{ "coverage",       genpafCoverage         },
	{ "covfrac",        genpafCoverageFrac     },
	{ "cov%",           genpafCoveragePct      },
	{ "continuity",     genpafContinuity       },
	{ "confrac",        genpafContinuityFrac   },
	{ "con%",           genpafContinuityPct    },
	{ "gaprate",        genpafGapRate          },
	{ genpafTNucsName,  genpafTargetNucs       },
	{ genpafTQualsName, genpafTargetQuals      },
	{ genpafQNucsName,  genpafQueryNucs        },
	{ genpafQQualsName, genpafQueryQuals       },
	{ "number",         genpafAlignmentNum     },
	{ "znumber",        genpafAlignmentNumZero },
	{ "chore",          genpafChoreId          },
	{ "NA",             genpafNA               },
	{ "~",              genpafMarker           },
	{ NULL,             0                      }  // (list terminator)
	};

global stringtokey genpafAliases[] =
	{
	{ "n1",          genpafName1            },
	{ "s1",          genpafStart1           },
	{ "z1",          genpafStart1Zero       },
	{ "e1",          genpafEnd1             },
	{ "l1",          genpafLength1          },
	{ "a1",          genpafAlign1           },
	{ "t1",          genpafText1            },
	{ "n2",          genpafName2            },
	{ "s2",          genpafStart2           },
	{ "z2",          genpafStart2Zero       },
	{ "s2+",         genpafStart2OnPlus     },
	{ "z2+",         genpafStart2ZeroOnPlus },
	{ "e2",          genpafEnd2             },
	{ "e2+",         genpafEnd2OnPlus       },
	{ "l2",          genpafLength2          },
	{ "a2",          genpafAlign2           },
	{ "t2",          genpafText2            },
	{ "d",           genpafDiagonal         },
	{ "diag",        genpafDiagonal         },
	{ "s",           genpafScore            },
	{ "id",          genpafIdentity         },
	{ "id%",         genpafIdentityPct      },
	{ "ident",       genpafIdentity         },
	{ "cov",         genpafCoverage         },
	{ "cov%",        genpafCoveragePct      },
	{ "con",         genpafContinuity       },
	{ "con%",        genpafContinuityPct    },
	{ "gap",         genpafGapRate          },
	{ NULL,          0                      }  // (list terminator)
	};
#else
global stringtokey genpafName[];
global stringtokey genpafAliases[];
#endif

//----------
//
// prototypes for routines in genpaf.c
//
//----------

void print_genpaf_job_header (FILE* f, char* keys);
void print_genpaf_job_footer (FILE* f);
void print_genpaf_header     (FILE* f, seq* seq1, seq* seq2);
void print_blast_job_header  (FILE* f);
void print_blast_job_footer  (FILE* f);
void print_blast_header      (FILE* f,
                              char* programName, char* args,
                              seq* seq1, seq* seq2);
void print_genpaf_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                              char* keys);
void print_genpaf_align_list_segments
                             (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                              char* keys, scoreset* scoring);
void print_genpaf_align      (FILE* f,
                              seq* seq1, unspos beg1, unspos end1,
                              seq* seq2, unspos beg2, unspos end2,
                              editscript* script, score s,
                              char* keys,
                              unspos idNumer, unspos idDenom,
                              unspos covNumer, unspos covDenom,
                              unspos conNumer, unspos conDenom,
                              unspos gapNumer, unspos gapDenom);
void print_genpaf_match      (FILE* f,
                              seq* seq1, unspos pos1,
                              seq* seq2, unspos pos2, unspos length,
                              score s,
                              char* keys);
char* parse_genpaf_keys      (char* s);

#undef global
#endif // genpaf_H
