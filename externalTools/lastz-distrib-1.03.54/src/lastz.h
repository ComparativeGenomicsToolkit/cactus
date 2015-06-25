//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: lastz.h
//
//----------

#ifndef lastz_H					// (prevent multiple inclusion)
#define lastz_H

// other files

#include "dna_utilities.h"		// dna/scoring stuff
#include "seeds.h"				// seed strategy stuff
#include "pos_table.h"			// position table stuff
#include "capsule.h"			// multi-process sharing stuff
#include "seed_search.h"		// seed hit search stuff
#include "gapped_extend.h"		// gapped alignment stuff
#include "masking.h"			// dynamic masking stuff
#include "diag_hash.h"			// diagonals hashing stuff

// establish ownership of global variables

#ifdef lastz_owner
#define global
#else
#define global extern
#endif

//----------
//
// data structures and types
//
//----------

global int debug;				// how much debug info to show the programmer
								//   0   => nothing
								//   100 => everything


typedef struct infcontrol
	{
	char*	inferFilename;		// Name of file to write inferred scores to.  If
								// .. this includes "%s", then each score set we
								// .. try during convergence will be written,
								// .. with %s replaced by a string of the form
								// .. sNNN or gNNN.  The final score set will
								// .. be written with %s replaced by an empty
								// .. string
	score	inferScale;			// The desired value for the maximum subsitution
								// .. score.  If this is zero, no scaling is
								// .. performed.
	int		writeAsInt;			// true => we should produce scores that are
								//         .. integers even if the score type
								//         .. is floating-point
	int		hspThresholdIsRatio;// non-zero => the hspThreshold and
	int		gappedThresholdIsRatio;//  .. gappedThreshold value(s) passed into
								//     .. drive_scoring_inference is to be used
								//     .. as a ratio, as per ratioXXX value
	int		gapOpenIsRatio;		// non-zero => the scoring->gapOpen and
	int		gapExtendIsRatio;	//     .. scoring->gapExtend value(s) passed
								//     .. into drive_scoring_inference is to be
								//     .. used as a ratio, as per ratioXXX value
	int		subIterations;		// max number of iterations for inferring
								// .. substitution scores
	int		gapIterations;		// max number of iterations for inferring gap
								// .. scores
	int		idIsPercentile;		// true => control.minIdentity and
								//         .. control.maxIdentity are
								//         .. percentiles
								// false => minIdentity and maxIdentity are
								//         .. the actual percentages
	} infcontrol;

enum
	{
	ratioNone = 0,
	ratioMaxSubScore,
	ratioMinSubScore,
	};

// control--
//	data structure for providing global information to all parts of the lastz
//	program;  most of these are control parameters, but some information about
//	the problem being solved is also passed here
//
//	items marked here with an (I) are also valid for control of the scoring
//	inference process;  items marked (I-only) are *only* valid for inference;
//	items marked (I-copy) are valid for inference but are copies of pointers
//	to memory that belongs to the regular control set

typedef struct control
	{
	// the sequences being processed

	char*		seq1Filename;	// target sequence
	struct seq*	seq1;

	char*		seq2Filename;	// query sequence
	struct seq*	seq2;

	u8*			rev1;			// (I) reverse of target and query sequences
	u8*			rev2;			//     .. (NOT reverse complement)

	//	control parameters

	int			inferScores;	//     true => infer scores from the sequences
	int			inferOnly;		//     true => prohibit alignment (other than
								//             .. what is performed for
								//             .. inference)
	infcontrol	ic;				// (I-only) Additional inference controls.

	int			selfCompare;	// (I) true  => align with the knowledge that
								//              .. seq1 and seq2 are the same
								//              .. sequence
	int			clonedQuery;	// (I) true  => seq1 and seq2 are cloned from
								//              .. the same sequence
	int			doSeedSearch;	// (I) false => prohibit the seed search (and
								//              .. prohibit alignment too)
	const s8*	charToBits;		//     tables to map sequence characters to
	const s8*	upperCharToBits;//     .. two-bit values, and illegal characters
								//     .. to -1;  indexed by a u8 value, 0..255;
								//     .. charToBits will normally consider
								//     .. upper and lower case the same;
								//     .. upperCharToBits will normally reject
								//     .. lower case
	int			whichStrand;	// (I) 0   => search + strand only
								//     < 0 => search - strand only
								//     > 0 => search both strands of target
	u32			step;			// (I) positional step size, indicating how
								//     .. often to store target word positions

	struct seed* hitSeed;		// (I) seeding strategy for hits
	int			maxIndexBits;	// (I) maximum number of index bits to use for
								//     .. the position table;  if the seed
								//     .. weight (in bits) is larger than this,
								//     .. the additional bits (effective hash
								//     .. collisions) are resolved using the
								//     .. sequences

	int			withTrans;		// (I) number of allowed transitions in a seed
								//     .. hit (0, 1 or 2)
	int			noHitFiltering;	//     true => just report every raw seed hit,
								//             .. with no filtering or
								//             .. processing
	u32			twinMinSpan;	// (I) span threshold for hits to be considered
	u32			twinMaxSpan;	// (I) .. twins;  twinMinSpan<=0 means we don't
								//     .. don't require twins (and we ignore
								//     .. twinMaxSpan);  see the twininfo struct
								//     .. in seed_search.h for details
	int			basicHitType;	// (I) the type of hit we require for a basic,
								//     .. single hit;  one of the hitXXX values
								//     .. below
	int			minMatches;		// (I) filter criteria for each seed hit;  we
	int			maxTransversions;//(I) .. require at least minMatches matches
	int			filterCaresOnly;// (I) .. and no more than maxTransversions
								//     .. transversions;  if minMatches<0 no
								//     .. filtering is performed;  if
								//     .. maxTransversions<0 there is no
								//     .. transversion limit; if
								//     .. filterCaresOnly is true the filter
								//     .. criteria only applies to the seed's
								//     .. "care" positions (and not to any
								//     .. don't-care positions)
#ifndef noSeedHitQueue
	int			seedHitQueueSize;//    number of entries to allocate for the
								//     .. seedHitQueue
#endif // not noSeedHitQueue

	int			readCapsule;	//     true => read capsule file
	int			writeCapsule;	//     true => write capsule file
	FILE*		capsuleFile;	//     file to write capsule to (can be NULL)
	char*		capsuleFilename;//     name of capsule file (can be NULL)
	capinfo*	capsule;		//     info about open capsule file
	u64			targetMem;		//     number of bytes to pre-allocate for
								//     .. the target sequence record (zero
								//     .. indicates no pre-allocation)
	u64			queryMem;		//     number of bytes to pre-allocate for
								//     .. the query sequence record (zero
								//     .. indicates no pre-allocation)

	FILE*		anchorsFile;	//     file to read anchors from instead of
								//     .. discovering them via seeding (can be
								//     .. NULL)
	char*		anchorsFilename;//     name of anchors file (can be NULL)
	char*		choresFilename;	//     name of chores file (can be NULL)

	int			gfExtend;		//     whether to extend seed hits into HSPs
								//     .. (one of gfexXXX, from seed_search.h)
	int			mergeAnchors;	//     true => after seed hits/HSPs are found,
								//             .. merge any overlapping segments
	int			chain;			// (I) true => perform chaining
	score		chainDiag;		// (I) diagonal chaining penalty
	score		chainAnti;		// (I) antidiagonal chaining penalty
	int			gappedExtend;	// (I) true => perform gapped extension

	scoreset*	scoring;		// (I) scoring set
	scoreset*	maskedScoring;	// (I) scoring set with lowercase penalized;  in
								//     .. general, we treat lower case as bad
								//     .. during the search for HSPs, then treat
								//     .. upper/lower case as equivalent during
								//     .. the gapped alignment stage
	score		xDrop;			// (I) threshold to stop extensions of
								//     .. *ungapped* matches;  if the score
								//     .. drops off by more than xDrop,
								//     .. extension stops
	score		yDrop;			// (I) threshold to stop extensions of *gapped*
								//     .. alignments (similar to xdrop)
	int			xDropUntrimmed;	// (I) true => xDrop extension does *not* trim
								//             .. to the scoring peak if it
								//             .. happens to run into the end
								//             .. of either sequence
	int			yDropUntrimmed; // (I) true => yDrop extension does *not* trim
								//             .. to the scoring peak if it
								//             .. happens to run into the end
								//             .. of either sequence
	sthresh		hspThreshold;	// (I) threshold for high scoring segment pairs
								//     .. (ungapped matches);  an HSP is
								//     .. discarded if its score is less than
								//     .. this threshold
	sthresh		gappedThreshold;// (I) threshold for high scoring alignments;  a
								//     .. gapped alignment is discarded if its
								//     .. score is less than this threshold
	int			entropicHsp;	// (I) true  => involve entropy in the decision
								//              .. to discard an HSP
	int			reportEntropy;	//     true  => report any HSPs that are
								//              .. discarded due to entropy
	int			gappedAllBounds;// (I) true  => bound gapped alignments by *all*
								//              .. gapped extensions of higher-
								//              .. scoring HSPs (a la blastz)
								//     false => bound gapped alignments only by
								//              .. gapped extensions that meet
								//              .. the score threshold
	int			mirrorHSP;		// (I) true  => each HSP or seed hit is output
								//              .. with its symmetric copy,
								//              .. reflected over the diagonal
	int			mirrorGapped;	// (I) true  => each gapped alignment is output
								//              .. with its symmetric copy,
								//              .. reflected over the diagonal
	int			inhibitTrivial;	// (I) true  => don't output the trivial self-
								//              .. alignment
	u32			tracebackMem;	//     number of bytes to allocate to track
								//     .. gapped alignment traceback
	tback*		traceback;		//     memory in which to track gapped alignment
								//     .. traceback
	int			nIsAmbiguous;	//     true  => N is an ambiguous nucleotide
								//     false => N is a sequence-splicer
	int			allowAmbiDNA;	//     true  => permit ambiguous DNA characters
								//              .. B,D,H,K,M,R,S,V,W,Y
								//     false => only A,C,G,T,N,X permitted
	score		ambiMatch;		//     (non-negative) penalty for matches
								//     .. among ambiguous DNA
	score		ambiMismatch;	//     (non-negative) penalty for mismatches
								//     .. among ambiguous DNA
	int			hspImmediate;	//     true => process HSPs immediately  (rather
								//             .. than collecting them in a
								//             .. table then processing en
								//             .. masse);  if gapped extension
								//             .. is to be performed, it is
								//             .. done immediately (rather than
								//             .. collecting them in a table
								//             .. then extending en masse)
	u32			searchLimit;	//     maximum number of "HSPs" allowed for a
								//     .. given query/strand;  zero indicates
								//     .. "no limit";  "HSPs" means gapped
								//     .. alignments if hspImmediate is true
	int			searchLimitWarn;//     true => warn user about reads that
								//             .. exceed searchLimit
	int			searchLimitKeep;//     true => report alignments for reads that
								//             .. exceed searchLimit (up to the
								//             .. limit)
	u32			numBestHsps;	//     maximum number of HSPs processed for a
								//     .. given query/strand;  zero indicates
								//     .. "no limit";  if this is non-zero, only
								//     .. the best N HSPs are kept after sorting
	float		maxPairedDepth;	//     maximum alignment "depth" we'll allow
								//     .. in one call to gapped_extend();  this
								//     .. relates to P/L, were P is the number
								//     .. of "paired bases" and L is the length
								//     .. of the sequence;  if this is zero it
								//     .. has no effect
	u64			maxPairedBases;	//     maximum number of "paired bases" we'll
								//     .. allow in one call to gapped_extend();
								//     .. if this is zero it has no effect,
								//     .. otherwise, it overrides maxPairedDepth
	int			overlyPairedWarn;//    true  => warn user about reads that
								//              .. exceed exceed the maximum
								//              .. alignment "depth"
	int			overlyPairedKeep;//    false => we discard all alignments for
								//              .. reads that exceed the
								//              .. maximum alignment "depth"
								//     true  => we output whatever alignments
								//              .. we happen to find prior to
								//              .. exceeding the limit

	float		wordCountKeep;	//     if non-zero, this specifies how to set
								//     .. wordCountLimit adaptively;  this value
								//     .. is between 0 and 1 and indicates a
								//     .. lower bound on the fraction of seed
								//     .. word positions that we will keep
	u32			wordCountLimit;	//     the maximum number of positions that a
								//     .. particular seed word can occur in the
								//     .. target;  words that occur more than
								//     .. this is deleted from the table;  if
								//     .. wordCountKeep is non-zero, this
								//     .. variable is adaptively set;  zero
								//     .. means no limit
	u32			maxWordCountChasm;//   if non-zero, this specifies the maximum
								//     .. length of an interval of discarded
								//     .. seed word positions (discarded due to
								//     .. wordCountKeep or wordCountLimit);  if
								//     .. there is an interval longer than this,
								//     .. we protect some seed word positions
								//     .. that would otherwise be discarded
	u32			dynamicMasking;	//     mask any position in target hit this many
								//     .. times;  zero indicates no masking
	FILE*		maskingFile;	//     file to write masked intervals to (can be
								//     .. NULL)
	char*		maskingFilename;//     name of masked intervals file (can be
								//     .. NULL)
	int			masking3Fields;	//     true  => write 3-field masked intervals
								//     false => write 2-field masked intervals
	FILE*		softMaskedFile;	//     file to write soft-masked intervals to
								//     .. (can be NULL)
	char*		softMaskedFilename;//  name of soft-masked intervals file (can
								//     .. be NULL)
	int			softMasked3Fields;//   true  => write 3-field soft-masked intervals
								//     false => write 2-field soft-masked intervals
	int			reportCensus;	//     true => report how many times each base
								//             .. in the target is aligned
	FILE*		censusFile;		//     file to write census to (can be NULL)
	char*		censusFilename;	//     name of census file (can be NULL)
	char		censusKind;		//     size of counting type to use for census
								//     .. (see the definition of the census type
								//     .. in masking.h).

	float		minIdentity;	// (I) the range of percent identity of HSPS we
	float		maxIdentity;	// (I) .. will keep (or alignment blocks if we
								//     .. are doing gapped extension);  these
								//     .. are values between 0.0 and 1.0
	float		minCoverage;	// (I) the range of query coverage of HSPS we
	float		maxCoverage;	// (I) .. will keep (or alignment blocks if we
								//     .. are doing gapped extension);  these
								//     .. are values between 0.0 and 1.0
	float		minContinuity;	// (I) the range of query continuity of HSPS we
	float		maxContinuity;	// (I) .. will keep (or alignment blocks if we
								//     .. are doing gapped extension);  these
								//     .. are values between 0.0 and 1.0
	double   minMatchCountRatio;// (I) ratio of match-count/query-length;  this
								//     .. is only valid if non-zero; and is used
								//     .. to reset minMatchCount as each query
								//     ... is loaded
	u32			minMatchCount;	// (I) the minimum match-count of HSPS we will
								//     .. keep (or alignment blocks if we are
								//     .. doing gapped extension)
	s32			maxMismatchCount;//(I) the maximum number of mismatches we'll
								//     .. allow in HSPS we keep (or alignment
								//     .. blocks if we are doing gapped
								//     .. extension);  -1 indicates we have no
								//     .. limit
	s32		maxSeparateGapsCount;//(I) the maximum number of gaps we'll
								//     .. allow, counting each run of gapped
								//     .. columns as one gap, in HSPS we keep
								//     .. (or alignment blocks if we are doing
								//     .. gapped extension);  -1 indicates we
								//     .. have no limit
	s32		maxGapColumnsCount;	// (I) the maximum number of gaps we'll
								//     .. allow, counting each gapped column
								//     .. as one gap, in HSPS we keep (or
								//     .. alignment blocks if we are doing
								//     .. gapped extension);  -1 indicates we
								//     .. have no limit

#ifdef densityFiltering
	float		maxDensity;		//     the maximum alignment density we will
								//     .. allow before discarding a query
								//     .. sequence;  zero means no limit
#endif // densityFiltering

	char*		outputFilename;	//     name of the file to write output to
	FILE*		outputFile;		//     file to write output to;  this may be
								//     .. stdout
	int			outputFormat;	//     format in which to write alignments;  one
								//     .. of fmtXXX values (in output.h)
	void*		outputInfo;		//     additional information for the particular
								//     .. output format chosen
	char*		readGroup;		//     additional information for SAM format
	char*		samRGTags;		//     additional information for SAM format
	int			endComment;		//     true => write a comment at the end of the
								//             .. output file, so that the user
								//             .. can tell we completed
	int			needTrueLengths;//     true => we need seq->trueLen to be
								//             .. correct (e.g. because we have
								//             .. an anchors file, or because
								//             .. outputFormat requires it)
	int			deGapifyOutput;	//     true => convert gapped alignments back to
								//             .. ungapped segments on output
	char*		dotplotFilename;//     name of the file to write dot plot to
	FILE*		dotplotFile;	//     file to write dot plot to
	char*		dotplotKeys;	//     genpaf keys for formatting dotplot

	// for inner alignment (interpolation)

	score		innerThreshold;	//     threshold for HSPs during interpolation
								//     .. stage (equivalent to hspThreshold);
								//     .. zero or negative indicates no
								//     .. interpolation
	struct seed* innerSeed;		//     seeding strategy for inner alignment hits
	u32			innerWindow;	//     windowSize interpolation parameter (see
								//     .. inner_interpolate() and inner.c)

	// for quantum DNA query

	int			targetIsQuantum;//     true => the target is quantum DNA
	int			queryIsQuantum;	//     true => the query  is quantum DNA
	score		ballScore;		//     minimum score required of a DNA word to
								//     .. be considered 'in' the ball
								//     .. surrounding a given quantum word;
								//     .. -1 indicates that this was never set

	// debug/info flags

	int			lajCompatible;	//     true => backward compatibility for LAJ
	u32			textContext;	//     the number of extra bp to print at the
								//     .. ends of matches for
								//     .. outputFormat=fmtText or fmtZeroText
	char*		args;			//     a copy of the argv options, less sequence
								//     .. names, suitable for printing
	int			verbosity;		// (I) how much info to bombard the user with
								//       0   => minimal
								//       10  => everything
	int			reportTiming;	// (I) true => report runtime to the output file
								//             .. (e.g as z records in a gfa
								//             ..  file)
	int			reportStats;	// (I) true => report search statistics to the
								//             .. output file (e.g as z records
								//             .. in a gfa file)
	int			showStats;		// (I) true => show search statistics
	FILE*		statsFile;		// (I-copy) file to write statistics to (can be
								//     .. NULL)
	char*		statsFilename;	//     name of stats file (can be NULL)
	int			showPosTable;	//     whether or not to show target positions
								//     .. table (one of spt_xxx values)
	} control;

// values for control parameters

enum
	{
	hitBad  = -1,
	hit_min = 0,
	hitSimple = hit_min,		// simple hit, no hash collision detection
	hitRecover,					// simple hit with hash collision recovery
	hit_max = hitRecover
	};

enum
	{
	spt_dont         = 0,		// don't show target positions table
	spt_table        = 1,		// show target positions table
	spt_countsonly   = 2,		// show target positions table counts only
	spt_withcounts   = 3,		// show target positions table and counts
	spt_distribution = 4		// show position counts distribution
	};

//----------
//
// prototypes for routines in lastz.c
//
//----------

void set_up_hit_processor (control* params, int collectingCensus,
                           hitprocessor* hitProc, void** hitProcInfo);

int  start_one_strand     (seq* target, postable* targPositions, seq* query,
                           int emptyAnchors, u32 prevAnchorsCount,
                           hitprocessor hitProc, void* hitProcInfo);
void finish_one_strand    (seq* target, u8* targetRev,
                           postable* targPositions,
                           seq* query,  u8* queryRev,
                           tback* traceback, census* targCensus);
void split_anchors        (int id);
void swap_anchor_sets     (void);
void print_job_header     (void);

//----------
//
// statistics for events in this module
//
//----------

#ifdef collect_stats

global struct
	{
	float runTime;
	} lastzStats;

// stats macros

#define lastz_count_stat(field)   ++lastzStats.field
#define lastz_uncount_stat(field) --lastzStats.field
#define lastz_set_stat(field,val) (lastzStats.field = val)
#else
#define lastz_count_stat(field)
#define lastz_uncount_stat(field)
#define lastz_set_stat(field,val)
#endif // collect_stats

#undef global
#endif // lastz_H
