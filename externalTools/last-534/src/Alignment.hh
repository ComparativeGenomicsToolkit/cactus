// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// This struct holds a gapped, pair-wise alignment.

#ifndef ALIGNMENT_HH
#define ALIGNMENT_HH
#include "ScoreMatrixRow.hh"
#include "SegmentPair.hh"
#include <stddef.h>  // size_t
#include <string>
#include <vector>
#include <iosfwd>

namespace cbrc{

typedef unsigned char uchar;

class GappedXdropAligner;
class GeneralizedAffineGapCosts;
class MultiSequence;
class Alphabet;
class Centroid;
class TwoQualityScoreMatrix;

struct AlignmentExtras {
  // Optional (probabilistic) attributes of an alignment.
  // To save memory, these are outside the main Alignment struct.
  std::vector<uchar> columnAmbiguityCodes;  // char or uchar?
  std::vector<double> expectedCounts;  // expected emission & transition counts
  double fullScore;  // a.k.a. forward score, sum-of-paths score
  AlignmentExtras() : fullScore(0) {}
};

struct Alignment{
  // make a single-block alignment:
  void fromSegmentPair( const SegmentPair& sp );

  // Make an Alignment by doing gapped X-drop extension in both
  // directions starting from a seed SegmentPair.  The resulting
  // Alignment might not be "optimal" (see below).
  // If outputType > 3: calculates match probabilities.
  // If outputType > 4: does gamma-centroid alignment.
  void makeXdrop( GappedXdropAligner& aligner, Centroid& centroid,
		  const uchar* seq1, const uchar* seq2, int globality,
		  const ScoreMatrixRow* scoreMatrix, int smMax,
		  const GeneralizedAffineGapCosts& gap, int maxDrop,
		  int frameshiftCost, size_t frameSize,
		  const ScoreMatrixRow* pssm2,
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2,
		  const Alphabet& alph, AlignmentExtras& extras,
		  double gamma = 0, int outputType = 0 );

  // Check that the Alignment has no prefix with score <= 0, no suffix
  // with score <= 0, and no sub-segment with score < -maxDrop.
  // Alignments that pass this test may be non-optimal in other ways.
  // If "globality" is non-zero, skip the prefix and suffix checks.
  bool isOptimal( const uchar* seq1, const uchar* seq2, int globality,
                  const ScoreMatrixRow* scoreMatrix, int maxDrop,
                  const GeneralizedAffineGapCosts& gap,
		  int frameshiftCost, size_t frameSize,
		  const ScoreMatrixRow* pssm2,
                  const TwoQualityScoreMatrix& sm2qual,
                  const uchar* qual1, const uchar* qual2 );

  void write( const MultiSequence& seq1, const MultiSequence& seq2,
	      char strand, bool isTranslated, const Alphabet& alph,
	      int format, std::ostream& os,
	      const AlignmentExtras& extras = AlignmentExtras() ) const;

  // data:
  std::vector<SegmentPair> blocks;  // the gapless blocks of the alignment
  int score;
  SegmentPair seed;  // the alignment remembers its seed

  size_t beg1() const{ return blocks.front().beg1(); }
  size_t beg2() const{ return blocks.front().beg2(); }
  size_t end1() const{ return blocks.back().end1(); }
  size_t end2() const{ return blocks.back().end2(); }

  void extend( std::vector< SegmentPair >& chunks,
	       std::vector< uchar >& ambiguityCodes,
	       GappedXdropAligner& aligner, Centroid& centroid,
	       const uchar* seq1, const uchar* seq2,
	       size_t start1, size_t start2,
	       bool isForward, int globality,
	       const ScoreMatrixRow* sm, int smMax, int maxDrop,
	       const GeneralizedAffineGapCosts& gap,
	       int frameshiftCost, size_t frameSize,
	       const ScoreMatrixRow* pssm2,
               const TwoQualityScoreMatrix& sm2qual,
               const uchar* qual1, const uchar* qual2,
	       const Alphabet& alph, AlignmentExtras& extras,
	       double gamma, int outputType );

  void writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, bool isTranslated, std::ostream& os,
		 const AlignmentExtras& extras ) const;

  void writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
		 char strand, bool isTranslated, const Alphabet& alph,
		 std::ostream& os, const AlignmentExtras& extras ) const;

  size_t numColumns( size_t frameSize ) const;

  char* writeTopSeq( const uchar* seq, const Alphabet& alph,
		     size_t frameSize, char* dest ) const;

  char* writeBotSeq( const uchar* seq, const Alphabet& alph,
		     size_t frameSize, char* dest ) const;

  char* writeTopQual( const uchar* qualities,
		      size_t qualsPerBase, char* dest ) const;

  char* writeBotQual( const uchar* qualities,
		      size_t qualsPerBase, char* dest ) const;
};

}  // end namespace cbrc
#endif  // ALIGNMENT_HH
