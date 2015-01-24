// Copyright 2008, 2009, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "Alphabet.hh"
#include "Centroid.hh"
#include "GappedXdropAligner.hh"
#include "GeneticCode.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "TwoQualityScoreMatrix.hh"
#include <cassert>

// make C++ tolerable:
#define IT(type) std::vector<type>::iterator
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

void Alignment::fromSegmentPair( const SegmentPair& sp ){
  blocks.assign( 1, sp );
  score = sp.score;
}

static void addExpectedCounts( double* expectedCounts,
			       const ExpectedCount& ec,
			       const Alphabet& alph ){
  for( unsigned i = 0; i < scoreMatrixRowSize; ++i ){
    unsigned x = alph.numbersToUppercase[i];
    if( x >= alph.size ) continue;
    for( unsigned j = 0; j < scoreMatrixRowSize; ++j ){
      unsigned y = alph.numbersToUppercase[j];
      if( y >= alph.size ) continue;
      expectedCounts[ x * alph.size + y ] += ec.emit[i][j];
    }
  }

  const int numEmissionCounts = alph.size * alph.size;
  double* transitionCounts = &expectedCounts[ numEmissionCounts ];

  transitionCounts[0] += ec.MM + ec.DM + ec.IM + ec.PM;  // match count
  transitionCounts[1] += ec.DD + ec.DM + ec.DI;  // deleted letter count
  transitionCounts[2] += ec.II + ec.IM;  // inserted letter count
  transitionCounts[3] += ec.DM + ec.DI;  // deletion open/close count
  transitionCounts[4] += ec.IM;  // insertion open/close count
  transitionCounts[5] += ec.DI;  // adjacent insertion & deletion count
  transitionCounts[7] += ec.PP + ec.MP;  // unaligned letter pair count
  transitionCounts[6] += ec.MP;  // pair-gap open/close count
  transitionCounts[8] += ec.PD;
  transitionCounts[9] += ec.PI;
  // MD = DM + DI - PD
  // MI = IM - DI - PI
  // PM = MP - PD - PI
  // DM + IM + PM = MD + MI + MP
  // SM, SD, SP, SI seem to be always zero.
  // MQ, SQ ?
}

static void countSeedMatches( double* expectedCounts,
			      const uchar* seq1beg, const uchar* seq1end,
			      const uchar* seq2beg, const Alphabet& alph ){
  while( seq1beg < seq1end ){
    unsigned x1 = alph.numbersToUppercase[ *seq1beg++ ];
    unsigned x2 = alph.numbersToUppercase[ *seq2beg++ ];
    if( x1 < alph.size && x2 < alph.size )
      ++expectedCounts[ x1 * alph.size + x2 ];
  }
}

// Does x precede and touch y in both sequences?
static bool isNext( const SegmentPair& x, const SegmentPair& y ){
  return x.end1() == y.beg1() && x.end2() == y.beg2();
}

void Alignment::makeXdrop( GappedXdropAligner& aligner, Centroid& centroid,
			   const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix, int smMax,
			   const GeneralizedAffineGapCosts& gap, int maxDrop,
			   int frameshiftCost, size_t frameSize,
			   const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2,
			   const Alphabet& alph, AlignmentExtras& extras,
			   double gamma, int outputType ){
  score = seed.score;
  if( outputType > 3 ) extras.fullScore = seed.score;

  if( outputType == 7 ){
    assert( seed.size > 0 );  // makes things easier to understand
    const int numEmissionCounts = alph.size * alph.size;
    const int numTransitionCounts = 10;
    std::vector<double>& expectedCounts = extras.expectedCounts;
    expectedCounts.resize( numEmissionCounts + numTransitionCounts );
    countSeedMatches( &expectedCounts[0],
		      seq1 + seed.beg1(), seq1 + seed.end1(),
		      seq2 + seed.beg2(), alph );
    expectedCounts[ numEmissionCounts ] += seed.size;  // match count
  }

  // extend a gapped alignment in the left/reverse direction from the seed:
  std::vector<uchar>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  extend( blocks, columnAmbiguityCodes, aligner, centroid, seq1, seq2,
	  seed.beg1(), seed.beg2(), false, globality,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // avoid the bizarre-seed assert

  // convert left-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedBeg1 = seed.beg1();
  SegmentPair::indexT seedBeg2 = aaToDna( seed.beg2(), frameSize );
  for( IT(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    i->start1 = seedBeg1 - i->start1 - i->size;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedBeg2 - i->start2, frameSize ) - i->size;
  }

  // extend a gapped alignment in the right/forward direction from the seed:
  std::vector<SegmentPair> forwardBlocks;
  std::vector<uchar> forwardAmbiguities;
  extend( forwardBlocks, forwardAmbiguities, aligner, centroid, seq1, seq2,
	  seed.end1(), seed.end2(), true, globality,
	  scoreMatrix, smMax, maxDrop, gap, frameshiftCost,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // avoid the bizarre-seed assert

  // convert right-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedEnd1 = seed.end1();
  SegmentPair::indexT seedEnd2 = aaToDna( seed.end2(), frameSize );
  for( IT(SegmentPair) i = forwardBlocks.begin(); i < forwardBlocks.end();
       ++i ){
    i->start1 = seedEnd1 + i->start1;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedEnd2 + i->start2, frameSize );
  }

  bool isMergeSeedReverse = !blocks.empty() && isNext( blocks.back(), seed );
  bool isMergeSeedForward =
    !forwardBlocks.empty() && isNext( seed, forwardBlocks.back() );

  // check that the seed isn't very bizarre and dubious:
  assert( seed.size > 0 || isMergeSeedReverse || isMergeSeedForward );

  // splice together the two extensions and the seed (a bit messy):

  blocks.reserve( blocks.size() + forwardBlocks.size() +
		  1 - isMergeSeedReverse - isMergeSeedForward );

  if( isMergeSeedReverse ) blocks.back().size += seed.size;
  else                     blocks.push_back(seed);

  if( isMergeSeedForward ){
    blocks.back().size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );

  if( outputType > 3 ){  // set the un-ambiguity of the core to a max value:
    columnAmbiguityCodes.insert( columnAmbiguityCodes.end(), seed.size, 126 );
  }

  columnAmbiguityCodes.insert( columnAmbiguityCodes.end(),
                               forwardAmbiguities.rbegin(),
                               forwardAmbiguities.rend() );
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix, int maxDrop,
			   const GeneralizedAffineGapCosts& gap,
			   int frameshiftCost, size_t frameSize,
			   const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2 ){
  int maxScore = 0;
  int runningScore = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      size_t gapBeg1 = (i-1)->end1();
      size_t gapEnd1 = i->beg1();
      size_t gapSize1 = gapEnd1 - gapBeg1;

      size_t gapBeg2 = (i-1)->end2();
      size_t gapEnd2 = i->beg2();
      size_t gapSize2, frameshift2;
      sizeAndFrameshift( gapBeg2, gapEnd2, frameSize, gapSize2, frameshift2 );
      if( frameshift2 ) runningScore -= frameshiftCost;

      runningScore -= gap.cost( gapSize1, gapSize2 );
      if( !globality && runningScore <= 0 ) return false;
      if( runningScore < maxScore - maxDrop ) return false;
    }

    const uchar* s1 = seq1 + i->beg1();
    const uchar* s2 = seq2 + i->beg2();
    const uchar* e1 = seq1 + i->end1();
    const ScoreMatrixRow* p2 = pssm2 ? pssm2 + i->beg2() : 0;
    const uchar* q1 = qual1 ? qual1 + i->beg1() : 0;
    const uchar* q2 = qual2 ? qual2 + i->beg2() : 0;

    while( s1 < e1 ){
      /**/ if( sm2qual ) runningScore += sm2qual( *s1++, *s2++, *q1++, *q2++ );
      else if( pssm2 )   runningScore += ( *p2++ )[ *s1++ ];
      else               runningScore += scoreMatrix[ *s1++ ][ *s2++ ];

      if( runningScore > maxScore ) maxScore = runningScore;
      else if( !globality && runningScore <= 0 ) return false;
      else if( !globality && (s1 == e1 && i+1 == blocks.end()) ) return false;
      else if( runningScore < maxScore - maxDrop ) return false;
    }
  }

  return true;
}

void Alignment::extend( std::vector< SegmentPair >& chunks,
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
			double gamma, int outputType ){
  if( frameSize ){
    assert( outputType < 4 );
    assert( !globality );
    assert( !pssm2 );
    assert( !sm2qual );
    assert( gap.isSymmetric() );

    size_t f = aaToDna( start2, frameSize ) + 1;
    size_t r = aaToDna( start2, frameSize ) - 1;

    const uchar* frame0 = seq2 + start2;
    const uchar* frame1 = seq2 + dnaToAa( isForward ? f : r, frameSize );
    const uchar* frame2 = seq2 + dnaToAa( isForward ? r : f, frameSize );

    score += aligner.align3( seq1 + start1, frame0, frame1, frame2, isForward,
                             sm, gap.delExist, gap.delExtend, gap.pairExtend,
			     frameshiftCost, maxDrop, smMax );

    size_t end1, end2, size;
    // This should be OK even if end2 < size * 3:
    while( aligner.getNextChunk3( end1, end2, size,
				  gap.delExist, gap.delExtend, gap.pairExtend,
				  frameshiftCost ) )
      chunks.push_back( SegmentPair( end1 - size, end2 - size * 3, size ) );

    return;
  }

  int extensionScore =
    sm2qual ? aligner.align2qual( seq1 + start1, qual1 + start1,
				  seq2 + start2, qual2 + start2,
				  isForward, globality, sm2qual,
				  gap.delExist, gap.delExtend,
				  gap.insExist, gap.insExtend,
				  gap.pairExtend, maxDrop, smMax )
    : pssm2 ? aligner.alignPssm( seq1 + start1, pssm2 + start2,
				 isForward, globality,
				 gap.delExist, gap.delExtend,
				 gap.insExist, gap.insExtend,
				 gap.pairExtend, maxDrop, smMax )
    :         aligner.align( seq1 + start1, seq2 + start2,
			     isForward, globality, sm,
			     gap.delExist, gap.delExtend,
			     gap.insExist, gap.insExtend,
			     gap.pairExtend, maxDrop, smMax );

  if( extensionScore == -INF ){
    score = -INF;  // avoid score overflow
    return;  // avoid ill-defined probabilistic alignment
  }

  score += extensionScore;

  if( outputType < 5 || outputType == 7 ){  // ordinary max-score alignment
    size_t end1, end2, size;
    while( aligner.getNextChunk( end1, end2, size,
				 gap.delExist, gap.delExtend,
				 gap.insExist, gap.insExtend,
				 gap.pairExtend ) )
      chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
  }

  if( outputType > 3 ){  // calculate match probabilities
    assert( !sm2qual );
    centroid.reset();
    centroid.forward( seq1, seq2, start1, start2, isForward, globality, gap );
    centroid.backward( seq1, seq2, start1, start2, isForward, globality, gap );

    if( outputType > 4 && outputType < 7 ){  // gamma-centroid / LAMA alignment
      centroid.dp( gamma );
      centroid.traceback( chunks, gamma );
    }

    centroid.getColumnAmbiguities( ambiguityCodes, chunks, isForward );
    extras.fullScore += centroid.logPartitionFunction();

    if( outputType == 7 ){
      ExpectedCount ec;
      centroid.computeExpectedCounts( seq1, seq2, start1, start2,
				      isForward, gap, ec );
      addExpectedCounts( &extras.expectedCounts[0], ec, alph );
    }
  }
}
