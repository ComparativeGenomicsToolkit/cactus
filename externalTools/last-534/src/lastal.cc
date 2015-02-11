// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "LastalArguments.hh"
#include "QualityPssmMaker.hh"
#include "OneQualityScoreMatrix.hh"
#include "TwoQualityScoreMatrix.hh"
#include "qualityScoreUtil.hh"
#include "LambdaCalculator.hh"
#include "GeneticCode.hh"
#include "SubsetSuffixArray.hh"
#include "Centroid.hh"
#include "GappedXdropAligner.hh"
#include "AlignmentPot.hh"
#include "Alignment.hh"
#include "SegmentPairPot.hh"
#include "SegmentPair.hh"
#include "ScoreMatrix.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "DiagonalTable.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "gaplessXdrop.hh"
#include "gaplessPssmXdrop.hh"
#include "gaplessTwoQualityXdrop.hh"
#include "io.hh"
#include "stringify.hh"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastal: " << x << '\n'

using namespace cbrc;

namespace {
  typedef MultiSequence::indexT indexT;
  typedef unsigned long long countT;

  LastalArguments args;
  Alphabet alph;
  Alphabet queryAlph;  // for translated alignment
  GeneticCode geneticCode;
  const unsigned maxNumOfIndexes = 16;
  SubsetSuffixArray suffixArrays[maxNumOfIndexes];
  ScoreMatrix scoreMatrix;
  GeneralizedAffineGapCosts gapCosts;
  GappedXdropAligner gappedXdropAligner;
  Centroid centroid( gappedXdropAligner );
  LambdaCalculator lambdaCalculator;
  MultiSequence query;  // sequence that hasn't been indexed by lastdb
  MultiSequence text;  // sequence that has been indexed by lastdb
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  OneQualityScoreMatrix oneQualityScoreMatrix;
  OneQualityScoreMatrix oneQualityScoreMatrixMasked;
  OneQualityExpMatrix oneQualityExpMatrix;
  QualityPssmMaker qualityPssmMaker;
  sequenceFormat::Enum referenceFormat;  // defaults to 0
  TwoQualityScoreMatrix twoQualityScoreMatrix;
  TwoQualityScoreMatrix twoQualityScoreMatrixMasked;
  int minScoreGapless;
  int isCaseSensitiveSeeds = -1;  // initialize it to an "error" value
  unsigned numOfIndexes = 1;  // assume this value, if unspecified
}

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( const std::string& matrixFile ){
  if( !matrixFile.empty() ){
    scoreMatrix.fromString( matrixFile );
  }
  else if( args.matchScore < 0 && args.mismatchCost < 0 && alph.isProtein() ){
    const char* b = args.isTranslated() ? "BLOSUM80" : "BLOSUM62";
    scoreMatrix.fromString( ScoreMatrix::stringFromName( b ) );
  }
  else{
    scoreMatrix.matchMismatch( args.matchScore, args.mismatchCost,
			       alph.letters );
  }

  scoreMatrix.init( alph.encode );

  // If the input is a PSSM, the score matrix is not used, and its
  // maximum score should not be used.  Here, we try to set it to a
  // high enough value that it has no effect.  This is a kludge - it
  // would be nice to use the maximum PSSM score.
  if( args.inputFormat == sequenceFormat::pssm ) scoreMatrix.maxScore = 10000;
  // This would work, except the maxDrops aren't finalized yet:
  // maxScore = std::max(args.maxDropGapped, args.maxDropFinal) + 1;
}

void makeQualityScorers(){
  const ScoreMatrixRow* m = scoreMatrix.caseSensitive;  // case isn't relevant
  double lambda = lambdaCalculator.lambda();
  const std::vector<double>& lp1 = lambdaCalculator.letterProbs1();
  bool isPhred1 = isPhred( referenceFormat );
  int offset1 = qualityOffset( referenceFormat );
  const std::vector<double>& lp2 = lambdaCalculator.letterProbs2();
  bool isPhred2 = isPhred( args.inputFormat );
  int offset2 = qualityOffset( args.inputFormat );

  if( referenceFormat == sequenceFormat::fasta ){
    if( isFastq( args.inputFormat ) ){
      if( args.maskLowercase > 0 )
        oneQualityScoreMatrixMasked.init( m, alph.size, lambda,
                                          &lp2[0], isPhred2, offset2,
                                          alph.numbersToUppercase, true );
      if( args.maskLowercase < 3 )
        oneQualityScoreMatrix.init( m, alph.size, lambda,
                                    &lp2[0], isPhred2, offset2,
                                    alph.numbersToUppercase, false );
      if( args.outputType > 3 ){
        const OneQualityScoreMatrix &m = (args.maskLowercase < 3) ?
            oneQualityScoreMatrix : oneQualityScoreMatrixMasked;
        oneQualityExpMatrix.init( m, args.temperature );
      }
    }
    else if( args.inputFormat == sequenceFormat::prb ){
      bool isMatchMismatch = (args.matrixFile.empty() && args.matchScore > 0);
      qualityPssmMaker.init( m, alph.size, lambda, isMatchMismatch,
                             args.matchScore, -args.mismatchCost,
                             offset2, alph.numbersToUppercase );
    }
  }
  else{
    if( isFastq( args.inputFormat ) ){
      if( args.maskLowercase > 0 )
        twoQualityScoreMatrixMasked.init( m, lambda, &lp1[0], &lp2[0],
                                          isPhred1, offset1, isPhred2, offset2,
                                          alph.numbersToUppercase, true);
      if( args.maskLowercase < 3 )
        twoQualityScoreMatrix.init( m, lambda, &lp1[0], &lp2[0],
                                    isPhred1, offset1, isPhred2, offset2,
                                    alph.numbersToUppercase, false );
      if( args.outputType > 3 ){
        ERR( "fastq-versus-fastq column probabilities not implemented" );
      }
    }
    else{
      ERR( "when the reference is fastq, the query must also be fastq" );
    }
  }
}

// Calculate statistical parameters for the alignment scoring scheme
// Meaningless for PSSMs, unless they have the same scale as the score matrix
void calculateScoreStatistics(){
  LOG( "calculating matrix probabilities..." );
  // the case-sensitivity of the matrix makes no difference here
  lambdaCalculator.calculate( scoreMatrix.caseSensitive, alph.size );
  if( lambdaCalculator.isBad() ){
    if( isQuality( args.inputFormat ) ||
        (args.temperature < 0 && args.outputType > 3) )
      ERR( "can't get probabilities for this score matrix" );
    else
      LOG( "can't get probabilities for this score matrix" );
  }else{
    LOG( "lambda=" << lambdaCalculator.lambda() );
  }
}

// Read the .prj file for the whole database
void readOuterPrj( const std::string& fileName, unsigned& volumes,
                   indexT& minSeedLimit,
                   countT& refSequences, countT& refLetters ){
  std::ifstream f( fileName.c_str() );
  if( !f ) ERR( "can't open file: " + fileName );
  unsigned version = 0;

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "version" ) iss >> version;
    if( word == "alphabet" ) iss >> alph;
    if( word == "numofsequences" ) iss >> refSequences;
    if( word == "numofletters" ) iss >> refLetters;
    if( word == "maxunsortedinterval" ) iss >> minSeedLimit;
    if( word == "masklowercase" ) iss >> isCaseSensitiveSeeds;
    if( word == "sequenceformat" ) iss >> referenceFormat;
    if( word == "volumes" ) iss >> volumes;
    if( word == "numofindexes" ) iss >> numOfIndexes;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( alph.letters.empty() || refSequences+1 == 0 || refLetters+1 == 0 ||
      isCaseSensitiveSeeds < 0 || referenceFormat >= sequenceFormat::prb ||
      numOfIndexes > maxNumOfIndexes ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) ERR( "can't read file: " + fileName );
  if( version < 294 && version > 0)
    ERR( "the lastdb files are old: please re-run lastdb" );
}

// Read a per-volume .prj file, with info about a database volume
void readInnerPrj( const std::string& fileName,
		   indexT& seqCount, indexT& seqLen ){
  std::ifstream f( fileName.c_str() );
  if( !f ) ERR( "can't open file: " + fileName );

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "numofsequences" ) iss >> seqCount;
    if( word == "numofletters" ) iss >> seqLen;
    if( word == "numofindexes" ) iss >> numOfIndexes;
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( seqCount+1 == 0 || seqLen+1 == 0 || numOfIndexes > maxNumOfIndexes ){
    f.setstate( std::ios::failbit );
  }
  if( !f ) ERR( "can't read file: " + fileName );
}

// Write match counts for each query sequence
void writeCounts( std::ostream& out ){
  LOG( "writing..." );

  for( indexT i = 0; i < matchCounts.size(); ++i ){
    out << query.seqName(i) << '\n';

    for( indexT j = args.minHitDepth; j < matchCounts[i].size(); ++j ){
      out << j << '\t' << matchCounts[i][j] << '\n';
    }

    out << '\n';  // blank line afterwards
  }
}

// Count all matches, of all sizes, of a query batch against a suffix array
void countMatches( char strand ){
  LOG( "counting..." );
  indexT seqNum = strand == '+' ? 0 : query.finishedSequences() - 1;

  for( indexT i = 0; i < query.finishedSize(); i += args.queryStep ){
    if( strand == '+' ){
      for( ;; ){
	if( seqNum == query.finishedSequences() ) return;
	if( query.seqEnd(seqNum) > i ) break;
	++seqNum;
      }
      // speed-up:
      if( args.minHitDepth > query.seqEnd(seqNum) - i ) continue;
    }
    else{
      indexT j = query.finishedSize() - i;
      for( ;; ){
	if( seqNum+1 == 0 ) return;
	if( query.seqBeg(seqNum) < j ) break;
	--seqNum;
      }
      // speed-up:
      if( args.minHitDepth > j - query.seqBeg(seqNum) ) continue;
    }

    for( unsigned x = 0; x < numOfIndexes; ++x )
      suffixArrays[x].countMatches( matchCounts[seqNum], query.seqReader() + i,
				    text.seqReader(), args.maxHitDepth );
  }
}

namespace Phase{ enum Enum{ gapless, gapped, final }; }

struct Dispatcher{
  const uchar* a;  // the reference sequence
  const uchar* b;  // the query sequence
  const uchar* i;  // the reference quality data
  const uchar* j;  // the query quality data
  const ScoreMatrixRow* p;  // the query PSSM
  const ScoreMatrixRow* m;  // the score matrix
  const TwoQualityScoreMatrix& t;
  int d;  // the maximum score drop
  int z;

  Dispatcher( Phase::Enum e ) :
      a( text.seqReader() ),
      b( query.seqReader() ),
      i( text.qualityReader() ),
      j( query.qualityReader() ),
      p( query.pssmReader() ),
      m( (e < args.maskLowercase) ?
         scoreMatrix.caseSensitive : scoreMatrix.caseInsensitive ),
      t( (e < args.maskLowercase) ?
         twoQualityScoreMatrixMasked : twoQualityScoreMatrix ),
      d( (e == Phase::gapless) ? args.maxDropGapless :
         (e == Phase::gapped ) ? args.maxDropGapped : args.maxDropFinal ),
      z( (args.inputFormat == sequenceFormat::fasta) ? 0 :
         (referenceFormat  == sequenceFormat::fasta) ? 1 : 2 ){}

  int forwardGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return forwardGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return forwardGaplessPssmXdropScore( a+x, p+y, d );
    return forwardGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  int reverseGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return reverseGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return reverseGaplessPssmXdropScore( a+x, p+y, d );
    return reverseGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  indexT forwardGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return forwardGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return forwardGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return forwardGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  indexT reverseGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return reverseGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return reverseGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return reverseGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  bool isOptimalGapless( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return isOptimalGaplessXdrop( a+x, a+e, b+y, m, d );
    if( z==1 ) return isOptimalGaplessPssmXdrop( a+x, a+e, p+y, d );
    return isOptimalGaplessTwoQualityXdrop( a+x, a+e, i+x, b+y, j+y, t, d );
  }

  int gaplessScore( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return gaplessAlignmentScore( a+x, a+e, b+y, m );
    if( z==1 ) return gaplessPssmAlignmentScore( a+x, a+e, p+y );
    return gaplessTwoQualityAlignmentScore( a+x, a+e, i+x, b+y, j+y, t );
  }
};

// Find query matches to the suffix array, and do gapless extensions
void alignGapless( SegmentPairPot& gaplessAlns,
		   char strand, std::ostream& out ){
  Dispatcher dis( Phase::gapless );
  DiagonalTable dt;  // record already-covered positions on each diagonal
  countT matchCount = 0, gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

  for( indexT i = 0; i < query.finishedSize(); i += args.queryStep ){
    for( unsigned x = 0; x < numOfIndexes; ++x ){
      const indexT* beg;
      const indexT* end;
      suffixArrays[x].match( beg, end, dis.b + i, dis.a,
			     args.oneHitMultiplicity,
			     args.minHitDepth, args.maxHitDepth );
      matchCount += end - beg;

      // Tried: if we hit a delimiter when using contiguous seeds, then
      // increase "i" to the delimiter position.  This gave a speed-up
      // of only 3%, with 34-nt tags.

      indexT gaplessAlignmentsPerQueryPosition = 0;

      for( /* noop */; beg < end; ++beg ){  // loop over suffix-array matches
	if( gaplessAlignmentsPerQueryPosition ==
	    args.maxGaplessAlignmentsPerQueryPosition ) break;

	indexT j = *beg;  // coordinate in the reference sequence

	if( dt.isCovered( i, j ) ) continue;

	int fs = dis.forwardGaplessScore( j, i );
	int rs = dis.reverseGaplessScore( j, i );
	int score = fs + rs;
	++gaplessExtensionCount;

	// Tried checking the score after isOptimal & addEndpoint, but
	// the number of extensions decreased by < 10%, and it was
	// slower overall.
	if( score < minScoreGapless ) continue;

	indexT tEnd = dis.forwardGaplessEnd( j, i, fs );
	indexT tBeg = dis.reverseGaplessEnd( j, i, rs );
	indexT qBeg = i - (j - tBeg);
	if( !dis.isOptimalGapless( tBeg, tEnd, qBeg ) ) continue;
	SegmentPair sp( tBeg, qBeg, tEnd - tBeg, score );

	if( args.outputType == 1 ){  // we just want gapless alignments
	  Alignment aln;
	  aln.fromSegmentPair(sp);
	  aln.write( text, query, strand, args.isTranslated(),
		     alph, args.outputFormat, out );
	}
	else{
	  gaplessAlns.add(sp);  // add the gapless alignment to the pot
	}

	++gaplessAlignmentsPerQueryPosition;
	++gaplessAlignmentCount;
	dt.addEndpoint( sp.end2(), sp.end1() );
      }
    }
  }

  LOG( "initial matches=" << matchCount );
  LOG( "gapless extensions=" << gaplessExtensionCount );
  LOG( "gapless alignments=" << gaplessAlignmentCount );
}

// Shrink the SegmentPair to its longest run of identical matches.
// This trims off possibly unreliable parts of the gapless alignment.
// It may not be the best strategy for protein alignment with subset
// seeds: there could be few or no identical matches...
void shrinkToLongestIdenticalRun( SegmentPair& sp, const Dispatcher& dis ){
  sp.maxIdenticalRun( dis.a, dis.b, alph.numbersToUppercase );
  sp.score = dis.gaplessScore( sp.beg1(), sp.end1(), sp.beg2() );
}

// Do gapped extensions of the gapless alignments
void alignGapped( AlignmentPot& gappedAlns, SegmentPairPot& gaplessAlns,
                  Phase::Enum phase ){
  Dispatcher dis(phase);
  indexT frameSize = args.isTranslated() ? (query.finishedSize() / 3) : 0;
  countT gappedExtensionCount = 0, gappedAlignmentCount = 0;

  // Redo the gapless extensions, using gapped score parameters.
  // Without this, if we self-compare a huge sequence, we risk getting
  // huge gapped extensions.
  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.items[i];

    int fs = dis.forwardGaplessScore( sp.beg1(), sp.beg2() );
    int rs = dis.reverseGaplessScore( sp.beg1(), sp.beg2() );
    indexT tEnd = dis.forwardGaplessEnd( sp.beg1(), sp.beg2(), fs );
    indexT tBeg = dis.reverseGaplessEnd( sp.beg1(), sp.beg2(), rs );
    indexT qBeg = sp.beg2() - (sp.beg1() - tBeg);
    sp = SegmentPair( tBeg, qBeg, tEnd - tBeg, fs + rs );

    if( !dis.isOptimalGapless( tBeg, tEnd, qBeg ) ){
      SegmentPairPot::mark(sp);
    }
  }

  erase_if( gaplessAlns.items, SegmentPairPot::isMarked );

  gaplessAlns.cull( args.cullingLimitForGaplessAlignments );
  gaplessAlns.sort();  // sort by score descending, and remove duplicates

  LOG( "redone gapless alignments=" << gaplessAlns.size() );

  for( size_t i = 0; i < gaplessAlns.size(); ++i ){
    SegmentPair& sp = gaplessAlns.get(i);

    if( SegmentPairPot::isMarked(sp) ) continue;

    Alignment aln;
    AlignmentExtras extras;  // not used
    aln.seed = sp;

    shrinkToLongestIdenticalRun( aln.seed, dis );

    // do gapped extension from each end of the seed:
    aln.makeXdrop( gappedXdropAligner, centroid, dis.a, dis.b, args.globality,
		   dis.m, scoreMatrix.maxScore, gapCosts, dis.d,
                   args.frameshiftCost, frameSize, dis.p,
                   dis.t, dis.i, dis.j, alph, extras );
    ++gappedExtensionCount;

    if( aln.score < args.minScoreGapped ) continue;

    if( !aln.isOptimal( dis.a, dis.b, args.globality, dis.m, dis.d, gapCosts,
			args.frameshiftCost, frameSize, dis.p,
                        dis.t, dis.i, dis.j ) ){
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-redundantization.
      continue;
    }

    gaplessAlns.markAllOverlaps( aln.blocks );
    gaplessAlns.markTandemRepeats( aln.seed, args.maxRepeatDistance );

    if( phase == Phase::final ) gappedAlns.add(aln);
    else SegmentPairPot::markAsGood(sp);

    ++gappedAlignmentCount;
  }

  LOG( "gapped extensions=" << gappedExtensionCount );
  LOG( "gapped alignments=" << gappedAlignmentCount );
}

// Print the gapped alignments, after optionally calculating match
// probabilities and re-aligning using the gamma-centroid algorithm
void alignFinish( const AlignmentPot& gappedAlns,
		  char strand, std::ostream& out ){
  Dispatcher dis( Phase::final );
  indexT frameSize = args.isTranslated() ? (query.finishedSize() / 3) : 0;

  if( args.outputType > 3 ){
    if( dis.p ){
      LOG( "exponentiating PSSM..." );
      centroid.setPssm( dis.p, query.finishedSize(), args.temperature,
                        oneQualityExpMatrix, dis.b, dis.j );
    }
    else{
      centroid.setScoreMatrix( dis.m, args.temperature );
    }
    centroid.setOutputType( args.outputType );
  }

  LOG( "finishing..." );

  for( size_t i = 0; i < gappedAlns.size(); ++i ){
    const Alignment& aln = gappedAlns.items[i];
    if( args.outputType < 4 ){
      aln.write( text, query, strand, args.isTranslated(),
		 alph, args.outputFormat, out );
    }
    else{  // calculate match probabilities:
      Alignment probAln;
      AlignmentExtras extras;
      probAln.seed = aln.seed;
      probAln.makeXdrop( gappedXdropAligner, centroid,
			 dis.a, dis.b, args.globality,
			 dis.m, scoreMatrix.maxScore, gapCosts, dis.d,
                         args.frameshiftCost, frameSize, dis.p, dis.t,
			 dis.i, dis.j, alph, extras,
			 args.gamma, args.outputType );
      probAln.write( text, query, strand, args.isTranslated(),
		     alph, args.outputFormat, out, extras );
    }
  }
}

void makeQualityPssm( bool isApplyMasking ){
  if( !isQuality( args.inputFormat ) || isQuality( referenceFormat ) ) return;

  LOG( "making PSSM..." );
  query.resizePssm();

  const uchar *seqBeg = query.seqReader();
  const uchar *seqEnd = seqBeg + query.finishedSize();
  const uchar *q = query.qualityReader();
  int *pssm = *query.pssmWriter();

  if( args.inputFormat == sequenceFormat::prb ){
    qualityPssmMaker.make( seqBeg, seqEnd, q, pssm, isApplyMasking );
  }
  else {
    const OneQualityScoreMatrix &m =
        isApplyMasking ? oneQualityScoreMatrixMasked : oneQualityScoreMatrix;
    makePositionSpecificScoreMatrix( m, seqBeg, seqEnd, q, pssm );
  }
}

// Scan one batch of query sequences against one database volume
void scan( char strand, std::ostream& out ){
  if( args.outputType == 0 ){  // we just want match counts
    countMatches( strand );
    return;
  }

  bool isApplyMasking = (args.maskLowercase > 0);
  makeQualityPssm(isApplyMasking);

  LOG( "scanning..." );

  SegmentPairPot gaplessAlns;
  alignGapless( gaplessAlns, strand, out );
  if( args.outputType == 1 ) return;  // we just want gapless alignments

  if( args.maskLowercase == 1 ) makeQualityPssm(false);

  AlignmentPot gappedAlns;

  if( args.maskLowercase == 2 || args.maxDropFinal != args.maxDropGapped ){
    alignGapped( gappedAlns, gaplessAlns, Phase::gapped );
    erase_if( gaplessAlns.items, SegmentPairPot::isNotMarkedAsGood );
  }

  if( args.maskLowercase == 2 ) makeQualityPssm(false);

  alignGapped( gappedAlns, gaplessAlns, Phase::final );

  if( args.outputType > 2 ){  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG( "nonredundant gapped alignments=" << gappedAlns.size() );
  }

  gappedAlns.sort();  // sort by score
  alignFinish( gappedAlns, strand, out );
}

// Scan one batch of query sequences against one database volume,
// after optionally translating the query
void translateAndScan( char strand, std::ostream& out ){
  if( args.isTranslated() ){
    LOG( "translating..." );
    const uchar* seq = query.seqReader();
    size_t size = query.finishedSize();
    std::vector<uchar> translation( size );
    geneticCode.translate( seq, seq + size, &translation[0] );
    query.swapSeq(translation);
    scan( strand, out );
    query.swapSeq(translation);
  }else{
    scan( strand, out );
  }
}

void readIndex( const std::string& baseName, indexT seqCount ) {
  LOG( "reading " << baseName << "..." );
  text.fromFiles( baseName, seqCount, isFastq( referenceFormat ) );
  for( unsigned x = 0; x < numOfIndexes; ++x ){
    if( numOfIndexes > 1 ){
      suffixArrays[x].fromFiles( baseName + char('a' + x),
				 isCaseSensitiveSeeds, alph.encode );
    }else{
      suffixArrays[x].fromFiles( baseName, isCaseSensitiveSeeds, alph.encode );
    }
  }
}

// Read one database volume
void readVolume( unsigned volumeNumber ){
  std::string baseName = args.lastdbName + stringify(volumeNumber);
  indexT seqCount = indexT(-1);
  indexT seqLen = indexT(-1);
  readInnerPrj( baseName + ".prj", seqCount, seqLen );
  minScoreGapless = args.calcMinScoreGapless( seqLen, numOfIndexes );
  readIndex( baseName, seqCount );
}

void reverseComplementPssm(){
  ScoreMatrixRow* beg = query.pssmWriter();
  ScoreMatrixRow* end = beg + query.finishedSize();

  while( beg < end ){
    --end;
    for( unsigned i = 0; i < scoreMatrixRowSize; ++i ){
      unsigned j = queryAlph.complement[i];
      if( beg < end || i < j ) std::swap( (*beg)[i], (*end)[j] );
    }
    ++beg;
  }
}

void reverseComplementQuery(){
  LOG( "reverse complementing..." );
  queryAlph.rc( query.seqWriter(), query.seqWriter() + query.finishedSize() );
  if( isQuality( args.inputFormat ) ){
    std::reverse( query.qualityWriter(),
		  query.qualityWriter() +
                  query.finishedSize() * query.qualsPerLetter() );
  }else if( args.inputFormat == sequenceFormat::pssm ){
    reverseComplementPssm();
  }
}

// Scan one batch of query sequences against all database volumes
void scanAllVolumes( unsigned volumes, std::ostream& out ){
  if( args.outputType == 0 ){
    matchCounts.clear();
    matchCounts.resize( query.finishedSequences() );
  }

  if( volumes+1 == 0 ) volumes = 1;

  for( unsigned i = 0; i < volumes; ++i ){
    if( text.unfinishedSize() == 0 || volumes > 1 ) readVolume( i );

    if( args.strand == 2 && i > 0 ) reverseComplementQuery();

    if( args.strand != 0 ) translateAndScan( '+', out );

    if( args.strand == 2 || (args.strand == 0 && i == 0) )
      reverseComplementQuery();

    if( args.strand != 1 ) translateAndScan( '-', out );
  }

  if( args.outputType == 0 ) writeCounts( out );

  LOG( "query batch done!" );
}

void writeHeader( countT refSequences, countT refLetters, std::ostream& out ){
  out << "# LAST version " <<
#include "version.hh"
      << "\n";
  out << "#\n";
  args.writeCommented( out );
  out << "# Reference sequences=" << refSequences
      << " normal letters=" << refLetters << "\n";
  out << "#\n";

  if( args.outputType == 0 ){  // we just want hit counts
    out << "# length\tcount\n";
  }
  else{  // we want alignments
    if( args.inputFormat != sequenceFormat::pssm || !args.matrixFile.empty() ){
      // we're not reading PSSMs, or we bothered to specify a matrix file
      scoreMatrix.writeCommented( out );
      // Write lambda?
      out << "#\n";
    }
    out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
    out << "# in the reverse complement of the 2nd sequence are used.\n";
    out << "#\n";

    if( args.outputFormat == 0 ){  // tabular format
      out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
	  << "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
    }
    else{  // MAF format
      out << "# name start alnSize strand seqSize alignment\n";
    }
  }

  out << "#\n";
}

// Read the next sequence, adding it to the MultiSequence
std::istream& appendFromFasta( std::istream& in ){
  indexT maxSeqLen = args.batchSize;
  if( maxSeqLen < args.batchSize ) maxSeqLen = indexT(-1);
  if( query.finishedSequences() == 0 ) maxSeqLen = indexT(-1);

  size_t oldSize = query.unfinishedSize();

  /**/ if( args.inputFormat == sequenceFormat::fasta )
    query.appendFromFasta( in, maxSeqLen );
  else if( args.inputFormat == sequenceFormat::prb )
    query.appendFromPrb( in, maxSeqLen, queryAlph.size, queryAlph.decode );
  else if( args.inputFormat == sequenceFormat::pssm )
    query.appendFromPssm( in, maxSeqLen, queryAlph.encode,
                          args.maskLowercase > 1 );
  else
    query.appendFromFastq( in, maxSeqLen );

  if( !query.isFinished() && query.finishedSequences() == 0 )
    ERR( "encountered a sequence that's too long" );

  // encode the newly-read sequence
  uchar* seq = query.seqWriter();
  size_t newSize = query.unfinishedSize();
  queryAlph.tr( seq + oldSize, seq + newSize );

  if( isPhred( args.inputFormat ) )  // assumes one quality code per letter:
    checkQualityCodes( query.qualityReader() + oldSize,
                       query.qualityReader() + newSize,
                       qualityOffset( args.inputFormat ) );

  return in;
}

void lastal( int argc, char** argv ){
  args.fromArgs( argc, argv );

  std::string matrixFile;
  if( !args.matrixFile.empty() ){
    matrixFile = ScoreMatrix::stringFromName( args.matrixFile );
    args.fromString( matrixFile );  // read options from the matrix file
    args.fromArgs( argc, argv );  // command line overrides matrix file
  }

  unsigned volumes = unsigned(-1);
  indexT minSeedLimit = 0;
  countT refSequences = -1;
  countT refLetters = -1;
  readOuterPrj( args.lastdbName + ".prj", volumes, minSeedLimit,
		refSequences, refLetters );

  if( minSeedLimit > 1 ){
    if( args.outputType == 0 )
      ERR( "can't use option -j 0: need to re-run lastdb with i <= 1" );
    if( minSeedLimit > args.oneHitMultiplicity )
      ERR( "can't use option -m < " + stringify(minSeedLimit) +
	   ": need to re-run lastdb with i <= " +
	   stringify(args.oneHitMultiplicity) );
    if( args.minHitDepth > 1 )
      ERR( "can't use option -l > 1: need to re-run lastdb with i <= 1" );
  }

  bool isMultiVolume = (volumes+1 > 0 && volumes > 1);
  args.setDefaultsFromAlphabet( alph.letters == alph.dna, alph.isProtein(),
                                isCaseSensitiveSeeds, isMultiVolume );
  makeScoreMatrix( matrixFile );
  gapCosts.assign( args.gapExistCost, args.gapExtendCost,
		   args.insExistCost, args.insExtendCost, args.gapPairCost );
  if( args.outputType > 0 ) calculateScoreStatistics();
  args.setDefaultsFromMatrix( lambdaCalculator.lambda() );
  minScoreGapless = args.calcMinScoreGapless( refLetters, numOfIndexes );
  if( !isMultiVolume ) args.minScoreGapless = minScoreGapless;
  if( args.outputType > 0 ) makeQualityScorers();

  if( args.isTranslated() ){
    if( alph.letters == alph.dna )  // allow user-defined alphabet
      ERR( "expected protein database, but got DNA" );
    queryAlph.fromString( queryAlph.dna );
    if( args.geneticCodeFile.empty() )
      geneticCode.fromString( geneticCode.standard );
    else
      geneticCode.fromFile( args.geneticCodeFile );
    geneticCode.codeTableSet( alph, queryAlph );
    query.initForAppending(3);
  }
  else{
    queryAlph = alph;
    query.initForAppending(1);
  }

  queryAlph.tr( query.seqWriter(),
                query.seqWriter() + query.unfinishedSize() );

  if( volumes+1 == 0 ) readIndex( args.lastdbName, refSequences );

  std::ostream& out = std::cout;
  writeHeader( refSequences, refLetters, out );
  out.precision(3);  // print non-integers more compactly
  countT queryBatchCount = 0;
  countT sequenceCount = 0;

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;

  for( char** i = *inputBegin ? inputBegin : defaultInput; *i; ++i ){
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while( appendFromFasta( in ) ){
      if( query.isFinished() ){
	++sequenceCount;
      }else{
        // this enables downstream parsers to read one batch at a time:
        out << "# batch " << queryBatchCount++ << "\n";
	scanAllVolumes( volumes, out );
	query.reinitForAppending();
      }
    }
  }

  if( query.finishedSequences() > 0 ){
    out << "# batch " << queryBatchCount << "\n";
    scanAllVolumes( volumes, out );
  }

  out << "# Query sequences=" << sequenceCount << "\n";
}

int main( int argc, char** argv )
try{
  lastal( argc, argv );
  if (!flush(std::cout)) ERR( "write error" );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << "lastal: out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << "lastal: " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
