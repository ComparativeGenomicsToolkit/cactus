// Copyright 2008, 2009, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#ifndef CENTROID_HH
#define CENTROID_HH
#include "GappedXdropAligner.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "SegmentPair.hh"
#include "OneQualityScoreMatrix.hh"
#include <cassert>
#include <stddef.h>  // size_t
#include <vector>
#include <iostream> // for debug

namespace cbrc{

  struct ExpectedCount{
  public:
    double emit[scoreMatrixRowSize][scoreMatrixRowSize];
    double MM, MD, MP, MI, MQ;
    double DD, DM, DI;
    double PP, PM, PD, PI;
    double II, IM;
    double SM, SD, SP, SI, SQ;
  public:
    ExpectedCount ();
    std::ostream& write (std::ostream& os) const;
  };
  /**
   * (1) Forward and backward algorithm on the DP region given by Xdrop algorithm
   * (2) \gamma-centroid decoding
   */
  class Centroid{
  public:
    Centroid( const GappedXdropAligner& xa_ );

    // Setters
    void setScoreMatrix( const ScoreMatrixRow* sm, double T );
    void setPssm ( const ScoreMatrixRow* pssm, size_t qsize, double T,
                   const OneQualityExpMatrix& oqem,
                   const uchar* sequenceBeg, const uchar* qualityBeg );
    void setOutputType( int m ) { outputType = m; }

    void reset( ) {
      numAntidiagonals = xa.numAntidiagonals();
      bestScore = 0;
      bestAntiDiagonal = 0;
      bestPos1 =0;
    }

    void forward( const uchar* seq1, const uchar* seq2,
		  size_t start1, size_t start2, bool isForward,
		  int globality, const GeneralizedAffineGapCosts& gap );

    void backward( const uchar* seq1, const uchar* seq2,
		   size_t start1, size_t start2, bool isForward,
		   int globality, const GeneralizedAffineGapCosts& gap );

    double dp( double gamma );
    void traceback( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_centroid( double gamma );
    void traceback_centroid( std::vector< SegmentPair >& chunks, double gamma ) const;

    double dp_ama( double gamma );
    void traceback_ama( std::vector< SegmentPair >& chunks, double gamma ) const;

    // Added by MCF: get the probability of each column in the alignment:
    void getColumnAmbiguities( std::vector< uchar >& ambiguityCodes,
                               const std::vector< SegmentPair >& chunks,
                               bool isForward );

    double logPartitionFunction() const;  // a.k.a. full score, forward score

    // Added by MH (2008/10/10) : compute expected counts for transitions and emissions
    void computeExpectedCounts ( const uchar* seq1, const uchar* seq2,
				 size_t start1, size_t start2, bool isForward,
				 const GeneralizedAffineGapCosts& gap,
				 ExpectedCount& count ) const;

  private:
    typedef double ExpMatrixRow[scoreMatrixRowSize];

    const GappedXdropAligner& xa;
    double T; // temperature
    size_t numAntidiagonals;
    double match_score[scoreMatrixRowSize][scoreMatrixRowSize];
    bool isPssm;
    std::vector<double> pssmExp; //
    ExpMatrixRow* pssmExp2; // pre-computed pssm for prob align
    int outputType;

    typedef std::vector< double > dvec_t;

    dvec_t fM; // f^M(i,j)
    dvec_t fD; // f^D(i,j) Ix
    dvec_t fI; // f^I(i,j) Iy
    dvec_t fP; // f^P(i,j)

    dvec_t bM; // b^M(i,j)
    dvec_t bD; // b^D(i,j)
    dvec_t bI; // b^I(i,j)
    dvec_t bP; // b^P(i,j)

    dvec_t pp; // posterior match probabilities

    dvec_t mD;
    dvec_t mI;
    dvec_t mX1;
    dvec_t mX2;

    dvec_t X; // DP tables for $gamma$-decoding

    dvec_t scale; // scale[n] is a scaling factor for the n-th anti-diagonal

    double bestScore;
    size_t bestAntiDiagonal;
    size_t bestPos1;

    void initForwardMatrix();
    void initBackwardMatrix();
    void initDecodingMatrix();

    void updateScore( double score, size_t antiDiagonal, size_t cur );

    // get DP matrix value at the given position
    double cellx( const dvec_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      return matrix[ xa.scoreEndIndex( antiDiagonal ) + seq1pos - xa.seq1start( antiDiagonal ) ];
    }

    // get DP matrix value "left of" the given position
    double horix( const dvec_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      assert( antiDiagonal > 0 );
      return cellx( matrix, antiDiagonal-1, seq1pos );
    }

    // get DP matrix value "above" the given position
    double vertx( const dvec_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      assert( antiDiagonal > 0 );
      return cellx( matrix, antiDiagonal-1, seq1pos+1 );
    }

    // get DP matrix value "diagonal from" the given position
    double diagx( const dvec_t& matrix,
                  size_t antiDiagonal, size_t seq1pos ) const{
      return cellx( matrix, antiDiagonal-2, seq1pos );
    }

    // get a pointer into a sequence, taking start and direction into account
    template < typename T >
    static const T* seqPtr( const T* seq, bool isForward, size_t pos ){
      if( isForward ) return seq + pos;
      else            return seq - pos - 1;
    }
  };

}  // end namespace cbrc
#endif  // CENTROID_HH
