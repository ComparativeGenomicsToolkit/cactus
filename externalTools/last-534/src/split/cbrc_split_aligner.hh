// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

#ifndef CBRC_SPLIT_ALIGNER_HH
#define CBRC_SPLIT_ALIGNER_HH

#include "cbrc_unsplit_alignment.hh"
#include "cbrc_int_exponentiator.hh"

#include "Alphabet.hh"
#include "MultiSequence.hh"

#include <string>
#include <vector>
#include <cmath>
#include <climits>
#include <stddef.h>  // size_t
#include <map>

namespace cbrc {

class SplitAligner {
public:
    SplitAligner() {
      setParams(-7, -1, -7, -1, -40, INT_MIN/2, 1.0, 33);  // xxx ???
      setSpliceParams(0.0, 7.0, 1.7);
    }

    // A gap of length k scores: gapExistenceScore + k * gapExtensionScore.

    // We allow for the possibility that insertions get different scores:
    // insExistenceScore + k * insExtensionScore.

    // "jumpScore" is the (negative) score for a trans-splice.

    // "restartScore" is the (negative) score for re-starting an
    // alignment, in the repeated matches algorithm from chapter 2 of
    // Durbin, Eddy, et al.  In that book, it is called "-T".

    // "qualityOffset" is 33 for fastq-sanger or 64 for fastq-illumina

    void setParams(int gapExistenceScoreIn, int gapExtensionScoreIn,
		   int insExistenceScoreIn, int insExtensionScoreIn,
		   int jumpScoreIn, int restartScoreIn, double scaleIn,
		   int qualityOffsetIn);

    void setSpliceParams(double splicePriorIn,
			 double meanLogDistIn, double sdevLogDistIn);

    void setScoreMat(const std::vector< std::vector<int> >& matrix,
		     const std::string& rowNames,
		     const std::string& colNames);

    void readGenome(const std::string& baseName);

    // XXX this should allow us to specify scores for gt-ag, at-ac, etc.
    void setSpliceSignals();

    // Outputs some algorithm parameters on lines starting with "#"
    void printParameters() const;

    // Prepares to analyze some candidate alignments for one query sequence
    void initForOneQuery(std::vector<UnsplitAlignment>::const_iterator beg,
			 std::vector<UnsplitAlignment>::const_iterator end);

    long viterbi();  // returns the optimal split-alignment score

    // Gets the chunks of an optimal split alignment.
    // For each chunk, it gets:
    // 1. The index of the candidate alignment that the chunk comes from
    // 2. The chunk's start coordinate in the query sequence
    // 3. The chunk's end coordinate in the query sequence
    // It gets the chunks in reverse order, from query end to query start.
    void traceBack(long viterbiScore,
		   std::vector<unsigned>& alnNums,
		   std::vector<unsigned>& queryBegs,
		   std::vector<unsigned>& queryEnds) const;

    // Calculates the alignment score for a segment of an alignment
    int segmentScore(unsigned alnNum,
		     unsigned queryBeg, unsigned queryEnd) const;

    void forward();

    void backward();

    // Returns one probability per column, for a segment of an alignment
    std::vector<double> marginalProbs(unsigned queryBeg, unsigned alnNum,
				      unsigned alnBeg, unsigned alnEnd) const;

    // Toggles between forward and reverse-complement splice signals
    void flipSpliceSignals();

    // The probability that the query uses splice signals in the
    // orientation currently set by flipSpliceSignals()
    double spliceSignalStrandProb() const;

private:
    static const int numQualCodes = 64;
    static int score_mat[64][64][numQualCodes];
    int maxMatchScore;
    int qualityOffset;
    int gapExistenceScore;
    int gapExtensionScore;
    int insExistenceScore;
    int insExtensionScore;
    int jumpScore;
    int restartScore;
    double jumpProb;
    double restartProb;
    double scale;
    IntExponentiator scaledExp;  // for fast calculation of exp(x / scale)
    unsigned numAlns;  // the number of candidate alignments (for 1 query)
    std::vector<UnsplitAlignment>::const_iterator alns;  // the candidates
    unsigned minBeg;  // the minimum query start coordinate of any candidate
    unsigned maxEnd;  // the maximum query end coordinate of any candidate
    std::vector<unsigned> dpBegs;  // dynamic programming begin coords
    std::vector<unsigned> dpEnds;  // dynamic programming end coords
    std::vector<size_t> matrixRowOrigins;  // layout of ragged matrices
    std::vector<int> Amat;  // scores at query bases, for each candidate
    std::vector<int> Dmat;  // scores between query bases, for each candidate
    std::vector<long> Vmat;  // DP matrix for Viterbi algorithm
    std::vector<long> Vvec;  // DP vector for Viterbi algorithm
    std::vector<double> Aexp;
    std::vector<double> Dexp;
    std::vector<double> Fmat;  // DP matrix for Forward algorithm
    std::vector<double> Bmat;  // DP matrix for Backward algorithm
    std::vector<double> rescales;  // the usual scaling for numerical stability

    std::vector<long> VmatRev;
    std::vector<long> VvecRev;
    std::vector<double> FmatRev;
    std::vector<double> BmatRev;
    std::vector<double> rescalesRev;

    std::vector<unsigned> sortedAlnIndices;
    std::vector<unsigned> oldInplayAlnIndices;
    std::vector<unsigned> newInplayAlnIndices;

    double splicePrior;
    double meanLogDist;
    double sdevLogDist;
    double spliceTerm1;
    double spliceTerm2;
    unsigned maxSpliceDist;
    std::vector<unsigned> spliceBegCoords;
    std::vector<unsigned> spliceEndCoords;
    std::vector<unsigned char> spliceBegSignals;
    std::vector<unsigned char> spliceEndSignals;
    std::vector<unsigned> rBegs;  // genomic beg coordinate of each candidate
    std::vector<unsigned> rEnds;  // genomic end coordinate of each candidate
    std::vector<unsigned> rnameAndStrandIds;
    std::vector<int> spliceScoreTable;  // lookup table
    std::vector<double> spliceProbTable;  // lookup table
    unsigned spliceTableSize;
    MultiSequence genome;
    Alphabet alphabet;
    typedef std::map<std::string, unsigned> StringNumMap;
    StringNumMap chromosomeIndex;
    int spliceBegScores[4 * 4 + 1];  // donor score for any dinucleotide
    int spliceEndScores[4 * 4 + 1];  // acceptor score for any dinucleotide
    double spliceBegProbs[4 * 4 + 1];
    double spliceEndProbs[4 * 4 + 1];
    unsigned spliceBegSignal(unsigned coordinate, char strand) const;
    unsigned spliceEndSignal(unsigned coordinate, char strand) const;
    int spliceBegScore(unsigned i, unsigned j) const {
      if (chromosomeIndex.empty()) return 0;
      return spliceBegScores[cell(spliceBegSignals, i, j)];
    }
    int spliceEndScore(unsigned i, unsigned j) const {
      if (chromosomeIndex.empty()) return 0;
      return spliceEndScores[cell(spliceEndSignals, i, j)];
    }
    double spliceBegProb(unsigned i, unsigned j) const {
      if (chromosomeIndex.empty()) return 1;
      return spliceBegProbs[cell(spliceBegSignals, i, j)];
    }
    double spliceEndProb(unsigned i, unsigned j) const {
      if (chromosomeIndex.empty()) return 1;
      return spliceEndProbs[cell(spliceEndSignals, i, j)];
    }
    int calcSpliceScore(double dist) const;
    int spliceScore(unsigned d) const
    { return d < spliceTableSize ? spliceScoreTable[d] : calcSpliceScore(d); }
    double calcSpliceProb(double dist) const
    { return scaledExp(calcSpliceScore(dist)); }
    double spliceProb(unsigned d) const
    { return d < spliceTableSize ? spliceProbTable[d] : calcSpliceProb(d); }
    void initSpliceCoords();
    void initSpliceSignals();
    void initRnameAndStrandIds();

    int maxJumpScore() const;

    void updateInplayAlnIndicesF(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    void updateInplayAlnIndicesB(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    unsigned findScore(unsigned j, long score) const;
    unsigned findSpliceScore(unsigned i, unsigned j, long score) const;
    long scoreFromSplice(unsigned i, unsigned j, unsigned oldNumInplay,
			 unsigned& oldInplayPos) const;
    long endScore() const;
    unsigned findEndScore(long score) const;

    // "dp" means "dynamic programming":
    unsigned dpBeg(unsigned i) const { return dpBegs[i]; }
    unsigned dpEnd(unsigned i) const { return dpEnds[i]; }

    template<typename T> T&
    cell(std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    template<typename T> const T&
    cell(const std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    // cell j in row i of a ragged matrix
    template<typename T> T&
    cell(std::vector<T>& v, unsigned i, unsigned j) const
    { return v[matrixRowOrigins[i] + j]; }

    // cell j in row i of a ragged matrix
    template<typename T> const T&
    cell(const std::vector<T>& v, unsigned i, unsigned j) const
    { return v[matrixRowOrigins[i] + j]; }

    template<typename T>
    void resizeVector(T& v) const
    { v.resize(maxEnd - minBeg + 1); }

    template<typename T>
    void resizeMatrix(T& m) const {
      // This reserves size for a ragged matrix, which is actually
      // stored in a flat vector.  There are numAlns rows, and row i
      // has dpEnd(i) - dpBeg(i) + 1 cells.  (The final cell per row
      // is used in some matrices but not others.)
      m.resize(matrixRowOrigins[numAlns-1] + dpEnd(numAlns-1) + 1);
    }

    double probFromSpliceF(unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    double probFromSpliceB(unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    void calcBaseScores(unsigned i);
    void calcInsScores(unsigned i);
    void calcDelScores(unsigned i);
    void calcScoreMatrices();
    void initForwardBackward();
    void initDpBounds();

    long scoreIndel(unsigned i, unsigned j) const {
      return cell(Vmat, i, j) + cell(Dmat, i, j);
    }
};

}

#endif
