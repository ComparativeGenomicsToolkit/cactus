// Copyright 2011 Martin C. Frith

// Score matrix for aligning letter1 to letter2, where letter2 has an
// associated quality code, but letter1 does not.

// Quality codes are transformed to quality scores like this:
// q  =  qualityCode - qualityOffset

// Quality scores are transformed to error probabilities like this:
// For phred qualities:        e(q)  =  10^(-q/10)
// Else for solexa qualities:  e(q)  =  1 / (1 + 10^(q/10))

// Error probabilities are transformed to "uncertainties" like this:
// u = e / (1 - backgroundProbability[letter2])

// The score matrix Sxy is adjusted by the quality data like this:
// Rxy    =  exp(lambda * Sxy)
// R'xyq  =  (1-u) * Rxy + u
// S'xyq  =  nearestInt[ ln(R'xyq) / lambda ]

// Some letters may be considered to be "masked" versions of other
// letters (e.g. indicated by lowercase in the original sequence).  If
// masking is not applied, masked letters get the same scores as their
// unmasked versions.  If masking is applied:
// maskedScore = min(score, 0).

// If either letter is ambiguous (e.g. 'N' for DNA), then quality data
// is ignored: S'xyq = Sxy.

#ifndef ONE_QUALITY_SCORE_MATRIX_HH
#define ONE_QUALITY_SCORE_MATRIX_HH

#include "ScoreMatrixRow.hh"

#include <vector>

namespace cbrc {

typedef unsigned char uchar;

inline int oneQualityMatrixIndex(int letter1, int letter2, int quality2) {
  return
      quality2 * scoreMatrixRowSize * scoreMatrixRowSize +
      letter2  * scoreMatrixRowSize +
      letter1;
}

class OneQualityScoreMatrix {
 public:
  void init(const ScoreMatrixRow *scoreMatrix,
            int numNormalLetters,  // typically 4 (ACGT)
            double lambda,  // scale factor for scoreMatrix
            const double *letterProbs2,  // scoreMatrix probs for 2nd sequence
            bool isPhred,  // phred or solexa qualities?
            int qualityOffset,  // typically 33 or 64
            const uchar *toUnmasked,  // maps letters to unmasked letters
            bool isApplyMasking);

  // Tests whether init has been called:
  operator const void *() const { return data.empty() ? 0 : this; }

  int operator()(int letter1, int letter2, int quality2) const
  { return data[oneQualityMatrixIndex(letter1, letter2, quality2)]; }

 private:
  std::vector<int> data;
};

// This class stores: exp(score(i, j, q) / temperature)
class OneQualityExpMatrix {
 public:
  void init(const OneQualityScoreMatrix &m, double temperature);

  // Tests whether init has been called:
  operator const void *() const { return data.empty() ? 0 : this; }

  double operator()(int letter1, int letter2, int quality2) const
  { return data[oneQualityMatrixIndex(letter1, letter2, quality2)]; }

 private:
  std::vector<double> data;
};

void makePositionSpecificScoreMatrix(const OneQualityScoreMatrix &m,
                                     const uchar *sequenceBeg,
                                     const uchar *sequenceEnd,
                                     const uchar *qualityBeg,
                                     int *destinationBeg);

void makePositionSpecificExpMatrix(const OneQualityExpMatrix &m,
                                   const uchar *sequenceBeg,
                                   const uchar *sequenceEnd,
                                   const uchar *qualityBeg,
                                   double *destinationBeg);

}

#endif
