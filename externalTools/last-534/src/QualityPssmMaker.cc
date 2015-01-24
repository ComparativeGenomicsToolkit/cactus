// Copyright 2011 Martin C. Frith

#include "QualityPssmMaker.hh"

#include "qualityScoreUtil.hh"

#include <algorithm>  // min
#include <cassert>
#include <cmath>
//#include <iostream>  // for debugging
#include <numeric>  // inner_product

namespace cbrc {

void QualityPssmMaker::init(const ScoreMatrixRow *scoreMatrix,
                            int numNormalLetters,
                            double lambda,
                            bool isMatchMismatchMatrix,
                            int matchScore,
                            int mismatchScore,
                            int qualityOffset,
                            const uchar* toUnmasked) {
  assert(numNormalLetters <= scoreMatrixRowSize);
  assert(lambda > 0);

  this->scoreMatrix = scoreMatrix;
  this->numNormalLetters = numNormalLetters;
  this->lambda = lambda;
  this->isMatchMismatchMatrix = isMatchMismatchMatrix;
  this->toUnmasked = toUnmasked;

  double matchExp = std::exp(lambda * matchScore);
  double mismatchExp = std::exp(lambda * mismatchScore);

  for (int q = 0; q < qualityCapacity; ++q) {
    double e = qualityUncertainty(q, qualityOffset, false, 0);
    double p = 1 - e;
    qualityToProbCorrect[q] = p;

    double expScore = p * matchExp + e * mismatchExp;
    assert(p > 0);
    assert(expScore > 0);
    qualityToMatchScore[q] = nearestInt(std::log(expScore) / lambda);
  }

  for (int i = 0; i < scoreMatrixRowSize; ++i)
    for (int j = 0; j < scoreMatrixRowSize; ++j)
      expMatrix[i][j] = std::exp(lambda * scoreMatrix[i][j]);  // can be 0
}

// This function is a wee bit slow, even if isMatchMismatchMatrix and
// !isApplyMasking.
void QualityPssmMaker::make(const uchar *sequenceBeg,
                            const uchar *sequenceEnd,
                            const uchar *qualityBeg,
                            int *pssmBeg,
                            bool isApplyMasking) const {
  while (sequenceBeg < sequenceEnd) {
    double letterProbs[scoreMatrixRowSize];
    double *letterProbsEnd = letterProbs + numNormalLetters;

    for (int j = 0; j < numNormalLetters; ++j) {
      uchar q = qualityBeg[j];
      letterProbs[j] = qualityToProbCorrect[q];
    }

    int letter2 = *sequenceBeg++;
    int unmasked2 = toUnmasked[letter2];
    bool isDelimiter2 = (unmasked2 == numNormalLetters);

    for (int letter1 = 0; letter1 < scoreMatrixRowSize; ++letter1) {
      int unmasked1 = toUnmasked[letter1];
      bool isNormal1 = (unmasked1 < numNormalLetters);
      bool isUseQuality = (isNormal1 && !isDelimiter2);

      int score = scoreMatrix[unmasked1][unmasked2];

      if (isUseQuality) {
        if (isMatchMismatchMatrix) {
          // This special case is unnecessary, but faster.
          // Unfortunately, it can give slightly different results (if
          // the letter probabilities don't add up to exactly 1).
          uchar q = qualityBeg[unmasked1];
          score = qualityToMatchScore[q];
        } else {
          // I think this is the right way round for an asymmetric
          // score matrix:
          double expScore = std::inner_product(letterProbs, letterProbsEnd,
                                               expMatrix[unmasked1], 0.0);
          assert(expScore > 0);
          score = nearestInt(std::log(expScore) / lambda);
        }
      }

      if (isApplyMasking)
        if (unmasked1 != letter1 || unmasked2 != letter2)
          score = std::min(score, 0);

      *pssmBeg++ = score;
    }

    qualityBeg += numNormalLetters;
  }
}

}
