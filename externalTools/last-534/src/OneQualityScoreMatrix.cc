// Copyright 2011 Martin C. Frith

#include "OneQualityScoreMatrix.hh"

#include "qualityScoreUtil.hh"

#include <algorithm>  // min
#include <cassert>
#include <cmath>
//#include <iostream>  // for debugging

namespace cbrc {

const int qualityCapacity = 128;

void OneQualityScoreMatrix::init(const ScoreMatrixRow *scoreMatrix,
                                 int numNormalLetters,
                                 double lambda,
                                 const double *letterProbs2,
                                 bool isPhred,
                                 int qualityOffset,
                                 const uchar *toUnmasked,
                                 bool isApplyMasking) {
  data.resize(qualityCapacity * scoreMatrixRowSize * scoreMatrixRowSize);

  for (int letter1 = 0; letter1 < scoreMatrixRowSize; ++letter1) {
    for (int letter2 = 0; letter2 < scoreMatrixRowSize; ++letter2) {
      int unmasked1 = toUnmasked[letter1];
      int unmasked2 = toUnmasked[letter2];

      bool isMasked = (unmasked1 != letter1 || unmasked2 != letter2);
      bool isMask = (isApplyMasking && isMasked);

      bool isNormal1 = (unmasked1 < numNormalLetters);
      bool isNormal2 = (unmasked2 < numNormalLetters);
      bool isUseQuality = (isNormal1 && isNormal2);

      int score = scoreMatrix[unmasked1][unmasked2];
      double expScore = std::exp(lambda * score);

      for (int q2 = 0; q2 < qualityCapacity; ++q2) {
        if (isUseQuality) {
          double p2 = letterProbs2[unmasked2];
          double u2 = qualityUncertainty(q2, qualityOffset, isPhred, p2);
          score = qualityPairScore(expScore, 0, u2, lambda);
        }

        if (isMask) score = std::min(score, 0);

        data[oneQualityMatrixIndex(letter1, letter2, q2)] = score;
      }
    }
  }
}

void OneQualityExpMatrix::init(const OneQualityScoreMatrix &m,
                               double temperature) {
  assert(temperature > 0);
  data.resize(qualityCapacity * scoreMatrixRowSize * scoreMatrixRowSize);

  for (int i = 0; i < scoreMatrixRowSize; ++i)
    for (int j = 0; j < scoreMatrixRowSize; ++j)
      for (int q = 0; q < qualityCapacity; ++q)
        data[oneQualityMatrixIndex(i, j, q)] =
            std::exp(m(i, j, q) / temperature);
}

void makePositionSpecificScoreMatrix(const OneQualityScoreMatrix &m,
                                     const uchar *sequenceBeg,
                                     const uchar *sequenceEnd,
                                     const uchar *qualityBeg,
                                     int *destinationBeg) {
  while (sequenceBeg < sequenceEnd) {
    int letter2 = *sequenceBeg++;
    int quality2 = *qualityBeg++;
    for (int letter1 = 0; letter1 < scoreMatrixRowSize; ++letter1)
      *destinationBeg++ = m(letter1, letter2, quality2);
  }
}

void makePositionSpecificExpMatrix(const OneQualityExpMatrix &m,
                                   const uchar *sequenceBeg,
                                   const uchar *sequenceEnd,
                                   const uchar *qualityBeg,
                                   double *destinationBeg) {
  while (sequenceBeg < sequenceEnd) {
    int letter2 = *sequenceBeg++;
    int quality2 = *qualityBeg++;
    for (int letter1 = 0; letter1 < scoreMatrixRowSize; ++letter1)
      *destinationBeg++ = m(letter1, letter2, quality2);
  }
}

}
