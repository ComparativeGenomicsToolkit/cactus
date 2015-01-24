// Copyright 2011 Martin C. Frith

#include "TwoQualityScoreMatrix.hh"

#include "qualityScoreUtil.hh"

#include <algorithm>  // min
#include <cassert>
#include <cmath>

namespace cbrc {

void TwoQualityMatrixIndexer::init(const uchar *toUnmasked) {
  indexMap.resize(qualityCapacity * scoreMatrixRowSize);

  for (int quality = 0; quality < qualityCapacity; ++quality) {
    int normalStart = quality * numQualityLetters;
    int maskedStart = normalStart + numNormalLetters;
    int abnormalPos = 0;

    for (int letter = 0; letter < scoreMatrixRowSize; ++letter) {
      int unmasked = toUnmasked[letter];
      int i = subindex(letter, quality);

      if (letter < numNormalLetters)
        indexMap[i] = normalStart + letter;
      else if (unmasked < numNormalLetters)
        indexMap[i] = maskedStart + unmasked;
      else
        indexMap[i] = abnormalPos++;
    }

    assert(abnormalPos <= minQuality * numQualityLetters);
  }
}

static int qualityEnd(const TwoQualityMatrixIndexer &indexer, int letter) {
  if (indexer(0, letter, 0, 0) == indexer(0, letter, 0, 1))
    return indexer.minQuality + 1;
  else
    return indexer.qualityCapacity;
}

void TwoQualityScoreMatrix::init(const ScoreMatrixRow *scoreMatrix,
                                 double lambda,
                                 const double *letterProbs1,
                                 const double *letterProbs2,
                                 bool isPhred1,
                                 int qualityOffset1,
                                 bool isPhred2,
                                 int qualityOffset2,
                                 const uchar *toUnmasked,
                                 bool isApplyMasking) {
  indexer.init(toUnmasked);
  data.resize(indexer.numSymbols * indexer.numSymbols);

  for (int letter1 = 0; letter1 < scoreMatrixRowSize; ++letter1) {
    for (int letter2 = 0; letter2 < scoreMatrixRowSize; ++letter2) {
      int unmasked1 = toUnmasked[letter1];
      int unmasked2 = toUnmasked[letter2];

      bool isMasked = (unmasked1 != letter1 || unmasked2 != letter2);
      bool isMask = (isApplyMasking && isMasked);

      bool isNormal1 = (unmasked1 < indexer.numNormalLetters);
      bool isNormal2 = (unmasked2 < indexer.numNormalLetters);
      bool isUseQuality = (isNormal1 && isNormal2);

      int score = scoreMatrix[unmasked1][unmasked2];
      double expScore = std::exp(lambda * score);

      int end1 = qualityEnd(indexer, letter1);
      int end2 = qualityEnd(indexer, letter2);

      for (int q1 = indexer.minQuality; q1 < end1; ++q1) {
        for (int q2 = indexer.minQuality; q2 < end2; ++q2) {
          if (isUseQuality) {
            double p1 = letterProbs1[unmasked1];
            double u1 = qualityUncertainty(q1, qualityOffset1, isPhred1, p1);
            double p2 = letterProbs2[unmasked2];
            double u2 = qualityUncertainty(q2, qualityOffset2, isPhred2, p2);
            score = qualityPairScore(expScore, u1, u2, lambda);
          }

          if (isMask) score = std::min(score, 0);

          data[indexer(letter1, letter2, q1, q2)] = score;
        }
      }
    }
  }
}

void TwoQualityExpMatrix::init(const TwoQualityScoreMatrix &m,
                               double temperature) {
  assert(temperature > 0);
  indexer = m.indexer;
  data.resize(indexer.numSymbols * indexer.numSymbols);

  for (int i1 = 0; i1 < scoreMatrixRowSize; ++i1) {
    for (int i2 = 0; i2 < scoreMatrixRowSize; ++i2) {
      int end1 = qualityEnd(indexer, i1);
      int end2 = qualityEnd(indexer, i2);

      for (int q1 = indexer.minQuality; q1 < end1; ++q1)
        for (int q2 = indexer.minQuality; q2 < end2; ++q2)
          data[indexer(i1, i2, q1, q2)] =
              std::exp(m(i1, i2, q1, q2) / temperature);
    }
  }
}

}
