// Copyright 2011 Martin C. Frith

// This class calculates a Position-Specific-Scoring-Matrix from
// a sequence with quality data.

// The sequence "letters" are actually small integers.  The first
// numNormalLetters letters are considered to be "normal" (e.g. ACGT
// for DNA).  The remaining letters are "abnormal" (e.g. ambiguous,
// delimiter).

// There must be numNormalLetters quality codes per letter (e.g. from
// PRB format).

// Quality codes are transformed to quality scores like this:
// q  =  qualityCode - qualityOffset

// Quality scores are transformed to error probabilities like this:
// e(q)  =  1 / (1 + 10^(q/10))

// The score matrix Sij is adjusted by the quality data like this:
// Rij    =  exp(lambda * Sij)
// P      =  1 - e
// R'ijq  =  sum(k)[ Pk * Rik ]
// S'ijq  =  nearestInt[ ln(R'ijq) / lambda ]

// The numNormalLetters-th letter is considered to be a delimiter.
// When the letter in this sequence is a delimiter, or when the letter
// in the other sequence is abnormal, the quality data is ignored:
// S'ijq = Sij.

// Some letters may be considered to be "masked" versions of other
// letters (e.g. indicated by lowercase in the original sequence).  If
// masking is not applied, masked letters get the same scores as their
// unmasked versions.  If masking is applied: score <- min(score, 0).

#ifndef QUALITY_PSSM_MAKER_HH
#define QUALITY_PSSM_MAKER_HH

#include "ScoreMatrixRow.hh"

namespace cbrc {

typedef unsigned char uchar;

class QualityPssmMaker {
 public:
  void init(const ScoreMatrixRow *scoreMatrix,
            int numNormalLetters,  // 4 for DNA
            double lambda,  // scale factor for scoreMatrix
            bool isMatchMismatchMatrix,
            int matchScore,
            int mismatchScore,  // score not cost!
            int qualityOffset,  // typically 64
            const uchar* toUnmasked);  // maps letters to unmasked letters

  void make(const uchar *sequenceBeg,
            const uchar *sequenceEnd,
            const uchar *qualityBeg,
            int *pssmBeg,
            bool isApplyMasking) const;

 private:
  static const int qualityCapacity = 128;

  const ScoreMatrixRow *scoreMatrix;
  int numNormalLetters;
  double lambda;
  bool isMatchMismatchMatrix;
  const uchar* toUnmasked;

  double qualityToProbCorrect[qualityCapacity];
  int qualityToMatchScore[qualityCapacity];
  double expMatrix[scoreMatrixRowSize][scoreMatrixRowSize];
};

}

#endif
