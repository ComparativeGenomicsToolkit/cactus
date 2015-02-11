// Copyright 2011 Martin C. Frith

#include "GappedXdropAligner.hh"
#include "GappedXdropAlignerInl.hh"

namespace cbrc {

int GappedXdropAligner::align3pssm(const uchar *seq,
                                   const ScoreMatrixRow *pssmFrame0,
                                   const ScoreMatrixRow *pssmFrame1,
                                   const ScoreMatrixRow *pssmFrame2,
                                   bool isForward,
                                   int gapExistenceCost,
                                   int gapExtensionCost,
                                   int gapUnalignedCost,
                                   int frameshiftCost,
                                   int maxScoreDrop,
                                   int maxMatchScore) {
  bool isAffine = gapUnalignedCost >= gapExistenceCost + 2 * gapExtensionCost;

  std::size_t maxSeq1begs[] = { 9, 9, 0, 9, 9, 9, 9 };
  std::size_t minSeq1ends[] = { 0, 0, 1, 0, 0, 0, 0 };

  int bestScore = 0;

  init3();

  for (std::size_t antidiagonal = 7; /* noop */; ++antidiagonal) {
    std::size_t seq1beg = arrayMin(maxSeq1begs);
    std::size_t seq1end = arrayMax(minSeq1ends);

    if (seq1beg >= seq1end) break;

    std::size_t scoreEnd = scoreEnds.back();
    std::size_t numCells = seq1end - seq1beg;

    initAntidiagonal3(seq1beg, scoreEnd, numCells);

    const ScoreMatrixRow *pssm =
        whichFrame(antidiagonal, pssmFrame0, pssmFrame1, pssmFrame2);

    std::size_t seq2pos = (antidiagonal - 7) / 3 - seq1beg;

    const uchar *s1 = isForward ? seq + seq1beg : seq - seq1beg - 1;
    const ScoreMatrixRow *s2 = isForward ? pssm + seq2pos : pssm - seq2pos - 1;

    if (isDelimiter(0, *s2)) {
      // prevent forward frameshifts from jumping over delimiters:
      if (maxSeq1begs[1] == seq1beg) ++maxSeq1begs[1];
      // Update maxScoreDrop in some clever way?
      // But be careful if the -1 frame starts in an initial delimiter.
    }

    int minScore = bestScore - maxScoreDrop;

    int *x0 = &xScores[scoreEnd];
    int *y0 = &yScores[scoreEnd];
    int *z0 = &zScores[scoreEnd];
    const int *y3 = &yScores[hori3(antidiagonal, seq1beg)];
    const int *z3 = &zScores[vert3(antidiagonal, seq1beg)];
    const int *x6 = &xScores[diag3(antidiagonal, seq1beg)];
    const int *x5 = &xScores[diag3(antidiagonal + 1, seq1beg)];
    const int *x7 = &xScores[diag3(antidiagonal - 1, seq1beg)];

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0last = x0 + numCells;

    *x0++ = *y0++ = *z0++ = -INF;  // add one pad cell

    const int *x0base = x0 - seq1beg;

    if (isAffine) {
      if (isForward)
        while (1) {
          int s = maxValue(*x5, *x7);
          int x = maxValue(*x6, s - frameshiftCost);
          int y = *y3 - gapExtensionCost;
          int z = *z3 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + (*s2)[*s1];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          ++s1;  --s2;  ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;
        }
      else
        while (1) {
          int s = maxValue(*x5, *x7);
          int x = maxValue(*x6, s - frameshiftCost);
          int y = *y3 - gapExtensionCost;
          int z = *z3 - gapExtensionCost;
          int b = maxValue(x, y, z);
          if (b >= minScore) {
            updateBest(bestScore, b, antidiagonal, x0, x0base);
            *x0 = b + (*s2)[*s1];
            int g = b - gapExistenceCost;
            *y0 = maxValue(g, y);
            *z0 = maxValue(g, z);
          }
          else *x0 = *y0 = *z0 = -INF;
          if (x0 == x0last) break;
          --s1;  ++s2;  ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;
        }
    } else {
      const int *y6 = &yScores[diag3(antidiagonal, seq1beg)];
      const int *z6 = &zScores[diag3(antidiagonal, seq1beg)];
      while (1) {
        int s = maxValue(*x5, *x7);
        int x = maxValue(*x6, s - frameshiftCost);
        int y = maxValue(*y3 - gapExtensionCost, *y6 - gapUnalignedCost);
        int z = maxValue(*z3 - gapExtensionCost, *z6 - gapUnalignedCost);
        int b = maxValue(x, y, z);
        if (b >= minScore) {
          updateBest(bestScore, b, antidiagonal, x0, x0base);
          *x0 = b + (*s2)[*s1];
          int g = b - gapExistenceCost;
          *y0 = maxValue(g, y);
          *z0 = maxValue(g, z);
        }
        else *x0 = *y0 = *z0 = -INF;
        if (x0 == x0last) break;
        ++x0;  ++y0;  ++z0;  ++y3;  ++z3;  ++x5;  ++x6;  ++x7;  ++y6;  ++z6;
        if (isForward) { ++s1;  --s2; }
        else           { --s1;  ++s2; }
      }
    }

    if (isDelimiter(*s1, *pssmFrame2) && isDelimiter(*s1, *pssmFrame0))
      updateMaxScoreDrop(maxScoreDrop, numCells, maxMatchScore);

    updateFiniteEdges3(maxSeq1begs, minSeq1ends, x0base, x0 + 1, numCells);
  }

  return bestScore;
}

}
