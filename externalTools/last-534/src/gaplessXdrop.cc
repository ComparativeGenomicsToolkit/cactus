// Copyright 2010, 2013 Martin C. Frith

#include "gaplessXdrop.hh"
#include <stdexcept>

static void err(const char *s) { throw std::overflow_error(s); }

namespace cbrc {

int forwardGaplessXdropScore(const uchar *seq1,
                             const uchar *seq2,
                             const ScoreMatrixRow *scorer,
                             int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += scorer[*seq1++][*seq2++];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in forward gapless extension");
  return score;
}

int reverseGaplessXdropScore(const uchar *seq1,
                             const uchar *seq2,
                             const ScoreMatrixRow *scorer,
                             int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += scorer[*--seq1][*--seq2];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in reverse gapless extension");
  return score;
}

const uchar *forwardGaplessXdropEnd(const uchar *seq1,
                                    const uchar *seq2,
                                    const ScoreMatrixRow *scorer,
                                    int score) {
  int s = 0;
  while (s < score) s += scorer[*seq1++][*seq2++];
  return seq1;
}

const uchar *reverseGaplessXdropEnd(const uchar *seq1,
                                    const uchar *seq2,
                                    const ScoreMatrixRow *scorer,
                                    int score) {
  int s = 0;
  while (s < score) s += scorer[*--seq1][*--seq2];
  return seq1;
}

bool isOptimalGaplessXdrop(const uchar *seq1,
                           const uchar *seq1end,
                           const uchar *seq2,
                           const ScoreMatrixRow *scorer,
                           int maxScoreDrop) {
  int score = 0;
  int maxScore = 0;
  while (seq1 < seq1end) {
    score += scorer[*seq1++][*seq2++];
    if (score > maxScore) maxScore = score;
    else if (score <= 0 ||                       // non-optimal prefix
             seq1 == seq1end ||                  // non-optimal suffix
             score < maxScore - maxScoreDrop) {  // excessive score drop
      return false;
    }
  }
  return true;
}

int gaplessAlignmentScore(const uchar *seq1,
                          const uchar *seq1end,
                          const uchar *seq2,
                          const ScoreMatrixRow *scorer) {
  int score = 0;
  while (seq1 < seq1end) score += scorer[*seq1++][*seq2++];
  return score;
}

}
