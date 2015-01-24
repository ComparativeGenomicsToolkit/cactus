// Copyright 2010, 2013 Martin C. Frith

#include "gaplessPssmXdrop.hh"
#include <stdexcept>

static void err(const char *s) { throw std::overflow_error(s); }

namespace cbrc {

int forwardGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*pssm++)[*seq++];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in forward gapless extension with PSSM");
  return score;
}

int reverseGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop) {
  int score = 0;
  int s = 0;
  while (true) {
    s += (*--pssm)[*--seq];  // overflow risk
    if (s < score - maxScoreDrop) break;
    if (s > score) score = s;
  }
  if (score - s < 0)
    err("score overflow in reverse gapless extension with PSSM");
  return score;
}

const uchar *forwardGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score) {
  int s = 0;
  while (s < score) s += (*pssm++)[*seq++];
  return seq;
}

const uchar *reverseGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score) {
  int s = 0;
  while (s < score) s += (*--pssm)[*--seq];
  return seq;
}

bool isOptimalGaplessPssmXdrop(const uchar *seq,
                               const uchar *seqEnd,
                               const ScoreMatrixRow *pssm,
                               int maxScoreDrop) {
  int score = 0;
  int maxScore = 0;
  while (seq < seqEnd) {
    score += (*pssm++)[*seq++];
    if (score > maxScore) maxScore = score;
    else if (score <= 0 ||                       // non-optimal prefix
             seq == seqEnd ||                    // non-optimal suffix
             score < maxScore - maxScoreDrop) {  // excessive score drop
      return false;
    }
  }
  return true;
}

int gaplessPssmAlignmentScore(const uchar *seq,
                              const uchar *seqEnd,
                              const ScoreMatrixRow *pssm) {
  int score = 0;
  while (seq < seqEnd) score += (*pssm++)[*seq++];
  return score;
}

}
