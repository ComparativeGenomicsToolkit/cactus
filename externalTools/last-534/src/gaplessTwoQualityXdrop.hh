// Copyright 2011 Martin C. Frith

// Functions that find gapless X-drop alignments between two sequences
// that both have quality scores.

// These functions are analogous to those described in
// gaplessXdrop.hh.

#ifndef GAPLESS_TWO_QUALITY_XDROP_HH
#define GAPLESS_TWO_QUALITY_XDROP_HH

namespace cbrc {

typedef unsigned char uchar;

class TwoQualityScoreMatrix;

int forwardGaplessTwoQualityXdropScore(const uchar *seq1,
                                       const uchar *qual1,
                                       const uchar *seq2,
                                       const uchar *qual2,
                                       const TwoQualityScoreMatrix &m,
                                       int maxScoreDrop);

int reverseGaplessTwoQualityXdropScore(const uchar *seq1,
                                       const uchar *qual1,
                                       const uchar *seq2,
                                       const uchar *qual2,
                                       const TwoQualityScoreMatrix &m,
                                       int maxScoreDrop);

const uchar *forwardGaplessTwoQualityXdropEnd(const uchar *seq1,
                                              const uchar *qual1,
                                              const uchar *seq2,
                                              const uchar *qual2,
                                              const TwoQualityScoreMatrix &m,
                                              int score);

const uchar *reverseGaplessTwoQualityXdropEnd(const uchar *seq1,
                                              const uchar *qual1,
                                              const uchar *seq2,
                                              const uchar *qual2,
                                              const TwoQualityScoreMatrix &m,
                                              int score);

bool isOptimalGaplessTwoQualityXdrop(const uchar *seq1,
                                     const uchar *seq1end,
                                     const uchar *qual1,
                                     const uchar *seq2,
                                     const uchar *qual2,
                                     const TwoQualityScoreMatrix &m,
                                     int maxScoreDrop);

int gaplessTwoQualityAlignmentScore(const uchar *seq1,
                                    const uchar *seq1end,
                                    const uchar *qual1,
                                    const uchar *seq2,
                                    const uchar *qual2,
                                    const TwoQualityScoreMatrix &m);

}

#endif
