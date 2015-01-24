// Copyright 2010 Martin C. Frith

// Functions that find gapless X-drop alignments between a sequence
// and a PSSM.

// These functions are analogous to those described in
// gaplessXdrop.hh.  Here, "pssm" replaces both "seq2" and "scorer".

#ifndef GAPLESS_PSSM_XDROP_HH
#define GAPLESS_PSSM_XDROP_HH

#include "ScoreMatrixRow.hh"

namespace cbrc {

typedef unsigned char uchar;

int forwardGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop);

int reverseGaplessPssmXdropScore(const uchar *seq,
                                 const ScoreMatrixRow *pssm,
                                 int maxScoreDrop);

const uchar *forwardGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score);

const uchar *reverseGaplessPssmXdropEnd(const uchar *seq,
                                        const ScoreMatrixRow *pssm,
                                        int score);

bool isOptimalGaplessPssmXdrop(const uchar *seq,
                               const uchar *seqEnd,
                               const ScoreMatrixRow *pssm,
                               int maxScoreDrop);

int gaplessPssmAlignmentScore(const uchar *seq,
                              const uchar *seqEnd,
                              const ScoreMatrixRow *pssm);

}

#endif
