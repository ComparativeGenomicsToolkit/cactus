// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

// This routine reads alignments of DNA reads to a genome, and
// estimates the probability that each alignment represents the
// genomic source of the read.  It assumes that the reads come in
// pairs, where each pair is from either end of a DNA fragment.

// The "rna" option makes it assume that the genomic fragment lengths
// follow a log-normal distribution (instead of a normal
// distribution).  In one test with human RNA, log-normal was a
// remarkably good fit, but not perfect.  The true distribution looked
// like a mixture of 2 log-normals: a dominant one for shorter
// introns, and a minor one for huge introns.  Thus, our use of a
// single log-normal fails to model rare, huge introns.  To compensate
// for that, the default value of "disjoint" can be increased when
// "rna" is used.

// (Should we try to estimate the prior probability of disjoint
// mapping from the data?  But maybe ignore low-scoring alignments for
// that?  Estimate disjoint maps to opposite strands of same
// chromosome = maps to same strand of same chromosome?)

#ifndef LAST_PAIR_PROBS_HH
#define LAST_PAIR_PROBS_HH

#include <string>
#include <vector>
#include <set>

struct LastPairProbsOptions {
  bool rna;
  bool estdist;
  double mismap;
  bool isFraglen;
  double fraglen;
  bool isSdev;
  double sdev;
  bool isDisjoint;
  double disjoint;
  std::set<std::string> circular;
  std::vector<std::string> inputFileNames;
  double outer;
  double inner;
  double disjointScore;
  double maxMissingScore1;
  double maxMissingScore2;
};

void lastPairProbs(LastPairProbsOptions& opts);

#endif
