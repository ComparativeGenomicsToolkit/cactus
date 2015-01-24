// Copyright 2010 Martin C. Frith

// This class calculates the expected number of alignments if we
// compare two completely random sequences.  It also works for two
// sets of sequences, e.g. two genomes.

#ifndef MCF_LOCAL_ALIGNMENT_EVALUER_HPP
#define MCF_LOCAL_ALIGNMENT_EVALUER_HPP

#include "sls_pvalues.hpp"

#include <vector>

namespace Mcf {

class LocalAlignmentEvaluer {
 public:
  typedef ncbi::blast::Sls::set_of_parameters GumbelParameters;

  LocalAlignmentEvaluer() { setBad(); }

  // Put us in a bad/undefined state.
  void setBad();

  // Are we in a bad/undefined state?
  bool isBad() const;

  // These routines might fail, putting us in the bad/undefined state.
  // letterProbs1 correspond to matrix rows (1st index), and
  // letterProbs2 correspond to matrix columns (2nd index).  The
  // letterProbs had better be normalized.

  void initGapless(const std::vector<double>& letterProbs1,
                   const std::vector<double>& letterProbs2,
                   const int *const *scoreMatrix);

  void initGapped(const std::vector<double>& letterProbs1,
                  const std::vector<double>& letterProbs2,
                  const int *const *scoreMatrix,
                  int gapExistCost, int gapExtendCost);

  // If we are in a bad/undefined state, the following functions have
  // undefined results!

  const GumbelParameters& gumbelParameters() const { return params; }

  // the expected number of alignments with score >= "score"
  double evalue(int score, double letterCount1, double letterCount2,
                double sequenceCount1, double sequenceCount2) /*const*/;

  // The minimum score with evalue <= maxEvalue
  // Never returns a score less than 1
  int minScore(double maxEvalue, double letterCount1, double letterCount2,
               double sequenceCount1, double sequenceCount2) /*const*/;

 private:
  GumbelParameters params;
};

}  // end namespace

#endif
