// Copyright 2010 Martin C. Frith

#include "mcf_local_alignment_evaluer.hpp"

#include "njn_localmaxstatmatrix.hpp"
#include "sls_alp_sim.hpp"

#include <cassert>
#include <cmath>

namespace Mcf {

void LocalAlignmentEvaluer::setBad() {
  params.lambda = -1;
}

bool LocalAlignmentEvaluer::isBad() const {
  return params.lambda < 0;
}

static void makeVector(std::vector<double>& vec, double val, double err) {
  vec.resize(2);
  vec[0] = val;
  vec[1] = val + err;
}

void LocalAlignmentEvaluer::initGapless(
    const std::vector<double>& letterProbs1,
    const std::vector<double>& letterProbs2,
    const int *const *scoreMatrix) {
  try {
    // don't limit the time, in order to get reproducible results:
    const double maxSeconds = 1e99;
    const double error = 1e-6;

    ncbi::blast::Njn::LocalMaxStatMatrix x(letterProbs1.size(), scoreMatrix,
                                           &letterProbs1[0], &letterProbs2[0],
                                           letterProbs2.size(), maxSeconds);

    if (x.getTerminated()) throw ncbi::blast::Sls::error("terminated", 1);

    params.lambda       = x.getLambda();
    params.lambda_error = error;
    params.C       = x.getC();
    params.C_error = error;
    params.K       = x.getK();
    params.K_error = error;

    params.a_I       = x.getA();
    params.a_I_error = error;
    params.a_J       = x.getA();
    params.a_J_error = error;

    params.sigma       = x.getAlpha();
    params.sigma_error = error;

    params.alpha_I       = x.getAlpha();
    params.alpha_I_error = error;
    params.alpha_J       = x.getAlpha();
    params.alpha_J_error = error;

    // a, a_error, alpha, alpha_error are never used?

    params.gapless_a       = x.getA();
    params.gapless_a_error = error;
    params.gapless_alpha       = x.getAlpha();
    params.gapless_alpha_error = error;

    params.G = 0;

    makeVector(params.m_LambdaSbs, params.lambda, error);
    makeVector(params.m_KSbs, params.K, error);
    makeVector(params.m_CSbs, params.C, error);
    makeVector(params.m_SigmaSbs, params.sigma, error);
    makeVector(params.m_AlphaISbs, params.alpha_I, error);
    makeVector(params.m_AlphaJSbs, params.alpha_J, error);
    makeVector(params.m_AISbs, params.a_I, error);
    makeVector(params.m_AJSbs, params.a_J, error);
  }
  catch (const ncbi::blast::Sls::error& e) {
    //std::cerr << "Sls::error: " << e.st << "\n";  // for debugging
    setBad();
  }
}

void LocalAlignmentEvaluer::initGapped(const std::vector<double>& letterProbs1,
                                       const std::vector<double>& letterProbs2,
                                       const int *const *scoreMatrix,
                                       int gapExistCost, int gapExtendCost) {
  initGapless(letterProbs1, letterProbs2, scoreMatrix);
  if (isBad()) return;

  try {
    const double lambdaTolerance = 0.01;  // ?
    const double kTolerance = 0.05;  // ?
    // don't limit the time, in order to get reproducible results:
    const double maxSeconds = 1e99;
    const double maxMegabytes = 200;
    const int randomSeed = 1;  // ?

    ncbi::blast::Sls::alp_data data(gapExistCost, gapExtendCost,
                                    lambdaTolerance, kTolerance,
                                    scoreMatrix, letterProbs1.size(),
                                    letterProbs1, letterProbs2,
                                    maxSeconds, maxMegabytes, randomSeed);

    ncbi::blast::Sls::alp_sim sim(&data);

    params.lambda       = sim.m_Lambda;
    params.lambda_error = sim.m_LambdaError;
    params.C       = sim.m_C;
    params.C_error = sim.m_CError;
    params.K       = sim.m_K;
    params.K_error = sim.m_KError;

    params.a_I       = sim.m_AI;
    params.a_I_error = sim.m_AIError;
    params.a_J       = sim.m_AJ;
    params.a_J_error = sim.m_AJError;

    params.sigma       = sim.m_Sigma;
    params.sigma_error = sim.m_SigmaError;

    params.alpha_I       = sim.m_AlphaI;
    params.alpha_I_error = sim.m_AlphaIError;
    params.alpha_J       = sim.m_AlphaJ;
    params.alpha_J_error = sim.m_AlphaJError;

    // a, a_error, alpha, alpha_error are never used?

    params.G = gapExistCost + gapExtendCost;

    params.m_LambdaSbs = sim.m_LambdaSbs;
    params.m_KSbs = sim.m_KSbs;
    params.m_CSbs = sim.m_CSbs;
    params.m_SigmaSbs = sim.m_SigmaSbs;
    params.m_AlphaISbs = sim.m_AlphaISbs;
    params.m_AlphaJSbs = sim.m_AlphaJSbs;
    params.m_AISbs = sim.m_AISbs;
    params.m_AJSbs = sim.m_AJSbs;
  }
  catch (const ncbi::blast::Sls::error& e) {
    //std::cerr << "Sls::error: " << e.st << "\n";  // for debugging
    setBad();
  }
}

double LocalAlignmentEvaluer::evalue(int score,
                                     double letterCount1,
                                     double letterCount2,
                                     double sequenceCount1,
                                     double sequenceCount2) /*const*/ {
  assert(!isBad());
  assert(sequenceCount1 > 0);
  assert(sequenceCount2 > 0);

  double averageLength1 = letterCount1 / sequenceCount1;
  double averageLength2 = letterCount2 / sequenceCount2;

  std::vector<double> out;
  std::vector<double> err;
  static ncbi::blast::Sls::pvalues x;
  x.calculate_P_values(score, score, averageLength1, averageLength2,
                       params, out, err);
  double evalue = out[0];
  return evalue * sequenceCount1 * sequenceCount2;
}

int LocalAlignmentEvaluer::minScore(double maxEvalue,
                                    double letterCount1,
                                    double letterCount2,
                                    double sequenceCount1,
                                    double sequenceCount2) /*const*/ {
  assert(!isBad());
  assert(maxEvalue >= 0);

  // first, get the minimum score without any edge correction
  double log_k_m_n = std::log(params.K * letterCount1 * letterCount2);
  double s = (log_k_m_n - std::log(maxEvalue)) / params.lambda;
  int score = static_cast<int>(std::ceil(s));
  if (score < 1) score = 1;
  // this shouldn't be necessary, but allow for floating-point weirdness:
  while (evalue(score, letterCount1, letterCount2,
                sequenceCount1, sequenceCount2) > maxEvalue) {
    ++score;
  }

  // now, take the edge correction into account
  do {
    --score;
    if (score == 0) break;
  } while (evalue(score, letterCount1, letterCount2,
                  sequenceCount1, sequenceCount2) <= maxEvalue);
  return score + 1;
}

}  // end namespace
