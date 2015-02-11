// Copyright 2011, 2013 Martin C. Frith

#ifndef QUALITY_SCORE_UTIL_HH
#define QUALITY_SCORE_UTIL_HH

#include "stringify.hh"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace cbrc {

typedef unsigned char uchar;

inline int nearestInt(double x) {
  if (x > 0) return  static_cast<int>(0.5 + x);
  else       return -static_cast<int>(0.5 - x);
}

inline double phredErrorProb(int qualityScore) {
  if (qualityScore < 0) qualityScore = 0;
  return std::pow(10.0, -0.1 * qualityScore);
}

inline double solexaErrorProb(int qualityScore) {
  return 1 / (1 + std::pow(10.0, 0.1 * qualityScore));
}

inline double qualityUncertainty(int qualityCode, int qualityOffset,
                                 bool isPhred, double letterProb) {
  int q = qualityCode - qualityOffset;
  double errorProb = isPhred ? phredErrorProb(q) : solexaErrorProb(q);
  // The next line of code is a kludge to avoid numerical instability.
  // An error probability of 1 is bizarre and probably shouldn't occur anyway.
  if (errorProb >= 1) errorProb = 0.999999;
  double otherProb = 1 - letterProb;
  assert(letterProb >= 0);
  assert(otherProb > 0);
  return errorProb / otherProb;
}

inline int qualityPairScore(double expScore, double uncertainty1,
                            double uncertainty2, double lambda) {
  double x = (1 - uncertainty1) * (1 - uncertainty2) * (expScore - 1) + 1;
  assert(lambda > 0);
  assert(x > 0);
  return nearestInt(std::log(x) / lambda);
}

inline void checkQualityCodes(const uchar *qualityBeg,
                              const uchar *qualityEnd, uchar minQuality) {
  while (qualityBeg < qualityEnd) {
    uchar q = *qualityBeg++;
    if (q < minQuality) {
      unsigned x = q;
      unsigned m = minQuality;
      throw std::runtime_error("invalid quality code '" + stringify(q) +
                               "' (ASCII " + stringify(x) +
                               "), minimum is '" + stringify(minQuality) +
                               "' (ASCII " + stringify(m) + ")");
    }
  }
}

}

#endif
