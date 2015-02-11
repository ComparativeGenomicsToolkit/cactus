// Copyright 2013 Martin C. Frith

#include "cbrc_linalg.hh"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cbrc {

void linalgSolve(double **m, double *v, unsigned s) {
  for (unsigned k = 0; k < s; ++k) {
    // partial pivoting:
    unsigned iMax = k;
    for (unsigned i = k; i < s; ++i)
      if (std::fabs(m[i][k]) > std::fabs(m[iMax][k]))
	iMax = i;
    if (iMax > k) {
      std::swap_ranges(m[k], m[k] + s, m[iMax]);
      std::swap(v[k], v[iMax]);
    }

    if (m[k][k] == 0.0)
      throw std::runtime_error("singular matrix");

    for (unsigned i = 0; i < s; ++i) {
      if (i == k) continue;
      double q = m[i][k] / m[k][k];
      for (unsigned j = k; j < s; ++j)
	m[i][j] -= q * m[k][j];
      v[i] -= q * v[k];
    }
  }

  for (unsigned k = 0; k < s; ++k)
    v[k] /= m[k][k];
}

}
