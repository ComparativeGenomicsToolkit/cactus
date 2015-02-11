// Copyright 2013 Martin C. Frith

// This routine solves simultaneous linear equations, such as:
//
//   m11 * x1  +  m12 * x2  +  m13 * x3   =   v1
//   m21 * x1  +  m22 * x2  +  m23 * x3   =   v2
//   m31 * x1  +  m32 * x2  +  m33 * x3   =   v3
//
// We know m and v, and we wish to determine x.
//
// Input: m should be a matrix of size s*s, and v should be a vector
// of length s.
//
// Output: the result is written into v.  The calculation modifies m.
//
// If the matrix is singular, an error is thrown.

#ifndef CBRC_LINALG_HH
#define CBRC_LINALG_HH

namespace cbrc {

void linalgSolve(double **m, double *v, unsigned s);

}

#endif
