// Copyright 2008, 2009, 2012 Martin C. Frith

// This struct holds parameters for so-called generalized affine gap
// costs (for pair-wise sequence alignment).  In this scheme, a "gap"
// may consist of unaligned regions in both sequences.  If these
// unaligned regions have sizes j and k, where j <= k, the cost is:

// a + b*(k-j) + c*j

// If c >= a + 2b, it reduces to standard affine gaps.  For more
// information, see: SF Altschul 1998 Proteins 32(1):88-96.

// In a further generalization, the costs for insertions and deletions
// may differ.  Typically, they would not differ:
// delExist = insExist = a
// delExtend = insExtend = b
// pairExtend = c

#ifndef GENERALIZEDAFFINEGAPCOSTS_HH
#define GENERALIZEDAFFINEGAPCOSTS_HH

namespace cbrc{

struct GeneralizedAffineGapCosts{
  int delExist;
  int delExtend;
  int insExist;
  int insExtend;
  int pairExtend;

  void assign( int a, int b, int A, int B, int c )
  { delExist = a; delExtend = b; insExist = A; insExtend = B; pairExtend = c; }

  bool isSymmetric() const
  { return insExist == delExist && insExtend == delExtend; }

  // Will standard affine gaps always suffice for maximal alignment scores?
  bool isAffine() const
  { return isSymmetric() && pairExtend >= delExist + 2 * delExtend; }

  // Return the score of a gap with the given sizes in a pair of
  // sequences, considering that it might be either one "generalized"
  // gap or two neighbouring "affine" gaps.
  // Here, gapSize2=0 would be a deletion, and gapSize1=0 would be an
  // insertion.
  int cost( int gapSize1, int gapSize2 ) const;
};

}

#endif
