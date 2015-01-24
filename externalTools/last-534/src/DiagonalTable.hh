// Copyright 2008, 2010 Martin C. Frith

// This struct records coverage of "diagonals" by gapless alignments,
// when comparing two sequences.  The diagonal is the coordinate in
// one sequence minus the coordinate in the other sequence.  We assume
// that one sequence is scanned sequentially, and the other is
// accessed randomly.  We record the furthest sequential position
// covered so far in each diagonal.  This lets us avoid triggering
// gapless alignments in places that are already covered.

// Since the number of diagonals may be huge, we map them to a smaller
// number of bins.  To keep the bin populations small, when checking
// if a position is covered, we discard information about earlier
// sequential positions.

#ifndef DIAGONALTABLE_HH
#define DIAGONALTABLE_HH
#include <utility>  // pair
#include <vector>

namespace cbrc{

struct DiagonalTable{
  typedef unsigned indexT;
  typedef std::pair<indexT, indexT> pairT;

  enum { BINS = 256 };  // use a power-of-two for faster modulus (maybe)
                        // 256 is much faster than 65536 in my tests

  // is this position on this diagonal already covered by an alignment?
  bool isCovered( indexT sequentialPos, indexT randomPos );

  // add an alignment endpoint to the table:
  void addEndpoint( indexT sequentialPos, indexT randomPos );

  std::vector<pairT> hits[BINS];
};

}  // end namespace cbrc
#endif  // DIAGONALTABLE_HH
