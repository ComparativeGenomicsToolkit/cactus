// Copyright 2008 Martin C. Frith

// This struct holds alignments, and includes a procedure to
// non-redundantize alignments that share endpoints.

#ifndef ALIGNMENTPOT_HH
#define ALIGNMENTPOT_HH
#include "Alignment.hh"
#include <algorithm>  // sort
#include <vector>

namespace cbrc{

struct AlignmentPot{
  typedef std::vector<Alignment>::iterator iterator;

  // add an alignment to the pot
  void add( const Alignment& aln ) { items.push_back(aln); }

  // the number of alignments in the pot
  size_t size() const { return items.size(); }

  // erase any alignment that shares an endpoint with a higher-scoring
  // alignment
  void eraseSuboptimal();

  // sort the alignments in descending order of score
  void sort() { std::sort( items.begin(), items.end(), moreScore ); }

  // data:
  std::vector<Alignment> items;

  static bool moreScore( const Alignment& x, const Alignment& y ){
    // Try to break ties, so that alignments come in a consistent
    // order.  This makes it easier to compare different results.
    return x.score != y.score ? x.score > y.score : lessBeg( x, y );
  }

  static bool lessBeg( const Alignment& x, const Alignment& y ){
    return
      x.beg1() != y.beg1() ? x.beg1() < y.beg1() :
      x.beg2() != y.beg2() ? x.beg2() < y.beg2() :
      x.score > y.score;
  }

  static bool lessEnd( const Alignment& x, const Alignment& y ){
    return
      x.end1() != y.end1() ? x.end1() < y.end1() :
      x.end2() != y.end2() ? x.end2() < y.end2() :
      x.score > y.score;
  }

  static void mark( Alignment& a ) { a.blocks[0].score = -1; }

  static bool isMarked( Alignment& a ) { return a.blocks[0].score == -1; }
};

}  // end namespace cbrc
#endif  // ALIGNMENTPOT_HH
