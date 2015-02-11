// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

// This struct holds segment-pairs, and allows us to find
// (near-)overlaps between sets of segment-pairs efficiently.  To find
// overlaps, we use binary search, after sorting the segment-pairs by
// position.  To mark a segment-pair that is overlapped, we set its
// score to zero.

#ifndef SEGMENT_PAIR_POT_HH
#define SEGMENT_PAIR_POT_HH
#include "SegmentPair.hh"
#include <algorithm>  // remove_if
#include <stddef.h>  // size_t
#include <vector>

namespace cbrc{

template< typename Container, typename Predicate >
void erase_if( Container& c, Predicate p ){
  c.erase( std::remove_if( c.begin(), c.end(), p ), c.end() );
}

struct SegmentPairPot{
  typedef SegmentPair::indexT indexT;
  typedef std::vector<SegmentPair>::iterator iterator;
  typedef std::vector<SegmentPair>::const_iterator const_iterator;

  // add a SegmentPair to the pot
  void add( const SegmentPair& sp ) { items.push_back(sp); }

  // the number of SegmentPairs in the pot
  size_t size() const { return items.size(); }

  // remove duplicate SegmentPairs, and then:
  // discard SegmentPairs whose query coordinates lie within limit or
  // more SegmentPairs that have higher score-per-length
  void cull( size_t limit );
  // does nothing if limit is 0

  // this must be called before using the following methods
  void sort();

  // get the i-th SegmentPair, sorted by score
  SegmentPair& get( size_t i ) { return *iters[i]; }

  // set the score of all items that overlap sp to zero
  void markOverlaps( const SegmentPair& sp );

  // set the score of all items that overlap anything in sps to zero
  void markAllOverlaps( const std::vector<SegmentPair>& sps );

  // mark (near-)tandem repeats within sp
  // to avoid death by dynamic programming when self-aligning a large sequence
  void markTandemRepeats( const SegmentPair& sp, indexT maxDistance );

  // data:
  std::vector<SegmentPair> items;
  std::vector<iterator> iters;

  // sort criterion for sorting by position
  static bool itemLess( const SegmentPair& x, const SegmentPair& y ){
    return x.diagonal() != y.diagonal() ? x.diagonal() < y.diagonal()
      : x.beg1() < y.beg1();
  }

  // sort criterion for sorting by score (in descending order)
  static bool iterLess( const iterator& x, const iterator& y ){
    // break ties in an arbitrary way, to make the results more reproducible:
    return (x->score  != y->score ) ? (x->score  > y->score )
      :    (x->start1 != y->start1) ? (x->start1 < y->start1)
      :                               (x->start2 < y->start2);
  }

  static void mark( SegmentPair& s ) { s.score = 0; }

  static bool isMarked( const SegmentPair& s ) { return s.score == 0; }

  // It turns out that we need a second kind of mark also:

  static void markAsGood( SegmentPair& s ) { s.size = 0; }

  static bool isNotMarkedAsGood( const SegmentPair& s ) { return s.size != 0; }
};

}

#endif
