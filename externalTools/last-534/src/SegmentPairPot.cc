// Copyright 2008, 2010, 2014 Martin C. Frith

#include "SegmentPairPot.hh"
#include <cassert>

// Check if n1/d1 < n2/d2, without overflow.
// This uses "continued fractions".
static bool lessFraction( size_t n1, size_t d1, size_t n2, size_t d2 ) {
  // assume that d1 > 0 and d2 > 0
  size_t q1 = n1 / d1;
  size_t q2 = n2 / d2;
  if( q1 < q2 ) return true;
  if( q1 > q2 ) return false;
  size_t r1 = n1 % d1;
  size_t r2 = n2 % d2;
  if( r2 == 0 ) return false;
  if( r1 == 0 ) return true;
  return lessFraction( d2, r2, d1, r1 );
}

namespace cbrc{

static bool cullingLess( const SegmentPair& x, const SegmentPair& y ){
  return x.beg2() != y.beg2() ? x.beg2() < y.beg2()
    :    x.score  != y.score  ? x.score  > y.score
    :                           x.beg1() < y.beg1();
}

void SegmentPairPot::cull( size_t limit ){
  if( !limit ) return;

  std::sort( items.begin(), items.end(), cullingLess );

  items.erase( unique( items.begin(), items.end() ), items.end() );
  // this is redundantly repeated in SegmentPairPot::sort()

  std::vector<iterator> stash;

  for( iterator i = items.begin(); i < items.end(); ++i ){
    size_t iBeg = i->beg2();
    size_t iEnd = i->end2();
    size_t iLen = i->size;
    int iScore = i->score;

    size_t numOfDominatingItems = 0;

    size_t x = 0;

    for( size_t y = 0; y < stash.size(); ++y ){
      iterator j = stash[ y ];
      size_t jEnd = j->end2();
      if( jEnd <= iBeg ) continue;
      stash[ x++ ] = j;
      if( jEnd < iEnd ) continue;
      size_t jLen = j->size;
      int jScore = j->score;
      if( lessFraction( iScore, iLen, jScore, jLen ) ) ++numOfDominatingItems;
    }

    stash.resize( x );

    if( numOfDominatingItems >= limit ){
      mark( *i );
    }else{
      stash.push_back( i );
    }
  }

  erase_if( items, isMarked );
}

void SegmentPairPot::sort(){
  std::sort( items.begin(), items.end(), itemLess );

  // Remove duplicates.  We assume that, after sorting, duplicates are
  // consecutive.  This will be true if non-duplicates never overlap,
  // which is true if all the SegmentPairs are "optimal".
  items.erase( unique( items.begin(), items.end() ), items.end() );

  iters.clear();
  iters.reserve( items.size() );

  for( iterator i = items.begin(); i < items.end(); ++i ){
    iters.push_back(i);
  }

  std::sort( iters.begin(), iters.end(), iterLess );
}

void SegmentPairPot::markOverlaps( const SegmentPair& sp ){
  iterator i = std::lower_bound( items.begin(), items.end(), sp, itemLess );

  // Assume we need to check just one previous item.  This assumption
  // will be true provided that the items never overlap each other.
  if( i > items.begin() &&
      (i-1)->diagonal() == sp.diagonal() &&
      (i-1)->end1() > sp.beg1() ){
    mark( *(i-1) );
  }

  while( i < items.end() &&
	 i->diagonal() == sp.diagonal() &&
	 i->beg1() < sp.end1() ){
    mark( *i );
    ++i;
  }
}

void SegmentPairPot::markAllOverlaps( const std::vector<SegmentPair>& sps ){
  for( const_iterator i = sps.begin(); i < sps.end(); ++i ){
    markOverlaps( *i );
  }
}

void SegmentPairPot::markTandemRepeats( const SegmentPair& sp,
					indexT maxDistance ){
  assert( !items.empty() );

  // Careful: if we are self-comparing lots of short sequences, there
  // may be many items on the same diagonal as sp.

  SegmentPair nextDiagonal( sp.beg1() + 1, sp.beg2(), 0, 0 );
  iterator n = std::lower_bound( items.begin(), items.end(),
				 nextDiagonal, itemLess );
  if( n == items.end() )  n = items.begin();

  iterator j = n;
  do{  // funny loop to deal with wrap-around
    if( j->diagonal() - sp.diagonal() > maxDistance )  break;
    if( j->beg2() >= sp.beg2() && j->end1() <= sp.end1() )  mark( *j );
    ++j;
    if( j == items.end() )  j = items.begin();
  }while( j != n );

  SegmentPair prevDiagonal( sp.end1() - 1, sp.end2(), 0, 0 );
  iterator p = std::lower_bound( items.begin(), items.end(),
				 prevDiagonal, itemLess );
  if( p == items.end() )  p = items.begin();

  iterator k = p;
  do{  // funny loop to deal with wrap-around
    if( k == items.begin() )  k = items.end();
    --k;
    if( sp.diagonal() - k->diagonal() > maxDistance )  break;
    if( k->beg1() >= sp.beg1() && k->end2() <= sp.end2() )  mark( *k );
  }while( k != p );
}

}
