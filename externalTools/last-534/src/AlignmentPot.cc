// Copyright 2008, 2010 Martin C. Frith

#include "AlignmentPot.hh"

namespace cbrc{

void AlignmentPot::eraseSuboptimal(){
  if( items.empty() ) return;
  int maxScore;

  std::sort( items.begin(), items.end(), lessBeg );
  maxScore = items[0].score;

  for( iterator i = items.begin() + 1; i < items.end(); ++i ){
    if( i->beg1() == (i-1)->beg1() && i->beg2() == (i-1)->beg2() ){
      if( i->score < maxScore ) mark( *i );
    }
    else maxScore = i->score;
  }

  std::sort( items.begin(), items.end(), lessEnd );
  maxScore = items[0].score;

  for( iterator i = items.begin() + 1; i < items.end(); ++i ){
    if( i->end1() == (i-1)->end1() && i->end2() == (i-1)->end2() ){
      if( i->score < maxScore ) mark( *i );
    }
    else maxScore = i->score;
  }

  items.erase( std::remove_if( items.begin(), items.end(), isMarked ),
	       items.end() );
}

}  // end namespace cbrc
