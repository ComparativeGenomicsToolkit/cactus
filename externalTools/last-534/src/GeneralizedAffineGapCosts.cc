// Copyright 2008, 2009, 2012 Martin C. Frith

#include "GeneralizedAffineGapCosts.hh"
#include <algorithm>
#include <cassert>

namespace cbrc{

int GeneralizedAffineGapCosts::cost( int gapSize1, int gapSize2 ) const{
  if( gapSize1 == 0 && gapSize2 == 0 ) return 0;

  int c = delExist + delExtend * gapSize1 + insExist + insExtend * gapSize2;

  if( gapSize1 >= gapSize2 ){
    int d =
      delExist + delExtend * (gapSize1 - gapSize2) + pairExtend * gapSize2;
    assert( d >= delExist );  // try to catch overflow errors
    c = std::min( c, d );
  }

  if( gapSize2 >= gapSize1 ){
    int d =
      insExist + insExtend * (gapSize2 - gapSize1) + pairExtend * gapSize1;
    assert( d >= insExist );  // try to catch overflow errors
    c = std::min( c, d );
  }

  return c;
}

}
