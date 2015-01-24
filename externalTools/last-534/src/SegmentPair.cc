// Copyright 2008, 2009, 2011, 2012 Martin C. Frith

#include "SegmentPair.hh"

namespace cbrc{

void SegmentPair::maxIdenticalRun( const uchar* seq1, const uchar* seq2,
				   const uchar* canonical ){
  const uchar* s1 = seq1 + beg1();
  const uchar* s2 = seq2 + beg2();
  const uchar* e1 = seq1 + end1();

  const uchar* bestEnd1 = s1;
  indexT bestSize = 0;
  indexT runSize = 0;

  while( s1 < e1 ){
    if( canonical[ *s1++ ] == canonical[ *s2++ ] ){
      ++runSize;
      if( runSize > bestSize ){
	bestSize = runSize;
	bestEnd1 = s1;
      }
    }
    else runSize = 0;
  }

  const uchar* bestBeg1 = bestEnd1 - bestSize;
  indexT offset = bestBeg1 - seq1 - beg1();
  start1 += offset;
  start2 += offset;
  size = bestSize;
}

}
