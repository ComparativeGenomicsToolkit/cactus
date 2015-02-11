// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

#include "SubsetSuffixArray.hh"
//#include <iostream>  // for debugging

using namespace cbrc;

// use past results to speed up long matches?
// could & probably should return the match depth
void SubsetSuffixArray::match( const indexT*& beg, const indexT*& end,
                               const uchar* queryPtr, const uchar* text,
                               indexT maxHits,
                               indexT minDepth, indexT maxDepth ) const{
  // the next line is unnecessary, but makes it faster in some cases:
  if( maxHits == 0 && minDepth < maxDepth ) minDepth = maxDepth;

  indexT depth = 0;
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  indexT startDepth = std::min( bucketDepth, minDepth );
  const indexT* bucketPtr = &buckets[0];

  while( depth < startDepth ){
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ){
      beg = end = &index[0];
      return;
    }
    ++depth;
    bucketPtr += subset * bucketSteps[depth];
    subsetMap = seed.nextMap( subsetMap );
  }

  indexT bucketBeg = *bucketPtr;
  indexT bucketEnd = depth ? *(bucketPtr + bucketSteps[depth]) : index.size();

  while( depth < bucketDepth ){
    if( bucketEnd - bucketBeg <= maxHits || depth >= maxDepth ){
      beg = &index[0] + bucketBeg;
      end = &index[0] + bucketEnd;
      return;
    }
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ){
      beg = end = &index[0];
      return;
    }
    ++depth;
    indexT step = bucketSteps[depth];
    bucketPtr += subset * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
    subsetMap = seed.nextMap( subsetMap );
  }

  // match using binary search:
  beg = &index[0] + bucketBeg;
  end = &index[0] + bucketEnd;

  if( depth < minDepth ){
    indexT d = depth;
    const uchar* s = subsetMap;
    while( depth < minDepth ){
      uchar subset = subsetMap[ queryPtr[depth] ];
      if( subset == CyclicSubsetSeed::DELIMITER ){
	beg = end;
	return;
      }
      ++depth;
      subsetMap = seed.nextMap( subsetMap );
    }
    equalRange2( beg, end, queryPtr + d, queryPtr + depth, text + d, s );
  }

  while( true ){
    if( indexT(end - beg) <= maxHits || depth >= maxDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ){
      beg = end;
      return;
    }
    equalRange( beg, end, text+depth, subsetMap, subset );
    ++depth;
    subsetMap = seed.nextMap( subsetMap );
  }
}

void SubsetSuffixArray::countMatches( std::vector<unsigned long long>& counts,
				      const uchar* queryPtr,
				      const uchar* text,
				      indexT maxDepth ) const{
  indexT depth = 0;
  const uchar* subsetMap = seed.firstMap();

  // match using buckets:
  indexT bucketDepth = maxBucketPrefix();
  const indexT* bucketPtr = &buckets[0];
  indexT bucketBeg = 0;
  indexT bucketEnd = index.size();

  while( depth < bucketDepth ){
    if( bucketBeg == bucketEnd ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += bucketEnd - bucketBeg;
    if( depth >= maxDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    ++depth;
    indexT step = bucketSteps[depth];
    bucketPtr += subset * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
    subsetMap = seed.nextMap( subsetMap );
  }

  // match using binary search:
  const indexT* beg = &index[0] + bucketBeg;
  const indexT* end = &index[0] + bucketEnd;

  while( true ){
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    if( depth >= maxDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    equalRange( beg, end, text+depth, subsetMap, subset );
    ++depth;
    subsetMap = seed.nextMap( subsetMap );
  }
}

void SubsetSuffixArray::equalRange( const indexT*& beg, const indexT*& end,
				    const uchar* textBase,
				    const uchar* subsetMap, uchar subset ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    uchar s = subsetMap[ textBase[ *mid ] ];
    if( s < subset ){
      beg = mid + 1;
    }else if( s > subset ){
      end = mid;
    }else{
      if( subset > 0 )  // this "if" is unnecessary, but makes it a bit faster
	beg = lowerBound( beg, mid, textBase, subsetMap, subset );
      end = upperBound( mid + 1, end, textBase, subsetMap, subset );
      return;
    }
  }
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::lowerBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    if( subsetMap[ textBase[ *mid ] ] < subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }
  return beg;
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::upperBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    if( subsetMap[ textBase[ *mid ] ] <= subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }
  return end;
}

void SubsetSuffixArray::equalRange2( const indexT*& beg, const indexT*& end,
				     const uchar* queryBeg,
				     const uchar* queryEnd,
				     const uchar* textBase,
				     const uchar* subsetMap ) const{
  const uchar* qBeg = queryBeg;
  const uchar* qEnd = qBeg;
  const uchar* tBeg = textBase;
  const uchar* tEnd = tBeg;
  const uchar* sBeg = subsetMap;
  const uchar* sEnd = sBeg;

  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    indexT offset = *mid;
    const uchar* q;
    const uchar* t;
    const uchar* s;
    if( qBeg < qEnd ){
      q = qBeg;
      t = tBeg + offset;
      s = sBeg;
    }else{
      q = qEnd;
      t = tEnd + offset;
      s = sEnd;
    }
    uchar x, y;
    for( ; ; ){  // loop over consecutive letters
      x = s[ *t ];  // this text letter's subset
      y = s[ *q ];  // this query letter's subset
      if( x != y ) break;
      ++q;  // next query letter
      if( q == queryEnd ){  // we found a full match to [queryBeg, queryEnd)
	beg = lowerBound2( beg, mid, qBeg, queryEnd, tBeg, sBeg );
	end = upperBound2( mid + 1, end, qEnd, queryEnd, tEnd, sEnd );
	return;
      }
      ++t;  // next text letter
      s = seed.nextMap( s );  // next mapping from letters to subsets
    }
    if( x < y ){
      beg = mid + 1;
      // the next 3 lines are unnecessary, but make it faster:
      qBeg = q;
      tBeg = t - offset;
      sBeg = s;
    }else{
      end = mid;
      // the next 3 lines are unnecessary, but make it faster:
      qEnd = q;
      tEnd = t - offset;
      sEnd = s;
    }
  }
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::lowerBound2( const indexT* beg, const indexT* end,
				const uchar* queryBeg, const uchar* queryEnd,
				const uchar* textBase,
				const uchar* subsetMap ) const{
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    indexT offset = *mid;
    const uchar* t = textBase + offset;
    const uchar* q = queryBeg;
    const uchar* s = subsetMap;
    for( ; ; ){  // loop over consecutive letters
      if( s[ *t ] < s[ *q ] ){
	beg = mid + 1;
	// the next 3 lines are unnecessary, but make it faster:
	queryBeg = q;
	textBase = t - offset;
	subsetMap = s;
	break;
      }
      ++q;  // next query letter
      if( q == queryEnd ){  // we found a full match to [queryBeg, queryEnd)
	end = mid;
	break;
      }
      ++t;  // next text letter
      s = seed.nextMap( s );  // next mapping from letters to subsets
    }
  }
  return beg;
}

const SubsetSuffixArray::indexT*
SubsetSuffixArray::upperBound2( const indexT* beg, const indexT* end,
				const uchar* queryBeg, const uchar* queryEnd,
				const uchar* textBase,
				const uchar* subsetMap ) const{
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    indexT offset = *mid;
    const uchar* t = textBase + offset;
    const uchar* q = queryBeg;
    const uchar* s = subsetMap;
    for( ; ; ){  // loop over consecutive letters
      if( s[ *t ] > s[ *q ] ){
	end = mid;
	// the next 3 lines are unnecessary, but make it faster:
	queryBeg = q;
	textBase = t - offset;
	subsetMap = s;
        break;
      }
      ++q;  // next query letter
      if( q == queryEnd ){  // we found a full match to [queryBeg, queryEnd)
	beg = mid + 1;
	break;
      }
      ++t;  // next text letter
      s = seed.nextMap( s );  // next mapping from letters to subsets
    }
  }
  return end;
}
