// Copyright 2008, 2009, 2010, 2011, 2013 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SubsetSuffixArray.hh"
#include <algorithm>  // iter_swap, min

using namespace cbrc;

namespace{
  typedef SubsetSuffixArray::indexT indexT;
  struct Stack{ indexT* beg; indexT* end; indexT depth; };
  Stack stack[1048576];  // big enough???
  Stack* sp = stack;
}

#define PUSH(b, e, d) sp->beg = b, sp->end = e, (sp++)->depth = d
#define  POP(b, e, d) b = (--sp)->beg, e = sp->end, d = sp->depth

static void insertionSort( const uchar* text, const CyclicSubsetSeed& seed,
			   indexT* beg, indexT* end, indexT depth ){
  const uchar* textBase = text + depth;
  const uchar* subsetMap = seed.subsetMap(depth);

  for( indexT* i = beg+1; i < end; ++i ){
    for( indexT* j = i; j > beg; --j ){
      const uchar* s = textBase + *(j-1);
      const uchar* t = textBase + *j;
      const uchar* m = subsetMap;

      while( m[ *s ] == m[ *t ] && m[ *s ] < CyclicSubsetSeed::DELIMITER ){
        ++s;
        ++t;
        m = seed.nextMap(m);
      }

      if( m[ *s ] <= m[ *t ] ) break;
      std::iter_swap( j, j-1 );
    }
  }
}

// Specialized sort for 1 symbol + 1 delimiter.
// E.g. wildcard positions in spaced seeds.
static void radixSort1( const uchar* text, const uchar* subsetMap,
			indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's

  while( end0 < end ){
    const indexT x = *end0;
    switch( subsetMap[ text[x] ] ){
    case 0:
      end0++;
      break;
    default:  // the delimiter subset
      *end0 = *--end;
      *end = x;
      break;
    }
  }

  PUSH( beg, end, depth );   // the '0's
}

// Specialized sort for 2 symbols + 1 delimiter.
// E.g. transition-constrained positions in subset seeds.
static void radixSort2( const uchar* text, const uchar* subsetMap,
                        indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's

  while( end1 < end ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      default:  // the delimiter subset
        *end1 = *--end;
        *end = x;
        break;
    }
  }

  PUSH( beg, end0, depth );  // the '0's
  PUSH( end0, end, depth );  // the '1's
}

// Specialized sort for 3 symbols + 1 delimiter.
// E.g. subset seeds for bisulfite-converted DNA.
static void radixSort3( const uchar* text, const uchar* subsetMap,
                        indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* beg2 = end;  // beginning of '2's

  while( end1 < beg2 ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      case 2:
        *end1 = *--beg2;
        *beg2 = x;
        break;
      default:  // the delimiter subset
        *end1 = *--beg2;
        *beg2 = *--end;
        *end = x;
        break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's
  PUSH( end0, beg2, depth );  // the '1's
  PUSH( beg2, end, depth );   // the '2's
}

// Specialized sort for 4 symbols + 1 delimiter.  E.g. DNA.
static void radixSort4( const uchar* text, const uchar* subsetMap,
			indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* end2 = beg;  // end of '2's
  indexT* beg3 = end;  // beginning of '3's

  while( end2 < beg3 ){
    const indexT x = *end2;
    switch( subsetMap[ text[x] ] ){
    case 0:
      *end2++ = *end1;
      *end1++ = *end0;
      *end0++ = x;
      break;
    case 1:
      *end2++ = *end1;
      *end1++ = x;
      break;
    case 2:
      end2++;
      break;
    case 3:
      *end2 = *--beg3;
      *beg3 = x;
      break;
    default:  // the delimiter subset
      *end2 = *--beg3;
      *beg3 = *--end;
      *end = x;
      break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's
  PUSH( end0, end1, depth );  // the '1's
  PUSH( end1, beg3, depth );  // the '2's
  PUSH( beg3, end, depth );  // the '3's
}

static void radixSortN( const uchar* text, const uchar* subsetMap,
			indexT* beg, indexT* end, indexT depth,
			unsigned subsetCount ){
  static indexT bucketSize[256];  // initialized to zero at startup
  /*  */ indexT* bucketEnd[256];  // "static" makes little difference to speed

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for( indexT* i = beg; i < end; /* noop */ ){
    uchar oracle[256];
    uchar* oracleEnd =
      oracle + std::min( sizeof(oracle), std::size_t(end - i) );
    for( uchar* j = oracle; j < oracleEnd; ++j )
      *j = subsetMap[ text[ *i++ ] ];
    for( uchar* j = oracle; j < oracleEnd; ++j )
      ++bucketSize[ *j ];
  }

  // get bucket ends, and put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  indexT* pos = beg;
  for( unsigned i = 0; i < subsetCount; ++i ){
    indexT* nextPos = pos + bucketSize[i];
    PUSH( pos, nextPos, depth );
    pos = nextPos;
    bucketEnd[i] = pos;
  }
  // don't sort within the delimiter bucket:
  bucketEnd[ CyclicSubsetSeed::DELIMITER ] = end;

  // permute items into the correct buckets:
  for( indexT* i = beg; i < end; /* noop */ ) {
    unsigned subset;  // unsigned is faster than uchar!
    indexT holdOut = *i;
    while( --bucketEnd[ subset = subsetMap[ text[holdOut] ] ] > i ){
      std::swap( *bucketEnd[subset], holdOut );
    }
    *i = holdOut;
    i += bucketSize[subset];
    bucketSize[subset] = 0;  // reset it so we can reuse it
  }
}

void SubsetSuffixArray::sortIndex( const uchar* text,
				   indexT maxUnsortedInterval ){
  const indexT minLength = 1;

  PUSH( &index.v.front(), &index.v.back() + 1, 0 );

  while( sp > stack ){
    indexT* beg;
    indexT* end;
    indexT depth;
    POP( beg, end, depth );

    if( end - beg <= maxUnsortedInterval && depth >= minLength ) continue;

    if( end - beg < 10 ){  // ???
      insertionSort( text, seed, beg, end, depth );
      continue;
    }

    const uchar* textBase = text + depth;
    const uchar* subsetMap = seed.subsetMap(depth);
    unsigned subsetCount = seed.subsetCount(depth);

    ++depth;

    switch( subsetCount ){
      case 1:  radixSort1( textBase, subsetMap, beg, end, depth );  break;
      case 2:  radixSort2( textBase, subsetMap, beg, end, depth );  break;
      case 3:  radixSort3( textBase, subsetMap, beg, end, depth );  break;
      case 4:  radixSort4( textBase, subsetMap, beg, end, depth );  break;
      default: radixSortN( textBase, subsetMap, beg, end, depth, subsetCount );
    }
  }
}
