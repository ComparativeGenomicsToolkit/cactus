// Copyright 2008 Martin C. Frith

#include "DiagonalTable.hh"

namespace cbrc{

bool DiagonalTable::isCovered( indexT sequentialPos, indexT randomPos ){

  indexT diagonal = sequentialPos - randomPos;  // wrap-around is OK
  indexT bin = diagonal % BINS;
  std::vector<pairT>& v = hits[bin];

  for( std::vector<pairT>::iterator i = v.begin(); i < v.end(); /* noop */ ){
    if( i->first >= sequentialPos ){
      if( i->second == diagonal ) return true;
      ++i;
    }else{
      i = v.erase(i);  // hopefully we rarely get here
    }
  }

  return false;
}

void DiagonalTable::addEndpoint( indexT sequentialPos, indexT randomPos ){

  indexT diagonal = sequentialPos - randomPos;  // wrap-around is OK
  indexT bin = diagonal % BINS;
  std::vector<pairT>& v = hits[bin];

  v.push_back( pairT( sequentialPos, diagonal ) );
}

}  // end namespace cbrc
