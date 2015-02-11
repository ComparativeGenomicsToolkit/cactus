// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "MultiSequence.hh"
#include "io.hh"
#include <sstream>
#include <algorithm>  // upper_bound
#include <cassert>
#include <cctype>  // isspace
#include <iterator>  // istreambuf_iterator

using namespace cbrc;

void MultiSequence::initForAppending( indexT padSizeIn ){
  padSize = padSizeIn;
  seq.v.assign( padSize, ' ' );
  ends.v.assign( 1, padSize );
  names.v.clear();
  nameEnds.v.assign( 1, 0 );
}

void MultiSequence::reinitForAppending(){
  seq.v.erase( seq.v.begin(), seq.v.begin() + ends.v.back() - padSize );
  names.v.erase( names.v.begin(),
                 names.v.begin() + nameEnds.v[ finishedSequences() ] );
  ends.v.resize(1);
  nameEnds.v.resize(1);
  if( !names.v.empty() ) nameEnds.v.push_back( names.v.size() );
}

void MultiSequence::fromFiles( const std::string& baseName, indexT seqCount,
                               std::size_t qualitiesPerLetter ){
  ends.m.open( baseName + ".ssp", seqCount + 1 );
  seq.m.open( baseName + ".tis", ends.m.back() );
  nameEnds.m.open( baseName + ".sds", seqCount + 1 );
  names.m.open( baseName + ".des", nameEnds.m.back() );
  padSize = ends.m[0];

  qualityScores.m.open( baseName + ".qua",
                        ends.m.back() * qualitiesPerLetter );
}

void MultiSequence::toFiles( const std::string& baseName ) const{
  memoryToBinaryFile( ends.begin(), ends.end(), baseName + ".ssp" );

  memoryToBinaryFile( seq.begin(), seq.begin() + ends.back(),
		      baseName + ".tis" );

  memoryToBinaryFile( nameEnds.begin(), nameEnds.begin() + ends.size(),
		      baseName + ".sds" );

  memoryToBinaryFile( names.begin(),
		      names.begin() + nameEnds[ finishedSequences() ],
		      baseName + ".des" );

  memoryToBinaryFile( qualityScores.begin(),
                      qualityScores.begin() + ends.back() * qualsPerLetter(),
                      baseName + ".qua" );
}

void MultiSequence::addName( std::string& name ){
  names.v.insert( names.v.end(), name.begin(), name.end() );
  nameEnds.v.push_back( names.v.size() );
  if( nameEnds.v.back() < names.v.size() )
    throw std::runtime_error("the sequence names are too long");
}

std::istream& MultiSequence::readFastaName( std::istream& stream ){
  std::string line, word;
  getline( stream, line );
  std::istringstream iss(line);
  iss >> word;
  if( !stream ) return stream;
  addName(word);
  return stream;
}

std::istream&
MultiSequence::appendFromFasta( std::istream& stream, indexT maxSeqLen ){
  if( isFinished() ){
    char c = '>';
    stream >> c;
    if( c != '>' )
      throw std::runtime_error("bad FASTA sequence data: missing '>'");
    readFastaName(stream);
    if( !stream ) return stream;
  }

  std::istreambuf_iterator<char> inpos(stream);
  std::istreambuf_iterator<char> endpos;
  while( inpos != endpos ){
    uchar c = *inpos;
    if( c == '>' ) break;  // we have hit the next FASTA sequence
    if( !std::isspace(c) ){
      if( seq.v.size() >= maxSeqLen ) break;
      seq.v.push_back(c);
    }
    ++inpos;
  }

  if( isFinishable(maxSeqLen) ) finish();
  return stream;
}

void MultiSequence::finish(){
  assert( !isFinished() );
  seq.v.insert( seq.v.end(), padSize, ' ' );
  ends.v.push_back( seq.v.size() );
  assert( ends.v.back() == seq.v.size() );
}

void MultiSequence::unfinish(){
  assert( isFinished() );
  ends.v.pop_back();
  seq.v.erase( seq.v.end() - padSize, seq.v.end() );
}

bool MultiSequence::isFinishable( indexT maxSeqLen ) const{
  return seq.v.size() + padSize <= maxSeqLen;
}

MultiSequence::indexT MultiSequence::whichSequence( indexT coordinate ) const{
  const indexT* u = std::upper_bound( ends.begin(), ends.end(), coordinate );
  assert( u != ends.begin() && u != ends.end() );
  return u - ends.begin() - 1;
}

std::string MultiSequence::seqName( indexT seqNum ) const{
  return std::string( names.begin() + nameEnds[ seqNum ],
		      names.begin() + nameEnds[ seqNum + 1 ] );
}
