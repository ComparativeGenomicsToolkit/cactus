// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

#include "CyclicSubsetSeed.hh"
#include "io.hh"
#include "stringify.hh"
#include <fstream>
#include <sstream>
#include <algorithm>  // sort
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

using namespace cbrc;

std::vector<std::string> CyclicSubsetSeed::fromName( const std::string& name ){
  const char* seedAlph[256] = {0};
  seedAlph['1'] = "A C G T";
  seedAlph['0'] = "ACGT";
  seedAlph['T'] = "AG CT";

  if( name == "BISF" ){
    seedAlph['1'] = "CT A G";
    return fromMask( seedAlph, "1111110101100" );
  }

  if( name == "BISR" ){
    seedAlph['1'] = "AG C T";
    return fromMask( seedAlph, "1111110101100" );
  }

  if( name == "MAM4" ){
    // From MC Frith & L Noe (2014) Nucleic Acids Research,
    // Supplementary Table 11, row 12.
    return fromMask( seedAlph,
		     "11100TT01T00T10TTTT,"
		     "TTTT110TT0T001T0T1T1,"
		     "11TT010T01TT0001T,"
		     "11TT10T1T101TT" );
  }

  if( name == "MAM8" ){
    // From MC Frith & L Noe (2014) Nucleic Acids Research,
    // Supplementary Table 12, second-last row.
    return fromMask( seedAlph,
		     "1101T1T0T1T00TT1TT,"
		     "1TTTTT010TT0TT01011TTT,"
		     "1TTTT10010T011T0TTTT1,"
		     "111T011T0T01T100,"
		     "1T10T100TT01000TT01TT11,"
		     "111T101TT000T0T10T00T1T,"
		     "111100T011TTT00T0TT01T,"
		     "1T1T10T1101101" );
  }

  std::string s = slurp( name );
  std::vector<std::string> v( 1, s );
  return v;
}

void CyclicSubsetSeed::fromString( const std::string& s,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  std::istringstream iss(s);
  fromStream( iss, isMaskLowercase, letterCode );
}

static bool isBlankOrComment( std::string line ){
  std::istringstream iss(line);
  char c;
  if( !(iss >> c) ) return true;
  if( c == '#' ) return true;
  return false;
}

void CyclicSubsetSeed::fromStream( std::istream& stream,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  clear();
  std::string line;
  while( std::getline( stream, line ) ){
    if( isBlankOrComment(line) ) continue;
    std::istringstream iss(line);
    appendPosition( iss, isMaskLowercase, letterCode );
  }
  if( span() && stream.eof() ) stream.clear( std::ios::eofbit );
}

static std::string exactSeed( const std::string& letters ){
  std::string result;
  for( unsigned i = 0; i < letters.size(); ++i ){
    if( i > 0 ) result += ' ';
    result += letters[i];
  }
  return result;
}

std::vector<std::string> CyclicSubsetSeed::fromMask( const std::string& alph,
						     const std::string& mask ){
  std::string es = exactSeed(alph);
  const char* seedAlph[256] = {0};
  seedAlph['1'] = es.c_str();
  seedAlph['#'] = es.c_str();
  seedAlph['0'] = alph.c_str();
  seedAlph['_'] = alph.c_str();
  seedAlph['-'] = alph.c_str();
  seedAlph['T'] = "AG CT";
  seedAlph['t'] = "AG CT";
  seedAlph['@'] = "AG CT";
  return fromMask( seedAlph, mask );
}

std::vector<std::string> CyclicSubsetSeed::fromMask( const char* seedAlph[],
						     const std::string& mask ){
  std::vector<std::string> v;
  int n = 0;

  for( unsigned i = 0; i < mask.size(); ++i ){
    uchar c = mask[i];
    if( c == ',' ){
      n = 0;
    }else{
      const char* x = seedAlph[c];
      if( !x ) ERR( "bad seed pattern: " + mask );
      if( !n++ ) v.push_back("");
      v.back() += x;
      v.back() += "\n";
    }
  }

  return v;
}

void CyclicSubsetSeed::addLetter( std::vector<uchar>& numbersToSubsets,
				  uchar letter, uchar subsetNum,
				  const uchar letterCode[] ){
  uchar number = letterCode[letter];
  if( number >= CyclicSubsetSeed::MAX_LETTERS )
    ERR( "bad symbol in subset-seed: " + stringify(letter) );
  if( numbersToSubsets[number] < CyclicSubsetSeed::DELIMITER )
    ERR( "repeated symbol in subset-seed: "  + stringify(letter) );
  numbersToSubsets[number] = subsetNum;
}

void CyclicSubsetSeed::appendPosition( std::istream& inputLine,
				       bool isMaskLowercase,
				       const uchar letterCode[] ){
  std::string inputWord;
  std::vector<std::string> subsetList;
  std::vector<uchar> numbersToSubsets( MAX_LETTERS, DELIMITER );

  for( unsigned subsetNum = 0; inputLine >> inputWord; ++subsetNum ){
    assert( subsetNum < DELIMITER );
    std::string subset;

    for( unsigned i = 0; i < inputWord.size(); ++i ){
      uchar upper = std::toupper( inputWord[i] );
      uchar lower = std::tolower( inputWord[i] );
      addLetter( numbersToSubsets, upper, subsetNum, letterCode );
      subset += upper;
      if( !isMaskLowercase && lower != upper ){
	addLetter( numbersToSubsets, lower, subsetNum, letterCode );
      }
    }

    std::sort( subset.begin(), subset.end() );  // canonicalize
    subsetList.push_back( subset );
  }

  subsetLists.push_back( subsetList );
  subsetMaps.insert( subsetMaps.end(),
		     numbersToSubsets.begin(), numbersToSubsets.end() );
}

void CyclicSubsetSeed::writePosition( std::ostream& out,
				      unsigned position ) const{
  assert( position < subsetLists.size() );
  for( unsigned i = 0; i < subsetLists[position].size(); ++i ){
    if( i > 0 ) out << ' ';
    out << subsetLists[position][i];
  }
}
