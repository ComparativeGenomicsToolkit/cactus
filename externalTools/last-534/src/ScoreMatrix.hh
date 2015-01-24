// Copyright 2008, 2009, 2011, 2014 Martin C. Frith

// This struct holds a score matrix for aligning pairs of residues,
// e.g. blosum62.  The delimiter symbol (space) aligned to anything
// gets a score of -INF.

// Maybe split this struct into two: ScoreMatrixEasy and ScoreMatrixFast?

#ifndef SCOREMATRIX_HH
#define SCOREMATRIX_HH
#include <string>
#include <vector>
#include <iosfwd>
#include <climits>  // INT_MAX

namespace cbrc{

struct ScoreMatrix{
  typedef unsigned char uchar;

  enum { INF = INT_MAX / 2 };  // big, but try to avoid overflow
  enum { MAT = 64 };           // matrix size = MAT x MAT
  enum { OUTPAD = 2 };         // cell-padding for output

  static std::string stringFromName( const std::string& name );

  void matchMismatch( int match, int mismatch, const std::string& letters );
  void fromString( const std::string& s );
  void init( const uchar encode[] );  // unspecified letters get minScore
  void writeCommented( std::ostream& stream ) const;  // write preceded by "#"

  std::string rows;                       // row headings (letters)
  std::string cols;                       // column headings (letters)
  std::vector< std::vector<int> > cells;  // scores
  int caseSensitive[MAT][MAT];
  int caseInsensitive[MAT][MAT];
  int minScore;
  int maxScore;
};

std::istream& operator>>( std::istream& stream, ScoreMatrix& mat );
std::ostream& operator<<( std::ostream& stream, const ScoreMatrix& mat );

}  // end namespace cbrc
#endif  // SCOREMATRIX_HH
