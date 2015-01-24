// Copyright 2010 Martin C. Frith

// This struct holds the command line arguments for lastex.

#ifndef LASTEX_ARGUMENTS_HH
#define LASTEX_ARGUMENTS_HH
#include <string>
#include <iosfwd>

namespace cbrc{

struct LastexArguments{
  // set the parameters to their default values:
  LastexArguments();

  // set parameters from a list of arguments:
  void fromArgs( int argc, char** argv, bool optionsOnly = false );

  // set parameters from a command line (by splitting it into arguments):
  void fromLine( const std::string& line, bool optionsOnly = true );

  // set parameters from lines beginning with "#last":
  void fromStream( std::istream& is, bool optionsOnly = true );

  // set parameters from lines beginning with "#last":
  void fromString( const std::string& s, bool optionsOnly = true );

  // set default option values that depend on input files:
  void setDefaultsFromAlphabet( bool isDna, bool isProtein );

  // options:
  int strand;
  int matchScore;
  int mismatchCost;
  int gapExistCost;
  int gapExtendCost;
  std::string matrixFile;
  bool isGapless;
  int score;
  double maxEvalue;
  int searchSpace;

  // positional arguments:
  std::string targetStatsFile;
  std::string queryStatsFile;
  std::string alignmentsFile;
};

}  // end namespace

#endif
