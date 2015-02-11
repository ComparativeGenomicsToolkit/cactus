// Copyright 2010, 2014 Martin C. Frith

#include "LastexArguments.hh"
#include "ScoreMatrix.hh"
#include "Alphabet.hh"
#include "io.hh"
#include "gumbel_params/mcf_local_alignment_evaluer.hpp"

#include <cctype>  // isdigit
#include <cmath>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>
#include <iomanip>  // setw
#include <iostream>
#include <numeric>  // accumulate
#include <sstream>
#include <stdexcept>

#define ERR(x) throw std::runtime_error(x)
#define WARN(x) std::cerr << "lastex: warning: " << x << '\n'

using namespace cbrc;

// Scale the elements of v so that they sum to 1
void normalize( std::vector<double>& v ){
  double t = std::accumulate( v.begin(), v.end(), 0.0 );
  // assume that t != 0

  for( std::vector<double>::iterator i = v.begin(); i < v.end(); ++i )
    *i /= t;
}

struct SequenceStatistics{
  double letterCount;
  double sequenceCount;
  std::string alphabet;
  std::vector<double> letterProbs;

  SequenceStatistics()
      : letterCount(0), sequenceCount(0), alphabet(""), letterProbs(0) {}
};

bool isBad( const SequenceStatistics& s ){
  if( s.letterCount <= 0 ) return true;
  if( s.sequenceCount <= 0 ) return true;
  if( s.alphabet.size() == 0 ) return true;
  if( s.letterProbs.size() != s.alphabet.size() ) return true;

  for( unsigned i = 0; i < s.letterProbs.size(); ++i ){
    if( s.letterProbs[i] < 0 ) return true;
  }

  for( unsigned i = 0; i < s.letterProbs.size(); ++i ){
    if( s.letterProbs[i] > 0 ) return false;
  }

  return true;
}

SequenceStatistics readStats( const std::string& fileName ){
  SequenceStatistics s;

  std::ifstream inFileStream;
  std::istream& f = openIn( fileName, inFileStream );

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream lineStream(line);
    getline( lineStream, word, '=' );
    if( word == "numofsequences" ) lineStream >> s.sequenceCount;
    if( word == "numofletters" ) lineStream >> s.letterCount;
    if( word == "alphabet" ) lineStream >> s.alphabet;
    if( word == "letterfreqs" ){
      double d;
      while( lineStream >> d ){
        s.letterProbs.push_back(d);
      }
    }
  }

  if( f.eof() && !f.bad() ) f.clear();
  if( isBad(s) ) f.setstate( std::ios::failbit );
  if( !f ) ERR( "can't read this counts file: " + fileName );

  normalize( s.letterProbs );
  return s;
}

namespace{
  LastexArguments args;
  Alphabet alph;
  ScoreMatrix scoreMatrix;
  SequenceStatistics stats1;
  SequenceStatistics stats2;
  Mcf::LocalAlignmentEvaluer evaluer;
}

// Set up a scoring matrix, based on the user options
void makeScoreMatrix( const std::string& matrixFile ){
  if( !matrixFile.empty() ){
    scoreMatrix.fromString( matrixFile );
  }
  else if( args.matchScore < 0 && args.mismatchCost < 0 && alph.isProtein() ){
    scoreMatrix.fromString( ScoreMatrix::stringFromName( "BL62" ) );
  }
  else{
    scoreMatrix.matchMismatch( args.matchScore, args.mismatchCost,
                               alph.letters );
  }

  scoreMatrix.init( alph.encode );
}

// Symmetrize the probabilities of letters and their complements
void addComplementProbs( std::vector<double>& letterProbs ){
  for( unsigned i = 0; i < letterProbs.size(); ++i ){
    unsigned j = alph.complement[i];
    if( j >= letterProbs.size() ) continue;  // avoid crashing on non-DNA
    if( j < i ) continue;
    double p = (letterProbs[i] + letterProbs[j]) / 2;
    letterProbs[i] = p;
    letterProbs[j] = p;
  }
}

void makeStrandStats(){
  if( args.strand < 2 ) return;
  addComplementProbs( stats1.letterProbs );
  addComplementProbs( stats2.letterProbs );
  stats2.letterCount *= 2;
  stats2.sequenceCount *= 2;
}

void makeEvaluer(){
  // we need to pass the score matrix as a pointer-to-pointers:
  std::vector<const int*> matrixPointers;
  for( unsigned i = 0; i < alph.size; ++i ){
    matrixPointers.push_back( &scoreMatrix.caseSensitive[i][0] );
  }

  if( args.isGapless)
    evaluer.initGapless( stats1.letterProbs, stats2.letterProbs,
                         &matrixPointers[0] );
  else
    evaluer.initGapped( stats1.letterProbs, stats2.letterProbs,
                        &matrixPointers[0],
                        args.gapExistCost, args.gapExtendCost );

  if( evaluer.isBad() )
    ERR( "can't evaluate scores for these score parameters and letter frequencies" );
}

// ******* Routines for printing a table of scores and E-values *******

void writeParameterWithError( const std::string& name, double p, double err ){
  std::cout << name << "\t" << p;
  if( !args.isGapless ) std::cout << "\t+/- " << err;
  std::cout << "\n";
}

void writeGumbelParameters(){
  const Mcf::LocalAlignmentEvaluer::GumbelParameters& g
      = evaluer.gumbelParameters();

  std::cout << "Parameters\n";
  writeParameterWithError( "lambda", g.lambda, g.lambda_error );
  writeParameterWithError( "k", g.K, g.K_error );
  writeParameterWithError( "c", g.C, g.C_error );
  writeParameterWithError( "aI", g.a_I, g.a_I_error );
  writeParameterWithError( "aJ", g.a_J, g.a_J_error );
  writeParameterWithError( "alphaI", g.alpha_I, g.alpha_I_error );
  writeParameterWithError( "alphaJ", g.alpha_J, g.alpha_J_error );
  writeParameterWithError( "sigma", g.sigma, g.sigma_error );
}

void writeLetterPercents(){
  std::cout.precision(2);
  std::cout << std::left;
  std::cout << "Letter percentages\n";

  for( unsigned i = 0; i < alph.size; ++i ){
    if( i > 0 ) std::cout << " ";
    std::cout << std::setw(3) << alph.letters[i];
  }
  std::cout << "\n";

  for( unsigned i = 0; i < alph.size; ++i ){
    if( i > 0 ) std::cout << " ";
    std::cout << std::setw(3) << 100 * stats1.letterProbs[i];
  }
  std::cout << "\n";

  for( unsigned i = 0; i < alph.size; ++i ){
    if( i > 0 ) std::cout << " ";
    std::cout << std::setw(3) << 100 * stats2.letterProbs[i];
  }
  std::cout << "\n";
}

void writeScoreAndEvalue( int score ){
  double evalue = evaluer.evalue( score,
                                  stats1.letterCount, stats2.letterCount,
                                  stats1.sequenceCount, stats2.sequenceCount );
  std::cout << score << "\t" << evalue << "\n";
}

void writeEvalues( const std::string& matrixString ){
  args.setDefaultsFromAlphabet( alph.letters == alph.dna, alph.isProtein() );
  makeScoreMatrix( matrixString );
  makeStrandStats();
  makeEvaluer();

  std::cout << "LAST version " <<
#include "version.hh"
            << "\n";
  std::cout << "\n";
  std::cout << "Score\tExpected number of alignments\n";

  if( args.score >= 0 ){
    writeScoreAndEvalue( args.score );
  }

  if( args.maxEvalue >= 0 ){
    int score = evaluer.minScore( args.maxEvalue,
                                  stats1.letterCount, stats2.letterCount,
                                  stats1.sequenceCount, stats2.sequenceCount );
    writeScoreAndEvalue( score );
  }

  if( args.score < 0 && args.maxEvalue < 0 ){
    int oldScore = 0;
    for( int i = -10; i <= 10; ++i ){
      double maxEvalue = std::pow( 10.0, i );  // the ".0" is needed!
      int s = evaluer.minScore( maxEvalue,
                                stats1.letterCount, stats2.letterCount,
                                stats1.sequenceCount, stats2.sequenceCount );
      if( s == oldScore ) continue;
      writeScoreAndEvalue( s );
      oldScore = s;
    }
  }

  std::cout << "\n";
  writeGumbelParameters();
  std::cout << "\n";
  writeLetterPercents();
}

// ******* Routines for reading alignments and assigning E-values *******

// Is the scoring scheme fully defined?
bool isDefinedScoringScheme(){
  if( args.matrixFile.empty() )
    if( args.matchScore < 0 || args.mismatchCost < 0)
      return false;
  if( !args.isGapless )
    if( args.gapExistCost < 0 || args.gapExtendCost < 0 )
      return false;
  return true;
}

char firstCharInLine( const std::string& line ){
  if( line.empty() ) return '\n';
  else return line[0];
}

// Return a copy of the line without the initial comment symbol, if any
std::string uncommented( const std::string& line ){
  if( firstCharInLine(line) == '#' ) return line.substr(1);
  else return line;
}

// Can this be the first line of a scoring matrix?
bool isMatrixHead( const std::string& line ){
  std::istringstream lineStream(line);
  std::string word;
  while( lineStream >> word ){
    if( word.size() > 1 ) return false;
  }
  if( word.empty() ) return false;
  return true;
}

// Can this be an interior line of a scoring matrix?
bool isMatrixBody( const std::string& line ){
  std::istringstream lineStream(line);
  std::string word;
  int score;
  lineStream >> word >> score;
  if( !lineStream || word.size() > 1 ) return false;
  return true;
}

// Read a value from a line, in the form "key=value"
// Return -1 if not found
int valueFromKey( const std::string& line, const std::string& key ){
  int value = -1;
  std::istringstream lineStream(line);
  std::string word, k;
  while( lineStream >> word ){
    std::istringstream wordStream(word);
    getline( wordStream, k, '=' );
    if( k == key ) wordStream >> value;
  }
  return value;
}

void parseParameter( const std::string& line, const std::string& key,
                     int& value ){
  int v = valueFromKey( line, key );
  if( v >= 0 ) value = v;
}

// Try to get alignment parameters from the header of an alignments file
void parseHeader( const std::vector<std::string>& lines,
                  const std::string& matrixString ){
  int a = -1, b = -1, c = -1, F = -1, j = -1, s = -1;
  std::string headerMatrix;

  for( unsigned i = 0; i < lines.size(); ++i ){
    const std::string line = uncommented( lines[i] );
    parseParameter( line, "a", a );
    parseParameter( line, "b", b );
    parseParameter( line, "c", c );
    parseParameter( line, "F", F );
    parseParameter( line, "j", j );
    parseParameter( line, "s", s );
    if( headerMatrix.empty() ? isMatrixHead(line) : isMatrixBody(line) ){
      headerMatrix += line + "\n";
    }
  }

  if( headerMatrix.empty() || a < 0 || b < 0 || j < 0 || s < 0 ){
    if( !isDefinedScoringScheme() )
      ERR( "please define all match, mismatch, and gap scores" );
    args.setDefaultsFromAlphabet( alph.letters == alph.dna, alph.isProtein() );
    makeScoreMatrix( matrixString );
    makeStrandStats();
  }
  else{  // we got the parameters from the header successfully
    if( args.matchScore < 0 && args.mismatchCost < 0 &&
        args.matrixFile.empty() && args.strand < 0 &&
        args.gapExistCost < 0 && args.gapExtendCost < 0 && !args.isGapless ){
      if( F > 0 ){
        WARN( "the E-values do not take frame-shifts into account" );
      }

      if( c > 0 && c < a + 2 * b ){
        WARN( "the E-values do not take generalized gaps into account" );
      }

      args.gapExistCost = a;
      args.gapExtendCost = b;
      args.isGapless = (j < 2);
      args.strand = s;
      makeScoreMatrix( headerMatrix );

      // F>0 means translated alignment.  In this case, the DNA counts
      // file should have counts for translated DNA (3 frames or 6
      // frames).
      if( F <= 0 ) makeStrandStats();
    }
    else{
      ERR( "can't redefine parameters that are in the alignments file" );
    }
  }

  makeEvaluer();
}

// Calculate E-values using either total lengths from count files, or
// individual lengths of the aligned sequences
double evalueForSequences( int score, double seqLen1, double seqLen2 ){
  // ugly, changes global stats1 and stats2 objects

  if( args.searchSpace / 2 ){
    stats1.letterCount = seqLen1;
    stats1.sequenceCount = 1;
  }

  if( args.searchSpace % 2 ){
    stats2.letterCount = seqLen2;
    stats2.sequenceCount = 1;

    if( args.strand > 1 ){
      stats2.letterCount *= 2;
      stats2.sequenceCount *= 2;
    }
  }

  return evaluer.evalue( score, stats1.letterCount, stats2.letterCount,
                         stats1.sequenceCount, stats2.sequenceCount );
}

// Parse the sequence length from a MAF alignment line
double lengthFromMafLine( const std::string& line ){
  double x;
  std::istringstream iss(line);
  std::string s;
  iss >> s >> s >> s >> s >> s >> x;
  if( !iss ) ERR( "bad MAF alignment line:\n" + line );
  return x;
}

// Parse the sequence lengths from MAF alignment lines
std::vector<double> lengthsFromMaf( const std::vector<std::string>& lines ){
  std::vector<double> lengths;
  for( unsigned i = 0; i < lines.size(); ++i ){
    if( firstCharInLine( lines[i] ) == 's' ){
      lengths.push_back( lengthFromMafLine( lines[i] ) );
    }
  }
  return lengths;
}

void parseMaf( const std::vector<std::string>& lines ){
  int score = valueFromKey( lines[0], "score" );
  if( score < 0 ) return;

  std::vector<double> seqLengths = lengthsFromMaf(lines);
  if( seqLengths.size() != 2 ) ERR( "non-pairwise MAF alignment" );

  double evalue = evalueForSequences( score, seqLengths[0], seqLengths[1] );
  if( args.maxEvalue >= 0 && evalue > args.maxEvalue ) return;

  std::cout << lines[0] << " expect=" << evalue << "\n";

  for( unsigned i = 1; i < lines.size(); ++i ){
    std::cout << lines[i] << "\n";
  }
}

void parseTabular( const std::vector<std::string>& lines ){
  std::istringstream iss( lines[0] );
  int score;
  double seqLen1, seqLen2;
  std::string s;
  iss >> score >> s >> s >> s >> s >> seqLen1 >> s >> s >> s >> s >> seqLen2;
  if( !iss ) ERR( "bad tabular alignment line:\n" + lines[0] );

  double evalue = evalueForSequences( score, seqLen1, seqLen2 );
  if( args.maxEvalue >= 0 && evalue > args.maxEvalue ) return;

  std::cout << lines[0] << "\t" << evalue << "\n";

  for( unsigned i = 1; i < lines.size(); ++i ){
    std::cout << lines[i] << "\n";
  }
}

void parseComment( const std::vector<std::string>& lines ){
  for( unsigned i = 0; i < lines.size(); ++i ){
    std::cout << lines[i] << "\n";
  }
}

// Process one chunk of an alignments file
void parseChunk( const std::vector<std::string>& lines,
                 const std::string& matrixString ){
  if( evaluer.isBad() ) parseHeader( lines, matrixString );
  if( lines.empty() ) return;
  char c = firstCharInLine( lines[0] );
  if( c == 'a' ) parseMaf( lines );
  else if( std::isdigit(c) ) parseTabular( lines );
  else parseComment( lines );
}

void evaluateAlignments( const std::string& matrixString ){
  std::ifstream inFileStream;
  std::istream& f = openIn( args.alignmentsFile, inFileStream );

  std::vector<std::string> chunk;
  std::string line;
  bool isComment = true;

  while( getline( f, line ) ){
    char c = firstCharInLine(line);

    if( c == 'a' || std::isdigit(c) || (c == '#' && !isComment) ){
      parseChunk( chunk, matrixString );
      chunk.clear();
      isComment = (c == '#');
    }

    chunk.push_back(line);
  }

  parseChunk( chunk, matrixString );
}

// ******* Main program *******

void lastex( int argc, char** argv ){
  args.fromArgs( argc, argv );

  std::string matrixString;
  if( !args.matrixFile.empty() ){
    matrixString = ScoreMatrix::stringFromName( args.matrixFile );
    args.fromString( matrixString );  // read options from the matrix file
    args.fromArgs( argc, argv );  // command line overrides matrix file
  }

  stats1 = readStats( args.targetStatsFile );
  stats2 = readStats( args.queryStatsFile );

  if( stats1.alphabet != stats2.alphabet )
    ERR( "the counts files have different alphabets" );

  alph.fromString( stats1.alphabet );

  std::cout.precision(3);

  if( args.alignmentsFile.empty() ) writeEvalues( matrixString );
  else evaluateAlignments( matrixString );
}

int main( int argc, char** argv )
try{
  lastex( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << "lastex: out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << "lastex: " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
