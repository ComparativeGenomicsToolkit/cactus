// Copyright 2010, 2014 Martin C. Frith

#include "LastexArguments.hh"
#include "stringify.hh"
#include <unistd.h>  // getopt
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cstring>  // strtok
#include <cstdlib>  // EXIT_SUCCESS
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

static void badopt( char opt, const char* arg ){
  ERR( std::string("bad option value: -") + opt + ' ' + arg );
}

namespace cbrc{

LastexArguments::LastexArguments() :
  strand(-1),  // depends on the alphabet
  matchScore(-1),  // depends on the alphabet
  mismatchCost(-1),  // depends on the alphabet
  gapExistCost(-1),  // depends on the alphabet
  gapExtendCost(-1),  // depends on the alphabet
  matrixFile(""),
  isGapless(false),
  score(-1),  // means: no score chosen
  maxEvalue(-1),  // means: no maximum E-value
  searchSpace(0){}

void LastexArguments::fromArgs( int argc, char** argv, bool optionsOnly ){
  std::string usage = "\
Usage:\n\
lastex [options] reference-counts-file query-counts-file [alignments-file]\n\
Calculate expected numbers of alignments for random sequences.\n\
\n\
Options (default settings):\n\
-h: show all options and their default settings\n\
-s: strands (2 for DNA, 1 for protein)\n\
-r: match score   (DNA: 1, protein: blosum62)\n\
-q: mismatch cost (DNA: 1, protein: blosum62)\n\
-p: match/mismatch score matrix\n\
-a: gap existence cost (DNA: 7, protein: 11)\n\
-b: gap extension cost (DNA: 1, protein:  2)\n\
-g: do calculations for gapless alignments\n\
-y: find the expected number of alignments with score >= this\n\
-E: maximum expected number of alignments\n\
-z: calculate the expected number of alignments per:\n\
    0 = reference counts file / query counts file\n\
    1 = reference counts file / each query sequence\n\
    2 = each reference sequence / query counts file\n\
    3 = each reference sequence / each query sequence\n\
    (0)";

  std::string help = usage + "\n\
\n\
Report bugs to: last-align (ATmark) googlegroups (dot) com\n\
LAST home page: http://last.cbrc.jp/\n\
";

  optind = 1;  // allows us to scan arguments more than once(???)
  int c;
  while( (c = getopt(argc, argv, "hs:r:q:p:a:b:gy:E:z:e:")) != -1 ){
    switch(c){
    case 'h':
      std::cout << help;
      throw EXIT_SUCCESS;
    case 's':
      unstringify( strand, optarg );
      if( strand < 1 || strand > 2 ) badopt( c, optarg );
      break;
    case 'r':
      unstringify( matchScore, optarg );
      if( matchScore <= 0 ) badopt( c, optarg );
      break;
    case 'q':
      unstringify( mismatchCost, optarg );
      if( mismatchCost <= 0 ) badopt( c, optarg );
      break;
    case 'p':
      matrixFile = optarg;
      break;
    case 'a':
      unstringify( gapExistCost, optarg );
      if( gapExistCost < 0 ) badopt( c, optarg );
      break;
    case 'b':
      unstringify( gapExtendCost, optarg );
      if( gapExtendCost <= 0 ) badopt( c, optarg );
      break;
    case 'g':
      isGapless = true;
      break;
    case 'y':
      unstringify( score, optarg );
      if( score <= 0 ) badopt( c, optarg );
      break;
    case 'E':
      unstringify( maxEvalue, optarg );
      if( maxEvalue < 0 ) badopt( c, optarg );
      break;
    case 'z':
      unstringify( searchSpace, optarg );
      if( searchSpace < 0 || searchSpace > 3 ) badopt( c, optarg );
      break;
    case 'e':
      // Allow but ignore this option, so we don't reject lastal
      // matrix files.  Maybe put other lastal options here too,
      // e.g. -d, -x
      break;
    case '?':
      ERR( "bad option" );
    }
  }

  if( optionsOnly ) return;
  if( optind + 2 > argc || optind + 3 < argc )
    ERR( "please give me two or three file names\n\n" + usage );
  targetStatsFile = argv[optind++];
  queryStatsFile = argv[optind++];
  if( optind < argc ) alignmentsFile = argv[optind];

  if( alignmentsFile.empty() && searchSpace > 0 )
    ERR( "can't use option -z > 0 without an alignments file" );

  if( !alignmentsFile.empty() && score >= 0 )
    ERR( "can't use option -y with an alignments file" );
}

void LastexArguments::fromLine( const std::string& line, bool optionsOnly ){
  std::vector<char> args( line.begin(), line.end() );
  args.push_back(0);  // don't forget the NUL terminator!
  std::vector<char*> argv;
  char* i = std::strtok( &args[0], " \t" );
  argv.push_back(i);
  while( i != NULL ){
    i = std::strtok( NULL, " \t" );
    argv.push_back(i);
  }
  fromArgs( argv.size()-1, &argv[0], optionsOnly );
}

void LastexArguments::fromStream( std::istream& is, bool optionsOnly ){
  std::string trigger = "#last";
  for( std::string line; std::getline( is, line ); /* noop */ )
    if( line.compare( 0, trigger.size(), trigger ) == 0 )
      fromLine( line, optionsOnly );
}

void LastexArguments::fromString( const std::string& s, bool optionsOnly ){
  std::istringstream iss(s);
  fromStream( iss, optionsOnly );
}

void LastexArguments::setDefaultsFromAlphabet( bool isDna, bool isProtein ){
  if( strand < 0 ) strand = (isDna ? 2 : 1);

  if( isProtein ){
    // default match & mismatch scores: Blosum62 matrix
    if( matchScore < 0 && mismatchCost >= 0 ) matchScore   = 1;  // idiot-proof
    if( mismatchCost < 0 && matchScore >= 0 ) mismatchCost = 1;  // idiot-proof
    if( gapExistCost   < 0 ) gapExistCost   =  11;
    if( gapExtendCost  < 0 ) gapExtendCost  =   2;
  }
  else{
    if( matchScore     < 0 ) matchScore     =   1;
    if( mismatchCost   < 0 ) mismatchCost   =   1;
    if( gapExistCost   < 0 ) gapExistCost   =   7;
    if( gapExtendCost  < 0 ) gapExtendCost  =   1;
  }
}

}  // end namespace
