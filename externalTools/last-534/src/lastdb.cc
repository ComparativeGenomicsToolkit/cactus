// Copyright 2008, 2009, 2010, 2011, 2013, 2014 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "LastdbArguments.hh"
#include "SubsetSuffixArray.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "io.hh"
#include "qualityScoreUtil.hh"
#include "stringify.hh"
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <numeric>  // accumulate

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastdb: " << x << '\n'

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, const LastdbArguments& args ){
  if( !args.userAlphabet.empty() )  alph.fromString( args.userAlphabet );
  else if( args.isProtein )         alph.fromString( alph.protein );
  else                              alph.fromString( alph.dna );
}

// Does the first sequence look like it isn't really DNA?
bool isDubiousDna( const Alphabet& alph, const MultiSequence& multi ){
  const uchar* seq = multi.seqReader() + multi.seqBeg(0);
  unsigned dnaCount = 0;

  for( indexT i = 0; i < 100; ++i ){  // look at the first 100 letters
    uchar c = alph.numbersToUppercase[ seq[i] ];
    if( c == alph.size ) return false;  // we hit the end of the sequence early
    if( c < alph.size || c == alph.encode[ (uchar)'N' ] ) ++dnaCount;
  }

  if( dnaCount < 90 ) return true;  // more than 10% unexpected letters
  else return false;
}

const unsigned maxNumOfIndexes = 16;

static void addSeeds( SubsetSuffixArray indexes[], unsigned& numOfIndexes,
		      const std::vector<std::string>& seedStrings,
		      const LastdbArguments& args, const Alphabet& alph ){
  for( unsigned x = 0; x < seedStrings.size(); ++x ){
    if( numOfIndexes >= maxNumOfIndexes ) ERR( "too many seed patterns" );
    CyclicSubsetSeed& seed = indexes[ numOfIndexes++ ].getSeed();
    seed.fromString( seedStrings[x], args.isCaseSensitive, alph.encode );
  }
}

// Set up the seed pattern(s), and return how many of them there are
unsigned makeSubsetSeeds( SubsetSuffixArray indexes[],
			  const LastdbArguments& args, const Alphabet& alph ){
  unsigned numOfIndexes = 0;
  const std::string& a = alph.letters;

  for( unsigned x = 0; x < args.subsetSeedFiles.size(); ++x ){
    const std::string& name = args.subsetSeedFiles[x];
    std::vector<std::string> s = CyclicSubsetSeed::fromName( name );
    addSeeds( indexes, numOfIndexes, s, args, alph );
  }

  for( unsigned x = 0; x < args.seedPatterns.size(); ++x ){
    const std::string& mask = args.seedPatterns[x];
    std::vector<std::string> s = CyclicSubsetSeed::fromMask( a, mask );
    addSeeds( indexes, numOfIndexes, s, args, alph );
  }

  if( numOfIndexes == 0 ){
    if( alph.letters == alph.dna ){
      const char* mask = "1T1001100101";  // YASS
      std::vector<std::string> s = CyclicSubsetSeed::fromMask( a, mask );
      addSeeds( indexes, numOfIndexes, s, args, alph );
    }
    else{
      std::vector<std::string> s = CyclicSubsetSeed::fromMask( a, "1" );
      addSeeds( indexes, numOfIndexes, s, args, alph );
    }
  }

  return numOfIndexes;
}

void writePrjFile( const std::string& fileName, const LastdbArguments& args,
		   const Alphabet& alph, countT sequenceCount,
		   const std::vector<countT>& letterCounts,
		   unsigned volumes, unsigned numOfIndexes ){
  countT letterTotal = std::accumulate( letterCounts.begin(),
                                        letterCounts.end(), countT(0) );

  std::ofstream f( fileName.c_str() );
  f << "version=" <<
#include "version.hh"
    << '\n';
  f << "alphabet=" << alph << '\n';
  f << "numofsequences=" << sequenceCount << '\n';
  f << "numofletters=" << letterTotal << '\n';
  f << "letterfreqs=";
  for( unsigned i = 0; i < letterCounts.size(); ++i ){
    if( i > 0 ) f << ' ';
    f << letterCounts[i];
  }
  f << '\n';

  if( !args.isCountsOnly ){
    f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
    f << "masklowercase=" << args.isCaseSensitive << '\n';
    if( args.inputFormat != sequenceFormat::fasta ){
      f << "sequenceformat=" << args.inputFormat << '\n';
    }
    if( volumes+1 > 0 ){
      f << "volumes=" << volumes << '\n';
    }
    else{
      f << "numofindexes=" << numOfIndexes << '\n';
    }
  }

  if( !f ) ERR( "can't write file: " + fileName );
}

// Make one database volume, from one batch of sequences
void makeVolume( SubsetSuffixArray indexes[], unsigned numOfIndexes,
		 const MultiSequence& multi, const LastdbArguments& args,
		 const Alphabet& alph, const std::vector<countT>& letterCounts,
		 const std::string& baseName ){
  LOG( "writing..." );
  writePrjFile( baseName + ".prj", args, alph, multi.finishedSequences(),
		letterCounts, -1, numOfIndexes );
  multi.toFiles( baseName );

  for( unsigned x = 0; x < numOfIndexes; ++x ){
    LOG( "sorting..." );
    indexes[x].sortIndex( multi.seqReader(), args.minSeedLimit );

    LOG( "bucketing..." );
    indexes[x].makeBuckets( multi.seqReader(), args.bucketDepth );

    LOG( "writing..." );
    indexT textLength = multi.finishedSize();
    if( numOfIndexes > 1 ){
      indexes[x].toFiles( baseName + char('a' + x), false, textLength );
    }
    else{
      indexes[x].toFiles( baseName, true, textLength );
    }

    indexes[x].clearPositions();
  }

  LOG( "done!" );
}

// The max number of sequence letters, such that the total volume size
// is likely to be less than volumeSize bytes.  (This is crude, it
// neglects memory for the sequence names, and the fact that
// lowercase-masked letters and DNA "N"s aren't indexed.)
static indexT maxLettersPerVolume( const LastdbArguments& args,
				   unsigned numOfIndexes ){
  size_t bytesPerLetter = isFastq( args.inputFormat ) ? 2 : 1;
  size_t maxIndexBytesPerPosition = sizeof(indexT) + 1;
  maxIndexBytesPerPosition *= numOfIndexes;
  size_t x = bytesPerLetter * args.indexStep + maxIndexBytesPerPosition;
  size_t y = args.volumeSize / x * args.indexStep;
  indexT z = y;
  if( z < y ) z = indexT(-1);
  return z;
}

// Read the next sequence, adding it to the MultiSequence and the SuffixArray
std::istream&
appendFromFasta( MultiSequence& multi, unsigned numOfIndexes,
		 const LastdbArguments& args, const Alphabet& alph,
		 std::istream& in ){
  indexT maxSeqLen = maxLettersPerVolume( args, numOfIndexes );
  if( multi.finishedSequences() == 0 ) maxSeqLen = indexT(-1);

  size_t oldSize = multi.unfinishedSize();

  if ( args.inputFormat == sequenceFormat::fasta )
    multi.appendFromFasta( in, maxSeqLen );
  else
    multi.appendFromFastq( in, maxSeqLen );

  if( !multi.isFinished() && multi.finishedSequences() == 0 )
    ERR( "encountered a sequence that's too long" );

  // encode the newly-read sequence
  uchar* seq = multi.seqWriter();
  size_t newSize = multi.unfinishedSize();
  alph.tr( seq + oldSize, seq + newSize );

  if( isPhred( args.inputFormat ) )  // assumes one quality code per letter:
    checkQualityCodes( multi.qualityReader() + oldSize,
                       multi.qualityReader() + newSize,
                       qualityOffset( args.inputFormat ) );

  return in;
}

void lastdb( int argc, char** argv ){
  LastdbArguments args;
  args.fromArgs( argc, argv );
  Alphabet alph;
  MultiSequence multi;
  SubsetSuffixArray indexes[maxNumOfIndexes];
  makeAlphabet( alph, args );
  unsigned numOfIndexes = makeSubsetSeeds( indexes, args, alph );
  multi.initForAppending(1);
  alph.tr( multi.seqWriter(), multi.seqWriter() + multi.unfinishedSize() );
  unsigned volumeNumber = 0;
  countT sequenceCount = 0;
  std::vector<countT> letterCounts( alph.size );
  std::vector<countT> letterTotals( alph.size );

  char defaultInputName[] = "-";
  char* defaultInput[] = { defaultInputName, 0 };
  char** inputBegin = argv + args.inputStart;

  for( char** i = *inputBegin ? inputBegin : defaultInput; *i; ++i ){
    std::ifstream inFileStream;
    std::istream& in = openIn( *i, inFileStream );
    LOG( "reading " << *i << "..." );

    while( appendFromFasta( multi, numOfIndexes, args, alph, in ) ){
      if( !args.isProtein && args.userAlphabet.empty() &&
          sequenceCount == 0 && isDubiousDna( alph, multi ) ){
        std::cerr << "lastdb: that's some funny-lookin DNA\n";
      }

      if( multi.isFinished() ){
        ++sequenceCount;
	const uchar* seq = multi.seqReader();
	size_t lastSeq = multi.finishedSequences() - 1;
	size_t beg = multi.seqBeg( lastSeq );
	size_t end = multi.seqEnd( lastSeq );
	alph.count( seq + beg, seq + end, &letterCounts[0] );
	if( args.isCountsOnly ){
	  // memory-saving, which seems to be important on 32-bit systems:
	  multi.reinitForAppending();
	}else{
	  for( unsigned x = 0; x < numOfIndexes; ++x ){
	    indexes[x].addPositions( seq, beg, end, args.indexStep );
	  }
	}
      }
      else{
	std::string baseName = args.lastdbName + stringify(volumeNumber++);
	makeVolume( indexes, numOfIndexes,
		    multi, args, alph, letterCounts, baseName );
	for( unsigned c = 0; c < alph.size; ++c )
	  letterTotals[c] += letterCounts[c];
	letterCounts.assign( alph.size, 0 );
	multi.reinitForAppending();
      }
    }
  }

  if( multi.finishedSequences() > 0 ){
    if( volumeNumber == 0 ){
      makeVolume( indexes, numOfIndexes,
		  multi, args, alph, letterCounts, args.lastdbName );
      return;
    }
    std::string baseName = args.lastdbName + stringify(volumeNumber++);
    makeVolume( indexes, numOfIndexes,
		multi, args, alph, letterCounts, baseName );
  }

  for( unsigned c = 0; c < alph.size; ++c ) letterTotals[c] += letterCounts[c];

  writePrjFile( args.lastdbName + ".prj", args, alph,
		sequenceCount, letterTotals, volumeNumber, numOfIndexes );
}

int main( int argc, char** argv )
try{
  lastdb( argc, argv );
  return EXIT_SUCCESS;
}
catch( const std::bad_alloc& e ) {  // bad_alloc::what() may be unfriendly
  std::cerr << "lastdb: out of memory\n";
  return EXIT_FAILURE;
}
catch( const std::exception& e ) {
  std::cerr << "lastdb: " << e.what() << '\n';
  return EXIT_FAILURE;
}
catch( int i ) {
  return i;
}
