// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

// A "subset seed" covers a range of sequence.  The size of this range
// is called its span.  At each position, it maps letters (encoded as
// numbers) to subsets of the letters (encoded as numbers).  The
// mapping may or may not be different at different positions.

// Subset seeds are described in:
// G Kucherov et al. J Bioinform Comput Biol. 2006 4(2):553-69.

// "Cyclic" means that the seed can cover arbitrary-size ranges, by
// cyclically repeating.

// The INPUT/OUTPUT FORMAT looks like this:
//   A C G T
//   AG CT
//   ACGT
// Each line represents one position, and each letter-group (separated
// by whitespace) is one subset.  In each position, any letters not
// specified will get mapped to the special DELIMITER subset.  Blank
// lines and comment lines starting with # are ignored.

// The input format is always case-insensitive.  If the
// isMaskLowercase argument of the reading routines is true, then all
// lowercase letters will get mapped to the DELIMITER subset,
// otherwise they will be treated like their uppercase equivalents.

#ifndef CYCLIC_SUBSET_SEED_HH
#define CYCLIC_SUBSET_SEED_HH

#include <string>
#include <vector>
#include <iosfwd>

namespace cbrc{

typedef unsigned char uchar;

class CyclicSubsetSeed{
public:
  enum { MAX_LETTERS = 64 };
  enum { DELIMITER = 255 };

  // This converts a seed name to a set of strings in the I/O format.
  // If the name isn't known, it assumes it's a file and tries to read it.
  static std::vector<std::string> fromName( const std::string& name );

  // This converts a mask to a set of strings in the I/O format.
  // The mask is something like: "1110TT,1001T1".  The "1"s are
  // must-match positions, the "0"s are don't care positions, and "T"
  // or "t" allows transitions but not transversions.  For consistency
  // with YASS/Iedera, you can also use "#" for match, "@" for
  // transition, and "_" or "-" for don't care.  You can have multiple
  // seed patterns separated by commas, which get returned as separate
  // strings.  The alph indicates the letters that can occur in "1"
  // and "0" positions.
  static std::vector<std::string> fromMask( const std::string& alph,
					    const std::string& mask );

  // This is a more general version, where seedAlph maps mask
  // characters to lines in the I/O format.
  static std::vector<std::string> fromMask( const char* seedAlph[],
					    const std::string& mask );

  void clear() { subsetLists.clear(); subsetMaps.clear(); }

  void fromString( const std::string& s,
		   bool isMaskLowercase, const uchar letterCode[] );

  void fromStream( std::istream& stream,
		   bool isMaskLowercase, const uchar letterCode[] );

  void appendPosition( std::istream& inputLine,
		       bool isMaskLowercase, const uchar letterCode[] );

  void writePosition( std::ostream& out, unsigned position ) const;

  unsigned span() const{
    return subsetLists.size();
  }

  const uchar* subsetMap( unsigned depth ) const{
    return &subsetMaps[0] + (depth % span()) * MAX_LETTERS;
  }

  unsigned subsetCount( unsigned depth ) const{
    return subsetLists[ depth % span() ].size();
  }

  const uchar* firstMap() const{
    return &subsetMaps[0];
  }

  const uchar* nextMap( const uchar* x ) const{
    const uchar* y = x + MAX_LETTERS;
    if( y == &subsetMaps[0] + subsetMaps.size() )
      y = &subsetMaps[0];
    return y;
  }

private:
  std::vector< std::vector<std::string> > subsetLists;
  std::vector<uchar> subsetMaps;

  static void addLetter( std::vector<uchar>& numbersToSubsets,
			 uchar letter, uchar subsetNum,
			 const uchar letterCode[] );
};

}  // end namespace
#endif
