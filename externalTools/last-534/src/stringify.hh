// Copyright 2008, 2009 Martin C. Frith

// Convert things (e.g. numbers) to and from strings.  These are
// replacements for boost::lexical_cast: by avoiding dependency on
// boost, we can distribute code more easily.

#ifndef STRINGIFY_HH
#define STRINGIFY_HH
#include <sstream>
#include <string>
#include <stdexcept>
#include <cassert>

namespace cbrc{

template<typename T>
std::string stringify( const T& x ){
  std::ostringstream oss;
  oss << x;
  assert(oss);
  return oss.str();
}

template<typename T>
void unstringify( T& x, const std::string& s ){
  std::istringstream iss(s);
  if( !(iss >> x) || !(iss >> std::ws).eof() ){
    throw std::runtime_error( "can't interpret: " + s );
  }
}

template<typename T>
void unstringifySize( T& x, const std::string& s ){
  std::istringstream iss(s);
  if( !(iss >> x) ){
    throw std::runtime_error( "can't interpret: " + s );
  }

  std::string suffix;
  if( iss >> suffix ){
    int multiplier;
    /**/ if( suffix == "K" ) multiplier = 1024;  // "KibiBytes"
    else if( suffix == "M" ) multiplier = 1024*1024;  // "MebiBytes"
    else if( suffix == "G" ) multiplier = 1024*1024*1024;  // "GibiBytes"
    else throw std::runtime_error( "can't interpret: " + s );
    if( (x * multiplier) / multiplier != x ){  // check for overflow
      throw std::runtime_error( "can't interpret (too big): " + s );
    }
    x *= multiplier;
  }

  if( !(iss >> std::ws).eof() ){
    throw std::runtime_error( "can't interpret: " + s );
  }
}

}  // end namespace cbrc
#endif  // STRINGIFY_HH
