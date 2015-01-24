// Copyright 2010 Martin C. Frith

// Container for a memory-mapped file.  Only reading is allowed, not
// writing.

// Since memory-mapping requires non-standard C++, the code that
// actually does that is isolated in another file.

#ifndef MMAP_HH
#define MMAP_HH

#include <algorithm>  // swap
#include <stdexcept>

#include "fileMap.hh"
#include "stringify.hh"

namespace cbrc{

template<typename T>
class Mmap{
public:
  // Make an empty Mmap.
  Mmap() : begin_(0), end_(0) {}

  // Make an Mmap with s items from a file.  Throws an exception if it
  // fails to read the file.  If s is zero, it doesn't try to read the
  // file.
  Mmap( const std::string& fileName, std::size_t s )
    : begin_(0), end_(0) { open( fileName, s ); };

  // Release the mapping, if not empty.
  ~Mmap() { close(); }

  // Map s items from a file.  Throws an exception if it fails to read
  // the file.  If s is zero, it doesn't try to read the file.
  // If a file is already being mapped, closes it first.
  void open( const std::string& fileName, std::size_t s );

  // Release the mapping, if not empty.  This makes the Mmap empty.
  void close();

  // Standard functions to access the stuff in the container.
  const T* begin() const { return begin_; }
  const T* end() const { return end_; }
  std::size_t size() const { return end_ - begin_; }
  bool empty() const { return end_ == 0; }
  const T& front() const { return *begin_; }
  const T& back() const { return *(end_ - 1); }
  const T& operator[]( std::size_t i ) const { return begin_[i]; }

  void swap( Mmap& m ){
    std::swap( begin_, m.begin_ );
    std::swap( end_, m.end_ );
  }

private:
  T* begin_;
  T* end_;

  // prevent copying:
  Mmap( const Mmap& );
  Mmap& operator=( const Mmap& );
};

template<typename T>
void Mmap<T>::open( const std::string& fileName, std::size_t s ){
  close();

  std::size_t bytes = s * sizeof(T);

  if( bytes / sizeof(T) < s )  // check for overflow
    throw std::runtime_error( "can't map " + stringify(s) +
                              " items of size " + stringify( sizeof(T) ) +
                              " (from file " + fileName + ")" );

  void* m = openFileMap( fileName, bytes );

  begin_ = static_cast<T*>(m);
  end_ = begin_ + s;
}

template<typename T>
void Mmap<T>::close(){
  std::size_t bytes = size() * sizeof(T);
  closeFileMap( begin_, bytes );
  begin_ = 0;
  end_ = 0;
}

}  // end namespace

#endif
