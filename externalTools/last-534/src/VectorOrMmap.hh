// Copyright 2010, 2011 Martin C. Frith

// Container that holds its data in either a std::vector or an Mmap.

#ifndef VECTOR_OR_MMAP_HH
#define VECTOR_OR_MMAP_HH

#include "Mmap.hh"

#include <vector>

namespace cbrc{

template<typename T>
struct VectorOrMmap{
  std::vector<T> v;
  Mmap<T> m;

  const T* begin() const { return v.empty() ? m.begin() : &v.front();    }
  const T* end()   const { return v.empty() ? m.end()   : &v.back() + 1; }

  std::size_t size() const { return v.empty() ? m.size() : v.size(); }

  bool empty() const { return v.empty() && m.empty(); }

  const T& front() const { return v.empty() ? m.front() : v.front(); }
  const T& back()  const { return v.empty() ? m.back()  : v.back();  }

  const T& operator[](std::size_t i) const { return v.empty() ? m[i] : v[i]; }
};

}  // end namespace

#endif
