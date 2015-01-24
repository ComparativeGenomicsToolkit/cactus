// Copyright 2011 Martin C. Frith

#ifndef SEQUENCE_FORMAT_HH
#define SEQUENCE_FORMAT_HH

#include <istream>

namespace cbrc {

namespace sequenceFormat {
enum Enum { fasta, fastqSanger, fastqSolexa, fastqIllumina, prb, pssm };
}

inline std::istream &operator>>(std::istream &s, sequenceFormat::Enum &f) {
  int i = 0;
  s >> i;
  if (i < 0 || i > sequenceFormat::pssm) s.setstate(std::ios::failbit);
  if (s) f = static_cast<sequenceFormat::Enum>(i);
  return s;
}

inline bool isQuality(sequenceFormat::Enum f) {
  return f != sequenceFormat::fasta && f != sequenceFormat::pssm;
}

inline bool isFastq(sequenceFormat::Enum f) {
  return
      f == sequenceFormat::fastqSanger ||
      f == sequenceFormat::fastqSolexa ||
      f == sequenceFormat::fastqIllumina;
}

inline int qualityOffset(sequenceFormat::Enum f) {
  return (f == sequenceFormat::fastqSanger) ? 33 : 64;
}  // The result is meaningless for non-quality formats.

inline bool isPhred(sequenceFormat::Enum f) {
  return
      f == sequenceFormat::fastqSanger || f == sequenceFormat::fastqIllumina;
}

}

#endif
