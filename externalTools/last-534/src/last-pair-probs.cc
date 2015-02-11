// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

#include "last-pair-probs.hh"

#include "io.hh"
#include "stringify.hh"

#include <algorithm>
#include <cctype>  // isalpha
#include <cerrno>
#include <cmath>
#include <cstdlib>  // atof
#include <cstring>  // strncmp
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits.h>
#include <cfloat>
#include <stddef.h>  // size_t

typedef const char *String;

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static bool isGraph(char c) {
  return c > ' ';  // faster than std::isgraph
}

static bool isSpace(char c) {
  return c > 0 && c <= ' ';  // faster than std::isspace
}

static bool isDigit(char c) {
  return c >= '0' && c <= '9';
}

static int wordCmp(const char* x, const char* y) {
  // Like strcmp, but stops at spaces.
  while (isGraph(*y)) {
    if (*x != *y) return *x - *y;
    ++x;
    ++y;
  }
  return isGraph(*x);
}

struct Alignment {
  const String *linesBeg;
  const String *linesEnd;
  const char *rName;
  long c;
  long rSize;
  double scaledScore;
  const char *qName;
  char strand;
  bool operator<( const Alignment& aln ) const {
    if (strand != aln.strand) return strand < aln.strand;
    int rNameCmp = wordCmp(rName, aln.rName);
    return rNameCmp < 0;
  }
};

static const char *readLong(const char *c, long &x) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  // this doesn't read negative numbers:
  if (!isDigit(*c)) return 0;
  long z = *c++ - '0';
  while (isDigit(*c)) {
    if (z > LONG_MAX / 10) return 0;
    z *= 10;
    long digit = *c++ - '0';
    if (z > LONG_MAX - digit) return 0;
    z += digit;
  }
  x = z;
  return c;
}

static const char *readDouble(const char *c, double &x) {
  if (!c) return 0;
  errno = 0;
  char *e;
  double z = std::strtod(c, &e);
  if (e == c || errno == ERANGE) return 0;
  x = z;
  return e;
}

static const char *readChar(const char *c, char &d) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  if (*c == 0) return 0;
  d = *c++;
  return c;
}

static const char *readWord(const char *c, String &s) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  const char *e = c;
  while (isGraph(*e)) ++e;
  if (e == c) return 0;
  s = c;
  return e;
}

static const char *skipWord(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  const char *e = c;
  while (isGraph(*e)) ++e;
  if (e == c) return 0;
  return e;
}

static const char *skipSpace(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  return c;
}

static double logSumExp(const double *beg, const double *end) {
  // Adds numbers, in log space, to avoid overflow.
  const double m = *std::max_element(beg, end);
  double s = 0.0;
  while (beg < end) {
    s += std::exp(*beg - m);
    ++beg;
  }
  return std::log(s) + m;
}

static double logSumExp(double x, double y) {
  // Special case of the preceding routine.
  return (x > y)
    ? std::log(1 + std::exp(y - x)) + x
    : std::log(std::exp(x - y) + 1) + y;
}

class AlignmentParameters {
  // Parses the score scale factor, minimum score, and genome size.
  double t;  // score scale factor
  double e;  // minimum score
  long   g;  // genome size

 public:
  AlignmentParameters() : t(-1), e(-1), g(-1) {}  // dummy values

  void update(const std::string& line) {
    std::istringstream ss(line);
    std::string i;
    while (ss >> i) {
      const char *c = i.c_str();
      if (t == -1.0 && i.substr(0,2) == "t=") {
	cbrc::unstringify(t, c + 2);
        if (t <= 0) err("t must be positive");
      }
      if (e == -1.0 && i.substr(0,2) == "e=") {
	cbrc::unstringify(e, c + 2);
        if (e <= 0) err("e must be positive");
      }
      if (g == -1.0 && i.substr(0,8) == "letters=") {
	cbrc::unstringify(g, c + 8);
        if (g <= 0) err("letters must be positive");
      }
    }
  }

  bool isValid() const { return t != -1 && e != -1 && g != -1; }

  void validate() const {
    if (t == -1) err("I need a header line with t=");
    if (e == -1) err("I need a header line with e=");
    if (g == -1) err("I need a header line with letters=");
  }

  double tGet() const { return t; }
  double eGet() const { return e; }
  long   gGet() const { return g; }
};

static bool isGoodQueryName(const char *nameEnd) {
  return nameEnd[-2] == '/' && (nameEnd[-1] == '1' || nameEnd[-1] == '2');
}

static void printAlignmentWithMismapProb(const Alignment& alignment,
                                         double prob, const char *suf) {
  const String *linesBeg = alignment.linesBeg;
  const String *linesEnd = alignment.linesEnd;
  const char *qName = alignment.qName;
  size_t qNameLen = skipWord(qName) - qName;
  if (isGoodQueryName(qName + qNameLen)) suf = "";
  char p[32];
  sprintf(p, "%.3g", prob);
  if (linesEnd - linesBeg == 1) {  // we have tabular format
    const char *c = *linesBeg;
    const char *d = c;
    for (int i = 0; i < 7; ++i) d = skipWord(d);
    std::cout.write(c, d - c);
    std::cout << suf << d << '\t' << p << '\n';
  } else {  // we have MAF format
    std::cout << *linesBeg << " mismap=" << p << '\n';
    const char *pad = *suf ? "  " : "";  // spacer to keep the alignment of MAF lines
    const char *rName = alignment.rName;
    size_t rNameLen = skipWord(rName) - rName;
    size_t rNameEnd = rNameLen + 2;  // where to insert the spacer
    size_t qNameEnd = qNameLen + 2;  // where to insert the suffix
    unsigned s = 0;
    for (const String *i = linesBeg + 1; i < linesEnd; ++i) {
      const char *c = *i;
      if (*c == 's' || *c == 'q') {
        if (*c == 's') s++;
        if (s == 1) {
	  std::cout.write(c, rNameEnd);
	  std::cout << pad << (c + rNameEnd) << '\n';
        } else {
	  std::cout.write(c, qNameEnd);
	  std::cout << suf << (c + qNameEnd) << '\n';
        }
      } else if (*c == 'p') {
        std::cout.write(c, 1);
        std::cout << pad << (c + 1) << '\n';
      } else {
        std::cout << c << '\n';
      }
    }
    std::cout << '\n';	// each MAF block should end with a blank line
  }
}

static long headToHeadDistance(const Alignment& x, const Alignment& y) {
  // The 5'-to-5' distance between 2 alignments on opposite strands.
  long length = x.c + y.c;
  if (length > x.rSize) length -= x.rSize;  // for circular chroms
  return length;
}

static double *conjointScores(const Alignment& aln1,
			      const Alignment *jBeg, const Alignment *jEnd,
			      double fraglen, double inner, bool isRna,
			      double *scores) {
  for (const Alignment *j = jBeg; j < jEnd; ++j) {
    const Alignment &aln2 = *j;
    if (aln1 < aln2) break;
    long length = headToHeadDistance(aln1, aln2);
    if (isRna) {  // use a log-normal distribution
      if (length <= 0) continue;
      double loglen = std::log((double)length);
      double x = loglen - fraglen;
      *scores++ = aln2.scaledScore + inner * (x * x) - loglen;
    } else {      // use a normal distribution
      if ((length > 0) != (fraglen > 0.0)) continue;  // ?
      double x = length - fraglen;
      *scores++ = aln2.scaledScore + inner * (x * x);
    }
  }
  return scores;
}

static void probForEachAlignment(const std::vector<Alignment>& alignments1,
				 const std::vector<Alignment>& alignments2,
				 const LastPairProbsOptions& opts,
				 double *zs) {
  size_t size2 = alignments2.size();
  std::vector<double> scoresVec(size2);
  double *scores = &scoresVec[0];
  for (size_t i = 0; i < size2; ++i)
    scores[i] = alignments2[i].scaledScore;
  double x = opts.disjointScore + logSumExp(scores, scores + size2);

  const Alignment *jBeg = &alignments2[0];
  const Alignment *jEnd = jBeg + size2;
  std::vector<Alignment>::const_iterator i;
  for (i = alignments1.begin(); i < alignments1.end(); ++i) {
    while (jBeg < jEnd && *jBeg < *i) ++jBeg;
    double *scoresEnd = conjointScores(*i, jBeg, jEnd, opts.fraglen,
				       opts.inner, opts.rna, scores);
    if (scoresEnd > scores) {
      double y = opts.outer + logSumExp(scores, scoresEnd);
      *zs++ = i->scaledScore + logSumExp(x, y);
    } else {
      *zs++ = i->scaledScore + x;
    }
  }
}

static void printAlnsForOneRead(const std::vector<Alignment>& alns1,
                                const std::vector<Alignment>& alns2,
                                const LastPairProbsOptions& opts,
                                double maxMissingScore, const char *suf) {
  size_t size1 = alns1.size();
  if (size1 == 0) return;
  size_t size2 = alns2.size();
  std::vector<double> zsVec(size1);
  double *zs = &zsVec[0];
  double w = maxMissingScore;

  if (size2 > 0) {
    probForEachAlignment(alns1, alns2, opts, zs);
    double w0 = -DBL_MAX;
    for (size_t i = 0; i < size2; ++i)
      w0 = std::max(w0, alns2[i].scaledScore);
    w += w0;
  } else {
    for (size_t i = 0; i < size1; ++i)
      zs[i] = alns1[i].scaledScore + opts.disjointScore;
  }

  double z = logSumExp(zs, zs + size1);
  double zw = logSumExp(z, w);

  for (size_t i = 0; i < size1; ++i) {
    double prob = 1.0 - std::exp(zs[i] - zw);
    if (prob <= opts.mismap) printAlignmentWithMismapProb(alns1[i], prob, suf);
  }
}

static void unambiguousFragmentLengths(const std::vector<Alignment>& alns1,
                                       const std::vector<Alignment>& alns2,
                                       std::vector<long>& lengths) {
  // Returns the fragment length implied by alignments of a pair of reads.
  long oldLen = LONG_MAX;
  std::vector<Alignment>::const_iterator i, j;
  std::vector<Alignment>::const_iterator jBeg = alns2.begin();
  std::vector<Alignment>::const_iterator jEnd = alns2.end();
  for (i = alns1.begin(); i < alns1.end(); ++i) {
    while (jBeg < jEnd && *jBeg < *i) ++jBeg;
    for (j = jBeg; j < jEnd; ++j) {
      if (*i < *j) break;
      long newLen = headToHeadDistance(*i, *j);
      if (oldLen == LONG_MAX) oldLen = newLen;
      else if (newLen != oldLen) return;  // the fragment length is ambiguous
    }
  }
  if (oldLen != LONG_MAX) lengths.push_back(oldLen);
}

static AlignmentParameters readHeaderOrDie(std::istream& lines) {
  std::string line;
  AlignmentParameters params;
  while (getline(lines, line)) {
    if (line.substr(0,1) == "#") {
      params.update(line);
      if (params.isValid())
        return params;
    } else if (line.find_first_not_of(" ") != std::string::npos) {
      break;
    }
  }
  params.validate();  // die
  return params;  // dummy
}

static Alignment parseAlignment(double score, const char *rName,
                                long rStart, long rSpan, long rSize,
                                const char *qName, char qStrand,
				const String *linesBeg, const String *linesEnd,
                                char strand, double scale,
                                const std::set<std::string>& circularChroms) {
  long c = -rStart;
  if (qStrand == '-') {
    c = rStart + rSpan;
    if (circularChroms.find(rName) != circularChroms.end() ||
        circularChroms.find(".") != circularChroms.end()) c += rSize;
  }

  Alignment parse;
  parse.linesBeg = linesBeg;
  parse.linesEnd = linesEnd;
  parse.rName = rName;
  parse.c = c;
  parse.rSize = rSize;
  parse.scaledScore = score / scale;  // needed in 2nd pass
  parse.qName = qName;
  parse.strand = (qStrand == strand) ? '+' : '-';
  return parse;
}

static double parseMafScore(const char *aLine) {
  const char *c = aLine;
  while ((c = skipWord(c))) {
    c = skipSpace(c);
    if (std::strncmp(c, "score=", 6) == 0) {
      double score;
      c = readDouble(c + 6, score);
      if (!c) err("bad score");
      return score;
    }
  }
  err("missing score");
  return 0.0;  // dummy;
}

static Alignment parseMaf(const String *linesBeg, const String *linesEnd,
			  char strand, double scale,
			  const std::set<std::string>& circularChroms) {
  double score = parseMafScore(*linesBeg);
  String rName, qName;
  char qStrand;
  long rStart, rSpan, rSize;
  unsigned n = 0;
  for (const String *i = linesBeg; i < linesEnd; ++i) {
    const char *c = *i;
    if (*c == 's') {
      if (n == 0) {
	c = skipWord(c);
        c = readWord(c, rName);
        c = readLong(c, rStart);
        c = readLong(c, rSpan);
        c = skipWord(c);
	c = readLong(c, rSize);
      } else if (n == 1) {
	c = skipWord(c);
	c = readWord(c, qName);
	c = skipWord(c);
	c = skipWord(c);
	c = readChar(c, qStrand);
      }
      n++;
    }
    if (!c) err("bad MAF line: " + std::string(*i));
  }
  if (n < 2) err("bad MAF");
  return parseAlignment(score, rName, rStart, rSpan, rSize, qName, qStrand,
                        linesBeg, linesEnd, strand, scale, circularChroms);
}

static Alignment parseTab(const String *linesBeg, const String *linesEnd,
			  char strand, double scale,
			  const std::set<std::string>& circularChroms) {
  double score;
  String rName, qName;
  char qStrand;
  long rStart, rSpan, rSize;
  const char *c = *linesBeg;
  c = readDouble(c, score);
  c = readWord(c, rName);
  c = readLong(c, rStart);
  c = readLong(c, rSpan);
  c = skipWord(c);
  c = readLong(c, rSize);
  c = readWord(c, qName);
  c = skipWord(c);
  c = skipWord(c);
  c = readChar(c, qStrand);
  if (!c) err("bad line: " + std::string(*linesBeg));
  return parseAlignment(score, rName, rStart, rSpan, rSize, qName, qStrand,
                        linesBeg, linesEnd, strand, scale, circularChroms);
}

static bool readBatch(std::istream& input,
		      char strand, const double scale,
		      const std::set<std::string>& circularChroms,
		      std::vector<char>& text, std::vector<String>& lines,
		      std::vector<Alignment>& alns) {
  // Yields alignment data from MAF or tabular format.
  text.clear();
  alns.clear();
  std::vector<size_t> lineStarts;
  std::string line;
  while (std::getline(input, line)) {
    const char *c = line.c_str();
    if (std::strncmp(c, "# batch ", 8) == 0) break;
    lineStarts.push_back(text.size());
    text.insert(text.end(), c, c + line.size() + 1);
  }

  size_t numOfLines = lineStarts.size();
  lines.resize(numOfLines);

  size_t mafStart = 0;
  for (size_t i = 0; i < numOfLines; ++i) {
    String *j = &lines[i];
    *j = &text[lineStarts[i]];
    const char *c = *j;
    if (isDigit(*c))
      alns.push_back(parseTab(j, j + 1, strand, scale, circularChroms));
    if (!std::isalpha(*c)) {
      if (mafStart < i)
	alns.push_back(parseMaf(&lines[mafStart], j,
				strand, scale, circularChroms));
      mafStart = i + 1;
    }
  }
  if (mafStart < numOfLines)
    alns.push_back(parseMaf(&lines[mafStart], &lines[0] + numOfLines,
			    strand, scale, circularChroms));

  size_t numOfAlns = alns.size();
  for (size_t i = 1; i < numOfAlns; ++i)
    if (wordCmp(alns[i - 1].qName, alns[i].qName) != 0)
      err("found 2 queries in 1 batch: did you forget lastal -i1?");

  stable_sort(alns.begin(), alns.end());

  return input;
}

static std::vector<long> readQueryPairs1pass(std::istream& in1, std::istream& in2,
                                             double scale1, double scale2,
                                             const std::set<std::string>& circularChroms) {
  std::vector<long> lengths;
  std::vector<char> text1, text2;
  std::vector<String> lines1, lines2;
  std::vector<Alignment> a1, a2;
  while (1) {
    bool ok1 = readBatch(in1, '+', scale1, circularChroms, text1, lines1, a1);
    bool ok2 = readBatch(in2, '-', scale2, circularChroms, text2, lines2, a2);
    unambiguousFragmentLengths(a1, a2, lengths);
    if (!ok1 || !ok2) break;
  }
  return lengths;
}

static void readQueryPairs2pass(std::istream& in1, std::istream& in2,
                                double scale1, double scale2,
                                const LastPairProbsOptions& opts) {
  std::vector<char> text1, text2;
  std::vector<String> lines1, lines2;
  std::vector<Alignment> a1, a2;
  while (1) {
    bool ok1 = readBatch(in1, '+', scale1, opts.circular, text1, lines1, a1);
    bool ok2 = readBatch(in2, '-', scale2, opts.circular, text2, lines2, a2);
    printAlnsForOneRead(a1, a2, opts, opts.maxMissingScore1, "/1");
    printAlnsForOneRead(a2, a1, opts, opts.maxMissingScore2, "/2");
    if (!ok1 || !ok2) break;
  }
}

static double myRound(double myFloat) {
  char buf[32];
  sprintf(buf, "%g", myFloat);
  return std::atof(buf);
}

static void estimateFragmentLengthDistribution(std::vector<long> lengths,
                                               LastPairProbsOptions& opts) {
  if (lengths.empty())
    err("can't estimate the distribution of distances");

  // Define quartiles in the most naive way possible:
  std::sort(lengths.begin(), lengths.end());
  const size_t sampleSize = lengths.size();
  const long quartile1 = lengths[sampleSize / 4];
  const long quartile2 = lengths[sampleSize / 2];
  const long quartile3 = lengths[sampleSize * 3 / 4];

  std::cout << "# distance sample size: " << sampleSize << "\n";
  std::cout << "# distance quartiles: "
	    << quartile1 << " " << quartile2 << " " << quartile3 << "\n";

  if (opts.rna && quartile1 <= 0)
    err("too many distances <= 0");

  const char* thing = (opts.rna) ? "ln[distance]" : "distance";

  if (!opts.isFraglen) {
    if (opts.rna) opts.fraglen = myRound(std::log((double)quartile2));
    else          opts.fraglen = double(quartile2);
    std::cout << "# estimated mean " << thing << ": " << opts.fraglen << "\n";
  }

  if (!opts.isSdev) {
    const double iqr =
      (opts.rna) ? std::log((double)quartile3) - std::log((double)quartile1)
      :            quartile3 - quartile1;
    // Normal Distribution: sdev = iqr / (2 * qnorm(0.75))
    opts.sdev = myRound(iqr / 1.34898);
    std::cout << "# estimated standard deviation of " << thing << ": "
	      << opts.sdev << "\n";
  }
}

static double safeLog(const double x) {
  if (x == 0.0) return -1.0e99;
  else          return std::log(x);
}

static void calculateScorePieces(LastPairProbsOptions& opts,
                                 const AlignmentParameters& params1,
                                 const AlignmentParameters& params2) {
  if (opts.sdev == 0.0) {
    opts.outer = opts.rna ? opts.fraglen : 0.0;
    opts.inner = -1.0e99;
  } else {  // parameters for a Normal Distribution (of fragment lengths):
    const double pi = atan(1.0) * 4.0;
    opts.outer = -std::log(opts.sdev * std::sqrt(2.0 * pi));
    opts.inner = -1.0 / (2.0 * (opts.sdev * opts.sdev));
  }
  opts.outer += safeLog(1.0 - opts.disjoint);

  if (params1.gGet() != params2.gGet()) err("unequal genome sizes");
  // Multiply genome size by 2, because it has 2 strands:
  opts.disjointScore = safeLog(opts.disjoint) - std::log(2.0 * params1.gGet());

  // Max possible influence of an alignment just below the score threshold:
  double maxLogPrior = opts.outer;
  if (opts.rna) maxLogPrior += (opts.sdev * opts.sdev) / 2.0 - opts.fraglen;
  opts.maxMissingScore1 = (params1.eGet() - 1) / params1.tGet() + maxLogPrior;
  opts.maxMissingScore2 = (params2.eGet() - 1) / params2.tGet() + maxLogPrior;
}

static void skipOneBatchMarker(std::istream& input) {
  std::string line;
  while (std::getline(input, line))
    if (line.find("# batch ") == 0) break;
}

void lastPairProbs(LastPairProbsOptions& opts) {
  const std::vector<std::string>& inputs = opts.inputFileNames;
  size_t n = inputs.size();

  if (!opts.isFraglen || !opts.isSdev) {
    std::ifstream inFile1, inFile2;
    std::istream& in1 = (n > 0) ? cbrc::openIn(inputs[0], inFile1) : std::cin;
    std::istream& in2 = (n > 1) ? cbrc::openIn(inputs[1], inFile2) : in1;
    if (n < 2) skipOneBatchMarker(in1);
    std::vector<long> lengths = readQueryPairs1pass(in1, in2, 1.0, 1.0,
						    opts.circular);
    estimateFragmentLengthDistribution(lengths, opts);
  }

  if (!opts.estdist) {
    std::ifstream inFile1, inFile2;
    std::istream& in1 = (n > 0) ? cbrc::openIn(inputs[0], inFile1) : std::cin;
    std::istream& in2 = (n > 1) ? cbrc::openIn(inputs[1], inFile2) : in1;
    AlignmentParameters params1 = readHeaderOrDie(in1);
    AlignmentParameters params2 = (n > 1) ? readHeaderOrDie(in2) : params1;
    calculateScorePieces(opts, params1, params2);
    std::cout << "# fraglen=" << opts.fraglen
              << " sdev=" << opts.sdev
              << " disjoint=" << opts.disjoint
              << " genome=" << params1.gGet() << "\n";
    if (n < 2) skipOneBatchMarker(in1);
    readQueryPairs2pass(in1, in2, params1.tGet(), params2.tGet(), opts);
  }
}
