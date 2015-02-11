// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

// Parse the command line and run last-pair-probs.

#include "last-pair-probs.hh"
#include "stringify.hh"

#include <getopt.h>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <new>  // bad_alloc

static void run(int argc, char* argv[]) {
  LastPairProbsOptions opts;

  opts.rna = false;
  opts.estdist = false;
  opts.mismap = 0.01;
  opts.isFraglen = false;
  opts.isSdev = false;
  opts.isDisjoint = false;

  const char *version = "last-pair-probs "
#include "version.hh"
"\n";

  std::string help = "\
Usage:\n\
  " + std::string(argv[0]) + " --help\n\
  " + std::string(argv[0]) + " [options] interleaved-alignments\n\
  " + std::string(argv[0]) + " [options] alignments1 alignments2\n\
\n\
Read alignments of paired DNA reads to a genome, and: (1) estimate the\n\
distribution of distances between paired reads, (2) estimate the probability\n\
that each alignment represents the genomic source of the read.\n\
\n\
Options:\n\
  -h, --help            show this help message and exit\n\
  -r, --rna             assume the reads are from potentially-spliced RNA\n\
  -e, --estdist         just estimate the distribution of distances\n\
  -m M, --mismap=M      don't write alignments with mismap probability > M\n\
                        (default: " + cbrc::stringify(opts.mismap) + ")\n\
  -f BP, --fraglen=BP   mean distance in bp\n\
  -s BP, --sdev=BP      standard deviation of distance\n\
  -d PROB, --disjoint=PROB\n\
                        prior probability of disjoint mapping (default: 0.02\n\
                        if -r, else 0.01)\n\
  -c CHROM, --circular=CHROM\n\
                        specifies that chromosome CHROM is circular (default:\n\
                        chrM)\n\
  -V, --version         show program's version number and exit\n\
";

  const char sOpts[] = "hrem:f:s:d:c:V";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "rna",      no_argument,       0, 'r' },
    { "estdist",  no_argument,       0, 'e' },
    { "mismap",   required_argument, 0, 'm' },
    { "fraglen",  required_argument, 0, 'f' },
    { "sdev",     required_argument, 0, 's' },
    { "disjoint", required_argument, 0, 'd' },
    { "circular", required_argument, 0, 'c' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return;
    case 'r':
      opts.rna = true;
      break;
    case 'e':
      opts.estdist = true;
      break;
    case 'm':
      cbrc::unstringify(opts.mismap, optarg);
      break;
    case 'f':
      opts.isFraglen = true;
      cbrc::unstringify(opts.fraglen, optarg);
      break;
    case 's':
      opts.isSdev = true;
      cbrc::unstringify(opts.sdev, optarg);
      if (opts.sdev < 0.0) {
        throw std::runtime_error("option -s: should be >= 0");
      }
      break;
    case 'd':
      opts.isDisjoint = true;
      cbrc::unstringify(opts.disjoint, optarg);
      if (opts.disjoint < 0.0) {
        throw std::runtime_error("option -d: should be >= 0");
      } else if (opts.disjoint > 1.0) {
        throw std::runtime_error("option -d: should be <= 1");
      }
      break;
    case 'c':
      opts.circular.insert(optarg);
      break;
    case 'V':
      std::cout << version;
      return;
    case '?':
      throw std::runtime_error("");
    }
  }

  if (optind == argc && !opts.estdist && (!opts.isFraglen || !opts.isSdev)) {
    std::cerr << help;
    throw std::runtime_error("");
  }

  if (!opts.isDisjoint) {
    opts.disjoint = opts.rna ? 0.02 : 0.01;
  }

  opts.inputFileNames.assign(argv + optind, argv + argc);

  if (opts.inputFileNames.size() > 2) {
    throw std::runtime_error("too many file names");
  }

  if (!opts.circular.size()) {
      opts.circular.insert("chrM");
  }
  std::ios_base::sync_with_stdio(false);  // makes std::cin much faster!!!

  lastPairProbs(opts);
}

int main(int argc, char* argv[]) {
  try {
    run(argc, argv);
    if (!flush(std::cout)) throw std::runtime_error("write error");
    return EXIT_SUCCESS;
  } catch (const std::bad_alloc& e) {  // bad_alloc::what() may be unfriendly
    std::cerr << argv[0] << ": out of memory\n";
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    const char *s = e.what();
    if (*s) std::cerr << argv[0] << ": " << s << '\n';
    return EXIT_FAILURE;
  }
}
