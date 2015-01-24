/* Copyright 2014 Martin C. Frith */

#include <getopt.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char progName[] = "last-merge-batches";

const char version[] =
#include "version.hh"
  ;

static void *mallocOrDie(size_t size) {
  void *p = malloc(size);
  if (!p && size) {
    fprintf(stderr, "%s: out of memory\n", progName);
    exit(EXIT_FAILURE);
  }
  return p;
}

static FILE *openOrDie(const char *fileName) {
  if (strcmp(fileName, "-") == 0) {
    return stdin;
  } else {
    FILE *f = fopen(fileName, "r");
    if (!f) {
      fprintf(stderr, "%s: can't open file: %s\n", progName, fileName);
      exit(EXIT_FAILURE);
    }
    return f;
  }
}

static char *readOrDie(char *buffer, int size, FILE *f, const char *fileName) {
  char *s = fgets(buffer, size, f);
  if (ferror(f)) {
    fprintf(stderr, "%s: can't read file: %s\n", progName, fileName);
    exit(EXIT_FAILURE);
  }
  return s;
}

static void writeOrDie(const char *s) {
  if (fputs(s, stdout) < 0) {
    fprintf(stderr, "%s: can't write the output\n", progName);
    exit(EXIT_FAILURE);
  }
}

static void flushOrDie(void) {
  if (fclose(stdout) < 0) {
    fprintf(stderr, "%s: can't finish writing the output\n", progName);
    exit(EXIT_FAILURE);
  }
}

static void lastMergeBatches(int fileNum, char **fileNames) {
  FILE **files = mallocOrDie(fileNum * sizeof *files);
  char buffer[4096];
  const int buflen = sizeof buffer;
  int isActive, i;

  for (i = 0; i < fileNum; ++i) {
    files[i] = openOrDie(fileNames[i]);
  }

  do {
    isActive = 0;
    for (i = 0; i < fileNum; ++i) {
      while (readOrDie(buffer, buflen, files[i], fileNames[i])) {
	if (memcmp(buffer, "# batch", 7) == 0) {
	  if (i+1 == fileNum) writeOrDie(buffer);
	  isActive = 1;
	  break;
	}
	writeOrDie(buffer);
      }
    }
  } while (isActive);
}

int main(int argc, char **argv) {
  const char sOpts[] = "hV";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      printf("\
Usage: %s files\n\
\n\
Read files of lastal output, merge corresponding batches, and write them.\n\
\n\
Options:\n\
 -h, --help         show this help message and exit\n\
 -V, --version      show version information and exit\n\
", argv[0]);
      return EXIT_SUCCESS;
    case 'V':
      printf("\
%s %s\n\
License GPLv3+: GNU GPL version 3 or later.\n\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.\n\
", progName, version);
      return EXIT_SUCCESS;
    case '?':
      return EXIT_FAILURE;
    }
  }

  lastMergeBatches(argc - optind, argv + optind);
  flushOrDie();
  return EXIT_SUCCESS;
}
