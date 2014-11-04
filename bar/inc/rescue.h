#ifndef RESCUE_H_
#define RESCUE_H_
#include "stPinchGraphs.h"

typedef struct {
    Name name;
    int64_t start;
    int64_t stop;
} bedRegion;

// reads the next bed line from the file. Assumes tab-delimited bed3 input.
int readNextBedLine(FILE *bedFile, bedRegion *curBedLine);

// Returns a hash of sequence Name => coverage array.
// Not thread-safe.
stHash *getIngroupCoverage(const char *ingroupCoverageDir, Flower *flower);

// Find any regions covered by outgroups that are in segments with no
// block, and "rescue" them into single-degree blocks.
void rescueCoveredRegions(stPinchThread *thread, FILE *bedFile,
                          bedRegion *curBedLine);
#endif // RESCUE_H_
