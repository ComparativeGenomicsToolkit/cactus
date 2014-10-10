#ifndef RESCUE_H_
#define RESCUE_H_
#include "stPinchGraphs.h"

// Returns a hash of sequence Name => coverage array.
// Not thread-safe.
stHash *getIngroupCoverage(const char *ingroupCoverageDir, Flower *flower);

// Find any regions covered by outgroups that are in segments with no
// block, and "rescue" them into single-degree blocks.
void rescueCoveredRegions(stPinchThread *thread, bool *coverageArray);
#endif // RESCUE_H_
