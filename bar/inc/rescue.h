#ifndef RESCUE_H_
#define RESCUE_H_
#include "stPinchGraphs.h"

typedef struct {
    Name name;
    int64_t start;
    int64_t stop;
} bedRegion;

bedRegion *bedRegion_construct(Name name, int64_t start, int64_t stop);

// Compare two bed regions in their little-endian format as mapped
// from the file. Returns 0 for any overlap.
int bedRegion_cmp(const bedRegion *region1, const bedRegion *region2);

// Find any regions covered by outgroups that are in segments with no
// block, and "rescue" them into single-degree blocks.
void rescueCoveredRegions(stPinchThread *thread, bedRegion *beds, size_t numBeds,
                          Name name, int64_t minSegmentLength, double coveredBasesThreshold);

#endif // RESCUE_H_
