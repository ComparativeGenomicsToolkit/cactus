// "bar rescue" -- ensures that any ingroup sequence that found an
// outgroup alignment in the blast stage still makes it into the
// ancestor after the bar stage.

#include "cactus.h"
#include "sonLib.h"
#include "stPinchGraphs.h"

typedef struct {
    Name name; // sequence Name, since the cap Name typically used
               // isn't easily accessible from flowers further down in
               // the hierarchy.
    int64_t start; // 0-based start, inclusive.
    int64_t stop; // 0-based end, exclusive.
} bedRegion;

// Compare two bed regions in their little-endian format as mapped
// from the file. Returns 0 for any overlap.
int bedRegion_cmp(const bedRegion *region1, const bedRegion *region2) {
    Name name1 = st_nativeInt64FromLittleEndian(region1->name);
    Name name2 = st_nativeInt64FromLittleEndian(region2->name);
    if (name1 < name2) {
        return -1;
    } else if (name1 == name2) {
        int64_t start1 = st_nativeInt64FromLittleEndian(region1->start);
        int64_t start2 = st_nativeInt64FromLittleEndian(region2->start);
        int64_t stop1 = st_nativeInt64FromLittleEndian(region1->stop);
        int64_t stop2 = st_nativeInt64FromLittleEndian(region2->stop);
        if (stop1 <= start2) {
            return -1;
        }
        if (stop2 <= start1) {
            return 1;
        }
        return 0;
    } else {
        assert(name1 > name2);
        return 1;
    }
}

bedRegion *bedRegion_construct(Name name, int64_t start, int64_t stop) {
    bedRegion *ret = st_malloc(sizeof(bedRegion));
    ret->name = st_nativeInt64ToLittleEndian(name);
    ret->start = st_nativeInt64ToLittleEndian(start);
    ret->stop = st_nativeInt64ToLittleEndian(stop);
    return ret;
}

static Name bedRegion_name(const bedRegion *region) {
    return st_nativeInt64FromLittleEndian(region->name);
}

static int64_t bedRegion_start(const bedRegion *region) {
    return st_nativeInt64FromLittleEndian(region->start);
}

static int64_t bedRegion_stop(const bedRegion *region) {
    return st_nativeInt64FromLittleEndian(region->stop);
}

// Find the closest bed region to the segment that is still less than
// the segment. This sets us up nicely to iterate forward along the
// bed array.
static bedRegion *seekToProperBedRegion(bedRegion *beds, size_t numBeds,
                                        stPinchSegment *segment,
                                        Name name) {
    bedRegion *targetRegion = bedRegion_construct(name,
                                                  stPinchSegment_getStart(segment),
                                                  stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment));

    size_t start = 0;
    size_t stop = numBeds;
    size_t pivot;
    bedRegion *pivotRegion;
    while ((pivot = start + (stop - start) / 2) != start) {
        pivotRegion = beds + pivot;
        int cmp = bedRegion_cmp(pivotRegion, targetRegion);
        if (cmp == -1) {
            // pivot less than target
            start = pivot;
        } else {
            // pivot greater than or equal to target
            stop = pivot;
        }
    }
    pivotRegion = beds + pivot;
    free(targetRegion);
    return pivotRegion;
}

// Searches forward linearly for the next bedRegion instead of doing a
// binary search. Has better cache behavior than seekToProperBedRegion
// once we are in the right area of the bed mmap.
static bedRegion *fastForwardToProperBedRegion(bedRegion *beds,
                                               size_t numBeds,
                                               bedRegion *begin,
                                               stPinchSegment *segment,
                                               Name name) {
    bedRegion *targetRegion = bedRegion_construct(name,
                                                  stPinchSegment_getStart(segment),
                                                  stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment));
    bedRegion *curRegion = begin;
    bedRegion *end = beds + numBeds;
    while (curRegion != end && bedRegion_cmp(curRegion, targetRegion) < 0) {
        curRegion++;
    }
    free(targetRegion);

    if (curRegion == end) {
        return NULL;
    }
    return curRegion;
}

// Find any regions in this thread covered by outgroups that are in
// segments with no block, and "rescue" them into single-degree
// blocks.
void rescueCoveredRegions(stPinchThread *thread, bedRegion *beds, size_t numBeds,
                          Name name) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    bedRegion *curRegion = seekToProperBedRegion(beds, numBeds, segment, name);
    while (segment != NULL && curRegion->name <= name) {
        if (stPinchSegment_getBlock(segment) == NULL) {
            // This is a potentially rescuable segment.
            curRegion = fastForwardToProperBedRegion(beds, numBeds, curRegion,
                                                     segment, name);
            if (curRegion == NULL || bedRegion_name(curRegion) > name) {
                // Reached end of the correct section of bed file.
                break;
            }
            int64_t start = stPinchSegment_getStart(segment);
            assert(start >= stPinchThread_getStart(thread));
            int64_t end = start + stPinchSegment_getLength(segment);
            assert(end <= stPinchThread_getStart(thread) + stPinchThread_getLength(thread));
            bool inCoveredRegion = false;
            for (int64_t i = start; i < end; i++) {
                if (i >= bedRegion_start(curRegion) && i < bedRegion_stop(curRegion)) {
                    if (!inCoveredRegion && i != start) {
                        // Wasn't covered before but now we are. Have
                        // to split up this block.
                        stPinchSegment_split(segment, i - 1);
                        break;
                    }
                    inCoveredRegion = true;
                } else if (inCoveredRegion) {
                    // Was covered before but not anymore. Have to
                    // split up this block.
                    stPinchSegment_split(segment, i - 1);
                    break;
                }
            }
            if (inCoveredRegion) {
                // Rescue this segment.
                stPinchBlock_construct2(segment);
            }
        }
        segment = stPinchSegment_get3Prime(segment);
    }
}
