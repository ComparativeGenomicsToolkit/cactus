// Only works on POSIX systems (*nix, Mac OS) due to opendir.
#define _POSIX_C_SOURCE 200809L

#include <dirent.h>
#include "cactus.h"
#include "sonLib.h"
#include "stPinchGraphs.h"

typedef struct {
    Name name;
    int64_t start;
    int64_t stop;
} bedRegion;

// reads the next bed line from the file. Assumes tab-delimited bed3 input.
bedRegion *readNextBedLine(FILE *bedFile) {
    Name name;
    int64_t start, stop;
    if (fscanf(bedFile, "%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n",
               &name, &start, &stop) != EOF) {
        assert(start >= 0);
        assert(stop >= start);
        bedRegion *ret = st_malloc(sizeof(bedRegion));
        ret->name = name;
        ret->start = start;
        ret->stop = stop;
        return ret;
    } else {
        return NULL;
    }
}

static void fastForwardToProperBedLine(FILE *bedFile, stPinchSegment *segment,
                                       bedRegion **curBedLine) {
    Name threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
    while (*curBedLine != NULL &&
           (((*curBedLine)->name < threadName) ||
            (((*curBedLine)->name == threadName) &&
             ((*curBedLine)->stop <= stPinchSegment_getStart(segment))))) {
        free(*curBedLine);
        *curBedLine = readNextBedLine(bedFile);
    }
}

// Find any regions in this thread covered by outgroups that are in
// segments with no block, and "rescue" them into single-degree
// blocks.
// The input bed file must be sorted by chromosome (numerically!), then
// start (numerically), then stop (numerically). No overlaps are allowed.
bedRegion *rescueCoveredRegions(stPinchThread *thread, FILE *bedFile,
                          bedRegion *curBedLine) {
    if (curBedLine != NULL && stPinchThread_getName(thread) < curBedLine->name) {
        return curBedLine;
    }
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    fastForwardToProperBedLine(bedFile, segment, &curBedLine);
    while (segment != NULL && curBedLine != NULL && curBedLine->name == stPinchThread_getName(thread)) {
        fastForwardToProperBedLine(bedFile, segment, &curBedLine);
        if (curBedLine == NULL || curBedLine->name != stPinchThread_getName(thread)) {
            // Reached end of the correct section of bed file.
            break;
        }
        if (stPinchSegment_getBlock(segment) == NULL) {
            // This is a potentially rescuable segment.
            int64_t start = stPinchSegment_getStart(segment);
            assert(start >= stPinchThread_getStart(thread));
            int64_t end = start + stPinchSegment_getLength(segment);
            assert(end <= stPinchThread_getStart(thread) + stPinchThread_getLength(thread));
            bool inCoveredRegion = false;
            for (int64_t i = start; i < end; i++) {
                if (i >= curBedLine->start && i < curBedLine->stop) {
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
    return curBedLine;
}
