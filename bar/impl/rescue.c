// "bar rescue" -- ensures that any ingroup sequence that found an
// outgroup alignment in the blast stage still makes it into the
// ancestor after the bar stage.

// Only works on POSIX systems (*nix, Mac OS) due to getline.
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

// reads the next bed line from the file. Assumes tab-delimited bed3
// or greater input.
// returns EOF on failure to read, returns 0 otherwise.
int readNextBedLine(FILE *bedFile, bedRegion *curBedLine) {
    char *line;
    size_t n = 0;
    if (getline(&line, &n, bedFile) == -1) {
        free(line);
        return EOF;
    }
    int ret = sscanf(line, "%" PRIi64 "\t%" PRIi64 "\t%" PRIi64,
                     &curBedLine->name, &curBedLine->start, &curBedLine->stop);
    if (ret == 3) {
        assert(curBedLine->start >= 0);
        assert(curBedLine->stop >= curBedLine->start);
        free(line);
        return 0;
    } else {
        free(line);
        st_errAbort("Improperly formatted bed file.");
    }
}

// Iterate through the bed file until we reach a line that's relevant
// to this segment.
static int fastForwardToProperBedLine(FILE *bedFile, stPinchSegment *segment,
                                      bedRegion *curBedLine) {
    Name threadName = stPinchThread_getName(stPinchSegment_getThread(segment));
    while (curBedLine != NULL &&
           ((curBedLine->name < threadName) ||
            ((curBedLine->name == threadName) &&
             (curBedLine->stop <= stPinchSegment_getStart(segment))))) {
        if (readNextBedLine(bedFile, curBedLine) == EOF) {
            return EOF;
        }
    }
    return 0;
}

// Find any regions in this thread covered by outgroups that are in
// segments with no block, and "rescue" them into single-degree
// blocks.
// The input bed file must be sorted by chromosome (numerically!), then
// start (numerically), then stop (numerically). No overlaps are allowed.
void rescueCoveredRegions(stPinchThread *thread, FILE *bedFile,
                          bedRegion *curBedLine) {
    if (stPinchThread_getName(thread) < curBedLine->name) {
        return;
    }
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    if (fastForwardToProperBedLine(bedFile, segment, curBedLine) == EOF) {
        // Reached end of bed file.
        return;
    }
    while (segment != NULL && curBedLine->name == stPinchThread_getName(thread)) {
        int bedStatus = fastForwardToProperBedLine(bedFile, segment, curBedLine);
        if (bedStatus == EOF || curBedLine->name != stPinchThread_getName(thread)) {
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
}
