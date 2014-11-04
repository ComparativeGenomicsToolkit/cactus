// Only works on POSIX systems (*nix, Mac OS) due to opendir.
#define _POSIX_C_SOURCE 200809L

#include <dirent.h>
#include "cactus.h"
#include "sonLib.h"
#include "stPinchGraphs.h"

// Returns a hash of sequence Name => coverage array.
// Not thread-safe.
stHash *getIngroupCoverage(const char *ingroupCoverageDir, Flower *flower) {
    stHash *ret = stHash_construct3((uint64_t (*)(const void *))stIntTuple_hashKey,
                                    (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct, free);
    // Look through the directory and read in all entries (we will
    // assume that the only thing in the ingroupCoverageDir is bed
    // files.)
    DIR *dir = opendir(ingroupCoverageDir);
    if (dir == NULL) {
        st_errnoAbort("Opening directory %s failed", ingroupCoverageDir);
    }
    struct dirent *dirent;
    while ((dirent = readdir(dir)) != NULL) {
        char *fullPath = stString_print("%s/%s", ingroupCoverageDir, dirent->d_name);
        FILE *bedFile = fopen(fullPath, "r");
        char *line = NULL;
        size_t n = 0;
        while (getline(&line, &n, bedFile) != -1) {
            stList *fields = stString_split(line);
            if (stList_length(fields) < 3) {
                st_errAbort("Invalid bed line (< 3 fields) in file %s: %s\n",
                            fullPath, line);
            }
            // Parse the info we need from the bed line.
            const char *capNameString = stList_get(fields, 0);
            Name capName;
            if (sscanf(capNameString, "%" PRIi64, &capName) != 1) {
                st_errAbort("Could not parse sequence Name %s in bed file %s",
                            capNameString, fullPath);
            }
            const char *startString = stList_get(fields, 1);
            int64_t start;
            if (sscanf(startString, "%" PRIi64, &start) != 1) {
                st_errAbort("Could not parse start pos %s in bed file %s",
                            startString, fullPath);
            }
            const char *endString = stList_get(fields, 2);
            int64_t end;
            if (sscanf(endString, "%" PRIi64, &end) != 1) {
                st_errAbort("Could not parse end pos %s in bed file %s",
                            endString, fullPath);
            }

            Cap *cap = flower_getCap(flower, capName);
            if (cap == NULL) {
                st_errAbort("In bed file %s: no cap loaded w/ Name %"
                            PRIi64, fullPath, capName);
            }

            Sequence *sequence = cap_getSequence(cap);
            assert(sequence != NULL);
            assert(end <= sequence_getLength(sequence) + 2);
            assert(start >= 0);
            Name sequenceName = sequence_getName(sequence);
            // OK, this is really, really inefficient. Should use char
            // and index into each bit. But this will all be replaced anyway.
            stIntTuple *queryTuple = stIntTuple_construct1(sequenceName);
            bool *coverageArray = stHash_search(ret, queryTuple);
            if (coverageArray == NULL) {
                // We haven't seen a bed entry with this sequence
                // yet. Initialize a coverage array.
                coverageArray = st_calloc(sequence_getLength(sequence),
                                          sizeof(bool));
                stHash_insert(ret, stIntTuple_construct1(sequenceName), coverageArray);
            }

            // Mark this region as covered.
            memset(coverageArray + start, 1, end - start);

            stIntTuple_destruct(queryTuple);
            stList_destruct(fields);
            free(line);
            line = NULL;
            n = 0;
        }
        fclose(bedFile);
    }
    closedir(dir);
    return ret;
}

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
