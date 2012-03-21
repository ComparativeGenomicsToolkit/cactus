/*
 * recursiveFileBuilder.c
 *
 *  Created on: 16 Mar 2012
 *      Author: benedictpaten
 */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"
#include "sonLib.h"
#include "recursiveFileBuilder.h"

typedef struct _recursiveFileBuilderEntry {
    Name capName;
    FILE *fileHandle;
    int64_t fileStartOffset;
    int64_t entryLength;
} RecursiveFileBuilderEntry;

static int recursiveFileBuilderEntry_cmpFn(RecursiveFileBuilderEntry *r1,
        RecursiveFileBuilderEntry *r2) {
    return cactusMisc_nameCompare(r1->capName, r2->capName);
}

struct _recursiveFileBuilder {
    stSortedSet *recursiveFileBuilderEntries;
    stList *childFilePointers;
    FILE *parentFileHandle;
    bool hasParent;
    int64_t entryLengthPointer;
    int64_t entryStartPointer;
};

RecursiveFileBuilder *recursiveFileBuilder_construct(const char *childDir,
        FILE *parentFileHandle, bool hasParent) {
    RecursiveFileBuilder *recursiveFileBuilder = st_malloc(
            sizeof(RecursiveFileBuilder));

    recursiveFileBuilder->recursiveFileBuilderEntries
            = stSortedSet_construct3(
                    (int(*)(const void *, const void *)) recursiveFileBuilderEntry_cmpFn,
                    free);

    recursiveFileBuilder->childFilePointers = stList_construct3(0,
            (void(*)(void *)) fclose);
    recursiveFileBuilder->parentFileHandle = parentFileHandle;
    recursiveFileBuilder->hasParent = hasParent;

    stList *childFileNames = childDir == NULL ? stList_construct() : stFile_getFileNamesInDirectory(childDir);

    for (int32_t j=0; j<stList_length(childFileNames); j++) { //For each child adjacency file
        char *childFileName = stFile_pathJoin(childDir, stList_get(childFileNames, j));
        FILE *childFileHandle = fopen(childFileName, "r");
        stList_append(recursiveFileBuilder->childFilePointers, childFileHandle);
        while (1) {
            char *line = stFile_getLineFromFile(childFileHandle);
            if (line == NULL) {
                break;
            }
            //st_uglyf("I am going to process the following line: #%s# for index %i %i\n", line, j, stSortedSet_size(recursiveFileBuilder->recursiveFileBuilderEntries));
            RecursiveFileBuilderEntry *recursiveFileBuilderEntry = st_malloc(
                    sizeof(RecursiveFileBuilderEntry));
            recursiveFileBuilderEntry->fileHandle = childFileHandle;
            recursiveFileBuilderEntry->fileStartOffset = ftell(childFileHandle);
            assert(recursiveFileBuilderEntry->fileStartOffset >= 0);
            int i = sscanf(line, "%" PRIi64 " %" PRIi64 "", &(recursiveFileBuilderEntry->capName), &(recursiveFileBuilderEntry->entryLength));
            assert(i == 2);
            assert(recursiveFileBuilderEntry->entryLength >= 0);
            int64_t currentPosition = ftell(childFileHandle);
            assert(currentPosition != -1);
            i = fseek(childFileHandle, currentPosition + recursiveFileBuilderEntry->entryLength, SEEK_SET);
            assert(i == 0);
            assert(
                    stSortedSet_search(
                            recursiveFileBuilder->recursiveFileBuilderEntries,
                            recursiveFileBuilderEntry) == NULL);
            stSortedSet_insert(
                    recursiveFileBuilder->recursiveFileBuilderEntries,
                    recursiveFileBuilderEntry);
            assert(
                                stSortedSet_search(
                                        recursiveFileBuilder->recursiveFileBuilderEntries,
                                        recursiveFileBuilderEntry) != NULL);
            //Now walk ahead by the length of the entry.
        }
        free(childFileName);
    }
    stList_destruct(childFileNames);

    return recursiveFileBuilder;
}

void recursiveFileBuilder_destruct(RecursiveFileBuilder *recursiveFileBuilder) {
    stList_destruct(recursiveFileBuilder->childFilePointers);
    stSortedSet_destruct(recursiveFileBuilder->recursiveFileBuilderEntries);
    free(recursiveFileBuilder);
}

static void recursiveFileBuilder_writeHead(RecursiveFileBuilder *recursiveFileBuilder, Cap *cap) {
    if (recursiveFileBuilder->hasParent) {
        int64_t i = fprintf(recursiveFileBuilder->parentFileHandle, "%20" PRIi64 " ", cap_getName(cap));
        assert(i > 0);
        recursiveFileBuilder->entryLengthPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
        assert(recursiveFileBuilder->entryLengthPointer != -1);
        i = fprintf(recursiveFileBuilder->parentFileHandle, "%20" PRIi64 "\n", INT64_MAX);
        assert(i > 0);
        recursiveFileBuilder->entryStartPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
        assert(recursiveFileBuilder->entryStartPointer != -1);
    }
}

static void recursiveFileBuilder_writeTail(RecursiveFileBuilder *recursiveFileBuilder) {
    if (recursiveFileBuilder->hasParent) {
        int64_t currentPosition = ftell(recursiveFileBuilder->parentFileHandle);
        int i = fseek(recursiveFileBuilder->parentFileHandle,
                recursiveFileBuilder->entryLengthPointer, SEEK_SET);
        assert(i == 0);
        fprintf(recursiveFileBuilder->parentFileHandle, "%20" PRIi64 "",
                currentPosition - recursiveFileBuilder->entryStartPointer);
        i = fseek(recursiveFileBuilder->parentFileHandle, currentPosition,
                SEEK_SET);
        assert(i == 0);
    }
}

static void recursiveFileBuilder_writeAdjacency(
        RecursiveFileBuilder *recursiveFileBuilder, Cap *cap,
        void (*terminalAdjacencyWriteFn)(FILE *, Cap *)) {
    static RecursiveFileBuilderEntry staticEntry;
    staticEntry.capName = cap_getName(cap);
    RecursiveFileBuilderEntry *a = stSortedSet_search(
            recursiveFileBuilder->recursiveFileBuilderEntries, &staticEntry);
    if(a != NULL) {
        fseek(a->fileHandle, a->fileStartOffset, SEEK_SET);
        for (int64_t i = 0; i < a->entryLength; i++) {
            putc(getc(a->fileHandle), recursiveFileBuilder->parentFileHandle);
        }
    }
    else if(terminalAdjacencyWriteFn != NULL) {
        terminalAdjacencyWriteFn(recursiveFileBuilder->parentFileHandle, cap);
    }
}

static void recursiveFileBuilder_writeSegment(
        RecursiveFileBuilder *recursiveFileBuilder, Segment *segment,
        void (*segmentWriteFn)(FILE *, Segment *)) {
    segmentWriteFn(recursiveFileBuilder->parentFileHandle, segment);
}

void recursiveFileBuilder_writeThread(RecursiveFileBuilder *recursiveFileBuilder,
        Cap *cap, void (*segmentWriteFn)(FILE *, Segment *), void (*terminalAdjacencyWriteFn)(FILE *, Cap *)) {
    /*
     * Iterate along thread.
     */
    recursiveFileBuilder_writeHead(recursiveFileBuilder, cap);
    while (1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) >= 1);
        if (cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) > 1) {
            recursiveFileBuilder_writeAdjacency(recursiveFileBuilder, cap, terminalAdjacencyWriteFn);
        }
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        recursiveFileBuilder_writeSegment(recursiveFileBuilder, cap_getSegment(adjacentCap),
                segmentWriteFn);
    }
    recursiveFileBuilder_writeTail(recursiveFileBuilder);
}

char *recursiveFileBuilder_getUniqueFileName(Flower *flower, const char *directory) {
    return stString_print("%s/%s.hal", directory,
            cactusMisc_nameToStringStatic(flower_getName(flower)));
}
