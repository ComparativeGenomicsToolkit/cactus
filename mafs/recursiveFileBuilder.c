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
            sizeof(recursiveFileBuilder));

    recursiveFileBuilder->recursiveFileBuilderEntries
            = stSortedSet_construct3(
                    (int(*)(const void *, const void *)) recursiveFileBuilderEntry_cmpFn,
                    free);
    recursiveFileBuilder->childFilePointers = stList_construct3(0,
            (void(*)(void *)) fclose);
    recursiveFileBuilder->parentFileHandle = parentFileHandle;
    recursiveFileBuilder->hasParent = hasParent;

    stList *childFileNames = stFile_getFileNamesInDirectory(childDir);

    for (int32_t j; j<stList_length(childFileNames); j++) { //For each child adjacency file
        char *childFileName = stFile_pathJoin(childDir, stList_get(childFileNames, j));
        FILE *childFileHandle = fopen(childFileName, "r");
        stList_append(recursiveFileBuilder->childFilePointers, childFileHandle);
        while (1) {
            char *line = stFile_getLineFromFile(childFileHandle);
            if (line == NULL) {
                break;
            }
            RecursiveFileBuilderEntry *recursiveFileBuilderEntry = st_malloc(
                    sizeof(RecursiveFileBuilderEntry));
            recursiveFileBuilderEntry->fileHandle = childFileHandle;
            recursiveFileBuilderEntry->fileStartOffset = ftell(childFileHandle);
            int i = sscanf(line, "%" PRIi64 " %" PRIi64 "", &recursiveFileBuilderEntry->capName, &recursiveFileBuilderEntry->entryLength);
            assert(i == 2);
            assert(recursiveFileBuilderEntry->entryLength >= 0);
            assert(
                    stSortedSet_search(
                            recursiveFileBuilder->recursiveFileBuilderEntries,
                            recursiveFileBuilderEntry) == NULL);
            stSortedSet_insert(
                    recursiveFileBuilder->recursiveFileBuilderEntries,
                    recursiveFileBuilderEntry);
        }
        free(childFileName);
        fclose(childFileHandle);
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
        fprintf(recursiveFileBuilder->parentFileHandle, "%" PRIi64 " ", cap_getName(cap));
        recursiveFileBuilder->entryLengthPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
        fprintf(recursiveFileBuilder->parentFileHandle, "%" PRIi64 "\n", INT64_MAX);
        recursiveFileBuilder->entryStartPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
    }
}

static void recursiveFileBuilder_writeTail(RecursiveFileBuilder *recursiveFileBuilder) {
    if (recursiveFileBuilder->hasParent) {
        int64_t currentPosition = ftell(recursiveFileBuilder->parentFileHandle);
        int i = fseek(recursiveFileBuilder->parentFileHandle,
                recursiveFileBuilder->entryLengthPointer, SEEK_SET);
        assert(i);
        fprintf(recursiveFileBuilder->parentFileHandle, "%" PRIi64 "",
                currentPosition - recursiveFileBuilder->entryStartPointer);
        i = fseek(recursiveFileBuilder->parentFileHandle, currentPosition,
                SEEK_SET);
        assert(i);
    }
}

static void recursiveFileBuilder_writeAdjacency(
        RecursiveFileBuilder *recursiveFileBuilder, Cap *cap) {
    stInt64Tuple *capName = stInt64Tuple_construct(1, cap_getName(cap));
    RecursiveFileBuilderEntry *a = stSortedSet_search(
            recursiveFileBuilder->recursiveFileBuilderEntries, capName);
    assert(a != NULL);
    fseek(a->fileHandle, a->fileStartOffset, SEEK_SET);
    for (int64_t i = 0; i < a->entryLength; i++) {
        putc(getc(a->fileHandle), recursiveFileBuilder->parentFileHandle);
    }
}

static void recursiveFileBuilder_writeSegment(
        RecursiveFileBuilder *recursiveFileBuilder, Segment *segment,
        void (*segmentWriteFn)(FILE *, Segment *)) {
    segmentWriteFn(recursiveFileBuilder->parentFileHandle, segment);
}

void recursiveFileBuilder_writeThread(RecursiveFileBuilder *recursiveFileBuilder,
        Cap *cap, void (*segmentWriteFn)(FILE *, Segment *)) {
    /*
     * Iterate along thread.
     */
    recursiveFileBuilder_writeHead(recursiveFileBuilder, cap);
    while (1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) >= 1);
        if (cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) > 1) {
            recursiveFileBuilder_writeAdjacency(recursiveFileBuilder, cap);
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
    return stString_print("%s/%s.maf", directory,
            cactusMisc_nameToStringStatic(flower_getName(flower)));
}
