/*
 * recursiveFileBuilder.c
 *
 *  Created on: 16 Mar 2012
 *      Author: benedictpaten
 */

#include <dirent.h>
#include <sys/stat.h>
#include <inttypes.h>

#include "cactus.h"
#include "sonLib.h"

typedef struct _recursiveFileBuilderEntry {
    Name capName;
    FILE *fileHandle;
    int64_t fileStartOffset;
    int64_t entryLength;
} RecursiveFileBuilderEntry;

int recursiveFileBuilderEntry_cmpFn(RecursiveFileBuilderEntry *r1,
        RecursiveFileBuilderEntry *r2) {
    return cactusMisc_nameCompare(r1->capName, r2->capName);
}

typedef struct _recursiveFileBuilder {
    stSortedSet *recursiveFileBuilderEntries;
    stList *childFilePointers;
    FILE *parentFileHandle;
    bool hasParent;
    int64_t entryLengthPointer;
    int64_t entryStartPointer;
} RecursiveFileBuilder;

#define ADJACENCY_INDEX_UNIQUE_STRING "ADJACENCY_INDEX"

stList *stFile_getFileNamesInDirectory(const char *dir) {
    stList *files = stList_construct3(0, free);
    DIR *dh = opendir(dir);
    struct dirent *file;//a 'directory entity' AKA file
    while ((file = readdir(dh)) != NULL) {
        if (file->d_name[0] != '.') {
            struct stat info;
            char *cA = pathJoin(dir, file->d_name);
            //ascertain if complete or not
            exitOnFailure(stat(cA, &info),
                    "Failed to get information about the file: %s\n",
                    file->d_name);
            if (!S_ISDIR(info.st_mode)) {
                st_logInfo("Processing file: %s\n", cA);
                stList_append(files, cA);
            } else {
                free(cA);
            }
        }
    }
    closedir(dh);
    return files;
}

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
        const char *childFileName = stList_get(childFileNames, j);
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
    }
    stList_destruct(childFileNames);

    return recursiveFileBuilder;
}

void recursiveFileBuilder_writeHead(RecursiveFileBuilder *recursiveFileBuilder, Cap *cap) {
    if (recursiveFileBuilder->hasParent) {
        fprintf(recursiveFileBuilder->parentFileHandle, "%" PRIi64 " ", cap_getName(cap));
        recursiveFileBuilder->entryLengthPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
        fprintf(recursiveFileBuilder->parentFileHandle, "%" PRIi64 "\n", INT64_MAX);
        recursiveFileBuilder->entryStartPointer = ftell(
                recursiveFileBuilder->parentFileHandle);
    }
}

void recursiveFileBuilder_writeTail(RecursiveFileBuilder *recursiveFileBuilder,
        Cap *cap) {
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

void recursiveFileBuilder_writeAdjacency(
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

void recursiveFileBuilder_writeSegment(
        RecursiveFileBuilder *recursiveFileBuilder, Segment *segment,
        void (*segmentWriteFn)(FILE *, Segment *)) {
    segmentWriteFn(recursiveFileBuilder->parentFileHandle, segment);
}
