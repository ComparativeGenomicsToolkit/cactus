#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "sonLib.h"
#include "cactus.h"
#include "flowerWriter.h"

FlowerWriter *flowerWriter_construct(FILE *fileHandle, int64_t maxFlowerGroupSize,
        int64_t maxFlowerSecondaryGroupSize) {
    FlowerWriter *flowerWriter = st_malloc(sizeof(FlowerWriter));
    flowerWriter->flowerNamesAndSizes = stList_construct3(0, free);
    flowerWriter->fileHandle = fileHandle;
    flowerWriter->maxFlowerGroupSize = maxFlowerGroupSize;
    flowerWriter->maxFlowerSecondaryGroupSize = maxFlowerSecondaryGroupSize;
    return flowerWriter;
}

typedef struct _flowerNameAndSize {
    Name flowerName;
    int64_t flowerSize;
} FlowerNameAndSize;

/*
 * Adds a given flower to the list to output.
 */
void flowerWriter_add(FlowerWriter *flowerWriter, Name flowerName, int64_t flowerSize) {
    FlowerNameAndSize *flowerNameAndSize = st_malloc(sizeof(FlowerNameAndSize));
    flowerNameAndSize->flowerName = flowerName;
    flowerNameAndSize->flowerSize = flowerSize;
    stList_append(flowerWriter->flowerNamesAndSizes, flowerNameAndSize);
}

/*
 * Dumps the flower names out to the stream,
 */

static int compareFlowersByName(const void *a, const void *b) {
    return cactusMisc_nameCompare(((FlowerNameAndSize *)a)->flowerName, ((FlowerNameAndSize *)b)->flowerName);
}

static void printName(FILE *fileHandle, Name name) {
    fprintf(fileHandle, "%" PRIi64 " ", name);
}

void encodeFlowersString(stList *flowerNamesAndSizes, FILE *fileHandle,
        int64_t maxFlowerSecondaryGroupSize) {
    fprintf(fileHandle, "%i ", stList_length(flowerNamesAndSizes));
    if (stList_length(flowerNamesAndSizes) > 0) {
        FlowerNameAndSize *flowerNameAndSize = stList_get(flowerNamesAndSizes, 0);
        Name name = flowerNameAndSize->flowerName;
        int64_t totalSize = flowerNameAndSize->flowerSize;
        if(flowerNameAndSize->flowerSize > maxFlowerSecondaryGroupSize) {
            fprintf(fileHandle, "b ");
        }
        printName(fileHandle, name);
        for (int32_t i = 1; i < stList_length(flowerNamesAndSizes); i++) {
            flowerNameAndSize = stList_get(flowerNamesAndSizes, i);
            if(totalSize + flowerNameAndSize->flowerSize > maxFlowerSecondaryGroupSize) {
                totalSize = 0;
                if(flowerNameAndSize->flowerSize > maxFlowerSecondaryGroupSize) {
                    fprintf(fileHandle, "b ");
                }
                else {
                    fprintf(fileHandle, "a ");
                }
            }
            printName(fileHandle, flowerNameAndSize->flowerName - name);
            name = flowerNameAndSize->flowerName;
            totalSize += flowerNameAndSize->flowerSize;
        }
    }
}

static void printFlowers(FlowerWriter *flowerWriter, stList *stack, int64_t totalSize) {
    fprintf(stdout, "%i", totalSize > flowerWriter->maxFlowerGroupSize);
    encodeFlowersString(stack, stdout, flowerWriter->maxFlowerSecondaryGroupSize);
    fprintf(stdout, "\n");
}

static void flowerWriter_flush(FlowerWriter *flowerWriter) {
    if (stList_length(flowerWriter->flowerNamesAndSizes) == 0) {
        return;
    }
    stList_sort(flowerWriter->flowerNamesAndSizes, compareFlowersByName);
    FlowerNameAndSize *flowerNameAndSize = stList_get(flowerWriter->flowerNamesAndSizes, 0);
    int64_t totalSize = flowerNameAndSize->flowerSize;
    stList *stack = stList_construct();
    stList_append(stack, flowerNameAndSize);
    for(int32_t i=1; i<stList_length(flowerWriter->flowerNamesAndSizes); i++) {
        flowerNameAndSize = stList_get(flowerWriter->flowerNamesAndSizes, i);
        if (flowerNameAndSize->flowerSize + totalSize > flowerWriter->maxFlowerGroupSize) { //  || stList_length(stack) > 5000) {
            printFlowers(flowerWriter, stack, totalSize);
            while(stList_length(stack) > 0) {
                stList_pop(stack);
            }
            totalSize = flowerNameAndSize->flowerSize;
            stList_append(stack, flowerNameAndSize);
        }
        else {
            totalSize += flowerNameAndSize->flowerSize;
            stList_append(stack, flowerNameAndSize);
        }
    }
    if(stList_length(stack) > 0) {
        printFlowers(flowerWriter, stack, totalSize);
    }
    stList_destruct(stack);
}

void flowerWriter_destruct(FlowerWriter *flowerWriter) {
    flowerWriter_flush(flowerWriter);
    stList_destruct(flowerWriter->flowerNamesAndSizes);
    free(flowerWriter);
}

static bool isSeperator(char *cA) {
    return cA[0] == 'a' || cA[0] == 'b';
}

static Name getName(FILE *fileHandle) {
    while(1) {
        char cA[100];
        int32_t j = fscanf(fileHandle, "%s", cA);
        (void) j;
        assert(j == 0);
        if(isSeperator(cA)) {
            continue;
        }
        Name name;
        j = sscanf(cA, NAME_STRING, &name);
        assert(j == 1);
        return name;
    }
}

static void addName(stList *flowerNamesList, Name name) {
    int64_t *iA = st_malloc(sizeof(int64_t));
    iA[0] = name;
    stList_append(flowerNamesList, iA);
}

stList *flowerWriter_decodeFlowersString(CactusDisk *cactusDisk, FILE *fileHandle) {
    int32_t flowerArgumentNumber;
    int32_t j = fscanf(fileHandle, "%i", &flowerArgumentNumber);
    (void) j;
    assert(j == 1);
    stList *flowerNamesList = stList_construct3(0, free);
    Name name = getName(fileHandle);
    addName(flowerNamesList, name);
    for (int32_t i = 1; i < flowerArgumentNumber; i++) {
        Name nName = getName(fileHandle) + name;
        addName(flowerNamesList, nName);
        name = nName;
    }
    stList *flowers = cactusDisk_getFlowers(cactusDisk, flowerNamesList);
    stList_destruct(flowerNamesList);
    return flowers;
}

stList *flowerWriter_parseFlowersFromStdin(CactusDisk *cactusDisk) {
    return flowerWriter_decodeFlowersString(cactusDisk, stdin);
}
