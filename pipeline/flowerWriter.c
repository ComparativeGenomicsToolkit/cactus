#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "sonLib.h"
#include "cactus.h"
#include "flowerWriter.h"

FlowerWriter *flowerWriter_construct(FILE *fileHandle, int64_t maxFlowerGroupSize) {
    FlowerWriter *flowerWriter = st_malloc(sizeof(FlowerWriter));
    flowerWriter->flowerNamesAndSizes = stList_construct3(0, free);
    flowerWriter->fileHandle = fileHandle;
    flowerWriter->maxFlowerGroupSize = maxFlowerGroupSize;
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

static void printFlowers(FlowerWriter *flowerWriter, stList *stack, int64_t totalSize) {
    fprintf(stdout, "%i", totalSize > flowerWriter->maxFlowerGroupSize);
    cactusMisc_encodeFlowersString(stack, stdout);
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
    stList_append(stack, &flowerNameAndSize->flowerName);
    for(int32_t i=1; i<stList_length(flowerWriter->flowerNamesAndSizes); i++) {
        flowerNameAndSize = stList_get(flowerWriter->flowerNamesAndSizes, i);
        if (flowerNameAndSize->flowerSize + totalSize > flowerWriter->maxFlowerGroupSize) { //  || stList_length(stack) > 5000) {
            printFlowers(flowerWriter, stack, totalSize);
            while(stList_length(stack) > 0) {
                stList_pop(stack);
            }
            totalSize = flowerNameAndSize->flowerSize;
            stList_append(stack, &flowerNameAndSize->flowerName);
        }
        else {
            totalSize += flowerNameAndSize->flowerSize;
            stList_append(stack, &flowerNameAndSize->flowerName);
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
