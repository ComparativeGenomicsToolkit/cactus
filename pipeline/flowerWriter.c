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
    flowerWriter->flowers = stList_construct();
    flowerWriter->fileHandle = fileHandle;
    flowerWriter->maxFlowerGroupSize = maxFlowerGroupSize;
    return flowerWriter;
}

/*
 * Adds a given flower to the list to output.
 */
void flowerWriter_add(FlowerWriter *flowerWriter, Flower *flower) {
    stList_append(flowerWriter->flowers, flower);
}

/*
 * Dumps the flower names out to the stream,
 */

static int compareFlowersByName(const void *a, const void *b) {
    return cactusMisc_nameCompare(flower_getName((Flower *)a), flower_getName((Flower *)b));
}

static void writeFlowers(FILE *fileHandle, Name firstFlowerName, Name lastFlowerName, int64_t totalSize) {
    fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %" PRIi64 "\n",
            firstFlowerName, lastFlowerName - firstFlowerName + 1, totalSize);
}

static void flowerWriter_flush(FlowerWriter *flowerWriter) {
    if(stList_length(flowerWriter->flowers) == 0) {
        return;
    }
    stList_sort(flowerWriter->flowers, compareFlowersByName);
    Flower *flower = stList_get(flowerWriter->flowers, 0);
    Name firstFlowerName = flower_getName(flower);
    Name lastFlowerName = firstFlowerName;
    int64_t totalSize = flower_getTotalBaseLength(flower);
    for(int32_t i=1; i<stList_length(flowerWriter->flowers); i++) {
        flower = stList_get(flowerWriter->flowers, i);
        int64_t flowerSize = flower_getTotalBaseLength(flower);
        assert(flower_getName(flower) - lastFlowerName >= 1);
        if(flower_getName(flower) - lastFlowerName > 1 || flowerSize + totalSize > flowerWriter->maxFlowerGroupSize) {
            writeFlowers(flowerWriter->fileHandle, firstFlowerName, lastFlowerName, totalSize);
            firstFlowerName = flower_getName(flower);
            lastFlowerName = firstFlowerName;
            totalSize = flowerSize;
        }
        else {
            totalSize += flowerSize;
            lastFlowerName = flower_getName(flower);
        }
    }
    writeFlowers(flowerWriter->fileHandle, firstFlowerName, lastFlowerName, totalSize);
}

void flowerWriter_destruct(FlowerWriter *flowerWriter) {
    flowerWriter_flush(flowerWriter);
    stList_destruct(flowerWriter->flowers);
    free(flowerWriter);
}

