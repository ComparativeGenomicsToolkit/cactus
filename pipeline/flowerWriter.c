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

static void writeFlowers(FILE *fileHandle, Name firstFlowerName, Name lastFlowerName, int64_t totalSize) {
    fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %" PRIi64 "\n",
        firstFlowerName, lastFlowerName - firstFlowerName + 1, totalSize);
}

static void flowerWriter_flush(FlowerWriter *flowerWriter) {
    if (stList_length(flowerWriter->flowerNamesAndSizes) == 0) {
        return;
    }
    stList_sort(flowerWriter->flowerNamesAndSizes, compareFlowersByName);
    FlowerNameAndSize *flowerNameAndSize = stList_get(flowerWriter->flowerNamesAndSizes, 0);
    Name firstFlowerName = flowerNameAndSize->flowerName;
    Name lastFlowerName = firstFlowerName;
    int64_t totalSize = flowerNameAndSize->flowerSize;
    for (int32_t i = 1; i < stList_length(flowerWriter->flowerNamesAndSizes); i++) {
        flowerNameAndSize = stList_get(flowerWriter->flowerNamesAndSizes, i);
        assert(flowerNameAndSize->flowerName - lastFlowerName >= 1);
        if (flowerNameAndSize->flowerName - lastFlowerName > 1 ||
            flowerNameAndSize->flowerSize + totalSize > flowerWriter->maxFlowerGroupSize) {
            writeFlowers(flowerWriter->fileHandle, firstFlowerName, lastFlowerName, totalSize);
            firstFlowerName = flowerNameAndSize->flowerName;
            lastFlowerName = firstFlowerName;
            totalSize = flowerNameAndSize->flowerSize;
        } else {
            totalSize += flowerNameAndSize->flowerSize;
            lastFlowerName = flowerNameAndSize->flowerName;
        }
    }
    writeFlowers(flowerWriter->fileHandle, firstFlowerName, lastFlowerName, totalSize);
}

void flowerWriter_destruct(FlowerWriter *flowerWriter) {
    flowerWriter_flush(flowerWriter);
    stList_destruct(flowerWriter->flowerNamesAndSizes);
    free(flowerWriter);
}
