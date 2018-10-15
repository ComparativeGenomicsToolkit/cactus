#include "sonLib.h"
#include "cactusGlobalsPrivate.h"

#define FLOWER_STREAM_BATCH_SIZE 50

struct _flowerWriter {
        stList *flowerNamesAndSizes;
        int64_t maxFlowerGroupSize;
        int64_t maxFlowerSecondaryGroupSize;
        FILE *fileHandle;
};

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

static void printSize(FILE *fileHandle, int64_t size) {
    fprintf(fileHandle, "%" PRIi64 " ", size);
}

static void flowerWriter_writeFlowersString(stList *flowerNamesAndSizes, FILE *fileHandle,
        int64_t maxFlowerSecondaryGroupSize) {
    fprintf(fileHandle, "%" PRIi64 " ", stList_length(flowerNamesAndSizes));
    if (stList_length(flowerNamesAndSizes) > 0) {
        FlowerNameAndSize *flowerNameAndSize = stList_get(flowerNamesAndSizes, 0);
        Name name = flowerNameAndSize->flowerName;
        int64_t totalSize = flowerNameAndSize->flowerSize;
        if(flowerNameAndSize->flowerSize > maxFlowerSecondaryGroupSize) {
            fprintf(fileHandle, "b ");
        }
        printName(fileHandle, name);
        printSize(fileHandle, flowerNameAndSize->flowerSize);
        for (int64_t i = 1; i < stList_length(flowerNamesAndSizes); i++) {
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
            printSize(fileHandle, flowerNameAndSize->flowerSize);
            name = flowerNameAndSize->flowerName;
            totalSize += flowerNameAndSize->flowerSize;
        }
    }
}

static void printFlowers(FlowerWriter *flowerWriter, stList *stack, int64_t totalSize) {
    fprintf(flowerWriter->fileHandle, "%i ", totalSize > flowerWriter->maxFlowerGroupSize);
    flowerWriter_writeFlowersString(stack, flowerWriter->fileHandle, flowerWriter->maxFlowerSecondaryGroupSize);
    fprintf(flowerWriter->fileHandle, "\n");
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
    for(int64_t i=1; i<stList_length(flowerWriter->flowerNamesAndSizes); i++) {
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
        int64_t j = fscanf(fileHandle, "%s", cA);
        (void) j;
        assert(j == 1);
        if(isSeperator(cA)) {
            continue;
        }
        Name name;
        j = sscanf(cA, NAME_STRING, &name);
        assert(j == 1);
        return name;
    }
    return 0;
}

static void addName(stList *flowerNamesList, Name name) {
    int64_t *iA = st_malloc(sizeof(int64_t));
    iA[0] = name;
    stList_append(flowerNamesList, iA);
}

stList *flowerWriter_parseNames(FILE *fileHandle) {
    int64_t flowerArgumentNumber;
    int64_t j = fscanf(fileHandle, "%" PRIi64 "", &flowerArgumentNumber);
    (void) j;
    assert(j == 1);
    stList *flowerNamesList = stList_construct3(0, free);
    Name name = getName(fileHandle);
    addName(flowerNamesList, name);
    for (int64_t i = 1; i < flowerArgumentNumber; i++) {
        Name nName = getName(fileHandle) + name;
        addName(flowerNamesList, nName);
        name = nName;
    }
    return flowerNamesList;
}

stList *flowerWriter_parseFlowersFromStdin(CactusDisk *cactusDisk) {
    stList *flowerNamesList = flowerWriter_parseNames(stdin);
    stList *flowers = cactusDisk_getFlowers(cactusDisk, flowerNamesList);
    stList_destruct(flowerNamesList);
    return flowers;
}

static FlowerStream *flowerStream_construct(stList *flowerNames, CactusDisk *cactusDisk) {
    FlowerStream *ret = malloc(sizeof(FlowerStream));
    ret->flowerNames = flowerNames;
    ret->flowerBatch = stList_construct();
    ret->curFlower = NULL;
    ret->nextIdx = 0;
    ret->cactusDisk = cactusDisk;
    return ret;
}

FlowerStream *flowerWriter_getFlowerStream(CactusDisk *cactusDisk, FILE *file) {
    stList *flowerNamesList = flowerWriter_parseNames(file);
    return flowerStream_construct(flowerNamesList, cactusDisk);
}

void flowerStream_destruct(FlowerStream *flowerStream) {
    if (flowerStream->curFlower != NULL) {
        flower_destruct(flowerStream->curFlower, false);
    }
    stList_destruct(flowerStream->flowerBatch);
    stList_destruct(flowerStream->flowerNames);
    free(flowerStream);
}

Flower *flowerStream_getNext(FlowerStream *flowerStream) {
    if (flowerStream->curFlower != NULL) {
        // Unload the previously loaded flower.
        flower_destruct(flowerStream->curFlower, false);
    }
    if (flowerStream->nextIdx >= stList_length(flowerStream->flowerNames)) {
        flowerStream->curFlower = NULL;
        return NULL;
    }
    if (stList_length(flowerStream->flowerBatch) == 0) {
        // Time to load the next batch of flowers from the DB.
        // Get the next batch of names.
        int64_t batchStart = flowerStream->nextIdx;
        int64_t batchEnd = flowerStream->nextIdx + FLOWER_STREAM_BATCH_SIZE;
        if (batchEnd > stList_length(flowerStream->flowerNames)) {
            batchEnd = stList_length(flowerStream->flowerNames);
        }
        stList *namesBatch = stList_construct2(batchEnd - batchStart);
        for (int64_t i = batchStart; i < batchEnd; i++) {
            stList_set(namesBatch, i - batchStart, stList_get(flowerStream->flowerNames, i));
        }
        // We want to be able to treat the batch like a stack and get
        // the same order, so we reverse it.
        stList_reverse(namesBatch);
        stList_destruct(flowerStream->flowerBatch);
        flowerStream->flowerBatch = cactusDisk_getFlowers(flowerStream->cactusDisk, namesBatch);
        stList_destruct(namesBatch);
    }
    flowerStream->curFlower = stList_pop(flowerStream->flowerBatch);
    flowerStream->nextIdx++;
    return flowerStream->curFlower;
}

int64_t flowerStream_size(const FlowerStream *flowerStream) {
    return stList_length(flowerStream->flowerNames);
}
