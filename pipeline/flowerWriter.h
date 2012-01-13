#ifndef FLOWER_WRITER_H_
#define FLOWER_WRITER_H_

/*
 * Functions for organising the communication of lists of flower as strings. Used to communicate
 * between getFlowers/extendFlowers and cactus_workflow.
 */

#include "sonLib.h"
#include "cactus.h"

typedef struct _flowerWriter {
        stList *flowers;
        int64_t maxFlowerGroupSize;
        FILE *fileHandle;
} FlowerWriter;

FlowerWriter *flowerWriter_construct(FILE *fileHandle, int64_t maxFlowerGroupSize);

void flowerWriter_destruct(FlowerWriter *flowerWriter);

/*
 * Adds a given flower to the list to output.
 */
void flowerWriter_add(FlowerWriter *flowerWriter, Flower *flower);


#endif FLOWER_WRITER_H_
