#ifndef FLOWER_WRITER_H_
#define FLOWER_WRITER_H_

/*
 * Functions for organising the communication of lists of flower as strings. Used to communicate
 * between getFlowers/extendFlowers and cactus_workflow.
 */

#include "sonLib.h"
#include "cactus.h"

FlowerWriter *flowerWriter_construct(FILE *fileHandle,
        int64_t maxFlowerGroupSize,
        int64_t maxFlowerSecondaryGroupSize);

void flowerWriter_destruct(FlowerWriter *flowerWriter);

/*
 * Adds a given flower to the list to output.
 */
void flowerWriter_add(FlowerWriter *flowerWriter, Name flowerName, int64_t flowerSize);

/*
 * Decodes a list of flower names and returns them from the filehandle.
 */
stList *flowerWriter_parseNames(FILE *fileHandle);

/*
 *  As parseFlowers, but reads the list from stdin.
 */
stList *flowerWriter_parseFlowersFromStdin(CactusDisk *cactusDisk);

#endif
