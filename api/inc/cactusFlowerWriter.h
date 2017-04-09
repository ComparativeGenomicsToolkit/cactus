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

typedef struct {
    stList *flowerNames;
    CactusDisk *cactusDisk;
    Flower *curFlower;
    size_t nextIdx;
} FlowerStream;

/*
 * Get a stream that will give you one flower at a time from
 * stdin. (It doesn't actually stream the flowers from the DB or the
 * names from stdin, but it's meant to be used to keep just one flower
 * in memory at a time.)
 *
 * Flower loading/unloading is managed for you, so keeping a reference
 * to one of the flowers around will cause problems.
 */
FlowerStream *flowerWriter_getFlowerStream(CactusDisk *cactusDisk, FILE *file);

/*
 * Free a flowerStream.
 */
void flowerStream_destruct(FlowerStream *flowerStream);

/*
 * Get the next flower, or NULL if there aren't any more in the stream.
 * NB: every call unloads the flower previously returned.
 */
Flower *flowerStream_getNext(FlowerStream *flowerStream);

/*
 * Get the total number of flowers in this flower stream. NB: not
 * affected by your current position within the stream.
 */
int64_t flowerStream_size(const FlowerStream *flowerStream);

#endif
