/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "sonLib.h"

/*
 * Prints info about a flower.
 */

int main(int argc, char *argv[]) {
    st_setLogLevelFromString(argv[1]);
    st_logDebug("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    st_logDebug("Set up the flower disk\n");

    Name flowerName = cactusMisc_stringToName(argv[3]);
    Flower *flower = cactusDisk_getFlower(cactusDisk, flowerName);

    int64_t totalBases = flower_getTotalBaseLength(flower);
    int32_t totalEnds = flower_getEndNumber(flower);
    int32_t totalCaps = flower_getCapNumber(flower);
    int32_t maxEndDegree = 0;
    int64_t maxAdjacencyLength = 0;

    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_getInstanceNumber(end) > maxEndDegree) {
            maxEndDegree = end_getInstanceNumber(end);
        }
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        Cap *cap;
        while((cap = end_getNext(capIt)) != NULL) {
            if(cap_getSequence(cap) != NULL) {
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(adjacentCap != NULL);
                int64_t adjacencyLength = cap_getCoordinate(cap) - cap_getCoordinate(adjacentCap);
                if(adjacencyLength < 0) {
                    adjacencyLength *= -1;
                }
                assert(adjacencyLength >= 1);
                if(adjacencyLength >= maxAdjacencyLength) {
                    maxAdjacencyLength = adjacencyLength;
                }
            }
        }
    }

    printf("total bases: %lli total-ends: %i total-caps: %i max-end-degree: %i max-adjacency-length: %lli\n",
            totalBases, totalEnds, totalCaps, maxEndDegree, maxAdjacencyLength);

    return 0;
}
