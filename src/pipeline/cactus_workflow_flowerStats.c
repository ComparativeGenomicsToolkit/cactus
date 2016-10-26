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
    int64_t totalEnds = flower_getEndNumber(flower);
    int64_t totalFreeEnds = flower_getFreeStubEndNumber(flower);
    int64_t totalAttachedEnds = flower_getAttachedStubEndNumber(flower);
    int64_t totalCaps = flower_getCapNumber(flower);
    int64_t totalBlocks = flower_getBlockNumber(flower);
    int64_t totalGroups = flower_getGroupNumber(flower);
    int64_t totalChains = flower_getChainNumber(flower);
    int64_t totalLinkGroups = 0;
    int64_t maxEndDegree = 0;
    int64_t maxAdjacencyLength = 0;
    int64_t totalEdges = 0;

    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while((end = flower_getNextEnd(endIt)) != NULL) {
        assert(end_getOrientation(end));
        if(end_getInstanceNumber(end) > maxEndDegree) {
            maxEndDegree = end_getInstanceNumber(end);
        }
        stSortedSet *ends = stSortedSet_construct();
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        Cap *cap;
        while((cap = end_getNext(capIt)) != NULL) {
            if(cap_getSequence(cap) != NULL) {
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(adjacentCap != NULL);
                End *adjacentEnd = end_getPositiveOrientation(cap_getEnd(adjacentCap));
                stSortedSet_insert(ends, adjacentEnd);
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
        end_destructInstanceIterator(capIt);
        totalEdges += stSortedSet_size(ends);
        if(stSortedSet_search(ends, end) != NULL) { //This ensures we count self edges twice, so that the division works.
            totalEdges += 1;
        }
        stSortedSet_destruct(ends);
    }
    assert(totalEdges % 2 == 0);
    flower_destructEndIterator(endIt);

    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(groupIt)) != NULL) {
        if(group_getLink(group) != NULL) {
            totalLinkGroups++;
        }
    }
    flower_destructGroupIterator(groupIt);

    printf("flower name: %" PRIi64 " total bases: %" PRIi64 " total-ends: %" PRIi64 " total-caps: %" PRIi64 " max-end-degree: %" PRIi64 " max-adjacency-length: %" PRIi64 " total-blocks: %" PRIi64 " total-groups: %" PRIi64 " total-edges: %" PRIi64 " total-free-ends: %" PRIi64 " total-attached-ends: %" PRIi64 " total-chains: %" PRIi64 " total-link groups: %" PRIi64 "\n",
            flower_getName(flower), totalBases, totalEnds, totalCaps, maxEndDegree, maxAdjacencyLength, totalBlocks, totalGroups, totalEdges/2, totalFreeEnds, totalAttachedEnds, totalChains, totalLinkGroups);

    return 0;
}
