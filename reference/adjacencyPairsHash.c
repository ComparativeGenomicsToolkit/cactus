/*
 * adjacencyPairsHash.c
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"

void adjacencyHash_add(stHash *adjacenciesHash, AdjacencyPair *adjacencyPair) {
#ifdef BEN_DEBUG
    assert(stHash_search(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair)) == NULL);
    assert(stHash_search(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair)) == NULL);
#endif
    stHash_insert(adjacenciesHash, adjacencyPair_getEnd1(adjacencyPair),
            adjacencyPair);
    stHash_insert(adjacenciesHash, adjacencyPair_getEnd2(adjacencyPair),
            adjacencyPair);
}

void adjacencyHash_remove(stHash *adjacencies,
        AdjacencyPair *adjacencyPair) {
#ifdef BEN_DEBUG
    assert(stHash_search(adjacencies, adjacencyPair_getEnd1(adjacencyPair)) == adjacencyPair);
    assert(stHash_search(adjacencies, adjacencyPair_getEnd2(adjacencyPair)) == adjacencyPair);
#endif
    stHash_remove(adjacencies, adjacencyPair_getEnd1(adjacencyPair));
    stHash_remove(adjacencies, adjacencyPair_getEnd2(adjacencyPair));
}

void adjacencyHash_cleanUp(stHash *adjacencies, Flower *flower) {
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        AdjacencyPair *adjacencyPair = stHash_search(adjacencies, end);
        if (adjacencyPair != NULL) {
            adjacencyHash_remove(adjacencies, adjacencyPair);
            adjacencyPair_destruct(adjacencyPair);
        }
    }
    flower_destructEndIterator(endIterator);
    assert(stHash_size(adjacencies) == 0);
    stHash_destruct(adjacencies);
}
