/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"

struct _adjacencyPair {
    /*
     * Struct to represent a pair of 'adjacent' ends.
     */
    End *end1;
    End *end2;
};

End *adjacencyPair_getEnd1(AdjacencyPair *adjacencyPair) {
    return adjacencyPair->end1;
}

End *adjacencyPair_getEnd2(AdjacencyPair *adjacencyPair) {
    return adjacencyPair->end2;
}

Group *adjacencyPair_getGroup(AdjacencyPair *adjacencyPair) {
    return end_getGroup(adjacencyPair_getEnd1(adjacencyPair));
}

AdjacencyPair *adjacencyPair_construct(End *end1, End *end2) {
#ifdef BEN_DEBUG
    assert(end_getOrientation(end1));
    assert(end_getOrientation(end2));
    assert(end1 != end2);
    assert(end_getName(end1) != end_getName(end2));
    assert(end_getFlower(end1) == end_getFlower(end2)); //same flower
    assert(end_getGroup(end1) == end_getGroup(end2));
#endif
    AdjacencyPair *adjacencyPair = st_malloc(sizeof(AdjacencyPair));
    if (cactusMisc_nameCompare(end_getName(end1), end_getName(end2)) <= 0) {
        adjacencyPair->end2 = end1;
        adjacencyPair->end1 = end2;
    } else {
        adjacencyPair->end1 = end1;
        adjacencyPair->end2 = end2;
    }
    assert(cactusMisc_nameCompare(end_getName(adjacencyPair_getEnd1(adjacencyPair)), end_getName(adjacencyPair_getEnd2(adjacencyPair))) >= 0);
    return adjacencyPair;
}

End *adjacencyPair_getOtherEnd(AdjacencyPair *adjacencyPair, End *end) {
#ifdef BEN_DEBUG
    assert(end_getOrientation(end));
    assert(adjacencyPair_getEnd1(adjacencyPair) == end || adjacencyPair_getEnd2(adjacencyPair) == end);
#endif
    return adjacencyPair_getEnd1(adjacencyPair) == end ? adjacencyPair_getEnd2(
            adjacencyPair) : adjacencyPair_getEnd1(adjacencyPair);
}

void adjacencyPair_destruct(AdjacencyPair *adjacencyPair) {
    free(adjacencyPair);
}

static int32_t childInstanceNumber(End *end) {
    End_InstanceIterator *iterator = end_getInstanceIterator(end);
    Cap *cap;
    int32_t i = 0;
    while ((cap = end_getNext(iterator)) != NULL) {
        if(cap_getChildNumber(cap) == 0) {
            i++;
        }
    }
    end_destructInstanceIterator(iterator);
    return i;
}

uint32_t adjacencyPair_getStrengthOfAdjacencyPair(AdjacencyPair *adjacencyPair) {
    End *end1 = adjacencyPair_getEnd1(adjacencyPair);
    End *end2 = adjacencyPair_getEnd2(adjacencyPair);
#ifdef BEN_DEBUG
    assert(end_getOrientation(end1));
    assert(end_getOrientation(end2));
#endif
    int32_t i = childInstanceNumber(end1) + childInstanceNumber(end2);

    if (i == 0) { //avoid divide by zero.
        return 0;
    }
    int32_t j = 0;
    //Walk through all the adjacencies and increment j by 2 if
    //the adjacency is between the two ends.
    End_InstanceIterator *iterator = end_getInstanceIterator(end1);
    Cap *cap;
    while ((cap = end_getNext(iterator)) != NULL) {
        Cap *cap2 = cap_getAdjacency(cap);
        if (cap2 != NULL && end_getPositiveOrientation(cap_getEnd(cap2))
                == end2) {
            j += 2;
        }
    }
    end_destructInstanceIterator(iterator);
    double strength = ((double) j) / i;
#ifdef BEN_DEBUG
    assert(strength <= 1.01);
    assert(strength >= -0.01);
#endif
    return strength > 0 ? strength * 10000 : 0;
}

uint32_t adjacencyPair_hashKey(AdjacencyPair *adjacencyPair) {
    return end_getName(adjacencyPair_getEnd1(adjacencyPair));
}

int adjacencyPair_hashEqual(AdjacencyPair *adjacencyPair1,
        AdjacencyPair *adjacencyPair2) {
    return adjacencyPair_getEnd1(adjacencyPair1) == adjacencyPair_getEnd1(
            adjacencyPair2) && adjacencyPair_getEnd2(adjacencyPair1)
            == adjacencyPair_getEnd2(adjacencyPair2);
}
