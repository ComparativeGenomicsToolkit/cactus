/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "adjacencySwitch.h"

struct _adjacencySwitch {
    AdjacencyPair *adjacencyPair1;
    AdjacencyPair *adjacencyPair2;
    bool switchState;
    uint32_t strength;
    int32_t pseudoAdjacencies;
};

static void adjacencySwitch_getAdjacencyPairs(AdjacencyPair *adjacencyPair1,
        AdjacencyPair *adjacencyPair2, AdjacencyPair **adjacencyPair3,
        AdjacencyPair **adjacencyPair4, bool switchState) {
    End *end1 = adjacencyPair_getEnd1(adjacencyPair1);
    End *end2 = adjacencyPair_getEnd2(adjacencyPair1);
    End *end3 = adjacencyPair_getEnd1(adjacencyPair2);
    End *end4 = adjacencyPair_getEnd2(adjacencyPair2);

    *adjacencyPair3 = adjacencyPair_construct(end1, switchState ? end3 : end4);
    *adjacencyPair4 = adjacencyPair_construct(end2, switchState ? end4 : end3);
}

AdjacencySwitch *adjacencySwitch_construct(AdjacencyPair *adjacencyPair1,
        AdjacencyPair *adjacencyPair2, bool switchState) {
#ifdef BEN_DEBUG
    assert(adjacencyPair_getGroup(adjacencyPair1) == adjacencyPair_getGroup(adjacencyPair2));
#endif
    AdjacencySwitch *adjacencySwitch = st_malloc(sizeof(AdjacencySwitch));
    adjacencySwitch->adjacencyPair1 = adjacencyPair1;
    adjacencySwitch->adjacencyPair2 = adjacencyPair2;
    adjacencySwitch->switchState = switchState;

    //Now calculate the strength and switch state..
    AdjacencyPair *adjacencyPair3, *adjacencyPair4;
    adjacencySwitch_getAdjacencyPairs(adjacencyPair1, adjacencyPair2,
            &adjacencyPair3, &adjacencyPair4, switchState);

    uint32_t i = adjacencyPair_getStrengthOfAdjacencyPair(adjacencyPair3);
    uint32_t j = adjacencyPair_getStrengthOfAdjacencyPair(adjacencyPair4);

    adjacencySwitch->strength = i + j; //The residualStrength
    adjacencySwitch->pseudoAdjacencies = (i == 0 ? 1 : 0) + (j == 0 ? 1 : 0);
    adjacencyPair_destruct(adjacencyPair3);
    adjacencyPair_destruct(adjacencyPair4);
    return adjacencySwitch;
}

void adjacencySwitch_destruct(AdjacencySwitch *adjacencySwitch) {
    free(adjacencySwitch);
}

AdjacencyPair *adjacencySwitch_getAdjacencyPair1(
        AdjacencySwitch *adjacencySwitch) {
    return adjacencySwitch->adjacencyPair1;
}

AdjacencyPair *adjacencySwitch_getAdjacencyPair2(
        AdjacencySwitch *adjacencySwitch) {
    return adjacencySwitch->adjacencyPair2;
}

int32_t adjacencySwitch_getNumberOfPseudoAdjacencies(AdjacencySwitch *adjacencySwitch) {
    return adjacencySwitch->pseudoAdjacencies;
}

uint32_t adjacencySwitch_getStrength(AdjacencySwitch *adjacencySwitch) {
    return adjacencySwitch->strength;
}

int adjacencySwitch_compareStrengthAndPseudoAdjacencies(AdjacencySwitch *adjacencySwitch1, AdjacencySwitch *adjacencySwitch2) {
    int32_t numberOfPsuedoAdjacencies1 = adjacencySwitch_getNumberOfPseudoAdjacencies(adjacencySwitch1);
    uint32_t adjacencySwitchStrength1 = adjacencySwitch_getStrength(adjacencySwitch1);
    int32_t numberOfPsuedoAdjacencies2 = adjacencySwitch_getNumberOfPseudoAdjacencies(adjacencySwitch2);
    uint32_t adjacencySwitchStrength2 = adjacencySwitch_getStrength(adjacencySwitch2);
    if(numberOfPsuedoAdjacencies1 < numberOfPsuedoAdjacencies2) {
        return 1;
    }
    if(numberOfPsuedoAdjacencies1 > numberOfPsuedoAdjacencies2) {
        return -1;
    }
    //equal number of pseudo adjacencies
    return (int32_t)adjacencySwitchStrength1 - adjacencySwitchStrength2;
}

void adjacencySwitch_switch(AdjacencySwitch *adjacencySwitch,
        stHash *adjacencies) {
    AdjacencyPair *adjacencyPair3, *adjacencyPair4;
    adjacencySwitch_getAdjacencyPairs(adjacencySwitch_getAdjacencyPair1(
            adjacencySwitch),
            adjacencySwitch_getAdjacencyPair2(adjacencySwitch),
            &adjacencyPair3, &adjacencyPair4, adjacencySwitch->switchState);

    adjacencyHash_remove(adjacencies, adjacencySwitch_getAdjacencyPair1(
            adjacencySwitch));
    adjacencyHash_remove(adjacencies, adjacencySwitch_getAdjacencyPair2(
            adjacencySwitch));
    adjacencyHash_add(adjacencies, adjacencyPair3);
    adjacencyHash_add(adjacencies, adjacencyPair4);
}

static stList *getAdjacencyPairs(stList *component, stHash *adjacencies) {
    stList *adjacencyPairs = stList_construct();
    for (int32_t i = 0; i < stList_length(component); i++) {
        End *end = stList_get(component, i);
        AdjacencyPair *adjacencyPair = stHash_search(adjacencies, end);
        assert(adjacencyPair != NULL);
        if(adjacencyPair_getEnd1(adjacencyPair) == end) {
            stList_append(adjacencyPairs, adjacencyPair);
        }
    }
    return adjacencyPairs;
}

AdjacencySwitch *adjacencySwitch_getStrongestAdjacencySwitch(stList *component,
        stList *component2, stHash *adjacencies) {
    AdjacencySwitch *adjacencySwitch = NULL;

    stList *adjacencyPairs1 = getAdjacencyPairs(component, adjacencies);
    stList *adjacencyPairs2 = getAdjacencyPairs(component2, adjacencies);

    for (int32_t i = 0; i < stList_length(adjacencyPairs1); i++) {
        AdjacencyPair *adjacencyPair = stList_get(adjacencyPairs1, i);
        Group *group = adjacencyPair_getGroup(adjacencyPair);
        for (int32_t j = 0; j < stList_length(adjacencyPairs2); j++) {
            AdjacencyPair *adjacencyPair2 = stList_get(adjacencyPairs2, j);
            Group *group2 = adjacencyPair_getGroup(adjacencyPair2);
            if(group == group2) {
                for (int32_t k = 0; k <= 1; k++) {
                    AdjacencySwitch *adjacencySwitch2 =
                            adjacencySwitch_construct(adjacencyPair,
                                    adjacencyPair2, k);
                    if(adjacencySwitch == NULL) {
                        adjacencySwitch = adjacencySwitch2;
                    }
                    else {
                        if(adjacencySwitch_compareStrengthAndPseudoAdjacencies(adjacencySwitch2, adjacencySwitch) > 0) {
                            adjacencySwitch_destruct(adjacencySwitch);
                            adjacencySwitch = adjacencySwitch2;
                        }
                        else {
                            adjacencySwitch_destruct(adjacencySwitch2);
                        }
                    }
                }
            }
        }
    }

    stList_destruct(adjacencyPairs1);
    stList_destruct(adjacencyPairs2);

    return adjacencySwitch;
}
