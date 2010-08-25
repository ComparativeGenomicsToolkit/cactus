#include "cactus.h"
#include "sonLib.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "adjacencySwitch.h"

struct _adjacencySwitch {
    AdjacencyPair *adjacencyPair1;
    AdjacencyPair *adjacencyPair2;
    bool switchState;
    double strength;
};

AdjacencySwitch *adjacencySwitch_construct(AdjacencyPair *adjacencyPair1,
        AdjacencyPair *adjacencyPair2) {
    assert(adjacencyPair_getGroup(adjacencyPair1) == adjacencyPair_getGroup(adjacencyPair2));
    AdjacencySwitch *adjacencySwitch = st_malloc(sizeof(AdjacencySwitch));
    adjacencySwitch->adjacencyPair1 = adjacencyPair1;
    adjacencySwitch->adjacencyPair2 = adjacencyPair2;
    //Now calculate the strength and switch state..
    double initialStrength = adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair1) + adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair2);
    AdjacencyPair *adjacencyPair3 = adjacencyPair_construct(
            adjacencyPair_getEnd1(adjacencyPair1), adjacencyPair_getEnd1(
                    adjacencyPair2));
    AdjacencyPair *adjacencyPair4 = adjacencyPair_construct(
            adjacencyPair_getEnd2(adjacencyPair1), adjacencyPair_getEnd2(
                    adjacencyPair2));
    double residualStrength1 = adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair3) + adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair4);
    adjacencyPair_destruct(adjacencyPair3);
    adjacencyPair_destruct(adjacencyPair4);
    adjacencyPair3 = adjacencyPair_construct(adjacencyPair_getEnd1(
            adjacencyPair1), adjacencyPair_getEnd2(adjacencyPair2));
    adjacencyPair4 = adjacencyPair_construct(adjacencyPair_getEnd2(
            adjacencyPair1), adjacencyPair_getEnd1(adjacencyPair2));
    double residualStrength2 = adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair3) + adjacencyPair_getStrengthOfAdjacencyPair(
            adjacencyPair4);
    adjacencyPair_destruct(adjacencyPair3);
    adjacencyPair_destruct(adjacencyPair4);
    if (residualStrength1 >= residualStrength2) {
        adjacencySwitch->strength = residualStrength1 - initialStrength;
        adjacencySwitch->switchState = 0;
    } else {
        adjacencySwitch->strength = residualStrength2 - initialStrength;
        adjacencySwitch->switchState = 1;
    }

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

double adjacencySwitch_getStrength(AdjacencySwitch *adjacencySwitch) {
    return adjacencySwitch->strength;
}

int32_t adjacencySwitch_cmpByStrength(AdjacencySwitch *adjacencySwitch,
        AdjacencySwitch *adjacencySwitch2) {
    return adjacencySwitch_getStrength(adjacencySwitch)
            > adjacencySwitch_getStrength(adjacencySwitch2) ? 1
            : (adjacencySwitch_getStrength(adjacencySwitch)
                    < adjacencySwitch_getStrength(adjacencySwitch2) ? -1 : 0);
}

void adjacencySwitch_switch(AdjacencySwitch *adjacencySwitch,
        stHash *adjacencies) {
    End *end1 = adjacencyPair_getEnd1(adjacencySwitch_getAdjacencyPair1(
            adjacencySwitch));
    End *end2 = adjacencyPair_getEnd2(adjacencySwitch_getAdjacencyPair1(
            adjacencySwitch));
    End *end3 = adjacencyPair_getEnd1(adjacencySwitch_getAdjacencyPair2(
            adjacencySwitch));
    End *end4 = adjacencyPair_getEnd2(adjacencySwitch_getAdjacencyPair2(
            adjacencySwitch));

    AdjacencyPair *adjacencyPair3 = adjacencyPair_construct(end1,
            adjacencySwitch->switchState ? end3 : end4);
    AdjacencyPair *adjacencyPair4 = adjacencyPair_construct(end2,
            adjacencySwitch->switchState ? end4 : end3);

    adjacencyHash_remove(adjacencies, adjacencySwitch_getAdjacencyPair1(
            adjacencySwitch));
    adjacencyHash_remove(adjacencies, adjacencySwitch_getAdjacencyPair2(
            adjacencySwitch));
    adjacencyHash_add(adjacencies, adjacencyPair3);
    adjacencyHash_add(adjacencies, adjacencyPair4);
}

AdjacencySwitch *adjacencySwitch_getStrongestAdjacencySwitch(stList *component,
        stList *component2) {
    AdjacencySwitch *adjacencySwitch = NULL;
    for (int32_t i = 0; i < stList_length(component); i++) {
        AdjacencyPair *adjacencyPair = stList_get(component, i);
        Group *group = adjacencyPair_getGroup(adjacencyPair);
        if (group_isTangle(group)) {
            for (int32_t j = 0; j < stList_length(component2); j++) {
                AdjacencyPair *adjacencyPair2 = stList_get(component2, j);
                Group *group2 = adjacencyPair_getGroup(adjacencyPair2);
                if (group_isTangle(group2)) {
                    AdjacencySwitch *adjacencySwitch2 =
                            adjacencySwitch_construct(adjacencyPair,
                                    adjacencyPair2);
                    if (adjacencySwitch2 != NULL) {
                        if (adjacencySwitch == NULL
                                || adjacencySwitch_cmpByStrength(
                                        adjacencySwitch2, adjacencySwitch) >= 0) {
                            adjacencySwitch = adjacencySwitch2;
                        }
                    }
                }
            }
        }
    }
    return adjacencySwitch;
}
