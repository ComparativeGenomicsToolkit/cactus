/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus_workflow_shared.h"

/*
 * Iterates through a tree an extends the cactus by creating new flowers for
 * non-terminal groups. Used during the core stage.
 */

//#ifdef BEN_DEBUG
static bool isOneAdjacencyComponent(Group *group) {
    /*
     * Checks a group is one adjacency component.
     */
    /*
     * Get the ends in a set.
     */
    stSortedSet *ends = stSortedSet_construct();
    Group_EndIterator *endIt = group_getEndIterator(group);
    End *end;
    while((end = group_getNextEnd(endIt)) != NULL) {
        assert(end_getPositiveOrientation(end) == end);
        stSortedSet_insert(ends, end);
    }
    group_destructEndIterator(endIt);

    /*
     * Now iterate through connections.
     */
    assert(stSortedSet_size(ends) > 0);
    end = stSortedSet_getFirst(ends);
    stList *stack = stList_construct();
    stSortedSet_remove(ends, end);
    stList_append(stack, end);


    while(stList_length(stack) > 0) {
        end = stList_pop(stack);
        Cap *cap;
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        while((cap = end_getNext(capIt)) != NULL) {
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            End *adjacentEnd = end_getPositiveOrientation(cap_getEnd(adjacentCap));
            assert(adjacentEnd != NULL);
            if(stSortedSet_search(ends, adjacentEnd) != NULL) {
                stSortedSet_remove(ends, adjacentEnd);
                stList_append(stack, adjacentEnd);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    bool b = stSortedSet_size(ends) == 0;

    /*
     * Cleanup
     */
    stSortedSet_destruct(ends);
    stList_destruct(stack);
    return b;
}
//#endif


static void extendFlowers(Flower *flower, bool createRedundantFlowerLinks) {
    Flower_GroupIterator *groupIterator;
    Group *group;
    if (flower_builtBlocks(flower)) {
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                int64_t size = group_getTotalBaseLength(group);
                assert(size >= 0);
                if (size >= minFlowerSize && size <= maxFlowerSize) {
                    Flower *nestedFlower = group_makeNestedFlower(group);
                    if(createRedundantFlowerLinks) {
                        flower_setBuiltBlocks(nestedFlower, 1);
                        Group *nestedGroup = flower_getFirstGroup(nestedFlower);
                        nestedFlower = group_makeNestedFlower(nestedGroup);
                        flowerWriter_add(flowerWriter, flower_getName(nestedFlower), size);
                    }
                    else {
                        flowerWriter_add(flowerWriter, flower_getName(nestedFlower), size);
                    }
//#ifdef BEN_DEBUG
                    assert(isOneAdjacencyComponent(flower_getFirstGroup(nestedFlower)));
                    assert(!flower_builtBlocks(nestedFlower));
                    assert(flower_isLeaf(nestedFlower));
                    assert(flower_getBlockNumber(nestedFlower) == 0);
                    assert(flower_getGroupNumber(nestedFlower) == 1);
                    assert(flower_isTerminal(nestedFlower));
//#endif
                }
            }
        }
        flower_destructGroupIterator(groupIterator);
    } else { //something went wrong last time, and the flower hasn't been filled in.. so we'll return it
        //again.
        flowerWriter_add(flowerWriter, flower_getName(flower), flower_getTotalBaseLength(flower));
    }
}

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    parseArgs(argc, argv);

    st_logDebug("Starting extending flowers\n");
    stList *flowers = parseFlowersFromStdin(cactusDisk);
    //stList *flowers = parseFlowers(argv + 6, argc - 6, cactusDisk);
    for (int32_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        extendFlowers(flower, 1); //flower_getTotalBaseLength(flower) > INT64_MAX);
        assert(!flower_isParentLoaded(flower)); //The parent should not be loaded.
    }
    stList_destruct(flowers);
    st_logDebug("Finish extending flowers\n");

    flowerWriter_destruct(flowerWriter);

    cactusDisk_write(cactusDisk);
    st_logDebug("Updated the cactus disk\n");

    cactusDisk_destruct(cactusDisk);

    st_logDebug("Am finished\n");
    return 0;
}
