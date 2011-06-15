#include "cactus.h"
#include "sonLib.h"

static int32_t intersectionSize(stSortedSet *set, stList *list) {
    /*
     * Returns the intersection size of the list and the set.
     */
    int32_t count = 0;
    for (int32_t j = 0; j < stList_length(list); j++) {
        if (stSortedSet_search(set, stList_get(list, j)) != NULL) {
            count++;
        }
    }
    return count;
}

static stList *filter(stList *list, stSortedSet *set) {
    /*
     * Returns a new list, identical to list, but with any elements containined in set removed.
     */
    stList *list2 = stList_construct();
    for (int32_t i = 0; i < stList_length(list); i++) {
        void *o = stList_get(list, i);
        if (stSortedSet_search(set, o) != NULL) {
            stList_append(list2, o);
        }
    }
    return list2;
}

static stList *filterLists(stList *listOfLists, stSortedSet *set) {
    /*
     * Takes a list of lists and returns a new list of lists whose elements are the product of applying
     * filter to each member of listOfLists in the same order.
     */
    stList *listOfLists2 = stList_construct3(0, (void (*)(void *))stList_destruct);
    for (int32_t i = 0; i < stList_length(listOfLists); i++) {
        stList_append(listOfLists2, filter(stList_get(listOfLists, i), set));
    }
    return listOfLists2;
}

static stSortedSet *getStubAndBlockEdgeSortedSet(stList *stubEdges,
        stList *blockEdges) {
    /*
     * Returns a sorted set containing the edges from the list of stub and block edges.
     */
    stList *stubAndBlockEdges = stList_copy(stubEdges, NULL);
    stList_appendAll(stubAndBlockEdges, blockEdges);
    stSortedSet *stubAndBlockEdgesSet = stList_getSortedSet(stubAndBlockEdges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn); //Set of only stub and block edges.
    stList_destruct(stubAndBlockEdges);
    return stubAndBlockEdgesSet;
}

stList *getComponents(stList *edges) {
    /*
     * Gets a list of connected components.
     */
    return NULL;
}

static stList *getComponents2(stList *adjacencyEdges, stList *stubEdges,
        stList *blockEdges) {
    /*
     * Gets a list of connected components for a set of adjacency, stub and block edges.
     * If stubEdges or blockEdges are NULL then they are ignored.
     */
    stList *allEdges = stList_construct(); //Build a concatenated list of all the block, stub and adjacency edges.
    if (stubEdges != NULL) {
        stList_appendAll(allEdges, adjacencyEdges);
    }
    if (stubEdges != NULL) {
        stList_appendAll(allEdges, stubEdges);
    }
    if (blockEdges != NULL) {
        stList_appendAll(allEdges, blockEdges);
    }
    stList *components = getComponents(allEdges); //Gets the graph components.
    stList_destruct(allEdges); //Cleanup the all edges.
    return components;
}

stList *mergeComponents(stList *components, stList *weightedEdges) {
    /*
     * Returns a single component, as a list of edges, by doing length(components)-1
     * pairs of edge switches.
     *
     * consisting of the edges
     * in components, and a set of substituted edges from weight edges that
     * are used to form the single merged component.
     */
    return NULL;
}

stList *doBestMergeOfTwoComponents(stList *components, stList *weightEdges) {
    /*
     * Returns a new list of components in which
     */
}


static stList *getStubAndBlockEdgeFreeComponents(stList *components,
        stList *stubEdges, stList *blockEdges) {
    stList *adjacencyOnlyComponents = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges,
            blockEdges);
    stList *adjacencyOnlyComponents = filterLists(components,
            stubAndBlockEdgesSet);
    stSortedSet_destruct(stubAndBlockEdgesSet);

    return adjacencyOnlyComponents;
}

static joinLists(stList *listOfLists) {
    stList *joinedList = stList_construct();
    for (int32_t i = 0; i < stList_length(joinedList); i++) {
        stList_appendAll(joinedList, stList_get(listOfLists, i));
    }
    return joinedList;
}

static stList *mergeCycles(stList *chosenEdges, stList *adjacencyEdges,
        stList *stubEdges, stList *blockEdges) {
    /*
     * Calculate components.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, blockEdges);

    /*
     * Calculate components that contain no stub edges.
     */
    stSortedSet *stubEdgesSet = stList_getSortedSet(stubEdges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);

    stList *stubContainingComponents = stList_construct();
    stList *stubFreeComponents = stList_construct();

    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_length(components, i);
        stList_append(
                intersectionSize(stubContainingComponents, component) > 0 ? stubContainingComponents
                        : stubFreeComponents, component);
    }
    assert(stList_length(stubContainingComponents) > 0);

    stSortedSet(stubEdgesSet);
    stList_destruct(components);

    /*
     * Merge the stub containing components into one 'global' component
     */
    stList *globalComponent = stList_construct();
    for (int32_t i = 0; i < stList_length(stubContainingComponents); i++) {
        stList_appendAll(globalComponent,
                stList_get(stubContainingComponents, i));
    }

    stList_destruct(stubContainingComponents);

    /*
     * Remove the stub/block edges from the components.
     */
    stList_append(stubFreeComponents, globalComponent);
    stList *adjacencyOnlyComponents = getStubAndBlockEdgeFreeComponents(
            stubFreeComponents, stubEdges, blockEdges);

    stList_destruct(stubFreeComponents);

    /*
     * Merge stub free components into the others.
     */
    stList *mergedComponents = mergeComponents(adjacencyOnlyComponents,
            adjacencyEdges);

    stList_destruct(adjacencyOnlyComponents);

    /*
     * Create an updated list of chosen edges.
     */
    stList *updatedChosenEdges = joinLists(mergedComponents);

    stList_destruct(mergedComponents);

    return updatedChosenEdges;
}

static stList *splitMultipleStubComponent(stList *component,
        stList *adjacencyEdges, stList *stubEdges, stList *blockEdges,
        stSortedSet *stubAndBlockEdgesSet) {
    //Filter out the the non-adjacency edges
    stList *adjacencyOnlyComponent = filter(component, stubAndBlockEdgesSet);
    stList *stubFreeComponents = getComponents2(adjacencyOnlyComponent, NULL,
            blockEdges);
    assert(stList_length(stubFreeComponents) >= 1);

    stList *splitComponents = stList_construct();

    if (stList_length(stubFreeComponents) > 1) {
        stList *l = filterLists(stubFreeComponents, stubAndBlockEdgesSet);
        stList *l2 = doBestMergeOfTwoComponents(l, adjacencyEdges);
        stList_destruct(l);
        stList *l = joinLists(l2);
        stList_destruct(l2);
        stList *l = getComponents(l2, adjacencyEdges, stubEdges, blockEdges);
        stList_destruct(l2);
        for (int32_t i = 0; i < stList_length(l); i++) {
            l2 = splitMultipleStubComponent(stList_get(l, i), adjacencyEdges,
                    stubEdges, blockEdges, stubAndBlockEdgesSet);
            stList_appendAll(splitComponents, l2);
            stList_destruct(l2);
        }
        stList_destruct(l);
    } else {
        stList_append(splitComponents, component);
    }

    stList_destruct(adjacencyOnlyComponent);
    stList_destruct(stubFreeComponents);

    return splitComponents;
}

static stList *splitMultipleStubCycles(stList *chosenEdges,
        stList *adjacencyEdges, stList *stubEdges, stList *blockEdges) {
    /*
     * Calculate components.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, blockEdges);

    /*
     * Find components with multiple stub edges.
     */
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges,
            blockEdges);

    stList *singleStubEdgeComponents = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList *splitComponents = splitMultipleStubComponent(component,
                adjacencyEdges, blockEdges, stubEdges, stubAndBlockEdgesSet);
        stList_appendAll(singleStubEdgeComponents, splitComponents);

        stList_destruct(splitComponents);
    }

    stList_destruct(components);

    /*
     * Remove the stub/block edges from the components.
     */
    stList *adjacencyOnlyComponents = filterLists(singleStubEdgeComponents,
            stubAndBlockEdgesSet);

    stList_destruct(singleStubEdgeComponents);

    /*
     * Merge the adjacency edges in the components into a single list.
     */
    stList *updatedChosenEdges = joinLists(adjacencyOnlyComponents);

    stList_destruct(adjacencyOnlyComponents);

    return updatedChosenEdges;
}

stList *chooseMatching(uint32_t nodeNumber, stList *adjacencyEdges,
        stList *stubEdges, stList *blockEdges, bool makeStubCyclesDisjoint,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    /*
     * Check the inputs.
     */
    assert(nodeNumber % 2 == 0);
    assert(stList_length(stubEdges) > 0);
    assert(
            stList_length(stubEdges) + stList_length(blockEdges) == nodeNumber
                    / 2);
    checkEdges(stubEdges, nodeNumber);
    checkEdges(blockEdges, nodeNumber);
    checkEdges(adjacencyEdges, nodeNumber);

    /*
     * First calculate the optimal matching.
     */
    stList *chosenEdges = matchingAlgorithm(adjacencyEdges, nodeNumber);

    /*
     * Merge in the stub free components.
     */
    stList *updatedChosenEdges = mergeCycles(chosenEdges, adjacencyEdges,
            stubEdges, blockEdges);
    stList_destruct(chosenEdges);
    updatedChosenEdges = chosenEdges;

    /*
     * Split stub components.
     */
    if (makeStubCyclesDisjoint) {
        updatedChosenEdges = splitMultipleStubCycles(chosenEdges,
                adjacencyEdges, stubEdges, blockEdges);
        stList_destruct(chosenEdges);
        updatedChosenEdges = chosenEdges;
    }

    return chosenEdges;
}
