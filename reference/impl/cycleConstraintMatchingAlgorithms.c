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
    stList *listOfLists2 = stList_construct3(0, (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(listOfLists); i++) {
        stList_append(listOfLists2, filter(stList_get(listOfLists, i), set));
    }
    return listOfLists2;
}

static stSortedSet *getStubAndBlockEdgeSortedSet(stList *stubEdges, stList *blockEdges) {
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

static void getComponentsP(stHash *nodesToEdges, stIntTuple *edge, int32_t position) {
    stIntTuple *node = stIntTuple_construct(1, stIntTuple_getPosition(edge, position));
    stList *edges;
    if ((edges = stHash_search(nodesToEdges, node)) == NULL) {
        edges = stList_construct();
        stHash_insert(nodesToEdges, node, edges);
    } else {
        stIntTuple_destruct(node);
    }
    stList_append(edges, edge);
}

static void getComponentsP2(stHash *nodesToEdges, int32_t node, stSortedSet *component) {
    stIntTuple *key = stIntTuple_construct(1, node);
    stList *edges = stHash_search(nodesToEdges, key);
    if (edges != NULL) {
        stHash_remove(nodesToEdges, key);
        for (int32_t i = 0; i < stList_length(edges); i++) {
            stIntTuple *edge = stList_get(edges, i);
            if(stSortedSet_search(component, edge) == NULL) {
                stSortedSet_insert(component, edge);
            }
            getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 0), component);
            getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 1), component);
        }
        stList_destruct(edges);
    }
    stIntTuple_destruct(key);
}

stList *getComponents(stList *edges) {
    /*
     * Gets a list of connected components, each connected component
     * being represented as a list of the edges, such that each edge is in exactly one
     * connected component.
     */
    stHash *nodesToEdges = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn, NULL, NULL);

    /*
     * Build a hash of nodes to edges.
     */
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        getComponentsP(nodesToEdges, edge, 0);
        getComponentsP(nodesToEdges, edge, 1);
    }

    /*
     * Traverse the edges greedily
     */
    stList *components = stList_construct3(0, (void(*)(void *))stList_destruct);
    stList *keys = stHash_getKeys(nodesToEdges);
    for (int32_t i = 0; i < stList_length(keys); i++) {
        stIntTuple *key = stList_get(keys, i);
        stList *edges = stHash_search(nodesToEdges, key);
        if (edges != NULL) { //We have a component to build
            stSortedSet *component = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, NULL);
            stHash_remove(nodesToEdges, key);
            for (int32_t i = 0; i < stList_length(edges); i++) {
                stIntTuple *edge = stList_get(edges, i);
                getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 0), component);
                getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 1), component);
            }
            stList_append(components, stSortedSet_getList(component));
            stList_destruct(edges);
            stSortedSet_destruct(component);
        }
    }
    assert(stHash_size(nodesToEdges) == 0);
    stHash_destruct(nodesToEdges);
    stList_setDestructor(keys, (void (*)(void *))stIntTuple_destruct);
    stList_destruct(keys);

    return components;
}

static stList *getComponents2(stList *adjacencyEdges, stList *stubEdges, stList *blockEdges) {
    /*
     * Gets a list of connected components for a set of adjacency, stub and block edges.
     * If adjacencyEdges, stubEdges or blockEdges are NULL then they are ignored.
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

void doBestMergeOfTwoComponents(stList *components, stList *weightedEdges) {
    /*
     * Modifies in place the list of components, merging together two components whose switch has lowest reduction in score.
     */
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
    stList_copy

    return NULL;
}

static stList *getStubAndBlockEdgeFreeComponents(stList *components, stList *stubEdges, stList *blockEdges) {
    /*
     * Returns a new list of components, such that each component has its stub and block edges removed.
     */
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges, blockEdges);
    stList *adjacencyOnlyComponents = filterLists(components, stubAndBlockEdgesSet);
    stSortedSet_destruct(stubAndBlockEdgesSet);

    return adjacencyOnlyComponents;
}

static stList *joinLists(stList *listOfLists) {
    /*
     * Returns new list which contains elements of the list of list concatenated in one list.
     */
    stList *joinedList = stList_construct();
    for (int32_t i = 0; i < stList_length(joinedList); i++) {
        stList_appendAll(joinedList, stList_get(listOfLists, i));
    }
    return joinedList;
}

static stList *mergeComponents2(stList *chosenEdges, stList *adjacencyEdges, stList *stubEdges, stList *blockEdges) {
    /*
     * Returns a new set of chosen edges, modified by adjacency switches such that every component
     * contains at least one stub edge.
     */

    /*
     * Calculate components.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, blockEdges);

    /*
     * Divide the components by the presence of one or more stub edges.
     */
    stSortedSet *stubEdgesSet = stList_getSortedSet(stubEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);

    stList *stubContainingComponents = stList_construct();
    stList *stubFreeComponents = stList_construct();

    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList_append(
                intersectionSize(stubEdgesSet, component) > 0 ? stubContainingComponents
                        : stubFreeComponents, component);
    }
    assert(stList_length(stubContainingComponents) > 0);
    stSortedSet_destruct(stubEdgesSet);

    /*
     * Merge the stub containing components into one 'global' component
     */
    stList *globalComponent = joinLists(stubContainingComponents);
    stList_destruct(stubContainingComponents);

    /*
     * Remove the stub/block edges from the components.
     */
    stList_append(stubFreeComponents, globalComponent);
    stList *adjacencyOnlyComponents = getStubAndBlockEdgeFreeComponents(stubFreeComponents, stubEdges, blockEdges);

    stList_destruct(stubFreeComponents);
    stList_destruct(globalComponent);
    stList_destruct(components); //We only clean this up now, as this frees the components it contains.

    /*
     * Merge stub free components into the others.
     */
    stList *mergedComponents = mergeComponents(adjacencyOnlyComponents, adjacencyEdges);
    stList_destruct(adjacencyOnlyComponents);

    /*
     * Create an updated list of chosen edges.
     */
    stList *updatedChosenEdges = joinLists(mergedComponents);
    stList_destruct(mergedComponents);

    return updatedChosenEdges;
}

static stList *splitMultipleStubComponent(stList *component, stList *adjacencyEdges, stList *stubEdges,
        stList *blockEdges, stSortedSet *stubAndBlockEdgesSet) {
    /*
     * Get sub-components containing only adjacency and block edges.
     */
    //Filter out the the non-adjacency edges
    stList *adjacencyOnlyComponent = filter(component, stubAndBlockEdgesSet);
    //Make it only the block edges present in the original component
    stList *stubFreeComponents = getComponents2(adjacencyOnlyComponent, NULL, blockEdges);
    stList_destruct(adjacencyOnlyComponent);
    assert(stList_length(stubFreeComponents) >= 1);

    stList *splitComponents = stList_construct(); //The list to return.

    if (stList_length(stubFreeComponents) > 1) {
        /*
         * Merge together the best two components.
         */
        stList *l = filterLists(stubFreeComponents, stubAndBlockEdgesSet);
        doBestMergeOfTwoComponents(l, adjacencyEdges); //This is inplace.
        stList *l2 = joinLists(l);
        stList_destruct(l);
        l = getComponents2(l2, stubEdges, blockEdges);
        stList_destruct(l2);

        /*
         * Call procedure recursively.
         */
        for (int32_t i = 0; i < stList_length(l); i++) {
            l2 = splitMultipleStubComponent(stList_get(l, i), adjacencyEdges, stubEdges, blockEdges,
                    stubAndBlockEdgesSet);
            stList_appendAll(splitComponents, l2);
            stList_setDestructor(l2, NULL);
            stList_destruct(l2);
        }
        stList_destruct(l);
    } else {
        stList_append(splitComponents, stList_copy(component, NULL));
    }

    stList_destruct(stubFreeComponents);

    return splitComponents;
}

static stList *splitMultipleStubCycles(stList *chosenEdges, stList *adjacencyEdges, stList *stubEdges,
        stList *blockEdges) {
    /*
     *  Returns an updated list of adjacency edges, such that each stub edge is a member of exactly one cycle.
     */

    /*
     * Calculate components.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, blockEdges);

    /*
     * Find components with multiple stub edges.
     */
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges, blockEdges);
    stList *singleStubEdgeComponents = stList_construct3(0, (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList *splitComponents = splitMultipleStubComponent(component, adjacencyEdges, blockEdges, stubEdges,
                stubAndBlockEdgesSet);
        stList_appendAll(singleStubEdgeComponents, splitComponents);
        stList_setDestructor(splitComponents, NULL); //Do this to avoid destroying the underlying lists
        stList_destruct(splitComponents);
    }
    stList_destruct(components);

    /*
     * Remove the stub/block edges from the components.
     */
    stList *adjacencyOnlyComponents = filterLists(singleStubEdgeComponents, stubAndBlockEdgesSet);
    stList_destruct(singleStubEdgeComponents);
    stSortedSet_destruct(stubAndBlockEdgesSet);

    /*
     * Merge the adjacency edges in the components into a single list.
     */
    stList *updatedChosenEdges = joinLists(adjacencyOnlyComponents);
    stList_destruct(adjacencyOnlyComponents);

    return updatedChosenEdges;
}

static void checkEdges(stList *edges, int32_t nodeNumber) {
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        assert(stIntTuple_getPosition(edge, 0) >= 0);
        assert(stIntTuple_getPosition(edge, 0) < nodeNumber);
        assert(stIntTuple_getPosition(edge, 1) >= 0);
        assert(stIntTuple_getPosition(edge, 1) < nodeNumber);
    }
}

stList *chooseMatching(uint32_t nodeNumber, stList *adjacencyEdges, stList *stubEdges, stList *blockEdges,
        bool makeStubCyclesDisjoint, stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    /*
     * Check the inputs.
     */
    assert(nodeNumber % 2 == 0);
    assert(stList_length(stubEdges) > 0);
    assert(stList_length(stubEdges) + stList_length(blockEdges) == nodeNumber / 2);
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
    stList *updatedChosenEdges = mergeComponents2(chosenEdges, adjacencyEdges, stubEdges, blockEdges);
    stList_destruct(chosenEdges);
    chosenEdges = updatedChosenEdges;

    /*
     * Split stub components.
     */
    if (makeStubCyclesDisjoint) {
        updatedChosenEdges = splitMultipleStubCycles(chosenEdges, adjacencyEdges, stubEdges, blockEdges);
        stList_destruct(chosenEdges);
        chosenEdges = updatedChosenEdges;
    }

    return chosenEdges;
}
