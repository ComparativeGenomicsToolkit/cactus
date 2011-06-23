#include "cactus.h"
#include "sonLib.h"

////////////////////////////////////
////////////////////////////////////
//Misc supporting functions
////////////////////////////////////
////////////////////////////////////

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

static stList *filter(stList *list, stSortedSet *set, bool include) {
    /*
     * Returns a new list, either containing the intersection with set if include is non-zero,
     * or containing the set difference if include is zero.
     */
    stList *list2 = stList_construct();
    for (int32_t i = 0; i < stList_length(list); i++) {
        void *o = stList_get(list, i);
        if ((stSortedSet_search(set, o) != NULL && include)
                || (stSortedSet_search(set, o) == NULL && !include)) {
            stList_append(list2, o);
        }
    }
    return list2;
}

static stList *filterToExclude(stList *list, stSortedSet *set) {
    /*
     * Returns a new list, identical to list, but with any elements contained in set removed.
     */
    return filter(list, set, 0);
}

stList *filterToInclude(stList *list, stSortedSet *set) {
    /*
     * Returns a new list, identical to list, but with any elements not contained in set removed.
     */
    return filter(list, set, 1);
}

static stList *filterListsToExclude(stList *listOfLists, stSortedSet *set) {
    /*
     * Takes a list of lists and returns a new list of lists whose elements are the product of applying
     * filterToExclude to each member of listOfLists in the same order.
     */
    stList *listOfLists2 = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(listOfLists); i++) {
        stList_append(listOfLists2,
                filterToExclude(stList_get(listOfLists, i), set));
    }
    return listOfLists2;
}

static void getNodesToEdgesHashP(stHash *nodesToEdges, stIntTuple *edge,
        int32_t position) {
    stIntTuple *node = stIntTuple_construct(1,
            stIntTuple_getPosition(edge, position));
    stList *edges;
    if ((edges = stHash_search(nodesToEdges, node)) == NULL) {
        edges = stList_construct();
        stHash_insert(nodesToEdges, node, edges);
    } else {
        stIntTuple_destruct(node);
    }
    stList_append(edges, edge);
}

stHash *getNodesToEdgesHash(stList *edges) {
    /*
     * Build a hash of nodes to edges.
     */

    stHash *nodesToEdges = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn,
            (void(*)(void *)) stIntTuple_destruct,
            (void(*)(void *)) stList_destruct);

    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        getNodesToEdgesHashP(nodesToEdges, edge, 0);
        getNodesToEdgesHashP(nodesToEdges, edge, 1);
    }

    return nodesToEdges;
}

static stList *getEdgesForNode(int32_t node, stHash *nodesToAdjacencyEdges) {
    /*
     * Get edges for given node, or NULL if none.
     */
    stIntTuple *node1 = stIntTuple_construct(1, node);
    stList *edges = stHash_search(nodesToAdjacencyEdges, node1);
    stIntTuple_destruct(node1);
    return edges;
}

static int32_t getOtherPosition(stIntTuple *edge, int32_t node) {
    /*
     * Get the other position to the given node in the edge.
     */
    if (stIntTuple_getPosition(edge, 0) == node) {
        return stIntTuple_getPosition(edge, 1);
    } else {
        assert(stIntTuple_getPosition(edge, 1) == node);
        return stIntTuple_getPosition(edge, 0);
    }
}

stIntTuple *getEdgeForNodes(int32_t node1, int32_t node2,
        stHash *nodesToAdjacencyEdges) {
    /*
     * Searches for any edge present in nodesToAdjacencyEdges that links node and node2.
     */
    stList *edges = getEdgesForNode(node1, nodesToAdjacencyEdges);
    if(edges == NULL) {
        return NULL;
    }
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        if (getOtherPosition(edge, node1) == node2) {
            return edge;
        }
    }
    return NULL;
}

static stSortedSet *getSetOfMergedLists(stList *list1, stList *list2) {
    /*
     * Returns a sorted set of two input lists.
     */
    stList *list3 = stList_copy(list1, NULL);
    stList_appendAll(list3, list2);
    stSortedSet *set = stList_getSortedSet(list3, NULL);
    stList_destruct(list3);
    return set;
}

stList *getStubAndChainEdgeFreeComponents(stList *listOfLists, stList *list1,
        stList *list2) {
    /*
     * Returns a new list of lists, with each sub list filtered to exclude members of list1 and list2.
     */
    stSortedSet *set = getSetOfMergedLists(list1, list2);
    stList *filteredListOfLists = filterListsToExclude(listOfLists, set);
    stSortedSet_destruct(set);

    return filteredListOfLists;
}

static stList *joinLists(stList *listOfLists) {
    /*
     * Returns new list which contains elements of the list of list concatenated in one list.
     */
    stList *joinedList = stList_construct();
    for (int32_t i = 0; i < stList_length(listOfLists); i++) {
        stList_appendAll(joinedList, stList_get(listOfLists, i));
    }
    return joinedList;
}

////////////////////////////////////
////////////////////////////////////
//Functions to compute connected components.
////////////////////////////////////
////////////////////////////////////

static void getComponentsP(stHash *nodesToEdges, int32_t node,
        stSortedSet *component) {
    stIntTuple *key = stIntTuple_construct(1, node);
    stList *edges = stHash_search(nodesToEdges, key);
    if (edges != NULL) {
        stHash_remove(nodesToEdges, key);
        for (int32_t i = 0; i < stList_length(edges); i++) {
            stIntTuple *edge = stList_get(edges, i);
            if (stSortedSet_search(component, edge) == NULL) {
                stSortedSet_insert(component, edge);
            }
            /*
             * Recursion on stack could equal the total number of nodes.
             */
            getComponentsP(nodesToEdges, stIntTuple_getPosition(edge, 0),
                    component);
            getComponentsP(nodesToEdges, stIntTuple_getPosition(edge, 1),
                    component);
        }
        stList_destruct(edges);
    }
    stIntTuple_destruct(key);
}

stList *getComponents(stList *edges) {
    /*
     * Gets a list of connected components, each connected component
     * being represented as a list of the edges, such that each edge is in exactly one
     * connected component. Allows for multi-graphs (multiple edges connecting two nodes).
     */

    stHash *nodesToEdges = getNodesToEdgesHash(edges);

    /*
     * Traverse the edges greedily
     */
    stList *components =
            stList_construct3(0, (void(*)(void *)) stList_destruct);
    stList *nodes = stHash_getKeys(nodesToEdges);
    while (stList_length(nodes) > 0) {
        stIntTuple *node = stList_pop(nodes);
        stList *edges = stHash_search(nodesToEdges, node);
        if (edges != NULL) { //We have a component to build
            stSortedSet *component = stSortedSet_construct();
            stHash_remove(nodesToEdges, node);
            for (int32_t i = 0; i < stList_length(edges); i++) {
                stIntTuple *edge = stList_get(edges, i);
                getComponentsP(nodesToEdges, stIntTuple_getPosition(edge, 0),
                        component);
                getComponentsP(nodesToEdges, stIntTuple_getPosition(edge, 1),
                        component);
            }
            stList_append(components, stSortedSet_getList(component));
            //Cleanup
            stSortedSet_destruct(component);
            stList_destruct(edges);
        }
        stIntTuple_destruct(node);
    }
    assert(stHash_size(nodesToEdges) == 0);
    stHash_destruct(nodesToEdges);
    stList_destruct(nodes);

    return components;
}

stList *getComponents2(stList *adjacencyEdges, stList *stubEdges,
        stList *chainEdges) {
    /*
     * Gets a list of connected components for a set of adjacency, stub and chain edges.
     * If adjacencyEdges, stubEdges or chainEdges are NULL then they are ignored.
     */
    stList *allEdges = stList_construct(); //Build a concatenated list of all the chain, stub and adjacency edges.
    if (adjacencyEdges != NULL) {
        stList_appendAll(allEdges, adjacencyEdges);
    }
    if (stubEdges != NULL) {
        stList_appendAll(allEdges, stubEdges);
    }
    if (chainEdges != NULL) {
        stList_appendAll(allEdges, chainEdges);
    }
    stList *components = getComponents(allEdges); //Gets the graph components.
    stList_destruct(allEdges); //Cleanup the all edges.
    return components;
}

////////////////////////////////////
////////////////////////////////////
//Functions to compute adjacency edge switches between cycles
////////////////////////////////////
////////////////////////////////////

typedef struct _AdjacencySwitch {
    /*
     * Structure  to represent an adjacency edge switch, in which two edges are removed and replaced with two others.
     */
    stIntTuple *oldEdge1;
    stIntTuple *oldEdge2;
    stIntTuple *newEdge1;
    stIntTuple *newEdge2;
    int32_t cost;
} AdjacencySwitch;

static AdjacencySwitch *adjacencySwitch_construct(stIntTuple *oldEdge1,
        stIntTuple *oldEdge2, stIntTuple *newEdge1, stIntTuple *newEdge2,
        int32_t cost) {
    AdjacencySwitch *adjacencySwitch = st_malloc(sizeof(AdjacencySwitch));
    adjacencySwitch->oldEdge1 = oldEdge1;
    adjacencySwitch->oldEdge2 = oldEdge2;
    adjacencySwitch->newEdge1 = newEdge1;
    adjacencySwitch->newEdge2 = newEdge2;
    assert(oldEdge1 != NULL);
    assert(oldEdge2 != NULL);
    assert(newEdge1 != NULL);
    assert(newEdge2 != NULL);
    adjacencySwitch->cost = cost;
    return adjacencySwitch;
}

static void adjacencySwitch_destruct(AdjacencySwitch *adjacencySwitch) {
    free(adjacencySwitch);
}

static AdjacencySwitch *adjacencySwitch_update(
        AdjacencySwitch *adjacencySwitch, stIntTuple *oldEdge1,
        stIntTuple *oldEdge2, stIntTuple *newEdge1, stIntTuple *newEdge2,
        int32_t cost) {
    /*
     * Convenience function, replaces one adjacency edge with a new one if the cost is lower.
     */
    if (adjacencySwitch == NULL) {
        return adjacencySwitch_construct(oldEdge1, oldEdge2, newEdge1,
                newEdge2, cost);
    }
    if (adjacencySwitch->cost > cost) {
        adjacencySwitch_destruct(adjacencySwitch);
        return adjacencySwitch_construct(oldEdge1, oldEdge2, newEdge1,
                newEdge2, cost);
    }
    return adjacencySwitch;
}

static stList *getValidEdges(int32_t node, stHash *nodesToAdjacencyEdges,
        stSortedSet *validEdges) {
    /*
     * Iterates through the edges for a node and returns those in the set of the valid edges.
     */
    stList *edges = getEdgesForNode(node, nodesToAdjacencyEdges);
    return filterToInclude(edges, validEdges);
}

static AdjacencySwitch *getBest4EdgeAdjacencySwitchP(stIntTuple *oldEdge1,
        int32_t node1, stHash *nodesToAdjacencyEdges,
        stSortedSet *allCurrentEdgesSet, stSortedSet *bridgingAdjacencyEdges) {
    /*
     * Returns the best adjacency switch for the given node and edge that
     * contains 4 existing edges.
     */
    int32_t node4 = getOtherPosition(oldEdge1, node1);
    AdjacencySwitch *minimumCostAdjacencySwitch = NULL;
    stList *validEdges = getValidEdges(node1, nodesToAdjacencyEdges,
            bridgingAdjacencyEdges);
    for (int32_t i = 0; i < stList_length(validEdges); i++) {
        stIntTuple *newEdge1 = stList_get(validEdges, i);
        int32_t node2 = getOtherPosition(newEdge1, node1);
        stList *validEdges2 = getValidEdges(node2, nodesToAdjacencyEdges,
                allCurrentEdgesSet);
        for (int32_t j = 0; j < stList_length(validEdges2); j++) {
            stIntTuple *oldEdge2 = stList_get(validEdges2, j);
            int32_t node3 = getOtherPosition(oldEdge2, node2);
            stIntTuple *newEdge2 = getEdgeForNodes(node3, node4,
                    nodesToAdjacencyEdges);
            assert(newEdge2 != NULL);
            int32_t cost = stIntTuple_getPosition(oldEdge1, 2)
                    + stIntTuple_getPosition(oldEdge2, 2)
                    - stIntTuple_getPosition(newEdge1, 2)
                    - stIntTuple_getPosition(newEdge2, 2);
            minimumCostAdjacencySwitch = adjacencySwitch_update(
                    minimumCostAdjacencySwitch, oldEdge1, oldEdge2, newEdge1,
                    newEdge2, cost);
        }
        stList_destruct(validEdges2);
    }
    stList_destruct(validEdges);
    assert(minimumCostAdjacencySwitch != NULL);
    return minimumCostAdjacencySwitch;
}

static AdjacencySwitch *getMinimumCostAdjacencySwitch(
        AdjacencySwitch *adjacencySwitch1, AdjacencySwitch *adjacencySwitch2) {
    /*
     * Convenience function that returns the adjacency switch of the two with lower cost,
     * allowing for NULL values and destroying
     * the higher cost adjacency switch in the process.
     */
    if (adjacencySwitch1 == NULL) {
        return adjacencySwitch2;
    }
    if (adjacencySwitch2 == NULL) {
        return adjacencySwitch1;
    }
    if (adjacencySwitch1->cost < adjacencySwitch2->cost) {
        adjacencySwitch_destruct(adjacencySwitch2);
        return adjacencySwitch1;
    }
    adjacencySwitch_destruct(adjacencySwitch1);
    return adjacencySwitch2;
}

static AdjacencySwitch *getBest4EdgeAdjacencySwitch2(stIntTuple *oldEdge1,
        stHash *nodesToAdjacencyEdges, stSortedSet *allCurrentEdgesSet,
        stSortedSet *bridgingAdjacencyEdges) {
    /*
     * Returns the best 3 or 4 edge switch (one including 3 or 4 edges) for the given existing edge, if they exist, or else NULL.
     */
    return getMinimumCostAdjacencySwitch(
            getBest4EdgeAdjacencySwitchP(oldEdge1,
                    stIntTuple_getPosition(oldEdge1, 0), nodesToAdjacencyEdges,
                    allCurrentEdgesSet, bridgingAdjacencyEdges),
            getBest4EdgeAdjacencySwitchP(oldEdge1,
                    stIntTuple_getPosition(oldEdge1, 1), nodesToAdjacencyEdges,
                    allCurrentEdgesSet, bridgingAdjacencyEdges));
}

static void getNodeSetOfEdgesP(stSortedSet *nodes, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    if (stSortedSet_search(nodes, i) == NULL) {
        stSortedSet_insert(nodes, i);
    } else {
        stIntTuple_destruct(i);
    }
}

stSortedSet *getNodeSetOfEdges(stList *edges) {
    /*
     * Returns a sorted set of the nodes in the given list of edges.
     */
    stSortedSet *nodes = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        getNodeSetOfEdgesP(nodes, stIntTuple_getPosition(edge, 0));
        getNodeSetOfEdgesP(nodes, stIntTuple_getPosition(edge, 1));
    }
    return nodes;
}

bool nodeInSet(stSortedSet *nodes, int32_t node) {
    /*
     * Returns non zero iff the node is in the set.
     */
    stIntTuple *i = stIntTuple_construct(1, node);
    bool b = stSortedSet_search(nodes, i) != NULL;
    stIntTuple_destruct(i);
    return b;
}

void addNodeToSet(stSortedSet *nodes, int32_t node) {
    /*
     * Adds a node to the set.
     */
    assert(!nodeInSet(nodes, node));
    stIntTuple *i = stIntTuple_construct(1, node);
    stSortedSet_insert(nodes, i);
}

static stSortedSet *getEdgesThatBridgeComponents(stList *components,
        stHash *nodesToAdjacencyEdges) {
    /*
     * Get set of adjacency edges that bridge between (have a node in two) components.
     */

    stSortedSet *bridgingAdjacencyEdges = stSortedSet_construct();

    for (int32_t i = 0; i < stList_length(components); i++) {
        stSortedSet *componentNodes = getNodeSetOfEdges(
                stList_get(components, i));
        stSortedSetIterator *it = stSortedSet_getIterator(componentNodes);
        stIntTuple *node;
        while ((node = stSortedSet_getNext(it)) != NULL) {
            stList *edges = stHash_search(nodesToAdjacencyEdges, node);
            assert(edges != NULL);
            for (int32_t j = 0; j < stList_length(edges); j++) {
                stIntTuple *edge = stList_get(edges, j);
                stIntTuple *node1 = stIntTuple_construct(1,
                        stIntTuple_getPosition(edge, 0));
                stIntTuple *node2 = stIntTuple_construct(1,
                        stIntTuple_getPosition(edge, 1));
                assert(
                        stSortedSet_search(componentNodes, node1) != NULL
                                || stSortedSet_search(componentNodes, node2)
                                        != NULL);
                if (stSortedSet_search(componentNodes, node1) == NULL
                        || stSortedSet_search(componentNodes, node2) == NULL) {
                    stSortedSet_insert(bridgingAdjacencyEdges, edge);
                }
                stIntTuple_destruct(node1);
                stIntTuple_destruct(node2);
            }
        }
        stSortedSet_destructIterator(it);
        stSortedSet_destruct(componentNodes);
    }

    return bridgingAdjacencyEdges;
}

static AdjacencySwitch *getBestAdjacencySwitch(stList *cycles,
        stHash *nodesToAdjacencyEdges) {
    /*
     * Returns the best 3 or 4 edge switch (one including 3 or 4 edges) for the given existing edge, if they exist, or else NULL.
     */
    assert(stList_length(cycles) > 0);

    stList *allComponentEdges = joinLists(cycles);
    assert(stList_length(allComponentEdges) > 0);
    stSortedSet *allComponentEdgesSet = stList_getSortedSet(allComponentEdges,
            NULL);
    assert(
            stList_length(allComponentEdges) == stSortedSet_size(
                    allComponentEdgesSet));

    /*
     * Get list of adjacency edges that bridge between (have a node in two) components.
     */
    stSortedSet *bridgingAdjacencyEdges = getEdgesThatBridgeComponents(cycles,
            nodesToAdjacencyEdges);
    assert(stList_length(allComponentEdges) > 0);

    /*
     * Look for the best 3 or 4 edge switch.
     */
    AdjacencySwitch *minimumCostAdjacencySwitch = NULL;
    for (int32_t i = 0; i < stList_length(allComponentEdges); i++) {
        minimumCostAdjacencySwitch = getMinimumCostAdjacencySwitch(
                minimumCostAdjacencySwitch,
                getBest4EdgeAdjacencySwitch2(stList_get(allComponentEdges, i),
                        nodesToAdjacencyEdges, allComponentEdgesSet,
                        bridgingAdjacencyEdges));
    }
    assert(minimumCostAdjacencySwitch != NULL);

    /*
     * Cleanup
     */
    stList_destruct(allComponentEdges);
    stSortedSet_destruct(allComponentEdgesSet);
    stSortedSet_destruct(bridgingAdjacencyEdges);

    return minimumCostAdjacencySwitch;
}

////////////////////////////////////
////////////////////////////////////
//Functions to merge disjoint cycles
////////////////////////////////////
////////////////////////////////////

static void doBestMergeOfTwoSimpleCycles(stList *cycles,
        stHash *nodesToAdjacencyEdges, stList *adjacencyEdges) {
    /*
     * Merge two simple cycles, using the best possible adjacency switch. Modifies components list in place,
     * destroying two old components and adding a new one. If new adjacency edges are needed then they are
     * added to the adjacency edges list.
     */
    assert(stList_length(cycles) > 1);

    /*
     * Get the best adjacency switch.
     */
    AdjacencySwitch *adjacencySwitch = getBestAdjacencySwitch(cycles,
            nodesToAdjacencyEdges);
    assert(adjacencySwitch != NULL);

    /*
     * Find the two components to merge.
     */
    stList *cyclesToMerge = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(cycles); i++) {
        stList *cycle = stList_get(cycles, i);
        if (stList_contains(cycle, adjacencySwitch->oldEdge1)) {
            assert(!stList_contains(cycle, adjacencySwitch->oldEdge2));
            stList_append(cyclesToMerge, cycle);
        } else if (stList_contains(cycle, adjacencySwitch->oldEdge2)) {
            stList_append(cyclesToMerge, cycle);
        }
    }

    /*
     * Now construct the new component and modify the list of components in place.
     */
    assert(stList_length(cyclesToMerge) == 2);
    stList *newComponent = joinLists(cyclesToMerge);
    assert(!stList_contains(newComponent, NULL));
    //Cleanup the old components
    assert(stList_contains(cycles, stList_get(cyclesToMerge, 0)));
    stList_removeItem(cycles, stList_get(cyclesToMerge, 0));
    assert(stList_contains(cycles, stList_get(cyclesToMerge, 1)));
    stList_removeItem(cycles, stList_get(cyclesToMerge, 1));
    stList_destruct(cyclesToMerge);
    //Now remove the old edges and add the new ones
    assert(stList_contains(newComponent, adjacencySwitch->oldEdge1));
    stList_removeItem(newComponent, adjacencySwitch->oldEdge1);
    assert(stList_contains(newComponent, adjacencySwitch->oldEdge2));
    stList_removeItem(newComponent, adjacencySwitch->oldEdge2);
    assert(!stList_contains(newComponent, adjacencySwitch->newEdge1));
    stList_append(newComponent, adjacencySwitch->newEdge1);
    assert(!stList_contains(newComponent, adjacencySwitch->newEdge2));
    stList_append(newComponent, adjacencySwitch->newEdge2);
    adjacencySwitch_destruct(adjacencySwitch); //Clean the adjacency switch.
    //Finally add the component to the list of components
    stList_append(cycles, newComponent);
}

stList *mergeSimpleCycles(stList *cycles, stList *adjacencyEdges) {
    /*
     * Takes a set of simple cycles (containing only the adjacency edges).
     * Returns a single simple cycle, as a list of edges, by doing length(components)-1
     * calls to doBestMergeOfTwoSimpleCycles.
     */

    /*
     * Build a hash of nodes to adjacency edges.
     */
    stHash *nodesToAdjacencyEdges = getNodesToEdgesHash(adjacencyEdges);

    cycles = stList_copy(cycles, NULL);
    for (int32_t i = 0; i < stList_length(cycles); i++) { //Clone the complete list
        assert(stList_length(stList_get(cycles, i)) > 0);
        assert(!stList_contains(stList_get(cycles, i), NULL));
        stList_set(cycles, i, stList_copy(stList_get(cycles, i), NULL));
    }
    while (stList_length(cycles) > 1) {
        doBestMergeOfTwoSimpleCycles(cycles, nodesToAdjacencyEdges,
                adjacencyEdges);
    }
    assert(stList_length(cycles) == 1);
    stList *mergedComponent = stList_get(cycles, 0);
    stList_destruct(cycles);
    stHash_destruct(nodesToAdjacencyEdges);
    return mergedComponent;
}

static stList *mergeSimpleCycles2(stList *chosenEdges, stList *adjacencyEdges,
        stList *stubEdges, stList *chainEdges) {
    /*
     * Returns a new set of chosen edges, modified by adjacency switches such that every simple cycle
     * contains at least one stub edge.
     */

    /*
     * Calculate components.
     */
    stList *components = getComponents2(chosenEdges, stubEdges, chainEdges);

    /*
     * Divide the components by the presence of one or more stub edges.
     */
    stSortedSet *stubEdgesSet = stList_getSortedSet(stubEdges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);

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
     * Remove the stub/chain edges from the components.
     */
    stList_append(stubFreeComponents, globalComponent);
    stList *adjacencyOnlyComponents = getStubAndChainEdgeFreeComponents(
            stubFreeComponents, stubEdges, chainEdges);

    stList_destruct(stubFreeComponents);
    stList_destruct(globalComponent);
    stList_destruct(components); //We only clean this up now, as this frees the components it contains.

    /*
     * Merge stub free components into the others.
     */
    stList *updatedChosenEdges = mergeSimpleCycles(adjacencyOnlyComponents,
            adjacencyEdges);
    stList_destruct(adjacencyOnlyComponents);

    return updatedChosenEdges;
}

////////////////////////////////////
////////////////////////////////////
//Functions to split cycles with multiple stubs in them.
////////////////////////////////////
////////////////////////////////////

static stSortedSet *getOddNodes(stList *cycle) {
    /*
     * Returns alternating nodes in a simple cycle.
     */

    //Set to return
    stSortedSet *nodes = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);

    stHash *nodesToEdges = getNodesToEdgesHash(cycle);
    int32_t node = stIntTuple_getPosition(stList_get(cycle, 0), 0);
    int32_t pNode = -1;
    int32_t counter = 0;
    bool b = 1;
    assert(stList_length(cycle) % 2 == 0);
    while (counter++ < stList_length(cycle)) {
        if (b) { //Make alternating
            addNodeToSet(nodes, node);
            b = 0;
        } else {
            b = 1;
        }
        stList *edges = getEdgesForNode(node, nodesToEdges);
        assert(stList_length(edges) == 2);
        stIntTuple *edge = stList_get(edges, 0);
        int32_t node2 = getOtherPosition(edge, node);
        if (node2 != pNode) {
            pNode = node;
            node = node2;
            continue;
        }
        edge = stList_get(edges, 1);
        node2 = getOtherPosition(edge, node);
        assert(node2 != pNode);
        pNode = node;
        node = node2;
    }
    stHash_destruct(nodesToEdges);

    assert(stList_length(cycle) / 2 == stSortedSet_size(nodes));

    return nodes;
}

static stList *getOddToEvenAdjacencyEdges(stSortedSet *oddNodes,
        stList *adjacencyEdges) {
    /*
     * Gets edges that include one node in the set of oddNodes, but not both.
     */
    stList *oddToEvenAdjacencyEdges = stList_construct();
    for (int32_t i = 0; i < stList_length(adjacencyEdges); i++) {
        stIntTuple *edge = stList_get(adjacencyEdges, i);
        if (nodeInSet(oddNodes, stIntTuple_getPosition(edge, 0)) ^ nodeInSet(
                oddNodes, stIntTuple_getPosition(edge, 1))) {
            stList_append(oddToEvenAdjacencyEdges, edge);
        }
    }
    return oddToEvenAdjacencyEdges;
}

static void splitIntoAdjacenciesStubsAndChains(stList *subCycle,
        stList *adjacencyEdges, stList *stubEdges, stList *chainEdges,
        stList **subAdjacencyEdges, stList **subStubEdges,
        stList **subChainEdges) {
    /*
     * Splits run into cycles and chains..
     */
    *subStubEdges = stList_construct();
    *subChainEdges = stList_construct();
    for (int32_t j = 0; j < stList_length(subCycle); j++) {
        stIntTuple *edge = stList_get(subCycle, j);
        if (stList_contains(stubEdges, edge)) {
            stList_append(*subStubEdges, edge);
        } else if (stList_contains(chainEdges, edge)) {
            stList_append(*subChainEdges, edge);
        }
    }
    *subAdjacencyEdges = stList_construct();
    stSortedSet *nodes = getNodeSetOfEdges(subCycle);
    for (int32_t j = 0; j < stList_length(adjacencyEdges); j++) {
        stIntTuple *edge = stList_get(adjacencyEdges, j);
        if (nodeInSet(nodes, stIntTuple_getPosition(edge, 0)) && nodeInSet(
                nodes, stIntTuple_getPosition(edge, 1))) {
            stList_append(*subAdjacencyEdges, edge);
        }
    }
    stSortedSet_destruct(nodes);
}

static stList *splitMultipleStubCycle(stList *cycle, stList *adjacencyEdges,
        stList *stubEdges, stList *chainEdges) {
    /*
     *  Takes a simple cycle containing k stub edges and splits into k cycles, each containing 1 stub edge.
     */

    /*
     * Get sub-components containing only adjacency and chain edges.
     */

    stSortedSet *stubAndChainEdgesSet = getSetOfMergedLists(stubEdges,
            chainEdges);
    stList *adjacencyEdgeMatching =
            filterToExclude(cycle, stubAndChainEdgesSet); //Filter out the the non-adjacency edges
    //Make it only the chain edges present in the original component
    stList *stubFreePaths = getComponents2(adjacencyEdgeMatching, NULL,
            chainEdges);
    stList_destruct(adjacencyEdgeMatching);
    assert(stList_length(stubFreePaths) >= 1);

    stList *splitCycles = stList_construct3(0,
            (void(*)(void *)) stList_destruct); //The list to return.


    if (stList_length(stubFreePaths) > 1) {
        /*
         * Build the list of adjacency edges acceptable in the merge
         */
        stSortedSet *oddNodes = getOddNodes(cycle);
        stList *oddToEvenAdjacencyEdges = getOddToEvenAdjacencyEdges(oddNodes,
                adjacencyEdges);
        stHash *nodesToOddToEvenAdjacencyEdges = getNodesToEdgesHash(
                oddToEvenAdjacencyEdges);

        /*
         * Merge together the best two components.
         */
        stList *l = filterListsToExclude(stubFreePaths, stubAndChainEdgesSet);
        doBestMergeOfTwoSimpleCycles(l, nodesToOddToEvenAdjacencyEdges,
                oddToEvenAdjacencyEdges); //This is inplace.
        stList *l2 = joinLists(l);
        stList_destruct(l);
        l = getComponents2(l2, stubEdges, chainEdges);
        assert(stList_length(l) == 2);
        stList_destruct(l2);

        /*
         * Cleanup
         */
        stSortedSet_destruct(oddNodes);
        stHash_destruct(nodesToOddToEvenAdjacencyEdges);
        stList_destruct(oddToEvenAdjacencyEdges);

        /*
         * Call procedure recursively.
         */
        for (int32_t i = 0; i < stList_length(l); i++) {
            /*
             * Split into adjacency edges, stub edges and chain edges.
             */
            stList *subCycle = stList_get(l, i);
            stList *subAdjacencyEdges;
            stList *subStubEdges;
            stList *subChainEdges;
            splitIntoAdjacenciesStubsAndChains(subCycle, adjacencyEdges,
                    stubEdges, chainEdges, &subAdjacencyEdges, &subStubEdges,
                    &subChainEdges);

            /*
             * Call recursively.
             */
            l2 = splitMultipleStubCycle(subCycle, subAdjacencyEdges,
                    subStubEdges, subChainEdges);
            stList_appendAll(splitCycles, l2);

            /*
             * Clean up
             */
            stList_setDestructor(l2, NULL);
            stList_destruct(l2);
            stList_destruct(subAdjacencyEdges);
            stList_destruct(subStubEdges);
            stList_destruct(subChainEdges);
        }
        stList_destruct(l);
    } else {
        stList_append(splitCycles, stList_copy(cycle, NULL));
    }

    stSortedSet_destruct(stubAndChainEdgesSet);
    stList_destruct(stubFreePaths);

    return splitCycles;
}

stList *splitMultipleStubCycles(stList *chosenEdges, stList *adjacencyEdges,
        stList *stubEdges, stList *chainEdges) {
    /*
     *  Returns an updated list of adjacency edges, such that each stub edge is a member of exactly one cycle.
     */

    /*
     * Calculate components.
     */
    stList *cycles = getComponents2(chosenEdges, stubEdges, chainEdges);

    /*
     * Find components with multiple stub edges.
     */
    stList *singleStubEdgeCycles = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(cycles); i++) {
        stList *subCycle = stList_get(cycles, i);
        stList *subAdjacencyEdges;
        stList *subStubEdges;
        stList *subChainEdges;
        splitIntoAdjacenciesStubsAndChains(subCycle, adjacencyEdges, stubEdges,
                chainEdges, &subAdjacencyEdges, &subStubEdges, &subChainEdges);
        stList *splitCycles = splitMultipleStubCycle(subCycle,
                subAdjacencyEdges, subStubEdges, subChainEdges);
        stList_appendAll(singleStubEdgeCycles, splitCycles);
        stList_setDestructor(splitCycles, NULL); //Do this to avoid destroying the underlying lists
        stList_destruct(splitCycles);
        stList_destruct(subAdjacencyEdges);
        stList_destruct(subStubEdges);
        stList_destruct(subChainEdges);
    }
    stList_destruct(cycles);

    /*
     * Remove the stub/chain edges from the components.
     */
    stSortedSet *stubAndChainEdgesSet = getSetOfMergedLists(stubEdges,
            chainEdges);
    stList *adjacencyOnlyComponents = filterListsToExclude(
            singleStubEdgeCycles, stubAndChainEdgesSet);
    stList_destruct(singleStubEdgeCycles);
    stSortedSet_destruct(stubAndChainEdgesSet);

    /*
     * Merge the adjacency edges in the components into a single list.
     */
    stList *updatedChosenEdges = joinLists(adjacencyOnlyComponents);
    stList_destruct(adjacencyOnlyComponents);

    return updatedChosenEdges;
}

////////////////////////////////////
////////////////////////////////////
//Top level function
////////////////////////////////////
////////////////////////////////////

static bool edgeInSet(stSortedSet *edges, int32_t node1, int32_t node2) {
    /*
     * Returns non-zero iff the edge linking node1 and node2 is in the set of edges.
     */
    if (node1 > node2) {
        return edgeInSet(edges, node2, node1);
    }
    stIntTuple *edge = stIntTuple_construct(2, node1, node2);
    bool b = stSortedSet_search(edges, edge) != NULL;
    stIntTuple_destruct(edge);
    return b;
}

static void addEdgeToSet(stSortedSet *edges, int32_t node1, int32_t node2) {
    /*
     * Adds an edge to the set.
     */
    if (node1 > node2) {
        addEdgeToSet(edges, node2, node1);
        return;
    }
    assert(!edgeInSet(edges, node1, node2));
    stIntTuple *edge = stIntTuple_construct(2, node1, node2);
    stSortedSet_insert(edges, edge);
    assert(edgeInSet(edges, node1, node2));
}

static void checkEdges(stList *edges, int32_t nodeNumber, bool coversAllNodes,
        bool isClique) {
    /*
     * Check the edges all refer to nodes between 0 and node number. If length three, check final argument
     * (which is a weight), is greater than or equal to zero. Check that if coversAllNodes is non-zero, all nodes are
     * covered. Checks edges from clique if isClique is non-zero.
     */
    stSortedSet *edgesSeen = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        /*
         * Check edges connect actual nodes.
         */
        assert(stIntTuple_getPosition(edge, 0) >= 0);
        assert(stIntTuple_getPosition(edge, 0) < nodeNumber);
        assert(stIntTuple_getPosition(edge, 1) >= 0);
        assert(stIntTuple_getPosition(edge, 1) < nodeNumber);
        /*
         * Check is not a multi-graph.
         */
        assert(
                !edgeInSet(edgesSeen, stIntTuple_getPosition(edge, 0),
                        stIntTuple_getPosition(edge, 1)));
        addEdgeToSet(edgesSeen, stIntTuple_getPosition(edge, 0),
                stIntTuple_getPosition(edge, 1));
        /*
         * Check weight, if weighted.
         */
        if (stIntTuple_length(edge) == 3) {
            assert(stIntTuple_getPosition(edge, 2) >= 0);
        } else {
            assert(stIntTuple_length(edge) == 2);
        }
    }
    stSortedSet_destruct(edgesSeen);
    /*
     * Check edges connect all nodes.
     */
    if (coversAllNodes) {
        stSortedSet *nodes = getNodeSetOfEdges(edges);
        assert(stSortedSet_size(nodes) == nodeNumber);
        stSortedSet_destruct(nodes);
    }
    /*
     * Check is clique.
     */
    if (isClique) {
        assert(
                stList_length(edges) == (nodeNumber * nodeNumber - nodeNumber)
                        / 2);
    }
}

static void makeMatchingPerfect(stList *chosenEdges, stList *adjacencyEdges,
        int32_t nodeNumber) {
    /*
     * While the the number of edges is less than a perfect matching add random edges.
     */
    stSortedSet *nodeSet = getNodeSetOfEdges(chosenEdges);
    stHash *nodesToAdjacencyEdges = getNodesToEdgesHash(adjacencyEdges);
    int32_t pNode = -1;
    for (int32_t i = 0; i < nodeNumber; i++) {
        if (!nodeInSet(nodeSet, i)) {
            if (pNode == -1) {
                pNode = i;
            } else {
                stList_append(chosenEdges,
                        getEdgeForNodes(pNode, i, nodesToAdjacencyEdges));
                pNode = -1;
            }
        }
    }
    assert(pNode == -1);
    stSortedSet_destruct(nodeSet);
    assert(stList_length(chosenEdges) * 2 == nodeNumber);
    stHash_destruct(nodesToAdjacencyEdges);
}

void checkInputs(uint32_t nodeNumber, stList *adjacencyEdges,
        stList *stubEdges, stList *chainEdges) {
    /*
     * Checks the inputs to the algorithm are as expected.
     */
    assert(nodeNumber % 2 == 0);
    if(nodeNumber > 0) {
        assert(stList_length(stubEdges) > 0);
    }
    assert(
            stList_length(stubEdges) + stList_length(chainEdges) == (nodeNumber
                    / 2));
    checkEdges(stubEdges, nodeNumber, 0, 0);
    checkEdges(chainEdges, nodeNumber, 0, 0);
    stList *stubsAndChainEdges = stList_copy(stubEdges, NULL);
    stList_appendAll(stubsAndChainEdges, chainEdges);
    checkEdges(stubsAndChainEdges, nodeNumber, 1, 0);
    stList_destruct(stubsAndChainEdges);
    checkEdges(adjacencyEdges, nodeNumber, 1, 1);
}

stList *getMatchingWithCyclicConstraints(uint32_t nodeNumber,
        stList *adjacencyEdges, stList *stubEdges, stList *chainEdges,
        bool makeStubCyclesDisjoint,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    /*
     * Check the inputs.
     */
    checkInputs(nodeNumber, adjacencyEdges, stubEdges, chainEdges);

    if(nodeNumber == 0) { //Some of the following functions assume there are at least 2 nodes.
        return stList_construct();
    }

    /*
     * First calculate the optimal matching.
     */
    stList *chosenEdges = matchingAlgorithm(adjacencyEdges, nodeNumber);
    makeMatchingPerfect(chosenEdges, adjacencyEdges, nodeNumber);

    /*
     * Merge in the stub free components.
     */
    stList *updatedChosenEdges = mergeSimpleCycles2(chosenEdges,
            adjacencyEdges, stubEdges, chainEdges);
    stList_destruct(chosenEdges);
    chosenEdges = updatedChosenEdges;

    /*
     * Split stub components.
     */
    if (makeStubCyclesDisjoint) {
        updatedChosenEdges = splitMultipleStubCycles(chosenEdges,
                adjacencyEdges, stubEdges, chainEdges);
        stList_destruct(chosenEdges);
        chosenEdges = updatedChosenEdges;
    }

    return chosenEdges;
}
