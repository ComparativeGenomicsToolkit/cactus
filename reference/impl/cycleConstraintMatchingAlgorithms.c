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

static stList *filter(stList *list, stSortedSet *set, bool include) {
    /*
     * Returns a new list, either containing the intersection with set if include is non-zero,
     * or containing the set difference if include is zero.
     */
    stList *list2 = stList_construct();
    for (int32_t i = 0; i < stList_length(list); i++) {
        void *o = stList_get(list, i);
        if ((stSortedSet_search(set, o) != NULL && include) || (stSortedSet_search(set, o) == NULL && !include)) {
            stList_append(list2, o);
        }
    }
    return list2;
}

static stList *filterToExclude(stList *list, stSortedSet *set) {
    /*
     * Returns a new list, identical to list, but with any elements containined in set removed.
     */
    return filter(list, set, 0);
}

static stList *filterToInclude(stList *list, stSortedSet *set) {
    return filter(list, set, 1);
}

static stList *filterLists(stList *listOfLists, stSortedSet *set) {
    /*
     * Takes a list of lists and returns a new list of lists whose elements are the product of applying
     * filter to each member of listOfLists in the same order.
     */
    stList *listOfLists2 = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(listOfLists); i++) {
        stList_append(listOfLists2, filterToExclude(stList_get(listOfLists, i), set));
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

static void getComponentsP(stHash *nodesToEdges, stIntTuple *edge,
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

static void getComponentsP2(stHash *nodesToEdges, int32_t node,
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
            getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 0),
                    component);
            getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 1),
                    component);
        }
        stList_destruct(edges);
    }
    stIntTuple_destruct(key);
}

static stList *getComponents(stList *edges) {
    /*
     * Gets a list of connected components, each connected component
     * being represented as a list of the edges, such that each edge is in exactly one
     * connected component.
     */
    stHash *nodesToEdges = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL,
            NULL);

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
    stList *components =
            stList_construct3(0, (void(*)(void *)) stList_destruct);
    stList *keys = stHash_getKeys(nodesToEdges);
    for (int32_t i = 0; i < stList_length(keys); i++) {
        stIntTuple *key = stList_get(keys, i);
        stList *edges = stHash_search(nodesToEdges, key);
        if (edges != NULL) { //We have a component to build
            stSortedSet
                    *component =
                            stSortedSet_construct3(
                                    (int(*)(const void *, const void *)) stIntTuple_cmpFn,
                                    NULL);
            stHash_remove(nodesToEdges, key);
            for (int32_t i = 0; i < stList_length(edges); i++) {
                stIntTuple *edge = stList_get(edges, i);
                getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 0),
                        component);
                getComponentsP2(nodesToEdges, stIntTuple_getPosition(edge, 1),
                        component);
            }
            stList_append(components, stSortedSet_getList(component));
            stList_destruct(edges);
            stSortedSet_destruct(component);
        }
    }
    assert(stHash_size(nodesToEdges) == 0);
    stHash_destruct(nodesToEdges);
    stList_setDestructor(keys, (void(*)(void *)) stIntTuple_destruct);
    stList_destruct(keys);

    return components;
}

static stList *getComponents2(stList *adjacencyEdges, stList *stubEdges,
        stList *blockEdges) {
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
    adjacencySwitch->newEdge1 = newEdge2;
    adjacencySwitch->cost = cost;
    return adjacencySwitch;
}

static void adjacencySwitch_destruct(AdjacencySwitch *adjacencySwitch) {
    free(adjacencySwitch);
}

static AdjacencySwitch *adjacencySwitch_update(AdjacencySwitch *adjacencySwitch, stIntTuple *oldEdge1,
        stIntTuple *oldEdge2, stIntTuple *newEdge1, stIntTuple *newEdge2,
        int32_t cost) {
    /*
     * Convenience function, replaces one adjacency edge with a new one if the cost is lower.
     */
    if(adjacencySwitch == NULL) {
        return adjacencySwitch_construct(oldEdge1,
                oldEdge2, newEdge1, newEdge2, cost);
    }
    if(adjacencySwitch->cost > cost) {
        adjacencySwitch_destruct(adjacencySwitch);
        return adjacencySwitch_construct(oldEdge1,
                        oldEdge2, newEdge1, newEdge2, cost);
    }
    return adjacencySwitch;
}

static stList *getValidEdges(int32_t node, stHash *nodesToWeightedEdges, stSortedSet *validEdges) {
    /*
     * Iterates through the edges for a node and returns those in the set of the valid edges.
     */
    stIntTuple *node1 = stIntTuple_construct(1, node);
    stList *edges = stHash_search(nodesToWeightedEdges, node1);
    assert(edges != NULL);
    return filterToInclude(edges, validEdges);
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

static AdjacencySwitch *getBest3Or4EdgeSwitchP(stIntTuple *oldEdge1, int32_t node1,
        stHash *nodesToWeightedEdges, stSortedSet *allCurrentEdgesSet,
        stSortedSet *bridgingSwitchEdges) {
    /*
     * Returns the best adjacency switch for the given node and edge that
     * contains at least three existing edges.
     */
    int32_t node4 = getOtherPosition(oldEdge1, node1);
    AdjacencySwitch *minimumCostAdjacencySwitch = NULL;
    stList *validEdges = getValidEdges(node1, nodesToWeightedEdges,
            bridgingSwitchEdges);
    for (int32_t i = 0; i < stList_length(validEdges); i++) {
        stIntTuple *newEdge1 = stList_get(validEdges, i);
        int32_t node2 = getOtherPosition(newEdge1, node1);
        stList *validEdges2 = getValidEdges(node2, nodesToWeightedEdges,
                allCurrentEdgesSet);
        for (int32_t j = 0; j < stList_length(validEdges2); j++) {
            /*
             * At this point we have a three edge switch.
             */
            stIntTuple *oldEdge2 = stList_get(validEdges2, j);
            int32_t _3EdgeCost = stIntTuple_getPosition(oldEdge1, 2) + + stIntTuple_getPosition(
                                oldEdge2, 2)
                                - stIntTuple_getPosition(newEdge1, 2);
                        minimumCostAdjacencySwitch = adjacencySwitch_update(minimumCostAdjacencySwitch,
                                oldEdge1, oldEdge2, newEdge1, NULL, _3EdgeCost);
            /*
             * Now search for a forth edge to complete the switch.
             */
            int32_t node3 = getOtherPosition(oldEdge2, node2);
            stList *validEdges3 = getValidEdges(node3, nodesToWeightedEdges,
                    bridgingSwitchEdges);
            for (int32_t k = 0; k < stList_length(validEdges3); k++) {
                stIntTuple *newEdge2 = stList_get(validEdges3, k);
                if (getOtherPosition(newEdge2, node3) == node4) { //Success we have a lovely 4 edge switch
                    int32_t _4EdgeCost = _3EdgeCost - stIntTuple_getPosition(
                            newEdge2, 2);
                    minimumCostAdjacencySwitch = adjacencySwitch_update(minimumCostAdjacencySwitch,
                                        oldEdge1, oldEdge2, newEdge1, newEdge2, _4EdgeCost);
                }
            }
            stList_destruct(validEdges3);
        }
        stList_destruct(validEdges2);
    }
    stList_destruct(validEdges);
    return minimumCostAdjacencySwitch;
}

static AdjacencySwitch *getBest3Or4EdgeSwitch2(stIntTuple *oldEdge1,
        stHash *nodesToWeightedEdges, stSortedSet *allCurrentEdgesSet,
        stSortedSet *bridgingSwitchEdges) {
    /*
     * Returns the best 3 or 4 edge switch (one including 3 or 4 edges) for the given existing edge, if they exist, or else NULL.
     */
    AdjacencySwitch *minimumCostAdjacencySwitch1 = getBest3Or4EdgeSwitchP(oldEdge1,
            stIntTuple_getPosition(oldEdge1, 0), nodesToWeightedEdges, allCurrentEdgesSet, bridgingSwitchEdges);
    AdjacencySwitch *minimumCostAdjacencySwitch2 = getBest3Or4EdgeSwitchP(oldEdge1,
                stIntTuple_getPosition(oldEdge1, 0), nodesToWeightedEdges, allCurrentEdgesSet, bridgingSwitchEdges);
    if(minimumCostAdjacencySwitch1->cost < minimumCostAdjacencySwitch2->cost) {
        free(minimumCostAdjacencySwitch2);
        return minimumCostAdjacencySwitch1;
    }
    free(minimumCostAdjacencySwitch1);
    return minimumCostAdjacencySwitch2;
}

static void getNodeSetOfEdgesP(stSortedSet *nodes, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    if (stSortedSet_search(nodes, i) == NULL) {
        stSortedSet_insert(nodes, i);
    } else {
        stIntTuple_destruct(i);
    }
}

static stSortedSet *getNodeSetOfEdges(stList *edges) {
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

static AdjacencySwitch *getBest3Or4EdgeSwitch(stList *components,
        stHash *nodesToWeightedEdges) {
    /*
     * Returns the best 3 or 4 edge switch (one including 3 or 4 edges) for the given existing edge, if they exist, or else NULL.
     */
    stList *allComponentEdges = joinLists(components);
    stSortedSet *allComponentEdgesSet = stList_getSortedSet(allComponentEdges,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    assert(
            stList_length(allComponentEdges) == stSortedSet_size(
                    allComponentEdgesSet));

    /*
     * Get list of adjacency edges that bridge two components.
     */
    stSortedSet *bridgingSwitchEdges = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stSortedSet *componentNodes = getNodeSetOfEdges(component);
        stSortedSetIterator *it = stSortedSet_getIterator(componentNodes);
        stIntTuple *node;
        while ((node = stSortedSet_getNext(it)) != NULL) {
            stList *edges = stHash_search(nodesToWeightedEdges, node);
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
                    stSortedSet_insert(bridgingSwitchEdges, edge);
                }
                stIntTuple_destruct(node1);
                stIntTuple_destruct(node2);
            }
        }
        stSortedSet_destructIterator(it);
        stSortedSet_destruct(componentNodes);
    }

    /*
     * Look for the best 3 or 4 edge switch.
     */
    AdjacencySwitch *minimumAdjacencySwitch = NULL;
    for (int32_t i = 0; i < stList_length(allComponentEdges); i++) {
        AdjacencySwitch *adjacencySwitch = getBest3Or4EdgeSwitch2(
                stList_get(allComponentEdges, i), nodesToWeightedEdges,
                allComponentEdgesSet, bridgingSwitchEdges);
        if (minimumAdjacencySwitch == NULL || adjacencySwitch->cost
                < minimumAdjacencySwitch->cost) {
            free(minimumAdjacencySwitch);
            minimumAdjacencySwitch = adjacencySwitch;
        }
    }

    stList_destruct(allComponentEdges);
    stSortedSet_destruct(allComponentEdgesSet);

    return minimumAdjacencySwitch;
}

static AdjacencySwitch *getBest2EdgeSwitch(stList *components,
        stHash *nodesToWeightedEdges) {
    /*
     * Look for the two lowest value adjacency edges in all current edges that are in a separate component and returns them as an adjacency switch
     * with now new adjacency edges.
     */
    stList *lowestScoringEdgeFromEachComponent = stList_construct();
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        assert(stList_length(component) > 0);
        stIntTuple *lowestScoringEdge = stList_get(component, 0);
        int32_t lowestScore = stIntTuple_getPosition(lowestScoringEdge, 2);
        for (int32_t j = 1; j < stList_length(component); j++) {
            stIntTuple *edge = stList_get(component, j);
            int32_t k = stIntTuple_getPosition(edge, 2);
            if (k < lowestScore) {
                lowestScore = k;
                lowestScoringEdge = edge;
            }
        }
        stList_append(lowestScoringEdgeFromEachComponent, lowestScoringEdge);
    }
    stIntTuple *lowestScoreEdge1 = stList_get(
            lowestScoringEdgeFromEachComponent, 0);
    int32_t edgeScore1 = stIntTuple_getPosition(lowestScoreEdge1, 2);
    stIntTuple *lowestScoreEdge2 = stList_get(
            lowestScoringEdgeFromEachComponent, 1);
    int32_t edgeScore2 = stIntTuple_getPosition(lowestScoreEdge2, 2);
    if (edgeScore1 > edgeScore2) {
        void *o = lowestScoreEdge1;
        int32_t i = edgeScore1;
        lowestScoreEdge1 = lowestScoreEdge2;
        edgeScore1 = edgeScore2;
        lowestScoreEdge2 = o;
        edgeScore2 = i;
    }
    for (int32_t i = 2; i < stList_length(lowestScoringEdgeFromEachComponent); i++) {
        stIntTuple *edge = stList_get(lowestScoringEdgeFromEachComponent, i);
        int32_t edgeScore = stIntTuple_getPosition(edge, 2);
        assert(edgeScore1 <= edgeScore2);
        if (edgeScore < edgeScore1) {
            lowestScoreEdge2 = lowestScoreEdge1;
            edgeScore2 = edgeScore1;
            lowestScoreEdge1 = edge;
            edgeScore1 = edgeScore;
        } else if (edgeScore < edgeScore2) {
            lowestScoreEdge2 = edge;
            edgeScore2 = edgeScore;
        }
    }
    return adjacencySwitch_construct(lowestScoreEdge1, lowestScoreEdge2, NULL,
            NULL, edgeScore1 + edgeScore2);
}

static AdjacencySwitch *getBestAdjacencySwitch(stList *components,
        stHash *nodesToWeightedEdges) {
    /*
     * Modifies in place the list of components, merging together two components whose
     * switch has lowest reduction in score
     * and destroying the two unmerged components.
     */
    assert(stList_length(components) > 1);

    /*
     * Look for the lowest cost switch with no supporting adjacencies.
     */
    AdjacencySwitch *minimumAdjacencySwitch = getBest2EdgeSwitch(
                components, nodesToWeightedEdges);

    /*
     * Look for the lowest cost switch with at least one supporting adjacency.
     */
    AdjacencySwitch *minimumAdjacencySwitch2 = getBest3Or4EdgeSwitch(components,
            nodesToWeightedEdges);

    if (minimumAdjacencySwitch2 != NULL) {
        if (minimumAdjacencySwitch2->cost <= minimumAdjacencySwitch->cost) {
            free(minimumAdjacencySwitch);
            return minimumAdjacencySwitch2;
        }
        free(minimumAdjacencySwitch2);
        return minimumAdjacencySwitch;
    }
    return minimumAdjacencySwitch;
}

static void doBestMergeOfTwoCycles(stList *components, stHash *nodesToWeightedEdges) {
    /*
     * Merge the best two cycles.
     */

    AdjacencySwitch *adjacencySwitch = getBestAdjacencySwitch(components, nodesToWeightedEdges);
    //Find the two components..

    //Merge the components..


}

stList *mergeCycles(stList *components, stList *weightedEdges) {
    /*
     * Returns a single component, as a list of edges, by doing length(components)-1
     * calls to doBestMergeOfTwoCycles.
     */
    components = stList_copy(components, NULL);
    for (int32_t i = 0; i < stList_length(components); i++) { //Clone the complete list
        stList_set(components, i, stList_copy(stList_get(components, i), NULL));
    }
    while (stList_length(components) > 1) {
        doBestMergeOfTwoCycles(components, weightedEdges);
    }
    assert(stList_length(components) == 1);
    stList *mergedComponent = stList_get(components, 0);
    stList_destruct(components);
    return mergedComponent;
}

static stList *getStubAndBlockEdgeFreeComponents(stList *components,
        stList *stubEdges, stList *blockEdges) {
    /*
     * Returns a new list of components, such that each component has its stub and block edges removed.
     */
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges,
            blockEdges);
    stList *adjacencyOnlyComponents = filterLists(components,
            stubAndBlockEdgesSet);
    stSortedSet_destruct(stubAndBlockEdgesSet);

    return adjacencyOnlyComponents;
}

static stList *mergeComponents2(stList *chosenEdges, stList *adjacencyEdges,
        stList *stubEdges, stList *blockEdges) {
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
     * Remove the stub/block edges from the components.
     */
    stList_append(stubFreeComponents, globalComponent);
    stList *adjacencyOnlyComponents = getStubAndBlockEdgeFreeComponents(
            stubFreeComponents, stubEdges, blockEdges);

    stList_destruct(stubFreeComponents);
    stList_destruct(globalComponent);
    stList_destruct(components); //We only clean this up now, as this frees the components it contains.

    /*
     * Merge stub free components into the others.
     */
    stList *mergedComponents = mergeCycles(adjacencyOnlyComponents,
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
    /*
     * Get sub-components containing only adjacency and block edges.
     */
    //Filter out the the non-adjacency edges
    stList *adjacencyOnlyComponent = filterToExclude(component, stubAndBlockEdgesSet);
    //Make it only the block edges present in the original component
    stList *stubFreeComponents = getComponents2(adjacencyOnlyComponent, NULL,
            blockEdges);
    stList_destruct(adjacencyOnlyComponent);
    assert(stList_length(stubFreeComponents) >= 1);

    stList *splitComponents = stList_construct(); //The list to return.

    if (stList_length(stubFreeComponents) > 1) {
        /*
         * Merge together the best two components.
         */
        stList *l = filterLists(stubFreeComponents, stubAndBlockEdgesSet);
        doBestMergeOfTwoCycles(l, adjacencyEdges); //This is inplace.
        stList *l2 = joinLists(l);
        stList_destruct(l);
        l = getComponents2(l2, stubEdges, blockEdges);
        stList_destruct(l2);

        /*
         * Call procedure recursively.
         */
        for (int32_t i = 0; i < stList_length(l); i++) {
            l2 = splitMultipleStubComponent(stList_get(l, i), adjacencyEdges,
                    stubEdges, blockEdges, stubAndBlockEdgesSet);
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

static stList *splitMultipleStubCycles(stList *chosenEdges,
        stList *adjacencyEdges, stList *stubEdges, stList *blockEdges) {
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
    stSortedSet *stubAndBlockEdgesSet = getStubAndBlockEdgeSortedSet(stubEdges,
            blockEdges);
    stList *singleStubEdgeComponents = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    for (int32_t i = 0; i < stList_length(components); i++) {
        stList *component = stList_get(components, i);
        stList *splitComponents = splitMultipleStubComponent(component,
                adjacencyEdges, blockEdges, stubEdges, stubAndBlockEdgesSet);
        stList_appendAll(singleStubEdgeComponents, splitComponents);
        stList_setDestructor(splitComponents, NULL); //Do this to avoid destroying the underlying lists
        stList_destruct(splitComponents);
    }
    stList_destruct(components);

    /*
     * Remove the stub/block edges from the components.
     */
    stList *adjacencyOnlyComponents = filterLists(singleStubEdgeComponents,
            stubAndBlockEdgesSet);
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
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        assert(stIntTuple_getPosition(edge, 0) >= 0);
        assert(stIntTuple_getPosition(edge, 0) < nodeNumber);
        assert(stIntTuple_getPosition(edge, 1) >= 0);
        assert(stIntTuple_getPosition(edge, 1) < nodeNumber);
    }
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
    stList *updatedChosenEdges = mergeComponents2(chosenEdges, adjacencyEdges,
            stubEdges, blockEdges);
    stList_destruct(chosenEdges);
    chosenEdges = updatedChosenEdges;

    /*
     * Split stub components.
     */
    if (makeStubCyclesDisjoint) {
        updatedChosenEdges = splitMultipleStubCycles(chosenEdges,
                adjacencyEdges, stubEdges, blockEdges);
        stList_destruct(chosenEdges);
        chosenEdges = updatedChosenEdges;
    }

    return chosenEdges;
}
