/*
 * massiveComponent.c
 *
 *  Created on: 22 Feb 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "stPinchGraphs.h"

static void *getValue(stHash *hash, int32_t node) {
    stIntTuple *nodeTuple = stIntTuple_construct(1, node);
    void *object = stHash_search(hash, nodeTuple);
    stIntTuple_destruct(nodeTuple);
    return object;
}

stList *stCaf_breakupComponentGreedily(stList *nodes, stList *edges, int32_t maxComponentSize) {
    /*
     * Make a component for each node in the graph
     */
    stHash *nodeToComponents = stHash_construct3((uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
    stListIterator *listIt = stList_getIterator(nodes);
    stIntTuple *node;
    while ((node = stList_getNext(listIt)) != NULL) {
        stSortedSet *component = stSortedSet_construct();
        stSortedSet_insert(component, node);
        assert(stHash_search(nodeToComponents, node) == NULL);
        stHash_insert(nodeToComponents, node, component);
    }
    stList_destructIterator(listIt);

    stList *sortedEdges = stList_copy(edges, NULL); //copy, to avoid messing input
    stList_sort(sortedEdges, (int(*)(const void *, const void *)) stIntTuple_cmpFn); //Sort in ascending order, so best edge first
    int32_t edgeScore = INT32_MAX;
    //While edges exist, try and put them into the graph.
    stList *edgesToDelete = stList_construct();
    int32_t totalComponents = stList_length(nodes);
    while (stList_length(sortedEdges) > 0) {
        stIntTuple *edge = stList_pop(sortedEdges);
        assert(edgeScore >= stIntTuple_getPosition(edge, 0));
        edgeScore = stIntTuple_getPosition(edge, 0);
        stSortedSet *component1 = getValue(nodeToComponents, stIntTuple_getPosition(edge, 1));
        stSortedSet *component2 = getValue(nodeToComponents, stIntTuple_getPosition(edge, 2));
        assert(component1 != NULL && component2 != NULL);
        if (component1 == component2) { //We're golden, as the edge is already contained within one component.
            continue;
        }
        if (stSortedSet_size(component1) + stSortedSet_size(component2) > maxComponentSize) { //This edge would make a too large component, so reject
            stList_append(edgesToDelete, edge);
            continue;
        }
        //Merge the components and replace references in the hash.
        if (stSortedSet_size(component1) < stSortedSet_size(component2)) {
            stSortedSet *component3 = component1;
            component1 = component2;
            component2 = component3;
        }
        assert(stSortedSet_size(component1) >= stSortedSet_size(component2));
        while (stSortedSet_size(component2) > 0) {
            node = stSortedSet_getLast(component2);
            stSortedSet_remove(component2, node);
            assert(stSortedSet_search(component1, node) == NULL);
            stSortedSet_insert(component1, node);
            stHash_insert(nodeToComponents, node, component1);
        }
        stSortedSet_destruct(component2);
        totalComponents -= 1;
    }

    st_logDebug(
            "We broke a graph with %i nodes and %i edges for a max component size of %i into %i distinct components with %i edges, discarding %i edges\n",
            stList_length(nodes), stList_length(edges), maxComponentSize, totalComponents,
            stList_length(edges) - stList_length(edgesToDelete), stList_length(edgesToDelete));

    //Cleanup
    stList_destruct(sortedEdges);
    stList *components = stHash_getValues(nodeToComponents);
    stSortedSet *componentsSet = stList_getSortedSet(components, NULL);
    stList_destruct(components);
    stSortedSet_setDestructor(componentsSet, (void(*)(void *)) stSortedSet_destruct);
    stSortedSet_destruct(componentsSet);
    stHash_destruct(nodeToComponents);

    return edgesToDelete;
}

static void convertToNodesAndEdges(stList *adjacencyComponent, stList **nodes, stList **edges) {

}

static void breakEdges(stPinchThreadSet *threadSet, stPinchEnd *pinchEnd1, stPinchEnd *pinchEnd2) {
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(stPinchEnd_getBlock(pinchEnd1));
    stPinchSegment *segment;
    while((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        stPinchSegment *segment2;
        if(stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(pinchEnd1), segment)) {
            segment2 = stPinchSegment_get5Prime(segment);
        } else {
            segment2 = stPinchSegment_get5Prime(segment);
        }
        if(stPinchSegment_getBlock(segment2) == NULL) {

        }
    }
}

void stCaf_breakupComponentsGreedily(stPinchThreadSet *threadSet, float maximumAdjacencyComponentSizeRatio) {
    double maximumAdjacencyComponentSize = maximumAdjacencyComponentSizeRatio * log(stPinchThreadSet_getTotalBlockNumber(threadSet) * 2);
    //Get adjacency components
    stHash *pinchEndsToAdjacencyComponents;
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents2(threadSet, &pinchEndsToAdjacencyComponents);
    for (int32_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (maximumAdjacencyComponentSize < stList_length(adjacencyComponent)) {
            //Get graph description
            stList *nodes, *edges;
            convertToNodesAndEdges(adjacencyComponent, &nodes, &edges);
            //Get the edges to remove
            stList *edgesToDelete = stCaf_breakupComponentGreedily(nodes, edges, maximumAdjacencyComponentSize);
            if (stList_length(edgesToDelete) > 0) {
                printf("Pinch graph component with %i nodes and %i edges is being split up by breaking %i edges\n", stList_length(nodes),
                        stList_length(edges), stList_length(edgesToDelete));
            }
            //Break edges;
            for(int32_t j=0; j<stList_length(edgesToDelete); j++) {
                stIntTuple *edge = stList_get(edgesToDelete, j);
                stPinchEnd *pinchEnd1 = stList_get(adjacencyComponent, stIntTuple_getPosition(edge, 0));
                stPinchEnd *pinchEnd2 = stList_get(adjacencyComponent, stIntTuple_getPosition(edge, 1));
                breakEdges(threadSet, pinchEnd1, pinchEnd2);
            }
            //Cleanup
            stList_destruct(edges);
            stList_destruct(nodes);
            stList_destruct(edgesToDelete);
        }
    }
}
