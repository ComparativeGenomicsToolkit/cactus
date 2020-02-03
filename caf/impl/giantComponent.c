/*
 * massiveComponent.c
 *
 *  Created on: 22 Feb 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "stPinchGraphs.h"
#include <math.h>
#include <stdlib.h>

static void *getValue(stHash *hash, int64_t node) {
    stIntTuple *nodeTuple = stIntTuple_construct1(node);
    void *object = stHash_search(hash, nodeTuple);
    stIntTuple_destruct(nodeTuple);
    return object;
}

stList *stCaf_breakupComponentGreedily(stList *nodes, stList *edges, int64_t maxComponentSize) {
    /*
     * Make a component for each node in the graph
     */
    stHash *nodeToComponents = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
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
    int64_t edgeScore = INT64_MAX;
    //While edges exist, try and put them into the graph.
    stList *edgesToDelete = stList_construct();
    int64_t totalComponents = stList_length(nodes);
    while (stList_length(sortedEdges) > 0) {
        stIntTuple *edge = stList_pop(sortedEdges);
        if (edgeScore < stIntTuple_get(edge, 0)) {
            st_errAbort("bad edgeScore");
        }
        edgeScore = stIntTuple_get(edge, 0);
        stSortedSet *component1 = getValue(nodeToComponents, stIntTuple_get(edge, 1));
        stSortedSet *component2 = getValue(nodeToComponents, stIntTuple_get(edge, 2));
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
            "We broke a graph with %" PRIi64 " nodes and %" PRIi64 " edges for a max component size of %" PRIi64 " into %" PRIi64 " distinct components with %" PRIi64 " edges, discarding %" PRIi64 " edges\n",
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
    //Make nodes
    *nodes = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    stHash *pinchEndsToNodesHash = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, NULL, NULL);
    for (int64_t i = 0; i < stList_length(adjacencyComponent); i++) {
        stIntTuple *node = stIntTuple_construct1(i);
        stList_append(*nodes, node);
        assert(stHash_search(pinchEndsToNodesHash, stList_get(adjacencyComponent, i)) == NULL);
        stHash_insert(pinchEndsToNodesHash, stList_get(adjacencyComponent, i), node);
    }
    //Make edges

    //First build a hash of edges to their multiplicity
    stHash *edgesToMultiplicityHash = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct,
            (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(adjacencyComponent); i++) {
        stPinchEnd *pinchEnd1 = stList_get(adjacencyComponent, i);
        int64_t node1 = stIntTuple_get(stHash_search(pinchEndsToNodesHash, pinchEnd1), 0);
        stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(stPinchEnd_getBlock(pinchEnd1));
        stPinchSegment *segment;
        while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
            bool traverse5Prime = stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(pinchEnd1), segment);
            stPinchSegment *segment2 = traverse5Prime ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
            while (segment2 != NULL) {
                if (stPinchSegment_getBlock(segment2) != NULL) {
                    stPinchEnd pinchEnd2 = stPinchEnd_constructStatic(stPinchSegment_getBlock(segment2),
                            stPinchEnd_endOrientation(traverse5Prime, segment2));
                    assert(stHash_search(pinchEndsToNodesHash, &pinchEnd2) != NULL);
                    int64_t node2 = stIntTuple_get(stHash_search(pinchEndsToNodesHash, &pinchEnd2), 0);
                    if (node1 != node2) { //Ignore self edges
                        stIntTuple *edge = node1 < node2 ? stIntTuple_construct2(node1, node2) : stIntTuple_construct2(node2, node1);
                        int64_t multiplicity = 1;
                        if (stHash_search(edgesToMultiplicityHash, edge) != NULL) {
                            stIntTuple *count = stHash_removeAndFreeKey(edgesToMultiplicityHash, edge);
                            multiplicity += stIntTuple_get(count, 0);
                            stIntTuple_destruct(count);
                            assert(multiplicity > 1);
                        }
                        stHash_insert(edgesToMultiplicityHash, edge, stIntTuple_construct1(multiplicity));
                    }
                    break;
                }
                segment2 = traverse5Prime ? stPinchSegment_get5Prime(segment2) : stPinchSegment_get3Prime(segment2);
            }
        }
    }
    //Now build edges, scoring them according to their multiplicity
    *edges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    stHashIterator *hashIt = stHash_getIterator(edgesToMultiplicityHash);
    stIntTuple *edge;
    while ((edge = stHash_getNext(hashIt)) != NULL) {
        stIntTuple *count = stHash_search(edgesToMultiplicityHash, edge);
        stList_append(*edges, stIntTuple_construct3(stIntTuple_get(count, 0), stIntTuple_get(edge, 0), stIntTuple_get(edge, 1)));
    }

    //Cleanup
    stHash_destruct(pinchEndsToNodesHash);
    stHash_destruct(edgesToMultiplicityHash);
}

static void breakEdges(stPinchThreadSet *threadSet, stPinchEnd *pinchEnd1, stPinchEnd *pinchEnd2) {
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(stPinchEnd_getBlock(pinchEnd1));
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        bool traverse5Prime = stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(pinchEnd1), segment);
        stPinchSegment *segment2 = traverse5Prime ? stPinchSegment_get5Prime(segment) : stPinchSegment_get3Prime(segment);
        int64_t start = stPinchSegment_getStart(segment) + (traverse5Prime ? 0 : stPinchSegment_getLength(segment) - 1);
        while (segment2 != NULL) {
            stPinchBlock *pinchBlock2 = stPinchSegment_getBlock(segment2);
            if (pinchBlock2 != NULL) {
                if (pinchBlock2 == stPinchEnd_getBlock(pinchEnd2) && stPinchEnd_endOrientation(traverse5Prime, segment2)
                        == stPinchEnd_getOrientation(pinchEnd2)) { //Have an edge
                    int64_t end = stPinchSegment_getStart(segment2) + (traverse5Prime ? stPinchSegment_getLength(segment2) - 1 : 0);
                    assert(end != start);
                    if (llabs(end - start) > 1) {
                        stPinchThread *thread = stPinchSegment_getThread(segment);
                        int64_t splitPoint = (end + start) / 2;
                        assert((splitPoint > start && splitPoint < end) || (splitPoint < start && splitPoint > end));
                        stPinchThread_split(thread, splitPoint);
                        stPinchThread_split(thread, splitPoint - 1);
                        stPinchSegment *segment3 = stPinchThread_getSegment(thread, splitPoint);
                        assert(stPinchSegment_getBlock(segment3) == NULL);
                        assert(stPinchSegment_getLength(segment3) == 1);
                        stPinchBlock_construct2(segment3);
                        st_logDebug("Split an edge in a giant component\n");
                    } else {
                        st_logInfo("Encountered an edge in a giant component which can not be broken due its short length\n");
                    }
                }
                break;
            }
            segment2 = traverse5Prime ? stPinchSegment_get5Prime(segment2) : stPinchSegment_get3Prime(segment2);
        }
    }
}

void stCaf_breakupComponentsGreedily(stPinchThreadSet *threadSet, float maximumAdjacencyComponentSizeRatio) {
    int64_t maximumAdjacencyComponentSize = maximumAdjacencyComponentSizeRatio * log(stPinchThreadSet_getTotalBlockNumber(threadSet) * 2);
    if (maximumAdjacencyComponentSize < 10) {
        maximumAdjacencyComponentSize = 10;
    }
    //Get adjacency components
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (maximumAdjacencyComponentSize < stList_length(adjacencyComponent)) {
            //Get graph description
            stList *nodes, *edges;
            convertToNodesAndEdges(adjacencyComponent, &nodes, &edges);
            //Get the edges to remove
            stList *edgesToDelete = stCaf_breakupComponentGreedily(nodes, edges, maximumAdjacencyComponentSize);
            //Break edges;
            int64_t unbrokenEdges = 0;
            for (int64_t j = 0; j < stList_length(edgesToDelete); j++) {
                stIntTuple *edge = stList_get(edgesToDelete, j);
                assert(stIntTuple_get(edge, 1) < stIntTuple_get(edge, 2));
                stPinchEnd *pinchEnd1 = stList_get(adjacencyComponent, stIntTuple_get(edge, 1));
                stPinchEnd *pinchEnd2 = stList_get(adjacencyComponent, stIntTuple_get(edge, 2));
                if (stPinchBlock_getDegree(stPinchEnd_getBlock(pinchEnd1)) > 1 && stPinchBlock_getDegree(stPinchEnd_getBlock(pinchEnd2))
                        > 1) {
                    breakEdges(threadSet, pinchEnd1, pinchEnd2);
                } else {
                    unbrokenEdges++;
                }
            }
            if (stList_length(edgesToDelete) > 0) {
                st_logInfo("Pinch graph component with %" PRIi64 " nodes and %" PRIi64 " edges is being split up by breaking %" PRIi64 " edges to reduce size to less than %" PRIi64 " max, but found %" PRIi64 " pointless edges \n",
                           stList_length(nodes), stList_length(edges), stList_length(edgesToDelete), maximumAdjacencyComponentSize, unbrokenEdges);
            }
            //Cleanup
            stList_destruct(edges);
            stList_destruct(nodes);
            stList_destruct(edgesToDelete);
        }
    }
    stList_destruct(adjacencyComponents);
}
