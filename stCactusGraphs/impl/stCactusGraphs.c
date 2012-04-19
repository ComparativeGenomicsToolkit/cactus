/*
 * stCactusGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include "sonLib.h"
#include "stCactusGraphs.h"
#include "3_Absorb3edge2x.h"

struct _stCactusNode {
    stCactusEdgeEnd *head;
    stCactusEdgeEnd *tail;
    void *nodeObject;
};

struct _stCactusEdgeEnd {
    stCactusEdgeEnd *otherEdgeEnd;
    stCactusNode *node;
    stCactusEdgeEnd *nEdgeEnd;
    void *endObject;
    stCactusEdgeEnd *link;
    bool linkOrientation;
    bool isChainEnd;
};

struct _stCactusGraph {
    stHash *objectToNodeHash;
};

//Node functions

stCactusNode *stCactusNode_construct(stCactusGraph *graph, void *nodeObject) {
    stCactusNode *node = st_calloc(1, sizeof(stCactusNode));
    node->nodeObject = nodeObject;
    stHash_insert(graph->objectToNodeHash, nodeObject, node);
    return node;
}

void *stCactusNode_getObject(stCactusNode *node) {
    return node->nodeObject;
}

stCactusNode_edgeEndIt stCactusNode_getEdgeEndIt(stCactusNode *node) {
    stCactusNode_edgeEndIt it;
    it.edgeEnd = node->head;
    return it;
}

stCactusEdgeEnd *stCactusNode_edgeEndIt_getNext(stCactusNode_edgeEndIt *it) {
    stCactusEdgeEnd *edgeEnd = it->edgeEnd;
    if (edgeEnd != NULL) {
        it->edgeEnd = it->edgeEnd->nEdgeEnd;
    }
    return edgeEnd;
}

//Private node functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edge);

static void stCactusNode_destruct(stCactusNode *node) {
    stCactusNode_edgeEndIt it = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&it)) != NULL) {
        stCactusEdgeEnd_destruct(edgeEnd);
    }
    free(node);
}

static void stCactusNode_mergeNodes(stCactusGraph *graph, stCactusNode *node1,
        stCactusNode *node2, void *(*mergeNodeObjects)(void *, void *)) {
    assert(
            stHash_search(graph->objectToNodeHash,
                    stCactusNode_getObject(node2)) == node2);
    stHash_remove(graph->objectToNodeHash, stCactusNode_getObject(node2));
    if (node1->head == NULL) {
        node1->head = node2->head;
    } else {
        assert(node1->tail->nEdgeEnd == NULL);
        node1->tail->nEdgeEnd = node1->head;
    }
    node1->tail = node2->tail;
    node1->nodeObject = mergeNodeObjects(node1->nodeObject, node2->nodeObject);
    free(node2);
}

//Edge functions

static void connectUpEdgeEnd(stCactusEdgeEnd *edgeEnd, stCactusNode *node,
        stCactusEdgeEnd *otherEdgeEnd, void *endObject) {
    edgeEnd->node = node;
    edgeEnd->otherEdgeEnd = otherEdgeEnd;
    edgeEnd->endObject = endObject;
    if (node->head == NULL) {
        node->head = edgeEnd;
    } else {
        node->tail->nEdgeEnd = edgeEnd;
    }
    node->tail = edgeEnd;
}

stCactusEdgeEnd *stCactusEdgeEnd_construct(stCactusGraph *graph,
        stCactusNode *node1, stCactusNode *node2, void *edgeEndObject1,
        void *edgeEndObject2) {
    stCactusEdgeEnd *edgeEnd1 = st_calloc(1, sizeof(stCactusEdgeEnd));
    stCactusEdgeEnd *edgeEnd2 = st_calloc(1, sizeof(stCactusEdgeEnd));

    connectUpEdgeEnd(edgeEnd1, node1, edgeEnd2, edgeEndObject1);
    connectUpEdgeEnd(edgeEnd2, node2, edgeEnd1, edgeEndObject2);
    return edgeEnd1;
}

void *stCactusEdgeEnd_getObject(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->endObject;
}

stCactusNode *stCactusEdgeEnd_getNode(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->node;
}

stCactusNode *stCactusEdgeEnd_getOtherNode(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->otherEdgeEnd->node;
}

stCactusEdgeEnd *stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->otherEdgeEnd;
}

stCactusEdgeEnd *stCactusEdgeEnd_getLink(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->link;
}

bool stCactusEdgeEnd_getLinkOrientation(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->linkOrientation;
}

bool stCactusEdgeEnd_isChainEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->isChainEnd;
}

//Private edge functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edge) {
    free(edge);
}

static void stCactusEdgeEnd_setLink(stCactusEdgeEnd *edgeEnd,
        stCactusEdgeEnd *otherEdgeEnd) {
    edgeEnd->link = otherEdgeEnd;
}

static void stCactusEdgeEnd_setLinkOrientation(stCactusEdgeEnd *edgeEnd,
        bool orientation) {
    edgeEnd->linkOrientation = orientation;
}

static void stCactusEdgeEnd_setIsChainEnd(stCactusEdgeEnd *edgeEnd,
        bool isChainEnd) {
    edgeEnd->isChainEnd = isChainEnd;
}

//Graph functions

stCactusGraph *stCactusGraph_construct(void *(*nodeIt)(),
        stCactusEdgeTuple(*edgeIt)(),
        void *(*mergeNodeObjects)(void *, void *), void *startNode) {
    stCactusGraph *cactusGraph = st_malloc(sizeof(stCactusGraph));
    cactusGraph->objectToNodeHash = stHash_construct2(NULL,
            (void(*)(void *)) stCactusNode_destruct);
    return cactusGraph;
}

void stCactusGraph_destruct(stCactusGraph *graph) {
    stHash_destruct(graph->objectToNodeHash);
    free(graph);
}

stCactusNode *stCactusGraph_getNode(stCactusGraph *node, void *nodeObject) {
    return stHash_search(node->objectToNodeHash, nodeObject);
}

void stCactusGraph_collapseToCactus(stCactusGraph *graph,
        void *(*mergeNodeObjects)(void *, void *), stCactusNode *startNode) {
    //Basic data structures
    stHash *nodesToPositions = stHash_construct();
    stHash *positionsToNodes = stHash_construct2(
            (void(*)(void *)) stIntTuple_destruct, NULL);
    stList *adjacencyList = stList_construct3(0,
            (void(*)(void *)) stList_destruct);

    //Set up and run the three edge connected code.
    stCactusGraphNodeIterator *nodeIt = stCactusGraphNodeIterator_construct(
            graph);
    stCactusNode *node;
    int32_t nodeIdCounter = 0;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stIntTuple *nodeId = stIntTuple_construct(1, nodeIdCounter++);
        stHash_insert(nodesToPositions, node, nodeId);
        stHash_insert(positionsToNodes, nodeId, node);
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    nodeIt = stCactusGraphNodeIterator_construct(graph);
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stList *edges = stList_construct();
        stList_append(adjacencyList, edges);
        stCactusNode_edgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stCactusNode *otherNode = stCactusEdgeEnd_getOtherNode(edgeEnd);
            stIntTuple *otherNodeId =
                    stHash_search(nodesToPositions, otherNode);
            assert(otherNodeId != NULL);
            stList_append(edges, otherNodeId);
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    //Now do the merging
    stList *_3EdgeConnectedComponents = computeThreeEdgeConnectedComponents(
            adjacencyList);
    for (int32_t i = 0; i < stList_length(_3EdgeConnectedComponents); i++) {
        stList *_3EdgeConnectedComponent = stList_get(
                _3EdgeConnectedComponents, i);
        stCactusNode *node = stList_get(_3EdgeConnectedComponent, 0);
        for (int32_t j = 1; j < stList_length(_3EdgeConnectedComponent); j++) {
            stCactusNode *otherNode = stList_get(_3EdgeConnectedComponent, j);
            if (otherNode == startNode) { //This prevents the start node from being destructed.
                otherNode = node;
                node = startNode;
            }
            stCactusNode_mergeNodes(graph, node, otherNode, mergeNodeObjects);
        }
    }
    //Cleanup
    stList_destruct(adjacencyList);
    stList_destruct(_3EdgeConnectedComponents);
    stHash_destruct(nodesToPositions);
    stHash_destruct(positionsToNodes);

    //Mark the cycles
    stCactusGraph_markCycles(graph, startNode);
}

stCactusGraphNodeIterator *stCactusGraphNodeIterator_construct(
        stCactusGraph *graph) {
    stCactusGraphNodeIterator *nodeIt = st_malloc(
            sizeof(stCactusGraphNodeIterator));
    nodeIt->it = stHash_getIterator(graph->objectToNodeHash);
    nodeIt->graph = graph;
    return nodeIt;
}

stCactusNode *stCactusGraphNodeIterator_getNext(
        stCactusGraphNodeIterator *nodeIt) {
    void *key = stHash_getNext(nodeIt->it);
    if (key != NULL) {
        return stHash_search(nodeIt->graph->objectToNodeHash, key);
    }
    return NULL;
}

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIterator *nodeIt) {
    stHash_destructIterator(nodeIt->it);
    free(nodeIt);
}

void stCactusGraph_unmarkCycles(stCactusGraph *graph) {
    stCactusGraphNodeIterator *nodeIt = stCactusGraphNodeIterator_construct(
            graph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNode_edgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stCactusEdgeEnd_setIsChainEnd(edgeEnd, 0);
            stCactusEdgeEnd_setLink(edgeEnd, NULL);
            stCactusEdgeEnd_setLinkOrientation(edgeEnd, 0);
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
}

void stCactusGraph_markCyclesP(stCactusEdgeEnd *edgeEnd,
        stHash *nodesOnChainToEdgeEnds, stList *chainPath) {
    stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
    stCactusEdgeEnd *edgeEnd2;
    if ((edgeEnd2 = stHash_search(nodesOnChainToEdgeEnds, node)) != NULL) { //We've traversed a cycle
        assert(edgeEnd != edgeEnd2);
        stCactusEdgeEnd_setLink(edgeEnd, edgeEnd2);
        stCactusEdgeEnd_setLink(edgeEnd2, edgeEnd);
        stCactusEdgeEnd_setIsChainEnd(edgeEnd, 1);
        stCactusEdgeEnd_setIsChainEnd(edgeEnd2, 1);
        stCactusEdgeEnd_setLinkOrientation(edgeEnd2, 1);
        //st_uglyf("I am making a chain %i:%i %i:%i\n",
        //        *((int32_t *)stCactusEdgeEnd_getObject(edgeEnd)),
        //        *((int32_t *)stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd))),
        //        *((int32_t *)stCactusEdgeEnd_getObject(edgeEnd2)),
        //        *((int32_t *)stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2))));
        for (int32_t j = stList_length(chainPath) - 1; j >= 0; j -= 2) {
            stCactusEdgeEnd *edgeEnd3 = stList_get(chainPath, j);
            stCactusEdgeEnd *edgeEnd4 = stList_get(chainPath, j - 1);
            assert(edgeEnd3 != edgeEnd);
            assert(edgeEnd4 != edgeEnd);
            assert(edgeEnd4 != edgeEnd2);
            if (edgeEnd3 == edgeEnd2) {
                break;
            }
            stCactusEdgeEnd_setLink(edgeEnd3, edgeEnd4);
            stCactusEdgeEnd_setLink(edgeEnd4, edgeEnd3);
            stCactusEdgeEnd_setLinkOrientation(edgeEnd3, 1);
        }
    } else {
        stList_append(chainPath, edgeEnd);
        stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
        while ((edgeEnd2 = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
            if (edgeEnd2 != edgeEnd && stCactusEdgeEnd_getLink(edgeEnd2)
                    == NULL) {
                stHash_insert(nodesOnChainToEdgeEnds, node, edgeEnd2);
                stList_append(chainPath, edgeEnd2);
                stCactusGraph_markCyclesP(
                        stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2),
                        nodesOnChainToEdgeEnds, chainPath);
                stHash_remove(nodesOnChainToEdgeEnds, node);
                stList_pop(chainPath);
            }
        }
        stList_pop(chainPath);
    }
}

void stCactusGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode) {
    //st_uglyf("Starting to mark the cycles\n");
    stHash *nodesOnChainToEdgeEnds = stHash_construct();
    stList *chainPath = stList_construct();
    stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stHash_insert(nodesOnChainToEdgeEnds, startNode, edgeEnd);
            stCactusGraph_markCyclesP(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd),
                    nodesOnChainToEdgeEnds, chainPath);
            stHash_remove(nodesOnChainToEdgeEnds, startNode);
        }
    }
    stHash_destruct(nodesOnChainToEdgeEnds);
    stList_destruct(chainPath);
    //st_uglyf("Finished marking the cycles\n");
}

void stCactusGraph_collapseBridgesP(stCactusNode *parentNode,
        stCactusEdgeEnd *edgeEnd, stList *nodesToMerge) {
    stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
    if (stCactusEdgeEnd_getLink(edgeEnd) == NULL) { //Is a bridge
        //Establish if this is a leaf
        stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd2;
        int32_t bridges = 0;
        while ((edgeEnd2 = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
            if (edgeEnd2 != edgeEnd) {
                if (stCactusEdgeEnd_getLink(edgeEnd2) == NULL) {
                    bridges++;
                }
                if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                    stCactusGraph_collapseBridgesP(parentNode,
                            stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2),
                            nodesToMerge);
                }
            }
        }
        if (bridges != 1) {
            stList_append(nodesToMerge, parentNode);
            stList_append(nodesToMerge, node);
        }
    } else {
        stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd2;
        while ((edgeEnd2 = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
            if (edgeEnd2 != edgeEnd && !stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                stCactusGraph_collapseBridgesP(node,
                        stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2), nodesToMerge);
            }
        }
    }
}

void stCactusGraph_collapseBridges(stCactusGraph *graph,
        stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *)) {
    stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    stList *bridgesToMerge = stList_construct();
    while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stCactusGraph_collapseBridgesP(startNode,
                    stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd), bridgesToMerge);
        }
    }
    //Now merge bridges
    for (int32_t i = 0; i < stList_length(bridgesToMerge); i += 2) {
        stCactusNode *parentNode = stList_get(bridgesToMerge, i);
        stCactusNode *nodeToMerge = stList_get(bridgesToMerge, i + 1);
        stCactusNode_mergeNodes(graph, parentNode, nodeToMerge,
                mergeNodeObjects);
    }
    stList_destruct(bridgesToMerge);
    stCactusGraph_unmarkCycles(graph);
    stCactusGraph_markCycles(graph, startNode);
}

