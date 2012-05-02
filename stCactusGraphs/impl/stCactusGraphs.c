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
    void (*destructNodeObjectFn)(void *);
    void (*destructEdgeEndObjectFn)(void *);
};

//Node functions

stCactusNode *stCactusNode_construct(stCactusGraph *graph, void *nodeObject) {
    stCactusNode *node = st_calloc(1, sizeof(stCactusNode));
    node->nodeObject = nodeObject;
    assert(stHash_search(graph->objectToNodeHash, nodeObject) == NULL);
    stHash_insert(graph->objectToNodeHash, nodeObject, node);
    return node;
}

void *stCactusNode_getObject(stCactusNode *node) {
    return node->nodeObject;
}

stCactusNodeEdgeEndIt stCactusNode_getEdgeEndIt(stCactusNode *node) {
    stCactusNodeEdgeEndIt it;
    it.edgeEnd = node->head;
    return it;
}

stCactusEdgeEnd *stCactusNodeEdgeEndIt_getNext(stCactusNodeEdgeEndIt *it) {
    stCactusEdgeEnd *edgeEnd = it->edgeEnd;
    if (edgeEnd != NULL) {
        it->edgeEnd = it->edgeEnd->nEdgeEnd;
    }
    return edgeEnd;
}

stCactusEdgeEnd *stCactusNode_getFirstEdgeEnd(stCactusNode *node) {
    return node->head;
}

//Private node functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edge, void(*destructEdgeEndObjectFn)(void *));

static void stCactusNode_destruct(stCactusNode *node, void(*destructNodeObjectFn)(void *),
        void(*destructEdgeEndObjectFn)(void *)) {
    stCactusNodeEdgeEndIt it = stCactusNode_getEdgeEndIt(node);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&it)) != NULL) {
        stCactusEdgeEnd_destruct(edgeEnd, destructEdgeEndObjectFn);
    }
    if (destructNodeObjectFn != NULL) {
        destructNodeObjectFn(node->nodeObject);
    }
    free(node);
}

static void stCactusNode_mergeNodes(stCactusGraph *graph, stCactusNode *node1, stCactusNode *node2,
        void *(*mergeNodeObjects)(void *, void *)) {
    if (node1 == node2) {
        return;
    }
    assert(stHash_search(graph->objectToNodeHash, stCactusNode_getObject(node1)) == node1);
    assert(stHash_search(graph->objectToNodeHash, stCactusNode_getObject(node2)) == node2);
    stHash_remove(graph->objectToNodeHash, stCactusNode_getObject(node1));
    stHash_remove(graph->objectToNodeHash, stCactusNode_getObject(node2));
    node1->nodeObject = mergeNodeObjects(node1->nodeObject, node2->nodeObject);
    stHash_insert(graph->objectToNodeHash, stCactusNode_getObject(node1), node1);
    if (node2->head != NULL) {
        if (node1->head == NULL) {
            node1->head = node2->head;
        } else {
            assert(node1->tail->nEdgeEnd == NULL);
            node1->tail->nEdgeEnd = node2->head;
        }
        stCactusEdgeEnd *edgeEnd = node2->head;
        while (edgeEnd != NULL) {
            edgeEnd->node = node1;
            edgeEnd = edgeEnd->nEdgeEnd;
        }
        node1->tail = node2->tail;
    }
    free(node2);
}

//Edge functions

static void connectUpEdgeEnd(stCactusEdgeEnd *edgeEnd, stCactusNode *node, stCactusEdgeEnd *otherEdgeEnd,
        void *endObject) {
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

stCactusEdgeEnd *stCactusEdgeEnd_construct(stCactusGraph *graph, stCactusNode *node1, stCactusNode *node2,
        void *edgeEndObject1, void *edgeEndObject2) {
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

stCactusEdgeEnd *stCactusEdgeEnd_getNextEdgeEnd(stCactusEdgeEnd *edgeEnd) {
    return edgeEnd->nEdgeEnd;
}

//Private edge functions

static void stCactusEdgeEnd_destruct(stCactusEdgeEnd *edgeEnd, void(*destructEdgeEndObjectFn)(void *)) {
    if (destructEdgeEndObjectFn != NULL) {
        destructEdgeEndObjectFn(edgeEnd->endObject);
    }
    free(edgeEnd);
}

static void stCactusEdgeEnd_setLink(stCactusEdgeEnd *edgeEnd, stCactusEdgeEnd *otherEdgeEnd) {
    edgeEnd->link = otherEdgeEnd;
}

static void stCactusEdgeEnd_setLinkOrientation(stCactusEdgeEnd *edgeEnd, bool orientation) {
    edgeEnd->linkOrientation = orientation;
}

static void stCactusEdgeEnd_setIsChainEnd(stCactusEdgeEnd *edgeEnd, bool isChainEnd) {
    edgeEnd->isChainEnd = isChainEnd;
}

//Graph functions

stCactusGraph *stCactusGraph_construct2(void(*destructNodeObjectFn)(void *), void(*destructEdgeEndObjectFn)(void *)) {
    stCactusGraph *cactusGraph = st_malloc(sizeof(stCactusGraph));
    cactusGraph->objectToNodeHash = stHash_construct();
    cactusGraph->destructNodeObjectFn = destructNodeObjectFn;
    cactusGraph->destructEdgeEndObjectFn = destructEdgeEndObjectFn;
    return cactusGraph;
}

stCactusGraph *stCactusGraph_construct() {
    return stCactusGraph_construct2(NULL, NULL);
}

void stCactusGraph_destruct(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(graph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNode_destruct(node, graph->destructNodeObjectFn, graph->destructEdgeEndObjectFn);
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    stHash_destruct(graph->objectToNodeHash);
    free(graph);
}

stCactusNode *stCactusGraph_getNode(stCactusGraph *node, void *nodeObject) {
    return stHash_search(node->objectToNodeHash, nodeObject);
}

void stCactusGraph_collapseToCactus(stCactusGraph *graph, void *(*mergeNodeObjects)(void *, void *),
        stCactusNode *startNode) {
    //Basic data structures
    stHash *nodesToPositions = stHash_construct();
    stHash *positionsToNodes = stHash_construct3((uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct, NULL);
    stList *adjacencyList = stList_construct3(0, (void(*)(void *)) stList_destruct);

    //Set up and run the three edge connected code.
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(graph);
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
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stCactusNode *otherNode = stCactusEdgeEnd_getOtherNode(edgeEnd);
            stIntTuple *otherNodeId = stHash_search(nodesToPositions, otherNode);
            assert(otherNodeId != NULL);
            stList_append(edges, otherNodeId);
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    //Now do the merging
    stList *_3EdgeConnectedComponents = computeThreeEdgeConnectedComponents(adjacencyList);
    for (int32_t i = 0; i < stList_length(_3EdgeConnectedComponents); i++) {
        stList *_3EdgeConnectedComponent = stList_get(_3EdgeConnectedComponents, i);
        stCactusNode *node = stHash_search(positionsToNodes, stList_get(_3EdgeConnectedComponent, 0));
        assert(node != NULL);
        for (int32_t j = 1; j < stList_length(_3EdgeConnectedComponent); j++) {
            stCactusNode *otherNode = stHash_search(positionsToNodes, stList_get(_3EdgeConnectedComponent, j));
            assert(otherNode != NULL);
            assert(node != otherNode);
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

stCactusGraphNodeIt *stCactusGraphNodeIterator_construct(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = st_malloc(sizeof(stCactusGraphNodeIt));
    nodeIt->it = stHash_getIterator(graph->objectToNodeHash);
    nodeIt->graph = graph;
    return nodeIt;
}

stCactusNode *stCactusGraphNodeIterator_getNext(stCactusGraphNodeIt *nodeIt) {
    void *key = stHash_getNext(nodeIt->it);
    if (key != NULL) {
        return stHash_search(nodeIt->graph->objectToNodeHash, key);
    }
    return NULL;
}

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIt *nodeIt) {
    stHash_destructIterator(nodeIt->it);
    free(nodeIt);
}

void stCactusGraph_unmarkCycles(stCactusGraph *graph) {
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(graph);
    stCactusNode *node;
    while ((node = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(node);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stCactusEdgeEnd_setIsChainEnd(edgeEnd, 0);
            stCactusEdgeEnd_setLink(edgeEnd, NULL);
            stCactusEdgeEnd_setLinkOrientation(edgeEnd, 0);
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
}

static void makeChain(stCactusEdgeEnd *edgeEnd, stCactusEdgeEnd *edgeEnd2, stList *chainPath) {
    assert(edgeEnd != edgeEnd2);
    stCactusEdgeEnd_setLink(edgeEnd, edgeEnd2);
    stCactusEdgeEnd_setLink(edgeEnd2, edgeEnd);
    stCactusEdgeEnd_setIsChainEnd(edgeEnd, 1);
    stCactusEdgeEnd_setIsChainEnd(edgeEnd2, 1);
    stCactusEdgeEnd_setLinkOrientation(edgeEnd2, 1);
    for (int32_t j = stList_length(chainPath) - 1; j >= 0; j -= 2) {
        stCactusEdgeEnd *edgeEnd3 = stList_get(chainPath, j);
        stCactusEdgeEnd *edgeEnd4 = stList_get(chainPath, j - 1);
        assert(edgeEnd3 != edgeEnd2);
        assert(edgeEnd4 != edgeEnd2);
        assert(edgeEnd4 != edgeEnd);
        if (edgeEnd3 == edgeEnd) {
            break;
        }
        stCactusEdgeEnd_setLink(edgeEnd3, edgeEnd4);
        stCactusEdgeEnd_setLink(edgeEnd4, edgeEnd3);
        stCactusEdgeEnd_setLinkOrientation(edgeEnd4, 1);
    }
}

void stCactusGraph_markCyclesP(stCactusEdgeEnd *edgeEnd, stHash *nodesOnChainToEdgeEnds, stList *chainPath) {
    stList_append(chainPath, edgeEnd);
    stList_append(chainPath, NULL);
    while (stList_length(chainPath) > 0) {
        stCactusEdgeEnd *edgeEnd2 = stList_pop(chainPath);
        edgeEnd = stList_peek(chainPath);
        stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
        if (edgeEnd2 == NULL) {
            edgeEnd2 = stCactusNode_getFirstEdgeEnd(node);
        } else {
            edgeEnd2 = stCactusEdgeEnd_getNextEdgeEnd(edgeEnd2);
        }
        if (edgeEnd2 == NULL) {
            stList_pop(chainPath);
            stHash_remove(nodesOnChainToEdgeEnds, node);
        } else {
            stList_append(chainPath, edgeEnd2);
            if (edgeEnd2 != edgeEnd && stCactusEdgeEnd_getLink(edgeEnd2) == NULL) {
                stHash_insert(nodesOnChainToEdgeEnds, node, edgeEnd2);
                stCactusEdgeEnd *edgeEnd3 = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2);
                stCactusEdgeEnd *edgeEnd4;
                if ((edgeEnd4 = stHash_search(nodesOnChainToEdgeEnds, stCactusEdgeEnd_getNode(edgeEnd3))) != NULL) { //We've traversed a cycle
                    makeChain(edgeEnd4, edgeEnd3, chainPath);
                } else {
                    stList_append(chainPath, edgeEnd3);
                    stList_append(chainPath, NULL);
                }
            }
        }
    }
}

void stCactusGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode) {
    stHash *nodesOnChainToEdgeEnds = stHash_construct();
    stList *chainPath = stList_construct();
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stCactusEdgeEnd *edgeEnd2 = stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd);
            if (stCactusEdgeEnd_getNode(edgeEnd2) == startNode) {
                makeChain(edgeEnd, edgeEnd2, chainPath);
            } else {
                stHash_insert(nodesOnChainToEdgeEnds, startNode, edgeEnd);
                stCactusGraph_markCyclesP(edgeEnd2, nodesOnChainToEdgeEnds, chainPath);
                stHash_remove(nodesOnChainToEdgeEnds, startNode);
            }
        }
    }
    stList_destruct(chainPath);
    assert(stHash_size(nodesOnChainToEdgeEnds) == 0);
    stHash_destruct(nodesOnChainToEdgeEnds);
}

void stCactusGraph_collapseBridgesP2(stList *stack, stCactusNode *parentNode, stCactusEdgeEnd *edgeEnd,
        stCactusEdgeEnd *edgeEnd2) {
    stList_append(stack, parentNode);
    stList_append(stack, edgeEnd);
    stList_append(stack, edgeEnd2);
}

void stCactusGraph_collapseBridgesP(stCactusNode *parentNode, stCactusEdgeEnd *edgeEnd, stList *nodesToMerge) {
    stList *stack = stList_construct();
    stList_append(stack, parentNode);
    stList_append(stack, edgeEnd);
    while (stList_length(stack) > 0) {
        edgeEnd = stList_pop(stack);
        parentNode = stList_pop(stack);
        stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
        if (stCactusEdgeEnd_getLink(edgeEnd) == NULL) { //Is a bridge
            //Establish if this is a leaf
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd2;
            int32_t bridges = 0;
            while ((edgeEnd2 = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                if (edgeEnd2 != edgeEnd) {
                    if (stCactusEdgeEnd_getLink(edgeEnd2) == NULL) {
                        bridges++;
                    }
                    if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                        stList_append(stack, parentNode);
                        stList_append(stack, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                    }
                }
            }
            if (bridges != 1) {
                stList_append(nodesToMerge, parentNode);
                stList_append(nodesToMerge, node);
            }
        } else if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
            stCactusEdgeEnd *edgeEnd2;
            while ((edgeEnd2 = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
                if (edgeEnd2 != edgeEnd && !stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                    stList_append(stack, node);
                    stList_append(stack, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2));
                }
            }
        }
    }
    stList_destruct(stack);
}

void stCactusGraph_collapseBridges(stCactusGraph *graph, stCactusNode *startNode,
        void *(*mergeNodeObjects)(void *, void *)) {
    stCactusNodeEdgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    stList *bridgesToMerge = stList_construct();
    while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stCactusGraph_collapseBridgesP(startNode, stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd), bridgesToMerge);
        }
    }
    //Now merge bridges
    for (int32_t i = 0; i < stList_length(bridgesToMerge); i += 2) {
        stCactusNode *parentNode = stList_get(bridgesToMerge, i);
        stCactusNode *nodeToMerge = stList_get(bridgesToMerge, i + 1);
        stCactusNode_mergeNodes(graph, parentNode, nodeToMerge, mergeNodeObjects);
    }
    stList_destruct(bridgesToMerge);
    stCactusGraph_unmarkCycles(graph);
    stCactusGraph_markCycles(graph, startNode);
}

