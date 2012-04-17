/*
 * stGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include "sonLib.h"
#include "stCactusGraphs.h"

struct _stCactusNode {
    stCactusEdgeEnd *head;
    stCactusEdgeEnd *tail;
    void *nodeObject;
};

struct _stCactusNode_edgeEndIt {
    stCactusEdgeEnd *edgeEnd;
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

struct _stCactusGraphNodeIterator {
    stHashIterator *it;
    stCactusGraph *graph;
};

//Node functions

stCactusNode *stCactusNode_construct(stCactusGraph *graph,
        void *nodeObject) {
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
    assert(stHash_search(graph->objectToNodeHash, stCactusNode_getObject(node2)) == node2);
    stHash_remove(graph->objectToNodeHash, stCactusNode_getObject(node2));
    if (node1->head == NULL) {
        node1->head = node2->head;
        node2->tail = node2->tail;
    } else {
        assert(node1->tail->nEdgeEnd == NULL);
        node1->tail->nEdgeEnd = node1->head;
        node1->tail = node2->tail;
    }
    node1->nodeObject = mergeNodeObjects(node1->nodeObject, node2->nodeObject);
    free(node2);
}

//Edge functions

stCactusEdgeEnd *stCactusEdgeEnd_construct(stCactusGraph *graph,
        stCactusNode *node1, stCactusNode *node2, void *edgeEndObject1,
        void *edgeEndObject2) {
    stCactusEdgeEnd *edgeEnd1 = st_calloc(1, sizeof(stCactusEdgeEnd));
    stCactusEdgeEnd *edgeEnd2 = st_calloc(1, sizeof(stCactusEdgeEnd));

    edgeEnd1->node = node1;
    edgeEnd1->endObject = edgeEndObject1;
    edgeEnd1->otherEdgeEnd = edgeEnd2;

    edgeEnd2->node = node2;
    edgeEnd2->endObject = edgeEndObject2;
    edgeEnd2->otherEdgeEnd = edgeEnd1;
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
    edgeEnd->linkOrientation = 1;
    otherEdgeEnd->linkOrientation = 0;
    edgeEnd->link = otherEdgeEnd;
    otherEdgeEnd->link = edgeEnd;
}

static void stCactusEdgeEnd_setIsChainEnd(stCactusEdgeEnd *edgeEnd,
        bool isChainEnd) {
    assert(edgeEnd->link != NULL);
    edgeEnd->isChainEnd = isChainEnd;
    edgeEnd->link->isChainEnd = isChainEnd;
}

//Graph functions

stCactusGraph *stCactusGraph_construct(void *(*nodeIt)(),
        stCactusEdgeTuple(*edgeIt)(),
        void *(*mergeNodeObjects)(void *, void *), void *startNode) {
    stCactusGraph *cactusGraph = st_malloc(sizeof(stCactusGraph));
    cactusGraph->objectToNodeHash = stHash_construct2(NULL, (void (*)(void *))stCactusNode_destruct);
    return cactusGraph;
}

stCactusNode *stCactusGraph_getNode(stCactusGraph *node, void *nodeObject) {
    return stHash_search(node->objectToNodeHash, nodeObject);
}

void stCactusGraph_collapseToCactus(
        stCactusGraph *graph, void *(*mergeNodeObjects)(void *, void *),
        stCactusNode *startNode) {
    stGraph_markCycles(graph, startNode);
}

void stCactusGraph_destruct(stCactusGraph *graph) {
    stHash_destruct(graph->objectToNodeHash);
    free(graph);
}

stCactusGraphNodeIterator stCactusGraphNodeIterator_construct(
        stCactusGraph *graph) {
    stCactusGraphNodeIterator nodeIt;
    nodeIt.it = stHash_getIterator(graph->objectToNodeHash);
    nodeIt.graph = graph;
    return nodeIt;
}

stCactusNode *stCactusGraphNodeIterator_getNext(
        stCactusGraphNodeIterator *nodeIt) {
    void *key = stHash_getNext(nodeIt->it);
    if(key != NULL) {
        return stHash_search(nodeIt->graph->objectToNodeHash, key);
    }
    return NULL;
}

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIterator nodeIt) {
    stHash_destructIterator(nodeIt.it);
}

void stGraph_markCyclesP(stCactusEdgeEnd *edgeEnd,
        stHash *nodesOnChainToEdgeEnds, stList *chainPath) {
    stCactusNode *node = stCactusEdgeEnd_getNode(edgeEnd);
    stCactusEdgeEnd *edgeEnd2;
    if ((edgeEnd2 = stHash_search(nodesOnChainToEdgeEnds, node)) != NULL) { //We've traversed a cycle
        stCactusEdgeEnd_setLink(edgeEnd, edgeEnd2);
        stCactusEdgeEnd_setIsChainEnd(edgeEnd, edgeEnd2);
        for (int32_t j = stList_length(chainPath) - 2; j > 0; j -= 2) {
            stCactusEdgeEnd *edgeEnd3 = stList_get(chainPath, j);
            if (edgeEnd3 == edgeEnd) {
                break;
            }
            stCactusEdgeEnd *edgeEnd4 = stList_get(chainPath, j - 1);
            stCactusEdgeEnd_setLink(edgeEnd3, edgeEnd4);
        }
    } else {
        stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(node);
        while ((edgeEnd2 = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
            if (edgeEnd2 != edgeEnd && stCactusEdgeEnd_getLink(edgeEnd2)
                    == NULL) {
                stHash_insert(nodesOnChainToEdgeEnds, node, edgeEnd2);
                stList_append(chainPath, edgeEnd2);
                stGraph_markCyclesP(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd),
                        nodesOnChainToEdgeEnds, chainPath);
                stHash_remove(nodesOnChainToEdgeEnds, node);
                stList_pop(chainPath);
            }
        }
    }
}

void stGraph_unmarkCycles(stCactusGraph *graph) {

}

void stGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode) {
    //Set the link edges to NULL.
    stHash *nodesOnChainToEdgeEnds = stHash_construct();
    stList *chainPath = stList_construct();
    stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
        if (!stCactusEdgeEnd_isChainEnd(edgeEnd)) {
            stHash_insert(nodesOnChainToEdgeEnds, startNode, edgeEnd);
            stList_append(chainPath, edgeEnd);
            stGraph_markCyclesP(stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd),
                    nodesOnChainToEdgeEnds, chainPath);
            stHash_remove(nodesOnChainToEdgeEnds, startNode);
            stList_pop(chainPath);
        }
    }
    stHash_destruct(nodesOnChainToEdgeEnds);
    stList_destruct(chainPath);
}

void stGraph_collapseBridgesP(stCactusNode *parentNode,
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
                    stGraph_collapseBridgesP(node,
                            stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2), nodesToMerge);
                    bridges++;
                } else if (stCactusEdgeEnd_getLinkOrientation(edgeEnd2)) {
                    stGraph_collapseBridgesP(parentNode,
                            stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2), nodesToMerge);
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
            if (edgeEnd2 != edgeEnd && (stCactusEdgeEnd_getLink(edgeEnd2)
                    == NULL || stCactusEdgeEnd_getLinkOrientation(edgeEnd2))) {
                stGraph_collapseBridgesP(node,
                        stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd2), nodesToMerge);
            }
        }
    }
}

void stGraph_collapseBridges(stCactusGraph *graph, stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *)) {
    stCactusNode_edgeEndIt edgeIterator = stCactusNode_getEdgeEndIt(startNode);
    stCactusEdgeEnd *edgeEnd;
    stList *bridgesToMerge = stList_construct();
    while ((edgeEnd = stCactusNode_edgeEndIt_getNext(&edgeIterator))) {
        if (stCactusEdgeEnd_getLink(edgeEnd) == NULL
                || stCactusEdgeEnd_getLinkOrientation(edgeEnd)) {
            stGraph_collapseBridgesP(startNode,
                    stCactusEdgeEnd_getOtherEdgeEnd(edgeEnd), bridgesToMerge);
        }
    }
    //Now merge bridges
    for (int32_t i = 0; i < stList_length(bridgesToMerge); i += 2) {
        stCactusNode *parentNode = stList_get(bridgesToMerge, i);
        stCactusNode *nodeToMerge = stList_get(bridgesToMerge, i + 1);
        stCactusNode_mergeNodes(graph, parentNode, nodeToMerge, mergeNodeObjects);
    }
    stList_destruct(bridgesToMerge);
    stGraph_unmarkCycles(graph);
    stGraph_markCycles(graph, startNode);
}

