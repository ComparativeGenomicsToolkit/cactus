/*
 * stCactusGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#ifndef ST_CACTUS_GRAPH_H_
#define ST_CACTUS_GRAPH_H_

typedef struct _stCactusNode stCactusNode;

typedef struct _stCactusEdgeEnd stCactusEdgeEnd;

typedef struct _stCactusNodeEdgeEndIt {
    stCactusEdgeEnd *edgeEnd;
} stCactusNodeEdgeEndIt;

typedef struct _stCactusGraph stCactusGraph;

typedef struct _stCactusGraphNodeIterator {
    stHashIterator *it;
    stCactusGraph *graph;
} stCactusGraphNodeIt;

//Node functions

stCactusNode *stCactusNode_construct(stCactusGraph *graph,
        void *nodeObject);

void *stCactusNode_getObject(stCactusNode *node);

stCactusNodeEdgeEndIt stCactusNode_getEdgeEndIt(stCactusNode *node);

stCactusEdgeEnd *stCactusNodeEdgeEndIt_getNext(stCactusNodeEdgeEndIt *it);

stCactusEdgeEnd *stCactusNode_getFirstEdgeEnd(stCactusNode *node);

//Edge functions

stCactusEdgeEnd *stCactusEdgeEnd_construct(stCactusGraph *graph,
        stCactusNode *node1, stCactusNode *node2, void *edgeEndObject1,
        void *edgeEndObject2);

void *stCactusEdgeEnd_getObject(stCactusEdgeEnd *edgeEnd);

stCactusNode *stCactusEdgeEnd_getNode(stCactusEdgeEnd *edgeEnd);

stCactusNode *stCactusEdgeEnd_getOtherNode(stCactusEdgeEnd *edgeEnd);

stCactusEdgeEnd *stCactusEdgeEnd_getOtherEdgeEnd(stCactusEdgeEnd *edgeEnd);

stCactusEdgeEnd *stCactusEdgeEnd_getLink(stCactusEdgeEnd *edgeEnd);

bool stCactusEdgeEnd_getLinkOrientation(stCactusEdgeEnd *edgeEnd);

bool stCactusEdgeEnd_isChainEnd(stCactusEdgeEnd *edgeEnd);

stCactusEdgeEnd *stCactusEdgeEnd_getNextEdgeEnd(stCactusEdgeEnd *edgeEnd);

//Graph functions

stCactusGraph *stCactusGraph_construct2(void (*destructNodeObjectFn)(void *), void (*destructEdgeEndObjectFn)(void *));

stCactusGraph *stCactusGraph_construct(void);

void stCactusGraph_collapseToCactus(
        stCactusGraph *graph, void *(*mergeNodeObjects)(void *, void *), stCactusNode *startNode);

stCactusNode *stCactusGraph_getNode(stCactusGraph *node, void *nodeObject);

void stCactusGraph_destruct(stCactusGraph *graph);

stCactusGraphNodeIt *stCactusGraphNodeIterator_construct(stCactusGraph *graph);

stCactusNode *stCactusGraphNodeIterator_getNext(stCactusGraphNodeIt *nodeIt);

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIt *);

void stCactusGraph_unmarkCycles(stCactusGraph *graph);

void stCactusGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode);

void stCactusGraph_collapseBridges(stCactusGraph *graph,
        stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *));

#endif /* ST_CACTUS_GRAPH_H_ */
