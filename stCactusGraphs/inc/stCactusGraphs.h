/*
 * stCactusGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#ifndef ST_CACTUS_GRAPH_H_
#define ST_CACTUS_GRAPH_H_

typedef struct _stCactusNode stCactusNode;

typedef struct _stCactusNode_edgeEndIt stCactusNode_edgeEndIt;

typedef struct _stCactusEdgeEnd stCactusEdgeEnd;

typedef struct _stCactusGraph stCactusGraph;

typedef struct _stCactusGraphNodeIterator stCactusGraphNodeIterator;

typedef struct _stCactusEdgeTuple stCactusEdgeTuple;

//Node functions

stCactusNode *stCactusNode_construct(stCactusGraph *graph,
        void *nodeObject);

void *stCactusNode_getObject(stCactusNode *node);

stCactusNode_edgeEndIt stCactusNode_getEdgeEndIt(stCactusNode *node);

stCactusEdgeEnd *stCactusNode_edgeEndIt_getNext(stCactusNode_edgeEndIt *it);

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

//Graph functions

stCactusGraph *stCactusGraph_construct();

void stCactusGraph_collapseToCactus(
        stCactusGraph *graph, void *(*mergeNodeObjects)(void *, void *), stCactusNode *startNode);

stCactusNode *stCactusGraph_getNode(stCactusGraph *node, void *nodeObject);

void stCactusGraph_destruct(stCactusGraph *graph);

stCactusGraphNodeIterator stCactusGraphNodeIterator_construct(stCactusGraph *graph);

stCactusNode *stCactusGraphNodeIterator_getNext(stCactusGraphNodeIterator *nodeIt);

void stCactusGraphNodeIterator_destruct(stCactusGraphNodeIterator);

void stGraph_unmarkCycles(stCactusGraph *graph);

void stGraph_markCycles(stCactusGraph *graph, stCactusNode *startNode);

void stCactusGraph_collapseBridges(stCactusGraph *graph,
        stCactusNode *startNode, void *(*mergeNodeObjects)(void *, void *));

#endif /* ST_CACTUS_GRAPH_H_ */
