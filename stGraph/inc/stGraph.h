/*
 * stGraph.h
 *
 *  Created on: 14 Apr 2012
 *      Author: benedictpaten
 */

#ifndef STGRAPH_H_
#define STGRAPH_H_

typedef struct _stNode {
        stEdgeEnd *head;
        stEdgeEnd *tail;
        void *nodeObject;
} stNode;

typedef struct _stNodeIt {
        stNode *node;
} stNodeIt;

typedef struct _stEdgeEnd {
        stEdgeEnd *otherEdgeEnd;
        stEdgeNode *node;
        stEdgeNode *nEdgeEnd;
        void *endObject;
} stEdgeEnd;

typedef struct _stGraph {
        stSortedSet *nodes;
} stGraph;

//Node functions

stNode *stNode_construct(stGraph *graph, void *nodeObject);

void stNode_destruct(stNode *node);

void *stNode_getObject(stNode *node);

stNodeIt stNode_getIt(stNode *node);

stEdgeEnd *stNodeIt_getNextEnd(stNodeIt *node);

stNode *stNodeIt_getNextNode(stNodeIt *node);

stNode *stNode_mergeNodes(stGraph *graph, stNode *node1, stNode *node2);

//Edge functions

stEdgeEnd *stEdgeEnd_construct(stGraph *graph, stNode *node1, stNode *node2, void *edgeEndObject1, void *edgeEndObject2);

void *stEdgeEnd_getObject(stEdgeEnd *edgeEnd);

stNode *stEdgeEnd_getNode(stEdgeEnd *edgeEnd);

stNode *stEdgeEnd_getOtherNode(stEdgeEnd *edgeEnd);

stEdgeEnd *stEdgeEnd_getOtherEnd(stEdgeEnd *edgeEnd);

stNode *stEdgeEnd_contractEdge(stGraph *graph, stEdgeEnd *edgeEnd);

void stEdgeEnd_reconnectEnd(stGraph *graph, stEdgeEnd *edgeEnd, stNode *newNode);

void stEdgeEnd_setLabelOrientation(stEdgeEnd *edgeEnd, bool inOrOut);

bool stEdgeEnd_getLabelOrientation(stEdgeEnd *edgeEnd);

//Graph functions

stGraph *stGraph_construct();

void stGraph_destruct(stGraph *graph);

void stGraph_collapseThreeEdgeConnectedComponents(stGraph *graph, void *(*mergeNodeObjects)(void *, void *));

/*
 * For each node n and unlabelled edge e_1, dp dfs, when cycle encountered label
 * internal node ends with 1 label source node ends with 0.
 */

void stGraph_labelEndsP(stEdgeEnd *edgeEnd, stHash *nodesOnChainToEdgeEnds) {
    stNode *node = stEdgeEnd_getNode(edgeEnd);
    stEdgeEnd *edgeEnd2;
    if((edgeEnd2 = stHash_search(nodesOnChainToEdgeEnds, node)) != NULL) { //We've traversed a cycle
        stEdgeEnd_setLabelOrientation(edgeEnd, 1);
        stEdgeEnd_setLabelOrientation(edgeEnd2, 2);
    }
    else {
        stNodeIt edgeIterator = stNode_getIt(node);
        while((edgeEnd2 = stNodeIt_getNextNode(&edgeIterator))) {
            if(edgeEnd2 != edgeEnd && stEdgeEnd_getLabelOrientation(edgeEnd2) == 0) {
                stHash_insert(nodesOnChainToEdgeEnds, node, edgeEnd2);
                stGraph_labelEndsP(stEdgeEnd_getOtherEnd(edgeEnd), nodesOnChainToEdgeEnds);
                stHash_remove(nodesOnChainToEdgeEnds, node);
            }
        }
    }
}

void stGraph_labelEnds(stGraph *graph, stNode *startNode) {
    //First reset graph by labelling all ends 0


    stHash *nodesOnChainToEdgeEnds = stHash_construct();
    stNodeIt edgeIterator = stNode_getIt(node);
    stEdgeEnd *edgeEnd;
    while((edgeEnd = stNodeIt_getNextNode(&edgeIterator))) {
        if(!stEdgeEnd_getLabelOrientation(edgeEnd)) {
            stHash_insert(nodesOnChainToEdgeEnds, startNode, edgeEnd);
            stGraph_labelEndsP(stEdgeEnd_getOtherEnd(edgeEnd), nodesOnChainToEdgeEnds);
            stHash_remove(nodesOnChainToEdgeEnds, startNode);
        }
    }
    stHash_destruct(nodesOnChainToEdgeEnds);
}

void stGraph_collapseBridgesAndLabelEndsP(stNode *rootNode, stNode *node) {

}

void stGraph_collapseBridgesAndLabelEnds(stGraph *graph, stNode *startNode) {
    //Label ends
    stGraph_labelEnds(graph, startNode);

    //Identify bridges
    stNodeIt edgeIterator = stNode_getIt(node);
    stEdgeEnd *edgeEnd;
    while((edgeEnd = stNodeIt_getNextNode(&edgeIterator))) {
        stGraph_collapseBridgesAndLabelEndsP(startNode, stEdgeEnd_getOtherEnd(edgeEnd));
    }
}


#endif /* STGRAPH_H_ */
