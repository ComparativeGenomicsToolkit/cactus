#ifndef CACTUS_GRAPH_H_
#define CACTUS_GRAPH_H_
#include "fastCMaths.h"
#include "commonC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
//#include "net.h"

#include "avl.h"
#include "cactus.h"
#include "commonC.h"
#include "sonLib.h"

/*
 * Data structures and methods used in building the cactus graph from
 * the pinch graph.
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus graph data structures for representing the
//basic cactus graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct CactusVertex {
    int32_t vertexID;
    struct List *edges;
};

struct CactusVertex *constructCactusVertex();

void destructCactusVertex(struct CactusVertex *vertex);

struct CactusEdge {
    struct CactusVertex *from;
    struct CactusVertex *to;
    struct List *pieces;
    struct CactusEdge *rEdge;
};

struct CactusEdge *constructCactusEdge(struct List *pieces);

void destructCactusEdge(struct CactusEdge *edge);

struct CactusGraph {
    struct List *vertices;
};

struct CactusGraph *constructCactusGraph(struct PinchGraph *pinchGraph,
        struct List *threeEdgeConnectedComponents);

void destructCactusGraph(struct CactusGraph *cactusGraph);

void checkCactusGraph(struct PinchGraph *pinchGraph,
        struct List *threeEdgeConnectedComponents,
        struct CactusGraph *cactusGraph);

void checkCactusContainsOnly2EdgeConnectedComponents(
        struct CactusGraph *cactusGraph);

int32_t cactusGraph_getEdgeNumber(struct CactusGraph *cactusGraph);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to calculate 2-edge connected components
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct List *computeBiConnectedComponents(struct CactusGraph *cactusGraph);

/*
 * Method sorts the bi-connected components into order.
 */
struct List
        *computeSortedBiConnectedComponents(struct CactusGraph *cactusGraph);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Method to compute DFS discovery time for vertices.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t *getDFSDiscoveryTimes(struct CactusGraph *cactusGraph);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//I/O Methods to interact with the 3-edge connected component
//code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

stList *writeOut3EdgeGraph(struct PinchGraph *pinchGraph,
        struct List *greyEdgeComponents);

struct List *readThreeEdgeComponents(struct PinchGraph *,
        struct List *greyEdgeComponents, stList *threeEdgeComponents);

void writeOutCactusGraph(struct CactusGraph *cactusGraph,
        struct PinchGraph *pinchGraph, FILE *fileHandle);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Script to construct the cactus graph
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void computeCactusGraph(struct PinchGraph *pinchGraph,
        struct CactusGraph **cactusGraph,
        struct List **threeEdgeConnectedComponents,
        int32_t excludeDegree1Edges);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to manipulate the cactus graph
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void circulariseStems(struct CactusGraph *cactusGraph, struct PinchGraph *pinchGraph);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which blocks in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

stSortedSet *filterBlocksByTreeCoverageAndLength(
        struct List *biConnectedComponents, Flower *net,
        float minimumTreeCoverage, /*Minimum tree coverage to be included */
        int32_t minimumBlockDegree, /*The minimum number of segments in a block to be included (>=)*/
        int32_t minimumBlockLength, /*The minimum length of an block to be included */
        int32_t minimumChainLength, /* Minimum chain length to be included */
        struct PinchGraph *pinchGraph);

void logTheChosenBlockSubset(struct List *biConnectedComponents,
        stSortedSet *chosenBlocks, struct PinchGraph *pinchGraph, Flower *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus graph misc functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *cactusEdgeToFirstPinchEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph);

/*
 * Returns non-zero if the edge is a stub (i.e. not a block end).
 */
int32_t isAStubCactusEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph);

/*
 * Returns non-zero if edge is stub end which is free (the dead end is not attached to the source vertex).
 */
int32_t isAFreeStubCactusEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph);

struct hashtable *createHashColouringPinchEdgesByChains(
        struct PinchGraph *pinchGraph, struct List *biConnectComponentsList);

#endif
